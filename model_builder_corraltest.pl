#!/usr/bin/perl
###Corral version
use warnings;
use strict;
use lib "/home/angus/Documents/work/text_analyzer";
use Util;
use Chart::Gnuplot;
use Time::HiRes qw(gettimeofday clock);

sub build_mod;
sub binning;
sub synthesize;

#Fetch arguments
#getting flag values 

my ($xdatas,$ydatas,$pdatas) = ([],[],[]); #freqs, D-values, p-values
#my @xdata;#general purpose debugging plotter;
#my @ydata;
#my @pdata;

my $PWD = `pwd`;
chomp $PWD;

if(scalar @ARGV ==0 || scalar @ARGV >2){
	print "Usage:model_builder.pl input [output]\n";
	exit;
}

my $ifname = shift @ARGV or die "Please enter the input file name\n" ;
my $ofname = (shift @ARGV);
if (not defined $ofname){$ofname = $ifname;}
open INPUT, $ifname or die "can't open the file\n";
open OUTPUT, (">".$ofname."log") or die "can't write to file.\n";

my $upperbound = 0;
my $multimodel=0;
my $cont = 0;
print ("Analyze with upper bound?[n]");
if(<STDIN> =~ /^[Yy][Ee]{0,1}[Ss]{0,1}\n$/){
	$upperbound = 1;
	print ("Build multiple model?[n]");
	$multimodel = <STDIN> =~ /^[Yy][Ee]{0,1}[Ss]{0,1}\n$/;
}
my $pthreshold=0.05;
print ("What's the threshold of p-value?[0.05]");
if(<STDIN> =~ /^(0.[\d]+)\n$/){
	$pthreshold = $1 +0;
}
my $earlytermination=0;
print ("Allow Early Termination?[n]");
$earlytermination = <STDIN> =~ /^[Yy][Ee]{0,1}[Ss]{0,1}\n$/;


#----------getting frequency file--------------
print OUTPUT "\n#----------getting frequency file--------------\n";
my %freq2count;
my @freqarr;
my @countarr;
my $wt0 = gettimeofday();
my $ct0=clock();

while(<INPUT>){
	chomp;
	next if $_ eq "";	
	/(.+)+:(.+)/;
	if ($1 eq "Frequency"){
		@freqarr = split(/\s/, $2);
	}
	elsif ($1 eq "# of words"){
		@countarr = split(/\s/, $2);
		for my $i(0..$#freqarr){
			$freq2count{$freqarr[$i]}=$countarr[$i];
		}	
	}
	else{die "Error: unknown line title";}
}
delete $freq2count{""};

@freqarr = sort {($a+0) <=> ($b+0)} keys %freq2count;
@countarr=();
for my $i (0..$#freqarr){
	push @countarr, $freq2count{$freqarr[$i]};
	if ($freqarr[$i]!=int($freqarr[$i])){$cont = 1;}
}

my $wt1 = gettimeofday();
my $ct1=clock();
print OUTPUT "Frequencies:\t", scalar @freqarr, "\n";
print OUTPUT "Used Time:\tWall=", $wt1-$wt0, ", CPU=", $ct1-$ct0, "\n";


#TODO --------------Binning---------------------might effectively reduce noise

sub binning{
	my $freqarr = shift;
	my $countarr = shift;
	my $bmin = shift;
	my $size = shift;
	my $barr = shift;
	my $bcarr = shift;
	my @data_b;
	my $lower = $bmin;
	my $upper = $lower*$size;
	my $i=0;
	my $current_bc=0;
	while($i<scalar @{$freqarr}){
		if($lower<=$freqarr->[$i] and $freqarr->[$i] < $upper){
			$current_bc+=$countarr->[$i];
			$i++;
		}
		else{
			push @{$barr}, $lower;
			push @{$bcarr}, $current_bc;
			$current_bc = 0;
			$lower = $upper;
			$upper*=$size;
		}
	}
	push @{$barr}, $lower;
	push @{$bcarr}, $current_bc;
	
}


#------------fitting---------------------------
print OUTPUT "\n#----------Constructng Models--------------\n";
$wt0 = gettimeofday();
$ct0 = clock();
my $modArr=[];
my $B=$#freqarr;
while(1){
	$upperbound = 1 if ($B != $#freqarr);
	my @para=build_mod(\@freqarr, \@countarr, $B);
	if($para[0] != -1){
		printf "\nrange= %u-%u, b/a= %.2f, sum=%u\nxmin= %f, xmax= %f, a= %.3f, KS=%.5f\n", 
			$para[4],$para[5],$para[1]/$para[0],$para[6],$para[0],$para[1]+0,$para[2],$para[3];
		printf OUTPUT "\nrange= %u-%u(%u), b/a= %.2f, sum=%u\nxmin= %f, xmax= %f, a= %.3f, KS=%.5f\n", 
			$para[4],$para[5],$para[5]-$para[4]+1,$para[1]/$para[0],$para[6],$para[0],$para[1]+0,$para[2],$para[3];
		push @{$modArr}, \@para;
	}
	$B = $para[4];
	if($para[0] == -1 && $earlytermination){$B = int($B*0.9);}
	$B--;
	last if $B<=1;
	last if (!$multimodel);
}
$wt1 = gettimeofday();
$ct1 = clock();
print OUTPUT scalar @{$modArr}, " models constructed.\n";
print OUTPUT "Used Time:\tWall=", $wt1-$wt0, ", CPU=", $ct1-$ct0, "\n";

my @checklist = sort {($b->[1]/$b->[0]) <=> ($a->[1]/$a->[0])} @{$modArr};

#--------------plotting------------------------------


my $N=0;
for my $i(reverse 0..$#freqarr){
	$N+=$countarr[$i];
}

my @Narr;
my @esums;
for (@{$modArr}){
	my $mod = $_;
	my $Xm = $mod->[0];
	my $XM = $mod->[1];
	push @Narr, 0;
	push @esums,0;
	for my $i(0..$#freqarr){
		if ($freqarr[$i]<$Xm) 	{next;}
		elsif($freqarr[$i]<=$XM){$Narr[-1]+=$countarr[$i];}
		else 					{$esums[-1]+=$countarr[$i];}
	}
}
my @tccdf;
my @tpdf;
my @tfreq;
#my $x = $freqarr[-1];
#my $pmod = 0;
#my $mod = $modArr->[$pmod];
#my $Xm = $mod->[0];
#my $XM = $mod->[1];
#my $alpha = $mod->[2];
	
for my $pmod(0..(scalar @{$modArr}-1)){
	my $mod = $modArr->[$pmod];
	my $Xm = $mod->[0];
	my $XM = $mod->[1];
	my $alpha = $mod->[2];
	my $l_limit = $mod->[4];
	my $u_limit = $mod->[5];
	if(!$cont){
		my @xs = $freqarr[$l_limit]..$freqarr[$u_limit];
		for my $x(reverse @xs){
			push @tfreq, log($x)/log(10);
			push @tccdf, log($Narr[$pmod]/$N*(zeta($x,$XM ,$alpha)/zeta($Xm,$XM,$alpha))+$esums[$pmod]/$N)/log(10);
			push @tpdf, $Narr[$pmod]/$N*($x**(-$alpha)/zeta($Xm,$XM ,$alpha));
		}	
	}
	else{
		for my $i(reverse $l_limit..$u_limit){
			my $x = $freqarr[$i];
			if ($x==$XM) {next;}#skip the last point, otherwise the logarithm plot can't be drawn.
			push @tfreq, log($x)/log(10);
			push @tccdf, log( $Narr[$pmod]/$N*($x**(1-$alpha)-$XM**(1-$alpha))/($Xm**(1-$alpha)-$XM**(1-$alpha))+$esums[$pmod]/$N )/log(10);
			push @tpdf, $Narr[$pmod]/$N*($x**(-$alpha)/($Xm**(1-$alpha)-$XM**(1-$alpha)) *($alpha-1));
		}
#		my @xs = @freqarr;
#		for my $x(reverse @xs){
#			push @tfreq, log($x)/log(10);
#			push @tccdf, log( $Narr[$pmod]/$N*($x**(1-$alpha)-$XM**(1-$alpha))/($Xm**(1-$alpha)-$XM**(1-$alpha))+$esums[$pmod]/$N )/log(10);
#			push @tpdf, $Narr[$pmod]/$N*($x**(-$alpha)/($Xm**(1-$alpha)-$XM**(1-$alpha)) *($alpha-1));
#		}	
	}
}
#for my $i(2..scalar @freqarr){
#	my $x = $freqarr[-$i];
#	push @tfreq, log($freqarr[-$i])/log(10);
#	push @tccdf, log($Narr[$pmod]/$N*(zeta($x,$XM ,$alpha)/zeta($Xm,$XM,$alpha))+$tsum/$N)/log(10);
#	push @tpdf, $Narr[$pmod]/$N*($x**(-$alpha)/zeta($Xm,$XM ,$alpha));
#	if ($x==$Xm){
#		$tsum += $Narr[$pmod];
#		$pmod++;
#		last if $pmod >scalar @{$modArr}-1;
#		$mod = $modArr->[$pmod];
#		$Xm = $mod->[0];
#		$XM = $mod->[1];
#		$alpha = $mod->[2];
#	}
#	
#}

#@tfreq = reverse @tfreq;
#@tccdf = reverse @tccdf;
#@tpdf = reverse @tpdf;

my @logtpdf;
for (@tpdf){
	push @logtpdf, log($_)/log(10);
}

my $tpdfdata = Chart::Gnuplot::DataSet->new(
	xdata	=> \@tfreq,
	ydata	=> \@logtpdf,
	style	=> "lines",
	color	=> "#007f00",
	linetype=> "solid",
	width	=> 3
);

my $sum_e=0;
for my $i(0..$#countarr){
	$sum_e+=$countarr[$i];
}

my @epdf;
for (@countarr){
	push @epdf, log($_/$sum_e)/log(10);
}


my @logfreq;
for my $i(0..$#freqarr){
	push @logfreq, log($freqarr[$i])/(log 10);
}

my $epdfdata= Chart::Gnuplot::DataSet->new(
	xdata => \@logfreq,
	ydata => \@epdf,
	color => "#000000",
	pointtype => 'plus',
	pointsize => 2
);

my $pdfchart = Chart::Gnuplot->new(
	output	=> "$ofname"."_pdf.png",
	title	=> "PDF",
	xlabel	=> "freqency",
	ylabel	=> "PDF",
	bg		=> "white",
	grid	=>{linetype=>'dash'},
);
if(scalar @tfreq !=0)	{$pdfchart->plot2d($epdfdata, $tpdfdata); }
else					{$pdfchart->plot2d($epdfdata); }

my @logxs;
my (@Ddatasets, @Pdatasets);
for (0.. scalar @{$xdatas}-1){
	my @logx;
	for my $i(0..scalar @{$xdatas->[$_]} -1){
		push @logx, (log $xdatas->[$_]->[$i])/(log 10);
	}
	push @logxs, \@logx;
	push @Ddatasets, Chart::Gnuplot::DataSet->new(
		xdata =>	\@logx,
		ydata =>	$ydatas->[$_]
	);
	push @Pdatasets, Chart::Gnuplot::DataSet->new(
		xdata =>	\@logx,
		ydata =>	$pdatas->[$_]
	);
}
#my @logx;
#for my $i (0..$#xdata){
#	push @logx, (log $xdata[$i])/(log 10);
#}

#my @logy;
#for my $i (0..$#ydata){
#	push @logy, (log $ydata[$i])/(log 10);
#}
#print scalar @ydata, scalar @logx, "\n";
#my $debugdata = Chart::Gnuplot::DataSet->new(
#	xdata => \@logx,
#	ydata => \@ydata,
#	color => "#ff0000",
#	pointtype => 'plus',
#	pointsize => 2
#	
#);
my $KSchart = Chart::Gnuplot->new(
	output	=> "$ofname"."_D.png",
	title	=> "D-value",
	xlabel	=> "frequency",
	ylabel	=> "D-value",
	bg		=> "white",
	grid	=>{linetype=>'dash'},
);
#$KSchart->plot2d($debugdata);
$KSchart->plot2d(@Ddatasets);

#my $pDataSet = Chart::Gnuplot::DataSet->new(
#	xdata => \@logx,
#	ydata => \@pdata,
#);
my $PChart = Chart::Gnuplot->new(
	output	=> "$ofname"."_P.png",
	title	=> "P-value",
	xlabel	=> "frequency",
	ylabel	=> "P-value",
	bg		=>	"white",
	grid	=>{linetype=>'dash'},
);
#$PChart->plot2d($pDataSet);
$PChart->plot2d(@Pdatasets);
my @eccdf;
my $s_e=0;
for my $i(1..scalar @freqarr){
	$s_e+=$countarr[-$i];
	push @eccdf, log($s_e/$sum_e)/log(10);
}
@eccdf = reverse @eccdf;
my $eccdfdata = Chart::Gnuplot::DataSet->new(
	xdata	=> \@logfreq,
	ydata	=> \@eccdf
);
my $tccdfdata = Chart::Gnuplot::DataSet->new(
	xdata	=> \@tfreq,
	ydata	=> \@tccdf,
	style	=> "lines",
	color	=> "#007f00",
	linetype=> "solid",
	width	=> 3
);

my $ccdfchart = Chart::Gnuplot->new(
	output	=> "$ofname"."_ccdf.png",
	title	=> "CCDF",
	xlabel	=> "frequency",
	ylabel	=> "CCDF",
	bg		=> "white",
	grid	=>{linetype=>'dash'},
);
if(scalar @tfreq !=0)	{$ccdfchart->plot2d($eccdfdata, $tccdfdata);}
else					{$ccdfchart->plot2d($eccdfdata);}

#my @KSindex = 0..$#xdata;
exit;
#---------goodness-of-fit----------------------
for (@checklist){
	print "===============Starting goodness-of-fit======================\n";
	print OUTPUT "\n-------------Starting goodness-of-fit---------------\n";
	$wt0 = gettimeofday();
	$ct0 = clock();
	my $model=$_;
	my $n = 2500;
	my $bt=0;#better than
	print "Xmin=",$model->[0], "\tXmax=",$model->[1],"\ta=",$model->[2],"\tKS=",$model->[3],"\n";
	print "Number of set of sythetic data generated = $n\n";
	print OUTPUT "Xmin=",$model->[0], "\tXmax=",$model->[1],"\ta=",$model->[2],"\tKS=",$model->[3],"\n";
	print OUTPUT "Number of set of sythetic data generated = $n\n";
	my @As;
	my @Ds;
	for my $i(1..$n){
		my @freqarr_s;
		my @countarr_s;
		synthesize(\@freqarr, \@countarr, $model, \@freqarr_s, \@countarr_s);
		print "$i: ";
		my @model_s=build_mod(\@freqarr_s, \@countarr_s, $#freqarr_s, $model->[3], $model->[2]);
		print "Xmin=",$model_s[0], "\tXmax=",$model_s[1],"\ta=",$model_s[2],"\tKS=",$model_s[3],"\n";
		push @As, $model_s[2];
		push @Ds, $model_s[3];
		if ($model->[3]<$model_s[3]){$bt++};
	}
	my $A_avg = 0;
	my $D_avg = 0;
	for my $i(0..$#As){
		$A_avg+=$As[$i];
		$D_avg+=$Ds[$i];
	}
	$A_avg/=$n;
	$D_avg/=$n;
	my $A_std=0;
	my $D_std=0;
	for my $i(0..$#As){
		$A_std+=($As[$i]-$A_avg)*($As[$i]-$A_avg);
		$D_std+=($Ds[$i]-$D_avg)*($Ds[$i]-$D_avg);
	}
	$A_std = ($A_std/($n-1))**0.5;
	$D_std = ($D_std/($n-1))**0.5;
	####These data are incorrect due to early termination.######
	print "alpha avg=$A_avg, std=$A_std\tKS avg=$D_avg, std=$D_std\n";
	print OUTPUT "alpha avg=$A_avg, std=$A_std\tKS avg=$D_avg, std=$D_std(data incorrect due to early termination)\n";

	
	
	
	my $p_value=$bt/$n;
	print "p-value = $p_value\n";
	print OUTPUT "\np-value = $p_value\n";
	$wt1 = gettimeofday();
	$ct1 = clock();
	print OUTPUT "Used Time:\tWall=", $wt1-$wt0, ", CPU=", $ct1-$ct0, "\n";
}


#===================subroutines========================
#TODO:The program is impractical if the speed of this method is not improved.
sub build_mod{
	my $freqarr = shift;
	my $countarr = shift;
	my $u_limit=shift;  #$freqarr[$u_limit] = Xmax
	#$u_limit=$#freqarr if (not defined $u_limit);
	my $sum =0;
	my $lnsum=0;
	my $alpha;
	my $xmin='inf';
	my $xmax=$freqarr->[$u_limit]; $xmax = 'inf' if (!$upperbound);
	my $Dmin='inf';
	my $P;
	my $summin=0;
	my $l_limit=$u_limit;
	#my @rev_CDF;
	my @_Aarr;	
	my @_sumarr;
	my ($zcount,$freepass)=(0,0);
	my ($xdata,$ydata,$pdata) = ([],[],[]); #freqs, D-values, p-values
	for my $n(0..$u_limit){
		my $_xmin=$freqarr->[$u_limit-$n];
		#print "_xmin=$_xmin\t";
		$sum+=$countarr->[$u_limit-$n];
		$lnsum+=$countarr->[$u_limit-$n]*log $_xmin;
		#next if ($n == 0 ||$n ==1);
		my $_alpha;
		if(!$cont){$_alpha=find_alpha($_xmin, $xmax, $sum, $lnsum);}  #TODO don't know the MLE of alpha with upper bound
		else {$_alpha = 1+$sum*(1/($lnsum-$sum*log($_xmin))) }
		if(!$cont && zeta($_xmin, $xmax, $_alpha)==0){ next;}
		if($cont && $_xmin**(1-$_alpha)-$xmax**(1-$_alpha) ==0) {next;}
		#next if ( zeta($_xmin,$xmax, $_alpha) ==0 );

		#TODO: is there a way to avoid this inner for-loop?
		my $D;
		if($cont){$D = KS_cont($freqarr, $countarr, $u_limit-$n,$u_limit,$_alpha,$sum,!$upperbound);}
		else{$D=KS($freqarr, $countarr, $u_limit-$n,$u_limit,$_alpha,$sum,!$upperbound);}
		
		
		push @{$xdata}, $_xmin;
		push @{$ydata}, $D;
		
		#calculating p-value
		my ($p, $numofdata)=(0,1);
		if ($earlytermination && $freepass>0){$freepass--;$p=0;}
		else{
			for (1..$numofdata){
				my $sdata = plrand_cont($_xmin, $xmax,$_alpha,$sum);
				my %f2c;
				for (@{$sdata}){
					if (exists $f2c{$_}){$f2c{$_}++;}
					else {$f2c{$_}=1;}
				}
				my @sfreqs = sort {($a+0)<=>($b+0)} keys(%f2c);
				my @scounts;
				my ($ssum,$slnsum)=(0,0);
				for (@sfreqs){
					push @scounts, $f2c{$_};
					$ssum+=$f2c{$_};
					$slnsum+=$f2c{$_}* log $_;
				}
				my $salpha;
				if (!$cont){$salpha=find_alpha($sfreqs[0], $xmax, $ssum, $slnsum);}
				else { $salpha = 1+$ssum*(1/($slnsum-$ssum*log($sfreqs[0]))); }
				my $sD=KS_cont(\@sfreqs,\@scounts,0,$#sfreqs,$salpha,$ssum,!$upperbound,$D);
				$p++ if ($D<$sD);
				#if($earlytermination && $p/$numofdata>$pthreshold){last;}
			}
			$p/=$numofdata;
			if ($p != 0){$zcount = 0; $freepass = 0;}
			else{$freepass = $zcount++;}
			print "xmax=$xmax\t_xmin=$_xmin\ta=$_alpha\tD=$D\tp=$p\n";
		}
		push @{$pdata},$p;
		#print "D=$D\tDmin=$Dmin\n";
		if($p>$pthreshold){
			$Dmin=$D;
			$xmin=$_xmin;
			$alpha=$_alpha;
			$l_limit=$u_limit-$n;
			$summin=$sum;
			$P=$p;
		}
	}
	if($l_limit != $u_limit || $u_limit == @{$freqarr}-1)#print the first model or whenever a probable model is found.
	{
		push @{$xdatas}, $xdata;
		push @{$ydatas}, $ydata;
		push @{$pdatas}, $pdata;
	}
	if($l_limit == $u_limit){$xmin = -1;}#indicating that no possible model found with this xmax
	
	return my @para=($xmin, $xmax, $alpha, $Dmin, $l_limit,$u_limit, $summin, $P);
}


sub synthesize{
	my $freqarr=shift;
	my $countarr=shift;
	my $mod=shift;
	my $freqarr_s=shift;
	my $countarr_s=shift;
	my $Xmin=$mod->[0];
	my $Xmax=$mod->[1];
	my $alpha=$mod->[2];
	my $ct_less=0;
	my $ct_more=0;
	my @headfreq;
	my @headcount;
	
	for my $i(0..scalar @{$freqarr}-1){
		if($freqarr->[$i]<$Xmin)		{$ct_less+=$countarr->[$i];
									push @headfreq, $freqarr->[$i];
									push @headcount, $countarr->[$i];}
		elsif($freqarr->[$i]<=$Xmax)	{$ct_more+=$countarr->[$i];}
		else{last;}
	}
	
	my $nless=0;
	my $nmore=0;
	for my $i(1..($ct_less+$ct_more) ){
		my $r=rand;
		if($r<=$ct_less/($ct_less+$ct_more))	{$nless++;}
		else									{$nmore++;}
	}
	my $rddata_1= sample(\@headfreq, \@headcount, $nless);
	my $rddata_2= plrand($Xmin, $Xmax, $alpha, $nmore);
	#--the two hashes should not be merged, so that their keys can be sorted faster
	my %freq2count_1;
	for(@{$rddata_1}){
		if(exists $freq2count_1{$_}){
			$freq2count_1{$_}++;
		}
		else{
			$freq2count_1{$_}=1;
		}
	}
	my @freqarr_1=sort {$a <=> $b} keys %freq2count_1;
	for (@freqarr_1){
		push @{$freqarr_s}, $_;
		push @{$countarr_s}, $freq2count_1{$_};
	}
	
	my %freq2count_2;
	for(@{$rddata_2}){
		if(exists $freq2count_2{$_}){
			$freq2count_2{$_}++;
		}
		else{
			$freq2count_2{$_}=1;
		}
	}
	my @freqarr_2=sort {$a <=> $b} keys %freq2count_2;
	for (@freqarr_2){
		push @{$freqarr_s}, $_;
		push @{$countarr_s}, $freq2count_2{$_};
	}
	
}
