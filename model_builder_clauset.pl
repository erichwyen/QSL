#!/usr/bin/perl
###Clauset version
use warnings;
use strict;
use lib "/home/angus/mylib";
use Util;
use Chart::Gnuplot;
use Time::HiRes qw(gettimeofday clock);

sub build_mod;
sub synthesize;

#Fetch arguments
#getting flag values 

my ($xdatas,$ydatas,$pdatas) = ([],[],[]); #freqs, D-values, p-values

my $PWD = `pwd`;
chomp $PWD;
my $ofname;
my $stdout;
my $outputdir = $PWD;
my $printSynthetic=0;
##getting parameters
my $i=0;
while($i < scalar @ARGV){
	if($ARGV[$i] eq '-o'){
		$ofname = $ARGV[$i+1];
		if(not defined $ofname) {die "Usage: -o <output name>\n";}
		splice(@ARGV, $i, 2);
	}
	elsif($ARGV[$i] eq '-d'){
		$outputdir = $ARGV[$i+1];
		if(not defined $outputdir){die "Usage: -d <output dirctory>\n";}
		splice(@ARGV, $i, 2);
		if(not -d $outputdir){`mkdir $outputdir`}
	}
	elsif($ARGV[$i] eq '--printSynthetic'){
		$printSynthetic = 1;
		splice(@ARGV, $i, 1);
	}
	elsif($ARGV[$i] eq '-O'){
		$stdout = $ARGV[$i+1];
		if(not defined $stdout) {die "Usage: -O <output stream name>\n";}
		splice(@ARGV, $i, 2);
	}	
	else{$i++;}
}
undef $i;

my $ifname = shift @ARGV or die "Please enter the input file name\n" ;
if (not defined $ofname){$ofname = $ifname; $ofname =~s/\..+//; $ofname=~s/.+\///;}
$ofname = "$outputdir/".$ofname;
open INPUT, $ifname or die "can't open the file\n";
open OUTPUT, (">".$ofname.".log") or die "can't write to file.\n";


my $upperbound = 0;
my $multimodel=0;
print ("Analyze with upper bound?[n]");
$upperbound = <STDIN> =~ /^[Yy][Ee]{0,1}[Ss]{0,1}\n$/;
print ("Build multiple models?[n]");
$multimodel = <STDIN> =~ /^[Yy][Ee]{0,1}[Ss]{0,1}\n$/;
my $pthreshold=0.05;
print ("What's the threshold of p-value?[0.05]");
if(<STDIN> =~ /^(0.[\d]+)\n$/){
	$pthreshold = $1 +0;
}
my $earlytermination=0;
print ("Allow Early Termination?[n]");
$earlytermination = <STDIN> =~ /^[Yy][Ee]{0,1}[Ss]{0,1}\n$/;
printf OUTPUT "#-------------Parameter Settings--------------\nupperbound=%s\nmultimodel=%s\npthreshold=%f\nearlytermination=%s\n",
$upperbound?'y':'n', $multimodel?'y':'n', $pthreshold, $earlytermination?'y':'n'; 

if(defined $stdout) {
	open MYSTDOUT, (">"."$outputdir/".$stdout) or die "can't write to file.\n";
	select MYSTDOUT;
}

#----------getting frequency file--------------
print OUTPUT "\n#----------getting frequency file--------------\n";
my @freqarr;
my @countarr;
my $wt0 = gettimeofday();
my $ct0=clock();
my $cont=0;
while(<INPUT>){
	chomp;
	next if (not $_ =~ /^[\d]+(\.[\d]+){0,1}\s+[\d]+/);
	my($freq, $count) = split(/\s+/, $_, 2);
	push @freqarr, $freq;
	push @countarr, $count;
	$cont = ($cont || (int($freq) != $freq));
}

my $wt1 = gettimeofday();
my $ct1=clock();
print OUTPUT "Frequencies:\t", scalar @freqarr, "\n";
print OUTPUT "Used Time:\tWall=", $wt1-$wt0, ", CPU=", $ct1-$ct0, "\n";



#------------fitting---------------------------
print OUTPUT "\n#----------Constructng Models--------------\n";
$wt0 = gettimeofday();
$ct0 = clock();
my $modArr=[];
my $B=$#freqarr;

while(1){
	$upperbound = 1 if ($B != $#freqarr);
	my @para=build_mod(\@freqarr, \@countarr, $B);
	printf "\nrange= %u-%u, b/a= %.2f, sum=%u\nxmin= %.5f, xmax= %f, a= %.3f, KS=%.5f\n", 
		$para[4],$para[5],$para[1]/$para[0],$para[6],$para[0],$para[1]+0,$para[2],$para[3];
	printf OUTPUT "\nrange= %u-%u(%u), b/a= %.2f, sum=%.5f\nxmin= %.5f, xmax= %f, a= %.3f, KS=%.5f\n", 
		$para[4],$para[5],$para[5]-$para[4]+1,$para[1]/$para[0],$para[6],$para[0],$para[1]+0,$para[2],$para[3];
	push @{$modArr}, \@para;
	$B = $para[4];
	last if $B<=1;
	last if (!$multimodel && scalar @{$modArr} != 0);
}
$wt1 = gettimeofday();
$ct1 = clock();
print OUTPUT "\n", scalar @{$modArr}, " models constructed.\n";
print OUTPUT "Used Time:\tWall=", $wt1-$wt0, ", CPU=", $ct1-$ct0, "\n";

#my @checklist = sort {($b->[1]/$b->[0]) <=> ($a->[1]/$a->[0])} @{$modArr};

#--------------plotting------------------------------


my $N=0;
for my $i(0..$#freqarr){
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
	}
}


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
my (@Ddatasets);
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
}

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


#---------goodness-of-fit----------------------
for (reverse @{$modArr}){
	print "===============Starting goodness-of-fit======================\n";
	print OUTPUT "\n-------------Starting goodness-of-fit---------------\n";
	$wt0 = gettimeofday();
	$ct0 = clock();
	my $model=$_;
	my $n = 100;
	my $bt=0;#better than
	printf "range= %u-%u, b/a= %.2f, sum=%u\nxmin= %.5f, xmax= %f, a= %.3f, KS=%.5f\n", 
		$model->[4],$model->[5],$model->[1]/$model->[0],$model->[6],$model->[0],$model->[1]+0,$model->[2],$model->[3];
	printf OUTPUT "\nrange= %u-%u(%u), b/a= %.2f, sum=%u\nxmin= %.5f, xmax= %f, a= %.3f, KS=%.5f\n", 
		$model->[4],$model->[5],$model->[5]-$model->[4]+1,$model->[1]/$model->[0],$model->[6],$model->[0],$model->[1]+0,$model->[2],$model->[3];
	#print "Xmin=",$model->[0], "\tXmax=",$model->[1],"\ta=",$model->[2],"\tKS=",$model->[3],"\n";
	print "Number of set of sythetic data generated = $n\n";
	#print OUTPUT "Xmin=",$model->[0], "\tXmax=",$model->[1],"\ta=",$model->[2],"\tKS=",$model->[3],"\n";
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
		
		next if (!$printSynthetic);
		my @ccdf_s;
		my $sum_s = 0;
		my ($sum_base, $sum_top) = (0,0);
		for my $j(reverse 0..$#freqarr_s){
			$sum_s+=$countarr_s[$j];
			push @ccdf_s, $sum_s;
			$sum_base+=$countarr_s[$j] if ($freqarr_s[$j] > $model_s[1]);
			$sum_top+=$countarr_s[$j] if ($freqarr_s[$j] >= $model_s[0]);
		}
		$sum_base/=$sum_s;
		$sum_top/=$sum_s;
		@ccdf_s = reverse @ccdf_s;
		my (@ccdf_smodel, @freqs_smodel);
		for my $j(0..$#freqarr_s){
			next if($freqarr_s[$j] == $model_s[1] && $sum_base ==0 && $cont);
			if($freqarr_s[$j]>=$model_s[0] && $freqarr_s[$j] <= $model_s[1]){
				push @freqs_smodel, $freqarr_s[$j];
				push @ccdf_smodel, $sum_base+($sum_top-$sum_base)*zeta($freqarr_s[$j], $model_s[1], $model_s[2])/zeta($model_s[0], $model_s[1], $model_s[2]) if(!$cont);
				push @ccdf_smodel, $sum_base+($sum_top-$sum_base)*zeta_cont($freqarr_s[$j], $model_s[1], $model_s[2])/zeta_cont($model_s[0], $model_s[1], $model_s[2]) if($cont);
			}
		}
		for my $j(0..$#freqarr_s){
			next if ($freqarr_s[$j]==0 || $ccdf_s[$j] ==0);
			$freqarr_s[$j] = log($freqarr_s[$j])/log(10);
			$ccdf_s[$j] = log($ccdf_s[$j]/$sum_s)/log(10);
		}
		for my $j(0..$#freqs_smodel){
			next if ($freqs_smodel[$j] ==0 || $ccdf_smodel[$j]==0);
			$freqs_smodel[$j] = log($freqs_smodel[$j])/log(10);
			$ccdf_smodel[$j] = log($ccdf_smodel[$j])/log(10);
		}
		my $sdataCCDF = Chart::Gnuplot::DataSet->new(
			xdata =>	\@freqarr_s,
			ydata =>	\@ccdf_s);
		my $smodelCCDF = Chart::Gnuplot::DataSet->new(
			xdata =>	\@freqs_smodel,
			ydata =>	\@ccdf_smodel,
			style	=> "lines",
			color	=> "#007f00",
			linetype=> "solid",
			width	=> 3);
		my $syntheticCCDFChart = Chart::Gnuplot->new(
			output	=> "$ofname"."_s$i"."_ccdf.png",
			title	=> "CCDF",
			xlabel	=> "frequency",
			ylabel	=> "CCDF",
			bg		=> "white",
			grid	=>{linetype=>'dash'},
		);
		if(scalar @freqs_smodel != 0){$syntheticCCDFChart->plot2d($sdataCCDF, $smodelCCDF);}
		else{$syntheticCCDFChart->plot2d($sdataCCDF);}
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
	print OUTPUT "alpha avg=$A_avg, std=$A_std\tKS avg=$D_avg, std=$D_std\n\n";

	
	
	
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
	my $targetD = shift;
	my $targetA = shift;
	#$u_limit=$#freqarr if (not defined $u_limit);
	my $sum =0;
	my $lnsum=0;
	my $alpha;
	my $xmin='inf';
	my $xmax=$freqarr->[$u_limit]; $xmax = 'inf' if (!$upperbound );
	my $Dmin='inf';
	my $summin=0;
	my $l_limit=$u_limit;
	#my @rev_CDF;
	my @_Aarr;	
	my @_sumarr;
	my ($xdata,$ydata) = ([],[]); #freqs, D-values, p-values
	for my $n(0..$u_limit){
		my $_xmin=$freqarr->[$u_limit-$n];
		#print "_xmin=$_xmin\t";
		$sum+=$countarr->[$u_limit-$n];
		$lnsum+=$countarr->[$u_limit-$n]*log $_xmin;
		next if ($n == 0);
		my $_alpha;
		if(!$cont){$_alpha=find_alpha($_xmin, $xmax, $sum, $lnsum);}
		else{$_alpha = 1+$sum/($lnsum-$sum*log($_xmin)) ;};  #TODO don't know the MLE of alpha with upper bound
		if(!$cont && zeta($_xmin, $xmax, $_alpha)==0){ next;} #The estimated alpha is too big and _xmin and xmax are too close
		if($cont && $_xmin**(1-$_alpha)-$xmax**(1-$_alpha) ==0) {next;}

		#TODO: is there a way to avoid this inner for-loop?
		my $D;
		if(!$cont){$D=KS($freqarr, $countarr, $u_limit-$n,$u_limit,$_alpha,$sum,!$upperbound);}
		else{$D=KS_cont($freqarr, $countarr, $u_limit-$n,$u_limit,$_alpha,$sum,!$upperbound);};
		next if($D==-1);
		push @$xdata, $_xmin;
		push @$ydata, $D;

		#print "D=$D\tDmin=$Dmin\n";
		if($D<$Dmin){
			$Dmin=$D;
			$xmin=$_xmin;
			$alpha=$_alpha;
			$l_limit=$u_limit-$n;
			$summin=$sum;
		}
	}
	if($xmin != 'inf'){
		push @$xdatas, $xdata;
		push @$ydatas, $ydata;
	}
	else{$xmin = -1};
	return my @para=($xmin, $xmax, $alpha, $Dmin, $l_limit,$u_limit, $summin);
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
	for my $i(0..($ct_less+$ct_more) ){
		my $r=rand;
		if($r<=$ct_less/($ct_less+$ct_more))	{$nless++;}
		else									{$nmore++;}
	}
	my $rddata_1= sample(\@headfreq, \@headcount, $nless);
	my $rddata_2;
	if(!$cont){$rddata_2= plrand($Xmin, $Xmax, $alpha, $nmore);}
	else{$rddata_2= plrand_cont($Xmin, $Xmax, $alpha, $nmore);}
	#--the two hashes should not be merged, so that their keys can be sorted faster
	my ($sfreqs1, $scounts1, $sfreqs2, $scounts2) = ([],[],[],[]);
	data_count($rddata_1,$sfreqs1,$scounts1);
	data_count($rddata_2,$sfreqs2,$scounts2);
	@{$freqarr_s} = (@{$sfreqs1}, @{$sfreqs2});
	@{$countarr_s} = (@{$scounts1}, @{$scounts2});
}
