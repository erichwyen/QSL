#!usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw(gettimeofday clock);

sub getline;
sub word2freq;

#getting flag values 
my ($ifname,$ofname);
#my $i=0;
#while($i < scalar @ARGV){
#	if($ARGV[$i] eq '-o'){
#		$ofname = $ARGV[$i+1];
#		if(not defined $ofname) {die "Usage: -o <output name>\n";}
#		splice (@ARGV, $i, 1);
#	}
#	else{$i++;}
#}
#undef $i;
$ifname = shift @ARGV or die "Please enter the input file name\n" ;
#if(not defined $ofname){$ofname = $ifname; $ofname =~s/\..+$//}#default
open INPUT, $ifname or die "Can't open the file:$!\n";
$ofname = shift;

my $wordlength = shift;
#print "Please enter the length of word with which to sample the text:";
#my $wordlength = <STDIN>;#in Chinese character
#chomp $wordlength;
#die "Word length must be an positive integer." if (int $wordlength != $wordlength or int $wordlength <= 0);
my $word2freq = {};

my $wtime0 = gettimeofday();
my $ctime0=clock();
my @wordlist = getline($word2freq);

#my @debug;
#TODO:Avoid overflow when dealing with large amount of data.
sub getline{
	my $word2freq = shift;
	my @wordque;
	for(1..$wordlength){push @wordque, [];}
	my $startflag = $wordlength-1;
	my $currenthead = 0;
	my $cit = 0; #character iterator
	my $intag =0; #within a tag
	#my $life = 700;
	while(<INPUT>){
		#removing tags and white spaces
		chomp;
		
		s|<qe>[^(<qe>)(</qe>)]*</qe>||g;
		s|<qe>.*||g;
		s|.*</qe>||g;
		s|<qn>[^(<qn>)(</qn>)]*</qn>||g;
		s|<qn>.*||g;
		s|.*</qn>||g;
		s|<[^<>]*>||g;
		s|。||g;
		s|，||g;
		s|「||g;
		s|」||g;
		s|、||g;
		s|△||g;
		s|；||g;
		s|：||g;
		s|…||g;
		s|？||g;
		s|—||g;
		s|─||g;
		s|囗||g;
		s|『||g;
		s|』||g;
		s|（||g;
		s|）||g;
		s|○||g;
		s|△||g;
		#constructing the word list one character at a time.
		my @arr = split("", $_); 
		my $composing = 0;
		#To prevent overflow, store word with a circular queue.
		for my $c(@arr){
			if($cit%3==0 && (ord($c)<224 || ord($c)>239)) {next;}#Chinese char-triplets always start within this range
			#print ord($c), " ";
			push @{$wordque[$currenthead]}, $c;
			$cit++;
			if($cit%3==0){
				if(not $composing){
					if(	ord($wordque[$currenthead]->[0]) == 239
					&& 	ord($wordque[$currenthead]->[1]) == 154
					&&	ord($wordque[$currenthead]->[2]) == 174){ $composing = 1;}
				}
				if($composing){
					if(	ord($wordque[$currenthead]->[-3]) == 239
					&& 	ord($wordque[$currenthead]->[-2]) == 154
					&&	ord($wordque[$currenthead]->[-1]) == 175){$composing = 0;}
				}
				if(not $composing){
					if($startflag<=0){#or equivalently, @wordque is full. ready to output the first word
						my $i = ($currenthead+1)%(scalar @wordque);
						my $newword="";
						for (@wordque){
							for (@{$wordque[$i]}){
								$newword.=$_;
							}
							$i++; $i%=(scalar @wordque);
						}
						if (exists $word2freq->{$newword}){$word2freq->{$newword}++;}
						else{$word2freq->{$newword}=1;}	
						#print "\t$newword\n";
					}
					$startflag--;			
					$currenthead++;	$currenthead%=(scalar @wordque);
					@{$wordque[$currenthead]} = ();
				}
			}
		}
		#$life--;
		#last if $life ==0;};
	}
}

#$mu->record('got data');

my @sortedlist = sort {$word2freq->{$b} <=> $word2freq->{$a};} keys(%{$word2freq});

#$mu->record('sorted');
#$mu->dump();

my %freq2count;
for (keys %{$word2freq}){
	if(exists $freq2count{$word2freq->{$_}} ){$freq2count{$word2freq->{$_}}+=1;}
	else {$freq2count{$word2freq->{$_}} =1; }
}
my @freqarr = sort {$a<=>$b} keys(%freq2count);
my @countarr = ();

my $wordcount = 0;
for (@freqarr){
	push @countarr, $freq2count{$_};
	$wordcount+=$_*$freq2count{$_};
}

#computing pdf, cdf, ccdf;
my @pdf;
my @ccdf;
my $sum = 0;
for my $i (1..scalar @countarr){
	push @pdf, $countarr[-$i]/(scalar @sortedlist);
	$sum+=$countarr[-$i];
	#print $sum, "\n";
	die "$sum goes over 1!!\n" if $sum >scalar @sortedlist;
	push @ccdf, $sum/(scalar @sortedlist);
	#print $i, ": " ,1-$cdf[$#cdf],"\n";
}

@pdf = reverse @pdf;
@ccdf = reverse @ccdf;

#---------------------output---------------------------------
open WRDFREQ, (">".$ofname.".wrdfreq") or die "Can't write to file!!:$!\n";
for(@sortedlist){
	print WRDFREQ $_, ": ", $word2freq->{$_}, "\n";
}

open FREQCOUNT, (">".$ofname.".freqcnt") or die "Can't write to file!!:$!\n";
print FREQCOUNT "Frequency\tCount\n";
for my $i(0..$#freqarr){
	print FREQCOUNT "$freqarr[$i]\t$countarr[$i]\n";
}

print "The number of total words: ", $wordcount, "\n";
print "The number of different words: ", scalar @sortedlist,"\n";
print "The number of different frequencies: ", scalar @freqarr,"\n"; 

my $wtime1 = gettimeofday();
my $ctime1 = clock();
print "Wall time = ", $wtime1-$wtime0, "\n";
print "CPU time = ", $ctime1-$ctime0, "\n";


