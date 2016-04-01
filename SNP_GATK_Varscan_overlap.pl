#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 <sample> <varscan> <gatk> <statfile>\n" if(@ARGV != 4);

my ($sample,$varscan,$gatk,$file)=@ARGV;

my %nsites;
my ($varsnp,$gatksnp,$gatkuniq,$varuniq,$overlap)=(0,0,0,0,0);
open IN1, ($varscan =~ /\.gz/? "gzip -dc $varscan |" : $varscan) || die $!;
<IN1>;
while(<IN1>){
	chomp;
	my @F = split /\t/;
	my ($chr,$pos,$ref,$var)=@F[0,1,2,3];
	my $id="$chr\t$pos\t$ref\t$var";
	my $type = $F[12];
	if($type eq "Germline"){
		$nsites{$id}=1;
		$varsnp++;
	}
}
close IN1;

open IN2, ($gatk =~ /\.gz/? "gzip -dc $gatk |" : $gatk) || die $!;
chomp(my $head=<IN2>);
print "$head\n";
while(<IN2>){
	chomp;
	my @F = split /\t/;
	my ($chr,$pos,$ref,$var)=@F[0,1,3,4];
	my $id="$chr\t$pos\t$ref\t$var";
	if(exists $nsites{$id}){
		print "$_\n";
		$overlap++;
	}
	else{
		$gatkuniq++;
	}
}
close IN2;

$gatksnp = $gatkuniq + $overlap;
$varuniq = $varsnp  - $overlap;
open OUT,">>$file" || die $!;
print OUT "Samples\tGATK\tVarscan\tGATK_uniq\tVarscan_uniq\tOverlap\n";
print OUT "$sample\t$gatksnp\t$varsnp\t$gatkuniq\t$varuniq\t$overlap\n";


