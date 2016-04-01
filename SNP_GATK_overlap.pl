#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 <normal> <tumor>\n" if(@ARGV != 2);

my ($normal,$tumor)=@ARGV;

my %nsites;
open IN1, ($normal =~ /\.gz/? "gzip -dc $normal |" : $normal) || die $!;
<IN1>;
while(<IN1>){
	chomp;
	my @F = split /\t/;
	my ($chr,$pos,$ref,$var)=@F[0,1,3,4];
	my $info = $F[71];
	my $id="$chr\t$pos\t$pos\t$ref\t$var";
	$nsites{$id}=$info;
}
close IN1;

#print "Chr\tStart\tEnd\tRef\tAlt\tGT_n\tRef_n\tAlt_n\tFre_n\tGT_t\tRef_t\tAlt_t\tFre_t\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tgenomicSuperDups\tavsnp142\tclinvar_20150629\n";

open IN2, ($tumor =~ /\.gz/? "gzip -dc $tumor |" : $tumor) || die $!;
chomp(my $head=<IN2>);
my @heads=(split /\t/,$head)[5..59];
print "Chr\tStart\tEnd\tRef\tAlt\tGT_n\tRef_n\tAlt_n\tFre_n\tGT_t\tRef_t\tAlt_t\tFre_t\t";
print join("\t",@heads),"\n";
while(<IN2>){
	chomp;
	my @F = split /\t/;
	my ($chr,$pos,$ref,$var)=@F[0,1,3,4];
	my $info = $F[71];
	my $id="$chr\t$pos\t$pos\t$ref\t$var";
	if(exists $nsites{$id}){
		my @ninfos=split /:/,$nsites{$id};
		my $ngt = $ninfos[0];
		my ($nref,$nvar)=(split /,/,$ninfos[1])[0,1];
		my $nfre = $nvar/($nref+$nvar);
		$nfre = sprintf("%.2f",$nfre);
		$nfre = 1 if($nfre==1);	
		my @tinfos=split /:/,$info;
		my $tgt = $tinfos[0];
		my ($tref,$tvar)=(split /,/,$tinfos[1])[0,1];
		my $tfre = $tvar/($tref+$tvar);
		$tfre = sprintf("%.2f",$tfre);	
		$tfre = 1 if($tfre==1);
		my @annos=@F[5..59];
		print "$id\t$ngt\t$nref\t$nvar\t$nfre\t$tgt\t$tref\t$tvar\t$tfre\t";
		print join("\t",@annos),"\n";
	}
}
close IN2;
