#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $usage=<<USAGE;
Description
        Statitics of overlap genes.
Parameter
        --in        [STR] input anno file list, required
        --out       [STR] output file, required
        --type      [STR] gene type, all or nonsynonymous gene, default all
        
Exmple          
        perl $0 --in filelist --out out.txt --type nonsynonymous
                
USAGE

my ($input,$output,$type);
my $genelist = "/ifs5/ST_ANNO/USER/songbin/2.project/breast_LR/gene/breast_gene.txt";

GetOptions(
	"in:s"    => \$input,
	"out:s"   => \$output,
	"type:s"  => \$type,
);
die "$usage" if(!$input || !$output);
$type ||= "all";

my (%files,@samples);
open LIST,"$input" || die $!;
open OUT,">$output" || die $!;
while(<LIST>){
	chomp;
	my ($name,$file)=(split)[0,1];
	$files{$name}=$file;
	push @samples,$name;
}
close LIST;

my (%fre,%site,%annos);
my $ids=0;
for my $sample(@samples){
	$ids++;
	&read_file($files{$sample},$ids,$type);
}
read_gene($genelist);

print OUT "Gene\tType\tDB_num\tDB_name\tSample_num\tSample_ratio\tSNP_num";
for my $sample(@samples){
	print OUT "\t$sample";
}
print OUT "\n";

foreach my $key (sort {$a cmp $b} keys %site){
	my ($sample_count,$site_count)=(0,0);
	my @infos;
	for my $i (1..$ids){
	        if(exists $site{$key}{$i}){
			$sample_count++;
			my $count1 = $#{$site{$key}{$i}} + 1;
			$site_count+=$count1;
			@infos=(@infos,$count1);
		}
		else{
			@infos=(@infos,0);
		}
	}
	my $sample_fre=$sample_count/$ids;
	$sample_fre = sprintf("%.2f",$sample_fre);
	print OUT "$key\t";
	if(exists $annos{$key}){
		print OUT "$annos{$key}";
	}
	else{
		print OUT "NA\t0\tNA";
	}
	print OUT "\t$sample_count\t$sample_fre\t$site_count\t",join("\t",@infos),"\n";
}
close OUT;

sub read_file{
    my ($file,$id,$type)=@_;
    open IN, ($file =~ /\.gz/? "gzip -dc $file |" : $file) || die $!;
    <IN>;
    while(<IN>){
	chomp;
	my @F = split /\t/;
	
	my ($chr,$pos,$ref,$var,$normal_fre,$tumor_fre,$gene,$func) = @F[0,1,3,4,8,13,14,16];

		my $flag=0;
		if($type eq "nonsynonymous"){
			if($func=~/nonsynonymous/ or $func=~/stop/){
				$flag=1;
			}
		}
		else{
			$flag=1;
		}
		next if($flag==0);

	my $mutation = $chr."-".$pos.":".$ref.">".$var;

	next if($gene=~/NONE/);	
	if($gene =~ /;/){
		$gene=(split /;/,$gene)[0];
		push @{$site{$gene}{$id}},$mutation;
	}	
	elsif($gene =~ /,/){
		my @genes = split /,/,$gene;
		foreach my $name (@genes){
			push @{$site{$name}{$id}},$mutation;
		}	
	}
	else{
		push @{$site{$gene}{$id}},$mutation;
	}
    }
}

sub read_gene{
	my ($file) = @_;
	open FILE,$file || die $!;
	while(<FILE>){
		chomp;
		next if(/^Gene/);
		my ($gene,$type,$num,$dbs)=(split /\t/)[0,1,2,3];
		$annos{$gene}="$type\t$num\t$dbs";
	}
}

sub mean{
	my ($fre_ref) = @_;
	my @fre = @{$fre_ref};
	my $total=0;
	for my $i (0..$#fre){
		$total+=$fre[$i];
	}
	$total/($#fre+1);
}
