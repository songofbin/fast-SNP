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
my (%genes,%sites,%annos,%DBs);
my $ids=0;
read_list($input);
read_gene($genelist);

print OUT "Chr\tPos\tRef\tAlt\t";
print OUT "Func.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tgenomicSuperDups\tavsnp142\tclinvar_20150629\t";
print OUT "Type\tDB_num\tDB_name\tSample_num\tSample_ratio\tSNP_fre_mean\t";
print OUT join("\t",@samples),"\n";

foreach my $site (sort {$a cmp $b} keys %sites){
	my ($sample_count,$snp_fre)=(0,0);
	my @infos;
	for my $i (1..$ids){
		my $fre_info = ".";
		if(exists $sites{$site}{$i}){
			$sample_count++;
			$fre_info = $sites{$site}{$i};
			my ($normal_fre,$tumor_fre)=(split /,/,$fre_info)[0,1];
			$snp_fre+=($normal_fre+$tumor_fre)/2;
		}
		@infos=(@infos,$fre_info);
	}

	$snp_fre = $snp_fre/$sample_count;
	$snp_fre = sprintf("%.2f",$snp_fre);
	my @anno_infos = @{$annos{$site}};
	my $gene = $genes{$site};	
	my $DB_info = (exists $DBs{$gene})?$DBs{$gene}:"NA\t0\tNA";
	my $sample_fre = $sample_count/$ids;
	$sample_fre = sprintf("%.2f",$sample_fre);
	print OUT "$site\t",join("\t",@anno_infos),"\t";
	print OUT "$DB_info\t$sample_count\t$sample_fre\t$snp_fre\t",join("\t",@infos),"\n";
}
close OUT;

sub read_list{
	my ($input) = @_;
	open LIST,"$input"  || die $!;
	open OUT,">$output" || die $!;
	while(<LIST>){
		chomp;
		my ($name,$file)=(split)[0,1];
		$files{$name}=$file;
		push @samples,$name;
	}
	close LIST;

	for my $sample(@samples){
		$ids++;
		&read_file($files{$sample},$ids,$type);
	}
}

sub read_file{
    my ($file,$id,$type)=@_;
    open IN, ($file =~ /\.gz/? "gzip -dc $file |" : $file) || die $!;
    <IN>;
    while(<IN>){
		chomp;
		my @F = split /\t/;
		my ($chr,$pos,$ref,$var,$normal_fre,$tumor_fre,$gene,$func) = @F[0,1,3,4,8,12,14,16];

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

		my $snp = "$chr\t$pos\t$ref\t$var";
		my @anno_infos=@F[13..21];
		$normal_fre = sprintf("%.2f",$normal_fre);
		$tumor_fre = sprintf("%.2f",$tumor_fre);	
		$normal_fre=1 if($normal_fre==1);
		$tumor_fre=1  if($tumor_fre==1);
		my $fre_info="$normal_fre,$tumor_fre";
		$sites{$snp}{$id}=$fre_info;
		if(!exists $genes{$snp}){$genes{$snp}=$gene;}
		if(!exists $annos{$snp}){
			@{$annos{$snp}} = @anno_infos; 
		}
	}
	close IN;
}

sub read_gene{
	my ($file) = @_;
	open FILE,$file || die $!;
	while(<FILE>){
		chomp;
		next if(/^Gene/);
		my ($gene,$type,$num,$dbs)=(split /\t/)[0,1,2,3];
		$DBs{$gene}="$type\t$num\t$dbs";
	}
	close FILE;
}
