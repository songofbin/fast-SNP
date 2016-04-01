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

my (%fre,%site);
my $ids=0;
for my $sample(@samples){
	$ids++;
	&read_file($files{$sample},$ids,$type);
}

print OUT "Gene\tSample_num";
for my $sample(@samples){
	print OUT "\t$sample\_num\t$sample\_fre\t$sample\_sites";
}
print OUT "\n";

foreach my $key (sort {$a cmp $b} keys %fre){
	my $count=0;
	my @infos;
	for my $i (1..$ids){
	        if(exists $fre{$key}{$i}){
			$count++;
			my $count1    = $#{$fre{$key}{$i}} + 1;
			my $fre_mean1 = &mean(\@{$fre{$key}{$i}});
			@infos=(@infos,$count1,$fre_mean1,join("; ",@{$site{$key}{$i}}));
		}
		else{
			@infos=(@infos,0,"-","-");
		}
	}		
	print OUT "$key\t$count\t",join("\t",@infos),"\n";
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
		push @{$fre{$gene}{$id}},$tumor_fre;
		push @{$site{$gene}{$id}},$mutation;
	}	
	elsif($gene =~ /,/){
		my @genes = split /,/,$gene;
		foreach my $name (@genes){
			push @{$fre{$name}{$id}},$tumor_fre;
			push @{$site{$name}{$id}},$mutation;
		}	
	}
	else{
		push @{$fre{$gene}{$id}},$tumor_fre;
		push @{$site{$gene}{$id}},$mutation;
	}
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
