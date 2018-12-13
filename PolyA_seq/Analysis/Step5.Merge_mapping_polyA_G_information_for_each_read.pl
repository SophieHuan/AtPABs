#!/usr/bin/perl
use strict;
use warnings;

# The input file
my $label="Col_1";
my $file1="../TAIR10/Sequence/TAIR10-cDNA_length_from_ensembel.txt";
my $file2="Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq_polyA_length.txt_Ginfor.txt";
my $file3="Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq_delA_cDNA_mapped.sam";

# Example of file1
#Gene|stable     ID
#AT3G11415|AT3G11415.1   1002
#AT5G24735|AT5G24735.1   585

# Example of file2
#index   polyA_length    seq     A_per   T_per   C_per   G_per   G_count NotA_count
#ST-E00243:521:HNFCLCCXY:3:1101:15514:1784       40      AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAACAAA        0.95    0       0.05    0       0       2
#ST-E00243:521:HNFCLCCXY:3:1101:10145:1836       25      AAAAAAAAAAAAAAAAAAAAAAAAA       1       0       0       0       0       0

# Example of file3
#ST-E00243:521:HNFCLCCXY:3:1101:15514:1784       0       AT1G20620|AT1G20620.4   2712    1       67M     *       0       0       ACTGTCTTTGTAATCTCTTTTTACAATAAATCAATGTTTCTTGAAGCAGTAAACAATTTCTAAGAAC     FFJJ<JJF7FAFJJJJFJ<JFAJJJFJJJJJJF-JJJJFJJJJFJJJFJFJJJJJ-<FJJJ-7A<FF     AS:i:-10        XS:i:-10        XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:65T0T0     YT:Z:UU
#ST-E00243:521:HNFCLCCXY:3:1101:10145:1836       0       AT2G01190|AT2G01190.1   2535    42      49M     *       0       0       CCATTTTCCTGGATCTTATATCATATTATTATCATTCATTTCATTTTGT       FFJFJ7JJFJJJFFFJJJJJJJJ-FJFJJAJJJJJJJJJJJJAFJJ<JJ AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:49 YT:Z:UU


# This inputs the length information of each cDNA
my %cDNA_end;
open (IN1,$file1)||die;
while (my $line=<IN1>) {
	chomp $line;
	my @lines=split("\t",$line);
	$cDNA_end{$lines[0]}=$lines[1];
}
close IN1;

# This inputs the poly(A) tail and G content information of each read
my %polyA;my %G;
open (IN2,$file2)||die;
while (my $line=<IN2>) {
	chomp $line;
	my @lines=split("\t",$line);
	$polyA{$lines[0]}=$lines[1];
	$G{$lines[0]}="$lines[3]\t$lines[4]\t$lines[5]\t$lines[6]\t$lines[7]\t$lines[8]";
}
close IN2;

open (IN3,$file3)||die;
open (OUT,">".$label."_mapping_polyA_G_infor.txt")||die;
print OUT "gene\tindex\tstart\tmapping\tend\tcDNA_length\tdistance_to_cDNA_end\tcDNA_percentage\tpolyA_length\tA_percentage\tT_percentage\tC_percentage\tG_percentage\tG_count\tNotA_count\n"; 
while (my $line=<IN3>) {
	chomp $line;
	my @lines=split("\t",$line);

	# This calculates the 3'-end position of each read
	if($lines[5]=~/(\d+)M$/ and exists $cDNA_end{$lines[2]}){
		my @map=split(/[A-Z]/,$lines[5]);
		my @num=split(/\d+/,$lines[5]);
		my $l=0;
		for (my $i=1;$i<@num;$i++){
			if ($num[$i]=~/^M/){
				$l=$l+$map[($i-1)];
			}
			elsif ($num[$i]=~/^I/){
				$l=$l;
			}
			elsif ($num[$i]=~/^D/){
				$l=$l+$map[($i-1)];
			}
		}
		my $end=$lines[3]+$l-1;
		my $per=$end/$cDNA_end{$lines[2]};
		my $dis=$cDNA_end{$lines[2]}-$end;
		if (exists $polyA{$lines[0]}) {
			print OUT "$lines[2]\t$lines[0]\t$lines[3]\t$lines[5]\t$end\t$cDNA_end{$lines[2]}\t$dis\t$per\t$polyA{$lines[0]}\t$G{$lines[0]}\n";
		}
	}
}
close IN3;
close OUT;
