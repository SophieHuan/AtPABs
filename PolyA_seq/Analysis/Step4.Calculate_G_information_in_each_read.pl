#!/usr/bin/perl
use strict;
use warnings;

# The input file
my $file="Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq_polyA_length.txt";

# Examples of the input file
# index	poly(A)_tail	poly(A)_sequences
#ST-E00243:521:HNFCLCCXY:3:1101:15514:1784	40	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAACAAA
#ST-E00243:521:HNFCLCCXY:3:1101:15432:1819	33	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#ST-E00243:521:HNFCLCCXY:3:1101:12824:1836	45	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#ST-E00243:521:HNFCLCCXY:3:1101:11941:2733	30	AAAGAAAAAAAAAAAAAAAAAAAAAAAAAA
#ST-E00243:521:HNFCLCCXY:3:1101:14915:2751	29	AAAAAAAGAAAAAAAAAAAAAAAAAAAAA
#ST-E00243:521:HNFCLCCXY:3:1101:8704:1907	33	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

open (IN,$file);
open (OUT,">".$file."_Ginfor.txt")||die;
print OUT "index\tpolyA_length\tseq\tA_per\tT_per\tC_per\tG_per\tG_count\tNotA_count\n";
while (my $line=<IN>) {
	chomp $line;
	my @lines=split("\t",$line);

	# This calculates the count of A/C/G/T in each poly(A) tail sequence
	my $A=0;my $C=0;my $G=0;my $T=0;my $notA=0;
	$A=$lines[2]=~ tr/A/A/;
	$C=$lines[2]=~ tr/C/C/;
	$G=$lines[2]=~ tr/G/G/;
	$T=$lines[2]=~ tr/T/T/;
	$notA=$C+$G+$T;
	
	# This calculate the percentage of A/C/G/T in each poly(A) tail sequence
	if (length($lines[2])>0){
		my $A_per=$A/length($lines[2]);
		my $C_per=$C/length($lines[2]);
		my $G_per=$G/length($lines[2]);
		my $T_per=$T/length($lines[2]);
		print OUT "$line\t$A_per\t$T_per\t$C_per\t$G_per\t$G\t$notA\n";
	}
}
close IN;
close OUT;
