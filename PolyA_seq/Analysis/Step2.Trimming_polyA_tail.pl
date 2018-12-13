#!/usr/bin/perl
use strict;
use warnings;

# The input file
my $file="Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq";

# Example of the input file
#@ST-E00243:521:HNFCLCCXY:3:1101:21227:2540 1:N:0:CGATGT
#TGTCTCTTACTCTATTATTTCATTCGAATACATCTTCTATAAATTATATATTCCACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#+
#FFJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJ-FJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJFJJJJ<FJJJJJJJJJJJJ-FAJJJJJJFJJFJJJJJA-
#@ST-E00243:521:HNFCLCCXY:3:1101:14580:2522 1:N:0:CGATGT
#TTTTTGATAGATTATATCAAATCCATGGATACTTTCTATAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#+
#FFJJJJJJJJJJJJJJJJJ-FJJJFJJJJJ<FJJJJJJJJFJJJJJJFAFJJJFJJJJJJJJJA-F<AAJJFAJJJJJJJJJF<7F<-F-<<<-7FJ<JFJ<-FA<7


# This defines the seperator of each read in .fastq document
$/="\@ST-";

open (IN,$file)||die;
open (OUT1,">".$file."_delA.fq")||die;
open (OUT2,">".$file."_polyA_length.txt")||die;
while (my $line=<IN>) {
	chomp $line;
	$line=~s/\n$// if ($line=~/\n$/);
	next if ($line eq "");
	my @lines=split(/\n/,$line);

	# This defines the poly(A) tail as constitutive A¡¯s at the 3¡ä-end of a read, allowing up to 5 interspersed G bases 
	if ($lines[1]=~/(A+)(A|G)(A+)(A|T|C|G)(A+)(A|G)(A+)(A|T|C|G)(A+)(A|G)$/) {

		# This calculates the length of poly(A) tail
		my $len=length($1)+length($2)+length($3)+length($4)+length($5)+length($6)+length($7)+length($8)+length($9)+length($10);

		# Trimmed reads with shorter than 18 nt in length are removed
		my @lables=split(" ",$lines[0]);
		if (length($lines[1])>= ($len+18)) {

			#print reads that trimmed ploy(A) tail
			print OUT1 "\@ST-$lines[0]\n".substr($lines[1],0,(length($lines[1])-$len))."\n$lines[2]\n".substr($lines[3],0,(length($lines[1])-$len))."\n";

			#print reads information with poly(A) tail length
			print OUT2 "ST-$lables[0]\t".$len."\t".$1.$2.$3.$4.$5.$6.$7.$8.$9.$10."\n";
		}
	}
	elsif ($lines[1]=~/(A+)(A|T|C|G)(A+)(A|G)(A+)(A|T|C|G)(A+)(A|G)(A+)(A|G)$/) {
		my $len=length($1)+length($2)+length($3)+length($4)+length($5)+length($6)+length($7)+length($8)+length($9)+length($10);
		my @lables=split(" ",$lines[0]);
		if (length($lines[1])>= ($len+18)) {
			print OUT1 "\@ST-$lines[0]\n".substr($lines[1],0,(length($lines[1])-$len))."\n$lines[2]\n".substr($lines[3],0,(length($lines[1])-$len))."\n";
			print OUT2 "ST-$lables[0]\t".$len."\t".$1.$2.$3.$4.$5.$6.$7.$8.$9.$10."\n";
		}
	}
	elsif ($lines[1]=~/(A+)(A|G)(A+)(A|T|C|G)(A+)(A|T|C|G)(A+)(A|G)$/) {
			my $len=length($1)+length($2)+length($3)+length($4)+length($5)+length($6)+length($7)+length($8);
			my @lables=split(" ",$lines[0]);
			if (length($lines[1])>= ($len+18)) {
				print OUT1 "\@ST-$lines[0]\n".substr($lines[1],0,(length($lines[1])-$len))."\n$lines[2]\n".substr($lines[3],0,(length($lines[1])-$len))."\n";
				print OUT2 "ST-$lables[0]\t".$len."\t".$1.$2.$3.$4.$5.$6.$7.$8."\n";
		}
	}
	elsif ($lines[1]=~/(A+)(A|T|C|G)(A+)(A|G)(A+)(A|T|C|G)(A+)(A|G)$/) {
			my $len=length($1)+length($2)+length($3)+length($4)+length($5)+length($6)+length($7)+length($8);
			my @lables=split(" ",$lines[0]);
			if (length($lines[1])>= ($len+18)) {
				print OUT1 "\@ST-$lines[0]\n".substr($lines[1],0,(length($lines[1])-$len))."\n$lines[2]\n".substr($lines[3],0,(length($lines[1])-$len))."\n";
				print OUT2 "ST-$lables[0]\t".$len."\t".$1.$2.$3.$4.$5.$6.$7.$8."\n";
		}
	}
	elsif ($lines[1]=~/(A+)(A|T|C|G)(A+)(A|T|C|G)(A+)(A|G)(A+)(A|G)$/) {
			my $len=length($1)+length($2)+length($3)+length($4)+length($5)+length($6)+length($7)+length($8);
			my @lables=split(" ",$lines[0]);
			if (length($lines[1])>= ($len+18)) {
				print OUT1 "\@ST-$lines[0]\n".substr($lines[1],0,(length($lines[1])-$len))."\n$lines[2]\n".substr($lines[3],0,(length($lines[1])-$len))."\n";
				print OUT2 "ST-$lables[0]\t".$len."\t".$1.$2.$3.$4.$5.$6.$7.$8."\n";
		}
	}
	elsif ($lines[1]=~/(A+)(A|T|C|G)(A+)(A|T|C|G)(A+)(A|G)$/){
			my $len=length($1)+length($2)+length($3)+length($4)+length($5)+length($6);
			my @lables=split(" ",$lines[0]);
			if (length($lines[1])>= ($len+18)) {
				print OUT1 "\@ST-$lines[0]\n".substr($lines[1],0,(length($lines[1])-$len))."\n$lines[2]\n".substr($lines[3],0,(length($lines[1])-$len))."\n";
				print OUT2 "ST-$lables[0]\t".$len."\t".$1.$2.$3.$4.$5.$6."\n";
		}
	}
}
close IN;
close OUT1;
close OUT2;