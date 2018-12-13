# Download the cDNA reference from Ensembl plants
# Build bowtie2 index of cDNA sequence
bowtie2-build ../TAIR10/Sequence/Bowtie2Index-cDNA/TAIR10-cDNA-Sequence.fa

# THis maps the trimmed reads to the reference (allowing up to 2 mismatches)
bowtie2 -x ../TAIR10/Sequence/Bowtie2Index-cDNA/TAIR10-cDNA-Sequence -U Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq_delA.fq -S Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq_delA_cDNA.sam

# This discards the unmapped reads
less Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq_delA_cDNA.sam |awk '{if ($3!~/^*/) {if ($1~/^ST/)print}}' >Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq_delA_cDNA_mapped.sam