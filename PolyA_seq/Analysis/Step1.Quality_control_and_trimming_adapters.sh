# This discards the poor-quality reads (>50% bases with a Phred score <20)
fastq_quality_filter -q 20 -p 50 -Q 33 -i Col_1_end1_raw.fq -o Col_1_end1_raw_q20p50Q33.fq

# This trims the 3'-adapter sequences and keeps only sequences which contained the adapter 
fastx_clipper -l 16 -a TGGAATTCTCGGG -c -Q 33 -i Col_1_end1_raw_q20p50Q33.fq -o Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter.fq

# This trims the 3-nt random sequence at the 5¡ä-end
fastx_trimmer -f 4 -Q 33 -i Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter.fq -o Col_1_end1_raw_q20p50Q33_l16aTGGAATTCTCGGG_adapter_del3.fq