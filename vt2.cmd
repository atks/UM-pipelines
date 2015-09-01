/net/fantasia/home/atks/dev/vt/pipeline/generate_vt2_calling_makefile.pl \
-o /net/fantasia/home/atks/topmed/20150901_freeze1_calling/run_3759samples \
-m 3759samples.topmed.vt.calling.mk \
-p 1000g \
-s 3759samples.index \
#window size of a chunk, use 0 if you want to process by sample
-w 10000000000 \
#sequence list with lengths, you can comment out sequences that you do not want in this list
-l hs37d5.seq.length.txt \
-r /net/topmed/working/atks/ref/hs37d5.fa
