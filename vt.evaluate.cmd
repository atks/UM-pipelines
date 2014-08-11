/net/fantasia/home/atks/dev/vt/pipeline/generate_variant_evaluation_pipeline_makefile.pl \
-g /net/fantasia/home/atks/bipolar/20140807_low_coverage_evaluation/bipolar.low_coverage.genotype.vcf.fileslist.txt \
-o /net/fantasia/home/atks/bipolar/20140807_low_coverage_evaluation \
-b /net/fantasia/home/atks/dev/vt \
-t /net/fantasia/home/atks/programs/cluster \
-m bipolar.low_coverage.gatk.ug.indels.evaluation.mk \
-c mini+ \
-d 25 \
#-s sample.sa \
#-p pedigree.ped \
-r /net/bipolar/atks/ref/hs37d5.fa \
-v INDEL
