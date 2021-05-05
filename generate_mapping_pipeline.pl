#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_mapping_pipeline.pl

=head1 SYNOPSIS

 generate_mapping_pipeline.pl [options]  
    
 --fastq  paired fastq file locations
  -r      reference genome
  -o 	  output directory
  -m      output make file
   
 example: ./generate.mapping.pipeline.pl --fastq /net/fantasia/home/atks/neptune/20120702_pilot/neptune.35samples.fastq.list -r /net/fantasia/home/atks/ref/genome/hs37d5.fa.gz -b /net/fantasia/home/atks/neptune/20120702_pilot/codes/bin -o /net/fantasia/home/atks/neptune/20120702_pilot/run -m neptune.35samples.mk --dbsnp /net/fantasia/home/atks/ref/dbsnp/dbsnp_135.b37.vcf -c 1000g
 
=head1 DESCRIPTION
 
=cut

#option variables
my $help;
my $verbose;
my $debug;

my %SAMPLE = ();

my $sampleFastQFile; 
my $refGenomeFile = "/net/fantasia/home/atks/ref/genome/hs37d5.fa.gz";
my $dbsnpFile = "/net/fantasia/home/atks/ref/dbsnp/dbsnp_135.b37.vcf.gz";
my $outputDir;
my $binDir;
my $makeFile = "mindel.mk";
my $cluster = "main";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug,
                'fastq:s'=>\$sampleFastQFile,
                'r:s'=>\$refGenomeFile,
                'dbsnpr:s'=>\$dbsnpFile,
                'o:s'=>\$outputDir,
                'b:s'=>\$binDir,
                'm:s'=>\$makeFile,
                'c:s'=>\$cluster) 
  || !defined($sampleFastQFile)
  || !defined($refGenomeFile)
  || !defined($binDir) 
  || !defined($outputDir) 
  || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(-verbose => 2);
    }
    else
    {
        pod2usage(1);
    }
}

my $refGenomeFAFile = $refGenomeFile;
if ($refGenomeFile =~ /(.+)\.gz$/)
{
	$refGenomeFAFile = $1;
}

#program directories
my $bwa = "/net/fantasia/home/atks/programs/bwa-0.6.2/bwa";
my $samtools = "/net/fantasia/home/atks/programs/samtools-0.1.17/samtools";
my $gatk = "/usr/bin/java -jar /net/fantasia/home/atks/programs/gatk-v1.6-13-g91f02df/dist/GenomeAnalysisTK.jar";
my $gatk4g = "/usr/bin/java -Xmx4g -jar /net/fantasia/home/atks/programs/gatk-v1.6-13-g91f02df/dist/GenomeAnalysisTK.jar";
my $mark_duplicates = "/usr/bin/java -Xmx16g  -jar /net/fantasia/home/atks/programs/picard-tools-1.72/jar/MarkDuplicates.jar";
my $bam_clipoverlap = "/usr/cluster/bin/bam clipOverlap";
my $bam_polishbam = "/usr/cluster/bin/bam polishBam";
my $verifybamid = "/net/fantasia/home/atks/programs/verifyBamID-20120620/verifyBamID";
my $verifybamidSiteVCF = "/net/fantasia/home/atks/programs/verifyBamID-20120620/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.vcf.gz";
my $qplot = "/net/fantasia/home/atks/programs/qplot/bin/qplot";
my $qplotDataDir = "/net/fantasia/home/atks/programs/qplot/data";

#read name of samples and fastQ file locations
open(SA,"$sampleFastQFile") || die "Cannot open $sampleFastQFile\n";
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $fastQFile1, $fastQFile2) = split('\t', $_);
        $SAMPLE{$sampleID}{FASTQ1} = $fastQFile1;
        $SAMPLE{$sampleID}{FASTQ2} = $fastQFile2;
    }
}
close(SA);

my @chrs = (1..22,"X", "Y", "MT");
my @tgts = ();
my @deps = ();
my @cmds = ();

my $noSamples = scalar(keys(%SAMPLE));
print STDERR "No of Samples: $noSamples\n";

mkdir($outputDir);

###############################
#1. Index Genome Reference file
###############################
push(@tgts,"$refGenomeFile.bwt.OK");
push(@deps,"");
my $cmd = "\t$bwa index -a bwtsw $refGenomeFile\n";
$cmd = "$cmd\ttouch $refGenomeFile.bwt.OK\n";
push(@cmds, $cmd);

##########################################
#2. Get SA Coordinates for each FASTQ file
##########################################
my $saiDir = "$outputDir/sai";
mkdir($saiDir);
for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$saiDir/$sampleID.1.sai.OK");
    push(@deps,"$refGenomeFile.bwt.OK");
    $cmd = "$bwa aln $refGenomeFile $SAMPLE{$sampleID}{FASTQ1} > $saiDir/$sampleID.1.sai";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $saiDir/$sampleID.1.sai.OK\n";
    push(@cmds, $cmd);
    
    push(@tgts,"$saiDir/$sampleID.2.sai.OK");
    push(@deps,"$refGenomeFile.bwt.OK");
    $cmd = "$bwa aln $refGenomeFile $SAMPLE{$sampleID}{FASTQ2} > $saiDir/$sampleID.2.sai";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $saiDir/$sampleID.2.sai.OK\n";
    push(@cmds, $cmd);
}

############################
#3. Perform paired alignment
############################
my $bamDir = "$outputDir/bam";
mkdir($bamDir);
for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$bamDir/$sampleID.bam.OK");
    push(@deps,"$refGenomeFile.bwt.OK $saiDir/$sampleID.1.sai.OK $saiDir/$sampleID.2.sai.OK");
    $cmd = "$bwa sampe $refGenomeFile $saiDir/$sampleID.1.sai $saiDir/$sampleID.2.sai $SAMPLE{$sampleID}{FASTQ1} $SAMPLE{$sampleID}{FASTQ2} | $samtools view -uhS - | $samtools sort -m 1000000000 - $bamDir/$sampleID";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.bam.OK\n";
    push(@cmds, $cmd);
}

#####################
#4. Remove Duplicates
#####################
for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$bamDir/$sampleID.dedup.bam.OK");
    push(@deps,"$bamDir/$sampleID.bam.OK");
    $cmd = "$mark_duplicates INPUT= $bamDir/$sampleID.bam OUTPUT=$bamDir/$sampleID.dedup.bam AS=TRUE VALIDATION_STRINGENCY=SILENT METRICS_FILE=$bamDir/$sampleID.dedup.bam.metrics";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.dedup.bam.OK\n";
    push(@cmds, $cmd);
}

#################
#5. Clip Overlap
#################
for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$bamDir/$sampleID.dedup.clip.bam.OK");
    push(@deps,"$bamDir/$sampleID.dedup.bam.OK");
    $cmd = "$bam_clipoverlap --in $bamDir/$sampleID.dedup.bam --out $bamDir/$sampleID.dedup.clip.bam --stats";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.dedup.clip.bam.OK\n";
    push(@cmds, $cmd);
}

###############
#6. Polish Bams
###############
for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$bamDir/$sampleID.dedup.clip.polished.bam.OK");
    push(@deps,"$bamDir/$sampleID.dedup.clip.bam.OK");
    $cmd = "$bam_polishbam -f $refGenomeFAFile --AS NCBI37 --UR file:$refGenomeFAFile --checkSQ --RG \"\@RG\tID:$sampleID\tSM:1\tLB:123\tCN:UM\tPL:ILLUMINA\" -i $bamDir/$sampleID.dedup.clip.bam -o $bamDir/$sampleID.dedup.clip.polished.bam -l $bamDir/$sampleID.dedup.clip.polished.bam.log";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\t" . makeMos("$samtools index $bamDir/$sampleID.dedup.clip.polished.bam") . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.dedup.clip.polished.bam.OK\n";
    push(@cmds, $cmd);
}

###############
#7. Recalibrate
###############
for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$bamDir/$sampleID.recal.csv.OK");
    push(@deps,"$bamDir/$sampleID.dedup.clip.polished.bam.OK");
    $cmd =  "\t" . makeMos("$gatk4g -l INFO -T CountCovariates -R $refGenomeFAFile -knownSites $dbsnpFile -OQ -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate --default_platform illumina -U ALLOW_UNSET_BAM_SORT_ORDER -S SILENT -I $bamDir/$sampleID.dedup.clip.polished.bam -recalFile $bamDir/$sampleID.recal.csv") . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.recal.csv.OK\n";
    push(@cmds, $cmd);
}

for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$bamDir/$sampleID.dedup.clip.polished.recal.bam.OK");
    push(@deps,"$bamDir/$sampleID.recal.csv.OK");
    $cmd =  "\t" . makeMos("$gatk4g -l INFO -T TableRecalibration -R $refGenomeFAFile -I $bamDir/$sampleID.dedup.clip.polished.bam -recalFile $bamDir/$sampleID.recal.csv --out $bamDir/$sampleID.dedup.clip.polished.recal.bam") . "\n";
    $cmd .= "\t" . makeMos("$samtools index $bamDir/$sampleID.dedup.clip.polished.recal.bam") . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.dedup.clip.polished.recal.bam.OK\n";
    push(@cmds, $cmd);
}

#########
#8. QPLOT
#########
my $qplotDir = "$outputDir/qplot";
my $refUMFAFile = $refGenomeFile;
$refUMFAFile =~ s/\.fa\.gz/.umfa/;
mkdir($qplotDir);
for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$qplotDir/$sampleID.qplot.OK");
    push(@deps,"$bamDir/$sampleID.dedup.clip.polished.recal.bam.OK");
    $cmd =  "\t" . makeMos("$qplot --reference $refUMFAFile --dbsnp $qplotDataDir/dbSNP130.hg19.tbl --gccontent $qplotDataDir/human.g1k.w100.gc --plot $qplotDir/$sampleID.pdf --stats $qplotDir/$sampleID.stats --Rcode $qplotDir/$sampleID.R --minMapQuality 0 --bamlabel $sampleID $bamDir/$sampleID.dedup.clip.polished.recal.bam") . "\n";
    $cmd .= "\ttouch $qplotDir/$sampleID.qplot.OK\n";
    push(@cmds, $cmd);
}
	
###############
#9. VerifybamID
###############
my $verifyBamIDDir = "$outputDir/verifybamid";
mkdir($verifyBamIDDir);
my $selfSMFiles = "";
my $verifybamidOKFiles = "";
for my $sampleID (keys(%SAMPLE))
{
    push(@tgts,"$verifyBamIDDir/$sampleID.verifybamid.OK");
    push(@deps,"$bamDir/$sampleID.dedup.clip.polished.recal.bam.OK");
    $cmd =  "\t" . makeMos("$verifybamid --site --vcf $verifybamidSiteVCF --genoError 0.005 --bam $bamDir/$sampleID.dedup.clip.polished.recal.bam --out $verifyBamIDDir/$sampleID.verifybamid") . "\n";
    $cmd .= "\ttouch $verifyBamIDDir/$sampleID.verifybamid.OK\n";
    push(@cmds, $cmd);
    $verifybamidOKFiles .= "$verifyBamIDDir/$sampleID.verifybamid.OK ";
	$selfSMFiles .= "$verifyBamIDDir/$sampleID.verifybamid.selfSM ";
}

#plot freemix values
my $plotDir = "$verifyBamIDDir/plot";
push(@tgts,"$plotDir/freemix.pdf.OK");
push(@deps, $verifybamidOKFiles);
$cmd =  "\t" . makeMos("$binDir/plot_verifybamidoutput -p $plotDir $selfSMFiles") . "\n";
$cmd .= "\ttouch $plotDir/freemix.pdf.OK\n";
push(@cmds, $cmd);

#BWA for alignment for paired end reads
#Picard for Dedupping
#BamclipOverlap
#GATK for base recalibration
#QPLOT
#verifybamID

#bwa index -a bwtsw ~/ref/genome/hs37d5.fa.gz 
#4.4 hrs

#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";
for(my $i=0; $i < @tgts; ++$i) 
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
##########
sub makeMos 
{
    my $cmd = shift;
    if ($cluster eq "1000g")
    {
    	return ("mosrun -E/tmp -i -t -j10,11,12,13 sh -c '$cmd'");
    }
    else
    {
    	return ("mosrun -E/tmp -j10-24,26,30-52,120-170 sh -c '$cmd'");
    	#return ("mosrun -E/tmp -j10-24,26,30-52,120-170 sh -c '$cmd'");
	}
}
