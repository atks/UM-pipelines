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
    
  -r      reference genome
  -o 	  output directory
  -m      output make file
  -c      cluster
   
 example: ./generate.mapping.pipeline.pl --fastq /net/fantasia/home/atks/neptune/20120702_pilot/neptune.35samples.fastq.list -r /net/fantasia/home/atks/ref/genome/hs37d5.fa.gz -b /net/fantasia/home/atks/neptune/20120702_pilot/codes/bin -o /net/fantasia/home/atks/neptune/20120702_pilot/run -m neptune.35samples.mk --dbsnp /net/fantasia/home/atks/ref/dbsnp/dbsnp_135.b37.vcf -c 1000g
 
=head1 DESCRIPTION
 
=cut

#option variables
my $help;
my $verbose;
my $debug;
my $mainNodes = "30,32-36,38-41,43,44,46-52,120,122-124,126,129,131-134,140,141,144-152,154,156-159,165-167";

my %SAMPLE = ();

#reference files
my $refGenomeFile = "/net/assembly/atks/ref_mem/hs37d5.fa.gz";
my $knownIndelsVCFFile1 = "/net/assembly/atks/ref/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz";
my $knownIndelsVCFFile2 = "/net/assembly/atks/ref/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz";
my $knownSNPsVCFFile = "/net/assembly/atks/ref/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz";
my $cluster = "mini";

#input/output options
my $rootDir = "/net/assembly/atks/neptune/20130529_run539mem";
my $sampleFastQFile = "$rootDir/neptune.24samples.run539.fastq.list.txt";
my $binDir = "$rootDir/codes/bin";
my $outputDir = "$rootDir/run";
my $makeFile = "neptune.run539.24samples.mk";

#program directories
my $bwa = "/net/fantasia/home/atks/programs/bwa-0.7.4/bwa";
my $samtools = "/net/fantasia/home/atks/programs/samtools-0.1.17/samtools";
my $gatk = "/usr/bin/java -jar /net/fantasia/home/atks/programs/gatk-v1.6-13-g91f02df/dist/GenomeAnalysisTK.jar";
my $mergeSamFiles = "/usr/bin/java -Xmx16g  -jar /net/fantasia/home/atks/programs/picard-tools-1.72/jar/MergeSamFiles.jar";
my $mark_duplicates = "/usr/bin/java -Xmx16g  -jar /net/fantasia/home/atks/programs/picard-tools-1.72/jar/MarkDuplicates.jar";
my $verifybamid = "/net/fantasia/home/atks/programs/verifyBamID-20120620/verifyBamID";
my $verifybamidSiteVCF = "/net/assembly/atks/ref/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.vcf.gz";
my $qplot = "/net/fantasia/home/atks/programs/qplot/bin/qplot";
my $qplotDataDir = "/net/assembly/atks/ref/qplot/data";
my $bamtools = "/net/fantasia/home/atks/programs/bamtools_1.07/bin/bam";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug,
                'fastq:s'=>\$sampleFastQFile,
                'r:s'=>\$refGenomeFile,
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


my $refGenomeFAFile = "/net/assembly/atks/ref/hs37d5.fa";

mkdir($outputDir);

#read name of samples and fastQ file locations
open(SA,"$sampleFastQFile") || die "Cannot open $sampleFastQFile\n";
my %SAMPLE_NAMES = ();
my %LANES = ();
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $RG, $fastQFile1, $fastQFile2) = split('\t', $_);
        if (!exists($SAMPLE_NAMES{$sampleID}))
        {
            $SAMPLE_NAMES{$sampleID} = ();
        }
        
        push(@{$SAMPLE_NAMES{$sampleID}}, $RG);
        
        $RG =~ /L([1-8])/;
        my $lane = $1;
        
        if (!exists($LANES{$lane}))
        {
            $LANES{$lane} = ();
        }
        
        push(@{$LANES{$lane}}, "$sampleID.$RG");
        
        
        if (!exists($SAMPLE{"$sampleID.$RG"}))
        {        
            $SAMPLE{"$sampleID.$RG"}{FASTQ1} = $fastQFile1;
            $SAMPLE{"$sampleID.$RG"}{FASTQ2} = $fastQFile2;
            $SAMPLE{"$sampleID.$RG"}{RG} = "$RG";
        }
    }
}
close(SA);

my @chrs = (1..22,"X", "Y", "MT");
my @tgts = ();
my @deps = ();
my @cmds = ();
my $cmd = "";

my $noSamples = scalar(keys(%SAMPLE));
print STDERR "No of Samples: $noSamples\n";

###############################
#1. Index Genome Reference file
###############################
#4h
push(@tgts,"$refGenomeFile.bwt.OK");
push(@deps,"");
$cmd = "\t$bwa index -a bwtsw $refGenomeFile\n";
$cmd = "$cmd\ttouch $refGenomeFile.bwt.OK\n";
push(@cmds, $cmd);

############################
#2. Perform paired alignment
############################
my $bamDir = "$outputDir/bam";
mkdir($bamDir);
for my $sampleID (keys(%SAMPLE))
{
    #?
    push(@tgts,"$bamDir/$sampleID.sam.OK");
    push(@deps,"$refGenomeFile.bwt.OK $SAMPLE{$sampleID}{FASTQ1} $SAMPLE{$sampleID}{FASTQ2}");
    $cmd = "$bwa mem -t 2 -M $refGenomeFile $SAMPLE{$sampleID}{FASTQ1} $SAMPLE{$sampleID}{FASTQ2} > $bamDir/$sampleID.sam"; 
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.sam.OK\n";
    push(@cmds, $cmd);
}

##############
#4. Polish bam
##############
for my $sampleID (keys(%SAMPLE))
{
    #8m
    push(@tgts,"$bamDir/$sampleID.bam.OK");
    push(@deps,"$bamDir/$sampleID.sam.OK");
    $cmd = "$samtools view -bSu $bamDir/$sampleID.sam | $samtools sort -n -o - $bamDir/$sampleID.nsort_tmp | $samtools fixmate /dev/stdin /dev/stdout | $samtools sort -o - $bamDir/$sampleID.csort_tmp | $samtools fillmd -u - $refGenomeFAFile.gz > $bamDir/$sampleID.bam";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.bam.OK\n";
    push(@cmds, $cmd);
}

########################
#5. Strip Tag and add RG
########################
for my $sampleID (keys(%SAMPLE))
{
    #$sampleID = Sample_\d+.<RG>
    
    #2m
    push(@tgts,"$bamDir/$sampleID.tagged.bam.OK");
    push(@deps,"$bamDir/$sampleID.bam.OK");
    $cmd = "\t" . makeMos("$samtools view -h $bamDir/$sampleID.bam | $binDir/strip_tags_and_addRG --rg $sampleID | $samtools view -bS - >  $bamDir/$sampleID.tagged.bam") . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.tagged.bam.OK\n";
    push(@cmds, $cmd);

    #2s
    push(@tgts,"$bamDir/$sampleID.tagged.bam.bai.OK");
    push(@deps,"$bamDir/$sampleID.tagged.bam.OK");
    $cmd = "$samtools index $bamDir/$sampleID.tagged.bam";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.tagged.bam.bai.OK\n";
    push(@cmds, $cmd);
}

#####################
#6. Local Realignment
#####################
#62m
push(@tgts,"$bamDir/indel.intervals.OK");
push(@deps,"$refGenomeFAFile $knownIndelsVCFFile1 $knownIndelsVCFFile2");
$cmd = "$gatk -T RealignerTargetCreator -R $refGenomeFAFile -o $bamDir/indel.intervals -known $knownIndelsVCFFile1 -known $knownIndelsVCFFile2";
$cmd = "\t" . makeMosBest($cmd) . "\n";
$cmd .= "\ttouch $bamDir/indel.intervals.OK\n";
push(@cmds, $cmd);

for my $sampleID (keys(%SAMPLE))
{
    #3m
    push(@tgts,"$bamDir/$sampleID.tagged.realigned.bam.OK");
    push(@deps,"$bamDir/indel.intervals.OK $bamDir/$sampleID.tagged.bam.bai.OK");
    $cmd = "$gatk -T IndelRealigner -rf NotPrimaryAlignment -R $refGenomeFAFile -I $bamDir/$sampleID.tagged.bam -o $bamDir/$sampleID.tagged.realigned.bam -targetIntervals $bamDir/indel.intervals -known $knownIndelsVCFFile1 -known $knownIndelsVCFFile2 -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing";
    $cmd = "\t" . makeMosBest($cmd) . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.tagged.realigned.bam.OK\n";
    push(@cmds, $cmd);  
}

#########################
#7. Merge BAMs by Samples
#########################
for my $sampleID (keys(%SAMPLE_NAMES))
{
    my $input = "";
    my $inputOK = "";
    for my $rgTag (@{$SAMPLE_NAMES{$sampleID}})
    {
        $input .= "INPUT=$bamDir/$sampleID.$rgTag.tagged.realigned.bam ";
        $inputOK .= "$bamDir/$sampleID.$rgTag.tagged.realigned.bam.OK ";
    }
    
    #m
    push(@tgts,"$bamDir/$sampleID.bam.OK");
    push(@deps,"$inputOK");
    $cmd = "\t" . makeMos("$mergeSamFiles $input OUTPUT=$bamDir/$sampleID.bam VALIDATION_STRINGENCY=SILENT") . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.bam.OK\n";
    push(@cmds, $cmd);
}


#################
#8. Recalibration
#################
for my $sampleID (keys(%SAMPLE_NAMES))
{
    #190m
    push(@tgts,"$bamDir/$sampleID.dedupped.recal.bam.OK");
    push(@deps,"$bamDir/$sampleID.bam.OK");
    $cmd = "\t" . makeMos("$bamtools dedup --recab --in $bamDir/$sampleID.bam --out $bamDir/$sampleID.dedupped.recal.bam --force --refFile $refGenomeFAFile  --dbsnp $knownSNPsVCFFile --storeQualTag OQ --maxBaseQual 40") . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.dedupped.recal.bam.OK\n";
    push(@cmds, $cmd);
    
    #4min
    push(@tgts,"$bamDir/$sampleID.dedupped.recal.bam.bai.OK");
    push(@deps,"$bamDir/$sampleID.dedupped.recal.bam.OK");
    $cmd = "$samtools index $bamDir/$sampleID.dedupped.recal.bam";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\ttouch $bamDir/$sampleID.dedupped.recal.bam.bai.OK\n";
    push(@cmds, $cmd);    
}  

#########
#9. QPLOT
#########
my $qplotDir = "$outputDir/qplot";
my $refUMFAFile = $refGenomeFAFile;
$refUMFAFile =~ s/\.fa\.gz/.umfa/;
mkdir($qplotDir);
my $qplotOKFiles = "";
my $qplotRFiles = "";
for my $sampleID (keys(%SAMPLE_NAMES))
{
    $qplotOKFiles .= " $qplotDir/$sampleID.qplot.OK";
    $qplotRFiles .= " $qplotDir/$sampleID.R";
    
    #63m
    push(@tgts,"$qplotDir/$sampleID.qplot.OK");
    push(@deps,"$bamDir/$sampleID.dedupped.recal.bam.OK");
    $cmd =  "\t" . makeMos("$qplot --reference $refUMFAFile --dbsnp $qplotDataDir/dbSNP135.hs37d5.tbl  --gccontent $qplotDataDir/human.g1k.w100.gc --plot $qplotDir/$sampleID.pdf --stats $qplotDir/$sampleID.stats --Rcode $qplotDir/$sampleID.R --minMapQuality 0 --bamlabel $sampleID.recal,$sampleID $bamDir/$sampleID.dedupped.recal.bam $bamDir/$sampleID.bam") . "\n";
    $cmd .= "\ttouch $qplotDir/$sampleID.qplot.OK\n";
    push(@cmds, $cmd);
}

#########################
#10. plot combined qplots
#########################
#1s
push(@tgts,"$qplotDir/combined.pdf.OK");
push(@deps, $qplotOKFiles);
$cmd =  "\t" . makeMos("/net/fantasia/home/atks/programs/python_2.7/bin/python $binDir/drawPool.py -o $qplotDir/combined $qplotRFiles") . "\n";
$cmd .= "\ttouch $qplotDir/combined.pdf.OK\n";
push(@cmds, $cmd);

#   Stage                           jobs#  time/job ptime   concurrent#       
#1. Index Genome Reference for bwa  1      4h       4h      1 
#2. Index fastq files for bwa       814    25m      10.8h   32
#3. bwa pair alignment              407    7m       1.5h    32
#4. Polish bam with Samtools        407    8m       1.8h    32
#5. Add RG                          407    2m       26m     32
#6. Index bam files                 407    2s       26s     32
#7. GATK Intervals for Realignment  1      62m      62m     1
#8. GATK local realignment          407    3m       78m     16
#9. Merge bams by samples           24     25m      50m     12
#10.Index merged bam files          24     4m       48m     12
#11.UM Dedupping and Recalibration  24     190m     6.4h    12
#12.Index bam file                  24     4m       48m     12 
#13.Qplot recalibrated bams         24     63m      126m    12
#14.VerifyBamID recalibrated bams   24     100m     200m    12
#                                       Total time: 35.1h

###############
#8. VerifybamID
###############
my $verifyBamIDDir = "$outputDir/verifybamid";
mkdir($verifyBamIDDir);
my $selfSMFiles = "";
my $verifybamidOKFiles = "";
for my $sampleID (keys(%SAMPLE_NAMES))
{
    #16m
    push(@tgts,"$verifyBamIDDir/$sampleID.verifybamid.OK");
    push(@deps,"$bamDir/$sampleID.dedupped.recal.bam.bai.OK");
    $cmd =  "\t" . makeMos("$verifybamid --site --vcf $verifybamidSiteVCF --genoError 0.005 --bam $bamDir/$sampleID.dedupped.recal.bam --out $verifyBamIDDir/$sampleID.verifybamid") . "\n";
    $cmd .= "\ttouch $verifyBamIDDir/$sampleID.verifybamid.OK\n";
    push(@cmds, $cmd);
    $verifybamidOKFiles .= "$verifyBamIDDir/$sampleID.verifybamid.OK ";
	$selfSMFiles .= "$verifyBamIDDir/$sampleID.verifybamid.selfSM ";
}

#plot freemix values
my $plotDir = "$verifyBamIDDir/plot_sample";
push(@tgts,"$plotDir/freemix.pdf.OK");
push(@deps, $verifybamidOKFiles);
$cmd =  "\t" . makeMos("$binDir/plot_verifybamidoutput -p $plotDir $selfSMFiles") . "\n";
$cmd .= "\ttouch $plotDir/freemix.pdf.OK\n";
push(@cmds, $cmd);

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
#my $mainNodes = "30,32-36,38-44,46-52,120,121-135,140-169";
sub makeMos 
{
    my $cmd = shift;
    if ($cluster eq "mini")
    {
    	return ("mosbatch -E/tmp -i -j10,11,12,13 sh -c \"$cmd\"");
    }
    else
    {
    	return ("mosbatch -E/tmp -i -j$mainNodes sh -c \"$cmd\"");
	}
}

sub makeMosBest 
{
    my $cmd = shift;
    if ($cluster eq "mini")
    {
    	return ("mosbatch -E/tmp -i -b sh -c \"$cmd\"");
    }
    else
    {
    	return ("mosbatch -E/tmp -i -j$mainNodes sh -c \"$cmd\"");
	}
}
