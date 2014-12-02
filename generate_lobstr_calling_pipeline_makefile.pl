#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_lobstr_calling_pipeline_makefile

=head1 SYNOPSIS

 generate_lobstr_pipeline_makefile [options]

  -f     fastq file list giving the location of each sample
         column 1: fastq file containing
         column 2: path of bam file
  -r     reference genome file
  -l     sequence length file
  -w     interval width
  -o     output directory
  -m     make file name

=head1 DESCRIPTION

This script generates the make file to discovery and genotype a set of individuals for lobSTR.

=cut

my $help;

my $outputDir;
my $vtDir;
my $clusterDir;
my $makeFile;
my $cluster;
my $sleep;
my $fastqListFile;
my $sequenceLengthFile;
my $intervalWidth = 1000000;
my $refGenomeFASTAFile;
my $variantType;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o:s'=>\$outputDir,
                'b:s'=>\$vtDir,
                't:s'=>\$clusterDir,
                'm:s'=>\$makeFile,
                'c:s'=>\$cluster,
                'd:s'=>\$sleep,
                'f:s'=>\$fastqListFile,
                'l:s'=>\$sequenceLengthFile,
                'i:s'=>\$intervalWidth,
                'r:s'=>\$refGenomeFASTAFile
                )
  || !defined($makeFile)
  || !defined($fastqListFile)
  || !defined($refGenomeFASTAFile))
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

#programs
#you can set the  maximum memory here to be whatever you want
my $samtools = "/net/fantasia/home/atks/dev/vt/comparisons/programs/samtools/samtools";
my $bcftools = "/net/fantasia/home/atks/dev/vt/comparisons/programs/bcftools/bcftools";
my $vt = "$vtDir/vt";

my $lobstr = "/net/fantasia/home/atks/dev/vt/comparisons/programs/lobSTR-3.0.2/bin/lobSTR";
my $lobstr_allelotype = "/net/fantasia/home/atks/dev/vt/comparisons/programs/lobSTR-3.0.2/bin/allelotype";
my $lobstrResourcePrefix = "/net/fantasia/home/atks/dev/vt/comparisons/programs/lobSTR-3.0.2/resource/lobstr_v3.0.2_hg19_ref/lobSTR_";
my $lobstrResourceChromSizeTabFile = "$lobstrResourcePrefix/chromsizes.tab";
my $lobstrSTRInfo = "/net/fantasia/home/atks/dev/vt/comparisons/programs/lobSTR-3.0.2/resource/lobstr_v3.0.2_hg19_strinfo.tab";
my $lobstrSTRPCRFreeModel = "/net/fantasia/home/atks/dev/vt/comparisons/programs/lobSTR-3.0.2/share/lobSTR/models/illumina_v3.pcrfree";

printf("generate_samtools_calling_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         vt path              %s\n", $vt);
printf("         cluster path         %s\n", $clusterDir);
printf("         make file            %s\n", $makeFile);
printf("         cluster              %s\n", $cluster);
printf("         sleep                %s\n", $sleep);
printf("         sample file          %s\n", $fastqListFile);
printf("         sequence length file %s\n", $sequenceLengthFile);
printf("         interval width       %s\n", $intervalWidth);
printf("         reference            %s\n", $refGenomeFASTAFile);
printf("\n");

my $bamDir = "$outputDir/bam";
mkpath($bamDir);
my $vcfOutDir = "$outputDir/vcf";
mkpath($vcfOutDir);
my $finalDir = "$outputDir/final";
mkpath($finalDir);
my $statsDir = "$outputDir/stats";
mkpath($statsDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $logFile = "$outputDir/run.log";

########################################
#Read file locations and name of samples
########################################
my %SAMPLE = ();
open(SA,"$fastqListFile") || die "Cannot open $fastqListFile\n";
my @samples = ();
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $readGroup, $fastq1Path, $fastq2Path) = split(/\s+/, $_);

        if (!exists($SAMPLE{$sampleID}))
        {
            $SAMPLE{$sampleID}{FASTQ1} = ();
            $SAMPLE{$sampleID}{FASTQ1} = ();
            $SAMPLE{$sampleID}{READGROUP} = ();
            $SAMPLE{$sampleID}{BAM_PREFIX} = ();
            $SAMPLE{$sampleID}{N} = 0;
            push(@samples, $sampleID);
        }

        my $no = $SAMPLE{$sampleID}{N}+1;

        push(@{$SAMPLE{$sampleID}{FASTQ1}}, $fastq1Path);
        push(@{$SAMPLE{$sampleID}{FASTQ2}}, $fastq2Path);
        push(@{$SAMPLE{$sampleID}{READGROUP}}, $readGroup);
        my ($name, $path, $suffix) = fileparse($fastq1Path, (".fastq.gz", ".fastq"));
        my $bamFilePrefix = $fastq1Path;
        $bamFilePrefix = "$bamDir/$name";
        push(@{$SAMPLE{$sampleID}{BAM_PREFIX}}, $bamFilePrefix);
        ++$SAMPLE{$sampleID}{N};
    }
}
close(SA);

#my $bamListFile = "$auxDir/bam.list";
#open(OUT,">$bamListFile") || die "Cannot open $bamListFile\n";
#print OUT $bamFiles;
#close(OUT);

print "read in " . scalar(keys(%SAMPLE)) . " samples\n";

###################
#Generate intervals
###################
my %intervalsByChrom = ();
my @intervals = ();
my @intervalNames = ();
my @intervalFiles = ();
my @CHROM = ();

open(SQ,"$sequenceLengthFile") || die "Cannot open $sequenceLengthFile\n";
while (<SQ>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $len) = split('\t', $_);

        print "processing $chrom\t$len ";

        push(@CHROM, $chrom);

        $intervalsByChrom{$chrom} = ();
        my $count = 0;
        for my $i (0 .. floor($len/$intervalWidth))
        {
            my $interval = "";
            my $intervalName = "";
            my $file = "";
            if ($i<floor($len/$intervalWidth))
            {
                $interval = $chrom . ":" . ($intervalWidth*$i+1) . "-" . ($intervalWidth*($i+1));
                $intervalName = $chrom . "_" . ($intervalWidth*$i+1) . "_" . ($intervalWidth*($i+1));
            }
            else
            {
                $interval = $chrom . ":" . ($intervalWidth*$i+1) . "-" . $len;
                $intervalName = $chrom . "_" . ($intervalWidth*$i+1) . "_" . $len;
            }

            push(@{$intervalsByChrom{$chrom}}, "$intervalName");
            push(@intervals, $interval);
            push(@intervalNames, $intervalName);
            push(@intervalFiles, $file);

            $count++;
        }

        print "added $count intervals\n";
    }
}
close(SQ);

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;
my $inputVCFFile;
my $outputVCFFile;

##########
#Alignment
##########

#****************************
#log start time for alignment
#****************************
$tgt = "$logDir/start.alignment.OK";
$dep = "";
@cmd = ("date | awk '{print \"lobSTR variant calling pipeline\\n\\nstart calling: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

my $bamFiles = "";
my $bamFilesOK = "";
my $sampleBAMFileIndicesOK = "";

for my $sampleID (@samples)
{
    my @fastq1 = @{$SAMPLE{$sampleID}{FASTQ1}};
    my @fastq2 = @{$SAMPLE{$sampleID}{FASTQ2}};
    my @readGroups = @{$SAMPLE{$sampleID}{READGROUP}};
    my @bam_prefixes = @{$SAMPLE{$sampleID}{BAM_PREFIX}};

    for my $i (0 .. $#fastq1)
    {
        my $outputBAMFilePrefix = "$bam_prefixes[$i]";
        my $outputBAMFile = "$bam_prefixes[$i].aligned.bam";
        $tgt = "$outputBAMFile.OK";
        $dep = "$fastq1[$i] $fastq2[$i]";
        @cmd = ("$lobstr --gzip -q --p1 $fastq1[$i] --p2 $fastq2[$i] --index-prefix $lobstrResourcePrefix --rg-sample NA12878 --rg-lib $readGroups[$i] -o $outputBAMFilePrefix");
        makeStep($tgt, $dep, @cmd);

        $bamFiles .= " $outputBAMFile";
        $bamFilesOK .= " $outputBAMFile.OK";
    }

    #generate new header  - write a seprate script to do this, the SQ order must follow the BAM file!
    my $newHeaderSAMFile = "$bamDir/$sampleID.sam";
    open(OUT,">$newHeaderSAMFile") || die "Cannot open $newHeaderSAMFile\n";
    open(CHROM,"$lobstrResourceChromSizeTabFile") || die "Cannot open $lobstrResourceChromSizeTabFile\n";
    
    print OUT "@HD\tVN:1.3\n";
    while(<CHROM>)
    {
        chomp;
        my ($chrom, $len) = split("\t");
        print OUT "@SQ\tSN:$chrom\tLN:$len\n";
    }
    close(CHROM);    
    for my $rg(@readGroups)
    {
        print OUT "@RG\tID:lobSTR;NA12878;$rg\tLB:$rg\tSM:NA12878\n";
    }
    close(OUT);

    #concatenate and sort
    my $outputBAMFilePrefix = "$bamDir/$sampleID";
    my $outputBAMFile = "$bamDir/$sampleID.bam";
    $tgt = "$outputBAMFile.OK";
    $dep = "$bamFilesOK";
    @cmd = ("$samtools cat -h $newHeaderSAMFile -o - $bamFiles | $samtools sort - $outputBAMFilePrefix"),
    makeStep($tgt, $dep, @cmd);
    
    #index
    my $inputBAMFile = "$bamDir/$sampleID.bam";
    $tgt = "$inputBAMFile.bai.OK";
    $dep = "$inputBAMFile.OK";
    @cmd = ("$samtools index $inputBAMFile"),
    makeStep($tgt, $dep, @cmd);
    
    $sampleBAMFileIndicesOK = " $inputBAMFile.bai.OK";
}

#**************************
#log end time for alignment
#**************************
$tgt = "$logDir/end.alignment.OK";
$dep = "$sampleBAMFileIndicesOK";
@cmd = ("date | awk '{print \"end: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

#############
#Allelotyping
#############

#*******************************
#log start time for allelotyping
#*******************************
$tgt = "$logDir/start.allelotyping.OK";
$dep = "$logDir/end.alignment.OK";
@cmd = ("date | awk '{print \"start allelotyping: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

my $inputBAMFiles = join(",", map {"$bamDir/$_.bam"} @samples);
my $outputVCFFilePrefix = "$finalDir/all";
$tgt = "$outputVCFFilePrefix.vcf.gz.OK";
$dep = join(" ", map {"$bamDir/$_.bam.OK"} @samples);
@cmd = ("$lobstr_allelotype --command classify --bam $inputBAMFiles --index-prefix $lobstrResourcePrefix --strinfo $lobstrSTRInfo --noise_model $lobstrSTRPCRFreeModel --out $outputVCFFilePrefix");
makeStep($tgt, $dep, @cmd);

#index
$inputVCFFile = "$finalDir/all.vcf.gz";
$tgt = "$inputVCFFile.tbi.OK";
$dep = "$inputVCFFile.OK";
@cmd = ("$vt index $inputVCFFile"),
makeStep($tgt, $dep, @cmd);

#add post processing to add sequences, fix chromosome names (should be fixed based on resource files), REF not to be empty and sorting of sequences.


#*****************************
#log end time for allelotyping
#*****************************
$tgt = "$logDir/end.allelotyping.OK";
$dep = "$bamFilesOK";
@cmd = ("date | awk '{print \"end: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/*.* $vcfOutDir/*.* $vcfOutDir/*/*.* $finalDir/*.* $statsDir/* $logDir/* $outputDir/intervals/*.*");

for(my $i=0; $i < @tgts; ++$i) {
    print MAK "$tgts[$i] : $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#Functions
##########
sub makeMos
{
    my $cmd = shift;

    if ($cluster eq "main")
    {
        return ("mosbatch -E/tmp -i -r`$clusterDir/pick_main_node $sleep` /bin/bash -c 'set pipefail; $cmd'");
    }
    elsif ($cluster eq "mini")
    {
        return ("mosbatch -E/tmp -i -r`$clusterDir/pick_mini_node $sleep` /bin/bash -c 'set pipefail; $cmd'");
    }
    elsif ($cluster eq "mini+")
    {
        return ("mosbatch -E/tmp -i -r`$clusterDir/pick_mini+_node $sleep` /bin/bash -c 'set pipefail; $cmd'");
    }
    else
    {
        print STDERR "$cluster not supported\n";
        exit(1);
    }
}

sub makeStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . makeMos($c) . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}