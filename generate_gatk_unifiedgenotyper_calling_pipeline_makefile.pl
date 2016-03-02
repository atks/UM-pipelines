#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_gatk_unifiedgenotyper_calling_pipeline_makefile

=head1 SYNOPSIS

 generate_gatk_unifiedgenotyper_pipeline_makefile [options]

  -s     sample file list giving the location of each sample
         column 1: sample name
         column 2: path of bam file
  -r     reference genome file
  -l     sequence length file
  -w     interval width
  -o     output directory
  -m     make file name

=head1 DESCRIPTION

This script generates the make file to discovery and genotype a set of individuals.

=cut

my $help;
my $outputDir = "run";
my $vtDir = "";
my $clusterDir = "";
my $makeFile = "Makefile";
my $cluster = "main";
my $sleep = 0;
my $sampleFile = "";
my $intervals = "";
my $sequenceLengthFile = "/net/fantasia/home/atks/dev/vt/pipeline/seq_length.txt";
my $intervalWidth = 1000000;
my $refGenomeFASTAFile = "";
my $jvmMemory = "2g";
my $variantType = "BOTH";
my $rawCopy = 0;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o:s'=>\$outputDir,
                'b:s'=>\$vtDir,
                't:s'=>\$clusterDir,
                'm:s'=>\$makeFile,
                'c:s'=>\$cluster,
                'd:s'=>\$sleep,
                's:s'=>\$sampleFile,
                'l:s'=>\$sequenceLengthFile,
                'i:s'=>\$intervalWidth,
                'r:s'=>\$refGenomeFASTAFile,
                'j:s'=>\$jvmMemory,
                'v:s'=>\$variantType,
                'x'=>\$rawCopy
                )
  || !defined($workDir)
  || !defined($makeFile)
  || !defined($sampleFile)
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
my $gatk = "/net/fantasia/home/atks/dev/vt/comparisons/programs/jdk1.7.0_25/bin/java -jar -Xmx$jvmMemory /net/fantasia/home/atks/dev/vt/comparisons/programs/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar";
my $vt = "$vtDir/vt";

printf("generate_gatk_unifiedgenotyper_calling_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         vt path              %s\n", $vt);
printf("         cluster path         %s\n", $clusterDir);
printf("         make file            %s\n", $makeFile);
printf("         cluster              %s\n", $cluster);
printf("         sleep                %s\n", $sleep);
printf("         sample file          %s\n", $sampleFile);
printf("         sequence length file %s\n", $sequenceLengthFile);
printf("         interval width       %s\n", $intervalWidth);
printf("         reference            %s\n", $refGenomeFASTAFile);
printf("         JVM Memory           %s\n", $jvmMemory);
printf("         variant types        %s\n", $variantType);
printf("         raw Copy             %s\n", $rawCopy);
printf("\n");

my $vcfOutDir = "$outputDir/vcf";
mkpath($vcfOutDir);
my $finalVCFOutDir = "$outputDir/final";
mkpath($finalVCFOutDir);
my $statsDir = "$outputDir/stats";
mkpath($statsDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $rawCopyDir = "$outputDir/raw";
if ($rawCopy)
{
    mkpath($rawCopyDir);
}
my $logFile = "$outputDir/run.log";

########################################
#Read file locations and name of samples
########################################
my %SAMPLE = ();
open(SA,"$sampleFile") || die "Cannot open $sampleFile\n";
my $bamFiles = "";
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $bamPath) = split(/\s+/, $_);
        $SAMPLE{$sampleID} = $bamPath;
        $bamFiles .= "$bamPath\n";
    }
}
close(SA);

my $bamListFile = "$auxDir/bam.list";
open(OUT,">$bamListFile") || die "Cannot open $bamListFile\n";
print OUT $bamFiles;
close(OUT);

print "read in " . scalar(keys(%SAMPLE)) . " samples\n";

###################
#Generate intervals
###################
my %intervalsByChrom = ();
my @intervals = ();
my @intervalFiles = ();
my @CHROM = ();

my $writeIntervals = 1;

if (-e "$outputDir/intervals/$intervalWidth.OK")
{
    print "$outputDir/intervals/$intervalWidth.OK exists, intervals wil not be generated.\n";
    $writeIntervals = 0;
}

mkpath("$outputDir/intervals/");
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
            my $file = "";
            if ($i<floor($len/$intervalWidth))
            {
                $interval = $chrom . "_" . ($intervalWidth*$i+1) . "_" . ($intervalWidth*($i+1));
                $file = "$outputDir/intervals/$interval.interval_list";
                if ($writeIntervals)
                {
                    open(INTERVAL, ">$file") || die "Cannot open $file\n";
                    print INTERVAL "$chrom:" . ($intervalWidth*$i+1) . "-" . ($intervalWidth*($i+1)) . "\n";
                    close(INTERVAL);
                }
            }
            elsif ($i*$intervalWidth!=$len)
            {
                $interval = $chrom . "_" . ($intervalWidth*$i+1) . "_" . $len;
                $file = "$outputDir/intervals/$interval.interval_list";
                if ($writeIntervals)
                {
                    open(INTERVAL, ">$file") || die "Cannot open $file\n";
                    print INTERVAL "$chrom:" . ($intervalWidth*$i+1) . "-" . $len . "\n";
                    close(INTERVAL);
                }
            }
            else
            {
                last;
            }
            
            push(@{$intervalsByChrom{$chrom}}, "$interval");
            push(@intervals, $interval);
            push(@intervalFiles, $file);

            $count++;
        }

        print "added $count intervals\n";
    }
}
close(SQ);

if ($writeIntervals)
{
    print `touch $outputDir/intervals/$intervalWidth.OK`;
}

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;
my $inputVCFFile;
my $outputVCFFile;

########
#Calling
########

#**************
#log start time
#**************
$tgt = "$logDir/start.calling.OK";
$dep = "";
@cmd = ("date | awk '{print \"gatk unifiedgenotyper variant calling pipeline\\n\\nstart calling: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

if ($intervalWidth!=0)
{
    my $intervalVCFFilesOK = "";
    for my $i (0 .. $#intervals)
    {
        #nct - number of computing threads
        #interval_padding ensures that you capture Indels that lie across a boundary. Note that UnifiedGenotyper uses locuswalker.
        #--max_alternate_alleles is set at 6 by default
        $outputVCFFile = "$vcfOutDir/$intervals[$i].genotypes.vcf";
        $tgt = "$outputVCFFile.OK";
        $dep = "";
        @cmd = ("$gatk -T UnifiedGenotyper -R $refGenomeFASTAFile -glm $variantType --interval_padding 100 -I $bamListFile --genotyping_mode DISCOVERY -o $outputVCFFile --output_mode EMIT_VARIANTS_ONLY -L $intervalFiles[$i]"),
        makeStep($tgt, $dep, @cmd);

        $intervalVCFFilesOK .= " $outputVCFFile.OK";
    }

    #************
    #log end time
    #************
    $tgt = "$logDir/end.calling.OK";
    $dep = "$intervalVCFFilesOK";
    @cmd = ("date | awk '{print \"end calling: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}
else
{
    $outputVCFFile = "$vcfOutDir/all.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "";
    @cmd = ("$gatk -T UnifiedGenotyper -R $refGenomeFASTAFile -glm $variantType -I $bamListFile --genotyping_mode DISCOVERY -o $outputVCFFile --output_mode EMIT_VARIANTS_ONLY");
    makeStep($tgt, $dep, @cmd);

    #************
    #log end time
    #************
    $tgt = "$logDir/end.calling.OK";
    $dep = "$outputVCFFile.OK";
    @cmd = ("date | awk '{print \"end calling: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}

###########################################
#Concatenate, normalize and drop duplicates
###########################################

if ($intervalWidth!=0)
{
    #**************
    #log start time
    #**************
    $tgt = "$logDir/start.concatenating.normalizing.OK";
    $dep = "$logDir/end.calling.OK";
    @cmd = ("date | awk '{print \"start concatenating and normalizing: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);

    my $chromSiteVCFFiles = "";
    my $chromSiteVCFFilesOK = "";
    my $chromSiteVCFIndicesOK = "";
    for my $chrom (@CHROM)
    {
        my $vcfListFile = "$auxDir/$chrom.vcf.list";
        open(OUT,">$vcfListFile") || die "Cannot open $vcfListFile\n";
        for my $interval (@{$intervalsByChrom{$chrom}})
        {
            print OUT "$vcfOutDir/$interval.genotypes.vcf\n";
        }
        close(OUT);
        
        #genotypes VCFs
        $outputVCFFile = "$finalVCFOutDir/$chrom.genotypes.vcf.gz";
        $tgt = "$outputVCFFile.OK";
        $dep = "$logDir/end.calling.OK";
        @cmd = ("$vt cat -L $vcfListFile -o + -w 1000 | $vt normalize + -o + -r $refGenomeFASTAFile 2> $statsDir/$chrom.normalize.log | $vt uniq + -o $outputVCFFile 2> $statsDir/$chrom.uniq.log");
        makeStep($tgt, $dep, @cmd);

        $inputVCFFile = "$finalVCFOutDir/$chrom.genotypes.vcf.gz";
        $tgt = "$inputVCFFile.tbi.OK";
        $dep = "$inputVCFFile.OK";
        @cmd = ("$vt index $inputVCFFile");
        makeStep($tgt, $dep, @cmd);

        #sites VCFs
        $inputVCFFile = "$finalVCFOutDir/$chrom.genotypes.vcf.gz";
        $outputVCFFile = "$finalVCFOutDir/$chrom.sites.vcf.gz";
        $tgt = "$outputVCFFile.OK";
        $dep = "$inputVCFFile.OK";
        @cmd = ("$vt view -s $inputVCFFile -o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);

        $inputVCFFile = "$finalVCFOutDir/$chrom.sites.vcf.gz";
        $tgt = "$inputVCFFile.tbi.OK";
        $dep = "$inputVCFFile.OK";
        @cmd = ("$vt index $inputVCFFile");
        makeStep($tgt, $dep, @cmd);
    }

    my $inputVCFFiles = join(" ", map {"$finalVCFOutDir/$_.sites.vcf.gz"} @CHROM);
    my $inputVCFFilesOK = join(" ", map {"$finalVCFOutDir/$_.sites.vcf.gz.OK"} @CHROM);
    $outputVCFFile = "$finalVCFOutDir/all.sites.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "$inputVCFFilesOK";
    @cmd = ("$vt cat $inputVCFFiles -o $outputVCFFile");
    makeStep($tgt, $dep, @cmd);

    $inputVCFFile = "$finalVCFOutDir/all.sites.vcf.gz";
    $tgt = "$inputVCFFile.tbi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeStep($tgt, $dep, @cmd);

    #************
    #log end time
    #************
    my $inputVCFFileIndicesOK = join(" ", map {"$finalVCFOutDir/$_.sites.vcf.gz.tbi.OK"} @CHROM);
    $tgt = "$logDir/end.concatenating.normalizing.OK";
    $dep = $inputVCFFileIndicesOK;
    @cmd = ("date | awk '{print \"end concatenating and normalizing: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
    
    if ($rawCopy)
    {
        for my $chrom (@CHROM)
        {
            my $vcfListFile = "$auxDir/$chrom.vcf.list";
            $outputVCFFile = "$rawCopyDir/$chrom.genotypes.vcf.gz";
            $tgt = "$outputVCFFile.OK";
            $dep = "$logDir/end.calling.OK";
            @cmd = ("$vt cat -L $vcfListFile -o + -w 1000 | $vt uniq + -o $outputVCFFile 2> $statsDir/$chrom.raw.uniq.log");
            makeStep($tgt, $dep, @cmd);

            $tgt = "$outputVCFFile.tbi.OK";
            $dep = "$outputVCFFile.OK";
            @cmd = ("$vt index $outputVCFFile");
            makeStep($tgt, $dep, @cmd);

            $inputVCFFile = "$rawCopyDir/$chrom.genotypes.vcf.gz";
            $outputVCFFile = "$rawCopyDir/$chrom.sites.vcf.gz";
            $tgt = "$outputVCFFile.OK";
            $dep = "$inputVCFFile.OK";
            @cmd = ("$vt view -s $inputVCFFile -o $outputVCFFile");
            makeStep($tgt, $dep, @cmd);

            $inputVCFFile = "$rawCopyDir/$chrom.sites.vcf.gz";
            $outputVCFFile = "$rawCopyDir/$chrom.sites.vcf.gz.tbi";
            $tgt = "$outputVCFFile.OK";
            $dep = "$inputVCFFile.OK";
            @cmd = ("$vt index $inputVCFFile");
            makeStep($tgt, $dep, @cmd);
        }

        my $inputVCFFiles = join(" ", map {"$rawCopyDir/$_.sites.vcf.gz"} @CHROM);
        my $inputVCFFilesOK = join(" ", map {"$rawCopyDir/$_.sites.vcf.gz.OK"} @CHROM);
        $outputVCFFile = "$rawCopyDir/all.sites.vcf.gz";
        $tgt = "$outputVCFFile.OK";
        $dep = "$inputVCFFilesOK";
        @cmd = ("$vt cat $inputVCFFiles -o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);

        $inputVCFFile = "$rawCopyDir/all.sites.vcf.gz";
        $tgt = "$inputVCFFile.tbi.OK";
        $dep = "$inputVCFFile.OK";
        @cmd = ("$vt index $inputVCFFile");
        makeStep($tgt, $dep, @cmd);
    }
}
else
{
    #**********************************
    #log start time for normalizing VCF
    #**********************************
    $tgt = "$logDir/start.normalization.OK";
    $dep = "$logDir/end.calling.OK";
    @cmd = ("date | awk '{print \"start normalization: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);

    $inputVCFFile = "$vcfOutDir/all.vcf";
    $outputVCFFile = "$finalVCFOutDir/all.genotypes.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "$logDir/end.genotyping.OK";
    @cmd = ("$vt normalize -r $refGenomeFASTAFile $inputVCFFile -o + 2> $statsDir/all.normalize.log | $vt uniq + -o $outputVCFFile 2> $statsDir/all.uniq.log");
    makeStep($tgt, $dep, @cmd);

    $inputVCFFile = "$finalVCFOutDir/all.genotypes.vcf.gz";
    $tgt = "$inputVCFFile.tbi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeStep($tgt, $dep, @cmd);

    $inputVCFFile = "$finalVCFOutDir/all.genotypes.vcf.gz";
    $outputVCFFile = "$finalVCFOutDir/all.sites.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt view -s $inputVCFFile -o $outputVCFFile");
    makeStep($tgt, $dep, @cmd);

    $inputVCFFile = "$finalVCFOutDir/all.sites.vcf.gz";
    $tgt = "$inputVCFFile.tbi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeStep($tgt, $dep, @cmd);

    #********************************
    #log end time for normalizing VCF
    #********************************
    $tgt = "$logDir/end.normalization.OK";
    $dep = "$inputVCFFile.tbi.OK";
    @cmd = ("date | awk '{print \"end normalization: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}

#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/*.* $vcfOutDir/*.* $vcfOutDir/*/*.* $finalVCFOutDir/*.* $statsDir/* $logDir/* $outputDir/intervals/*.*");

for(my $i=0; $i < @tgts; ++$i) {
    print MAK "$tgts[$i] : $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#Functions
##########

#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, @cmd) = @_;

    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    else
    {
        makeSlurm($partition, $tgt, $dep, @cmd);
    }
}

#run slurm jobs
sub makeSlurm
{
    my ($partition, $tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        #contains pipe
        if ($c=~/\|/)
        {
            ++$slurmScriptNo;
            my $slurmScriptFile = "$slurmScriptsDir/$slurmScriptNo.sh";
            open(IN, ">$slurmScriptFile");
            print IN "#!/bin/bash\n";
            print IN "set -o pipefail; $c";
            close(IN);
            chmod(0755, $slurmScriptFile);

            $cmd .= "\techo '" . $c . "'\n";
            $cmd .= "\tsrun -p $partition $slurmScriptFile\n";
        }
        else
        {
            $cmd .= "\tsrun -p $partition " . $c . "\n";
        }
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\tset -o pipefail; " . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

#run a local phony job
sub makePhonyJob
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    push(@cmds, $cmd);
}