#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_samtools_calling_pipeline_makefile

=head1 SYNOPSIS

 generate_samtools_pipeline_makefile [options]

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

my $outputDir;
my $vtDir;
my $clusterDir;
my $makeFile;
my $cluster;
my $sleep;
my $sampleFile;
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
                's:s'=>\$sampleFile,
                'l:s'=>\$sequenceLengthFile,
                'i:s'=>\$intervalWidth,
                'r:s'=>\$refGenomeFASTAFile
                )
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
my $samtools = "/net/fantasia/home/atks/dev/vt/comparisons/programs/samtools/samtools";
my $bcftools = "/net/fantasia/home/atks/dev/vt/comparisons/programs/bcftools/bcftools";
my $vt = "$vtDir/vt";

printf("generate_samtools_calling_pipeline_makefile.pl\n");
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
            elsif ($i*$intervalWidth!=$len)
            {
                $interval = $chrom . ":" . ($intervalWidth*$i+1) . "-" . $len;
                $intervalName = $chrom . "_" . ($intervalWidth*$i+1) . "_" . $len;
            }
            else
            {
                last;
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

########
#Calling
########

#**************
#log start time
#**************
$tgt = "$logDir/start.calling.OK";
$dep = "";
@cmd = ("date | awk '{print \"samtools variant calling pipeline\\n\\nstart calling: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

if ($intervalWidth!=0)
{
    my $intervalVCFFilesOK = "";
    for my $i (0 .. $#intervals)
    {
        $outputVCFFile = "$vcfOutDir/$intervalNames[$i].vcf.gz";
        $tgt = "$outputVCFFile.OK";
        $dep = "";
        @cmd = ("$samtools  mpileup -ugf $refGenomeFASTAFile -b $bamListFile -r $intervals[$i] | $bcftools call -vmO z -o $outputVCFFile"),
        makeStep($tgt, $dep, @cmd);

        $intervalVCFFilesOK .= " $outputVCFFile.OK";
    }

    #************
    #log end time
    #************
    $tgt = "$logDir/end.calling.OK";
    $dep = "$intervalVCFFilesOK";
    @cmd = ("date | awk '{print \"end: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}
else
{
    $outputVCFFile = "$vcfOutDir/all.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "";
    @cmd = ("$samtools  mpileup -ugf $refGenomeFASTAFile -b $bamListFile | $bcftools call -vmO z -o $outputVCFFile");
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
        $chromSiteVCFFiles .= " $finalVCFOutDir/$chrom.sites.vcf.gz";
        $chromSiteVCFFilesOK .= " $finalVCFOutDir/$chrom.sites.vcf.gz.OK";
        $chromSiteVCFIndicesOK .= " $finalVCFOutDir/$chrom.sites.vcf.gz.tbi.OK";

        my $inputChromosomeIntervalVCFFiles = "";
        for my $interval (@{$intervalsByChrom{$chrom}})
        {
            $inputChromosomeIntervalVCFFiles .= " $vcfOutDir/$interval.vcf.gz";
        }

        #genotypes VCFs
        $outputVCFFile = "$finalVCFOutDir/$chrom.genotypes.vcf.gz";
        $tgt = "$outputVCFFile.OK";
        $dep = "$logDir/end.calling.OK";
        @cmd = ("$vt cat $inputChromosomeIntervalVCFFiles -o + | $vt normalize + -o + -r $refGenomeFASTAFile 2> $statsDir/$chrom.normalize.log | $vt uniq + -o $outputVCFFile 2> $statsDir/$chrom.uniq.log");
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

    $outputVCFFile = "$finalVCFOutDir/all.sites.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "$chromSiteVCFFilesOK";
    @cmd = ("$vt cat $chromSiteVCFFiles -o $outputVCFFile");
    makeStep($tgt, $dep, @cmd);

    $inputVCFFile = "$finalVCFOutDir/all.sites.vcf.gz";
    $tgt = "$inputVCFFile.tbi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeStep($tgt, $dep, @cmd);

    #************
    #log end time
    #************
    $tgt = "$logDir/end.concatenating.normalizing.OK";
    $dep = "$chromSiteVCFIndicesOK";
    @cmd = ("date | awk '{print \"end concatenating and normalizing: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
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
    @cmd = ("$vt normalize -r $refGenomeFASTAFile $inputVCFFile -o + | $vt mergedups + -o $outputVCFFile ");
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