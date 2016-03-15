#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_freebayes_calling_pipeline_makefile

=head1 SYNOPSIS

 generate_freebayes_pipeline_makefile [options]

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
my $makeFile = "Makefile";
my $partition = "local";
my $sampleFile = "";
my $intervalWidth = 1000000;
my $refGenomeFASTAFile = "";
my $slurmScriptsSubDir = "";
my $slurmScriptNo = 0;


#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o:s'=>\$outputDir,
                'm:s'=>\$makeFile,
                'p:s'=>\$partition,
                's:s'=>\$sampleFile,
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
my $freebayes = "/net/fantasia/home/atks/dev/vt/comparisons/programs/freebayes/bin/freebayes-1.0.2-6-g3ce827d";
my $vt = "/net/fantasia/home/atks/dev/vt/comparisons/programs/vt/vt";
my $injectContigs = "/net/fantasia/home/atks/dev/vt/comparisons/programs/scripts/inject_contigs";
my $contigsFile = "/net/fantasia/home/atks/dev/vt/comparisons/NA12878/indices/contigs.txt";

printf("generate_freebayes_calling_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         vt path              %s\n", $vt);
printf("         make file            %s\n", $makeFile);
printf("         sample file          %s\n", $sampleFile);
printf("         interval width       %s\n", $intervalWidth);
printf("         reference            %s\n", $refGenomeFASTAFile);
printf("\n");

my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $finalDir = "$outputDir/final";
mkpath($finalDir);
my $statsDir = "$outputDir/stats";
mkpath($statsDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $slurmScriptsDir = "$outputDir/slurm_scripts/$slurmScriptsSubDir";
mkpath($slurmScriptsDir);

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

open(SQ,"$refGenomeFASTAFile.fai") || die "Cannot open $refGenomeFASTAFile.fai\n";
while (<SQ>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $len, $offset, $others) = split('\t');

        last if ($chrom=~/^GL/);
            
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
my $inputVCFFiles;
my $outputVCFFile;

########
#Calling
########

#**************
#log start time
#**************
$tgt = "$logDir/start.calling.OK";
$dep = "";
@cmd = ("date | awk '{print \"freebayes variant calling pipeline\\n\\nstart calling: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

for my $i (0 .. $#intervals)
{
    #call
    $outputVCFFile = "$auxDir/$intervalNames[$i].genotypes.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "";
    @cmd = ("$freebayes -f $refGenomeFASTAFile -L $bamListFile -r $intervals[$i] -v $outputVCFFile"),
    makeJob($partition, $tgt, $dep, @cmd);

    #inject contigs
    $inputVCFFile = "$auxDir/$intervalNames[$i].genotypes.vcf.gz";
    $tgt = "$inputVCFFile.hdr.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$injectContigs -v $inputVCFFile -c $contigsFile");
    makeJob($partition, $tgt, $dep, @cmd);
}

#************
#log end time
#************
$tgt = "$logDir/end.calling.OK";
$dep = join(" ", map {"$auxDir/$_.genotypes.vcf.gz.hdr.OK"} @intervalNames);
@cmd = ("date | awk '{print \"end: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

################################
#Concatenate and drop duplicates
################################

#**************
#log start time
#**************
$tgt = "$logDir/start.concatenatingOK";
$dep = "$logDir/end.calling.OK";
@cmd = ("date | awk '{print \"start concatenating and normalizing: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

for my $chrom (@CHROM)
{
    $inputVCFFiles = join(" ", map {"$auxDir/$_.genotypes.vcf.gz"} @{$intervalsByChrom{$chrom}});
    $outputVCFFile = "$finalDir/$chrom.genotypes.vcf.gz";
    $tgt = "$outputVCFFile.OK";
    $dep = "$logDir/end.calling.OK";
    @cmd = ("$vt cat $inputVCFFiles -o + -w 1000 | $vt uniq + -o $outputVCFFile 2> $statsDir/$chrom.uniq.log");
    makeJob($partition, $tgt, $dep, @cmd);

    $inputVCFFile = "$finalDir/$chrom.genotypes.vcf.gz";
    $tgt = "$inputVCFFile.tbi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeJob($partition, $tgt, $dep, @cmd);
}

$inputVCFFiles = join(" ", map {"$finalDir/$_.genotypes.vcf.gz"} (1..22,"X","Y","MT"));
$outputVCFFile = "$finalDir/all.genotypes.vcf.gz";
$tgt = "$outputVCFFile.OK";
$dep = join(" ", map {"$finalDir/$_.genotypes.vcf.gz.OK"} (1..22,"X","Y","MT"));
@cmd = ("$vt cat $inputVCFFiles -o $outputVCFFile");
makeJob($partition, $tgt, $dep, @cmd);

$inputVCFFile = "$finalDir/all.genotypes.vcf.gz";
$tgt = "$inputVCFFile.tbi.OK";
$dep = "$inputVCFFile.OK";
@cmd = ("$vt index $inputVCFFile");
makeJob($partition, $tgt, $dep, @cmd);

$inputVCFFile = "$finalDir/all.genotypes.vcf.gz";
$outputVCFFile = "$finalDir/all.sites.vcf.gz";
$tgt = "$outputVCFFile.OK";
$dep = "$inputVCFFile.OK";
@cmd = ("$vt view -s  $inputVCFFile -o $outputVCFFile");
makeJob($partition, $tgt, $dep, @cmd);

$inputVCFFile = "$finalDir/all.sites.vcf.gz";
$tgt = "$inputVCFFile.tbi.OK";
$dep = "$inputVCFFile.OK";
@cmd = ("$vt index $inputVCFFile");
makeJob($partition, $tgt, $dep, @cmd);

#************
#log end time
#************
$tgt = "$logDir/end.concatenating.OK";
$dep = "$finalDir/all.sites.vcf.gz.tbi.OK";
@cmd = ("date | awk '{print \"end concatenating and normalizing: \"\$\$0}' >> $logFile");
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
push(@cmds, "\t-rm -rf $outputDir/*.* $auxDir/*.* $auxDir/*/*.* $finalDir/*.* $statsDir/* $logDir/* $outputDir/intervals/*.*");

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