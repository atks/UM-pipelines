#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_trf_calling_pipeline_makefile

=head1 SYNOPSIS

 generate_trf_calling_pipeline_makefile [options]

  -f     fastq file list giving the location of each sample
         column 1: fastq file containing
         column 2: path of bam file
  -r     reference genome file
  -l     sequence length file
  -w     interval width
  -o     output directory
  -m     make file name

=head1 DESCRIPTION

This script generates the make file to discover STRs with Tandem Repeat Finder.

=cut

my $help;

my $outputDir;
my $vtDir;
my $clusterDir;
my $makeFile;
my $cluster;
my $sleep;
my $sequenceLengthFile;
my $refGenomeFASTAFile;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o:s'=>\$outputDir,
                'b:s'=>\$vtDir,
                't:s'=>\$clusterDir,
                'm:s'=>\$makeFile,
                'c:s'=>\$cluster,
                'd:s'=>\$sleep,
                'l:s'=>\$sequenceLengthFile,
                'r:s'=>\$refGenomeFASTAFile
                )
  || !defined($makeFile)
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
my $trf = "/net/fantasia/home/atks/dev/vt/comparisons/programs/trf-4.07b/trf";
my $trfScriptsDir = "/net/fantasia/home/atks/dev/vt/comparisons/programs/trf-4.07b/scripts/";
my $vt = "$vtDir/vt";

printf("generate_trf_calling_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         vt path              %s\n", $vt);
printf("         cluster path         %s\n", $clusterDir);
printf("         make file            %s\n", $makeFile);
printf("         cluster              %s\n", $cluster);
printf("         sleep                %s\n", $sleep);
printf("         sequence length file %s\n", $sequenceLengthFile);
printf("         reference            %s\n", $refGenomeFASTAFile);
printf("\n");

my $fastaDir = "$outputDir/fasta";
mkpath($fastaDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $finalDir = "$outputDir/final";
mkpath($finalDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $logFile = "$outputDir/run.log";

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;

my @CHROMS = ();

open(SQ,"$sequenceLengthFile") || die "Cannot open $sequenceLengthFile\n";
while (<SQ>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $len) = split('\t', $_);

        push(@CHROMS, $chrom);
    }
}
close(SQ);

##########
#Discovery
##########

#****************************
#log start time for discovery
#****************************
$tgt = "$logDir/start.discovery.OK";
$dep = "";
@cmd = ("date | awk '{print \"Tandem Repeat Finder STR calling pipeline\\n\\nstart calling: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);


#FASTA File preparation
$tgt = "$fastaDir/all.OK";
$dep = "$refGenomeFASTAFile";
@cmd = ("$trfScriptsDir/split_fasta_file $refGenomeFASTAFile -o $fastaDir");
makeLocalStep($tgt, $dep, @cmd);

my $matchWeight = 2; 
my $mismatch = 7;
my $delta = 7;
my $pm = 80;
my $pi = 10;
my $minScore = 14;
my $maxPeriod = 500;

#trf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]
#Where: (all weights, penalties, and scores are positive)
#  File = sequences input file
#  Match  = matching weight
#  Mismatch  = mismatching penalty
#  Delta = indel penalty
#  PM = match probability (whole number)
#  PI = indel probability (whole number)
#  Minscore = minimum alignment score to report
#  MaxPeriod = maximum period size to report
#  [options] = one or more of the following :
#               -m    masked sequence file
#               -f    flanking sequence
#               -d    data file
#               -h    suppress html output
#               -r    no redundancy elimination
#               -ngs  more compact .dat output on multisequence files, returns 0 on success. You may pipe input in with this option using - for file name. Short 50 flanks are appended to .dat output. See more information on TRF Unix Help web page.

#Call STRs with TRF
my $chromDatFiles = "";
my $chromDatFilesOK = "";
for my $chrom (@CHROMS)
{
    my $inputFASTAFile = "$fastaDir/$chrom.fa";
    my $outputDir = "$auxDir";
    
    my $chromDatFile = "$auxDir/$chrom.fa.$matchWeight.$mismatch.$delta.$pm.$pi.$minScore.$maxPeriod.dat";
    $tgt = "$chromDatFile.OK";
    $dep = "$fastaDir/all.OK";
    @cmd = ("cd $outputDir; $trf $inputFASTAFile $matchWeight $mismatch $delta $pm $pi $minScore $maxPeriod -f -d -h 2> $auxDir/$chrom.log > /dev/null; echo done > /dev/null");
    makeStep($tgt, $dep, @cmd);   

    $chromDatFiles .= " $chromDatFile";
    $chromDatFilesOK .= " $chromDatFile.OK";
}

#Concatenate into VCF file
my $outputVCFFile = "$finalDir/all.sites.vcf.gz";
$tgt = "$outputVCFFile.OK";
$dep = $chromDatFilesOK;
@cmd = ("$trfScriptsDir/cat_trf_dat_files $chromDatFiles -o $outputVCFFile");
makeStep($tgt, $dep, @cmd);   

#**************************
#log end time for discovery
#**************************
$tgt = "$logDir/end.discovery.OK";
$dep = "$outputVCFFile.OK";
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
push(@cmds, "\t-rm -rf $outputDir/*.* $fastaDir/*.* $finalDir/*.* $logDir/* ");

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