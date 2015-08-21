#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_vt2_calling_makefile

=head1 SYNOPSIS

 generate_vt2_calling_makefile [options]

  -s     sample file list giving the location of each sample
         column 1: sample name
         column 2: path of bam file
  -r     reference sequence fasta file
  -b     binaries directory : location of binaries required for this pipeline
  -o     output directory : location of all output files
  -m     output make file

 example:

=head1 DESCRIPTION

=cut

#option variables
my $help;

#
my $outputDir;
my $vtDir = "/net/fantasia/home/atks/dev/vt";
my $clusterDir = "/net/fantasia/home/atks/programs/cluster";
my $makeFile = "Makefile";
my $partition = "nomosix";
my $sleep = 0;
my $sampleFile = "";
my $intervalWidth = 20000000;
my $sequenceLengthFile =
my $refGenomeFASTAFile = "";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o:s'=>\$outputDir,
                'b:s'=>\$vtDir,
                'm:s'=>\$makeFile,
                'p:s'=>\$partition,
                's:s'=>\$sampleFile,
                'w:s'=>\$intervalWidth,
                'l:s'=>\$sequenceLengthFile,
                'r:s'=>\$refGenomeFASTAFile
                )
  || !defined($outputDir)
  || !defined($sampleFile)
  || !defined($sequenceLengthFile)
  || !defined($refGenomeFASTAFile)
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

#programs
my $vt = "$vtDir/vt";
my $samtools = "/net/fantasia/home/atks/programs/samtools/samtools";
my $bam = "/usr/cluster/bin/bam";

printf("generate_vt2_calling_makefile.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         vt path              %s\n", $vtDir);
printf("         make file            %s\n", $makeFile);
printf("         partition            %s\n", $partition);
printf("         sample file          %s\n", $sampleFile);
printf("         interval width       %s\n", $intervalWidth);
printf("         sequence length file %s\n", $sequenceLengthFile);
printf("         reference            %s\n", $refGenomeFASTAFile);
printf("\n");

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;
my $inputVCFFile;
my $outputVCFFile;

mkpath($outputDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $finalDir = "$outputDir/final";
mkpath($finalDir);

###########################
#Read samples and BAM paths
###########################
my %BAMFILE = ();
my @SAMPLE = ();
open(INDEX,"$sampleFile") || die "Cannot open $sampleFile\n";
while (<INDEX>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $bamPath) = split('\t', $_);
        $BAMFILE{$sampleID} = $bamPath;
        push(@SAMPLE, $sampleID);
    }
}
close(INDEX);

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

###############
#log start time
###############
my $logFile = "$logDir/run.log";
$tgt = "$logFile.start.OK";
$dep = "";
@cmd = "date | awk '{print \"vt calling pipeline\\n\\nstart: \"\$\$0}' > $logFile";
makeLocalStep($tgt, $dep, @cmd);

#############
#1. Discovery
#############
my $candidateSitesVCFFiles = "";
my $candidateSitesVCFOKFiles = "";

if ($intervalWidth!=0)
{
    my $intervalVCFFilesOK = "";
    for my $sampleID (@SAMPLE)
    {
        mkdir("$auxDir/$sampleID/");
        for my $i (0 .. $#intervals)
        {
            $outputVCFFile = "$auxDir/$sampleID/$intervals[$i].genotypes.bcf";
            $tgt = "$outputVCFFile.OK";
            $dep = "";
            #@cmd = ("$samtools view -h $BAMFILE{$sampleID} 20 -u | $bam clipoverlap --in -.ubam --out -.ubam | $vt discover2 -z -q 20 -b + -r $refGenomeFASTAFile -s $sampleID -I $intervalFiles[$i] -o $outputVCFFile 2> $auxDir/$sampleID/$intervals[$i].discover2.log");
            @cmd = ("$vt discover2 -z -q 20 -b $BAMFILE{$sampleID} -r $refGenomeFASTAFile -s $sampleID -I $intervalFiles[$i] -o $outputVCFFile 2> $auxDir/$sampleID/$intervals[$i].discover2.log");
            makeJob($partition, $tgt, $dep, @cmd);
        }
    }
    
    #************
    #log end time
    #************
    $tgt = "$logDir/end.calling.OK";
    $dep = "$intervalVCFFilesOK";
    @cmd = ("date | awk '{print \"end calling: \"\$\$0}' >> $logFile");
    makeJob("local", $tgt, $dep, @cmd);
}
else
{
#    $outputVCFFile = "$vcfOutDir/all.vcf.gz";
#    $tgt = "$outputVCFFile.OK";
#    $dep = "";
#    @cmd = ("$gatk -T UnifiedGenotyper -R $refGenomeFASTAFile -glm $variantType -I $bamListFile --genotyping_mode DISCOVERY -o $outputVCFFile --output_mode EMIT_VARIANTS_ONLY");
#    makeStep($tgt, $dep, @cmd);
#
#    #************
#    #log end time
#    #************
#    $tgt = "$logDir/end.calling.OK";
#    $dep = "$outputVCFFile.OK";
#    @cmd = ("date | awk '{print \"end calling: \"\$\$0}' >> $logFile");
#    makeLocalStep($tgt, $dep, @cmd);
}

if ($intervalWidth!=0)
{
    for my $i (0 .. $#intervals)
    {
        open(IN, ">$auxDir/$intervals[$i]_vcf_file.list");
        my @files = map {"$auxDir/$_/$intervals[$i].genotypes.bcf"} @SAMPLE;
        print IN join("\n", @files); 
        close(IN);
        
#        $outputVCFFile = "$auxDir/$sampleID/$intervals[$i].genotypes.bcf";
#        $tgt = "$outputVCFFile.OK";
#        $dep = "";
##            @cmd = ("$samtools view -h $BAMFILE{$sampleID} 20 | $bam clipoverlap --in - --out - | $vt discover2 -l -z -q 20 -b $BAMFILE{$sampleID} -r $refGenomeFASTAFile -s $sampleID -I $intervalFiles[$i] -o $outputVCFFile 2> $auxDir/$sampleID/$intervals[$i].discover2.log");
#        @cmd = ("$vt discover2 -l -q 20 -b $BAMFILE{$sampleID} -r $refGenomeFASTAFile -s $sampleID -I $intervalFiles[$i] -o $outputVCFFile 2> $auxDir/$sampleID/$intervals[$i].discover2.log");
#        makeJob($partition, $tgt, $dep, @cmd);
        
    }
    
}

##write candidate VCF files into a list
#my $candidateVCFFileList = "$auxDir/candidate_vcf_files.txt";
#open(IN, ">$candidateVCFFileList") ||die "Cannot open $candidateVCFFileList\n";
#my $temp = $candidateSitesVCFFiles;
#$temp =~ s/ /\n/g;
#$temp =~ s/^\s+//;
#print IN "$temp\n";
#close(IN);
#
##merging and filtering of initial sites using likelhood ratio cut off of 2.
#my $candidateVCFFile = "$auxDir/all.sites.$ext";
#$tgt = "$candidateVCFFile.OK";
#$dep =  $candidateSitesVCFOKFiles;
#@cmd = ("$vt merge_candidate_variants -L $candidateVCFFileList -o $candidateVCFFile 2> $candidateVCFFile.log");
#makeStep($tgt, $dep, @cmd);
#
##index candidate sites
#$tgt = "$candidateVCFFile.$indexExt.OK";
#$dep = "$candidateVCFFile.OK";
#@cmd = ("$vt index $candidateVCFFile 2> $candidateVCFFile.$indexExt.log");
#makeStep($tgt, $dep, @cmd);

###############
##2. Genotyping
###############

##construct probes for candidate sites
#my $probesVCFFile = "$auxDir/probes.sites.$ext";
#$tgt = "$probesVCFFile.OK";
#$dep = "$candidateVCFFile.OK";
#@cmd = ("$vt construct_probes $candidateVCFFile -r $refGenomeFASTAFile -o $probesVCFFile 2> $auxDir/probes.log");
#makeStep($tgt, $dep, @cmd);
#
##index probes sites
#$tgt = "$probesVCFFile.$indexExt.OK";
#$dep = "$probesVCFFile.OK";
#@cmd = ("$vt index $probesVCFFile 2> $probesVCFFile.$indexExt.log");
#makeStep($tgt, $dep, @cmd);
#
##per sample discovery of sites
#my $candidateSitesVCFOutDir = "$outputDir/vcf";
#my $sampleGenotypesVCFFiles;
#my $sampleGenotypesVCFIndexOKFiles;
#mkdir($candidateSitesVCFOutDir);
#for my $sampleID (keys(%SAMPLE))
#{
#    mkdir("$candidateSitesVCFOutDir/$sampleID/");
#
#    my $sampleDir = "$candidateSitesVCFOutDir/$sampleID";
#    my $sampleVCFFile = "$sampleDir/$sampleID.genotypes.$ext";
#    my $sampleVCFFileIndex = "$sampleDir/$sampleID.genotypes.$ext.$indexExt";
#
#    #genotype
#    $tgt = "$sampleVCFFile.OK";
#    $dep = "$probesVCFFile.OK";
#    @cmd = ("$vt genotype -b $SAMPLE{$sampleID}{HC} $probesVCFFile -r $refGenomeFASTAFile -s $sampleID -o $sampleVCFFile  2> $sampleDir/genotype.log");
#    makeStep($tgt, $dep, @cmd);
#
#    #index
#    $tgt = "$sampleVCFFileIndex.OK";
#    $dep = "$sampleVCFFile.OK";
#    @cmd = ("$vt index $sampleVCFFile");
#    makeStep($tgt, $dep, @cmd);
#
#    $sampleGenotypesVCFFiles .= "$sampleVCFFile\n";
#    $sampleGenotypesVCFIndexOKFiles .= " $sampleVCFFileIndex.OK";
#}

#######################
##3. Merge and Annotate
#######################
#
##make merge list
#my $mergedVCFFileList = "$auxDir/merge.vcf.list.txt";
#chomp($sampleGenotypesVCFFiles);
#open(IN, ">$mergedVCFFileList") || die "Cannot open $mergedVCFFileList";
#print IN $sampleGenotypesVCFFiles;
#close(IN);
#
##merge
#my $mergedVCFFile = "$finalDir/all.genotypes.$ext";
#$tgt = "$mergedVCFFile.OK";
#$dep = "$sampleGenotypesVCFIndexOKFiles";
#@cmd = ("$vt paste -L $mergedVCFFileList -o + | $vt compute_features + -o + 2> $finalDir/compute_features.log | $vt remove_overlap + -o $mergedVCFFile 2> $finalDir/remove_overlap.log");
#makeStep($tgt, $dep, @cmd);
#
##index
#$tgt = "$mergedVCFFile.$indexExt.OK";
#$dep = "$mergedVCFFile.OK";
#@cmd = ("$vt index $mergedVCFFile");
#makeStep($tgt, $dep, @cmd);

##############
##log end time
##############
#$tgt = "$logFile.end.OK";
#$dep = "$probesVCFFile.$indexExt.OK";
#@cmd = ("\tdate | awk '{print \"end: \"\$\$0}' >> $logFile");
#makeLocalStep($tgt, $dep, @cmd);

####################
#Write out make file
####################

open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK ".PHONY: clean\n\n";
print MAK "all: @tgts\n\n";

#clean
$tgt = "clean";
$dep = "";
@cmd = ("-rm -rf $finalDir/*.OK $auxDir/*.OK $logDir/*.OK");
makePhonyJob($tgt, $dep, @cmd);

for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
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
        $cmd .= "\tsrun -p $partition " . $c . "\n";
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
        $cmd .= "\t" . $c . "\n";
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