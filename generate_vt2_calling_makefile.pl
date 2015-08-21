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
my $vt = "/net/fantasia/home/atks/programs/vt/vt";
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
my $processBySample = ($intervalWidth==0);

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
        mkpath("$auxDir/$sampleID");
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

#############
#1. Discovery
#############

if ($processBySample)
{
    for my $sampleID (@SAMPLE)
    {
        $outputVCFFile = "$auxDir/$sampleID/all.genotypes.bcf";
        $tgt = "$outputVCFFile.OK";
        $dep = "";
        @cmd = ("$samtools view -h $BAMFILE{$sampleID} 20 -u | $bam clipoverlap --in -.ubam --out -.ubam | $vt discover2 -z -q 20 -b + -r $refGenomeFASTAFile -s $sampleID -o $outputVCFFile 2> $auxDir/$sampleID/all.discover2.log");
        #@cmd = ("$vt discover2 -z -q 20 -b $BAMFILE{$sampleID} -r $refGenomeFASTAFile -s $sampleID -I $intervalFiles[$i] -o $outputVCFFile 2> $auxDir/$sampleID/$intervals[$i].discover2.log");
        makeJob($partition, $tgt, $dep, @cmd);
    }    
}
else
{
    #log start time
    my $logFile = "$logDir/run.log";
    $tgt = "$logFile.start.OK";
    $dep = "";
    @cmd = "date | awk '{print \"vt calling pipeline\\n\\nstart: \"\$\$0}' > $logFile";
    makeLocalStep($tgt, $dep, @cmd);

    my @intervalSampleDiscoveryVCFFilesOK = ();
    
    #mine variants from aligned reads
    for my $i (0 .. $#intervals)
    {
        my $intervalVCFFilesOK = "";
        for my $sampleID (@SAMPLE)
        {
            $outputVCFFile = "$auxDir/$sampleID/$intervals[$i].genotypes.bcf";
            $tgt = "$outputVCFFile.OK";
            $dep = "";
            @cmd = ("$samtools view -h $BAMFILE{$sampleID} 20 -u | $bam clipoverlap --in -.ubam --out -.ubam | $vt discover2 -z -q 20 -b + -r $refGenomeFASTAFile -s $sampleID -I $intervalFiles[$i] -o $outputVCFFile 2> $auxDir/$sampleID/$intervals[$i].discover2.log");
            #@cmd = ("$vt discover2 -z -q 20 -b $BAMFILE{$sampleID} -r $refGenomeFASTAFile -s $sampleID -I $intervalFiles[$i] -o $outputVCFFile 2> $auxDir/$sampleID/$intervals[$i].discover2.log");
            makeJob($partition, $tgt, $dep, @cmd);
            
            $intervalVCFFilesOK .= " $outputVCFFile.OK";
        }
        
        push(@intervalSampleDiscoveryVCFFilesOK, $intervalVCFFilesOK);
    }
    
    #merge variants
    for my $i (0 .. $#intervals)
    {
        my $vcfFileList = "$auxDir/$intervals[$i]_vcf_file.list";
        open(IN, ">$vcfFileList");
        my @files = map {"$auxDir/$_/$intervals[$i].genotypes.bcf"} @SAMPLE;
        print IN join("\n", @files); 
        close(IN);
        
        $outputVCFFile = "$auxDir/$intervals[$i].sites.bcf";
        $tgt = "$outputVCFFile.OK";
        $dep = $intervalSampleDiscoveryVCFFilesOK[$i];
        @cmd = ("$vt merge_candidate_variants -L $vcfFileList -o $outputVCFFile 2> $auxDir/$intervals[$i].merge_candidate_variants.log");
        makeJob($partition, $tgt, $dep, @cmd);
    }
    
    #concatenate variants by chromosome
    my $chromVCFFilesOK = "";
    for my $chrom (@CHROM)
    {
        my @intervals =  @{$intervalsByChrom{$chrom}};
        my $vcfFileList = "$auxDir/$chrom" . "_vcf_file.list";
        open(IN, ">$vcfFileList");
        my @files = map {"$auxDir/$_.sites.bcf"} @intervals;
        print IN join("\n", @files); 
        close(IN);
                
        $outputVCFFile = "$finalDir/$chrom.sites.bcf";
        $tgt = "$outputVCFFile.OK";
        my @filesOK = map {"$auxDir/$_.sites.bcf.OK"} @intervals;
        $dep = join(" ", @filesOK);
        @cmd = ("$vt cat -L $vcfFileList -o $outputVCFFile");
        makeJob($partition, $tgt, $dep, @cmd);  
        
        $chromVCFFilesOK .= " $outputVCFFile.OK";    
    }    
    
    #log end time
    $tgt = "$logDir/end.discovery.OK";
    $dep = "$chromVCFFilesOK";
    @cmd = ("date | awk '{print \"end discovery: \"\$\$0}' >> $logFile");
    makeJob("local", $tgt, $dep, @cmd);
}



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