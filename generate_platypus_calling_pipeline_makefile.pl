#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_platypus_calling_pipeline_makefile

=head1 SYNOPSIS

 generate_platypus_calling_pipeline_makefile [options]

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

my $outputDir = "";
my $vtDir = "";
my $clusterDir = "";
my $makeFile = "Makefile";
my $cluster = "main";
my $sleep = 0;
my $sampleFile = "";
my $intervals = "";
my $sequenceLengthFile = "";
my $intervalWidth = "";
my $refGenomeFASTAFile = "";
my $jvmMemory = "2g";
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
                'x'=>\$rawCopy
                )
  || !defined($outputDir)
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
my $platypus = "/net/fantasia/home/atks/dev/vt/comparisons/programs/python_2.7.3/python /net/fantasia/home/atks/dev/vt/comparisons/programs/Platypus_0.7.8/Platypus.py";
my $injectContigs = "/net/fantasia/home/atks/dev/vt/comparisons/programs/scripts/inject_contigs";
my $vt = "$vtDir/vt";

printf("generate_platypus_calling_pipeline_makefile.pl\n");
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
printf("         Raw Copy             %s\n", $rawCopy);
printf("\n");

my $vcfOutDir = "$outputDir/vcf";
mkpath($vcfOutDir);
my $finalVCFOutDir = "$outputDir/final";
mkpath($finalVCFOutDir);
my $statsDir = "$outputDir/stats";
mkpath($statsDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $logFile = "$outputDir/run.log";
my $rawCopyDir = "$outputDir/raw";
if ($rawCopy)
{
    mkpath($rawCopyDir);
}

########################################
#Read file locations and name of samples
########################################
my %SAMPLE = ();
my @sample = ();
open(SA,"$sampleFile") || die "Cannot open $sampleFile\n";
my $bamFiles = "";
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        chomp;
        my ($sampleID, $bamPath) = split(/\s+/);
        $SAMPLE{$sampleID} = $bamPath;
        push(@sample, $sampleID);
        $bamFiles .= ($bamFiles eq "" ? "" : ",") . "$bamPath";
    }
}
close(SA);

print "read in " . scalar(keys(%SAMPLE)) . " samples\n";

###################
#Generate intervals
###################
my %intervalsByChrom = ();
my @intervals = ();
my @intervalNames = ();
my @CHROM = ();

if ($intervalWidth!=0)
{
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
                
                $count++;
            }

            print "added $count intervals\n";
        }
    }
    close(SQ);
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

#**************************
#log start time for calling
#**************************
$tgt = "$logDir/start.calling.OK";
$dep = "";
@cmd = ("date | awk '{print \"platypus variant calling pipeline\\n\\nstart calling: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

my $intervalVCFFiles = "";
my $intervalVCFFileHdrsOK = "";
my $contigsFile = "/net/fantasia/home/atks/dev/vt/comparisons/na12878/contigs.txt";

if ($intervalWidth!=0)
{
    for my $i (0..$#intervals)
    {
        my $interval = $intervals[$i];
        my $intervalName = $intervalNames[$i];        
        $outputVCFFile = "$vcfOutDir/$intervalName.vcf";
        $tgt = "$outputVCFFile.OK";
        $dep = "";
        @cmd = ("$platypus callVariants --bamFiles=$bamFiles --refFile=$refGenomeFASTAFile --output=$outputVCFFile --regions=$interval");
        makeStep($tgt, $dep, @cmd);
        
        $tgt = "$outputVCFFile.hdr.OK";
        $dep = "$outputVCFFile.OK";
        @cmd = ("$injectContigs -v $outputVCFFile -c $contigsFile");
        makeStep($tgt, $dep, @cmd);
            
        $intervalVCFFiles .= " $outputVCFFile";
        $intervalVCFFileHdrsOK .= " $outputVCFFile.hdr.OK";
    }
    
    #************************
    #log end time for calling
    #************************
    $tgt = "$logDir/end.calling.OK";
    $dep = "$intervalVCFFileHdrsOK";
    @cmd = ("date | awk '{print \"end calling: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}
else
{
    $outputVCFFile = "$vcfOutDir/all.vcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "";
    @cmd = ("$platypus callVariants --bamFiles=$bamFiles --refFile=$refGenomeFASTAFile --output=$outputVCFFile");
    makeStep($tgt, $dep, @cmd);

    $tgt = "$outputVCFFile.hdr.OK";
    $dep = "$outputVCFFile.OK";
    @cmd = ("$injectContigs -v $outputVCFFile -c $contigsFile");
    makeStep($tgt, $dep, @cmd);
    
    #************************
    #log end time for calling
    #************************
    $tgt = "$logDir/end.calling.OK";
    $dep = "$outputVCFFile.hdr.OK";
    @cmd = ("date | awk '{print \"end calling: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}

###########################################
#Concatenate, normalize and drop duplicates
###########################################
if ($intervalWidth!=0)
{
    #*************************************************
    #log start time for concating and normalizing VCFs
    #*************************************************
    $tgt = "$logDir/start.concatenation.normalization.OK";
    $dep = "$logDir/end.calling.OK";
    @cmd = ("date | awk '{print \"start concatenating and normalizing: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);

    my $chromGenotypeVCFFilesOK = "";
    my $chromSiteVCFFiles = "";
    my $chromSiteVCFFilesOK = "";

    for my $chrom (@CHROM)
    {
        my $inputChromosomeIntervalVCFFiles = "";
        my $inputChromosomeIntervalVCFFileHdrsOK = "";
        $outputVCFFile = "$finalVCFOutDir/$chrom.genotypes.vcf.gz";
        $chromGenotypeVCFFilesOK .= " $outputVCFFile.OK";
        $chromSiteVCFFiles .= " $finalVCFOutDir/$chrom.sites.vcf.gz";
        $chromSiteVCFFilesOK .= " $finalVCFOutDir/$chrom.sites.vcf.gz.OK";

        my $chromVCFFileListFile = "$auxDir/$chrom.list";
        open(OUT, ">$chromVCFFileListFile") || die "Cannot open $chromVCFFileListFile";
        for my $interval (@{$intervalsByChrom{$chrom}})
        {
            print OUT "$vcfOutDir/$interval.vcf\n";
            $inputChromosomeIntervalVCFFiles .= " $vcfOutDir/$interval.vcf";
            $inputChromosomeIntervalVCFFileHdrsOK .= " $vcfOutDir/$interval.vcf.hdr.OK";
        }
        close(OUT);

        #genotypes VCFs
        $tgt = "$outputVCFFile.OK";
        $dep = "$logDir/end.calling.OK";
        @cmd = ("$vt cat -L $chromVCFFileListFile -o + | $vt normalize + -r $refGenomeFASTAFile - 2> $statsDir/$chrom.normalize.log | $vt uniq - -o $outputVCFFile 2> $statsDir/$chrom.uniq.log");
        makeStep($tgt, $dep, @cmd);

        $tgt = "$outputVCFFile.tbi.OK";
        $dep = "$outputVCFFile.OK";
        @cmd = ("$vt index $outputVCFFile");
        makeStep($tgt, $dep, @cmd);

        #sites VCFs
        $inputVCFFile = "$finalVCFOutDir/$chrom.genotypes.vcf.gz";
        $outputVCFFile = "$finalVCFOutDir/$chrom.sites.vcf.gz";
        $tgt = "$outputVCFFile.OK";
        $dep = "$inputVCFFile.OK";
        @cmd = ("$vt view -s $inputVCFFile -o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);

        $tgt = "$outputVCFFile.tbi.OK";
        $dep = "$outputVCFFile.OK";
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
    
    if ($rawCopy)
    {
        for my $chrom (@CHROM)
        {
            my $inputGenotypeVCFFiles = "";
            my $inputGenotypeVCFFilesOK = "";
            $chromGenotypeVCFFilesOK .= " $rawCopyDir/$chrom.genotypes.vcf.gz.OK";
    
            for my $interval (@{$intervalsByChrom{$chrom}})
            {
                $inputGenotypeVCFFiles .= " $vcfOutDir/$interval.vcf";
                $inputGenotypeVCFFilesOK .= " $vcfOutDir/$interval.vcf.OK";
            }
    
            $outputVCFFile = "$rawCopyDir/$chrom.genotypes.vcf.gz";
            $tgt = "$outputVCFFile.OK";
            $dep = "$inputGenotypeVCFFilesOK";
            @cmd = ("$vt cat $inputGenotypeVCFFiles -o + | $vt uniq + -o $outputVCFFile 2> $statsDir/$chrom.raw.uniq.log");
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
    
    #***************************************************
    #log end time for concatenating and normalizing VCFs
    #***************************************************
    $tgt = "$logDir/end.concatenation.normalization.OK";
    $dep = "$finalVCFOutDir/all.sites.vcf.gz.tbi.OK";
    @cmd = ("date | awk '{print \"end concat and normalize: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}
else
{
    #**********************************
    #log start time for normalizing VCF
    #**********************************
    $tgt = "$logDir/start.normalization.OK";
    $dep = "$logDir/end.genotyping.OK";
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
    print MAK "$tgts[$i]: $deps[$i]\n";
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