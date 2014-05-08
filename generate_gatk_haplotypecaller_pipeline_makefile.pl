#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_gatk_haplotypecaller_pipeline_makefile

=head1 SYNOPSIS

 generate_gatk_haplotypecaller_calling_pipeline_makefile [options]

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
                'j:s'=>\$jvmMemory
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
#you can set the  maximum memory here to be whatever you want
my $gatk = "/net/fantasia/home/atks/programs/jdk1.7.0_25/bin/java -jar -Xmx$jvmMemory /net/fantasia/home/atks/programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar";
my $gatk64g = "/net/fantasia/home/atks/programs/jdk1.7.0_25/bin/java -jar -Xmx64g /net/fantasia/home/atks/programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar";
my $vt = "$vtDir/vt";

printf("generate_gatk_ug_calling_makefile.pl\n");
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
printf("\n");

mkpath($outputDir);

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
        my ($sampleID, $bamPath) = split(/\s+/, $_);
        $SAMPLE{$sampleID} = $bamPath;
        push(@sample, $sampleID);
        $bamFiles .= "$bamPath\n";
    }
}
close(SA);

my $bamListFile = "$outputDir/bam.list";
open(OUT,">$bamListFile") || die "Cannot open $bamListFile\n";
print OUT $bamFiles;
close(OUT);

print "read in " . scalar(keys(%SAMPLE)) . " samples\n";
my $vcfOutDir = "$outputDir/vcf";
mkpath($vcfOutDir);
my $finalVCFOutDir = "$outputDir/final";
mkpath($finalVCFOutDir);

###################
#Generate intervals
###################
my %intervalsByChrom = ();
my @intervals = ();
my @CHROM = ();

my $writeIntervals = 1;

if ($intervalWidth!=0)
{
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
                else
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
    
                push(@{$intervalsByChrom{$chrom}}, "$interval");
                push(@intervals, $interval);
    
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
}

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;

###############
#log start time
###############
my $logFile = "$outputDir/run.log";
$tgt = "$logFile.start.OK";
$dep = "";
@cmd = ("date | awk '{print \"gatk unified genotyper calling pipeline\\n\\nstart: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

######################
#Per sample processing
######################
my %GVCFFilesByInterval = ();
my $GVCFFiles = "";
my $GVCFFilesOK = "";

if ($intervalWidth!=0)
{
    for my $interval (@intervals)
    {
        $GVCFFilesByInterval{$interval} = " ";
        
        for my $sampleID (@sample)
        {
            my $outputVCFFile = "$vcfOutDir/$sampleID.$interval.vcf";
            $tgt = "$outputVCFFile.OK";
            $dep = "";
            @cmd = ("$gatk -T HaplotypeCaller -R $refGenomeFASTAFile -I $SAMPLE{$sampleID} " .
                    "--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 " .
                    "-L $outputDir/intervals/$interval.interval_list " .
                    "-o $outputVCFFile");
            makeStep($tgt, $dep, @cmd);
        
            $GVCFFilesByInterval{$interval} .= " --variant $outputVCFFile";
            $GVCFFilesOK .= " $outputVCFFile.OK";
        }
    }
}
else
{
    for my $sampleID (@sample)
    {
        my $outputVCFFile = "$vcfOutDir/$sampleID.vcf";
        $tgt = "$outputVCFFile.OK";
        $dep = "";
        @cmd = ("$gatk -T HaplotypeCaller -R $refGenomeFASTAFile -I $SAMPLE{$sampleID} " .
                "--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 " .
                "-o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);
    
        $GVCFFiles .= " --variant $outputVCFFile";
        $GVCFFilesOK .= " $outputVCFFile.OK";
    }
}


#############
#log end time
#############
$tgt = "$logFile.end.gvcfs_creation.OK";
$dep = "$GVCFFilesOK";
@cmd = ("date | awk '{print \"end gvcfs creation: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

############
#Combine VCF
############
$tgt = "$logFile.start.merge.OK";
$dep = "$logFile.end.gvcfs_creation.OK";
@cmd = ("date | awk '{print \"start merge: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

my $mergedVCFFilesOK = "";

if ($intervalWidth!=0)
{
    for my $interval (@intervals)
    {
        my $outputVCFFile = "$vcfOutDir/$interval.vcf";
        $tgt = "$outputVCFFile.OK";
        $dep = "$logFile.start.merge.OK";;
        @cmd = ("$gatk -T CombineGVCFs -R $refGenomeFASTAFile $GVCFFilesByInterval{$interval} -o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);
        
        $mergedVCFFilesOK .= " $outputVCFFile.OK";
    }
}
else
{
    my $outputVCFFile = "$vcfOutDir/all.vcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "$logFile.start.merge.OK";;
    @cmd = ("$gatk -T CombineGVCFs -R $refGenomeFASTAFile $GVCFFiles -o $outputVCFFile");
    makeStep($tgt, $dep, @cmd);
    
    $mergedVCFFilesOK = "$outputVCFFile.OK";
}

$tgt = "$logFile.end.merge.OK";
$dep = "$mergedVCFFilesOK";
@cmd = ("date | awk '{print \"end merge: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

#########
#Genotype
#########
$tgt = "$logFile.start.genotype.OK";
$dep = "$logFile.end.merge.OK";
@cmd = ("date | awk '{print \"start genotype: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

my $genotypeVCFFilesOK = "";

if ($intervalWidth!=0)
{
    for my $interval (@intervals)
    {
        my $inputVCFFile = "$vcfOutDir/$interval.vcf";
        my $outputVCFFile = "$vcfOutDir/$interval.genotypes.vcf";
        $tgt = "$outputVCFFile.OK";
        $dep = "";
        @cmd = ("$gatk -T GenotypeGVCFs -R $refGenomeFASTAFile --variant $inputVCFFile -o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);
        
        $genotypeVCFFilesOK .= " $outputVCFFile.OK";
    }
}
else
{
    my $inputVCFFile = "$vcfOutDir/all.vcf";
    my $outputVCFFile = "$finalVCFOutDir/all.genotypes.vcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "";
    @cmd = ("$gatk -T CombineGVCFs -R $refGenomeFASTAFile --variant $inputVCFFile -o $outputVCFFile");
    makeStep($tgt, $dep, @cmd);
    
    $genotypeVCFFilesOK = "$outputVCFFile.OK";
}

$tgt = "$logFile.end.genotype.OK";
$dep = "$genotypeVCFFilesOK";
@cmd = ("date | awk '{print \"end genotype: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

############################################
##Concatenate, normalize and drop duplicates
############################################
if ($intervalWidth!=0)
{
    $tgt = "$logFile.start.concat.normalize.OK";
    $dep = "$logFile.end.genotype.OK";
    @cmd = ("date | awk '{print \"start concat and normalize: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
    
    my $chromGenotypeVCFFilesOK = "";

    for my $chrom (@CHROM)
    {
        my $inputGenotypeVCFFiles = "";
        my $inputGenotypeVCFFilesOK = "";
        $chromGenotypeVCFFilesOK .= " $finalVCFOutDir/$chrom.genotypes.vcf.gz.OK";
        
        for my $interval (@{$intervalsByChrom{$chrom}})
        {
            $inputGenotypeVCFFiles .= " $vcfOutDir/$interval.genotypes.vcf";
            $inputGenotypeVCFFilesOK .= " $vcfOutDir/$interval.genotypes.vcf.OK";
        }
    
        my $outputVCFFile = "$finalVCFOutDir/$chrom.genotypes.vcf.gz";
        $tgt = "$outputVCFFile.OK";
        $dep = "$inputGenotypeVCFFilesOK";
        @cmd = ("$vt concat $inputGenotypeVCFFiles -o + | $vt normalize -r $refGenomeFASTAFile + -o + | $vt mergedups + -o $outputVCFFile ");
        makeStep($tgt, $dep, @cmd);
        
        $tgt = "$finalVCFOutDir/$chrom.genotypes.vcf.gz.tbi.OK";
        $dep = "$finalVCFOutDir/$chrom.genotypes.vcf.gz.OK";
        @cmd = ("$vt index $outputVCFFile");
        makeStep($tgt, $dep, @cmd);
        
        my $inputVCFFile = "$finalVCFOutDir/$chrom.genotypes.vcf.gz";
        $outputVCFFile = "$finalVCFOutDir/$chrom.sites.vcf.gz";
        $tgt = "$outputVCFFile.OK";
        $dep = "$inputVCFFile.OK";
        @cmd = ("$vt view -s $inputVCFFile -o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);
        
        $inputVCFFile = "$finalVCFOutDir/$chrom.sites.vcf.gz";
        $outputVCFFile = "$finalVCFOutDir/$chrom.sites.vcf.gz.tbi";
        $tgt = "$outputVCFFile.OK";
        $dep = "$inputVCFFile.OK";
        @cmd = ("$vt index $inputVCFFile");
        makeStep($tgt, $dep, @cmd);
        
        
    }

    $tgt = "$logFile.end.start.concat.normalize.OK";
    $dep = "$chromGenotypeVCFFilesOK";
    @cmd = ("date | awk '{print \"end concat and normalize: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}



#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

######
#clean
######
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm $outputDir/*.OK $vcfOutDir");

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