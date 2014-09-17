#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_vt_v0_calling_makefile

=head1 SYNOPSIS

 generate_vt_v0_calling_makefile [options]

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
my $vtDir;
my $clusterDir;
my $ext;
my $makeFile;
my $cluster;
my $sleep;
my $sampleFile;
my $sequenceLengthFile;
my $intervalWidth;
my $variantTypes;
my $refGenomeFASTAFile;
my $dbsnpVCFFile;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'o:s'=>\$outputDir,
                'b:s'=>\$vtDir,
                't:s'=>\$clusterDir,
                'e:s'=>\$ext,
                'm:s'=>\$makeFile,
                'c:s'=>\$cluster,
                'd:s'=>\$sleep,
                's:s'=>\$sampleFile,
                'l:s'=>\$sequenceLengthFile,
                'i:s'=>\$intervalWidth,
                'v:s'=>\$variantTypes,
                'r:s'=>\$refGenomeFASTAFile,
                'u:s'=>\$dbsnpVCFFile
                )
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

printf("generate_vt_v0_calling_makefile.pl\n");
printf("\n");
printf("options: output dir      %s\n", $outputDir);
printf("         vt path         %s\n", $vtDir);
printf("         cluster path    %s\n", $clusterDir);
printf("         ext             %s\n", $ext);
printf("         make file       %s\n", $makeFile);
printf("         cluster         %s\n", $cluster);
printf("         sleep           %s\n", $sleep);
printf("         sample file     %s\n", $sampleFile);
printf("         sequence length file %s\n", $sequenceLengthFile);
printf("         variant types   %s\n", $variantTypes);
#printf("         dbsnp           %s\n", $dbsnpVCFFile);
printf("         reference       %s\n", $refGenomeFASTAFile);
printf("\n");

my $indexExt = $ext eq "vcf.gz" ? "tbi" : "csi";
my $logDir = "$outputDir/log";
mkpath($logDir);
my $vcfOutputDir = "$outputDir/vcf";
mkpath($vcfOutputDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $finalDir = "$outputDir/final";
mkpath($finalDir);
my $statsDir = "$outputDir/stats";
mkpath($statsDir);

#read file locations and name of samples
my %SAMPLE = ();
my @sample = ();
open(INDEX,"$sampleFile") || die "Cannot open $sampleFile\n";
while (<INDEX>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $bamPath) = split(/\s+/);
        $SAMPLE{$sampleID} = $bamPath;
        push(@sample, $sampleID);
    }
}
close(INDEX);

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;
my $inputVCFFile;
my $outputVCFFile;

###################
#Generate intervals
###################
my %intervalsByChrom = ();
my %intervalsByChromOK = ();
my @intervals = ();
my @intervalNames = ();
my @intervalFiles = ();

open(SQ,"$sequenceLengthFile") || die "Cannot open $sequenceLengthFile\n";
while (<SQ>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $len) = split(/\s+/);

        print "processing $chrom\t$len ";

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

            push(@{$intervalsByChrom{$chrom}}, "$vcfOutputDir/$interval.vcf");
            push(@{$intervalsByChromOK{$chrom}}, "$vcfOutputDir/$interval.vcf.OK");
            push(@intervals, $interval);
            push(@intervalNames, $intervalName);
            push(@intervalFiles, $file);

            $count++;
        }

        print "added $count intervals\n";
    }
}
close(SQ);

#############
#1. Discovery
#############

#**************
#log start time
#**************
my $logFile = "$outputDir/run.log";
$tgt = "$logDir/start.discovery.OK";
$dep = "";
@cmd = "date | awk '{print \"vt calling pipeline\\n\\nstart discovery: \"\$\$0}' > $logFile";
makeLocalStep($tgt, $dep, @cmd);

my %GVCFFilesByInterval = ();
my $GVCFFiles = "";
my $GVCFFilesOK = "";

if ($intervalWidth!=0)
{
    for my $sampleID (@sample)
    {
        my $sampleDir = "$vcfOutputDir/$sampleID";
        mkpath("$sampleDir");
        
        my $sampleCandidateIntervalVCFFiles = "";
        my $sampleCandidateIntervalVCFFilesOK = "";
        for my $i (0..$#intervals)
        {
            $outputVCFFile = "$sampleDir/$intervalNames[$i].$ext";
            $tgt = "$outputVCFFile.OK";
            $dep = "";
            @cmd = "$vt discover -b $SAMPLE{$sampleID} -o + -v $variantTypes -r $refGenomeFASTAFile -s $sampleID -i $intervals[$i] 2> $sampleDir/discover.log | $vt normalize + -r $refGenomeFASTAFile -o $outputVCFFile 2> $sampleDir/normalize.log";
            makeStep($tgt, $dep, @cmd);

            $sampleCandidateIntervalVCFFiles .= " $outputVCFFile";
            $sampleCandidateIntervalVCFFilesOK .= " $outputVCFFile.OK";
        }

        #write candidate VCF files into a list
        my $sampleCandidateIntervalVCFFileList = "$vcfOutputDir/$sampleID/interval_vcf_file_list.txt";
        open(IN, ">$sampleCandidateIntervalVCFFileList") ||die "Cannot open $sampleCandidateIntervalVCFFileList\n";
        my $temp = $sampleCandidateIntervalVCFFiles;
        $temp =~ s/ /\n/g;
        $temp =~ s/^\s+//;
        print IN "$temp\n";
        close(IN);

        #concatenate
        $outputVCFFile = "$vcfOutputDir/$sampleID/$sampleID.sites.$ext";
        $tgt = "$outputVCFFile.OK";
        $dep = "$sampleCandidateIntervalVCFFilesOK";
        @cmd = ("$vt concat -L $sampleCandidateIntervalVCFFileList -o + | $vt view + -w 1000 -o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);
        
        $inputVCFFile = "$vcfOutputDir/$sampleID/$sampleID.sites.$ext";
        $tgt = "$inputVCFFile.$indexExt.OK";
        $dep = "$inputVCFFile.OK";
        @cmd = ("$vt index $inputVCFFile");
        makeStep($tgt, $dep, @cmd);
    }

    #write candidate VCF files into a list
    my $sampleVCFFileList = "$vcfOutputDir/sample_vcf_file_list.txt";
    my @sampleVCFFiles = map {"$vcfOutputDir/$_/$_.sites.$ext"} @sample;
    open(IN, ">$sampleVCFFileList") ||die "Cannot open  $sampleVCFFileList\n";
    my $temp = join("\n", @sampleVCFFiles);
    print IN "$temp\n";
    close(IN);

    #merging and filtering of initial sites using likelhood ratio cut off of 2.
    my $outputVCFFile = "$finalDir/all.sites.$ext";
    my @sampleVCFFileIndicesOK = map {"$vcfOutputDir/$_/$_.sites.$ext.$indexExt.OK"} @sample;
    my $sampleVCFFileIndicesOK = join(" ", @sampleVCFFileIndicesOK);
    $tgt = "$outputVCFFile.OK";
    $dep =  $sampleVCFFileIndicesOK;
    @cmd = ("$vt merge_candidate_variants -L $sampleVCFFileList -o $outputVCFFile 2> $outputVCFFile.log");
    makeStep($tgt, $dep, @cmd);

    #index candidate sites
    $tgt = "$outputVCFFile.$indexExt.OK";
    $dep = "$outputVCFFile.OK";
    @cmd = ("$vt index $outputVCFFile 2> $outputVCFFile.$indexExt.log");
    makeStep($tgt, $dep, @cmd);
}
else
{
#    for my $sampleID (@sample)
#    {
#        my $outputVCFFile = "$vcfOutputDir/$sampleID/$sampleID.vcf";
#        $tgt = "$outputVCFFile.OK";
#        $dep = "";
#        @cmd = "$vt discover -b $SAMPLE{$sampleID} -o + -v $variantTypes -r $refGenomeFASTAFile -s $sampleID $intervalsOption 2> $sampleDir/discover.log | $vt normalize + -r $refGenomeFASTAFile -o + 2> $sampleDir/normalize.log | $vt mergedups + -o $sampleVCFFile  2> $sampleDir/mergedups.log";
#        makeStep($tgt, $dep, @cmd);
#
#        $GVCFFiles .= " $outputVCFFile";
#        $GVCFFilesOK .= " $outputVCFFile.OK";
#    }
#
#    #write candidate VCF files into a list
#    my $candidateVCFFileList = "$auxDir/candidate_vcf_files.txt";
#    open(IN, ">$candidateVCFFileList") ||die "Cannot open $candidateVCFFileList\n";
#    my $temp = $candidateSitesVCFFiles;
#    $temp =~ s/ /\n/g;
#    $temp =~ s/^\s+//;
#    print IN "$temp\n";
#    close(IN);
#
#    #merging and filtering of initial sites using likelhood ratio cut off of 2.
#    my $candidateVCFFile = "$auxDir/all.sites.$ext";
#    $tgt = "$candidateVCFFile.OK";
#    $dep =  $candidateSitesVCFOKFiles;
#    @cmd = ("$vt merge_candidate_variants -L $candidateVCFFileList -o $candidateVCFFile 2> $candidateVCFFile.log");
#    makeStep($tgt, $dep, @cmd);
#
#    #index candidate sites
#    $tgt = "$candidateVCFFile.$indexExt.OK";
#    $dep = "$candidateVCFFile.OK";
#    @cmd = ("$vt index $candidateVCFFile 2> $candidateVCFFile.$indexExt.log");
#    makeStep($tgt, $dep, @cmd);
}

###############
##2. Genotyping
###############
#
##construct probes for candidate sites
#my $probesVCFFile = "$auxDir/probes.sites.$ext";
#$tgt = "$probesVCFFile.OK";
#$dep = "$candidateVCFFile.$indexExt.OK";
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

######################
#3. Merge and Annotate
######################

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
#@cmd = ("$vt paste -L $mergedVCFFileList -o + | $vt compute_features + -o + 2> $finalDir/compute_features.log | $vt annotate_dbsnp_rsid - -d $dbsnpVCFFile  2>$finalDir/annotate_dbsnp_rsid.log | $vt remove_overlap - -o $mergedVCFFile 2> $finalDir/remove_overlap.log");
#makeStep($tgt, $dep, @cmd);
#
##index
#$tgt = "$mergedVCFFile.$indexExt.OK";
#$dep = "$mergedVCFFile.OK";
#@cmd = ("$vt index $mergedVCFFile");
#makeStep($tgt, $dep, @cmd);
#
##generate site file
#my $mergedVCFSitesFile = "$finalDir/all.sites.$ext";
#$tgt = "$mergedVCFSitesFile.OK";
#$dep = "$mergedVCFFile.OK";
#@cmd = ("$vt view -s $mergedVCFFile -o $mergedVCFSitesFile");
#makeStep($tgt, $dep, @cmd);
#
##index
#$tgt = "$mergedVCFSitesFile.$indexExt.OK";
#$dep = "$mergedVCFSitesFile.OK";
#@cmd = ("$vt index $mergedVCFSitesFile");
#makeStep($tgt, $dep, @cmd);

#############
#log end time
#############
#$tgt = "$logFile.end.OK";
#$dep = "$probesVCFFile.$indexExt.OK";
#@cmd = ("\tdate | awk '{print \"end: \"\$\$0}' >> $logFile");
#makeLocalStep($tgt, $dep, @cmd);

####################
#Write out make file
####################

open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/*.* $vcfOutputDir/*.* $vcfOutputDir/*/*.* $finalDir/*.* $statsDir/* $logDir/* $outputDir/intervals/*.*");

for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
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