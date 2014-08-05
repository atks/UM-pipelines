#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_vt_calling_makefile

=head1 SYNOPSIS

 generate_vt_calling_makefile [options]

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
my $ext = "bcf";
my $makeFile = "Makefile";
my $cluster = "mini";
my $sleep = 0;
my $sampleFile = "";
my $intervals = "";
my $variantTypes = "indels";
my $refGenomeFASTAFile;
my $dbsnpVCFFile;

#programs



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
                'v:s'=>\$variantTypes,
                'i:s'=>\$intervals,
                'r:s'=>\$refGenomeFASTAFile,
                'u:s'=>\$dbsnpVCFFile
                )
  || !defined($outputDir)
  || !defined($makeFile)
  || !defined($sampleFile)
  || !defined($refGenomeFASTAFile)
  || !defined($dbsnpVCFFile)
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

my $indexExt = $ext eq "vcf.gz" ? "tbi" : "csi";
my $vt = "$vtDir/vt";

my $intervalsOption = ($intervals ne "") ? ((-e $intervals) ? "-I $intervals" : "-i $intervals") : "";

printf("generate_vt_calling_makefile.pl\n");
printf("\n");
printf("options: output dir      %s\n", $outputDir);
printf("         vt path         %s\n", $vtDir);
printf("         cluster path    %s\n", $clusterDir);
printf("         ext             %s\n", $ext);
printf("         make file       %s\n", $makeFile);
printf("         cluster         %s\n", $cluster);
printf("         sleep           %s\n", $sleep);
printf("         sample file     %s\n", $sampleFile);
printf("         variant types   %s\n", $variantTypes);
printf("         intervals       %s\n", $intervalsOption eq "" ? "none" : $intervalsOption);
printf("         dbsnp           %s\n", $dbsnpVCFFile);
printf("         reference       %s\n", $refGenomeFASTAFile);
printf("\n");

#read file locations and name of samples
my %SAMPLE = ();
open(INDEX,"$sampleFile") || die "Cannot open $sampleFile\n";
while (<INDEX>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $hcBamPath) = split('\t', $_);
        $SAMPLE{$sampleID}{HC} = $hcBamPath;
    }
}
close(INDEX);

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;
 
mkpath($outputDir);
my $logDir = "$outputDir/log";
mkpath($logDir);
my $auxDir = "$outputDir/aux";
mkpath($auxDir);
my $finalDir = "$outputDir/final";
mkpath($finalDir);

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

#per sample discovery of sites
my $candidateSitesVCFOutputDir = "$outputDir/vcf";
my $candidateSitesVCFFiles = "";
my $candidateSitesVCFOKFiles = "";
mkdir($candidateSitesVCFOutputDir);
for my $sampleID (keys(%SAMPLE))
{
    mkdir("$candidateSitesVCFOutputDir/$sampleID/");

    my $sampleVCFFile = "$candidateSitesVCFOutputDir/$sampleID/$sampleID.sites.$ext";
    my $sampleDir = "$candidateSitesVCFOutputDir/$sampleID";

    #mine, left align and merge duplicate variants
    $tgt = "$sampleVCFFile.OK";
    $dep = "";
    @cmd = "$vt discover -b $SAMPLE{$sampleID}{HC} -o + -v $variantTypes -r $refGenomeFASTAFile -s $sampleID $intervalsOption 2> $sampleDir/discover.log | $vt normalize + -r $refGenomeFASTAFile -o + 2> $sampleDir/normalize.log | $vt mergedups + -o $sampleVCFFile  2> $sampleDir/mergedups.log";
    makeStep($tgt, $dep, @cmd);
    
    #index
    $tgt = "$sampleVCFFile.$indexExt.OK";
    $dep = "$sampleVCFFile.OK";
    @cmd = "$vt index $sampleVCFFile";
    makeStep($tgt, $dep, @cmd);

    $candidateSitesVCFFiles .= " $sampleVCFFile";
    $candidateSitesVCFOKFiles .= " $sampleVCFFile.$indexExt.OK";
}

#write candidate VCF files into a list
my $candidateVCFFileList = "$auxDir/candidate_vcf_files.txt";
open(IN, ">$candidateVCFFileList") ||die "Cannot open $candidateVCFFileList\n";
my $temp = $candidateSitesVCFFiles;
$temp =~ s/ /\n/g;
$temp =~ s/^\s+//;
print IN "$temp\n";
close(IN);

#merging and filtering of initial sites using likelhood ratio cut off of 2.
my $candidateVCFFile = "$auxDir/all.sites.$ext";
$tgt = "$candidateVCFFile.OK";
$dep =  $candidateSitesVCFOKFiles;
@cmd = ("$vt merge_candidate_variants -L $candidateVCFFileList -o $candidateVCFFile 2> $candidateVCFFile.log");
makeStep($tgt, $dep, @cmd);

#index candidate sites
$tgt = "$candidateVCFFile.$indexExt.OK";
$dep = "$candidateVCFFile.OK";
@cmd = ("$vt index $candidateVCFFile 2> $candidateVCFFile.$indexExt.log");
makeStep($tgt, $dep, @cmd);

##############
#2. Genotyping
##############

#construct probes for candidate sites
my $probesVCFFile = "$auxDir/probes.sites.$ext";
$tgt = "$probesVCFFile.OK";
$dep = "$candidateVCFFile.$indexExt.OK";
@cmd = ("$vt construct_probes $candidateVCFFile -r $refGenomeFASTAFile -o $probesVCFFile 2> $auxDir/probes.log");
makeStep($tgt, $dep, @cmd);

#index probes sites
$tgt = "$probesVCFFile.$indexExt.OK";
$dep = "$probesVCFFile.OK";
@cmd = ("$vt index $probesVCFFile 2> $probesVCFFile.$indexExt.log");
makeStep($tgt, $dep, @cmd);

#per sample discovery of sites
my $candidateSitesVCFOutDir = "$outputDir/vcf";
my $sampleGenotypesVCFFiles;
my $sampleGenotypesVCFIndexOKFiles;
mkdir($candidateSitesVCFOutDir);
for my $sampleID (keys(%SAMPLE))
{
    mkdir("$candidateSitesVCFOutDir/$sampleID/");

    my $sampleDir = "$candidateSitesVCFOutDir/$sampleID";
    my $sampleVCFFile = "$sampleDir/$sampleID.genotypes.$ext";
    my $sampleVCFFileIndex = "$sampleDir/$sampleID.genotypes.$ext.$indexExt";
    
    #genotype
    $tgt = "$sampleVCFFile.OK";
    $dep = "$probesVCFFile.OK";
    @cmd = ("$vt genotype -b $SAMPLE{$sampleID}{HC} $probesVCFFile -r $refGenomeFASTAFile -s $sampleID -o $sampleVCFFile  2> $sampleDir/genotype.log");
    makeStep($tgt, $dep, @cmd);
    
    #index
    $tgt = "$sampleVCFFileIndex.OK";
    $dep = "$sampleVCFFile.OK";
    @cmd = ("$vt index $sampleVCFFile");
    makeStep($tgt, $dep, @cmd);

    $sampleGenotypesVCFFiles .= "$sampleVCFFile\n";
    $sampleGenotypesVCFIndexOKFiles .= " $sampleVCFFileIndex.OK";
}

######################
#3. Merge and Annotate
######################

#make merge list
my $mergedVCFFileList = "$auxDir/merge.vcf.list.txt";
chomp($sampleGenotypesVCFFiles);
open(IN, ">$mergedVCFFileList") || die "Cannot open $mergedVCFFileList";
print IN $sampleGenotypesVCFFiles;
close(IN);

#merge
my $mergedVCFFile = "$finalDir/all.genotypes.$ext";
$tgt = "$mergedVCFFile.OK";
$dep = "$sampleGenotypesVCFIndexOKFiles";
@cmd = ("$vt paste -L $mergedVCFFileList -o + | $vt compute_features + -o + 2> $finalDir/compute_features.log | $vt annotate_dbsnp_rsid - -d $dbsnpVCFFile  2>$finalDir/annotate_dbsnp_rsid.log | $vt remove_overlap - -o $mergedVCFFile 2> $finalDir/remove_overlap.log");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$mergedVCFFile.$indexExt.OK";
$dep = "$mergedVCFFile.OK";
@cmd = ("$vt index $mergedVCFFile");
makeStep($tgt, $dep, @cmd);

#generate site file
my $mergedVCFSitesFile = "$finalDir/all.sites.$ext";
$tgt = "$mergedVCFSitesFile.OK";
$dep = "$mergedVCFFile.OK";
@cmd = ("$vt view -s $mergedVCFFile -o $mergedVCFSitesFile");
makeStep($tgt, $dep, @cmd);

#index
$tgt = "$mergedVCFSitesFile.$indexExt.OK";
$dep = "$mergedVCFSitesFile.OK";
@cmd = ("$vt index $mergedVCFSitesFile");
makeStep($tgt, $dep, @cmd);

#############
#log end time
#############
$tgt = "$logFile.end.OK";
$dep = "$probesVCFFile.$indexExt.OK";
@cmd = ("\tdate | awk '{print \"end: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

####################
#Write out make file
####################

open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
$tgt = "clean";
$dep = "";
@cmd = ("-rm -rf $finalDir/*.OK $auxDir/*.OK $logDir/*.OK $candidateSitesVCFOutputDir/*/*.OK");
makeStep($tgt, $dep, @cmd);

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