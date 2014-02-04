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

my $workDir = "";
my $outputDir = "run";
my $vtDir = "codes/vt";
my $clusterDir = "";
my $ext = "bcf";
my $makeFile = "Makefile";
my $cluster = "main";
my $sleep = 0;
my $sampleFile = "";
my $intervals = "";
my $variantTypes = "indels";
my $refGenomeFASTAFile = "";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'w:s'=>\$workDir,
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
                'r:s'=>\$refGenomeFASTAFile
                )
  || !defined($workDir)
  || !defined($makeFile)
  || !defined($sampleFile)
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

my $outDir = "$workDir/$outputDir";
my $indexExt = $ext eq "vcf.gz" ? "tbi" : "csi";
my $vt = "$vtDir/vt";

my $intervalsOption = ($intervals ne "") ? ((-e $intervals) ? "-I $intervals" : "-i $intervals") : "";

printf("generate_vt_calling_makefile.pl\n");
printf("\n");
printf("options: work dir        %s\n", $workDir);
printf("         out dir         %s\n", $outDir);
printf("         vt path         %s\n", $vt);
printf("         cluster path    %s\n", $clusterDir);
printf("         ext             %s\n", $ext);
printf("         make file       %s\n", $makeFile);
printf("         cluster         %s\n", $cluster);
printf("         sleep           %s\n", $sleep);
printf("         sample file     %s\n", $sampleFile);
printf("         variant types   %s\n", $variantTypes);
printf("         intervals       %s\n", $intervalsOption eq "" ? "none" : $intervalsOption);
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
my $cmd;

mkdir($outDir);

###############
#log start time
###############
my $logFile = "$outDir/run.log";
push(@tgts,"$logFile.start.OK");
push(@deps, "");
$cmd = "\tdate | awk '{print \"vt calling pipeline\\n\\nstart: \"\$\$0}' > $logFile\n";
$cmd = $cmd . "\ttouch $logFile.start.OK\n";
push(@cmds, $cmd);

#############
#1. Discovery
#############

#per sample discovery of sites
my $candidateSitesVCFOutDir = "$outDir/vcf";
my $candidateSitesVCFFiles = "";
my $candidateSitesVCFOKFiles = "";
mkdir($candidateSitesVCFOutDir);
for my $sampleID (keys(%SAMPLE))
{
    mkdir("$candidateSitesVCFOutDir/$sampleID/");

    my $sampleVCFFile = "$candidateSitesVCFOutDir/$sampleID/$sampleID.sites.$ext";
    my $sampleDir = "$candidateSitesVCFOutDir/$sampleID";

    #mine, left align and merge duplicate variants
    push(@tgts,"$sampleVCFFile.OK");
    push(@deps,"");
    $cmd = "$vt discover -b $SAMPLE{$sampleID}{HC} -o + -v $variantTypes -r $refGenomeFASTAFile -s $sampleID $intervalsOption 2> $sampleDir/discover.log | $vt normalize + -r $refGenomeFASTAFile -o + 2> $sampleDir/normalize.log | $vt mergedups + -o $sampleVCFFile  2> $sampleDir/mergedups.log";
    $cmd = "\t" . &makeMos($cmd) . "\n";
    $cmd .= "\ttouch $sampleVCFFile.OK\n";
    push(@cmds, $cmd);

    #index
    push(@tgts,"$sampleVCFFile.$indexExt.OK");
    push(@deps,"$sampleVCFFile.OK");
    $cmd = "$vt index $sampleVCFFile";
    $cmd = "\t" . &makeMos($cmd) . "\n";
    $cmd .= "\ttouch $sampleVCFFile.$indexExt.OK\n";
    push(@cmds, $cmd);

    $candidateSitesVCFFiles .= " $sampleVCFFile";
    $candidateSitesVCFOKFiles .= " $sampleVCFFile.$indexExt.OK";
}

#write candidate VCF files into a list
my $candidateVCFFileList = "$outDir/candidate_vcf_files.txt";
open(IN, ">$candidateVCFFileList") ||die "Cannot open $candidateVCFFileList\n";
my $temp = $candidateSitesVCFFiles;
$temp =~ s/ /\n/g;
$temp =~ s/^\s+//;
print IN "$temp\n";
close(IN);

#merging and filtering of initial sites using likelhood ratio cut off of 2.
my $candidateVCFFile = "$outDir/all.sites.$ext";
push(@tgts,"$candidateVCFFile.OK");
push(@deps, $candidateSitesVCFOKFiles);
$cmd = "$vt merge_candidate_variants -L $candidateVCFFileList -o $candidateVCFFile 2> $candidateSitesVCFOutDir/candidates.log";
$cmd = "\t" . &makeMos($cmd) . "\n";
$cmd .= "\ttouch $candidateVCFFile.OK\n";
push(@cmds, $cmd);

#index candidate sites
push(@tgts,"$candidateVCFFile.$indexExt.OK");
push(@deps, "$candidateVCFFile.OK");
$cmd = "$vt index $candidateVCFFile 2> $candidateVCFFile.$indexExt.log";
$cmd = "\t" . &makeMos($cmd) . "\n";
$cmd .= "\ttouch $candidateVCFFile.$indexExt.OK\n";
push(@cmds, $cmd);

##############
#2. Genotyping
##############

#construct probes for candidate sites
my $probesVCFFile = "$outDir/probes.sites.$ext";
push(@tgts, "$probesVCFFile.OK");
push(@deps, "$candidateVCFFile.OK");
$cmd = "$vt construct_probes $candidateVCFFile -r $refGenomeFASTAFile -o $probesVCFFile 2> $outDir/probes.log";
$cmd = "\t" . &makeMos($cmd) . "\n";
$cmd .= "\ttouch $probesVCFFile.OK\n";
push(@cmds, $cmd);

#index probes sites
push(@tgts,"$probesVCFFile.$indexExt.OK");
push(@deps, "$probesVCFFile.OK");
$cmd = "$vt index $probesVCFFile 2> $probesVCFFile.$indexExt.log";
$cmd = "\t" . &makeMos($cmd) . "\n";
$cmd .= "\ttouch $probesVCFFile.$indexExt.OK\n";
push(@cmds, $cmd);

#############
#log end time
#############
push(@tgts,"$logFile.end.OK");
push(@deps, "$probesVCFFile.$indexExt.OK");
$cmd = "\tdate | awk '{print \"end: \"\$\$0}' >> $logFile\n";
$cmd = $cmd . "\ttouch $logFile.end.OK\n";
push(@cmds, $cmd);

#per sample genotyping of sites
#for my $sampleID (keys(%SAMPLE))
#{
#    my $sampleVCFFile = "$candidateSitesVCFOutDir/$sampleID/$sampleID.genotypes.$ext";
#
#    #genotype
#    push(@tgts,"$sampleVCFFile.OK");
#    push(@deps,"$probesVCFFile.OK");
#    $cmd = "$vt genotype -b $SAMPLE{$sampleID}{HC} -r $refGenomeFASTAFile -s $sampleID $intervals -o $sampleVCFFile";
#    $cmd = "\t" . &makeMos($cmd) . "\n";
#    $cmd .= "\ttouch $sampleVCFFile.OK\n";
#    push(@cmds, $cmd);
#
#    #index
#    push(@tgts,"$sampleVCFFile.$indexExt.OK");
#    push(@deps,"$sampleVCFFile.OK");
#    $cmd = "$vt index $sampleVCFFile";
#    $cmd = "\t" . &makeMos($cmd) . "\n";
#    $cmd .= "\ttouch $sampleVCFFile.$indexExt.OK\n";
#    push(@cmds, $cmd);
#}

####################
#Write out make file
####################

open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts,"clean");
push(@deps, "");
$cmd = "\t-rm -rf $outDir/*.OK $candidateSitesVCFOutDir/*/*.OK\n";
push(@cmds, $cmd);

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
