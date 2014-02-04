#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_vt_genotyping_makefile

=head1 SYNOPSIS

 generate_vt_genotyping_makefile [options]

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
my $ext = "bcf";
my $makeFile = "Makefile";
my $cluster = "main";
my $sleep = 0;
my $sampleFile = "";
my $intervals = "";
my $refGenomeFASTAFile = "";
my $candidateSitesVCFFile;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'w:s'=>\$workDir,
                'o:s'=>\$outputDir,
                'b:s'=>\$vtDir,
                'e:s'=>\$ext,
                'm:s'=>\$makeFile,
                'c:s'=>\$cluster,
                'd:s'=>\$sleep,
                's:s'=>\$sampleFile,
                'g:s'=>\$candidateSitesVCFFile,
                'i:s'=>\$intervals,
                'r:s'=>\$refGenomeFASTAFile
                )
  || !defined($workDir)
  || !defined($makeFile)
  || !defined($sampleFile)
  || !defined($candidateSitesVCFFile)
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
my $vt = "$workDir/$vtDir/vt";
my $vtx =  "/net/fantasia/home/atks/1000g/20130605_phase3_discovery/codes/bin/vtx";

my $intervalsOption = ($intervals ne "") ? ((-e $intervals) ? "-I $intervals" : "-i $intervals") : "";

printf("generate_vt_genotyping_makefile.pl\n");
printf("\n");
printf("options: work dir        %s\n", $workDir);
printf("         out dir         %s\n", $outDir);
printf("         vt path         %s\n", $vt);
printf("         ext             %s\n", $ext);
printf("         make file       %s\n", $makeFile);
printf("         cluster         %s\n", $cluster);
printf("         sleep           %s\n", $sleep);
printf("         sample file     %s\n", $sampleFile);
printf("         candidates VCF  %s\n", $candidateSitesVCFFile);
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
my $logFile = "$outDir/run.genotype.log";
push(@tgts,"$logFile.start.OK");
push(@deps, "");
$cmd = "\tdate | awk '{print \"vt calling pipeline\\n\\nstart: \"\$\$0}' > $logFile\n";
$cmd = $cmd . "\ttouch $logFile.start.OK\n";
push(@cmds, $cmd);


####################
#1. Construct Probes
####################

#construct probes for candidate sites
my $probesVCFFile = "$outDir/probes.sites.vcf.gz";
push(@tgts, "$probesVCFFile.OK");
push(@deps, "");
$cmd = "$vt construct_probes $candidateSitesVCFFile -r $refGenomeFASTAFile -o $probesVCFFile 2> $outDir/probes.log";
$cmd = "\t" . &makeMos($cmd) . "\n";
$cmd .= "\ttouch $probesVCFFile.OK\n";
push(@cmds, $cmd);

#index probes sites
push(@tgts,"$probesVCFFile.tbi.OK");
push(@deps, "$probesVCFFile.OK");
$cmd = "$vt index $probesVCFFile 2> $probesVCFFile.tbi.log";
$cmd = "\t" . &makeMos($cmd) . "\n";
$cmd .= "\ttouch $probesVCFFile.tbi.OK\n";
push(@cmds, $cmd);

##############
#2. Genotyping
##############

#per sample discovery of sites
my $candidateSitesVCFOutDir = "$outDir/vcf";
my $sampleGenotypesVCFFiles;
my $sampleGenotypesVCFIndexOKFiles;
mkdir($candidateSitesVCFOutDir);
for my $sampleID (keys(%SAMPLE))
{
    mkdir("$candidateSitesVCFOutDir/$sampleID/");

    my $sampleDir = "$candidateSitesVCFOutDir/$sampleID";
    my $sampleVCFFile = "$sampleDir/$sampleID.genotypes";
    
    #genotype
    push(@tgts,"$sampleVCFFile.vcf.OK");
    push(@deps,"$probesVCFFile.OK");
    $cmd = "$vtx genotype -b $SAMPLE{$sampleID}{HC} -i $probesVCFFile -g $refGenomeFASTAFile -s $sampleID -o $sampleVCFFile.vcf";
    $cmd = "\t" . &makeMos($cmd) . "\n";
    $cmd .= "\ttouch $sampleVCFFile.vcf.OK\n";
    push(@cmds, $cmd);

    #compress
    push(@tgts,"$sampleVCFFile.$ext.OK");
    push(@deps,"$sampleVCFFile.vcf.OK");
    $cmd = "$vt view $sampleVCFFile.vcf -o $sampleVCFFile.$ext";
    $cmd = "\t" . &makeMos($cmd) . "\n";
    $cmd .= "\ttouch $sampleVCFFile.$ext.OK\n";
    push(@cmds, $cmd);

    #index
    push(@tgts,"$sampleVCFFile.$ext.$indexExt.OK");
    push(@deps,"$sampleVCFFile.$ext.OK");
    $cmd = "$vt index $sampleVCFFile.$ext";
    $cmd = "\t" . &makeMos($cmd) . "\n";
    $cmd .= "\ttouch $sampleVCFFile.$ext.$indexExt.OK\n";
    push(@cmds, $cmd);

    $sampleGenotypesVCFFiles .= "$sampleVCFFile.$ext\n";
    $sampleGenotypesVCFIndexOKFiles .= " $sampleVCFFile.$ext.$indexExt.OK";
}

###################
#3. Merge genotypes
###################

#make merge list
my $mergedVCFFileList = "$outDir/merge.vcf.list.txt";
chomp($sampleGenotypesVCFFiles);
open(IN, ">$mergedVCFFileList") || die "Cannot open $mergedVCFFileList";
print IN $sampleGenotypesVCFFiles;
close(IN);

#merge
my $mergedVCFFile = "$outDir/all.genotypes.$ext";
push(@tgts,"$mergedVCFFile.OK");
push(@deps,"$sampleGenotypesVCFIndexOKFiles");
$cmd = "$vt merge -L $mergedVCFFileList -o $mergedVCFFile";
$cmd = "\t" . &makeMos($cmd) . "\n";
$cmd .= "\ttouch $mergedVCFFile.OK\n";
push(@cmds, $cmd);

#index
push(@tgts,"$mergedVCFFile.$indexExt.OK");
push(@deps,"$mergedVCFFile.OK");
$cmd = "$vt index $mergedVCFFile";
$cmd = "\t" . &makeMos($cmd) . "\n";
$cmd .= "\ttouch $mergedVCFFile.$indexExt.OK\n";
push(@cmds, $cmd);

##############
#log end time
#############
push(@tgts,"$logFile.end.OK");
push(@deps, "$mergedVCFFile.$indexExt.OK");
$cmd = "\tdate | awk '{print \"end: \"\$\$0}' >> $logFile\n";
$cmd = $cmd . "\ttouch $logFile.end.OK\n";
push(@cmds, $cmd);


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
        return ("mosbatch -E/tmp -i -r`/net/fantasia/home/atks/programs/cluster/pick_main_node $sleep` /bin/bash -c 'set pipefail; $cmd'");
    }
    elsif ($cluster eq "mini")
    {
        return ("mosbatch -E/tmp -i -r`/net/fantasia/home/atks/programs/cluster/pick_mini_node $sleep` /bin/bash -c 'set pipefail; $cmd'");
    }
    elsif ($cluster eq "mini+")
    {
        return ("mosbatch -E/tmp -i -r`/net/fantasia/home/atks/programs/cluster/pick_mini+_node $sleep` /bin/bash -c 'set pipefail; $cmd'");
    }    
    else
    {
        print STDERR "$cluster not supported\n";
        exit(1);
    }
}
