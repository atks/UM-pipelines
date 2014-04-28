#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_gatk_vqsr_pipeline_makefile

=head1 SYNOPSIS

 generate_gatk_vqsr_pipeline_makefile [options]

  -s     sample file list giving the location of each sample
         column 1: sample name
         column 2: path of bam file
  -r     reference genome file
  -l     sequence length file
  -o     output directory
  -m     make file name

=head1 DESCRIPTION

This script generates the make file to filter a set of variants using GATKs VQSR.

=cut

my $help;

my $inputDir = "";
my $outputDir = "";
my $vtDir = "";
my $clusterDir = "";
my $makeFile = "Makefile";
my $cluster = "main";
my $sleep = 0;
my $sequenceFile = "/net/fantasia/home/atks/dev/vt/pipeline/seq_length.txt";
my $refGenomeFASTAFile = "";
my $jvmMemory = "2g";
my $variantType = "BOTH";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'i:s'=>\$inputDir,
                'o:s'=>\$outputDir,
                'b:s'=>\$vtDir,
                't:s'=>\$clusterDir,
                'm:s'=>\$makeFile,
                'c:s'=>\$cluster,
                'd:s'=>\$sleep,
                'l:s'=>\$sequenceFile,
                'r:s'=>\$refGenomeFASTAFile,
                'j:s'=>\$jvmMemory,
                'v:s'=>\$variantType
                )
  || !defined($inputDir)
  || !defined($outputDir)
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
#my $gatk = "/net/fantasia/home/atks/programs/jdk1.7.0_25/bin/java -jar -Xmx$jvmMemory /net/fantasia/home/atks/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar";
my $gatk = "/net/fantasia/home/atks/programs/jdk1.7.0_25/bin/java -jar -Xmx$jvmMemory /net/fantasia/home/atks/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar";
my $gatk64g = "/net/fantasia/home/atks/programs/jdk1.7.0_25/bin/java -jar -Xmx64g /net/fantasia/home/atks/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar";
my $vt = "$vtDir/vt";

printf("generate_gatk_ug_calling_makefile.pl\n");
printf("\n");
printf("options: input dir            %s\n", $inputDir);
printf("         out dir              %s\n", $outputDir);
printf("         vt path              %s\n", $vt);
printf("         cluster path         %s\n", $clusterDir);
printf("         make file            %s\n", $makeFile);
printf("         cluster              %s\n", $cluster);
printf("         sleep                %s\n", $sleep);
printf("         sequence file        %s\n", $sequenceFile);
printf("         reference            %s\n", $refGenomeFASTAFile);
printf("         JVM Memory           %s\n", $jvmMemory);
printf("\n");

mkpath($outputDir);

####################################
#Read sequences and make directories
####################################
my $vqsrDir = "$outputDir/vqsr";
mkpath($vqsrDir);
my $auxVQSRDir = "$outputDir/vqsr_aux";
mkpath($auxVQSRDir);

my @SEQ = ();
open(SQ,"$sequenceFile") || die "Cannot open $sequenceFile\n";
while (<SQ>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $len) = split('\t', $_);

        push(@SEQ, $chrom);

    }
}
close(SQ);
print "added " . @SEQ . " sequences\n";

my @tgts = ();
my @deps = ();
my @cmds = ();
my $cmd;

###############
#log start time
###############
my $logFile = "$auxVQSRDir/run.vqsr.log";
push(@tgts,"$logFile.start.OK");
push(@deps, "");
$cmd = "\tdate | awk '{print \"gatk vqsr pipeline\\n\\nstart: \"\$\$0}' > $logFile\n";
$cmd = $cmd . "\ttouch $logFile.start.OK\n";
push(@cmds, $cmd);

#########################
#Merge variant sites list
#########################
my $mergedSitesVCFFile = "$auxVQSRDir/merged_sites_vcf_list.txt";
my $vqsrSitesVCFFilesIndicesOK = "";
open(MERGED_SITES_VCF, ">$mergedSitesVCFFile") || die "Cannot open $mergedSitesVCFFile";
for my $chrom (@SEQ)
{
    print MERGED_SITES_VCF "$inputDir/$chrom.sites.vcf.gz\n";
    $vqsrSitesVCFFilesIndicesOK .= " $vqsrDir/$chrom.vqsr.sites.vcf.gz.tbi.OK";
}

#######################
#Concatenate sites list
#######################
my $allSitesVCFFile = "$auxVQSRDir/all.sites.vcf.gz";
push(@tgts,"$allSitesVCFFile.OK");
push(@deps, $mergedSitesVCFFile);
$cmd = "\t$vt concat -L $mergedSitesVCFFile -o $allSitesVCFFile\n";
$cmd .= "\ttouch $allSitesVCFFile.OK\n";
push(@cmds, $cmd);

push(@tgts,"$allSitesVCFFile.tbi.OK");
push(@deps, "$allSitesVCFFile.OK");
$cmd = "\t$vt index $allSitesVCFFile\n";
$cmd .= "\ttouch $allSitesVCFFile.tbi.OK\n";
push(@cmds, $cmd);

###################
##VQSR SNP Training
###################
#if ($variantType eq "BOTH" || $variantType eq "SNP")
#{
#    push(@tgts,"$auxVQSRDir/recalibrate_SNP.recal.OK");
#    push(@deps," $allSitesVCFFile.OK");
#    $cmd = "\t$gatk64g -T VariantRecalibrator " .
#                   "-R $refGenomeFASTAFile " .
#                   "-input $allSitesVCFFile " .
#                   "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/hapmap_3.3.b37.vcf.gz " .
#                   "-resource:omni,known=false,training=true,truth=true,prior=12.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/1000G_omni2.5.b37.vcf.gz " .
#                   "-resource:1000G,known=false,training=true,truth=false,prior=10.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/1000G_phase1.snps.high_confidence.b37.vcf.gz " .
#                   "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz " .
#                   "-an DP " .
#                   "-an QD " .
#                   "-an FS " .
#                   "-an MQRankSum " .
#                   "-an ReadPosRankSum " .
#                   "-mode SNP " .
#                   "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 " .
#                   "-recalFile $auxVQSRDir/recalibrate_SNP.recal " .
#                   "-tranchesFile $auxVQSRDir/recalibrate_SNP.tranches " .
#                   "-rscriptFile $auxVQSRDir/recalibrate_SNP_plots.R " .
#                   "\n";
#    $cmd = "\t" . $cmd . "\n";
#    $cmd .= "\n\ttouch $auxVQSRDir/recalibrate_SNP.recal.OK\n";
#    push(@cmds, $cmd);
#}

########################
##VQSR SNP Recalibration
########################
my @sortedChromosomes = sort {if ($a=~/^\d+$/ && $b=~/^\d+$/){$a<=>$b} else { if ($a eq "MT") {return 1} elsif($b eq "MT") {return -1}else{$a cmp $b} }} @SEQ;

#if ($variantType eq "BOTH" || $variantType eq "SNP")
#{
#    for my $chrom (@sortedChromosomes)
#    {
#         #recalibrate
#        my $outputVQSRVCFFile = ($variantType eq "BOTH") ? "$auxVQSRDir/$chrom.snps.vcf.gz" : "$vqsrDir/$chrom.vqsr.genotypes.vcf.gz";
#        my $dependencyVCFFile = "$inputDir/$chrom.vcf.gz";
#                
#        push(@tgts,"$outputVQSRVCFFile.OK");
#        push(@deps,"$vqsrDir/recalibrate_SNP.recal.OK $dependencyVCFFile.OK");
#        $cmd = "\t$gatk -T ApplyRecalibration " .
#                       "-R $refGenomeFASTAFile " .
#                       "-input $dependencyVCFFile " .
#                       "-mode SNP " .
#                       "-ts_filter_level 99.0 " .
#                       "-recalFile $vqsrDir/recalibrate_SNP.recal " .
#                       "-tranchesFile $vqsrDir/recalibrate_SNP.tranches " .
#                       "-o $outputVQSRVCFFile";
#        $cmd = "\t" . makeMos($cmd) . "\n";
#        $cmd .= "\n\ttouch $outputVQSRVCFFile.OK\n";
#        push(@cmds, $cmd);
#    
#        if ($variantType eq "SNP")
#        {
#            my $outputVQSRSitesVCFFile = $outputVQSRVCFFile;
#            $outputVQSRSitesVCFFile =~ s/genotypes/sites/;
#        
#            push(@tgts,"$outputVQSRSitesVCFFile.OK");
#            push(@deps,"$outputVQSRVCFFile.OK");
#            $cmd = "$vt view -s $outputVQSRVCFFile -o $outputVQSRSitesVCFFile";
#            $cmd = "\t" . makeMos($cmd) . "\n";
#            $cmd .= "\n\ttouch $outputVQSRSitesVCFFile.OK\n";
#            push(@cmds, $cmd);
#        }
#    }
#}

####################
#VQSR INDEL Training
####################
if ($variantType eq "BOTH" || $variantType eq "INDEL")
{
    push(@tgts,"$auxVQSRDir/recalibrate_INDEL.recal.OK");
    push(@deps," $allSitesVCFFile.tbi.OK");
    $cmd = "\t$gatk64g -T VariantRecalibrator " .
                   "-R $refGenomeFASTAFile " .
                   "-input $allSitesVCFFile " .
                   "-resource:hapmap,known=false,training=true,truth=true,prior=12.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz " .
                   "-an DP " .
                   "-an FS " .
                   "-an MQRankSum " .
                   "-an ReadPosRankSum " .
                   "-mode INDEL " .
                   "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 " .
                   "--maxGaussians 4 ".
                   "-recalFile $auxVQSRDir/recalibrate_INDEL.recal " .
                   "-tranchesFile $auxVQSRDir/recalibrate_INDEL.tranches " .
                   "-rscriptFile $auxVQSRDir/recalibrate_INDEL_plots.R " .
                   "\n";
    $cmd .= "\ttouch $auxVQSRDir/recalibrate_INDEL.recal.OK\n";
    push(@cmds, $cmd);
}

#########################
#VQSR INDEL Recalibration
#########################
if ($variantType eq "BOTH" || $variantType eq "INDEL")
{
    for my $chrom (@sortedChromosomes)
    {
        #recalibrate
        my $outputVQSRVCFFile = "$vqsrDir/$chrom.vqsr.genotypes.vcf.gz";
        my $dependencyVCFFile = ($variantType eq "BOTH") ? "$auxVQSRDir/$chrom.snps.vcf.gz" : "$inputDir/$chrom.vcf.gz";
        
        push(@tgts,"$outputVQSRVCFFile.OK");
        push(@deps,"$dependencyVCFFile.OK");
        $cmd = "\t$gatk -T ApplyRecalibration " .
                       "-R $refGenomeFASTAFile " .
                       "-input $dependencyVCFFile " .
                       "-ts_filter_level 99.0 " .
                       "-mode INDEL " .
                       "-recalFile $auxVQSRDir/recalibrate_INDEL.recal " .
                       "-tranchesFile $auxVQSRDir/recalibrate_INDEL.tranches " .
                       "-o $outputVQSRVCFFile";
        $cmd = "\t" . makeMos($cmd) . "\n";
        $cmd .= "\n\ttouch $outputVQSRVCFFile.OK\n";
        push(@cmds, $cmd);
    
        #index
        push(@tgts,"$outputVQSRVCFFile.tbi.OK");
        push(@deps,"$outputVQSRVCFFile.OK");
        $cmd = "$vt index $outputVQSRVCFFile";
        $cmd = "\t" . makeMos($cmd) . "\n";
        $cmd .= "\n\ttouch $outputVQSRVCFFile.tbi.OK\n";
        push(@cmds, $cmd);    
        
        my $outputVQSRSitesVCFFile = $outputVQSRVCFFile;
        $outputVQSRSitesVCFFile =~ s/genotypes/sites/;
    
        #extract sites
        push(@tgts,"$outputVQSRSitesVCFFile.OK");
        push(@deps,"$outputVQSRVCFFile.OK");
        $cmd = "$vt view -s $outputVQSRVCFFile -o $outputVQSRSitesVCFFile";
        $cmd = "\t" . makeMos($cmd) . "\n";
        $cmd .= "\n\ttouch $outputVQSRSitesVCFFile.OK\n";
        push(@cmds, $cmd);
        
        #index sites
        push(@tgts,"$outputVQSRSitesVCFFile.tbi.OK");
        push(@deps,"$outputVQSRSitesVCFFile.OK");
        $cmd = "$vt index $outputVQSRSitesVCFFile";
        $cmd = "\t" . makeMos($cmd) . "\n";
        $cmd .= "\n\ttouch $outputVQSRSitesVCFFile.tbi.OK\n";
        push(@cmds, $cmd);   
        
    }
}

#############
#log end time
#############
push(@tgts,"$logFile.end.OK");
push(@deps, "$vqsrSitesVCFFilesIndicesOK");
$cmd = "\tdate | awk '{print \"end: \"\$\$0}' >> $logFile\n";
$cmd = $cmd . "\ttouch $logFile.end.OK\n";
push(@cmds, $cmd);

#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts,"clean");
push(@deps, "");
$cmd = "\t-rm -rf $outputDir/*.OK\n";
push(@cmds, $cmd);

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
