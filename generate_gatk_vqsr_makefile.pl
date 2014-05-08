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
printf("         variant type         %s\n", $variantType);
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
my $tgt;
my $dep;
my @cmd;

###############
#log start time
###############
my $logFile = "$auxVQSRDir/run.vqsr.log";

$tgt = "$logFile.start.OK";
$dep = "";
@cmd = ("date | awk '{print \"gatk vqsr pipeline\\n\\nstart: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

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
$tgt = "$allSitesVCFFile.OK";
$dep = "";
@cmd = ("$vt concat -L $mergedSitesVCFFile -o $allSitesVCFFile");
makeStep($tgt, $dep, @cmd);

$tgt = "$allSitesVCFFile.tbi.OK";
$dep = "$allSitesVCFFile.OK";
@cmd = ("$vt index $allSitesVCFFile");
makeStep($tgt, $dep, @cmd);

###################
##VQSR SNP Training
###################
if ($variantType eq "BOTH" || $variantType eq "SNP")
{
    $tgt = "$logFile.start.snp.recalibration.OK";
    $dep = "$allSitesVCFFile.tbi.OK";
    @cmd = ("date | awk '{print \"start snp recalibration: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);

    $tgt = "$auxVQSRDir/recalibrate_SNP.recal.OK";
    $dep = "$logFile.start.snp.recalibration.OK";
    @cmd = ("$gatk64g -T VariantRecalibrator " .
                   "-R $refGenomeFASTAFile " .
                   "-input $allSitesVCFFile " .
                   "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/hapmap_3.3.b37.vcf.gz " .
                   "-resource:omni,known=false,training=true,truth=true,prior=12.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/1000G_omni2.5.b37.vcf.gz " .
                   "-resource:1000G,known=false,training=true,truth=false,prior=10.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/1000G_phase1.snps.high_confidence.b37.vcf.gz " .
                   "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz " .
                   "-an DP " .
                   "-an QD " .
                   "-an FS " .
                   "-an MQRankSum " .
                   "-an ReadPosRankSum " .
                   "-mode SNP " .
                   "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 " .
                   "-recalFile $auxVQSRDir/recalibrate_SNP.recal " .
                   "-tranchesFile $auxVQSRDir/recalibrate_SNP.tranches " .
                   "-rscriptFile $auxVQSRDir/recalibrate_SNP_plots.R ");
    makeLocalStep($tgt, $dep, @cmd);

    $tgt = "$logFile.end.snp.recalibration.OK";
    $dep = "$auxVQSRDir/recalibrate_SNP.recal.OK";
    @cmd = ("date | awk '{print \"end snp recalibration: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}

####################
#VQSR INDEL Training
####################
if ($variantType eq "BOTH" || $variantType eq "INDEL")
{
    $tgt = "$logFile.start.indel.recalibration.OK";
    $dep = "$allSitesVCFFile.tbi.OK";
    @cmd = ("date | awk '{print \"start indel recalibration: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);

    $tgt = "$auxVQSRDir/recalibrate_INDEL.recal.OK";
    $dep = "$logFile.start.indel.recalibration.OK";
    @cmd = ("$gatk64g -T VariantRecalibrator " .
                   "-R $refGenomeFASTAFile " .
                   "-input $allSitesVCFFile " .
                   "-resource:mills,known=false,training=true,truth=true,prior=12.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz " .
                   "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /net/fantasia/home/atks/ref/gatk/2.8_b37/dbsnp_138.b37.vcf.gz " .
                   "-an QD " .
                   "-an DP " .
                   "-an FS " .
                   "-an ReadPosRankSum " .
                   "-an MQRankSum " .
                   "-an InbreedingCoeff " .
                   "-mode INDEL " .
                   "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 " .
                   "--maxGaussians 4 ".
                   "-recalFile $auxVQSRDir/recalibrate_INDEL.recal " .
                   "-tranchesFile $auxVQSRDir/recalibrate_INDEL.tranches " .
                   "-rscriptFile $auxVQSRDir/recalibrate_INDEL_plots.R "
              );
    makeStep($tgt, $dep, @cmd);

    $tgt = "$logFile.end.indel.recalibration.OK";
    $dep = "$auxVQSRDir/recalibrate_INDEL.recal.OK";
    @cmd = ("date | awk '{print \"end indel recalibration: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}

####################
##VQSR Recalibration
####################
my @sortedChromosomes = sort {if ($a=~/^\d+$/ && $b=~/^\d+$/){$a<=>$b} else { if ($a eq "MT") {return 1} elsif($b eq "MT") {return -1}else{$a cmp $b} }} @SEQ;

if ($variantType eq "BOTH")
{
    $tgt = "$logFile.start.snp.indel.filtering.OK ";
    $dep = "$logFile.end.snp.recalibration.OK $logFile.end.indel.recalibration.OK";
    @cmd = ("date | awk '{print \"start filtering: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);

    my $VQSRSiteVCFFileIndicesOK = "";

    for my $chrom (@sortedChromosomes)
    {
        #recalibrate SNPs
        my $outputVCFFile = "$auxVQSRDir/$chrom.snps.vcf.gz";
        my $inputVCFFile = "$inputDir/$chrom.genotypes.vcf.gz";

        $tgt = "$outputVCFFile.OK";
        $dep = "$auxVQSRDir/recalibrate_SNP.recal.OK $inputVCFFile";
        @cmd = ("$gatk -T ApplyRecalibration " .
                       "-R $refGenomeFASTAFile " .
                       "-input $inputVCFFile " .
                       "-mode SNP " .
                       "-ts_filter_level 99.0 " .
                       "-recalFile $auxVQSRDir/recalibrate_SNP.recal " .
                       "-tranchesFile $auxVQSRDir/recalibrate_SNP.tranches " .
                       "-o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);

        #recalibrate Indels
        $outputVCFFile = "$vqsrDir/$chrom.vqsr.genotypes.vcf.gz";
        $inputVCFFile = "$auxVQSRDir/$chrom.snps.vcf.gz";

        $tgt = "$outputVCFFile.OK";
        $dep = "$auxVQSRDir/recalibrate_INDEL.recal.OK $inputVCFFile.OK";
        @cmd = ("$gatk -T ApplyRecalibration " .
                       "-R $refGenomeFASTAFile " .
                       "-input $inputVCFFile " .
                       "-ts_filter_level 99.0 " .
                       "-mode INDEL " .
                       "-recalFile $auxVQSRDir/recalibrate_INDEL.recal " .
                       "-tranchesFile $auxVQSRDir/recalibrate_INDEL.tranches " .
                       "-o $outputVCFFile");
        makeStep($tgt, $dep, @cmd);

        #index
        my $VQSRVCFFile = "$vqsrDir/$chrom.vqsr.genotypes.vcf.gz";
        $tgt = "$VQSRVCFFile.tbi.OK";
        $dep = "$VQSRVCFFile.OK";
        @cmd = ("$vt index $VQSRVCFFile");
        makeStep($tgt, $dep, @cmd);

        my $VQSRSitesVCFFile = "$vqsrDir/$chrom.vqsr.sites.vcf.gz";

        #extract sites
        $tgt = "$VQSRSitesVCFFile.OK";
        $dep = "$VQSRVCFFile.OK";
        @cmd = ("$vt view -s $VQSRVCFFile -o $VQSRSitesVCFFile");
        makeStep($tgt, $dep, @cmd);

        #index sites
        $tgt = "$VQSRSitesVCFFile.tbi.OK";
        $dep = "$VQSRSitesVCFFile.OK";
        @cmd = ("$vt index $VQSRSitesVCFFile");
        makeStep($tgt, $dep, @cmd);

        $VQSRSiteVCFFileIndicesOK .= " $VQSRSitesVCFFile.tbi.OK";
    }

    $tgt = "$logFile.end.snp.indel.filtering.OK ";
    $dep = "$VQSRSiteVCFFileIndicesOK";
    @cmd = ("date | awk '{print \"end filtering: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}
else
{
    my $variant = $variantType eq "SNP" ? "snp" : "indel";
    $tgt = "$logFile.start.filtering.OK";
    $dep = "$logFile.end.$variant.recalibration.OK";
    @cmd = ("date | awk '{print \"start filtering: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);

    my $VQSRSiteVCFFileIndicesOK = "";

    for my $chrom (@sortedChromosomes)
    {
        my $outputVCFFile = "$vqsrDir/$chrom.vqsr.genotypes.vcf.gz";
        my $inputVCFFile = "$inputDir/$chrom.genotypes.vcf.gz";

        if ($variantType eq "SNP")
        {

            $tgt = "$outputVCFFile.OK";
            $dep = "$vqsrDir/recalibrate_SNP.recal.OK $inputVCFFile.OK";
            @cmd = ("$gatk -T ApplyRecalibration " .
                           "-R $refGenomeFASTAFile " .
                           "-input $inputVCFFile " .
                           "-mode SNP " .
                           "-ts_filter_level 99.0 " .
                           "-recalFile $vqsrDir/recalibrate_SNP.recal " .
                           "-tranchesFile $vqsrDir/recalibrate_SNP.tranches " .
                           "-o $outputVCFFile");
            makeStep($tgt, $dep, @cmd);
        }
        else
        {
            $tgt = "$outputVCFFile.OK";
            $dep = "$auxVQSRDir/recalibrate_INDEL.recal.OK $inputVCFFile.OK";
            @cmd = ("$gatk -T ApplyRecalibration " .
                           "-R $refGenomeFASTAFile " .
                           "-input $inputVCFFile " .
                           "-ts_filter_level 99.0 " .
                           "-mode INDEL " .
                           "-recalFile $auxVQSRDir/recalibrate_INDEL.recal " .
                           "-tranchesFile $auxVQSRDir/recalibrate_INDEL.tranches " .
                           "-o $outputVCFFile"
                    );
            makeStep($tgt, $dep, @cmd);
        }

        my $VQSRVCFFile = "$vqsrDir/$chrom.vqsr.genotypes.vcf.gz";
        my $VQSRSitesVCFFile = " $vqsrDir/$chrom.vqsr.sites.vcf.gz";

        #index
        $tgt = "$VQSRVCFFile.tbi.OK";
        $dep = "$VQSRVCFFile.OK";
        @cmd = ("$vt index $VQSRVCFFile");
        makeStep($tgt, $dep, @cmd);

        #extract sites
        $tgt = "$VQSRSitesVCFFile.OK";
        $dep = "$VQSRVCFFile.OK";
        @cmd = ("$vt view -s $outputVCFFile -o $VQSRSitesVCFFile");
        makeStep($tgt, $dep, @cmd);

        #index sites
        $tgt = "$VQSRSitesVCFFile.tbi.OK";
        $dep = "$VQSRSitesVCFFile.OK";
        @cmd = ("$vt index $VQSRSitesVCFFile");
        makeStep($tgt, $dep, @cmd);

        $VQSRSiteVCFFileIndicesOK .= " $VQSRSitesVCFFile.tbi.OK";
    }

    $tgt = "$logFile.end.OK";
    $dep = "$vqsrSitesVCFFilesIndicesOK";
    @cmd = ("date | awk '{print \"end recalibration: \"\$\$0}' >> $logFile");
    makeLocalStep($tgt, $dep, @cmd);
}

####################
#Write out make file
####################
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

######
#Clean
######
$tgt = "clean";
$dep = "";
@cmd = ("-rm -rf $outputDir/*.OK");
makeLocalStep($tgt, $dep, @cmd);

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