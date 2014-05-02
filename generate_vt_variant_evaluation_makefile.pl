#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_sardinia_indel_analysis_makefile

=head1 SYNOPSIS

 generate_gatk_calling_pipeline_makefile [options]  
    
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

#programs
my $vt = "/net/fantasia/home/atks/sardinia/20140421_sardinia_indel_analysis/codes/vt/vt";
my $clusterDir = "/net/fantasia/home/atks/programs/cluster";

#options
my $help;

#cluster options
my $makeFile = "indel.analysis.sardinia.mk";
my $cluster = "mini+";
my $sleep = 60;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help))
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

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;


my @chrom = (1..22, "X", "Y", "MT");

my $refGenomeFASTAFile = "/net/fantasia/home/atks/ref/genome/hs37d5.fa";

my $workDir = "/net/fantasia/home/atks/sardinia/20140421_sardinia_indel_analysis";
my $vcfDir = "$workDir/analyse_indels/vcf";
my $plotDir = "$workDir/analyse_indels/plot";
my $logDir = "$workDir/analyse_indels/log";
mkpath($vcfDir);
mkpath($plotDir);
mkpath($logDir);

###################
#Subset variant set
###################
for my $chrom (@chrom)
{   
    my $sampleFile = "$workDir/2210samples.sardinia.low_coverage.sa";
    my $inputVCFFile = "/net/fantasia/home/atks/sardinia/20140121_variant_calling/run_gatk_ug_indels/vqsr/$chrom.vqsr.genotypes.vcf.gz";
    my $outputVCFFile = "$vcfDir/$chrom.genotypes.bcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "";
    @cmd = ("$vt subset -s $sampleFile  $inputVCFFile -o $outputVCFFile 2> $logDir/$chrom.subset.log");
    makeStep($tgt, $dep, @cmd);
}

#################
#Compute features
#################
my $annotatedSiteVCFFilesOK = "";
my $annotatedSiteVCFFiles = "";
for my $chrom (@chrom)
{  
    $annotatedSiteVCFFilesOK .= " $vcfDir/$chrom.annotated.sites.bcf.OK";
    $annotatedSiteVCFFiles .= " $vcfDir/$chrom.annotated.sites.bcf";
     
    my $inputVCFFile = "$vcfDir/$chrom.genotypes.bcf";
    my $outputVCFFile = "$vcfDir/$chrom.annotated.genotypes.bcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt compute_features $inputVCFFile -o $outputVCFFile 2> $logDir/$chrom.compute_features.log");
    makeStep($tgt, $dep, @cmd);
    
    $inputVCFFile = "$vcfDir/$chrom.annotated.genotypes.bcf";
    $outputVCFFile = "$vcfDir/$chrom.annotated.sites.bcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt view -s $inputVCFFile -o $outputVCFFile");
    makeStep($tgt, $dep, @cmd);
      
    $inputVCFFile = "$vcfDir/$chrom.annotated.genotypes.bcf";
    $tgt = "$inputVCFFile.csi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeStep($tgt, $dep, @cmd);
    
    $inputVCFFile = "$vcfDir/$chrom.annotated.sites.bcf";
    $tgt = "$inputVCFFile.csi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeStep($tgt, $dep, @cmd);
    
}

my $outputVCFFile = "$vcfDir/all.annotated.sites.bcf";
$tgt = "$outputVCFFile.OK";
$dep = $annotatedSiteVCFFilesOK;
@cmd = ("$vt concat $annotatedSiteVCFFiles -o + 2> $logDir/concat.log | $vt annotate_variants -r $refGenomeFASTAFile -g /net/fantasia/home/atks/dev/vt/bundle/public/grch37/gencode.v19.annotation.gtf.gz  + -o $outputVCFFile 2> $logDir/annotate_variants.log");
makeStep($tgt, $dep, @cmd);

$outputVCFFile = "$vcfDir/all.annotated.sites.bcf";
my $outputVCFFileIndex = "$vcfDir/all.annotated.sites.bcf.csi";
$tgt = "$outputVCFFileIndex.OK";
$dep = $annotatedSiteVCFFilesOK;
@cmd = ("$vt index $outputVCFFile");
makeStep($tgt, $dep, @cmd);

###############
#Profile Indels
###############

#Summary
my $plotFile = "$plotDir/summary.pdf";
$tgt = "$plotFile.OK";
$dep = "$vcfDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt peek $vcfDir/all.annotated.sites.bcf -x $plotDir/tabulate_summary -y $plotFile  2> $logDir/summary.log");
makeLocalStep($tgt, $dep, @cmd);

#AFS
$plotFile = "$plotDir/afs.pdf";
$tgt = "$plotFile.OK";
$dep = "$vcfDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt profile_afs $vcfDir/all.annotated.sites.bcf -c VT_AC -n VT_AN -x $plotDir/plot_afs -y $plotFile  2> $logDir/profile_afs.log");
makeLocalStep($tgt, $dep, @cmd);

#HWE
$plotFile = "$plotDir/hwe.pdf";
$tgt = "$plotFile.OK";
$dep = "$vcfDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt profile_hwe $vcfDir/all.annotated.sites.bcf -h VT_HWE_LPVAL -a VT_MLEAF -x $plotDir/plot_hwe -y $plotFile  2> $logDir/profile_hwe.log");
makeLocalStep($tgt, $dep, @cmd);

#Mendelian
my $tabulateFile = "$plotDir/mendelian.pdf";
my $pedigreeFile = "/net/fantasia/home/atks/sardinia/20140421_sardinia_indel_analysis/freeze_120306.sequenced.ped";
$tgt = "$tabulateFile.OK";
$dep = "$vcfDir/20.annotated.genotypes.bcf.OK $vcfDir/20.annotated.genotypes.bcf.csi.OK ";
@cmd = ("$vt profile_mendelian $vcfDir/20.genotypes.bcf -f \"PASS\" -p $pedigreeFile -x $plotDir/tabulate_mendelian -y $tabulateFile  2> $logDir/profile_mendelian.log");
makeLocalStep($tgt, $dep, @cmd);
	
#Overlap
$tabulateFile = "$plotDir/indels.pdf";
$tgt = "$tabulateFile.OK";
$dep = "$vcfDir/all.annotated.sites.bcf.OK $vcfDir/all.annotated.sites.bcf.csi.OK ";
@cmd = ("$vt profile_indels $vcfDir/all.annotated.sites.bcf -f \"PASS&&N_ALLELE==2&&VTYPE==INDEL\" -r $refGenomeFASTAFile -g /net/fantasia/home/atks/ref/vt/grch37/indel.reference.txt -x $plotDir/tabulate_indels -y $tabulateFile  2> $logDir/profile_indels.log");
makeLocalStep($tgt, $dep, @cmd);	
	
#length analysis
$plotFile = "$plotDir/len.pdf";
$tgt = "$plotFile.OK";
$dep = "$vcfDir/all.annotated.sites.bcf.OK $vcfDir/all.annotated.sites.bcf.csi.OK ";
@cmd = ("$vt profile_len $vcfDir/all.annotated.sites.bcf -a VT_HWEAF -b VT_AB -x $plotDir/plot_len -y $plotFile  2> $logDir/profile_len.log");
makeLocalStep($tgt, $dep, @cmd);	
	
	
#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";
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