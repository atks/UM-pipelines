#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_variant_evaluation_pipeline_makefile.pl

=head1 SYNOPSIS

 generate_variant_evaluation_pipeline_makefile [options]  
    
  -s     sample file list to be analysed
  -r     reference genome file
  -l     sequence length file
  -w     interval width
  -o     output directory
  -m     make file name

=head1 DESCRIPTION

This script generates the make file to discovery and genotype a set of individuals.
 
=cut

#programs
my $vtDir = "/net/fantasia/home/atks/programs/vt";
my $clusterDir = "/net/fantasia/home/atks/programs/cluster";

#options
my $help;
my $genotypeVCFFileList;
my $pedigreeFile;
my $refGenomeFASTAFile = "/net/fantasia/home/atks/ref/genome/hs37d5.fa";
my $gencodeAnnotationGTFFile = "/net/fantasia/home/atks/dev/vt/bundle/public/grch37/gencode.v19.annotation.gtf.gz";
my $gencodeCodingRegionsBEDFile = "/net/fantasia/home/atks/dev/vt/bundle/public/grch37/gencode.cds.bed.gz";
my $mdustRegionsBEDFile = "/net/fantasia/home/atks/dev/vt/bundle/public/grch37/mdust.bed.gz";
my $snpsReferenceFile = "/net/fantasia/home/atks/dev/vt/bundle/public/grch37/snp.reference.txt";
my $indelReferenceFile = "/net/fantasia/home/atks/dev/vt/bundle/public/grch37/indel.reference.txt";
my $makeFile;
my $cluster = "mini+";
my $sampleFile;
my $variantType = "BOTH",
my $sleep = 60;
my $outputDir;

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'g:s'=>\$genotypeVCFFileList,
                's:s'=>\$sampleFile,
                'p:s'=>\$pedigreeFile,
                'v:s'=>\$variantType,
                'o:s'=>\$outputDir,
                'm:s'=>\$makeFile,
                'b:s'=>\$vtDir,
                'r:s'=>\$refGenomeFASTAFile,
                't:s'=>\$clusterDir,
                'd:s'=>\$sleep,
                'c:s'=>\$cluster)
  || !defined($outputDir)
  || !defined($makeFile)
  || !defined($genotypeVCFFileList)
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

my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my @cmd;

my $vt = "$vtDir/vt";

my $auxDir = "$outputDir/aux";
my $plotDir = "$outputDir/plot";
my $logDir = "$outputDir/log";
my $logFile = "$outputDir/log/run.log"; 
mkpath($auxDir);
mkpath($plotDir);
mkpath($logDir);

printf("generate_variant_evaluation_pipeline_makefile.pl\n");
printf("\n");
printf("options: output dir               %s\n", $outputDir);
printf("         vt path                  %s\n", $vt);
printf("         cluster path             %s\n", $clusterDir);
printf("         make file                %s\n", $makeFile);
printf("         cluster                  %s\n", $cluster);
printf("         sleep                    %s\n", $sleep);
printf("         genotype VCF file list   %s\n", $genotypeVCFFileList);
printf("         sample file              %s\n", !defined($sampleFile)?"None":$sampleFile);
printf("         pedigree file            %s\n", !defined($pedigreeFile)?"None":$pedigreeFile);
printf("         variant types            %s\n", $variantType);
printf("         reference                %s\n", $refGenomeFASTAFile);
printf("         GENCODE GTF              %s\n", $gencodeAnnotationGTFFile);
printf("         GENCODE CDS              %s\n", $gencodeCodingRegionsBEDFile);
printf("         mdust regions            %s\n", $mdustRegionsBEDFile);
printf("\n");

#read file locations and name of samples
my @CHROM = ();
my %CHROM = ();
open(INDEX,"$genotypeVCFFileList") || die "Cannot open $genotypeVCFFileList\n";
while (<INDEX>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $VCFFile) = split('\t');
        $CHROM{$chrom} = $VCFFile;
        push(@CHROM, $chrom);
    }
}
close(INDEX);

###################
#Subset variant set
###################
my $subsetOption = !defined($sampleFile)?"view":"subset -s $sampleFile";

$tgt = "$logFile.start.subset.OK";
$dep = "";
@cmd = ("date | awk '{print \"start subset: \"\$\$0}' > $logFile");
makeLocalStep($tgt, $dep, @cmd);

my $chromSubsetVCFFilesOK = "";

for my $chrom (@CHROM)
{
    my $inputVCFFile = $CHROM{$chrom};
    my $outputVCFFile = "$auxDir/$chrom.genotypes.bcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "$logFile.start.subset.OK";
    @cmd = ("$vt $subsetOption $inputVCFFile -o $outputVCFFile 2> $logDir/$chrom.subset.log");
    makeStep($tgt, $dep, @cmd);
    
    $chromSubsetVCFFilesOK .= " $outputVCFFile.OK";
}

$tgt = "$logFile.end.subset.OK";
$dep = "$chromSubsetVCFFilesOK";
@cmd = ("date | awk '{print \"end subset: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

##################
##Compute features
##################
$tgt = "$logFile.start.annotation.OK";
$dep = "$logFile.end.subset.OK";
@cmd = ("date | awk '{print \"start annotation: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

my $annotatedSiteVCFFiles = "";
my $annotatedSiteVCFFilesOK = "";
for my $chrom (@CHROM)
{  
    $annotatedSiteVCFFilesOK .= " $auxDir/$chrom.annotated.sites.bcf.OK";
    $annotatedSiteVCFFiles .= " $auxDir/$chrom.annotated.sites.bcf";
     
    my $inputVCFFile = "$auxDir/$chrom.genotypes.bcf";
    my $outputVCFFile = "$auxDir/$chrom.annotated.genotypes.bcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "$logFile.start.annotation.OK";
    @cmd = ("$vt compute_features $inputVCFFile -o + 2> $logDir/$chrom.compute_features.log | " .
      "$vt annotate_variants + -r $refGenomeFASTAFile -g $gencodeCodingRegionsBEDFile -m $mdustRegionsBEDFile -o $outputVCFFile 2>$logDir/$chrom.annotation.log");
    makeStep($tgt, $dep, @cmd);
    
    $inputVCFFile = "$auxDir/$chrom.annotated.genotypes.bcf";
    $outputVCFFile = "$auxDir/$chrom.annotated.sites.bcf";
    $tgt = "$outputVCFFile.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt view -s $inputVCFFile -o $outputVCFFile");
    makeStep($tgt, $dep, @cmd);
      
    $inputVCFFile = "$auxDir/$chrom.annotated.genotypes.bcf";
    $tgt = "$inputVCFFile.csi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeStep($tgt, $dep, @cmd);
    
    $inputVCFFile = "$auxDir/$chrom.annotated.sites.bcf";
    $tgt = "$inputVCFFile.csi.OK";
    $dep = "$inputVCFFile.OK";
    @cmd = ("$vt index $inputVCFFile");
    makeStep($tgt, $dep, @cmd);
    
}

$tgt = "$logFile.end.annotation.OK";
$dep = "$annotatedSiteVCFFilesOK";
@cmd = ("date | awk '{print \"s end annotation: \"\$\$0}' >> $logFile");
makeLocalStep($tgt, $dep, @cmd);

$tgt = "$auxDir/all.annotated.sites.bcf.OK";
$dep = $annotatedSiteVCFFilesOK;
@cmd = ("$vt concat $annotatedSiteVCFFiles -o $auxDir/all.annotated.sites.bcf 2> $logDir/concat.log");
makeStep($tgt, $dep, @cmd);

$tgt = "$auxDir/all.annotated.sites.bcf.csi.OK";
$dep = "$auxDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt index $auxDir/all.annotated.sites.bcf");
makeStep($tgt, $dep, @cmd);

#########
##Profile
#########
my $pdfFiles = "";
my $pdfFilesOK = "";

#Summary
my $plotFile = "$plotDir/summary.pass.pdf";
$pdfFiles .= " $plotFile";
$pdfFilesOK .= " $plotFile.OK";
$tgt = "$plotFile.OK";
$dep = "$auxDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt peek $auxDir/all.annotated.sites.bcf -x $plotDir/tabulate_summary -y $plotFile -f \"PASS\"  2> $logDir/summary.pass.log");
makeLocalStep($tgt, $dep, @cmd);

$plotFile = "$plotDir/summary.fail.pdf";
$pdfFiles .= " $plotFile";
$pdfFilesOK .= " $plotFile.OK";
$tgt = "$plotFile.OK";
$dep = "$auxDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt peek $auxDir/all.annotated.sites.bcf -x $plotDir/tabulate_summary -y $plotFile -f \"~PASS\"  2> $logDir/summary.fail.log");
makeLocalStep($tgt, $dep, @cmd);

##Chromosome
$plotFile = "$plotDir/chromosome.pdf";
$pdfFiles .= " $plotFile";
$pdfFilesOK .= " $plotFile.OK";
$tgt = "$plotFile.OK";
$dep = "$auxDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt profile_chrom $auxDir/all.annotated.sites.bcf -x $plotDir/plot_chrom -y $plotFile  2> $logDir/profile_chrom.log");
makeLocalStep($tgt, $dep, @cmd);

##AFS
$plotFile = "$plotDir/afs.pdf";
$pdfFiles .= " $plotFile";
$pdfFilesOK .= " $plotFile.OK";
$tgt = "$plotFile.OK";
$dep = "$auxDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt profile_afs $auxDir/all.annotated.sites.bcf -c VT_AC -n VT_AN -x $plotDir/plot_afs -y $plotFile  2> $logDir/profile_afs.log");
makeLocalStep($tgt, $dep, @cmd);

##HWE
$plotFile = "$plotDir/hwe.pdf";
$pdfFiles .= " $plotFile";
$pdfFilesOK .= " $plotFile.OK";
$tgt = "$plotFile.OK";
$dep = "$auxDir/all.annotated.sites.bcf.OK";
@cmd = ("$vt profile_hwe $auxDir/all.annotated.sites.bcf -h VT_HWE_LPVAL -a VT_MLEAF -x $plotDir/plot_hwe -y $plotFile  2> $logDir/profile_hwe.log");
makeLocalStep($tgt, $dep, @cmd);

#Mendelian
my $tabulateFile = "$plotDir/mendelian.pdf";
$pdfFiles .= " $tabulateFile";
$pdfFilesOK .= " $tabulateFile.OK";
$tgt = "$tabulateFile.OK";
$dep = "$auxDir/20.annotated.genotypes.bcf.OK $auxDir/20.annotated.genotypes.bcf.csi.OK ";
@cmd = ("$vt profile_mendelian $auxDir/20.genotypes.bcf -f \"PASS\" -p $pedigreeFile -x $plotDir/tabulate_mendelian -y $tabulateFile  2> $logDir/profile_mendelian.log");
makeLocalStep($tgt, $dep, @cmd);
	
#Overlap
$tabulateFile = "$plotDir/overlap.indels.pdf";
$pdfFiles .= " $tabulateFile";
$pdfFilesOK .= " $tabulateFile.OK";
$tgt = "$tabulateFile.OK";
$dep = "$auxDir/all.annotated.sites.bcf.OK $auxDir/all.annotated.sites.bcf.csi.OK ";
@cmd = ("$vt profile_indels $auxDir/all.annotated.sites.bcf -f \"PASS&&N_ALLELE==2&&VTYPE==INDEL\" -r $refGenomeFASTAFile -g $indelReferenceFile -x $plotDir/tabulate_indels -y $tabulateFile  2> $logDir/profile_indels.log");
makeLocalStep($tgt, $dep, @cmd);	

###Overlap
#$tabulateFile = "$plotDir/overlap.snps.pdf";
#$tgt = "$tabulateFile.OK";
#$dep = "$auxDir/all.annotated.sites.bcf.OK $auxDir/all.annotated.sites.bcf.csi.OK ";
#@cmd = ("$vt profile_indels $auxDir/all.annotated.sites.bcf -f \"PASS&&N_ALLELE==2&&VTYPE==INDEL\" -r $refGenomeFASTAFile -g $snpReferenceFile -x $plotDir/tabulate_snps -y $tabulateFile  2> $logDir/profile_snps.log");
#makeLocalStep($tgt, $dep, @cmd);		
	
#length analysis
$plotFile = "$plotDir/len.pdf";
$pdfFiles .= " $plotFile";
$pdfFilesOK .= " $plotFile.OK";
$tgt = "$plotFile.OK";
$dep = "$auxDir/all.annotated.sites.bcf.OK $auxDir/all.annotated.sites.bcf.csi.OK ";
@cmd = ("$vt profile_len $auxDir/all.annotated.sites.bcf -a VT_HWEAF -b VT_AB -x $plotDir/plot_len -y $plotFile  2> $logDir/profile_len.log");
makeLocalStep($tgt, $dep, @cmd);	

#combine reports
$plotFile = "$plotDir/report.pdf";
$tgt = "$plotFile.OK";
$dep = "$pdfFilesOK";
@cmd = ("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$plotFile $pdfFiles");
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