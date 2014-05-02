#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_gatk_calling_pipeline_makefile

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

my $help;

my $workDir = "";
my $outputDir = "run";
my $vtDir = "";
my $clusterDir = "";
my $makeFile = "Makefile";
my $cluster = "main";
my $sleep = 0;
my $sampleFile = "";
my $intervals = "";
my $sequenceLengthFile = "/net/fantasia/home/atks/dev/vt/pipeline/seq_length.txt";
my $intervalWidth = 1000000;
my $refGenomeFASTAFile = "";
my $jvmMemory = "2g";
my $variantType = "BOTH";
my $caller = "UnifiedGenotyper";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                'w:s'=>\$workDir,
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
                'x:s'=>\$caller,
                'j:s'=>\$jvmMemory,
                'v:s'=>\$variantType
                )
  || !defined($workDir)
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
my $gatk = "/net/fantasia/home/atks/programs/jdk1.7.0_25/bin/java -jar -Xmx$jvmMemory /net/fantasia/home/atks/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar";   
my $gatk64g = "/net/fantasia/home/atks/programs/jdk1.7.0_25/bin/java -jar -Xmx64g /net/fantasia/home/atks/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar";   
my $vt = "$vtDir/vt";
my $outDir = "$workDir/$outputDir";

printf("generate_gatk_ug_calling_makefile.pl\n");
printf("\n");
printf("options: work dir             %s\n", $workDir);
printf("         out dir              %s\n", $outDir);
printf("         vt path              %s\n", $vt);
printf("         cluster path         %s\n", $clusterDir);
printf("         make file            %s\n", $makeFile);
printf("         cluster              %s\n", $cluster);
printf("         sleep                %s\n", $sleep);
printf("         sample file          %s\n", $sampleFile);
printf("         sequence length file %s\n", $sequenceLengthFile);
printf("         interval width       %s\n", $intervalWidth);
printf("         reference            %s\n", $refGenomeFASTAFile);
printf("         caller               %s\n", $caller);
printf("         JVM Memory           %s\n", $jvmMemory);
printf("         variant types        %s\n", $variantType);
printf("\n");

mkpath($outDir);

########################################
#Read file locations and name of samples
########################################
my %SAMPLE = ();
open(SA,"$sampleFile") || die "Cannot open $sampleFile\n";
my $bamFiles = "";
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $bamPath) = split(/\s+/, $_);
        $SAMPLE{$sampleID} = $bamPath;
        $bamFiles .= "$bamPath\n";
    }
}
close(SA);

my $bamListFile = "$outDir/bam.list";
open(OUT,">$bamListFile") || die "Cannot open $bamListFile\n";
print OUT $bamFiles;
close(OUT);

print "read in " . scalar(keys(%SAMPLE)) . " samples\n";
my $vcfOutDir = "$outDir/vcf";
mkpath($vcfOutDir);
my $finalVCFOutDir = "$outDir/final";
mkpath($finalVCFOutDir);

###################
#Generate intervals
###################
my %intervalsByChrom = ();
my %intervalsByChromOK = ();
my @intervals = ();
my @intervalFiles = ();

my $writeIntervals = 1;

if (-e "$outDir/intervals/$intervalWidth.OK")
{
    print "$outDir/intervals/$intervalWidth.OK exists, intervals wil not be generated.\n";
    $writeIntervals = 0;
}

mkpath("$outDir/intervals/");
open(SQ,"$sequenceLengthFile") || die "Cannot open $sequenceLengthFile\n";
while (<SQ>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($chrom, $len) = split('\t', $_);
                
        print "processing $chrom\t$len ";
        
        $intervalsByChrom{$chrom} = ();
        my $count = 0;
        for my $i (0 .. floor($len/$intervalWidth))
        {   
            my $interval = "";
            my $file = "";
            if ($i<floor($len/$intervalWidth))
            {
                $interval = $chrom . "_" . ($intervalWidth*$i+1) . "_" . ($intervalWidth*($i+1));
                $file = "$outDir/intervals/$interval.interval_list";
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
                $file = "$outDir/intervals/$interval.interval_list";
                if ($writeIntervals) 
                {
                    open(INTERVAL, ">$file") || die "Cannot open $file\n";
                    print INTERVAL "$chrom:" . ($intervalWidth*$i+1) . "-" . $len . "\n";
                    close(INTERVAL);
                }
            }
            
            push(@{$intervalsByChrom{$chrom}}, "$vcfOutDir/$interval.vcf");
            push(@{$intervalsByChromOK{$chrom}}, "$vcfOutDir/$interval.vcf.OK");
            push(@intervals, $interval);
            push(@intervalFiles, $file);
            
            $count++;
        } 
        
        print "added $count intervals\n";
    }
}
close(SQ);

if ($writeIntervals) 
{
    print `touch $outDir/intervals/$intervalWidth.OK`;
}

my @tgts = ();
my @deps = ();
my @cmds = ();
my $cmd;

###############
#log start time
###############
my $logFile = "$outDir/run.log";
push(@tgts,"$logFile.start.OK");
push(@deps, "");
$cmd = "\tdate | awk '{print \"gatk unified genotyper calling pipeline\\n\\nstart: \"\$\$0}' > $logFile\n";
$cmd = $cmd . "\ttouch $logFile.start.OK\n";
push(@cmds, $cmd);

#########################
#Discovery and Genotyping
#########################
for my $i (0 .. $#intervals)
{   
    push(@tgts,"$vcfOutDir/$intervals[$i].vcf.OK");
    push(@deps,"");
    #nct - number of computing threads
    #interval_padding ensures that you capture Indels that lie across a boundary. Note that UnifiedGenotyper uses locuswalker.
    #--max_alternate_alleles is set at 6 by default
    
    if ($caller eq "UnifiedGenotyper")
    {
        $cmd = "\t$gatk -T UnifiedGenotyper -R $refGenomeFASTAFile -glm $variantType --interval_padding 100 -I $bamListFile --genotyping_mode DISCOVERY -o $vcfOutDir/$intervals[$i].vcf --output_mode EMIT_VARIANTS_ONLY -L $intervalFiles[$i]";
    }
    else
    {
        $cmd = "\t$gatk -T HaplotypeCaller -R $refGenomeFASTAFile --interval_padding 100 -I $bamListFile --genotyping_mode DISCOVERY -o $vcfOutDir/$intervals[$i].vcf -L $intervalFiles[$i]";
    }
    
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\n\ttouch $vcfOutDir/$intervals[$i].vcf.OK\n";
    push(@cmds, $cmd);
}

###########################################
#Concatenate, normalize and drop duplicates
###########################################
open(IN, ">$finalVCFOutDir/merge_vcf_list.txt") || die "Cannot open merge_vcf_list.txt";
my @sortedChromosomes = sort {if ($a=~/^\d+$/ && $b=~/^\d+$/){$a<=>$b} else { if ($a eq "MT") {return 1} elsif($b eq "MT") {return -1}else{$a cmp $b} }} keys(%intervalsByChrom);
for my $chrom (@sortedChromosomes)
{
    my $intervalFilesOK = join(' ', @{$intervalsByChromOK{$chrom}});
    
    push(@tgts,"$finalVCFOutDir/$chrom.vcf.gz.OK");
    push(@deps,"$intervalFilesOK");
    $cmd = "$vt concat " . join(' ', @{$intervalsByChrom{$chrom}})  . " -o + | $vt normalize + -o + -r $refGenomeFASTAFile | $vt mergedups + -o $finalVCFOutDir/$chrom.vcf.gz";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\n\ttouch $finalVCFOutDir/$chrom.vcf.gz.OK\n";
    push(@cmds, $cmd);
    
    print IN "$finalVCFOutDir/$chrom.vcf.gz\n";    
}
close(IN);

my $chromVCFIndicesOK = "";
my $chromSitesVCFIndicesOK = "";
for my $chrom (@sortedChromosomes)
{
    #index main file
    push(@tgts,"$finalVCFOutDir/$chrom.vcf.gz.tbi.OK");
    push(@deps,"$finalVCFOutDir/$chrom.vcf.gz.OK");
    $cmd = "$vt index $finalVCFOutDir/$chrom.vcf.gz";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\n\ttouch $finalVCFOutDir/$chrom.vcf.gz.tbi.OK\n";
    push(@cmds, $cmd);
    
    $chromVCFIndicesOK .= " $finalVCFOutDir/$chrom.vcf.gz.tbi.OK"; 

    #sites    
    push(@tgts,"$finalVCFOutDir/$chrom.sites.vcf.gz.OK");
    push(@deps,"$finalVCFOutDir/$chrom.vcf.gz.tbi.OK");
    $cmd = "$vt view -s $finalVCFOutDir/$chrom.vcf.gz -o $finalVCFOutDir/$chrom.sites.vcf.gz";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\n\ttouch $finalVCFOutDir/$chrom.sites.vcf.gz.OK\n";
    push(@cmds, $cmd);
    
    #index sites
    push(@tgts,"$finalVCFOutDir/$chrom.sites.vcf.gz.tbi.OK");
    push(@deps,"$finalVCFOutDir/$chrom.sites.vcf.gz.OK");
    $cmd = "$vt index $finalVCFOutDir/$chrom.sites.vcf.gz";
    $cmd = "\t" . makeMos($cmd) . "\n";
    $cmd .= "\n\ttouch $finalVCFOutDir/$chrom.sites.vcf.gz.tbi.OK\n";
    push(@cmds, $cmd);
       
    $chromSitesVCFIndicesOK .= " $finalVCFOutDir/$chrom.sites.vcf.gz.tbi.OK";    
}

#############
#log end time
#############
push(@tgts,"$logFile.end.OK");
push(@deps, "$chromVCFIndicesOK");
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
$cmd = "\t-rm -rf $outDir/*.OK $vcfOutDir/*.OK $finalVCFOutDir/*.OK\n";
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
