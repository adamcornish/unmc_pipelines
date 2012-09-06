#!/usr/bin/perl
use warnings;
use strict;

usage() unless $#ARGV == 0;
my $in     = shift;
my $config = `cat $in`;
usage() if $config !~ /READS_DIR/;

######## Start Variables ########

my $bin       = ( $config =~ /BIN\s+(\S+)/ )       ? $1 : "/opt/var_calling";
my $snpEff    = ( $config =~ /SNPEFF\s+(\S+)/ )    ? $1 : "/opt/snpeff";
my $ref_dir   = ( $config =~ /REF_DIR\s+(\S+)/ )   ? $1 : "/safer/genomes/Homo_sapiens/UCSC/hg19";
my $bt2_idx   = ( $config =~ /BT2\s+(\S+)/ )       ? $1 : "$ref_dir/Sequence/BowtieIndex/ucsc.hg19";
my $threads   = ( $config =~ /THREADS\s+(\S+)/ )   ? $1 : "24";
my $reads_dir = ( $config =~ /READS_DIR\s+(\S+)/ ) ? $1 : ".";
my $step      = ( $config =~ /STEP\s+(\S+)/ )      ? $1 : "0";
my $exp_name  = ( $config =~ /NAME\s+(\S+)/ )      ? $1 : "experiment";
my $gatk      = "$bin/GenomeAnalysisTK.jar";
my $dbsnp     = "$ref_dir/Annotation/Variation/dbsnp.vcf";
my $omni      = "$ref_dir/Annotation/Variation/omni.vcf";
my $hapmap    = "$ref_dir/Annotation/Variation/hapmap.vcf";
my $mills     = "$ref_dir/Annotation/Variation/indels.vcf";
my $ref       = "$ref_dir/Sequence/WholeGenomeFasta/ucsc.hg19.fasta";
die "There are no read 1 fastq reads in $reads_dir. The read 1 reads must be formatted as follows: *_R1.fastq.\n" unless ( `ls $reads_dir/*_R1*.fastq` );
die "There are no read 2 fastq reads in $reads_dir. The read 2 reads must be formatted as follows: *_R2.fastq.\n" unless ( `ls $reads_dir/*_R2*.fastq` );
chomp ( my @reads  = `ls $reads_dir/*fastq` );
chomp ( my $time   = `date +%T` );

#TODO: add options for (bwa|bowtie), (cancer|noncancer), (snpeff|annovar), (GATK|mpileup), (exome|rna)

print "Options used     :\n",
      "\tBIN      : $bin\n",
      "\tSNPEFF   : $snpEff\n",
      "\tREF_DIR  : $ref_dir\n",
      "\tBT2_IDX  : $bt2_idx\n",
      "\tTHREADS  : $threads\n",
      "\tREADS_DIR: $reads_dir\n",
      "\tSTEP     : $step\n",
      "\tNAME     : $exp_name\n";

######## End Variables ########

for ( my $i = 0; $i < @reads; $i += 2 )
{
    my ($name)      = $reads[$i] =~ /^.+\/(.+?)_/;
    my $R1          = $reads[$i];
    my $R2          = $reads[$i+1];
    my $div_threads = sprintf ( "%d", $threads / (@reads / 2)); # divide threads by the number of samples so as not to overtax the system
    system ( "nohup run_alignments.pl -1 $R1 -2 $R2 -n $name -b $bin -p $div_threads -x $bt2_idx > run_alignments_$name.nohup &" );
}

my $skip = 1;
while ( $skip )
{
    my @done = `ls *done 2>/dev/null`;
    if ( @done != @reads/2 ) { sleep 10 } 
    else { $skip = 0 };
}

chomp ( my @names = `ls $reads_dir/*fastq* | sed 's/_R.*fastq//' | uniq` );
my $JAVA_pre = "java -jar";
my $GATK_pre = "$JAVA_pre $gatk -T";
my $filters  = "-filter 'QD < 2.0' -filterName 'QD' ".
               "-filter 'DP < 8' -filterName 'DP' ".
               "-filter 'MQ < 35.0' -filterName 'MQ' ".
               "-filter 'FS > 60.0' -filterName 'FS' ".
               "-filter 'HaplotypeScore > 13.0' -filterName 'HaplotypeScore' ".
               "-filter 'MQRankSum < -12.5' -filterName 'MQRankSum' ".
               "-filter 'ReadPosRankSum < -8.0' -filterName 'ReadPosRankSum'".
               "-filter 'InbreedingCoeff < -0.8' -filterName 'InbreedingCoeff'";
my $fixed_RG = "";
foreach my $name ( @names ) { $fixed_RG .= "-I $name.fixed_RG.bam "; }

my @gatk = (
             "$GATK_pre BaseRecalibrator -R $ref -knownSites $dbsnp -o recal_data.grp $fixed_RG",
             "$GATK_pre PrintReads -R $ref -BQSR recal_data.grp -o $exp_name.BQSR.bam $fixed_RG",
            #"$GATK_pre ReduceReads -R $ref -I BQSR.bam -o reduced.bam", # only use this if you're using UnifiedGenotyper; it doesn't work well with HaplotypeCaller
             "$GATK_pre RealignerTargetCreator -R $ref -I $exp_name.BQSR.bam -known $dbsnp -o $exp_name.indel_realigner.intervals",
             "$GATK_pre IndelRealigner -R $ref $fixed_RG -known $dbsnp -o $exp_name.indels_realigned.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals $exp_name.indel_realigner.intervals",
            #"samtools index $exp_name.recalibrated.bam",
            #"$GATK_pre UnifiedGenotyper -nt $threads -R $ref recalibrated.bam -o raw.vcf -glm BOTH -D $dbsnp",
             "$GATK_pre HaplotypeCaller -nt $threads -R $ref $exp_name.recalibrated.bam -o $exp_name.raw.vcf -D $dbsnp -stand_call_conf 50.0 -stand_emit_conf 10.0",
            "$GATK_pre VariantRecalibrator -R $ref -nt $threads -input $exp_name.raw.vcf -mG 6 -mode BOTH -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 $omni -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $dbsnp -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $mills -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ -an FS -an DP -an InbreedingCoeff -recalFile $exp_name.recal.out -tranchesFile $exp_name.tranches.out", #only used with HaplotypeCaller
            #"$GATK_pre VariantRecalibrator -R $ref -nt $threads -input raw.vcf -mG 6 -mode SNP -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 $omni -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $dbsnp -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ -an FS -an DP -recalFile snvs.recal -tranchesFile snvs.model", #only used with UnifiedGenotyper
            #"$GATK_pre VariantRecalibrator -R $ref -nt $threads -input raw.vcf -mG 6 -mode INDEL -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $mills -an QD -an ReadPosRankSum -an FS -recalFile indels.recal -tranchesFile indel.model",#only used with UnifiedGenotyper
             "$GATK_pre ApplyRecalibration -R $ref -input $exp_name.raw.vcf -ts_filter_level 99.0 -tranchesFile $exp_name.tranches.out -recalFile $exp_name.recal.out -o $exp_name.recalibrated.vcf",
             "$GATK_pre VariantFiltration -R $ref -V $exp_name.recalibrated.vcf -o $exp_name.all.vcf $filters",
             "grep -P '\\sTruth' $exp_name.all.vcf > $exp_name.hard.vcf",
             "grep -P '^#'       $exp_name.all.vcf > $exp_name.pass.vcf",
             "grep PASS          $exp_name.all.vcf >> $exp_name.pass.vcf",
             "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$exp_name.pass.html -v -i vcf -o txt hg19 $exp_name.pass.vcf > $exp_name.pass.txt",
             "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$exp_name.hard.html -v -i vcf -o txt hg19 $exp_name.hard.vcf > $exp_name.hard.txt",
           );

chomp ( $time = `date +%T` );
print "[$time] Working on all samples.\n";
for ( my $i = $step; $i < @gatk; $i++ )
{
    my $current_step = $gatk[$i];
    my $nom = sprintf ( "%02d", $i );
    chomp ( $time = `date +%T` );
    my ($clean_step) = $current_step;
    $clean_step =~ s/ -/\n                  -/g if length ($clean_step) > 256;
    print "[$time][$nom/$#gatk] Running this step: \n\n", " "x18, "$clean_step\n\n";
    system ( $current_step );
}

sub usage
{
    die <<USAGE;

    Usage: perl $0 <configuration_file.txt>

    Your configuration file MUST have a READS_DIR specified.

    Configuration options available

      OPTION    Default                                  Description
      BIN       /opt/var_calling                         Absolute location of the Picard Tools and GATK jar files
      SNPEFF    /opt/snpeff                              Absolute location of snpEff and its requisite files
      REF_DIR   /data/genomes/Homo_sapiens/UCSC/hg19/    Absolute location of the reference directory
      BT2       REF_DIR/Sequence/BowtieIndex/ucsc.hg19   Absolute location of the Bowtie2 index
      THREADS   24                                       Number of threads to use in parallelizable modules
      READS_DIR N/A                                      Absolute location of the reads that are going to be used
      STEP      0                                        The step to start at in the pipeline (0-indexed).     
      NAME      experiment                               The name you want to give to this experiment.

USAGE
}
