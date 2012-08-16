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
my $ref_dir   = ( $config =~ /REF_DIR\s+(\S+)/ )   ? $1 : "/data/genomes/Homo_sapiens/UCSC/hg19";
my $bt2_idx   = ( $config =~ /BT2\s+(\S+)/ )       ? $1 : "$ref_dir/Sequence/BowtieIndex/ucsc.hg19";
my $threads   = ( $config =~ /THREADS\s+(\S+)/ )   ? $1 : "24";
my $memory    = ( $config =~ /MEMORY\s+(\S+)/ )    ? $1 : "48";
my $reads_dir = ( $config =~ /READS_DIR\s+(\S+)/ ) ? $1 : ".";
my $step      = ( $config =~ /STEP\s+(\S+)/ )      ? $1 : "0";
my $gatk      = "$bin/GenomeAnalysisTK.jar";
my $dbsnp     = "$ref_dir/Annotation/Variation/dbsnp.vcf";
my $omni      = "$ref_dir/Annotation/Variation/omni.vcf";
my $hapmap    = "$ref_dir/Annotation/Variation/hapmap.vcf";
my $mills     = "$ref_dir/Annotation/Variation/indels.vcf";
my $ref       = "$ref_dir/Sequence/WholeGenomeFasta/ucsc.hg19.fasta";
die "There are no read 1 fastq reads in $reads_dir. The read 1 reads must be formatted as follows: *_R1.fastq.\n" unless ( `ls $reads_dir/*_R1.fastq*` );
die "There are no read 2 fastq reads in $reads_dir. The read 2 reads must be formatted as follows: *_R2.fastq.\n" unless ( `ls $reads_dir/*_R2.fastq*` );
chomp ( my @reads  = `ls $reads_dir/*fastq*` );
chomp ( my $time   = `date +%T` );

print "Options used     :\n",
      "\tBIN      : $bin\n",
      "\tSNPEFF   : $snpEff\n",
      "\tREF_DIR  : $ref_dir\n",
      "\tBT2_IDX  : $bt2_idx\n",
      "\tTHREADS  : $threads\n",
      "\tMEMORY   : $memory\n",
      "\tREADS_DIR: $reads_dir\n",
      "\tSTEP     : $step\n";

######## End Variables ########

for ( my $i = 0; $i < @reads; $i += 2 )
{
    my ($name)     = $reads[$i] =~ /^.+\/(.+?)_/;
    my $R1         = $reads[$i];
    my $R2         = $reads[$i+1];
    my $JAVA_pre   = "java -Xmx${memory}g -jar";
    my $GATK_pre   = "$JAVA_pre $gatk -T";
    my $snpFilters = "-filter 'QD < 2.0' -filterName 'QD' ".
                     "-filter 'DP < 8' -filterName 'DP' ".
                     "-filter 'MQ < 40.0' -filterName 'MQ' ".
                     "-filter 'FS > 60.0' -filterName 'FS' ".
                     "-filter 'HaplotypeScore > 13.0' -filterName 'HaplotypeScore' ".
                     "-filter 'MQRankSum < -12.5' -filterName 'MQRankSum' ".
                     "-filter 'ReadPosRankSum < -8.0' -filterName 'ReadPosRankSum'";
    my $indelFilts = "-filter 'QD < 2.0' -filterName 'QD' ".
                     "-filter 'DP < 8' -filterName 'DP' ".
                     "-filter 'ReadPosRankSum < -20.0' -filterName 'ReadPosRankSum' ".
                     "-filter 'FS > 200.0' -filterName 'FS' ".
                     "-filter 'InbreedingCoeff < -0.8' -filterName 'InbreedingCoeff'";
    my @steps      = (
                       "bowtie2 --very-sensitive-local -x $bt2_idx -p $threads -1 $R1 -2 $R2 -S $name.sam",
                       "samtools view -bS $name.sam -o $name.bam",
                       "$JAVA_pre $bin/SortSam.jar INPUT=$name.bam OUTPUT=$name.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT",
                       "$JAVA_pre $bin/MarkDuplicates.jar I=$name.sorted.bam O=$name.dup_removed.bam REMOVE_DUPLICATES=true M=mark_dups_metrics_file",
                       "$JAVA_pre $bin/AddOrReplaceReadGroups.jar I=$name.sorted.bam O=$name.fixed_RG.bam SO=coordinate RGID=$name RGLB=$name RGPL=illumina RGPU=$name RGSM=$name VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
                       "$GATK_pre RealignerTargetCreator -R $ref -I $name.fixed_RG.bam -known $dbsnp -o $name.indel_realigner.intervals",
                       "$GATK_pre IndelRealigner -R $ref -I $name.fixed_RG.bam -known $dbsnp -o $name.indels_realigned.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals $name.indel_realigner.intervals",
                       "$GATK_pre CountCovariates -nt $threads -R $ref --knownSites $dbsnp -I $name.indels_realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -dP illumina -recalFile $name.recal.csv",
                       "$GATK_pre TableRecalibration -R $ref -I $name.indels_realigned.bam --out $name.recalibrated.bam -recalFile $name.recal.csv",
                       "samtools index $name.recalibrated.bam",
                       "$GATK_pre UnifiedGenotyper -nt $threads -R $ref -I $name.recalibrated.bam -o $name.raw.snvs.vcf   -glm SNP   -D $dbsnp",
                       "$GATK_pre UnifiedGenotyper -nt $threads -R $ref -I $name.recalibrated.bam -o $name.raw.indels.vcf -glm INDEL -D $mills",
                       "$GATK_pre VariantRecalibrator -R $ref -nt $threads -input $name.raw.snvs.vcf   -mG 6 -mode SNP   -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 $omni -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $dbsnp -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ -an FS -an DP -recalFile $name.snvs.recal.out -tranchesFile $name.snvs.tranches.out",
                       "$GATK_pre VariantRecalibrator -R $ref -nt $threads -input $name.raw.indels.vcf -mG 6 -mode INDEL -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $mills -an QD -an ReadPosRankSum -an FS -recalFile $name.indels.recal.out -tranchesFile $name.indels.tranches.out",
                       "$GATK_pre ApplyRecalibration -R $ref -input $name.raw.snvs.vcf   -ts_filter_level 99.0 -tranchesFile $name.snvs.tranches.out   -recalFile $name.snvs.recal.out   -o $name.recalibrated.snvs.vcf",
                       "$GATK_pre ApplyRecalibration -R $ref -input $name.raw.indels.vcf -ts_filter_level 99.0 -tranchesFile $name.indels.tranches.out -recalFile $name.indels.recal.out -o $name.recalibrated.indels.vcf",
                       "$GATK_pre VariantFiltration -R $ref -V $name.recalibrated.snvs.vcf   -o $name.all.snvs.vcf   $snpFilters",
                       "$GATK_pre VariantFiltration -R $ref -V $name.recalibrated.indels.vcf -o $name.all.indels.vcf $indelFilts",
                       "grep -P '\sTruth' $name.all.snvs.vcf > $name.hard.snvs.vcf",
                       "grep -P '^#'      $name.all.snvs.vcf > $name.pass.snvs.vcf",
                       "grep PASS         $name.all.snvs.vcf >> $name.pass.snvs.vcf",
                       "grep -P '\sTruth' $name.all.indels.vcf> $name.hard.indels.vcf",
                       "grep -P '^#'      $name.all.indels.vcf> $name.pass.indels.vcf",
                       "grep PASS         $name.all.indels.vcf>> $name.pass.indels.vcf",
                       "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$name.pass.snvs.html   -v -i vcf -o txt hg19 $name.pass.snvs.vcf   > $name.pass.snvs.txt",
                       "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$name.hard.snvs.html   -v -i vcf -o txt hg19 $name.hard.snvs.vcf   > $name.hard.snvs.txt",
                       "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$name.pass.indels.html -v -i vcf -o txt hg19 $name.pass.indels.vcf > $name.pass.indels.txt",
                       "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$name.hard.indels.html -v -i vcf -o txt hg19 $name.hard.indels.vcf > $name.hard.indels.txt",
                     );

    chomp ( $time = `date +%T` );
    print "[$time][- / -] Working on sample $name.\n";
    my $nom;
    for ( my $i = $step; $i < @steps; $i++ )
    {
        $nom = sprintf ( "%02d", $i );
        my $current_step = $steps[$i];
        chomp ( $time = `date +%T` );
        my ($clean_step) = $current_step;
        $clean_step =~ s/ -/\n                  -/g if length ($clean_step) > 256;
        print "[$time][$nom/$#steps] Running this step: \n\n", " "x18, "$clean_step\n\n";
        print "current_step: $current_step\n";
        system ( $current_step );
    }
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
      MEMORY    48                                       Amount of memory, in gigabytes, to use
      READS_DIR N/A                                      Absolute location of the reads that are going to be used
      STEP      0                                        The step to start at in the pipeline (0-indexed).     

USAGE
}
