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
#my $indels    = "$ref_dir/Annotation/Variation/indels.vcf";
my $indels    = $dbsnp;
my $exome_bed = "$ref_dir/Annotation/Genes/refSeq_genes.bed";
my $ref       = "$ref_dir/Sequence/WholeGenomeFasta/ucsc.hg19.fasta";
my $log       = "run.log";
die "There are no read 1 fastq reads in $reads_dir. The read 1 reads must be formatted as follows: *_R1.fastq.\n" unless ( `ls $reads_dir/*_R1.fastq` );
die "There are no read 2 fastq reads in $reads_dir. The read 2 reads must be formatted as follows: *_R2.fastq.\n" unless ( `ls $reads_dir/*_R2.fastq` );
chomp ( my @reads  = `ls $reads_dir/*fastq` );
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
    my ($name) = $reads[$i] =~ /^.+\/(.+?)_/;
    my $R1     = $reads[$i];
    my $R2     = $reads[$i+1];
    my $JAVA_pre = "java -Xmx${memory}g -jar";
    my $GATK_pre   = "$JAVA_pre $gatk -T";
    my @steps      = (
                       "bowtie2 --local -x $bt2_idx -p $threads -1 $R1 -2 $R2 -S $name.sam >> $log.stdout 2>> $log.stderr",
                       "samtools view -bS $name.sam -o $name.bam >> $log.stdout 2>> $log.stderr",
                       "$JAVA_pre $bin/SortSam.jar INPUT=$name.bam OUTPUT=$name.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT >> $log.stdout 2>> $log.stderr",
                       "$JAVA_pre $bin/AddOrReplaceReadGroups.jar I=$name.sorted.bam O=$name.fixed_RG.bam SO=coordinate RGID=$name RGLB=$name RGPL=illumina RGPU=$name RGSM=$name VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre RealignerTargetCreator -R $ref -I $name.fixed_RG.bam -known $indels -o $name.indel_realigner.intervals", # >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre IndelRealigner -R $ref -I $name.fixed_RG.bam -known $indels -o $name.indels_realigned.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals $name.indel_realigner.intervals", # >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre CountCovariates -nt $threads -R $ref --knownSites $dbsnp -I $name.indels_realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -dP illumina -recalFile $name.recal.csv", # >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre TableRecalibration -R $ref -I $name.indels_realigned.bam --out $name.recalibrated.bam -recalFile $name.recal.csv", #>> $log.stdout 2>> $log.stderr",
                       "samtools index $name.recalibrated.bam", # >> $log.stdout 2>> $log.stderr",
                       "java -Xmx8g -jar $gatk -T UnifiedGenotyper -nt $threads -R $ref -I $name.recalibrated.bam -o $name.raw.snvs.vcf   -glm SNP   -mbq 20 -stand_emit_conf 10.0 -D $dbsnp", #  >> $log.stdout 2>> $log.stderr",
                       "java -Xmx8g -jar $gatk -T UnifiedGenotyper -nt $threads -R $ref -I $name.recalibrated.bam -o $name.raw.indels.vcf -glm INDEL -mbq 20 -stand_emit_conf 10.0 -D $indels", # >> $log.stdout 2>> $log.stderr",
                      #"$GATK_pre VariantFiltration -R $ref -V $name.raw_snvs.vcf   -o $name.filtered_snvs.vcf   -filter 'DP <= 20 || HRun > 8 || QD < 5.0' -filterName 'standard_filters' -filter 'MQ0 >= 4 && ((MQ0 / (1.0  * DP)) > 0.1 )' -filterName 'hard_to_validate' >> $log.stdout 2>> $log.stderr",
                      #"$GATK_pre VariantFiltration -R $ref -V $name.raw_indels.vcf -o $name.filtered_indels.vcf -filter 'DP <= 20 || HRun > 8 || QD < 5.0' -filterName 'standard_filters' -filter 'MQ0 >= 4 && ((MQ0 / (1.0  * DP)) > 0.1 )' -filterName 'hard_to_validate' >> $log.stdout 2>> $log.stderr",
                      #"cat $name.filter_snvs.vcf   | grep -P '^#' > $name.pass.snvs.vcf",
                      #"cat $name.filter_indels.vcf | grep -P '^#' > $name.pass.indels.vcf",
                      #"cat $name.filter_snvs.vcf   | grep PASS >> $name.pass.snvs.vcf",
                      #"cat $name.filter_indels.vcf | grep PASS >> $name.pass.indels.vcf",
                       "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$name.snvs.html   -minC 15 -v -i vcf -o txt hg19 $name.raw.snvs.vcf   > $name.snvs.txt",
                       "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$name.indels.html -minC 15 -v -i vcf -o txt hg19 $name.raw.indels.vcf > $name.indels.txt",
                       "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$name.snvs.html   -minC 15 -v -i vcf -o vcf hg19 $name.raw.snvs.vcf   > $name.snvs.vcf",
                       "$JAVA_pre $snpEff/snpEff.jar eff -c $snpEff/snpEff.config -s ./$name.indels.html -minC 15 -v -i vcf -o vcf hg19 $name.raw.indels.vcf > $name.indels.vcf",
                       "$GATK_pre -T VariantRecalibrator -R $ref -nt $threads -input $name.snvs.vcf -mG 6 -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 $omni -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $dbsnp -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ -recalFile recal.out -tranchesFile tranches.out",
                       "$GATK_pre -T ApplyRecalibration -R $ref -input $name.snvs.vcf -ts_filter_level 99.0 -tranchesFile tranches.out -recalFile recal.out -o $name.recalibrated.filtered.snvs.vcf",
                     );

    my $nom    = "00";
    chomp ( $time = `date +%T` );
    print "[$time][- / -] Working on sample $name.\n";
    for ( my $i = $step; $i < @steps; $i++ )
    {
        my $current_step = $steps[$i];
        chomp ( $time = `date +%T` );
        my ($clean_step) = $current_step;
        $clean_step =~ s/ -/\n                  -/g if length ($clean_step) > 256;
        print "[$time][$nom/$#steps] Running this step: \n\n", " "x18, "$clean_step\n\n";
        system ( $current_step );
        $nom = sprintf ( "%02d", ++$nom);
    }
   #system ( "mkdir $name; mv $name.* $name;" ); 
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
