#!/usr/bin/perl
use warnings;
use strict;

#
# Configuration options available
#
#     OPTION    Default                                                 Description
#     BIN       /opt/var_calling                                        Absolute location of the Picard Tools and GATK jar files
#     REF_DIR   /safer/genomes/Homo_sapiens/UCSC/hg19                   Absolute location of the reference directory
#     BT2       /safer/genomes/Homo_sapiens/Sequence/BowtieIndex/hg19   Absolute location of the Bowtie2 index
#     THREADS   64                                                      Number of threads to use in parallelizable modules
#     MEMORY    64                                                      Amount of memory, in gigabytes, to use
#     READS_DIR N/A                                                     Absolute location of the reads that are going to be used
#

usage() unless $#ARGV == 0;
my $in     = shift;
my $config = `cat $in`;
usage() if $config !~ /READS_DIR/;

######## Start Variables ########

my $bin       = ( $config =~ /^BIN\s+(\S+)/ )       ? $1 : "/opt/var_calling";
my $ref_dir   = ( $config =~ /^REF_DIR\s+(\S+)/ )   ? $1 : "/safer/genomes/Homo_sapiens/UCSC/hg19";
my $bt2_idx   = ( $config =~ /^BT2\s+(\S+)/ )       ? $1 : "$ref_dir/Sequence/BowtieIndex/hg19";
my $threads   = ( $config =~ /^THREADS\s+(\S+)/ )   ? $1 : "64";
my $memory    = ( $config =~ /^MEMORY\s+(\S+)/ )    ? $1 : "64";
my $reads_dir = ( $config =~ /^READS_DIR\s+(\S+)/ ) ? $1 : ".";
my $gatk      = "$bin/GenomeAnalysisTK.jar";
my $dbsnp     = "$ref_dir/Annotation/Variation/dbsnp.vcf";
my $indels    = "$ref_dir/Annotation/Variation/indels.vcf";
my $exome_bed = "$ref_dir/Annotation/Genes/refSeq_genes.bed";
my $ref       = "$ref_dir/Sequence/WholeGenomeFasta/ref.fa";
my $log       = "run.log";
die "There are no read 1 fastq reads in $reads_dir. The read 1 reads must be formatted as follows: *_R1.fastq.\n" unless ( `ls $reads_dir/*_R1.fastq` );
die "There are no read 2 fastq reads in $reads_dir. The read 2 reads must be formatted as follows: *_R2.fastq.\n" unless ( `ls $reads_dir/*_R2.fastq` );
chomp ( my @reads  = `ls $reads_dir/*fastq` );
chomp ( my $time   = `date +%T` );

######## End Variables ########

for ( my $i = 0; $i < @reads; $i += 2 )
{
    my ($name) = $reads[$i] =~ /^(.+?)_/;
    my $R1     = $reads[$i];
    my $R2     = $reads[$i+1];
    my $JAVA_pre = "java -Xmx${memory}g -jar";
    my $GATK_pre   = "$JAVA_pre $gatk -T";
    my @steps      = (
                       "bowtie2 --local -x $bt2_idx -p $threads -1 $R1 -2 $R2 -S $name.sam >> $log.stdout 2>> $log.stderr",
                       "samtools view -bS $name.sam -o $name.bam >> $log.stdout 2>> $log.stderr",
                       "$JAVA_pre $bin/SortSam.jar INPUT=$name.bam OUTPUT=$name.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT >> $log.stdout 2>> $log.stderr",
                       "$JAVA_pre $bin/AddOrReplaceReadGroups.jar I=$name.sorted.bam O=$name.fixed_RG.bam SO=coordinate RGID=$name RGLB=$name RGPL=illumina RGPU=$name RGSM=$name VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre RealignerTargetCreator -R $ref -I $name.fixed_RG.bam -known $indels -o $name.indel_realigner.intervals >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre IndelRealigner -R $ref -I $name.fixed_RG.bam -known $indels -o $name.indels_realigned.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals $name.indel_realigner.intervals >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre CountCovariates -nt $threads -R $ref --knownSites $dbsnp -I $name.indels_realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -dP illumina -recalFile $name.recal.csv >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre TableRecalibration -R $ref -I $name.indels_realigned.bam --out $name.recalibrated.bam -recalFile $name.recal.csv >> $log.stdout 2>> $log.stderr",
                       "samtools index $name.recalibrated.bam >> $log.stdout 2>> $log.stderr",
                       "java -Xmx8g -jar $gatk -T UnifiedGenotyper -nt $threads -R $ref -I $name.recalibrated.bam -o $name.raw_snvs.vcf              -mbq 20 -stand_call_conf 40.0 -stand_emit_conf 10.0 -dcov 1000 -D $dbsnp  -L $exome_bed >> $log.stdout 2>> $log.stderr",
                       "java -Xmx8g -jar $gatk -T UnifiedGenotyper -nt $threads -R $ref -I $name.recalibrated.bam -o $name.raw_indels.vcf -glm INDEL -mbq 20 -stand_call_conf 40.0 -stand_emit_conf 10.0 -dcov 1000 -D $indels -L $exome_bed >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre VariantFiltration -R $ref -V $name.raw_snvs.vcf   -o $name.filtered_snvs.vcf   -filter 'DP <= 20 || HRun > 8 || QD < 5.0' -filterName 'standard_filters' -filter 'MQ0 >= 4 && ((MQ0 / (1.0  * DP)) > 0.1 )' -filterName 'hard_to_validate' >> $log.stdout 2>> $log.stderr",
                       "$GATK_pre VariantFiltration -R $ref -V $name.raw_indels.vcf -o $name.filtered_indels.vcf -filter 'DP <= 20 || HRun > 8 || QD < 5.0' -filterName 'standard_filters' -filter 'MQ0 >= 4 && ((MQ0 / (1.0  * DP)) > 0.1 )' -filterName 'hard_to_validate' >> $log.stdout 2>> $log.stderr",
                     );

    my $nom    = "00";
    chomp ( $time = `date +%T` );
    print "[$time][- / -] Working on sample $name.\n";
    foreach my $step ( @steps )
    {
        chomp ( $time = `date +%T` );
        my ($clean_step) = $step;
        $clean_step =~ s/ -/\n                  -/g;
        print "[$time][$nom/$#steps] Running this step: \n", " "x18, "$clean_step\n";
       #system ( $step );
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

         OPTION    Default                                                 Description
         BIN       /opt/var_calling                                        Absolute location of the Picard Tools and GATK jar files
         REF_DIR   /safer/genomes/Homo_sapiens/UCSC/hg19                   Absolute location of the reference directory
         BT2       /safer/genomes/Homo_sapiens/Sequence/BowtieIndex/hg19   Absolute location of the Bowtie2 index
         THREADS   64                                                      Number of threads to use in parallelizable modules
         MEMORY    64                                                      Amount of memory, in gigabytes, to use
         READS_DIR N/A                                                     Absolute location of the reads that are going to be used
     
USAGE
}
