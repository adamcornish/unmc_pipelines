#!/usr/bin/perl
use warnings;
use strict;

######## Start Variables ########

my $bin       = "/data/tools/tmp_bin";
my $gatk      = "$bin/GenomeAnalysisTK.jar";
my $dbsnp     = "/data/hg19/dbsnp134.vcf";
my $indels    = "/data/hg19/dbsnp134_indels.vcf";
my $ref       = "/data/hg19/ref.fa";
my $exome_bed = "/data/hg19/refSeq_genes.bed";
my $name      = "smith_exome";
my $read_1    = "smith_R1.fastq";
my $read_2    = "smith_R2.fastq";
my $log       = "run.log";
my $threads   = "24";
my $memory    = "42";
my $i         = 1;
my $total     = 18;
chomp ( my $time   = `date +%T` );

######## End Variables ########

print "[$time][",$i++,"/$total] Aligning read 1 using bwa aln.\n";
system ( "$bin/bwa aln -q 1 -t $threads $ref $read_1 > R1.sai 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Aligning read 2 using bwa aln.\n";
system ( "$bin/bwa aln -q 1 -t $threads $ref $read_2 > R2.sai 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Aligning both reads using bwa sampe.\n";
system ( "$bin/bwa sampe $ref R1.sai R2.sai $read_1 $read_2 > a.sam 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Converting a.sam to a.bam using samtools view -bS.\n";
system ( "samtools view -bS a.sam -o a.bam >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Sorting a.bam using PicardTools SortSam.jar.\n";
system ( "java -Xmx${memory}g -jar $bin/SortSam.jar INPUT=a.bam OUTPUT=b.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Marking Duplicates using PicardTools MarkDuplicates.jar.\n";
system ( "java -Xmx${memory}g -jar $bin/MarkDuplicates.jar INPUT=b.bam OUTPUT=b_dup.bam M=metric VALIDATION_STRINGENCY=SILENT >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Indexing b.bam using samtools index.\n";
system ( "samtools index b_dup.bam >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Creating realignment targets - needed for indel realignment - using GATK RealignerTargetCreator.\n";
system ( "java -Xmx${memory}g -jar $gatk -T RealignerTargetCreator -R $ref -I b_dup.bam -known $indels -o indel_realigner.intervals >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Realigning around indels using GATK IndelRealigner.\n";
system ( "java -Xmx${memory}g -jar $gatk -T IndelRealigner -R $ref -I b_dup.bam -known $indels -o c.bam -targetIntervals indel_realigner.intervals >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Sorting c.bam using PicardTools SortSam.jar.\n";
system ( "java -Xmx${memory}g -jar $bin/SortSam.jar INPUT=c.bam OUTPUT=d.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Indexing d.bam using samtools index.\n";
system ( "samtools index d.bam >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Fixing Read Group in d.bam using PicardTools AddOrReplaceReadGroups.jar.\n";
system ( "java -Xmx${memory}g -jar $bin/AddOrReplaceReadGroups.jar I=d.bam O=e.bam SO=coordinate RGID=$name RGLB=$name RGPL=illumina RGPU=$name RGSM=$name VALIDATION_STRINGENCY=SILENT >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Indexing e.bam using samtools index.\n";
system ( "samtools index e.bam >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Counting covariates using GATK CountCovariates.\n";
system ( "java -Xmx${memory}g -jar $gatk -T CountCovariates -nt $threads -R $ref --knownSites $dbsnp -I e.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -dP illumina -recalFile recal.csv >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Recalibrating read quality scores using GATK TableRecalibration.\n";
system ( "java -Xmx${memory}g -jar $gatk -T TableRecalibration -R $ref -I e.bam --out f.bam -recalFile recal.csv >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Indexing f.bam using samtools index.\n";
system ( "samtools index f.bam >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Performing variant calling using GATK UnifiedGenotyper - SNVs.\n";
system ( "java -Xmx8g -jar $gatk -T UnifiedGenotyper -nt $threads -R $ref -I f.bam -o raw_variants.snvs.vcf              -mbq 20 -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 -D $dbsnp  -L $exome_bed >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Performing variant calling using GATK UnifiedGenotyper - INDELs.\n";
system ( "java -Xmx8g -jar $gatk -T UnifiedGenotyper -nt $threads -R $ref -I f.bam -o raw_variants.indels.vcf -glm INDEL -mbq 20 -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 -D $indels -L $exome_bed >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Filtering variants using GATK VariantFiltration - SNVs.\n";
system ( "java -Xmx${memory}g -jar $gatk -T VariantFiltration -R $ref -V raw_variants.snvs.vcf -o filtered_variants.snvs.vcf -filter 'DP <= 20 || HRun > 8 || QD < 5.0' -filterName 'standard_filters' -filter 'MQ0 >= 4 && ((MQ0 / (1.0  * DP)) > 0.1 )' -filterName 'hard_to_validate' >> $log.stdout 2>> $log.stderr" );
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] Filtering variants using GATK VariantFiltration - INDELs.\n";
system ( "java -Xmx${memory}g -jar $gatk -T VariantFiltration -R $ref -V raw_variants.indels.vcf -o filtered_variants.indels.vcf -filter 'DP <= 20 || HRun > 8 || QD < 5.0' -filterName 'standard_filters' -filter 'MQ0 >= 4 && ((MQ0 / (1.0  * DP)) > 0.1 )' -filterName 'hard_to_validate' >> $log.stdout 2>> $log.stderr" );
=cut
chomp ( $time = `date +%T` );
print "[$time][",$i++,"/$total] .\n";
system ( "java -jar /data/tools/GATK/dist/GenomeAnalysisTK.jar -T VariantRecalibrator -R /data/hg19_ref/ref.fa -input filtered_variants.vcf -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -tranchesFile tranches_file.txt -mode BOTH -Rscript /usr/local/bin/ -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 /data/hg19_ref/dbsnp134.vcf -recalFile recal_file.txt -rscriptFile r_script.R -resource:omni,known=false,training=true,truth=false,prior=12.0 /data/hg19_ref/1000G_omni2.5.hg19.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/hg19_ref/hapmap_3.3.hg19.vcf \ >> $log.stdout 2>> $log.stderr" );
