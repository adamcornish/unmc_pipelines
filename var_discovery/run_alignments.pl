#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opts;
getopts ( '1:2:n:b:p:x:s:G:R:D:', \%opts );

my $R1       = $opts{1};
my $R2       = $opts{2};
my $name     = $opts{n};
my $bin      = $opts{b};
my $threads  = $opts{p};
my $bt2_idx  = $opts{x};
my $step     = $opts{s};
my $gatk     = $opts{G};
my $ref      = $opts{R};
my $dbsnp    = $opts{D};
my $JAVA_pre = "java -Xms8G -Xmx8G -XX:MaxPermSize=1G -jar";
my $GATK_pre = "$JAVA_pre $gatk -T";

chomp ( my $time = `date +%T` );
print "[$time] Working on sample $name.\n";

my @aln = (
            "bowtie2 --very-sensitive-local -x $bt2_idx -p $threads -1 $R1 -2 $R2 -S tmp/$name.sam",
            "samtools view -bS tmp/$name.sam -o tmp/$name.bam",
            "$JAVA_pre $bin/SortSam.jar INPUT=tmp/$name.bam OUTPUT=tmp/$name.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT",
            "$JAVA_pre $bin/MarkDuplicates.jar I=tmp/$name.sorted.bam O=tmp/$name.dup_removed.bam REMOVE_DUPLICATES=true M=tmp/$name.mark_dups_metrics_file",
            "$JAVA_pre $bin/AddOrReplaceReadGroups.jar I=tmp/$name.sorted.bam O=tmp/$name.fixed_RG.bam SO=coordinate RGID=$name RGLB=$name RGPL=illumina RGPU=$name RGSM=$name VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
            "$GATK_pre RealignerTargetCreator -R $ref -I tmp/$name.fixed_RG.bam -known $dbsnp -o tmp/$name.indel.intervals", # -o must have .intervals extension
            "$GATK_pre IndelRealigner -R $ref -I tmp/$name.fixed_RG.bam -known $dbsnp -o tmp/$name.indels_realigned.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals tmp/$name.indel.intervals",
            "$GATK_pre BaseRecalibrator -R $ref -knownSites $dbsnp -I tmp/$name.indels_realigned.bam -o tmp/$name.recal_data.grp", # can't do -nt flag until v. 2.2
            "$GATK_pre PrintReads -R $ref -BQSR tmp/$name.recal_data.grp -I tmp/$name.indels_realigned.bam -o tmp/$name.BQSR.bam",
            "echo > $name.done",
          );

for ( my $i = $step; $i < @aln; $i++ )
{
    my $current_step = $aln[$i];
    my ($clean_step) = $current_step;
    my $nom          = sprintf ( "%02d/%02d", $i, $#aln );
    $clean_step      =~ s/ -/\n                  -/g if length ($clean_step) > 256;
    chomp ( $time = `date +%T` );
    print "[$time][$nom]  Running this step: \n\n", " "x18, "$clean_step\n\n";
    system ( $current_step );
}
