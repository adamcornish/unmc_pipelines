#!/usr/bin/perl
use warnings;
use strict;

print "Usage: perl $0 directory.\n\tSimply provide the directory that has the *finished.snvs* vcf files in them.\n\n";
my $dir = ( $ARGV[0] ) ? $ARGV[0] : "./";
chomp ( my @files = `ls *finished.snvs*` );

foreach my $file ( @files )
{
    print "Working on $file\n";
    my @db   = "";
    my @data = `grep -v '^#' $file`;
    my $curr = "";
    my ($name ) = $file =~ /(.+?)\./;
    open OUT, ">$name.snvs.txt";
    print OUT "Chromosome\tPosition\tdbSNP ID\tRef\tAlt\tMutation\tImpact\tType\tCodon Change\tAA Change\tGene\tAccession ID\tExon\tPhyloP\tSIFT\tPolyphen2\tLRT\tMutationTaster\n";

    foreach my $datum ( @data )
    {
        my ($chr)   = $datum =~ /^(\S+)/;
        if ( $curr ne $chr )
        {
            print "\tReading in $chr.\n";
            @db = `cat /data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/dbNSFP/$chr | grep -vP '^#'`;
            $curr = $chr;
        }
        my @split   = split "\t", $datum;
        my ($eff)   = $datum =~ /EFF=(.+?)\s/;
        my @effects = split ",", $eff;
        foreach my $effect ( @effects )
        {
            my ($loc, $parens) = $effect =~ /(.+?)\((.+?)\)/;
            my ($impact,$type,$codon,$aa,$gene,$rna,$coding,$acc,$exon) = split /\|/, $parens;
            print OUT "$split[0]\t$split[1]\t$split[2]\t$split[3]\t$split[4]\t$loc\t$impact\t$type\t$codon\t$aa\t$gene\t$acc\t$exon";
            if ( $effect =~ /NON_SYN/)
            {
                my ($ref_AA, $alt_AA) = $aa =~ /([A-Z])\d+([A-Z])/;
                my ($ref_nuc, $alt_nuc, $index) = ($split[3], $split[4], $split[1]);
                my $result = bin_search ( \@split, \@db, $ref_AA, $alt_AA );
                if ( $result )
                {
                    my @col = split "\t", $result;
                    my $out = "\t$col[7]\t$col[8]\t$col[9]\t$col[10]\t$col[12]";
                    print OUT $out;
                }
            }
            print OUT "\n";
        }
    }

    close OUT;
}

sub bin_search
{
    my ( $arr, $db, $ref_AA, $alt_AA ) = @_;
    my $max = $#$db;
    my $min = 0;

    while ( $max >= $min ) 
    {
        my $mid = int( ( $max + $min ) / 2 );
        my @s   = split "\t", $$db[$mid];
        if    ( $s[6] < $$arr[1] ) { $min = $mid + 1; }
        elsif ( $s[6] > $$arr[1] ) { $max = $mid - 1; }
        else
        { 
            my $regex = "$alt_AA\\s$$arr[1]";
            for ( my $i = ($mid - 4); $i < ($mid + 4); $i++ )
            {
                return $$db[$i] if $$db[$i] =~ /$regex/;
                return "" if ($i == ($mid + 3));
            }
        }
    }
}
