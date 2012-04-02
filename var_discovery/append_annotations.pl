#!/usr/bin/perl
use warnings;
use strict;

my @db   = "";
my @data = `cat $ARGV[0] | grep -v '^#'`;
my $curr = "";
open OUT, ">ddl";

foreach my $datum ( @data )
{
    my ($chr)   = $datum =~ /^(\S+)/;
    if ( $curr ne $chr )
    {
        print "Reading in $chr.\n";
        @db = `cat /data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/dbNSFP/$chr | grep -vP '^#'`;
        print "Sorting the db.\n";
        @db = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, /\S+\s\S+\s\S+\s\S+\s\S+\s\S+\s(\S+)/s] } @db;
        print "Done reading and sorting the db.\n";
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
                my $out = "\t$col[16]\t$col[17]\t$col[18]\t$col[19]\t$col[20]\t$col[21]\t$col[22]\t$col[23]\t$col[24]\t$col[25]";
                print OUT $out;
            }
        }
        print OUT "\n";
    }
}

close OUT;

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
