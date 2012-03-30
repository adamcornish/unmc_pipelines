#!/usr/bin/perl
use warnings;
use strict;

my @list = `cat $ARGV[0] | grep -v '^#'`;
my $curr = "";
my @db;
open OUT, ">wat";
foreach my $item ( @list )
{
    my ($chrom) = $item =~ /^(\S+)/;
    my @split  = split "\t", $item;
    if ( $curr ne $chrom )
    {
        @db = `cat /data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/dbNSFP/chr$chrom`;
        $curr = $chrom;
    }
    if ( $item =~ /NON_SYN/ )
    {
        my ($chr)  = $split[0] =~ /chr(\S+)/;
        my $search = "$chr\\s\\S+\\s$split[3]\\s$split[4]\\s\\S\\s\\S\\s$split[1]";
        print "search  : $search\n";
    }
    else
    {
        my ($eff)   = $item =~ /EFF=(.+?)\s/;
        my @effects = split ",", $eff;
        foreach my $effect ( @effects )
        {
            my ($loc, $parens) = $effect =~ /(.+?)\((.+?)\)/;
            my ($impact,$type,$codon,$aa,$gene,$rna,$coding,$acc,$exon) = split "|", $parens;
            print "$split[0]\t$split[1]\t$split[2]\t$split[3]\t$split[4]\t$loc\t$impact\t$type\t$codon\t$aa\t$gene\t$acc\t$exon\t";
            <>;
            print OUT "$split[0]\t$split[1]\t$split[2]\t$split[3]\t$split[4]\t$loc\t$impact\t$type\t$codon\t$aa\t$gene\t$acc\t$exon\t";
        }
    }
}

sub bin_search
{
    my ( $elem, $list ) = @_;
    my $max = $#$list;
    my $min = 0;

    while ( $max >= $min ) {
        my $mid    = int( ( $max + $min ) / 2 );
        my ($item) = $list[$mid] =~ /\S+\s\S+\s\S\s\S\s\S\s\S\s(\S+)/;
        if    ( $list->[$mid] < $elem ) { $min = $mid + 1; }
        elsif ( $list->[$mid] > $elem ) { $max = $mid - 1; }
        else                              { return $mid; }
    }
    return;
}
