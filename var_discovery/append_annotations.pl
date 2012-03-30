#!/usr/bin/perl
use warnings;
use strict;

my @list = `cat $ARGV[0] | grep -v '^#'`;
my $curr = "";
open OUT, ">y";
open DDL, ">z";
my ($result, $prev_search) = ("","");
foreach my $item ( @list )
{
    my ($chrom) = $item =~ /^(\S+)/;
    my @split   = split "\t", $item;
    my ($eff)   = $item =~ /EFF=(.+?)\s/;
    my @effects = split ",", $eff;
    my $db      = "/data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/dbNSFP/$chrom";

    foreach my $effect ( @effects )
    {
        my ($loc, $parens) = $effect =~ /(.+?)\((.+?)\)/;
        my ($impact,$type,$codon,$aa,$gene,$rna,$coding,$acc,$exon) = split /\|/, $parens;
        if ( $loc =~ /NON_SYNON/ )
        {
            my ($ref_AA, $new_AA) = $aa =~ /([A-Z])\d+([A-Z])/;
            my $search = "$split[3]\\s$split[4]\\s$ref_AA\\s$new_AA\\s$split[1]";
            if ( $search ne $prev_search )
            {
                $result = `grep -P '$search' $db`;
                $prev_search = $search;
            }
            print "eff: $effect\n\t\t$result\n";
            print DDL "$result\n";
        }
        else
        {
            print OUT "$split[0]\t$split[1]\t$split[2]\t$split[3]\t$split[4]\t$loc\t$impact\t$type\t$codon\t$aa\t$gene\t$acc\t$exon\t\n";
        }
        print "prev_search: $prev_search\n";
    }
}
close OUT;
close DDL;
