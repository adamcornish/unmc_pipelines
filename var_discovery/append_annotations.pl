#!/usr/bin/perl
use warnings;
use strict;

my @list = `cat $ARGV[0] | grep -v '^#'`;
my $curr = 0;
my @db;
foreach my $item ( @list )
{
    my ($chrom) = $item =~ /^(\S+)/;
    $curr = $chrom unless $curr;
    if ( $item =~ /NON_SYN/ )
    {
        my @split  = split "\t", $item;
        my ($chr)  = $split[0] =~ /chr(\S+)/;
        my $chunk  = `grep $split[1] /data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/dbNSFP/chr$chr`;
        my $search = "$chr\\s\\S+\\s$split[3]\\s$split[4]\\s\\S\\s\\S\\s$split[1]";
        print "search  : $search\n";
        my ($line) = $chunk =~ /($search.+?)\n/;
        print "line is : $line\n";
    }
}
