#!/usr/local/bin/perl -w
use strict;

chomp ( my @folders = `ls -d /work/unmc_ngs/acornish/parallel/mckeithan/*` );

foreach my $folder ( @folders) 
{
    my ($job_name ) = $folder =~ /.+\/(.+)/;
    open ( OUT, ">tmp.pbs" );
    print OUT <<END;
#!/bin/sh
#PBS -N $job_name
#PBS -l select=40
#PBS -l mem=64GB
#PBS -l walltime=72:00:00
#PBS -e $job_name.stderr
#PBS -o $job_name.stdout
cd $folder
gorrammit.pl /work/unmc_ngs/acornish/conf.$job_name
END
    close OUT;
    system ( "qsub tmp.pbs" );
}
