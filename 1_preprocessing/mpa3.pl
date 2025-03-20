#!/usr/bin/perl

foreach $ar (@ARGV)
{

@field = split (/_/, $ar);
$name = $ar;
$name =~ s/_1.fastq//;

$r1 = $ar; 
$r2 = $ar; 
$r2 =~ s/1.fastq/2.fastq/; 

system("metaphlan $r1,$r2 --bowtie2db bowtie_mpa3_v30 --bowtie2out $name.bowtie2.bz2 --nproc 4 --input_type fastq  -t rel_ab --index mpa_v30_CHOCOPhlAn_201901 --mpa3 -o ../profiled_MPA3/profiled_$name.txt");
}
