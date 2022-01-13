#!/usr/bin/perl

use strict;
use warnings;

my $gemma = $ARGV[0];
my $gene  = $ARGV[1];
my $out   = $ARGV[2];


## make a hash for the gene names
my %gene_hash;

## fill the hash using the cluster id
open(GENE,'<',$gene) or die $!;
while(<GENE>){
    chomp;
    my @linestuff = split(" ",$_);
    my ($cluster_id,$chr,$intron_start,$intron_end,$gene_id) = @linestuff[0..4];
    $gene_hash{$cluster_id} = $gene_id;
    # print "$cluster_id\t$gene_id\n";
}
close(GENE);


## open the GEMMA results, iterate over
open(GEMMA,"gunzip -c $gemma |") or die $!;
open(OUT,'>', $out) or die $!;
while(<GEMMA>){
    chomp;

    # extract relevant fields
    my @linestuff = split("\t", $_);
    my ($chr,$intron_start,$intron_end,$cluster_id) = split(":",$linestuff[0]);
    $cluster_id =~ s/-9$//;
    print OUT join("\t",@linestuff)."\t".$gene_hash{$cluster_id}."\n" if exists $gene_hash{$cluster_id};
}
close(GEMMA);
close(OUT);
