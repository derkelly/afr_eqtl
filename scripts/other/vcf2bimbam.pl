#!/usr/bin/perl

use strict;
use warnings;

my $vcf = $ARGV[0];
my $out = $ARGV[1];

open(VCF,"gunzip -c $vcf |") or die $!;
open(OUT,'>', $out) or die $!;
while(<VCF>){
    chomp;
    next if /^#/;

    # extract relevant fields
    my @linestuff = split("\t", $_);
    my ($chr,$pos,$id,$ref,$alt) = @linestuff[0..4];

    my @doses    = ();
    my $sum_alls = 0; my $foo = 0;

    my $numstuff = scalar(@linestuff) - 9;

    # extract the dosage information for all genotypes
    foreach my $geno (splice(@linestuff, 9)){
	my ($all1,$all2,$dose) = $geno =~ /(\d)\|(\d)\:(\d\.?\d*)/;
	push @doses, $dose;
	$sum_alls    += ($all1 + $all2);
    }

    my $dose_line;
    
    # if the alternative allele has a frequency above 50%, then the
    # ref is the minor allele and the doses must be subtracted from 2
    if ($sum_alls <= $numstuff){
	$dose_line = join("\t", @doses);
    	print OUT "$chr\t$pos\t$id\t$alt\t$ref\t$dose_line\n";
    }
    else{
	@doses = map { 2 - $_ } @doses;
	$dose_line = join("\t", @doses);
 	print OUT "$chr\t$pos\t$id\t$ref\t$alt\t$dose_line\n";
    }
}
close(VCF);
close(OUT);
