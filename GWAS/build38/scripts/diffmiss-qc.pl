#!/usr/bin/perl

use strict;

open IN, '<', $ARGV[0].".tMiss.missing"
        or die "Cannot open .missing file (".$ARGV[0].".tMiss.missing): $!\n";
open OUT, '>', $ARGV[1];
while(<IN>){
	s/^\s+//;
	my @fields = split /\s+/, $_;
	unless($fields[0] eq 'CHR'){
		if($fields[4] < 0.000001){
			print OUT "$fields[1]\n";
		}
	}
}
