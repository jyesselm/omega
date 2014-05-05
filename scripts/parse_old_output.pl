#!/usr/bin/perl

use strict;
use warnings;

open(FILE, $ARGV[0]);

my @FileContents = <FILE>;

close(FILE);

foreach my $Line (@FileContents) {
	
	my @spl = split /\s+/, $Line;
	
	next if $spl[9] == 0;
	
	print $Line;
	
}
