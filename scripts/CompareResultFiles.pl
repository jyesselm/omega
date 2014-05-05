#!/usr/bin/perl -w

=head1 NAME


=head1 SYNOPSIS



=head1 DESCRIPTION


=head2 EXPORT


=head1 AUTHOR

Joseph Yesselman, E<lt>jyesselm@umich.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Joseph Yesselman

This library is free software; you can redistribute it and/or modifycd
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=head1 FUNCTIONS

=cut

use strict;
use warnings;
use Carp;
use Math::Trig;

my $HeavyAtom = $ARGV[0];
my $Hydrogen  = $ARGV[1];
my $Filter;

my @dat_files = grep { $_ =~ /\.dat$/ } <*>;

my %HydrogenData;
my %HeavyAtomData;

my @hydrogen_included_files;
my @heavyatom_incuded_files;

my %categoriesinfile;
my %referenceinfo;

foreach my $file (@dat_files) {
	
	if($file =~ /$Hydrogen/) {
		
		open(FILE, $file);
		
		my @FileContents = <FILE>;
		
		close(FILE);
		
		shift @FileContents; pop @FileContents;
		
		my $categorylist = shift @FileContents;
		
		my @categories = split /\s+/, $categorylist;
		
		my @referencedata = splice(@categories, 0, 8);
		
		$referenceinfo{$file} = [$referencedata[0], $referencedata[1], $referencedata[3], $referencedata[4],$referencedata[5],$referencedata[6],$referencedata[7]];
		
		$categoriesinfile{$file} = \@categories;
		
		foreach my $Line (@FileContents) {
			
			my @spl = split /\s+/, $Line;
			
			my $key = $spl[0] . "-" . $spl[1] . "-" . $spl[3] . "-" . $spl[4]  . "-" . $spl[5] . "-" . $spl[6] . "-" . $spl[7];
			
			splice(@spl, 0, 8);
						
			$HydrogenData{$key} = [$file, @spl];
						
		}
		
		push @hydrogen_included_files, $file;
		
	}
	
	elsif($file =~ /$HeavyAtom/) {

		open(FILE, $file);

		my @FileContents = <FILE>;

		close(FILE);

		shift @FileContents; pop @FileContents;

		my $categorylist = shift @FileContents;

		my @categories = split /\s+/, $categorylist;

		splice(@categories, 0, 7);

		$categoriesinfile{$file} = \@categories;

		foreach my $Line (@FileContents) {

			my @spl = split /\s+/, $Line;
			
			my $key = $spl[0] . "-" . $spl[1] ."-" . $spl[2] . "-" . $spl[3]  . "-" . $spl[4] . "-" .	 $spl[5] . "-" . $spl[6];

			splice(@spl, 0, 7);
		
			$HeavyAtomData{$key} = [$file, @spl];

		}
		
		push @heavyatom_incuded_files, $file;

	}
	
} 

my @file_pairs;

#print scalar(@hydrogen_included_files) . "\n";

foreach my $ha_file (@heavyatom_incuded_files) {
		
	my @ha_members = split /\_/, substr($ha_file, 0, -4);
	
	foreach my $h_file (@hydrogen_included_files) {
		
		my $same = 0;
		
		my @h_members = split /\_/, substr($h_file, 0, -4);
		
		foreach my $ha_member (@ha_members) {
			
			next if $ha_member eq $HeavyAtom;
			
			foreach my $h_member (@h_members) {
				
				next if $h_member eq $Hydrogen;
				
				if($h_member eq $ha_member) {
					
					$same = 1; last;
					
				}
				
				last if $same;
			} 
			
			last if $same;
			
		}
				
		if($same) {
			
			push @file_pairs, [$ha_file, $h_file];
			
		}
		
	}
	
}

#foreach my $file_pair (@file_pairs) {
	
#	print join(" ", @$file_pair)  ."\n";
	
#}

#exit 1;


my %seen;

my @uniquekeys = grep { ! $seen{$_} ++ }(keys %HydrogenData, keys %HeavyAtomData);

foreach my $file_pair (@file_pairs) {

	my $ha_cats = $categoriesinfile{$file_pair->[0]};
	my $h_cats = $categoriesinfile{$file_pair->[1]};
	my $reference_info = $referenceinfo{$file_pair->[1]};
	
	print join(" ", @$reference_info) . " " . join(" ", @$h_cats) . " " . join(" ", @$ha_cats) . "\n";

	foreach my $key (@uniquekeys) {
	
		my $isCorrectData = 0;
		my $isDefined = 0;
	
		my $ha_data = $HeavyAtomData{$key};
		my $h_data = $HydrogenData{$key};
		
		my @split_key = split /\-/, $key;
		
		my $String = join(" ", @split_key) . " ";
		
		if(defined $h_data) {
			
			$isDefined++;
			
			$isCorrectData++ if $h_data->[0] eq $file_pair->[1];
			
			my @data = map { $_ } @$h_data;
			
			shift @data;
			
			$String .= join(" ", @data) . " ";
						
		}
		
		else {
			
			$String .= join(" ", map { "N/A"} (0 .. @$h_cats-1)) . " ";
			
		}
		
		if(defined $ha_data) {
			
			$isDefined++;

			$isCorrectData++ if $ha_data->[0] eq $file_pair->[0];
			
			my @data = map { $_ } @$ha_data;

			shift @data;

			$String .= join(" ", @data) . " ";
			
		}
		
		else {
			
			$String .= join(" ", map { "N/A"} (0 .. @$ha_cats-1)) . " ";		
			
		}
		
		if($isDefined == $isCorrectData && $isDefined > 0) {
			
			print $String . "\n";
			
		}
		
	
	}

}
