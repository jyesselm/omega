#!/usr/bin/perl

use strict;
use warnings;

my @pdbs = grep { $_ =~ /\.pdb$/} <*>;

my $ParameterFile;
my $type1;
my $type2;
my $RemoveWaters;

#Process cmd arguments need to put in hash!
foreach my $i (0 .. @ARGV-1) {

	next if $ARGV[$i] !~ /^-/;

	if(lc($ARGV[$i]) =~ /^-cutoffs/) {

		$ParameterFile = $ARGV[$i+1];

	}

	elsif(lc($ARGV[$i]) =~ /^-type1/) {

		$type1 = $ARGV[$i+1];

	}
	
	elsif(lc($ARGV[$i]) =~ /^-type2/) {

		$type2 = $ARGV[$i+1];

	}

	elsif(lc($ARGV[$i]) =~ /-removewaters/) {

		$RemoveWaters = 1;

	}

}


foreach my $pdb (@pdbs) {
	 
   my $command = "/Users/skullnite/Dropbox/Public/Work/hbonds/FinalProgram/test.pl -cutoffs $ParameterFile -pdb $pdb";

   if($RemoveWaters) { $command .= " -removewaters"}

   my $name = substr($pdb,0,-4);
	
   system($command);
   system("/Users/skullnite/Dropbox/Public/Work/hbonds/FinalProgram/CompareResultFiles.pl $type1 $type2 > $name.out");	

   exit(0);
	
}


