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

#It is important to install the MATCH toolset to run this script see: 
use lib $ENV{'MATCH'} . "/lib"; 
use lib "/Users/skullnite/Dropbox/Public/Work/hbonds/FinalProgram";

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use Type;
use MoleculeFileHandler;

use Constraints;

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;

my $ParameterFile;
my $PDB;
my $RemoveWaters = 0;
my $DEBUG = 0;
my $CreateDataFiles = 1;

#Process cmd arguments need to put in hash!
foreach my $i (0 .. @ARGV-1) {

	next if $ARGV[$i] !~ /^-/;

	if(lc($ARGV[$i]) =~ /^-cutoffs/) {

		$ParameterFile = $ARGV[$i+1];

	}

	elsif(lc($ARGV[$i]) =~ /^-pdb/) {

		$PDB = $ARGV[$i+1];

	}

	elsif(lc($ARGV[$i]) =~ /-removewaters/) {

		$RemoveWaters = 1;

	}

}

#Load cutoff file, -cutoffs at command line
my $ConstraintObjects = ConstructCutoffObjectFromCutoffFile($ParameterFile);

croak("No PDB file was specified please use -pdb filename.pdb commandline argument\n") unless defined $PDB;

#Build molecule from pdb -pdb at command
my $MoleculeFileHandler = MoleculeFileHandler->New($PDB);

my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

my @Atoms = map { @{$_->getAtoms} } @$Structures;

#To Handle Deutrium Atoms which is annoying to pattern match out since its element is D vs H, just changing them to Hs
#Really shouldn't make a difference
foreach my $Atom (@Atoms) {
	
	if($Atom->getElement eq "D") {
		
		$Atom->setElement("H");
		$Atom->setState("H.1");
		
	}
	
}

my %ElementsofTypes;
my %Types;

#Look at each Constraint object and collect all the types, and save all of the 
foreach my $Constraint (@$ConstraintObjects) {
	
	my $Types = $Constraint->{Types};

  foreach my $Type (@$Types) {
			
	 	$Type->getLookUpTable->SetAtomLevelsUsingNodeLevels();
	
		my $Atoms = $Type->getAtoms;
		
		my $HeadAtom;
		
		foreach my $Atom (@$Atoms) {
			
			if($Atom->getLevel == 0) { $HeadAtom = $Atom }
						
		}
		
		$ElementsofTypes{$HeadAtom->getElement} = 1;
		$Types{$Type->getName} = $Type;
	
	}
	
} 

my %AtomHash;

my @UniqueTypes = values %Types;

foreach my $Atom (@Atoms) {

  $Atom->setType("");

  next unless $ElementsofTypes{$Atom->getElement};

  if($RemoveWaters) {
	
		next if $Atom->getMolecule->getName =~ /HOH|DOD/;
	
	}

  my $AtomLookUpTable = LookUpTable->New;

  $AtomLookUpTable->Initiate($Atom);

  $AtomLookUpTable->Update();

	foreach my $Type (@UniqueTypes) {
								
		next unless $Type->getLookUpTable->AreAllNodesSharedBetween([$AtomLookUpTable]);

		if(! defined $AtomHash{$Type->getName}) { $AtomHash{$Type->getName} = [] }

		push @{$AtomHash{$Type->getName}}, $Atom;
		
		#print $Atom->getName . " " . $Atom->getMolecule->getName . " "  . $Type->getName . "\n";
		
	}

}

#exit 1;

my $OutputHandler = *STDOUT;

my $pdb_name = substr($PDB,0,-4);

foreach my $Constraint (@$ConstraintObjects) {
	
	my $Types = $Constraint->{Types};
	
	my @Type1Atoms = @{$AtomHash{$Types->[0]->getName}};
	my @Type2Atoms = @{$AtomHash{$Types->[1]->getName}};
	
	my ($T1, $T2) = ($Types->[0]->getName, $Types->[1]->getName);
	
	if($CreateDataFiles) {
		
		open(FILE, ">$T1\_$T2.dat");
		
		$OutputHandler = *FILE;
		
	}
	
	print $OutputHandler "Number of Type " . $Types->[0]->getName . " Atoms: " . scalar(@Type1Atoms) . " Number of Type " . $Types->[1]->getName . " Atoms: ". scalar(@Type2Atoms) . "\n";
		
	my $HitCount = 0;

	my %Seen;
	
	my $ConstraintGeometrics = $Constraint->{Geometrics};
	
	my @sorted_keys = map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [OrderGeometricType($_), $_] } keys %$ConstraintGeometrics;
	
	my @other_info = ("$T1-AtomName", "$T1-ResName", "$T1-ResNum");
	my @other_info2 = ("$T2-AtomName", "$T2-ResName", "$T2-ResNum");
	
	my @other_subs = (sub { return $_[0]->getName}, sub { return $_[0]->getMolecule->getName}, sub { return $_[0]->getMolecule->getNum });
	my @other_subs2 = (sub { return $_[0]->getName}, sub { return $_[0]->getMolecule->getName}, sub { return $_[0]->getMolecule->getNum });
	
	
	if($Type1Atoms[0]->getElement =~ /^H|D/) {
	
	  @other_info = ("$T1-HeavyAtom",@other_info);
		@other_subs = (sub { return $_[0]->getBondedAtoms->[0]->getName }, @other_subs);
		
	}
	
	elsif($Type1Atoms[1]->getElement =~ /^H|D/) {
		
	 @other_info2 = ("$T2-HeavyAtom",@other_info2);
   @other_subs2 = (sub { return $_[0]->getBondedAtoms->[0]->getName }, @other_subs2);
		
	}
	
	@other_info = ("PDBFile", @other_info);
	@other_subs = (sub { return $pdb_name }, @other_subs);
	
	print $OutputHandler join(" ", @other_info, @other_info2, @sorted_keys) . "\n";
	
	my @other_info_lengths = map { length($_) } @other_info;
	my @other_info2_lengths = map { length($_) } @other_info2;
	my @sorted_keys_lengths = map { length($_) } @sorted_keys;
	
	
	print scalar(@Type1Atoms) . " " . scalar(@Type2Atoms) . "\n";
		
	my $type_count = 0;
		
	foreach my $T1Atom (@Type1Atoms) {
		
		$type_count++;
		
		print $type_count . "\n";
		
		$T1Atom->setType($Types->[0]->getName);
		
		foreach my $T2Atom (@Type2Atoms) {
		
			next if $T1Atom eq $T2Atom;
			
			next if $T1Atom->IsBondedTo($T2Atom);
			
			$T2Atom->setType($Types->[1]->getName);
									
			my $Success = $Constraint->DoAtomsStatisifyConstraints([$T1Atom, $T2Atom]);
					
			next unless $Success;
			
			$HitCount++;
			
			$Seen{$T1Atom} = 1;
						
			my $Results = $Constraint->{Results};
			
			my $count = 0;
					
			foreach my $sub (@other_subs) {
				
				my $length = $other_info_lengths[$count];
				
				print $OutputHandler sprintf("%-$length" . "s", &$sub($T1Atom)) . " ";
				
				$count++;
				
			}
			
			$count = 0;
			
			foreach my $sub (@other_subs2) {

				my $length = $other_info2_lengths[$count];

				print $OutputHandler sprintf("%-$length" . "s", &$sub($T2Atom)) . " ";

				$count++;

			}
			
			$count = 0;
			
			foreach my $key (@sorted_keys) {

				my $length = $sorted_keys_lengths[$count];

				my $value = "N/A";

				if(defined $Results->{$key}) {
					
					$value = sprintf("%.2f", $Results->{$key});
					
				}

				print $OutputHandler sprintf("%-$length" . "s", $value) . " ";

				$count++;

			}
			
			print $OutputHandler "\n";
			
		
		}	
			
	}
	
	print $OutputHandler "Stats: Total Type1 Atoms = " . scalar(@Type1Atoms) . " Meet Constraints = " . $HitCount . " Percent = " . ($HitCount/scalar(@Type1Atoms)*100) . "%\n";
	
	if($CreateDataFiles) {
		
		close(FILE);
		
	}
	
}


=head2 ConstructCutoffObjectFromCutoffFile

Usage: 
  my $CutoffObjects = ConstructCutoffObjectFromCutoffFile($CutoffFilePath);

Arguments:
  $CutoffFilePath: the absolute path of the cutoff file path 

Returns:
  reference to an array of constraint objects from Constraints.pm

Member Of:
  Standalone

Synopsis:
  Builds constraint objects (Constraints.pm) from the type and constraints in cutoff file.

=cut

sub ConstructCutoffObjectFromCutoffFile {
	
	my $CutoffFilePath = shift;
	
	croak("no cutoff file specified please specify with -cutoffs commandline argument\n") unless defined $CutoffFilePath;
	
	croak("$CutoffFilePath file does not exist please double check to make sure you have the correct path\n") unless -e $CutoffFilePath;
 	
	open(FILE, $CutoffFilePath);
	
	my @FileContents = <FILE>;
	
	close(FILE);
	
	my @ConstraintObjects;
	
	my @TypeDeclartions = grep { $_ =~ /^[^\-]+\s*\=\s*\S+/ && $_ !~ /#/ } @FileContents;
	
	croak ("There are less then 2 specified types in $CutoffFilePath, please see examples for correct syntax\n") if scalar(@TypeDeclartions) < 2;
		
	my @Types;
	
	#Build Type objects from PerlChemistry module to use in chemoinformatic typing engine 
	foreach my $TypeDeclartion (@TypeDeclartions) {
		
		my @split = split /\s+/,$TypeDeclartion;
		
		#split[0] = Type Name
		#split[2] = Type String, a description of the chemical informatics	
		my $Type = Type->New($split[0], 1, $split[2]);
		
		$Type->Initiate(); 
				
		push @Types, $Type;
				
	}
	
	my @TypeNames = map { $_->getName } @Types;
 	
  my @ContraintDeclarations = map { my @spl = split/\s+/, $_; [@spl] } grep { $_ =~ /^\S+\-\S+\s*\=\s*\S+/ && $_ !~ /#/ } @FileContents;

	croak ("There are no constraints in $CutoffFilePath, please see examples for correct syntax\n") if scalar(@ContraintDeclarations) < 1;

  #Perform exhustive pairing of Types declared in file any constraint declaration that contains both will be included in the current constraint object
  foreach my $i (0 .. @TypeNames-1) {

		croak ("Types cannot start with the letter \"R\", R is a reserved letter to designate the use of an atom bound to a type\n") if $TypeNames[$i] =~ /^R/; 

    foreach my $j ($i+1 .. @TypeNames-1) {
		
	    my @Constraints;
		
			#Does the Constraint Declaration contain both types, if so add it to the list 
			foreach my $ContraintDeclaration (@ContraintDeclarations) {
		
			  push @Constraints, $ContraintDeclaration if isConstraintDeclartionPartofConstraintObject([$Types[$i], $Types[$j]], $ContraintDeclaration->[0]);
										
			}
			
			next unless @Constraints;
			
			#for DEBUGing purposes print all constraints of this Constraint object
			
			print $TypeNames[$i] . " " . $TypeNames[$j] . "\n" if $DEBUG;
			
			foreach my $Constraint (@Constraints) {
				
				print $Constraint->[0] . "\n" if $DEBUG;
				
			}
			
			#Build Contraints Object			
	    my $ConstraintObject = Constraints->New([$Types[$i], $Types[$j]], \@Constraints);
	
			$ConstraintObject->Initiate();
	
			push @ConstraintObjects, $ConstraintObject;

		}
	
  }

  return \@ConstraintObjects;

}

=head2 isConstraintDeclartionPartofConstraintObject

Usage: 
  my $Success = isConstraintDeclartionPartofConstraintObject($Types, $ConstraintDeclartion)

Arguments:
  $Types: the names of the two current types, example: ["A", "D"]
  $ConstraintDeclartion: the first element in a constraint declartion from the cutoff file, example D-A.MaxDistance

Member Of:
	subfunction of ConstructCutoffObjectFromCutoffFile

Synopsis:
  Checks to see if the current constraint declartion is part of the two current types this is done using string comparisions checking to see if the
type names are within the members of the constraint declaration

=cut

sub isConstraintDeclartionPartofConstraintObject {
	
	my ($Types, $ConstraintDeclartion) = @_;
	
	#Example: get RD-D-A substring of RD-D-A.MaxAngle
	my @SplitOverPeriod = split /\./, $ConstraintDeclartion;
	
	#Example: breaks RD-D-A into (RD,D,A) array
	my @ConstraintDeclartionMembers = split /\-/, $SplitOverPeriod[0];
	
	my $OverlapCount = 0;
		
	foreach my $Member (@ConstraintDeclartionMembers) {
		
		my $Success = 0;
		
		foreach my $Type (@$Types) {
						
			my $TypeName = $Type->getName;
			
			#To be related to a type it must either equal one of the members or be 
			if($TypeName eq $Member || $Member =~ /^(R|H)$TypeName$/ ) {
				
				$Success = 1; last;
				
			}
			
		}
		
		$OverlapCount++ if $Success;
		
	}
		
	return $OverlapCount == @ConstraintDeclartionMembers ? 1 : 0;
	
}

sub OrderGeometricType {
	
	my $Geometric = shift;
	
	if(lc($Geometric) =~ /plane/) {
		
		return 3;
		
	}
	
	elsif(lc($Geometric) =~ /distance/) {
		
		return 1;
		
	}
	
	else { 
		
	  return 2;
		
	}
	
}



