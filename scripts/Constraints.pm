package Constraints;

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

#use 5.010000;
use strict;
use warnings;
use Carp;

use VectorFunctions ':func';


our %ConstraintTypeHash = (
	
	maxdistance    => \&IsWithinMaxDistance,
	mindistance    => \&IsWithinMinDistance,
	maxangle  	 	 => \&IsAngleWithinMax,
	minangle  	 	 => \&IsAngleWithinMin,
	maxplaneangle => \&IsPlanarAngleWithinMax
	
);

our %ConstraintValues = (
	
	maxdistance    => 1,
	mindistance 	 => 1,
	maxangle 			 => 2,
	minangle 			 => 2,
	maxplaneangle  => 3,
	
);

our %ConstraintFuntionArgumentExpection = (

	CalculateVectorDifferenceBetweenAtoms => ['Atom', 'Atom'],

);

our %test = (


	
);


sub New { 

  my ($Class, $Types, $Constraints) = @_;

  my $Self = {
	
	  Types 			=> $Types,
	  Constraints => $Constraints,
	  Conditions  => [], 
	  DataHash    => { $Types->[0]->getName => undef, $Types->[1]->getName => undef},
	  Geometrics   => {}, 
	 
  };

	foreach my $i (1 .. 10) {
		
	  $Self->{DataHash}->{"Result$i"} = undef;
		
	}

  bless $Self, $Class;

  return $Self;	
	
}

sub Initiate { 

  my $Self = shift;

  my $Constraints = $Self->{Constraints};

  my @Conditions = $Self->ParseFileConstraintsIntoSuccessConditions($Constraints);
		
	$Self->{Conditions} = \@Conditions;
		
}


#Organize Constraints in conditions for success
sub ParseFileConstraintsIntoSuccessConditions {
	
	my ($Self, $Constraints) = @_;
	
	my @Conditions = ([]);

  my @InAllConditions;

  foreach my $Constraint (@$Constraints) {

	  my $Declaration = $Constraint->[0];

	  my @SplitOverPeriod = split /\./, $Declaration;

	  my $AtomInclusions = $SplitOverPeriod[0];
	  my $ConstraintType;
	  my $ConditionNum = 0;  

		#Specifies what condition it is apart of example: A-B.0.MaxDistance = 3.00, part of condition 0
	  if(scalar(@SplitOverPeriod) == 3) {

			$ConditionNum = $SplitOverPeriod[1];
			$ConstraintType = $SplitOverPeriod[2];	

			my $ProcessingSubroutine = LinkConstraintTypeToSubroutine($ConstraintType);

		  my $ThingsToCalculate = $Self->DeconstructInclusionsIntoRequirements($AtomInclusions, $ConstraintType);

			if(!$Conditions[$ConditionNum]) {

				$Conditions[$ConditionNum] = [];

			}

			#print $ConditionNum . " " . $ProcessingSubroutine . " " . $Constraint->[2] . "\n";

			push @{$Conditions[$ConditionNum]}, [$ThingsToCalculate, $ProcessingSubroutine, $Constraint->[2], "$AtomInclusions.$ConstraintType"];


		}

		#Does not have a condition number such as: A-B.MaxDistance = 3.00, this sets this constraint for all conditions
		else {

			$ConstraintType = $SplitOverPeriod[1];

			my $ProcessingSubroutine = LinkConstraintTypeToSubroutine($ConstraintType);

			my $ThingsToCalculate = $Self->DeconstructInclusionsIntoRequirements($AtomInclusions, $ConstraintType);

			push @InAllConditions, [$ThingsToCalculate, $ProcessingSubroutine, $Constraint->[2], "$AtomInclusions.$ConstraintType"];

		}

	}

		
	foreach my $i (0 .. @Conditions-1) {

		my %SeenConstraint = map { $_->[3] => 1} @{$Conditions[$i]};
		
		foreach my $Constraint (@InAllConditions) {
		
			next if $SeenConstraint{$Constraint->[3]};
			
			push @{$Conditions[$i]}, $Constraint;
	
		}
		
		my @SortedConstraints = @{$Conditions[$i]}; #sort { $ConstraintValues{lc($a->[3])} <=> $ConstraintValues{lc($b->[3])} } @{$Conditions[$i]};
		
		$Conditions[$i] = \@SortedConstraints;

	}
		
	return @Conditions;
	
}

sub DoAtomsStatisifyConstraints {
	
	my ($Self, $Atoms) = @_;
	
	my $Conditions = $Self->{Conditions};
	
	my $Success = 0;
	
  $Self->{Results} = { };
		
	foreach my $Atom (@$Atoms) {
		
		$Self->{DataHash}->{$Atom->getType} = $Atom;
		
	}

	my $ConditionCount = 0;

	foreach my $Condition (@$Conditions) {
		
		my $ConstraintSatisifedCount = 0;

		foreach my $Constraint (@$Condition) {

			my $ProcessingSubrountine = $Constraint->[1];

			my $ResultCount = 1;
			
			my $Failed = 0;

			foreach my $ThingToCalculate (@{$Constraint->[0]}) {

				my $Sub = $ThingToCalculate->[0];

				my @DeferencedArgs = map { ${$ThingToCalculate->[$_]} } (1 .. @$ThingToCalculate-1);

				$Self->{DataHash}->{"Result$ResultCount"} = getCalculatedInfo($Sub, @DeferencedArgs);

				if($Self->{DataHash}->{"Result$ResultCount"} == -1) { $Failed=1; last; }

				#print $ResultCount . " " . $Self->{DataHash}->{"Result$ResultCount"} . "\n";
			
				$ResultCount++;

			}
			
			if($Failed == 1) {
				
				$ConstraintSatisifedCount++; next;
			
			}

			$ResultCount--;

			#print $Self->{DataHash}->{"Result$ResultCount"} . " " . $Constraint->[2] . "\n";

			my $IsConstraintSatisifed = &$ProcessingSubrountine($Self->{DataHash}->{"Result$ResultCount"}, $Constraint->[2]);
			
			$Self->{Results}->{ $test{$Constraint->[0]->[-1]} } = $Self->{DataHash}->{"Result$ResultCount"};

			if($IsConstraintSatisifed) {

				$ConstraintSatisifedCount++;

			}

		}
		
		if($ConstraintSatisifedCount == @$Condition) {
			
			$Success = 1; last;
			
		}

	}
			
	return $Success;
		
}


sub DeconstructInclusionsIntoRequirements {
	
	my ($Self, $AtomInclusions, $ConstraintType) = @_;
	
	my @ThingsToCalculate;
	
	my @AtomNames = split /\-/, $AtomInclusions;
	
	my $Types = $Self->{Types};
		
	my $CalculateCount = 1;
	
	my $AtomNameCount = -1;
	
	my @AtomPositionValues;
		
	foreach my $AtomName (@AtomNames) {
		
		my $IsAType = 0;
					
		foreach my $Type (@$Types) {
			
			if($AtomName eq $Type->getName) { 
				
				$IsAType = 1; 
								
				#print scalar(@A) . " " . $AtomName . "\n";
				
				push @AtomPositionValues, \$Self->{DataHash}->{$AtomName};
			
			}
			
		}
		
		next if $IsAType; 
				
		if($AtomName =~ /^R(\w+)/) {
			
			push @AtomPositionValues, \$Self->{DataHash}->{"Result$CalculateCount"};
			
			push @ThingsToCalculate, [\&getRAtoms, \$Self->{DataHash}->{$1}];
						
			$CalculateCount++;
			
			next;
			
		}
		
		elsif($AtomName =~ /^H(\w+)/) {
			
			push @AtomPositionValues, \$Self->{DataHash}->{"Result$CalculateCount"};

			push @ThingsToCalculate, [\&getHAtoms, \$Self->{DataHash}->{$1}];

			$CalculateCount++;

			next;
			
		}
		
	}
		
  push @ThingsToCalculate, [\&CalculateVectorDifferenceBetweenAtoms, $AtomPositionValues[1], $AtomPositionValues[0]]; #Returns to Result1
	
	my $FinalCalculate = undef;
		
	if($ConstraintType =~ /Distance/) {
				
		$FinalCalculate = [\&Magnitude, \$Self->{DataHash}->{"Result$CalculateCount"}];
						
		push @ThingsToCalculate, $FinalCalculate;		
		
	}
	
	elsif($ConstraintType =~ /Angle/ && $ConstraintType !~ /Plane/) {
		
		$CalculateCount++;
		
		my $Current = $CalculateCount;
		my $Last = $CalculateCount-1;
		
		$FinalCalculate = [\&CalculateAngleBetweenVectors, \$Self->{DataHash}->{"Result$Current"}, \$Self->{DataHash}->{"Result$Last"}];
				
		push @ThingsToCalculate, [\&CalculateVectorDifferenceBetweenAtoms, $AtomPositionValues[1], $AtomPositionValues[2]];
		push @ThingsToCalculate, $FinalCalculate;
		
	}
	
	elsif($ConstraintType =~ /Plane/) {
		
		push @ThingsToCalculate, [\&getRAtoms, $AtomPositionValues[1]]; #Returns to Result2 (is OC atom)
		push @ThingsToCalculate, [\&getRAtoms, \$Self->{DataHash}->{"Result2"}]; #Returns to Result3 (is OCR atom)
		push @ThingsToCalculate, [\&CalculateVectorDifferenceBetweenAtoms, \$Self->{DataHash}->{"Result2"}, $AtomPositionValues[1]]; #Returns to Result4 (OC_O_Vector)
		push @ThingsToCalculate, [\&CalculateVectorDifferenceBetweenAtoms, \$Self->{DataHash}->{"Result3"}, \$Self->{DataHash}->{"Result2"}]; #Returns to Result5 (OCR_OC_Vector)
		$FinalCalculate = [\&CalculatePlanarAngle, \$Self->{DataHash}->{Result4}, \$Self->{DataHash}->{Result5}, \$Self->{DataHash}->{Result1}]; #Returns to Result6, Plane Angle?
		
		push @ThingsToCalculate, $FinalCalculate;
		
	}
	
	
	#print $FinalCalculate . "\n";
		
	my ($Property) = (lc($ConstraintType) =~ /^(?:min|max)(\S+)/);
			
	$test{$FinalCalculate} = join("-", @AtomNames) . ".$Property";
	
	$Self->{Geometrics}->{join("-", @AtomNames) . ".$Property"} = 1;
		
	return \@ThingsToCalculate;
}


sub LinkConstraintTypeToSubroutine {
	
	my $ConstraintType = shift;
	
	unless(defined $ConstraintTypeHash{lc($ConstraintType)}) {
		
		croak("$ConstraintType is a not a valid Constraint");
	  
		
	}
	
	return $ConstraintTypeHash{lc($ConstraintType)};
	
	
}


sub getRAtoms {
	
	my $Atom = shift;
	
	my $BondedAtoms = $Atom->getBondedAtoms;
	
	my @RAtoms = grep { $_->getElement ne "H" && $_->getElement ne "D" } @$BondedAtoms;
	
	return -1 if scalar(@RAtoms) == 0;
	
	return $RAtoms[0];
	
}

sub getHAtoms {
	
	my $Atom = shift;
	
	my $BondedAtoms = $Atom->getBondedAtoms;
	
	my @HAtoms = grep { $_->getElement eq "H" || $_->getElement eq "D" } @$BondedAtoms;
	
	return $HAtoms[0];
	
}


{

  my $HashedResults;

	sub getCalculatedInfo {
	
		my ($Sub, @Args) = @_;
		
		my $Key = join(" ", ($Sub,@Args));
		
		if(defined $HashedResults->{$Key}) {
			
			return $HashedResults->{$Key};
			
		}
		
		my $Result = &$Sub(@Args);
		
		$HashedResults->{$Key} = $Result;
		
		return $Result;
	
	}

}

1;
