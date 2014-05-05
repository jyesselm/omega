package VectorFunctions;

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
use Math::Trig;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BaseObject ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw(IsAngleWithinMax IsAngleWithinMin IsPlanarAngleWithinMax IsWithinMaxDistance IsWithinMinDistance CalculateAngleBetweenVectors CalculatePlanarAngle CalculateVectorDifference Magnitude CalculateVectorDifferenceBetweenAtoms)],
											 
										 'vars' => [ qw () ],
										
										 'func' => [ qw(IsAngleWithinMax IsAngleWithinMin IsPlanarAngleWithinMax IsWithinMaxDistance IsWithinMinDistance CalculateAngleBetweenVectors CalculatePlanarAngle CalculateVectorDifference Magnitude CalculateVectorDifferenceBetweenAtoms)]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';


sub IsAngleWithinMax {
	
	my ($Angle, $MaxAngle) = @_;
	
	return $MaxAngle >= $Angle ? 1 : 0;
	
}



sub IsAngleWithinMaxOld {
	
	my ($Info, $Angle, $Requirements, $MaxAngle) = @_;
	
	unless(defined $Angle) {

	  my @RequirementArray = map { $Info->{$_} } @$Requirements;

    $Angle = CalculateAngleBetweenVectors(@RequirementArray);

  }

	return $Angle if $MaxAngle >= $Angle;

	return -1;
 
}

sub IsAngleWithinMin {
	
	my ($Angle, $MinAngle) = @_;
	
	return $MinAngle <= $Angle ? 1 : 0;
	
}

sub IsPlanarAngleWithinMax {
	
	my ($PlaneAngle, $MaxPlaneAngle) = @_;
	
	return $MaxPlaneAngle >= $PlaneAngle ? 1 : 0;
	
}

sub IsWithinMaxDistance {
	
	my ($Distance, $MaxDistance) = @_;
	
	return $MaxDistance >= $Distance ? 1 : 0;
		
}

sub IsWithinMinDistance {
	
	my ($Distance, $MinDistance) = @_;
	
	return $MinDistance <= $Distance ? 1 : 0;
	
}

sub IsWithinMaxDistanceOld {
	
	my ($Info, $Distance, $Requirements, $MaxDistance) = @_;
	
	unless(defined $Distance) {
	
	  my @RequirementArray = map { $Info->{$_} } @$Requirements;
	
	  $Distance = Magnitude($RequirementArray[0]);
	
  }
				
	return $Distance if $MaxDistance > $Distance;
	
	return -1;
	
}

=head2 CalculateAngleBetweenVectors

Usage: CalculateAngleBetweenVectors($Vector1, $Vector2);

Arguments:
  Both Vectors are cartesian vectors in the form [X,Y,Z]

Synopsis:
  Calculates the Angle given two vectors AB and AC such that 

B
 \
  A -- C 

giving you the angle BAC, uses the trig relationship BAC = ACos(AB*AC/||AB||*||AC||) 

=cut

sub CalculateAngleBetweenVectors {
	
	my ($Vector1, $Vector2) = @_;
		
  my ($Mag1, $Mag2) = (Magnitude($Vector1), Magnitude($Vector2));
  
  my $DotProduct = DotProduct($Vector1, $Vector2);

	my $Angle = ArcCosine($DotProduct /($Mag1*$Mag2))*(180/3.14159265);
	
	return $Angle;
	
}

=head2 CalculatePlanarAngle

Usage: CalculatePlanarAngle($Vector1, $Vector2, $Vector3);

Arguments:
  All Vectors are cartesian vectors in the form [X,Y,Z]

Synopsis:
  Calculate The 3D plane that Vector1 and Vector2 share, used with O=C-R groups, then determines the angle in which Vector3 is out of the plane,
this is achieved by calculating the $Vector1 X $Vector2 to get the normal vector and then taking the arcsin of Vec1*Vec2 / ||Vec1||*||Vec2||;

=cut

sub CalculatePlanarAngle {
	
	my ($Vector1, $Vector2, $Vector3) = @_;
	
	my $NormalVector = CrossProduct($Vector1, $Vector2);
		
	my $DotProduct = DotProduct($NormalVector, $Vector3);

  return -1 if Magnitude($NormalVector) == 0;
		
	my $PlaneAngle = asin($DotProduct / (Magnitude($NormalVector)*Magnitude($Vector3)))*(180/3.14159265);
	
	#If $PlanarAngle 
	
	if($PlaneAngle < 0) {
		
	  my $NegNormalVector = [map { -$_ } @$NormalVector];
		
		my $DotProduct = DotProduct($NegNormalVector, $Vector3);

		my $NegAngle = asin($DotProduct / (Magnitude($NegNormalVector)*Magnitude($Vector3)))*(180/3.14159265);
			
    $PlaneAngle = $NegAngle;
		
	}
	
	return $PlaneAngle;
	
}

#Math Functions!

sub DotProduct {
	
  my ($Vector1, $Vector2) = @_;

  return $Vector1->[0]*$Vector2->[0] + $Vector1->[1]*$Vector2->[1] + $Vector1->[2]*$Vector2->[2];

}

sub CalculateVectorDifference {
	
	my ($Vector1, $Vector2) = @_;
	
	return [map { $Vector1->[$_] - $Vector2->[$_] } (0 .. scalar(@{$Vector1}) - 1)];
	
}

sub CalculateVectorDifferenceBetweenAtoms {
	
	my ($Atom1, $Atom2) = @_;
	
	my $Vector1 = $Atom1->getCartesianCoordinates;
	my $Vector2 = $Atom2->getCartesianCoordinates;
	
	return [map { $Vector1->[$_] - $Vector2->[$_] } (0 .. scalar(@{$Vector1}) - 1)];
	
}


sub Magnitude {
	
	my $Vector = shift;
	
	my $Total = 0;
	
	foreach my $Element (@$Vector) {
		
		$Total += $Element**2;
		
	}
	return sqrt($Total);
}

sub ArcCosine { 
	my $Num = shift;
	return atan2( sqrt(1 - $Num * $Num), $Num ); 
}

sub CrossProduct {
	
	my ($V1, $V2) = @_;
	
	return [$V1->[1]*$V2->[2] - $V1->[2]*$V2->[1],
	 				$V1->[2]*$V2->[0] - $V1->[0]*$V2->[2], 
					$V1->[0]*$V2->[1] - $V1->[1]*$V2->[0] ];
	
}

1;