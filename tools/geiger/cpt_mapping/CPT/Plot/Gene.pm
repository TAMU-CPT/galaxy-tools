package CPT::Plot::Gene;
use Moose;
use strict;
use warnings;

# ABSTRACT: Stupid representation of a gene. Does not handle joined genes
has 'start'  => (is => 'rw', isa => 'Int' );
has 'end'    => (is => 'rw', isa => 'Int' );
has 'tag'    => (is => 'rw', isa => 'Str' );
has 'label'  => (is => 'rw', isa => 'Str' );
has 'strand' => (is => 'rw', isa => 'Str' );
has 'color'  => (is => 'rw', isa => 'Any' );


sub getLocations {
	my ($self) = @_;
	return [ $self->start(), $self->end() ];
}

no Moose;
1;
