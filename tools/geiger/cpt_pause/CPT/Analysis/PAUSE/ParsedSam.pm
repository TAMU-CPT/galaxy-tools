package CPT::Analysis::PAUSE::ParsedSam;

# ABSTRACT: Library for use in PAUSE analysis
use strict;
use warnings;
use Moose;
use SVG;

has 'coverage_density'               => ( is => 'rw', isa => 'ArrayRef' );
has 'read_starts'                    => ( is => 'rw', isa => 'ArrayRef' );
has 'read_ends'                      => ( is => 'rw', isa => 'ArrayRef' );
has 'max'                            => ( is => 'rw', isa => 'Int' );
has 'stats_start_max'                => ( is => 'rw', isa => 'Num' );
has 'stats_end_max'                  => ( is => 'rw', isa => 'Num' );
has 'stats_start_mean'               => ( is => 'rw', isa => 'Num' );
has 'stats_end_mean'                 => ( is => 'rw', isa => 'Num' );
has 'stats_start_standard_deviation' => ( is => 'rw', isa => 'Num' );
has 'stats_end_standard_deviation'   => ( is => 'rw', isa => 'Num' );

no Moose;
1;
