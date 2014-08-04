use strict;
use warnings;
use Data::Dumper;
use Carp;
use Test::More tests => 98;
use Test::Exception;

require_ok('CPT::Util');
my $libCPT = CPT::Util->new();

my %colour_model = %{ $libCPT->JSONYAMLopts( file => 'data/colouring.yaml' ) };

my %regex_containers;
my %custom_regexes;
foreach my $key ( keys %colour_model ) {
	foreach ( @{ $colour_model{$key}{'members'} } ) {
		my $regex = qr/\b$_\b/i;
		$regex_containers{$regex} = $colour_model{$key}{'title'};

		#$regex_containers{$regex} = $colour_model{$key}{'color'};
	}
	foreach ( keys %{ $colour_model{$key}{'custom'} } ) {
		$custom_regexes{$_} = {
			'title'  => $colour_model{$key}{'title'},
			'colour' => $colour_model{$key}{'color'},
			'is'     => $colour_model{$key}{'custom'}{$_}{'is'},
			'isnot'  => $colour_model{$key}{'custom'}{$_}{'isnot'},
		};
	}
}

open( COLOUR_TESTS, '<', 't/color_tests.tsv' );
while (<COLOUR_TESTS>) {
	chomp $_;
	my ( $string, $should_match ) = split( /\t/, $_ );
	test_string( $string, $should_match );
}
close(COLOUR_TESTS);

sub test_string {
	my ( $string, $expected ) = @_;
	my $result = method($string);
	ok( $result eq $expected,
		"Tested [$string]  Got [$result] Expected [$expected]" );
}

sub method {
	my ($string) = @_;
	my $matched_result = '';
	foreach my $regex ( keys %regex_containers ) {
		if ( $string =~ $regex ) {
			$matched_result = $regex_containers{$regex};
		}
	}
	foreach my $custom ( keys %custom_regexes ) {

		# Do we possibly want to have another look at this one?
		my $care = 0;
		foreach ( @{ $custom_regexes{$custom}{'is'} } ) {
			if ( $string =~ /\b$_\b/i ) {
				$care = 1;
			}
		}
		foreach ( @{ $custom_regexes{$custom}{'isnot'} } ) {
			if ( $string =~ /\b$_\b/i ) {
				$care = 0;
			}
		}
		if ($care) {
			my $is_okay         = 1;
			my $ok_to_overwrite = 1
			  ; # Fixes a strange bug. If we have a match, and then we have a subpart of that match which hits a custom element, we want to make sure that it'll ONLY overwrite if we didn't specifically exclude items like the one we hit.
			foreach ( @{ $custom_regexes{$custom}{'is'} } ) {
				if ( $string !~ /\b$_\b/i ) {
					$is_okay = 0;
				}
			}
			foreach ( @{ $custom_regexes{$custom}{'isnot'} } ) {
				if ( $string =~ /\b$_\b/i ) {
					$ok_to_overwrite = 0;
					$is_okay         = 0;
				}
			}
			if ( $is_okay && $ok_to_overwrite ) {
				$matched_result =
				  $custom_regexes{$custom}{'title'};
			}
		}
	}
	return $matched_result;
}
