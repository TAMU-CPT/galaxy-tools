package CPT::Plot::Colours;


# ABSTRACT: Color transformation library.

our %artemis_colours = (
	0  => 'rgb(255,255,255)',
	1  => 'rgb(100,100,100)',
	2  => 'rgb(255,0,0)',
	3  => 'rgb(0,255,0)',
	4  => 'rgb(0,0,255)',
	5  => 'rgb(0,255,255)',
	6  => 'rgb(255,0,255)',
	7  => 'rgb(255,255,0)',
	8  => 'rgb(152,251,152)',
	9  => 'rgb(135,206,250)',
	10 => 'rgb(255,165,0)',
	11 => 'rgb(200,150,100)',
	12 => 'rgb(255,200,200)',
	13 => 'rgb(170,170,170)',
	14 => 'rgb(0,0,0)',
	15 => 'rgb(255,63,63)',
	16 => 'rgb(255,127,127)',
	17 => 'rgb(255,191,191)',
);

sub new {
	my $class = shift;
	my $self  = {@_};
	bless $self, $class;
	return $self;
}

sub getColour {
	my ( $self, $string ) = @_;
	if ($string) {
		my $colour_result;
		if ( $string =~ qr/^\s*(\d+)\s*$/ ) {
			$colour_result = $artemis_colours{$1};
		}
		elsif ( $string =~ qr/^\s*(\d+)\s+(\d+)\s+(\d+)\s*$/ ) {
			$colour_result = "rgb($1,$2,$3)";
		}
		else {
			warn "Bad Colour Specfication";
			return undef;
		}
		return $colour_result;
	}
	else {
		return undef;
	}
}

1;

