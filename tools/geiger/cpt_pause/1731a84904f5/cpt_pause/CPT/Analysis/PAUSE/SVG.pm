package CPT::Analysis::PAUSE::SVG;

# ABSTRACT: Library for use in PAUSE analysis
use strict;
use warnings;
use Moose;
use Data::Dumper;
use List::MoreUtils qw(each_array);
use SVG;

has 'svg'                => ( is => 'rw' );
has 'width'              => ( is => 'rw', isa => 'Int' );
has 'height'             => ( is => 'rw', isa => 'Int' );
has 'vertical_offset'    => ( is => 'rw', isa => 'Int' );
has 'start_end_max_num'  => ( is => 'rw', isa => 'Int' );
has 'num_rows'           => ( is => 'rw', isa => 'Int' );
has 'row_size'           => ( is => 'rw', isa => 'Int' );
has 'row_width'          => ( is => 'rw', isa => 'Int' );
has 'x_border'           => ( is => 'rw', isa => 'Int' );
has 'y_border'           => ( is => 'rw', isa => 'Int' );
has 'line_height'        => ( is => 'rw', isa => 'Int' );
has 'inter_line_spacing' => ( is => 'rw', isa => 'Int' );
has 'max'                => ( is => 'rw', isa => 'Int' );
has 'fasta_id'           => ( is => 'rw', isa => 'Any' );

sub setup {
	my ($self) = @_;
	$self->svg(
		SVG->new(
			width  => $self->width(),
			height => $self->height(),
		)
	);
}

sub add_header {
	my ($self, @refs) = @_;
	
	$self->plot_title( 'Plot of ' . $self->fasta_id() );

	my $i = 0;
	foreach(@refs){
		my @subrefs = @{$_};
		foreach(@subrefs){
			$i++;
			my %d = %{$_};
			$self->plot_key($d{name},$d{line},$d{fill}, $i);
		}
	}
	$self->vertical_offset( $self->vertical_offset() - ($i - 1) * 20);
}

my $global_pline_idx = 0;
sub plot_track {
	my ( $self, $points_ref, $stroke, $fill, $id) = @_;

	$global_pline_idx++;
	$self->svg()->polyline(
		%{$points_ref},
		id    => 'pline_' . $id . '-' . $global_pline_idx,
		style => {
			'fill-opacity' => .5,
			'stroke'       => $stroke,
			'fill'         => $fill,
		}
	  )

}

sub make_scale {
	my ( $self, $i, $start, $stop ) = @_;

	# Left axis label, must be rotated
	my $tmp_x = $self->fix_x_value(-30);       #$self->x_border()-30;
	my $tmp_y = $self->fix_y_value($i) + 50;
	$self->svg()->text(
		id            => 'left_side_label_row_' . $i,
		x             => $tmp_x - 20,
		y             => $tmp_y - 20,
		'font-family' => 'Helvetica, sans-serif',
		'transform'   => sprintf( 'rotate(-90 %s %s)', $tmp_x, $tmp_y ),
	)->cdata('Start/End Hit Count Scale');

	# Right axis label, must be rotated
	$tmp_x = $self->fix_x_value( $self->row_width() + 60 );
	$tmp_y = $self->fix_y_value($i);
	$self->svg()->text(
		id            => 'right_side_label_row_' . $i,
		x             => $tmp_x - 20,
		y             => $tmp_y,
		'font-family' => 'Helvetica, sans-serif',
		'transform'   => sprintf( 'rotate(-90 %s %s)', $tmp_x, $tmp_y ),
	)->cdata('Coverage Density');

	# Horizontal increments
	for ( my $k = -4 ; $k <= 4 ; $k++ ) {

		# Left side label
		my $y_position = $k / 4 * $self->line_height();
		$self->svg()->text(
			id => sprintf( 'label_left_side_row_%s_%s', $i, $k ),
			x  => $self->fix_x_value(-30),
			y             => $self->fix_y_value( $i, $y_position ),
			'font-family' => 'Helvetica, sans-serif',
		)->cdata( int( $self->start_end_max_num() * abs( $k / 4 ) ) );

		# Right side label
		$self->svg()->text(
			id => sprintf( 'label_right_side_row_%s_%s', $i, $k ),
			x => $self->fix_x_value( $self->row_width() + 10 ),
			y => $self->fix_y_value( $i, $y_position ),
			'font-family' => 'Helvetica, sans-serif',
		)->cdata( int( $self->max() * abs( $k / 4 ) ) );

		# Vertical lines
		$self->svg()->line(
			x1 => $self->fix_x_value( $self->row_width() ),
			x2 => $self->fix_x_value(0),
			y1 => $self->fix_y_value( $i, $y_position ),
			y2 => $self->fix_y_value( $i, $y_position ),
			id => sprintf( 'vertical_increment_row_%s_%s', $i, $k ),
			opacity        => .25,
			stroke         => 'rgb(0,0,0)',
			'stroke-width' => '2',
		);
	}

	# Vertical Increments
	my $number_of_increments = 10;
	for (
		my $k = 0 ;
		$k <= $self->row_width();
		$k += ( $self->row_width() / $number_of_increments )
	  )
	{
		# We get % of way across (k/num_inc) and we multiply by the width value, to get % of width which we adjust with start to get correct value
		my $b =   ($k / $self->row_width()) * ($stop-$start) + $start;
		my $kb = $b/1000;
		$self->svg()->text(
			id =>
			  sprintf( 'vertical_line_label_row_%s_%s', $i, $k ),
			x => $self->fix_x_value($k),
			y => $self->fix_y_value( $i, $self->line_height() + 20 ),
			'font-family' => 'Helvetica, sans-serif',
		)->cdata( ( $kb ) . ' kb' );

		$self->svg()->line(
			x1 => $self->fix_x_value($k),
			x2 => $self->fix_x_value($k),
			y1 => $self->fix_y_value($i, -$self->line_height() ),
			y2 => $self->fix_y_value($i, $self->line_height() ),
			id => sprintf( 'vertical_line_row_%s_%s', $i, $k ),
			opacity        => .5,
			stroke         => 'rgb(0,0,0)',
			'stroke-width' => '1',
		);
	}
}

sub plot_title {
	my ( $self, $string ) = @_;
	$self->svg()->text(
		id            => 'label_plot_title',
		x             => $self->x_border(),
		y             => 50 + $self->vertical_offset(),
		'font-family' => 'Helvetica, sans-serif',
		'font-size'   => '150%',
	)->cdata($string);
	$self->vertical_offset( $self->vertical_offset() + 25 );
}

sub plot_key {
	my ( $self, $text, $stroke, $colour, $i) = @_;

	$self->svg()->rectangle(
		x              => $self->x_border(),
		y              => 50 + $self->vertical_offset() - 15,
		width          => 15,
		height         => 15,
		id             => 'label_key_example'. $i,
		'fill-opacity' => .5,
		'stroke'       => $stroke,
		'fill'         => $colour,
	);
	$self->svg()->text(
		id            => 'label_key_string'. $i,
		x             => $self->x_border() + 20,
		y             => 50 + $self->vertical_offset(),
		'font-family' => 'Helvetica, sans-serif',
	)->cdata($text);
	$self->vertical_offset( $self->vertical_offset + 20 );
}

sub xmlify {
	my ($self) = @_;
	return $self->svg()->xmlify();
}

sub x_values_for_range_scaled {
	my ( $self, $start, $end, $pieces ) = @_;
	my @vals = ($start);
	my $by   = ( $end - $start ) / $pieces;
	## For all values from the x_border to xborder+row_width, add a value of row_width split/row_size (i.e., how far for EACH INDIVIDUAL value)
	for ( my $i = $start ; $i < $end ; $i += $by ) {
		push( @vals, $i );
	}
	push( @vals, $end );
	return @vals;
}

sub fix_x_values {
	my ( $self, @values ) = @_;
	return map { $self->fix_x_value($_) } @values;
}

sub fix_x_value {
	my ( $self, $val ) = @_;
	return $val + $self->x_border();
}

sub fix_y_value {
	my ( $self, $i, $val ) = @_;
	return (
		$self->vertical_offset() + $val - $self->line_height() + (
			( ( 2 + $i ) * $self->line_height() ) +
			  ( $i * $self->inter_line_spacing() ) +
			  $self->y_border()
		)
	);
}

sub fix_all_y_values {
	my ( $self, $i, @arrays_to_fix, ) = @_;
	for ( my $j = 0 ; $j < scalar @arrays_to_fix ; $j++ ) {

		# For each array in postive_y (AoA)
		#
		# we cast to array, then we map this, then we have this
		# in an anonymous array which means we can just
		# replace. This is probably not as efficient as looking
		# at every value directly and doing "in place"
		# replacement, but I don't know how that would be
		# written here...
		$arrays_to_fix[$j] =
		  [ map { $self->fix_y_value( $i, $_ ) }
			  @{ $arrays_to_fix[$j] } ];
	}
	return @arrays_to_fix;
}

sub copy_data {
	my ( $self, $start, $stop, $data_to_ref, $data_from_ref, $max) = @_;
	# Copy data from the original array to the new one, transforming out
	# the subset of interest. This is done across an AoA
	my @data_to   = @{$data_to_ref};
	my @data_from = @{$data_from_ref};
	for ( my $k = 0 ; $k < scalar @data_from ; $k++ ) {
		foreach ( my $j = $start ; $j < $stop ; $j++ )
		{    #1 to 10_000 in the genome
			if ( defined ${ $data_from[$k] }[$j] ) {
				push(
					@{ $data_to[$k] },
					-(
						$self->line_height() *
						  ${ $data_from[$k] }[$j] /
						  $max
					)
				);
			}
			else {
				push( @{ $data_to[$k] }, 0 );
			}
		}
	}
	return @data_to;
}

sub plot_individual_row {
	my ( $self, $start, $stop, $i, $regular_ref, $rescale_ref ) = @_;

	my @regular = map { ${$_}{data} } @{$regular_ref};
	my @rescale = map { ${$_}{data} } @{$rescale_ref};
	my @regular_y;
	my @rescale_y;
	## Ensure we duplicate the number of arrays.
	foreach (@regular) {
		push( @regular_y, [] );
	}
	foreach (@rescale) {
		push( @rescale_y, [] );
	}

	# Determine bounds of row
	$self->push_all( \@regular_y, 0 );
	$self->push_all( \@rescale_y, 0 );

	@regular_y = $self->copy_data( $start , $stop , \@regular_y , \@regular , $self->start_end_max_num());
	@rescale_y = $self->copy_data( $start , $stop , \@rescale_y , \@rescale , $self->max());
	#print @rescale_y;

	# Set up our X values
	my @x_values         = $self->fix_x_values( $self->x_values_for_range_scaled(0, $self->row_width(), ( $stop - $start )));
	#my @x_values_rescale = $self->fix_x_values( $self->x_values_for_range_scaled(0, $self->row_width(), ( $stop - $start )));

	$self->push_all( \@regular_y, 0);
	$self->push_all( \@rescale_y, 0);
	# Fix the ys
	@regular_y = $self->fix_all_y_values( $i, @regular_y );
	@rescale_y = $self->fix_all_y_values( $i, @rescale_y );
	# Prepare our styling
	my @regular_line = map { ${$_}{line} } @{$regular_ref};
	my @rescale_line = map { ${$_}{line} } @{$rescale_ref};
	my @regular_fill = map { ${$_}{fill} } @{$regular_ref};
	my @rescale_fill = map { ${$_}{fill} } @{$rescale_ref};

	# Add data to plot
	$self->svg_add_track( \@x_values, \@regular_y , \@regular_line , \@regular_fill, "$i-$start-$stop");
	$self->svg_add_track( \@x_values, \@rescale_y , \@rescale_line , \@rescale_fill ,"$i-$start-$stop");

	# scale
	$self->make_scale( $i, $start, $stop );
}

sub debug {
	my ( $self, $title, @arrs ) = @_;
	print "=" x 16 . "\n";
	foreach (@arrs) {
		my @arr = @{$_};
		printf "Array %s : %s\n", $title, scalar @arr;
		print "\t"
		  . join( ',',
			map { sprintf( '%-10d', int($_) ) } @arr[ 0 .. 10 ] )
		  . "\n";
		my $a = scalar(@arr) - 11;
		my $b = scalar(@arr) - 1;
		print "\t"
		  . join( ',',
			map { sprintf( '%-10d', int($_) ) } @arr[ $a .. $b ] )
		  . "\n";
	}
}

sub svg_add_track {
	my ( $self, $x_values_ref, $data_ref, $line_ref, $fill_ref, $base_track_id) = @_;
	my @x_values = @{$x_values_ref};
	my @data     = @{$data_ref};
	my @lines    = @{$line_ref};
	my @fills    = @{$fill_ref};

	my $it = each_array( @data, @lines, @fills );
	while ( my ( $pry, $prl, $prf ) = $it->() ) {
		my $plot_data = $self->svg()->get_path(
			x     => \@x_values,
			y     => $pry,
			-type => 'polyline',
			-closed => 'false' #specify that the polyline is closed.
		);
		$self->plot_track( $plot_data, $prl, $prf, "$base_track_id" );
	}
}

sub plot_data {
	my ( $self, %d ) = @_;
	$self->add_header($d{regular},$d{rescale});

	##loop through rows
	foreach ( my $i = 0 ; $i < $self->num_rows() ; $i++ ) {
		my ( $start, $stop ) =
		  ( $i * $self->row_size(), ( $i + 1 ) * $self->row_size() );
		$self->plot_individual_row( $start, $stop, $i,
			$d{regular},
			$d{rescale},
		);
	}
}

sub plot_data_subset {
	my ( $self, %d ) = @_;

	$self->add_header($d{regular},$d{rescale});
	##loop through rows
	$self->plot_individual_row( $d{from}, $d{to}, 0, $d{regular},
		$d{rescale}, );
}

sub push_all {
	my ( $self, $array_ref, @values ) = @_;
	foreach ( @{$array_ref} ) {
		push( @{$_}, @values );
	}
}

no Moose;
1;
