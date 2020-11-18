package CPT::Plot::Base;
use Data::Dumper;
use CPT::Plot::Label;
use CPT::Plot::Class;
use CPT::Plot::Gene;
use CPT::Plot::Colours;
use Bio::SeqIO;
use SVG;
use Moose;

# ABSTRACT: Main plotting class for genome mapper

has 'svg' => ( is => 'rw', isa =>'Any');
has 'line_count' => ( is => 'rw', isa =>'Num', default => 1);
has '_ft_count' => ( is => 'rw', isa =>'Num', default => 0);
has 'classes' => ( is => 'rw', isa =>'HashRef');

# Labels 
has 'label' =>( is=>'rw',isa=>'Bool');

has 'label_pos' => (is => 'rw', isa => 'Any');
has 'label_shrink_mode' => (is => 'rw', isa => 'Any');
has 'label_callouts' => (is => 'rw', isa => 'Any');
has 'label_from' => (is => 'rw', isa => 'Any');
has 'label_text_source' => (is => 'rw', isa => 'Any');
has 'label_numeric_features' => (is => 'rw', isa => 'Any');
has 'label_query' => (is => 'rw', isa => 'Any');
has 'label_numbering_count' => ( is => 'rw', isa => 'Any', default => 1);

has 'justified' =>( is=>'rw',isa=>'Str');
# CHanged to any b/c unpassed = undef
has 'separate_strands' =>( is=>'rw',isa=>'Any');
has 'double_line_for_overlap' =>( is=>'rw',isa=>'Any');
has 'opacity' =>( is=>'rw',isa=>'Str');
has 'view' =>( is=>'rw',isa=>'Str');

has 'color_scheme' => ( is => 'rw', isa =>'HashRef');
has 'wanted_tags' => (is => 'rw', isa => 'HashRef');
has 'genome_length' => (is => 'rw', isa => 'Int');
has 'features' => (is => 'rw',isa => 'ArrayRef');
has 'start' => (is => 'rw', isa => 'Int');
has 'end' => (is => 'rw', isa => 'Int');

has 'avgRowLength' => (is => 'rw', isa => 'Int');
has 'calc_height' => (is => 'rw', isa => 'Int');
has 'calc_width' => (is => 'rw', isa => 'Int');
has 'x_offset' => (is => 'rw', isa => 'Num');
has 'y_offset' => (is => 'rw', isa => 'Num');
has 'ils' => (is => 'rw', isa => 'Num');
has 'zoom' => (is => 'rw', isa => 'Num');
has 'rows' => (is => 'rw', isa => 'Num');
has 'split_factor' => (is => 'rw', isa => 'Num');

has 'rowdata' => (is => 'rw', isa => 'HashRef');
has '_internal_maxrowlength' => (is => 'rw', isa => 'Num');

my $color_spec = CPT::Plot::Colours->new( 'default' => '#000000' );
our ( $parser, $tree, $cb );

sub init {
	my ($self) = @_;
	my %classes;
	my %cs = %{$self->color_scheme()};
	foreach my $key(keys %cs){
		$classes{ $key } = CPT::Plot::Class->new(
			'key'      => $key,
			'color'    => $cs{$key}{color},
			'border'   => $cs{$key}{border},
			'plot'     => $cs{$key}{plot},
			'included' => 1,
		);
	}
	$self->classes(\%classes);
	$self->init_label_stuff();
	$self->filterFeatues();
}
sub init_label_stuff {
	my ($self) = @_;

	if ( $self->{'label_from'} eq 'custom' ) {
		use Parse::BooleanLogic;
		$parser = new Parse::BooleanLogic( operators => [ '', 'OR' ] );
		$tree = $parser->as_array( $self->label_query );
		print $parser;

		#foreach bio feature,
		#if solve == 1, then add to our return,
		#else doesn't match
		#endforeach
		#my $new_tree = $parser->solve($tree,$filter);
		$cb = sub {
			my $query   = $_[0]->{'operand'};
			my $feature = $_[1];

			my $negate = 0;
			if ( substr( $query, 0, 1 ) eq '!' ) {    #negate
				$negate = 1;
				$query = substr( $query, 1 );
			}
			if ( $query =~ m/([^:]+):["']{0,1}([^'"]*)["']{0,1}/ ) {
				my ( $k, $v ) = ( $1, $2 );
				my $result;
				if ( $k eq 'contains' ) {
					my $values = join(
						"\t",
						map {
							if ( $_ ne
								"translation" )
							{
								join(
									'',
									$feature
									  ->get_tag_values
									  (
										$_
									  )
								);
							}
						  } $feature->get_all_tags()
					);
					if ( $values =~ m/$v/i ) {
						$result = 1;
					}
					else {
						$result = 0;
					}
				}
				elsif ( $k eq 'key' ) {
					if ( $v =~ m/,/ ) {
						$result = 0;
						foreach ( split( /,/, $v ) ) {
							if ( $feature
								->primary_tag eq
								$_ )
							{
								$result = 1;
							}
						}
					}
					else {
						$result =
						  $feature->primary_tag eq $v;
					}
				}
				elsif ( $k eq 'tag' ) {
					if ( $v =~ m/([^=]+)=(.*)/ ) {
						my ( $tag_name, $tag_match ) =
						  ( $1, $2 );
						if ( $feature->has_tag($1) ) {
							if (
								join(
									'',
									$feature
									  ->get_tag_values
									  (
										$1
									  )
								) =~ /$2/i
							  )
							{
								$result = 1;
							}
							else {
								$result = 0;
							}
						}
						else {
							$result = 0;
						}
					}
					else {
						$result = $feature->has_tag($v);
					}
				}
				else {

					#error
					$result = 0;
				}
				return ( $negate xor $result );
			}
			else {

				#error
				return 0;
			}

			#error
			return 0;
		};
	}
}
sub filterFeatues {
	my ( $self ) = @_;

	#$self->{'wanted_tags'} = map { $_ => 1 } split(/,/,$self->{'q'});
	my %tags = map { $_ => 1 } split( /,/, "tRNA,CDS" );
	$self->wanted_tags(\%tags);
	my @feats = @{$self->features()};
	for my $feat_object ( @feats ) {
		my $should_add = 1;
		if ( $feat_object->primary_tag eq 'source' ) {
			$should_add = 0;
		}
		if(defined $self->start() &&
			$feat_object->start < $self->start()){
			$should_add = 0;
		}
		if(defined $self->end() &&
			$feat_object->end > $self->end()){
			$should_add = 0;
		}
		if($should_add){
			$self->addGene($feat_object);
		}
	}
}
sub addGene {
	my ( $self, $feat_object ) = @_;
	my $tag   = $feat_object->primary_tag;
	my $label = "";
	if ( $self->label() ) {

		#If it meets the criteria specified for labelling an object, set the label, else don't set a label
		if ( $self->label_from() eq 'custom') {
			if($parser->solve($tree,$cb,$feat_object)){
				if($feat_object->has_tag($self->label_text_source())){
					$label = join(' ', $feat_object->get_tag_values($self->label_text_source()));
				}else{
					$label = '[]';
				}
			}
			#if($feat_object->has_tag($self->label_text_source())){
				#$label = ' '.join(' ', $feat_object->get_tag_values($self->label_text_source()));
			#}
		}
		elsif ( $self->label_from() eq 'numeric') {
			if ( ${$self->wanted_tags()}{$tag} ) {
				$label = $self->label_numbering_count();
				$self->label_numbering_count($self->label_numbering_count()+1);
			}
		}else{
			die $self->label_from();
		}
	}
	my @color_arr;
	my $color;
	if ( $feat_object->has_tag('color') ) {
		push( @color_arr, $feat_object->get_tag_values('color') );
	}
	if ( $feat_object->has_tag('color') ) {
		push( @color_arr, $feat_object->get_tag_values('color') );
	}
	if ( scalar @color_arr ) {
		$color = $color_arr[0];
	}

	my $gene = CPT::Plot::Gene->new(
		'tag'    => $tag,
		'label'  => $label,
		'start'  => $feat_object->start,
		'end'    => $feat_object->end,
		'strand' => $feat_object->strand,
		'color' => $color,
	);

#This is a "failsafe" addition of classes, in case the user didn't specify a color
	if( ! defined ${$self->classes()}{$tag} ) {
		${$self->classes()}{$tag} = CPT::Plot::Class->new(
			'key'         => $tag,
			'color'      => '#000000',
			'border'      => 1,
			'plot'        => 1,
			'included' => 1,
		);
	}else{
		${$self->classes()}{$tag}->addObject($gene);
	}
}

sub partitionLines {
	my ($self) = @_;

# To use when I finally get partitioning.pm working
#sub partitionLines{
#	my ($self) = @_;
#
#	my $partioner = Partitioning->new(
#		genome_length => $self->genome_length(),
#		rows          => $self->rows(),
#		justified     => $self->justified(),
#	);
#
#	# Add data to it
#	foreach(keys %classes){
#		if($classes{$_}->isIncludedInPartioning()){
#			$partioner->add($classes{$_}->getItemList());
#		}
#	}
#	# Run &&  get Results
#	my %result = %{$partioner->run()};
#	# . . .
#	print Dupmer %results;
#	# Profit
#	exit 1;
#	# This is supposed to merge two hashes. [http://perldoc.perl.org/perlfaq4.html#How-do-I-merge-two-hashes%3f]
#	@self{keys %result} = values %result;

	my @items;

	$self->avgRowLength(
	  int( $self->genome_length() / $self->rows() * $self->split_factor() )
	); #TODO, allow adjusting or just re-calc? need to benchmark first I guess.
	$self->calc_height(
	  int( ( 1 + $self->rows() )
		  * $self->ils() ));
	$self->calc_width(
	  int( $self->avgRowLength()
		  / $self->zoom() ));

	my $fake_count = 100;
	if ($fake_count) {
		for ( my $i = 0 ; $i <= $fake_count ; $i++ ) {
			my $key =
			  int( $self->genome_length() * $i / $fake_count );
			push( @items, [ $key, $key, 1 ] );
		}
	}

	my %classes = %{$self->classes()};
	foreach ( keys %classes ) {
		if ( $classes{$_}->included() ) {
			push( @items, @{ $classes{$_}->getItemList() } );
		}
	}

	#Sort based on where each item starts
	@items = sort { ${$a}[0] <=> ${$b}[0] } @items;

	#my $z = '(' . join('),(',map { "${$_}[0],${$_}[1]" } @items ) . ')';
	#print join("\n",split(/(.{1,120})/msxog, $z)) . "\n";
	my %rowdata;

	my ( $longest_last_object, $thisRowEnd, $currentRow ) = ( 1, 1 + $self->avgRowLength(), 1 );
	$rowdata{1}{start} = 1;
	foreach my $item_ref (@items) {
		my ( $item_start, $item_end ) = @{$item_ref};

		#print "\t$item_start\t$item_end\t$thisRowEnd\n";
		if ( $item_start >= $thisRowEnd || $item_end > $thisRowEnd ) {

       # This was just cleaned up from the following commented out piece of code
			if (       $self->justified() eq 'justify'
				|| $item_start >=
				$rowdata{$currentRow}{end} )
			{
				$rowdata{$currentRow}{end} =
				  $thisRowEnd;
			}
			else {
				$rowdata{$currentRow}{end} =
				  max( $longest_last_object, $item_start );
			}

# There was a corner case here:
# O    represents the end of a gene,
# ---  represents a gene
# |    represents $thisRowEnd
#
#
# ------O     |  O---------
# In this case, the second end would be chosen as
# max($longest_last_object,$item_start), which is NOT what we
# want.  You want  | to be chosen, not O, so in the case that
# item_start is >= current row end (or should that be >?), we
# use this.
#
# ------O     |
#          O--+--------
# This case works fine
#
#
# ------O     |
#    O--------+--------
# This case also works fine
#
#
#				if($self->justified()){
#					$rowdata{$currentRow}end() = $thisRowEnd;
#				}else{
#					if($item_start <= $rowdata{$currentRow}end()){
#						$rowdata{$currentRow}end() = max($longest_last_object,$item_start);
#					}else{
#						$rowdata{$currentRow}end() = $thisRowEnd;
#					}
#				}
			$self->_internal_maxrowlength(max(
				$self->_internal_maxrowlength(),
				$rowdata{$currentRow}{end} -
				  $rowdata{$currentRow}{start}
			));
			$currentRow++;

		#print "$item_start $rowdata{$currentRow-1}{end}\n";
			if ( $item_start <=
				$rowdata{ $currentRow - 1 }{end} )
			{
				$rowdata{$currentRow}{start} =
				  $item_start;
			}
			else { #nonjustified never encounters the following line
				$rowdata{$currentRow}{start} =
				  $rowdata{ $currentRow - 1 }{end} +
				  1;
			}
			$thisRowEnd =
			  $self->avgRowLength() +
			  $rowdata{$currentRow}{start};
		}
	}

#	if($self->justified()){
#		foreach my $item_ref(@items){
#				my ($item_start, $item_end) = @{$item_ref};
#				# If the item starts OR ends after this row is supposed to end
#				# print "\t$item_start\t$item_end\t$thisrowend\n";
#				if($item_start >= $thisRowEnd || $item_end >  $thisRowEnd){
#					$rowdata{$currentRow}end() = $thisRowEnd;
#					#Internal max row length is the length of the longest row
#					$self->_internal_maxrowlength'} = max($self->{'_internal_maxrowlength'},$rowdata{$currentRow}{'end'}-$rowdata{$currentRow}{'start());
#					#Update which row we're on (so we aren't using +1s everywhere)
#					$currentRow++;
#					if($item_start <= $rowdata{$currentRow-1}end()){
#						$rowdata{$currentRow}start() = $item_start;
#					}else{
#						$rowdata{$currentRow}start'} = $rowdata{$currentRow-1}{'end() + 1;
#					}
#					#tracks where the current row ends
#					#print Dumper $rowdata;
#					#print ">>$thisRowEnd\t".$self->avgRowLength'}." + ".$rowdata{$currentRow}{'start()."\n";
#					$thisRowEnd = $self->avgRowLength'} + $rowdata{$currentRow}{'start();
#					#print ">>$thisRowEnd\t".$self->avgRowLength'}." + ".$rowdata{$currentRow}{'start()."\n";
#				}
#		}
#	}else{#Non justified, raggedright
#		foreach my $item_ref(@items){
#				my ($item_start, $item_end) = @{$item_ref};
#				#print "\t$item_start\t$item_end\t$thisrowend\n";
#				if($item_start >= $thisRowEnd || $item_end >  $thisRowEnd){
##					print "\t> $item_start\t$item_end\t$thisRowEnd\n";
##					print "Candidate for ending [" . ($item_start >= $thisRowEnd) ."]\t[" .($item_end >= $thisRowEnd) . "]\n";
##					# If we have ``justified'' rulers, they all need to the be the SAME length (until the last)
##					print "              -- $rowdata{$currentRow}end()$thisRowEnd\n";
#					$rowdata{$currentRow}end() = max($longest_last_object,$item_start);
#					#Internal max row length is the length of the longest row
#					$self->_internal_maxrowlength'} = max($self->{'_internal_maxrowlength'},$rowdata{$currentRow}{'end'}-$rowdata{$currentRow}{'start());
#					#Update which row we're on (so we aren't using +1s everywhere)
#					$currentRow++;
#					#if($item_start <= $rowdata{$currentRow-1}end()){
#						$rowdata{$currentRow}start() = $item_start;
#					#}
#					#tracks where the current row ends
#					$thisRowEnd = $self->avgRowLength'} + $rowdata{$currentRow}{'start();
#				}
#				$longest_last_object = max($longest_last_object,$item_end);
#		}
#	}
#make sure the final row length is set, in addition to the _int_max_rowlength
	$thisRowEnd = $rowdata{$currentRow}{end} =
	  $self->genome_length() + 1;    #Putative
	$self->_internal_maxrowlength(max(
		$self->_internal_maxrowlength(),
		$rowdata{$currentRow}{end} -
		  $rowdata{$currentRow}{start}
	));
	$rowdata{max} = $currentRow;

	if(defined $self->{start} && defined $self->{end}){
		%rowdata = (
			'1' => { 'end' => $self->{end}, 'start' => $self->{start} },
			'max' => 1,
		);
	}

	$self->rowdata(\%rowdata);


}
sub getSVG {
	my($self) = @_;
	return $self->svg();
}
# SVG
sub createSVG {
	my ($self) = @_;
	my %rowdata = %{$self->rowdata()};
	$self->calc_height( int( ( 1 + $rowdata{max} ) * $self->ils() ));
	$self->calc_width(int( $self->avgRowLength() /  $self->zoom()));

	$self->svg(SVG->new(
		width  => $self->calc_width() + 2 * $self->x_offset(),
		height => $self->calc_height() + 2 * $self->y_offset(),
	));
	#$self->svg()->title( id => 'documenfeatures from t-title' )->cdata("Genome Map of [$file_name]");

	my $ui_group = $self->svg()->tag(
		'g',
		id    => 'group_ui',
		style => {
			stroke         => '#000000',
			fill           => '#000000',
			'fill-opacity' => 1,
		}
	);

	foreach ( my $i = 1 ; $i <= $rowdata{max} ; $i++ ) {
		$self->_addRuler( $i, $ui_group );
	}

	my %classes = %{$self->classes()};
	foreach my $class_key ( keys %classes ) {
		#print "Adding features from $class_key\n";
		my $class = $classes{$class_key};
		if ( !$class->plot() ) {
			next;
		}
		my $group = $self->svg()->tag(
			'g',
			id    => 'group_' . $class->key(),
			style => {
				stroke => (
					$class->plot()
					? (
						$class->border()
						? "black"
						: "none"
					  )
					: 'none'
				),
				fill           => $class->color(),
				'fill-opacity' => $self->opacity(),
			}
		);
		my @data = @{ $class->getObjects() };
		foreach my $gene (@data) {
			my ( $start, $end ) =
			  ( $gene->start(), $gene->end() );
			my $row = calculateRow( $self, $start, $end );
			addFeature(
				$self,
				group    => $group,
				row      => $row,
				start    => $start,
				end      => $end,
				key      => $gene->tag(),
				strand   => $gene->strand(),
				label    => $gene->label(),
				ui_group => $ui_group,
				color    => $gene->color(),
			);

		}
	}

}
sub calculateRow {
	my ( $self, $start, $end ) = @_;
	my %rowdata = %{$self->rowdata()};
	for ( my $i = 1 ; $i <= $rowdata{max} ; $i++ ) {
		if (
			   $start > $rowdata{$i}{start} - 1
			&& $start < $rowdata{$i}{end} + 1
			&& $end > $rowdata{$i}{start} - 1
			&& $end < $rowdata{$i}{end} + 1

		  )
		{
			return $i;
		}
	}

#print "<b>$start,$end,".$self->rowdata'}{$i}{'start'}.",".$self->{'rowdata'}{$i}{'end()."<\/b>\n";
	return 1.5;
}
sub _addRuler {
	my ( $self, $row, $ui_group ) = @_;
	my $y_fix = $self->ils() * ( $row - 1 );

	#	my @d = (
	#		$self->calc_width(),
	#		$self->rowdata'}{$row}{'end(),
	#		$self->rowdata'}{$row}{'start(),
	#		($self->rowdata'}{$row}{'end'}-$self->{'rowdata'}{$row}{'start()),
	#		$self->_internal_maxrowlength(),
	#	);
	#	print join("\t",@d),"\n";
	my %rowdata = %{$self->rowdata()};
	my $line_width =
	  $self->calc_width() *
	  ( $rowdata{$row}{end} -
		  $rowdata{$row}{start} ) /
	  $self->_internal_maxrowlength();

#print "Adding ruler\t".$self->rowdata'}{$row}{'start'}."\t".$self->{'rowdata'}{$row}{'end'}."\t" . ($self->{'rowdata'}{$row}{'end'} - $self->{'rowdata'}{$row}{'start()) . "\n";

	$ui_group->line(
		id => 'ui_element_' . ( $self->line_count() + rand() ),
		x1 => 0 + $self->x_offset(),
		x2 => $line_width + $self->x_offset(),
		y1 => $y_fix + $self->y_offset(),
		y2 => $y_fix + $self->y_offset()
	);

	#	print "Ruler is being plotted from $y_fix to $line_width\n";
	if ( $self->separate_strands() ) {
		#$ui_group->rectangle(
			#id => 'ui_element_' . ( $self->line_count() + rand() ) . "_" . rand(1),
			#x      => 0 + $self->x_offset(),
			#y      => $y_fix - 2.5 + $self->y_offset(),
			#width  => $line_width,
			#height => 5
		#);

		#$y_fix += 100;
	}

	if ( $self->double_line_for_overlap() && $row > 1 )
	{    #This shows any duplicated part of the scale
		if ( $rowdata{ $row - 1 }{end} -
			$rowdata{$row}{start} >= 0 )
		{    #Equal to zero indicates ONE base of overlap
			$ui_group->line(
				id => 'ui_element_' . ( $self->line_count() + rand() ),
				y1 => $y_fix - 5 + $self->y_offset(),
				y2 => $y_fix - 5 + $self->y_offset(),
				x1 => 0 + $self->x_offset(),
				x2 => $self->calc_width() * (
					$rowdata{ $row - 1 }{end} -
					  $rowdata{$row}{start}
				  ) / $self->_internal_maxrowlength() +
				  $self->x_offset(),

#$calc_width*($rowdata{$row-1}end'}-$rowdata{$row}{'start'})/$self->{'_internal_maxrowlength'} + $self->{'x_offset(),
			);
		}
	}
	$ui_group->line(
		id => 'ui_element_' . ( $self->line_count() + rand() ),
		x1 => 0 + $self->x_offset(),
		x2 => $line_width + $self->x_offset(),
		y1 => $y_fix + $self->y_offset(),
		y2 => $y_fix + $self->y_offset()
	);
	foreach ( $rowdata{$row}{start}-1 .. $rowdata{$row}{end} )
	{
		if ( $_ % 1000 == 0 && $_ % 10000 != 0 ) {
			my $current_location =
			  $self->calc_width() *
			  ( $_ - $rowdata{$row}{start} ) /
			  $self->_internal_maxrowlength();
			$ui_group->line(
				id => 'ui_element_' . ( $self->line_count() + rand() ),
				x1 => $current_location + $self->x_offset(),
				x2 => $current_location + $self->x_offset(),
				y1 => $y_fix + $self->y_offset(),
				y2 => $y_fix + 5 + $self->y_offset(),
			);
		}
		if ( $_ % 10000 == 0 ) {
			my $current_location =
			  $self->calc_width() *
			  ( $_ - $rowdata{$row}{start} ) /
			  $self->_internal_maxrowlength();
			$ui_group->line(
				id => 'ui_element_' . ( $self->line_count() + rand() ),
				x1 => $current_location + $self->x_offset(),
				x2 => $current_location + $self->x_offset(),
				y1 => $y_fix + $self->y_offset(),
				y2 => $y_fix + 10 + $self->y_offset(),
			);
			$ui_group->text(
				id => 'ui_text' . ( $self->line_count() + rand() ),
				x => $current_location + 10 +
				  $self->x_offset(),
				y      => $y_fix + 20 + $self->y_offset(),
				-cdata => ( $_ / 1000 ) . " kb",
				'fill' => '#000000',
				'fill-opacity' => 1,
				'font-family'  => 'mono',
				'stroke'       => 'none'
			);
		}

		if ( ($_ == $rowdata{$row}{start}-1 || $_ == $rowdata{$row}{end}) && ( $_ % 10000 != 0) ) {
			my $current_location =
			  $self->calc_width() *
			  ( $_ - $rowdata{$row}{start} ) /
			  $self->_internal_maxrowlength();
			$ui_group->line(
				id => 'ui_element_' . ( $self->line_count() + rand() ),
				x1 => $current_location + $self->x_offset(),
				x2 => $current_location + $self->x_offset(),
				y1 => $y_fix + $self->y_offset(),
				y2 => $y_fix + 10 + $self->y_offset(),
			);
			$ui_group->text(
				id             => 'ui_text' . ( $self->line_count() + rand() ),
				x              => $current_location + $self->x_offset(),
				y              => $y_fix + 20 + $self->y_offset(),
				-cdata         => sprintf('%d kb', ( $_ / 1000 )),
				'fill'         => '#000000',
				'fill-opacity' => 1,
				'font-family'  => 'mono',
				'stroke'       => 'none'
			);
		}
	}
}
sub addFeature {
	my ( $self, %data ) = @_;
	my %rowdata = %{$self->rowdata()};
	my $x =
	  $self->calc_width() *
	  ( $data{'start'} - $rowdata{ $data{'row'} }{'start'} ) /
	  $self->_internal_maxrowlength() + $self->x_offset();
	my $w =
	  $self->calc_width() *
	  ( $data{'end'} - $data{'start'} ) /
	  $self->_internal_maxrowlength();
	my $h = 15;
	my $y =
	  ( $data{'row'} - 1 ) *
	  $self->ils() + $self->y_offset() - $h / 2;

	my $id = "$x$y$w$h" . rand();

#print "Item(".$data{'start'}.",".$data{'end'}.",".$data{'row'}.") =\t($x,$y,$w,$h)\n";

	if ( $self->separate_strands() ) {
		$y += -$data{'strand'} * 30;
	}

	if ( $self->view() eq 'alt_random' ) {    # Max add = 20
		$y += 4 * ( $x % 5 );
	}
	elsif ( $self->view() eq 'alt_every' ) {    # Max add = 10
		  # We (Sort of like a convolution?) multiply by strand This has
		  # the following effect; when on the top strand, we will only
		  # ever add a positive to the height of the item (moving it
		  # downward and closer to the ruler). On the bottom strand
		  # however, we only ever add a negative to the height of the
		  # item (moving it upwards towards the ruler). This allows the
		  # items on the top and bottom to stay balanced.
		$y +=
		  $data{'strand'} * 10 * ( ( $self->_ft_count($self->_ft_count() +1 ) ) % 2 ) -
		  10 * $data{'strand'};

		# However, This is imperfect, since we add items based on class,
		# not from left to right
	}
	elsif ( $self->view() eq 'alt_artemis' ) {    # Max add = 20?
		   # Muwahahahaha. Sorry. Determined coefficient and constant by
		   # trial and error, but this matches up with the artemis view
		   # without an if/else based on which strand. :D
		$y +=
		  10 * ( ( $data{'start'} - 2 * $data{'strand'} + 1 ) % 3 ) -
		  10 * $data{'strand'};
	}

	my $item_color = $color_spec->getColour( $data{'color'} );
	if ($item_color) {
		$data{'group'}->rectangle(
			x      => ($x),
			y      => $y,
			width  => $w,
			height => $h,
			id     => $id,
			fill   => $color_spec->getColour( $data{'color'} )
		);
	}
	else {
		$data{'group'}->rectangle(
			x      => ($x),
			y      => $y,
			width  => $w,
			height => $h,
			id     => $id,
		);
	}
	if ( $self->label() && $data{'label'} ) {

		my ( $lx, $ly );
		my @char_data = split( //, $data{label} );

		#Exit early if we don't even want to plot.
		my $is_too_small = ( scalar(@char_data) * 2 > $w );
		if ( $self->label_shrink_mode() eq 'cutoff' && $is_too_small ) {
			return;
		}

		#Font Scaling
		my $font_scaling = 100;
		if ( $self->label_shrink_mode() eq 'shrink') {
			$font_scaling *= $w / ( 8 * scalar(@char_data) );
		}

		# Horizontal positioning
		$lx =
		  $x +
		  $w / 2
		  ; #Horizontally center it, but this is by the leading edge of the text
		if ( scalar(@char_data) * 8 > $w && $self->label_shrink_mode() eq 'shrink') {
			$lx -=
			  scalar(@char_data) * 4 *
			  $font_scaling / 100
			  ; #Adjustment for scaled text. Determined by experiment
		}
		else {
			$lx -=
			  scalar(@char_data) * 4
			  ; #Move four pixels left for every character in the label
		}

		# Vertical positioning
		if ( $self->label_pos() eq "above" ) {    #Label is ABOVE
			if ($self->separate_strands() && $data{'strand'} == -1 )
			{
				$ly =
				  $y +
				  $h / 2 + 10 + 30
				  ; #Need to consider below strand, only one strand.
			}
			else {
				$ly =
				  $y +
				  $h / 2 - 30
				  ; #Need to consider below strand, only one strand.
			}
		}
		else {              #Label is ON
			$ly = $y + $h / 2 + 5;
		}

		$self->plot_label( $lx, $ly, $font_scaling, $data{'label'},
			$data{'ui_group'} );



		if (       $self->label_callouts()
			&& $self->label_pos() eq "above" )
		{
			$data{'ui_group'}->line(
				id => 'l' . "_" . rand(1),
				x1 => $x + ( $w / 2 ),
				x2 => $x + ( $w / 2 ),
				y1 => (
					$self->separate_strands()
					  && $data{'strand'} eq '-1' ? $y + $h
					: $y
				),
				y2 => (
					$self->separate_strands()
					  && $data{'strand'} eq '-1' ? $ly - 12
					: $ly
				)
			);
		}

	}
}
sub plot_label {
	my ( $self, $x, $y, $font_size, $label, $ui_group ) = @_;
	$ui_group->text(
		id             => 'text' . rand(1),
		x              => $x,
		y              => $y,
		-cdata         => $label,
		'fill'         => '#000000',
		'fill-opacity' => 1,
		'font-family'  => 'mono',
		'font-size'    => $font_size . '%',
		'stroke'       => 'none'
	);
}

sub max ($$) { $_[ $_[0] < $_[1] ] }
sub min ($$) { $_[ $_[0] > $_[1] ] }

no Moose;
1;
