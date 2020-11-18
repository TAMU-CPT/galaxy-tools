#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use File::Temp qw/ tempdir tempfile /;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
				file_format => ['genbank', 'embl', 'fasta', 'txt'],
				multiple => 1,
			}
		],
		[
			'label',
			'Label (per input file)',
			{
				required => 1,
				validate => 'String',
				multiple => 1,
			}
		],
		[
			'zoom',
			'How zoomed in the image is. Be careful with this option. It represents the number of bases to plot in a single pixel. For large genomes, this can mean very large images, and should be lowered appropriately. For a value of 50, 50 bases would be considered a single pixel in the output image. For 1Mbp of genomes totaly (say 5 x 200 kb phages), this would result in a 20,000 pixel image.',
			{
				required => 1,
				validate => 'Int',
				min => 20,
				default => 50,
			}
		],
		[],
		[
			'java_7_bin',
			'Location of the Java binary',
			{
				required => 1,
				validate => 'File/Input',
				default =>
				  '/usr/lib/jvm/java-7-openjdk-amd64/bin/java',
				hidden           => 1,
				_galaxy_specific => 1,
				_show_in_galaxy  => 0,

			}
		],
		['local', 'Running locally under test environment. Only affects paths', { validate => 'Flag', hidden => 1, _show_in_galaxy => 0}],
	],
	'outputs' => [
		[
			'dotplot',
			'MIST Plot',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'mist',
				data_format    => 'text/html',
				default_format => 'HTML',
			}
		],
	],
	'defaults' => [
		'appname' => 'MIST',
		'appid'   => 'MIST',
		'appvers' => '1.94',
		'appdesc' =>
'Multiple Interrelated Sequence doT plotter. Uses a stripped down vesion of Gepard (Dr. Jan Krumsiek/HelmholtzZentrum/IBIS) for dot plotting.',
	],
	'tests' => [
	],
);

# Process GBKs to Fasta Files
my @files = map { { file => $_ } } @{ $options->file };
my @labels = @{ $options->label };
# Some gbks/FA have multiple entries. In that case we'll be appending multiple file/label objects.
my @new_files;

# Create a temp dir for this run, keeps our data cleaner.
my $dir = tempdir('cpt.mist.XXXXXXX',CLEANUP => 1 );

# Check that our args are OK
if( scalar(@files) != scalar(@labels) ){
	die 'please provide as many labels as files';
}


# Loop over input files
for ( my $i = 0 ; $i < scalar(@files) ; $i++ ) {
	my %file = %{ $files[$i] };
	#if ( $labels[$i] ne '__none__' && $labels[$i] ne 'None' ) {
		#$file{custom_label} = $labels[$i];
	#}

	# Parse out data, possibly will return multiple "files"
	my @parsed_new_files = extract_info_from_file($file{file}, $labels[$i]);
	# Copy to new_files
	push(@new_files, @parsed_new_files);
}

sub extract_info_from_file {
	my ($file, $default_label) = @_;
	use CPT::Bio;
	my $bio = CPT::Bio->new();
	my %args = (
		'file'      => $file,
		'header'    => 1,
		'subset'    => 'whole',
	);

	my @array_of_new_file_info;
	foreach my $arrayref($bio->parseFile(%args)){
		foreach my $subseqref(@{$arrayref}){
			my %local_file_info;

			my ($header,$seq) = @{$subseqref};
			my ( $fh, $path ) = tempfile(
				'XXXXXXXXX',
				UNLINK => 1,
				DIR    => $dir,
				SUFFIX => '.fa'
			);

			# /tmp/YYYYYYYYY/XXXXXXXXX.fasta
			my $id = substr( $path, rindex( $path, '/' ) + 1, -6 );

			# Store the info
			$local_file_info{fasta_path}   = $path;
			$local_file_info{id}           = $id;
			$local_file_info{header}       = substr($header,1), #strip >
			# We'll just go ahead and set this here. That way we
			# can remove logic to check header vs custom later
			$local_file_info{custom_label} = substr($header,1), #strip >

			# Create fasta file
			print $fh "$header\n$seq\n";

			push(@array_of_new_file_info, \%local_file_info);
		}
	}

	# If there was only one entry in this file, we'll set custom label to
	# whatever the user asked for. Otherwise we had multiple entries in the
	# file, and so we take the sequence's ID instead.
	if(scalar(@array_of_new_file_info) == 1){
		$array_of_new_file_info[0]{custom_label} = $default_label;
	}

	# Return AoH
	return @array_of_new_file_info;
}


# Prepare a reverse index because we'll eventually want to lookup by ID
my %id_to_file_ridx;
foreach (@new_files) {
	my %f = %{$_};
	$id_to_file_ridx{ $f{id} } = \%f;
}

# Process the FASTA Files
use IPC::Run3;
my ( $gin, $gout, $gerr );

# Create the directories we'll need
use File::Spec::Functions qw/catfile catdir/;
mkdir( catdir( $dir, 'png' ) );
mkdir( catdir( $dir, 'thumb' ) );
mkdir( catdir( $dir, 'prefinal' ) );
mkdir( catdir( $dir, 'final' ) );

my @img_array;

# Get gepard jar location
use File::ShareDir;
use File::Spec;
my $data_dir = File::ShareDir::dist_dir('CPT-CLI');
my $gepard_location= File::Spec->catfile(
	$data_dir,
	'Gepard_Stripped.jar'
);

# Prevent stupid errors.
if($options->{local}){
	$gepard_location = '../data/Gepard_Stripped.jar';
}
if(!-e $gepard_location){
	die 'Could not find specified Gepard JAR file';
}

sub run_gepard{
	my ($seq_from, $seq_to, $zoom, $output) = @_;

	my ($gin,$gout,$gerr);
	my @command = (
		$options->{'java_7_bin'},
		'-jar',
		$gepard_location,
		'-seq1',
		$seq_from,
		'-seq2',
		$seq_to,
		'-matrix',
		'edna.mat',
		'-zoom', $zoom,
		'-outfile',
		$output,
	);
	run3 \@command, \$gin, \$gout, \$gerr;
}

sub resize_image {
	my ($scaling, $from, $to) = @_;
	my @command = ( "convert", '-resize', $scaling, $from, $to );

	my ($gin,$gout,$gerr);
	run3 \@command, \$gin, \$gout, \$gerr;
}

# Loop through and down
for ( my $i = 0 ; $i < scalar(@new_files) ; $i++ ) {

	# Each row
	$img_array[$i] = [];
	for ( my $j = 0 ; $j < scalar(@new_files) ; $j++ ) {

		# Collect our objects
		my %file0 = %{ $new_files[$i] };
		my %file1 = %{ $new_files[$j] };

		run_gepard(
			$file0{'fasta_path'},
			$file1{'fasta_path'},
			$options->{zoom},
			catfile(
				$dir,
				'png',
				sprintf( '%s-%s.png', $file0{id}, $file1{id}, ),
			)
		);

		#print join( ' ', @command ), "\n";
		my $png_loc =
		  catfile( $dir, 'png',
			sprintf( '%s-%s.png', $file0{id}, $file1{id} ) );
		my $thumb_loc =
		  catfile( $dir, 'thumb',
			sprintf( '%s-%s.png', $file0{id}, $file1{id} ) );

		# Resize our image
		resize_image('30%',$png_loc,$thumb_loc);

		$img_array[$i][$j] = {
			from  => $file0{id},
			to    => $file1{id},
			orig  => $png_loc,
			thumb => $thumb_loc,
		};
	}
}

my $CARDINALITY         = scalar(@new_files);
my $INTER_IMAGE_BORDERS = 2;
my $IMAGE_BORDER        = 50;

# Create montage
my @transposed;
for my $row (@img_array) {
	for my $column ( 0 .. $#{$row} ) {
		push( @{ $transposed[$column] }, $row->[$column] );
	}
}

my @complete_image_list;
for my $row (@transposed) {
	my @r = @{$row};
	push( @complete_image_list, map { ${$_}{thumb} } @r );
}
my @command = (
	"montage",
	@complete_image_list,
	"-tile",
	"${CARDINALITY}x${CARDINALITY}",
	'-geometry',
	'+0+0',
	'-border',
	$INTER_IMAGE_BORDERS,
	'-bordercolor',
	'purple',
	catfile( $dir, 'prefinal', 'montage_0.png' )
);

run3 \@command, \$gin, \$gout, \$gerr;

# Add border
@command = (
	'convert',
	catfile( $dir, 'prefinal', 'montage_0.png' ),
	'-bordercolor',
	'purple',
	'-border',
	sprintf( "%dx%d", 1, 1 ),
	catfile( $dir, 'prefinal', 'montage_1.png' ),
);
run3 \@command, \$gin, \$gout, \$gerr;
@command = (
	'convert',
	catfile( $dir, 'prefinal', 'montage_1.png' ),
	'-bordercolor',
	'gray',
	'-border',
	sprintf( "%dx%d", $IMAGE_BORDER, $IMAGE_BORDER ),
	catfile( $dir, 'prefinal', 'montage_2.png' ),
);
run3 \@command, \$gin, \$gout, \$gerr;

for ( my $i = 0 ; $i < scalar(@new_files) ; $i++ ) {
	for ( my $j = 0 ; $j < scalar(@new_files) ; $j++ ) {
		my %current_image = %{ $img_array[$i][$j] };

		#from/to/orig/thumb
		my ( $in, $out, $err );
		@command = ( 'identify', $current_image{thumb} );
		run3 \@command, \$in, \$out, \$err;

#/tmp/Q0Gk8V53PV/thumb/Db5GU0Osc-Db5GU0Osc.png PNG 291x291 291x291+0+0 8-bit PseudoClass 256c 21KB 0.000u 0:00.010
#/tmp/Q0Gk8V53PV/thumb/9LuA65WNY-9LuA65WNY.png PNG 1342x1342 1342x1342+0+0 8-bit PseudoClass 256c 841KB 0.000u 0:00.000
		if ( $out =~ /(\d+)x(\d+)\+0\+0/ ) {
			$current_image{height} = $1;
			$current_image{width}  = $2;
			$img_array[$i][$j]     = \%current_image;
		}
		else {
			printf STDERR (
				"In: %s\nOut:%s\nErr: %s\nCmd: %s\n",
				$in, $out, $err, join( " ", @command )
			);
			die 'Could not get size for ' . $current_image{thumb};
		}
	}
}

#print Dumper @img_array;

my @tmp = map { ${$_}{width} } @{ $img_array[0] };
use List::Util qw(sum);
my $img_size_without_borders = sum(@tmp);
my $border_size              = 2 * $INTER_IMAGE_BORDERS * $CARDINALITY;

# The +1 and +2 are as a result of adding a 1 width purple border, so the border is consistent everywhere.
my $total_size_two_border =
  $img_size_without_borders + $border_size + $IMAGE_BORDER + $IMAGE_BORDER + 2;
my $total_size_one_border =
  $img_size_without_borders + $border_size + $IMAGE_BORDER + 1;

# Creating labels

my $current_sum = $IMAGE_BORDER + $INTER_IMAGE_BORDERS;
my @convert_arguments_top;
my @convert_arguments_left;
my $counter     = 0;
my $left_offset = $total_size_one_border + 20;

for ( my $i = 0 ; $i < scalar(@new_files) ; $i++ ) {
	my %current_image = %{ $img_array[$i][0] };
	if ( ${ $id_to_file_ridx{ $current_image{from} } }{custom_label} ) {
		$a =
		  ${ $id_to_file_ridx{ $current_image{from} } }{custom_label};
	}
	else {
		$a = ${ $id_to_file_ridx{ $current_image{from} } }{header};
	}

	push( @convert_arguments_top,
		'-fill', 'black', '-annotate',
		sprintf( '+%s+30', $current_sum ), $a );
	push( @convert_arguments_left,
		'-fill', 'black', '-annotate',
		sprintf( '+%s+%s', $current_sum, $left_offset ), $a );
	$current_sum +=
	  $current_image{width} + $INTER_IMAGE_BORDERS + $INTER_IMAGE_BORDERS;
}

# Applying labels
@command = (
	'convert',
	catfile( $dir, 'prefinal', 'montage_2.png' ),
	'-rotate',
	'-90',
	'-pointsize',
	'24',
	'-font',
	'Ubuntu-Mono-Regular',
	@convert_arguments_left,
	'-rotate',
	'90',
	@convert_arguments_top,
	'-pointsize',
	'14',
	'-annotate',
	sprintf( '+%s+%s', $IMAGE_BORDER, $left_offset ),
"Produced by the CPT's MIST (Multiple Interrelated Sequence doT plotter). Written by Eric Rasche <rasche.eric\@yandex.ru>.\nDot plots produced by the Gepard Dot Plotter by Dr. Jan Krumsiek",
	catfile( $dir, 'final', 'large.png' ),
);
run3 \@command, \$gin, \$gout, \$gerr;


resize_image('50%',catfile( $dir, 'final', 'large.png' ), catfile( $dir, 'final', 'small.png' ));

use File::Copy;
my @files_to_move;
push(@files_to_move,
	[catfile( $dir, 'final', 'large.png' ),'large'],
	[catfile( $dir, 'final', 'small.png' ),'small'],
);

# Now we'll need to move all the files to a hierarchy specified by galaxy. Yuck!
my $mistmap = "";

my $scaling = 0.5;
# Vertical Loop
my $cur_y   = $scaling * $IMAGE_BORDER + 1;
for ( my $i = 0 ; $i < scalar(@new_files) ; $i++ ) {
	# Horizontal loop
	my $cur_x = $scaling * $IMAGE_BORDER + 1;
	for ( my $j = 0 ; $j < scalar(@new_files) ; $j++ ) {
		my %current_image = %{ $img_array[$i][$j] };
		my $new_filename  = sprintf( "%s-%s", $current_image{from}, $current_image{to} );
		push(@files_to_move, [$current_image{orig}, $new_filename ]);


		my $current_image_width = $current_image{width};
		my $current_image_height = $current_image{height};

		$mistmap .= sprintf('<area shape="rect" coords="%s,%s,%s,%s" alt="%s" href="%s" />',
			$cur_x,
			$cur_y,
			int(
				
				  ( $cur_x + $scaling *$current_image{width} + 1 )
			),
			int(
				  ( $cur_y + $scaling *$current_image{height} + 1 )
			),
			$current_image{from} . $current_image{to},
			sprintf( "%s-%s.png",
				$current_image{to}, $current_image{from} ),
		) . "\n";
		$cur_x += $scaling * ( $current_image{width} + 2 );
	}
	my %tmp = %{ $img_array[$i][0] };
	$cur_y += $scaling * ( $tmp{height} + 2 );
}

my $html = sprintf( '
<!DOCTYPE html>
<html>
<body>
<h1>Mist Results</h1>
<p>Each section of mist output is now clickable to view a higher resolution version of that subsection</p>
<p>Additionally, a high resolution version is available <a href="large.png">here</a></p>
<img src="small.png" usemap="#mistmap">
<map name="mistmap">
%s
</map>
</body>
</html>
',
	$mistmap,
);

use CPT::OutputFiles;
my $mist_output = CPT::OutputFiles->new(
	name => 'dotplot',
	GGO => $ggo,
);
$mist_output->CRR(data => $html);

#We prepared an array of files to move with their original location and a new
#filename. Extensions are assumed to be PNG in this case.
foreach my $move(@files_to_move){
	my ($original, $new) = @{$move};
	my @produced_files = $mist_output->subCRR(data_format => "Dummy", format_as => "Dummy", filename=>$new, extension => 'png');
	move($original,$produced_files[0]);
}

=head1 DESCRIPTION

MIST produces an N-by-N dot matrix of a given set of input genomes (fasta or genbank).

=cut
