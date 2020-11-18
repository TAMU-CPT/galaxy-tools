#!/usr/bin/perl
#
#	   Code written by Eric Rasche
#			   mailto:rasche.eric@yandex.ru
#			   tel:404.692.2048
#			   http://eric.rasche.co.uk
#	   for
#			   Center for Phage Technology
#

use strict;
use warnings;

use CPT;
use Data::Dumper;
use HTML::Entities;    # encode_entities

my %color_hash = (
	'red'    => 'red',
	'green'  => 'green',
	'orange' => 'orange',
	'yellow' => 'yellow',
	'blue'   => 'blue',
	'black'  => 'black',
	'white'  => 'white',
);

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
		[
			'tag' => "Analyse by a specific tag",
			{
				validate => 'Genomic/Tag',
				required => 1,
				default  => 'all'
			}
		],
		[
			'number_residues', 'Number the residues (every 10)',
			{ validate => 'Flag', }
		],
		[
			'residues_per_line', '# of residues per line, defaults to 70',
			{
				required => 1,
				min => 10,
				max => 1000,
				default  => 70,
				validate => 'Int',
			}
		],
		[],
		[
			'group', 'Name for a group of residues, free text',
			{
				required => 1,
				validate => 'String',
				multiple => 1,
				default =>
				  [ 'polar', 'hydrophobic', 'charged' ],
			}
		],
		[
			'match',
			'String of letters that should be matched',
			{
				required => 1,
				validate => 'String',
				multiple => 1,
				default  => [ 'HSTQNCY', 'AVLIMPFWG', 'ERDK' ]
			}
		],
		[
			'fg',
			'Foreground colour for a group',
			{
				required => 1,
				validate => 'Option',
				multiple => 1,
				default  => [ 'white', 'black', 'white' ],
				options  => \%color_hash,
			}
		],
		[
			'bg',
			'Background colour for a group',
			{
				required => 1,
				validate => 'Option',
				multiple => 1,
				default  => [ 'red', 'green', 'blue' ],
				options  => \%color_hash,
			}
		],
		['nocli', 'Do not output using the terminal, only save to the file' ],
	],

	'outputs' => [
		[
			'charges',
			'Charges Output',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'charges',
				data_format    => 'text/html',
				default_format => 'HTML'
			}
		],
	],
	'defaults' => [
		'appid'   => 'Charges',
		'appname' => 'Charges',
		'appdesc' => 'colour sequences based on rules.',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => 'Default',
			params => {
				'file' => 'test-data/inputs/multi.gbk',
				'tag' => 'CDS',
				'nocli' => '',
			},
			outputs      => {
				'charges' => ['charges.html', 'test-data/outputs/charges.gbk.html'],
			}
		},
		{
			test_name => 'from FA',
			params => {
				'file' => 'test-data/inputs/multi.cds.fa',
				'tag' => 'CDS',
				'nocli' => '',
			},
			outputs      => {
				'charges' => ['charges.html', 'test-data/outputs/charges.fa.html'],
			}
		},
	],
);

# The colours are actually more complicated than that.
# This is the standard colour attributes (0-7)
# There are also 'bright_' versions of each.

# However, for 256 color terminals, you can also use
# ansi0-ansi15, grey0-grey23
# rgbRGB where R,G,B \in [0,5]

# We can do some fun things here, (finally!).  We'll always produce an HTML
# file for posterity However, if we're running under CLI (!opts->galaxy) we can
# use coloured terminal output! :3
my %grouped;
#my @group = @{$options->{group}};
my @match = @{ $options->{match} };
my @fg    = @{ $options->{fg} };
my @bg    = @{ $options->{bg} };

if ( scalar(@match) == scalar(@fg) && scalar(@fg) == scalar(@bg) ) {
	use List::MoreUtils qw(each_array);
	#use Digest::Crc32;
	#my $crc = Digest::Crc32->new();
	my $ea = each_array( @match, @fg, @bg );
	while ( my ( $a, $b, $c ) = $ea->() ) {
		$a =~ s/[^A-Za-z]*//g;
		my $id = 'group_' . $a; #$crc->strcrc32($a);
		$grouped{$a} = {
			match  => qr/[$a]/i,
			fg     => $b,
			bg     => $c,
			css_id => $id,
			css    => sprintf(
				".%s{\nbackground: %s; color: %s;\n}\n",
				$id, $c, $b
			),
		};
	}
	my %args = (
		'file'      => $options->file,
		'callback'  => \&func,
		'translate' => 1,
		'header'    => 1,
		'subset'    => $options->{'tag'}
	);
	use CPT::Bio;
	my $cptbio = CPT::Bio->new();
	$cptbio->parseFile(%args);
}
else {
	die
'The size of match, foreground, and background were not identical';
}

my $html_output;

sub func {
	my $response_ref = shift;
	my @response     = @{$response_ref};

	$html_output .= sprintf( '
<style type="text/css">
	%s
</style>
<ul class="list">
	%s 
</ul>',
		join( "\n", map { $grouped{$_}{'css'} } keys(%grouped) )
		,    # Join the pregen'd css
		join(
			"\n",
			map {
"<li style=\"list-style:none\"> [$_]<span style=\"float:left;width:20px;color:"
				  . $grouped{$_}{fg}
				  . "\">&nbsp;</span><span style=\"float:left;width:20px;background:"
				  . $grouped{$_}{bg}
				  . "\">&nbsp;</span></li>"
			  } keys(%grouped)
		),    #create the items/coloured spans
	);

	# Running under the galaxy environment
	if ( $options->{galaxy} ) {
		# Nothing special
	}
	else {
		require Term::ANSIColor;
	}
	foreach my $row (@response) {
		my ( $header, $sequence ) = @{$row};

		# Escape the fasta leader sequence.
		$html_output .= '<h3>' . encode_entities($header) . '</h3>';
		$html_output .=
		  '<pre>' . chargedAAnew( $sequence, 'html' ) . '</pre>';
		if ( !$options->{galaxy} && !$options->{nocli} ) {
			printf( "%s: \n", $header );
			print chargedAAcli($sequence);
			print "\n";
		}
	}

	use CPT::OutputFiles;
	my $output = CPT::OutputFiles->new(
		name => 'charges',
		libCPT => $libCPT,
	);
	$output->CRR(data => $html_output);

}

sub wrapAAcli {
	use Term::ANSIColor;
	my ( $character, $type ) = @_;
	foreach ( keys(%grouped) ) {
		if ( $character =~ $grouped{$_}{match} ) {
			print color sprintf( "%s on_%s",
				$grouped{$_}{fg}, $grouped{$_}{bg} );
			last;
		}
	}
	print $character;
	print color 'reset';
}

sub wrapAA {
	my ( $character, $type ) = @_;
	my $return;
	foreach ( keys(%grouped) ) {
		if ( $character =~ $grouped{$_}{match} ) {
			$return = sprintf( '<span class="%s">%s</span>',
				$grouped{$_}{css_id}, $character );
			last;
		}
	}
	if(!defined($return)){
		$return = $character;
	}
	return $return;
}

sub chargedAAnew {
	my ( $orig_seq, $type ) = @_;

	#Strip all the whitespace
	$orig_seq =~ s/\s*//g;
	my @aa_seq = split( //, $orig_seq );

	my @charges = ();
	foreach my $char (@aa_seq) {
		if ( $char =~ /[kr]/i ) {
			push( @charges, '+' );
		}
		elsif ( $char =~ /[de]/i ) {
			push( @charges, '-' );
		}
		else {
			push( @charges, ' ' );
		}
	}

	# Calculate the number of lines needed
	my $num_of_lines =
	  int( ( scalar @aa_seq ) / $options->{'residues_per_line'} ) + 1;

	my $response = '';
	for ( my $j = 0 ; $j < $num_of_lines ; $j++ ) {

		#print out the charge
		my $line_charges;
		my $line_residues;
		my $line_numbers = "";
		for (
			my $k = ( $j * $options->{residues_per_line} ) ;
			$k < ( ( $j + 1 ) * $options->{residues_per_line} )
			&& $k < scalar(@aa_seq) ;
			$k++
		  )
		{
			$line_charges .= $charges[$k];
			$line_residues .= wrapAA( $aa_seq[$k], $type );

		}
		if ( $options->{number_residues} ) {
			for (
				my $k = ( $j * $options->{residues_per_line} ) ;
				$k <
				( ( $j + 1 ) * $options->{residues_per_line} )
				&& $k < scalar(@aa_seq) ;
				$k += 10
			  )
			{
				$line_numbers .= sprintf( "%10s", $k + 10);
			}
			$line_numbers .= "\n";
		}
		if(!$line_charges){$line_charges = "";}
		if(!$line_residues){$line_residues = "";}
		$response .= "$line_charges\n$line_residues\n$line_numbers\n\n";
	}
	return $response;
}

sub chargedAAcli {
	my ( $orig_seq, $type ) = @_;

	#Strip all the whitespace
	$orig_seq =~ s/\s*//g;
	my @aa_seq = split( //, $orig_seq );

	my @charges = ();
	foreach my $char (@aa_seq) {
		if ( $char =~ /[kr]/i ) {
			push( @charges, '+' );
		}
		elsif ( $char =~ /[de]/i ) {
			push( @charges, '-' );
		}
		else {
			push( @charges, ' ' );
		}
	}

	# Calculate the number of lines needed
	my $num_of_lines =
	  int( ( scalar @aa_seq ) / $options->{'residues_per_line'} ) + 1;

	my $response = '';
	for ( my $j = 0 ; $j < $num_of_lines ; $j++ ) {

		#print out the charge
		my $line_residues;
		my $line_numbers = "";
		print "\t";
		for (
			my $k = ( $j * $options->{residues_per_line} ) ;
			$k < ( ( $j + 1 ) * $options->{residues_per_line} )
			&& $k < scalar(@aa_seq) ;
			$k++
		  )
		{
			print $charges[$k];
		}
		print "\n\t";
		for (
			my $k = ( $j * $options->{residues_per_line} ) ;
			$k < ( ( $j + 1 ) * $options->{residues_per_line} )
			&& $k < scalar(@aa_seq) ;
			$k++
		  )
		{
			wrapAAcli( $aa_seq[$k], $type );
		}
		print "\n";
		if ( $options->{number_residues} ) {
			print "\t";
			for (
				my $k = ( $j * $options->{residues_per_line} ) ;
				$k <
				( ( $j + 1 ) * $options->{residues_per_line} )
				&& $k < scalar(@aa_seq) ;
				$k += 10
			  )
			{
				printf( "%10s", $k + 10);
			}
			print "\n";
		}
	}
}

=head1 DESCRIPTION

Charges simply pretty-prints fasta sequences with their charges labelled above the sequence

=cut
