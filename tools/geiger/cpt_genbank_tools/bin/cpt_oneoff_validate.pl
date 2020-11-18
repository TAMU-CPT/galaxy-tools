#!/usr/bin/perl

use strict;
use warnings;
use autodie;
use CPT::GalaxyGetOpt;
use Data::Dumper;
use Text::CSV;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[ 'file'       , 'Input file Spanin CSV Database'                  , { required => 1 , validate => 'File/Input' } ] ,
		[ 'rz'         , 'Column of the Rz sequence (numbers start at 1)'  , { required => 1 , validate => 'Int' , default => 1,} ]  ,
		[ 'rz1'        , 'Column of the Rz1 sequence (numbers start at 1)' , { required => 1 , validate => 'Int' , default => 2,} ]  ,
		[ 'has_header' , 'Should we ignore the first line' ]               ,

	],
	'outputs' => [
		[
			'results',
			'Report of Valid/Invalid spanins',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'spanin_report',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'OneOff.Rohit.Spanins',
		'appname' => 'Spanin validator',
		'appdesc' => 'Validates spanins from Rohit\'s Excel database',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

my ($rows_ref,$header_ref) = read_csv_file($options->{file}, headers => $options->{has_header});

my @rows = @{$rows_ref};
my @headers = @{$header_ref};

sub read_csv_file {
	my ($file, %attr) = @_;
	my @headers;
	my @rows;
	my $csv = Text::CSV->new( { binary => 1 } )    # should set binary attribute.
	  or die "Cannot use CSV: " . Text::CSV->error_diag();
	open my $fh, "<:encoding(utf8)", $file;
	if ( $attr{headers} ) {
		push( @headers, $csv->getline($fh) );
	}
	while ( my $row = $csv->getline($fh) ) {
		push @rows, $row;
	}
	$csv->eof or $csv->error_diag();
	close $fh;
	return (\@rows,\@headers);
}

use CPT::External::TMHMM;
use CPT::External::LipoP;

push( @{ $headers[0] }, 'TMDs', 'Lipo' );

for ( my $i = 0 ; $i < scalar @rows ; $i++ ) {
	# pass I/O through TMHMM/LipoP respectively, report fails.
	my @row = @{ $rows[$i] };
	my ( $spanin_a, $spanin_b ) = ( $row[ $options->{rz} -1 ], $row[ $options->{rz1} -1 ] );
	my $tmhmm = CPT::External::TMHMM->new();
	my $lipop = CPT::External::LipoP->new();
	$tmhmm->analyze_seq($spanin_a);
	$lipop->analyze_seq($spanin_b);
	my $tmds     = $tmhmm->num_predicted();
	my $cleavage = $lipop->cleavage();
	if ( $tmds > 0 ) {
		push @row, $tmds;
	}
	else {
		push @row, 'NO TMDS FOUND';
	}
	if ( defined $cleavage ) {
		push @row, join( "-", @{$cleavage} );
	}
	else {
		push @row, 'NO CLEAVAGE';
	}
	$rows[$i] = \@row;
}

my %results = (
	'Sheet1' => {
		'header' => \@{ $headers[0] },
		'data'   => \@rows,
	}
);

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'results',
	GGO => $ggo,
);
$crr_output->CRR(data => \%results);

#foreach(my $row, @rows){
# DB Columns of I/O/E, validate against
# Check if overlapping.
# the script would check the beginning nt of the o-spanin and see if
# > it is INSIDE the orf of the upstream i-spanin.  If so, it is either Embedded
# > or Overlapped.  In this case, if the last nt of the o-spanin is INSIDE  the
# > i-spanin, it is EMBEDDED.  Else it is OVERLAPPED.  If the start nt of the
# > o-spanin is outside the ORF of the i-spanin, it is SEPARATED.

# Report on
#> Actually, it would be good to generate a report that gave numbers about the
#> relative position of the start and stop codons of the o-spanin are compared
#> to the start and stop codons of the i-spanin.   Might make a sentence or a
#> table or figure in the paper.
#}
