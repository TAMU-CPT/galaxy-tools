#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:   404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#

use strict;
use warnings;

use CPT;
use File::Temp qw/ tempfile /;
use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();
use IPC::Run qw(run timeout);

my $libCPT = CPT->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file|f', 'Input file',
			{ required => 1, validate => 'File/Input' }
		],
		[
			'minimum_confidence',
			'Minimum confidence level to include',
			{ required => 1, validate => 'Float', default => 80.0 }
		],
		[],
		['Transterm: Scoring Parameters'],
		[
			'gc',
			'Score of a G-C pair',
			{ required => 1, validate => 'Float', default => -2.3 }
		],
		[
			'au',
			'Score of an A-U pair',
			{ required => 1, validate => 'Float', default => -0.9 }
		],
		[
			'gu',
			'Score of a G-U pair',
			{ required => 1, validate => 'Float', default => 1.3 }
		],
		[
			'mm',
			'Score of any other pair',
			{ required => 1, validate => 'Float', default => 3.5 }
		],
		[
			'gap',
			'Score of a gap in the hairpin',
			{ required => 1, validate => 'Int', default => 6 }
		],
		[
			'max_hp_score',
			'Maximum allowable hairpin score',
			{ required => 1, validate => 'Float', default => -2 }
		],
		[
			'max_tail_score',
			'Maximum allowable tail score',
			{ required => 1, validate => 'Float', default => -2.5 }
		],
		[
			'loop_penalty',
'The cost of loops of various lengths can be set using --loop_penalty=f1,f2,f3,f4,f5,...fn, where f1 is the cost of a loop of length --min_loop, f2 is the cost of a loop of length --min_loop+1, as so on. If there are too few terms to cover up to max_loop, the last term is repeated.',
			{
				required => 1,
				validate => 'String',
				default  => '1,2,3,4,5,6,7,8,9,10,11'
			}
		],
		['Transterm: Length Parameters'],
		[
			'max_len',
			'Total extent of hairpin <= n NT long',
			{ required => 1, validate => 'Int', default => 59 }
		],
		[
			'min_stem',
			'Stem must be n nucleotides long',
			{ required => 1, validate => 'Int', default => 4 }
		],
		[
			'max_loop',
			'The loop portion can be no longer than n',
			{ required => 1, validate => 'Int', default => 13 }
		],
		[
			'min_loop',
			'Loop portion of the hairpin must be at least n long',
			{ required => 1, validate => 'Int', default => 3 }
		],

		['Transterm: U-rich region parameters'],
		[
			'uwin_require',
'Number of "U" nucleotides in the --uwin_length long region.',
			{ required => 1, validate => 'Int', default => 3 }
		],
	],
	'outputs' => [
		[
			'gff3',
			'GFF3 File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'transterm_gff',
				data_format    => 'text/plain',
				default_format => 'TXT',
			}
		],
		[
			'csv_report',
			'Report of Transterm HP run',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'transterm',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			}
		],
	],
	'defaults' => [
		'appid'   => 'Transterm',
		'appname' => 'Transterm HP',
		'appvers' => '1.94',
		'appdesc' =>
		  'Add terminators to your annotated genomes via Transterm HP',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/single.gbk',
			},
			outputs => {
				'genbank' => ['transterm.gbk', 'test-data/outputs/transterm.gbk'],
				'csv_report' => ['transterm.Sheet1.csv', 'test-data/outputs/transterm.csv'],
			}
		},
	],
);

#  ['start-cut=i'                   , '???'                                                 , { required => 1 , validate => 'Int'} ] ,
#  ['end-cut=i'                     , '???'                                                 , { required => 1 , validate => 'Int'} ] ,

use CPT::Bio;
my $bio = CPT::Bio->new();
my %args   = ( 'file' => $options->file, );
my $seqobj = ${ $bio->requestCopy(%args) };

# Export the genome as FASTA
my $genome_export_results_ref = $bio->parseFile(
	'file'   => $options->file,
	'header' => 1,
	'subset' => 'whole',
);

my @results = @{$genome_export_results_ref};
my ( $fh, $fasta_path ) = tempfile(
	'transterm_tmp_XXXXXX',
	UNLINK => 1,
	DIR    => '/tmp/',
	SUFFIX => '.fasta'
);
my $sequence_header_id = substr( $results[0][0], 1 );
print $fh ">$sequence_header_id\n";
my $seq = $results[0][1];
$seq =~ s/(.{60})/$1\n/g;
print $fh $seq . "\n";
close($fh);

# Create the coordinates file
my ( $coords_fh, $coords_path ) = tempfile(
	'transterm_tmp_XXXXXX',
	UNLINK => 1,
	DIR    => '/tmp/',
	SUFFIX => '.coords'
);
foreach my $feat($seqobj->get_SeqFeatures) {
	if ( $feat->primary_tag eq 'CDS' ) {
		print $coords_fh join( " ",
			'CDS', $feat->start, $feat->end, $sequence_header_id )
		  . "\n";
	}
}
close($coords_fh);
# Run Transterm
my @command = (
	"transterm",                          '-p',
	'/usr/share/transtermhp/expterm.dat', '--all-context'
);
foreach (
	qw(gc au gu mm gap max_hp_score max_tail_score loop_penalty max_len min_stem max_loop min_loop uwin_require)
  )
{
	my $transterm_flag = $_;
	$transterm_flag =~ s/_/-/g;

	if ( $options->{$_} ) {
		push( @command, "--$transterm_flag=" . $options->{$_} );
	}
}

push( @command, $fasta_path, $coords_path );

#print join(" ",@command),"\n";

my ( $in, $out, $err );
run \@command, \$in, \$out, \$err, timeout(10) or warn "Transterm: $?";
my @lines = split( "\n", $out );

my $has_started_results = 0;

my $regex_term =
qr/  TERM \d+ \s* (\d+) - (\d+)\s*(\+|-) [^ ]* \s* (\d+)\s* ([0-9.-]+)\s* -([0-9.]+)\s*\|\s*(.*)/;
my $parts_term = qr/  ([^ ]*)\s+([^ ]*)\s+([^ ]*)\s+([^ ]*)\s+([^ ]*)/;

my %terminators = (
	'Sheet1' => {
		'header' => [
			"Start",    "End",
			"Strand",   "Confidence",
			"HP Score", "Tail Score",
			"Notes",    "5' Tail",
			"5' Stem",  "Loop",
			"3' Stem",  "3' Loop"
		],
		'data' => [],
	},
);
my @hits;

for ( my $i = 0 ; $i < scalar @lines ; $i++ ) {
	my $line = $lines[$i];
	if ( substr( $line, 0, 8 ) eq 'SEQUENCE' ) {
		$has_started_results = 1;
	}
	if ( $has_started_results && substr( $line, 0, 2 ) eq '  ' ) {
		my @hit_item;
		if ( $line =~ $regex_term ) {
			push( @hit_item, $1, $2, $3, $4, $5, $6, $7 );
		}
		if ( $lines[ $i + 1 ] =~ $parts_term ) {
			push( @hit_item, $1, $2, $3, $4, $5 );
		}
		push( @{ $terminators{'Sheet1'}{'data'} }, \@hit_item );
		push( @hits,                               \@hit_item );
		$i++;
	}
}

# TODO This should be moved to an OOP module but w/e.
my @gff3 = ('##gff-version 3');
my $i = 0;
foreach (@hits) {
	my @hit = @{$_};
	print join( ",", @hit ), "\n" if $options->{verbose};

#5363  , 5382 , +      , 91         , -11.4    , 4.03204    , opp_overlap 5362 5357 , AGATGTAAACACACT , GGCTCACC , TTAG , GGTGGGCC , TTTCTGCGTTTAATA
#start , end  , strand , confidence , hp score , tail score , notes                 , 5' tail         , 5'stem   , loop , 3' stem  , 3'loop
# 0    , 1    , 2      , 3          , 4        , 5          , 6                     , 7               , 8        , 9    , 10       , 11
    my @gff3_line = (
        $sequence_header_id,
        'TransTermHP',
        'terminator',
        $hit[0],
        $hit[1],
        $hit[3],
        $hit[2],
        '.',
        sprintf('ID="terminator_%s";Note="%s"', $i, $hit[6]),
    );
    $i++;
    use Data::Dumper;
    push(@gff3, join("\t", @gff3_line));
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'csv_report',
	GGO => $ggo,
);
$crr_output->CRR(data => \%terminators);

my $crr_output2 = CPT::OutputFiles->new(
	name => 'gff3',
	GGO => $ggo,
);
$crr_output2->CRR(data => join("\n",@gff3));
