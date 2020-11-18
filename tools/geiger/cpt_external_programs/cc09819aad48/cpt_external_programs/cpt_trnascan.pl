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
use File::Temp;
use CPT;
use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my %code = (
	'universal'    => "Universal",
	'ciliate'      => "Ciliate Nuclear",
	'vertebrate'   => "Vertebrate Mito",
	'invertebrate' => "Invertebrate Mito",
	'yeast'        => "Yeast Mito",
	'echinoderm'   => "Echinoderm Mito",
	'mold'         => "Mold and Protozoan Mito"
);
my %eutrna = (
	'deafult' => "Default",
	'relaxed' => "Relaxed",
	'strict'  => "Strict (Pavesi params)",
);
my %search_mode = (
	'default'     => "Default",
	'cove'        => "Cove only (very slow)",
	'trna_only'   => "tRNAscan only",
	'eufind_only' => "EufindtRNA only",
	'eufind_cove' => "EufindtRNA then Cove",
	'trna_cove'   => "tRNAscan then Cove"
);
my %source = (
	'mixed'       => "Mixed (general tRNA model)",
	'euk'         => "Eukaryotic",
	'bact'        => "Bacterial",
	'arch'        => "Archaeal",
	'mito_chloro' => "Mito and Chloroplast"
);

my $options = $ggo->getOptions(
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
			'disable_pseudo',
			'Disable pseudo gene checking',
			{ validate => 'Flag' }
		],
		[ 'cove_cutoff', 'Cove score cutoff', { validate => 'Float' } ],
		[
			'intermediate_cutoff',
			'Intermediate score cutoff',
			{ validate => 'Float' }
		],
		[
			'eutrna',
			'EufindtRNA search parameters',
			{
				validate => 'Option',
				options  => \%eutrna,
				required => 1,
				default  => 'default',
			}
		],
		[
			'genetic_code',
			'Genetic Code for tRNA Isotype Prediction',
			{
				validate => 'Option',
				options  => \%code,
				required => 1,
				default  => 'universal',
			}
		],
		[
			'search_mode',
			'Search Mode',
			{
				validate => 'Option',
				options  => \%search_mode,
				required => 1,
				default  => 'default',
			}
		],
		[
			'source', 'Source',
			{
				validate => 'Option',
				options  => \%source,
				required => 1,
				default  => 'bact',
			}
		],
	],
	'outputs' => [
		[
			'genbank',
			'Genbank file with annotated tRNAs',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'added_trnas',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
		[
			'report',
			'tRNAscanSE report',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'transcan',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'tRNAscanSE',
		'appname' => 'tRNAscanSE',
		'appvers' => '1.94',
		'appdesc' => 'adds tRNAs to your GBK files',
	],
	'tests' => [
	]
);

my $temp_fasta = File::Temp->new(
	TEMPLATE => 'galaxy.trnascan.XXXXX',
	DIR      => '/tmp/',
	SUFFIX   => '.fa'
);
my %args = (
	'file'     => $options->file,
	'callback' => \&func,
	'subset'   => 'whole',
	'header'   => 1,
);
use CPT::Bio;
my $cptbio = CPT::Bio->new();
$cptbio->parseFile(%args);

sub func {
	my $response_ref      = shift;
	my @response          = @{$response_ref};
	my $complete_response = "";
	foreach (@response) {
		my ( $header, $sequence, $other ) = @{$_};
		$sequence =~ s/(.{80})/$1\n/g;
		$complete_response .= "$header\n$sequence\n";
	}
	print $temp_fasta $complete_response;
}

my %opts = (
	'pseudo'  => '',
	'covecut' => '',
	'inter'   => '',
);
my %source_lookup = (
	"Mixed (general tRNA model)" => "-G",
	"Eukaryotic"                 => "",
	"Bacterial"                  => "-B",
	"Archaeal"                   => "-A",
	"Mito/Chloroplast"           => "-O"
);
my %search_mode_lookup = (
	"Default"               => "",
	"Cove only (very slow)" => "-C",
	"tRNAscan only"         => "-T",
	"EufindtRNA only"       => "-E",
	"EufindtRNA -> Cove"    => "-E -C",
	"tRNAscan -> Cove"      => "-T -C",
);

my @command = ();
if ( $options->{disable_pseudo} ) {
	push( @command, $opts{'pseudo'} );
}
if ( $options->{cove_cutoff} ) {
	push( @command, $opts{'covecut'}, $options->{cove_cutoff} );
}
if ( $options->{intermediate_cutoff} ) {
	push( @command, $opts{'inter'}, $options->{intermediate_cutoff} );
}
if ( $options->{search_mode} ) {
	my $param = $search_mode_lookup{ $options->{search_mode} };
	if ($param) {
		push( @command, $param );
	}
}
if ( $options->{source} ) {
	my $param = $source_lookup{ $options->{source} };
	if ($param) {
		push( @command, $param );
	}
}

push( @command, '-b' );
push( @command, $temp_fasta );
use IPC::System::Simple qw(capture);
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };
$ENV{'PATH'} = '/bin:/usr/bin:/usr/local/bin';
my $path     = $ENV{'PATH'};                         # $path now NOT tainted
my $captured = capture( "tRNAscan-SE", @command );
my @tRNA     = ();
my @data     = ();
foreach ( split( /\n/, $captured ) ) {
	chomp;
	my @row = split( /\t/, $_ );

	#  tRNAs found on the reverse (lower) strand are indicated by having
	#  the Begin (5') bound greater than the End (3') bound.
	my $c = $row[5];
	$c =~ tr/T/U/;
	my $strand;
	if ( $row[2] > $row[3] ) {
		$strand = -1;
	}
	else {
		$strand = 1;
	}
	push(
		@tRNA,
		{
			start            => $row[3],
			end              => $row[2],
			strand           => $strand,
			trna_type        => $row[4],
			codon_recognized => $c,
			cove_score       => $row[8],
		}
	);
	push( @data, [ $row[3], $row[2], $strand, $row[4], $c, $row[8] ] );

	# The  'tRNA  Type'  is the predicted amino acid charged to the tRNA
	# molecule based on the predicted anticodon (written 5'->3') displayed
	# in the next column.   tRNAs that fit criteria for potential
	# pseudogenes (poor primary or secondary structure), will be marked
	# with "Pseudo" inthe 'tRNA Type' column (pseudogene checking is
	# further discussed in theMethods section of the program manual).  If
	# there is a predicted intronin the tRNA, the next  two columns
	# indicate the nucleotide bounds.  Ifthere is no predicted intron, both
	# of these columns contain zero.

	# The  final  column  is  the  Cove score for the tRNA in bits of
	# information.  Specifically, it is a log-odds score: the log of the
	# ratio of the probability of the sequence given the tRNA covariance
	# model used (developed fromhand-alignment of 1415 tRNAs), and the
	# probability of the sequence given a simple random sequence model.
	# tRNAscan-SE counts any sequence that attains a score of 20.0 bits or
	# larger as a tRNA (based on empirical studies con‐ducted by Eddy &
	# Durbin in ref #2).

}
my %results = (
	'tRNAs' => {
		header =>
		  [ 'Start', 'End', 'Strand', 'tRNA Type', 'Cove Score' ],
		data => \@data,
	}
);
use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqobj = ${ $bio->requestCopy(%args) };
foreach (@tRNA) {
	my %tRNA = %{$_};
	my $feat = new Bio::SeqFeature::Generic(
		-start       => $tRNA{start},
		-end         => $tRNA{end},
		-strand      => 0,
		-primary_tag => 'tRNA',
		-tag         => {
			product    => "tRNA-" . $tRNA{trna_type},
			"evidence" => 'tRNAscan-SE',
			"note"     => [
				"codon_recognized=" . $tRNA{codon_recognized},
				"Cove Score: " . $tRNA{cove_score},
			],
		}
	);
	$seqobj->add_SeqFeature($feat);
}
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'genbank',
	GGO => $ggo,
);
$crr_output->CRR(data => $seqobj);

$crr_output = CPT::OutputFiles->new(
	name => 'report',
	GGO => $ggo,
);
$crr_output->CRR(data => \%results);
