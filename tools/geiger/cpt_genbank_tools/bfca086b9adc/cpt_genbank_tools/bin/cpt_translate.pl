#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#

use strict;
use warnings;

use CPT;
use Data::Dumper;

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file',
			'Input Fasta file',
			{
				required => 1,
				validate => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			'tn_table',
			'Translation Table',
			{
				required => 1,
				validate => 'Option',
				options => {
					"1" => "The Standard Code",
					"2" => "The Vertebrate Mitochondrial Code",
					"3" => "The Yeast Mitochondrial Code",
					"4" => "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
					"5" => "The Invertebrate Mitochondrial Code ",
					"6" => "The Ciliate, Dasycladacean and Hexamita Nuclear Code",
					"9" => "The Echinoderm and Flatworm Mitochondrial Code",
					"10" => "The Euplotid Nuclear Code",
					"11" => "The Bacterial, Archaeal and Plant Plastid Code",
					"12" => "The Alternative Yeast Nuclear Code",
					"13" => "The Ascidian Mitochondrial Code",
					"14" => "The Alternative Flatworm Mitochondrial Code",
					"15" => "Blepharisma Nuclear Code",
					"16" => "Chlorophycean Mitochondrial Code",
					"21" => "Trematode Mitochondrial Code",
					"22" => "Scenedesmus Obliquus Mitochondrial Code",
					"23" => "Thraustochytrium Mitochondrial Code",
					"24" => "Pterobranchia Mitochondrial Code",
					"25" => "Candidate Division SR1 and Gracilibacteria Code",
				},
				default => '11',
			},
		],
		[
			'strip_stops',
			'Remove stop codons from translated sequence',
			{ validate => 'Flag' },
		],
	],
	'outputs' => [
		[
			'fasta',
			'Translated fasta file',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'translated',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
	],
	'defaults' => [
		'appid'   => 'FastaTranslator',
		'appname' => 'Translate Fasta Sequeces',
		'appvers' => '1.94',
		'appdesc' => 'dna to protein with selectable coding table',
	],
	'tests' => [
	],
);
use CPT::Bio;
my $bio = CPT::Bio->new();
$bio->set_codon_table($options->{tn_table});

my %args = (
	'file'     => $options->file,
	'callback' => \&func,
	'header'   => 1,
	'subset'   => 'CDS',
);

$bio->parseFile(%args);

my $results;
sub func {
	my ($response_ref) = @_;
	my @response = @{$response_ref};
	foreach(@response){
		my ($header, $seq) = @{$_};
		$seq = $bio->translate($seq);
		if($options->{strip_stops}){
			$seq =~ s/\*//g;
			$seq =~ s/\+//g;
			$seq =~ s/#//g;
		}
		$results .= "$header\n$seq\n";
	}
}

use CPT::OutputFiles;
my $psmout = CPT::OutputFiles->new(
	name => 'fasta',
	libCPT => $libCPT,
);
$psmout->CRR(data => $results);

=head1 DESCRIPTION

Translate sequences from DNA to AAs.

=cut
