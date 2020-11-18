#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
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
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/multi.cds2.fa',
			},
			outputs => {
				'my_output_data' => ["translated.fa", 'test-data/outputs/translated.fa' ],
			},
		},
	],
);
use CPT::BioData;
use CPT::Bio;
my $bd = CPT::BioData->new();
my $codonTable = $bd->getTranslationTable($options->{tn_table});
my $bio = CPT::Bio->new(codonTable => $codonTable);

my $seqio = Bio::SeqIO->new(
    -file   => $options->file,
    -format => 'fasta',
);

my $output = '';
while ( my $seqobj = $seqio->next_seq() ) {
    my $seq_id = $bio->_getIdentifier($seqobj);
    my $seq_prop = $seqobj->description();
    my $seq = $bio->intelligent_get_seq($seqobj);
	$seq = $bio->translate($seq);
    if($options->{strip_stops}){
        $seq =~ s/\*.*//g;
        $seq =~ s/\+.*//g;
        $seq =~ s/#.*//g;
    }else{
        $seq =~ s/\*.*/*/g;
        $seq =~ s/\+.*/+/g;
        $seq =~ s/#.*/#/g;
    }
    if(length($seq) > 0 ){
        $output .= sprintf(">%s %s\n%s\n", $seq_id, $seq_prop, $seq);
    }
}

use CPT::OutputFiles;
my $psmout = CPT::OutputFiles->new(
    name => 'fasta',
    GGO => $ggo,
);
$psmout->CRR(data => $output);

=head1 DESCRIPTION

Translate sequences from DNA to AAs.

=cut
