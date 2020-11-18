#!/usr/bin/env perl
use CPT::GalaxyGetOpt;
use Data::Dumper;
use CPT::GBK2GFF3;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
		[ 'is_circular' => 'Set this flag if the genome is circular' ],
		[
			'seqid',
			'Seqid used in GFF3 file',
			{ required => 1, validate => 'String' }
		],
		[
			'id_prefix',
'While the ID must be unique only within the scope of the GFF3 file, by the GFF3 spec, Chado treats them as globally unique',
			{ required => 1, validate => 'String' }
		],
		[
			'cpt_modifications',
'Coerce all transformation specialised for the CPT to Genbank (mostly specification of programmes we use, e.g., transterm, etc)',
		],
	],
	'outputs' => [
		[
			'gff3',
			'Converted GFF3 File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'converted',
				data_format    => 'genomic/interval',
				default_format => 'GFF3',
			}
		],
	],
	'defaults' => [
		'appid'   => 'gbk2gff3',
		'appname' => 'Genbank To GFF3',
		'appvers' => '1.94',
		'appdesc' => 'conversion utility. ',
	],
	'tests' => [
	],
);

my $seqid = $options->{'seqid'};
$seqid =~ s/ /_/g;
my $gff = CPT::GBK2GFF3->new(
	seqid               => $seqid,
	id_prefix           => $options->{'id_prefix'},
	annotation_software => {
		'source'     => 'Newbler',
		'CDS'        => 'GenemarkS',
		'tRNA'       => 'tRNAScan-SE',
		'terminator' => 'TranstermHP',
	},
	override_source => 'Genbank',
);
$gff->header(['##gff-version 3']);

use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqio = $bio->getSeqIO($options->{file});

while ( my $seqobj = $seqio->next_seq() ) {
	$gff->sequence( $seqobj->seq() );
	my @features = $seqobj->get_SeqFeatures();
	foreach my $feat (@features) {
		$gff->add_feature($feat);
	}
}

my $output = $gff->get_gff3_file();
use CPT::OutputFiles;
my $output2 = CPT::OutputFiles->new(
        name => 'gff3',
        GGO => $ggo,
);
$output2->CRR(data => $output);

=head1 DESCRIPTION

Convert GBK files to GFF3 files. This is a B<ALPHA> quality tool, so beware!

=cut
