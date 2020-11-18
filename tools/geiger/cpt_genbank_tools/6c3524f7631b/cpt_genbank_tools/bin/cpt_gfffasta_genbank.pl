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
use Data::Dumper;

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'fasta',
			'Input Fasta file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			'gff3',
			'Input GFF3 file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['gff3'],
			},
		],
	],
	'outputs' => [
		[
			'genbank',
			'Merged Genbank File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'merged.gbk',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank'
			}
		],
	],
	'defaults' => [
		'appid'   => 'GfftoGenbank',
		'appname' => 'GFF3 + Fasta to Genbank',
		'appvers' => '1.94',
		'appdesc' => '',
	],
	'tests' => [
	],
);
use Bio::Tools::GFF;
use Bio::SeqIO;

#Valid keys were extracted by grepping the genbank manual page basically for everything that went /[^=] and then sorting/uniqing, etc until we got this list.
my $tags = 'allele anticodon artificial_location
bio_material bound_moiety cell_line cell_type chromosome citation
clone clone_lib codon_start collected_by collection_date compare
country cultivar culture_collection db_xref dev_stage direction
EC_number ecotype environmental_sample estimated_length exception
experiment focus frequency function gap_type gene gene_synonym
germline haplogroup haplotype host identified_by inference isolate
isolation_source lab_host lat_lon linkage_evidence locus_tag
macronuclear map mating_type mobile_element_type mod_base mol_type
ncRNA_class note number old_locus_tag operon organelle organism
partial PCR_conditions PCR_primers phenotype plasmid pop_variant
product protein_id proviral pseudo rearranged replace
ribosomal_slippage rpt_family rpt_type rpt_unit_range rpt_unit_seq
satellite segment serotype serovar sex specimen_voucher
standard_name strain sub_clone sub_species sub_strain tag_peptide
tissue_lib tissue_type transgenic translation transl_except
transl_table trans_splicing variety';
my %valid_tags = map { $_ => 1 } split( /\s+/, $tags );

my $keys = "-10_signal -35_signal 3'UTR 5'UTR
CAAT_signal CDS C_region D-loop D_segment GC_signal J_segment LTR
N_region RBS STS S_region TATA_signal V_region V_segment
assembly_gap attenuator enhancer exon gap gene iDNA intron mRNA
mat_peptide misc_RNA misc_binding misc_difference misc_feature
misc_recomb misc_signal misc_structure mobile_element
modified_base ncRNA old_sequence operon oriT polyA_signal
polyA_site precursor_RNA prim_transcript primer_bind promoter
protein_bind rRNA rep_origin repeat_region sig_peptide source
stem_loop tRNA terminator tmRNA transit_peptide unsure variation";
my %valid_keys = map { $_ => 1 } split( /\s+/, $keys );

#These tags have equivalents that should just be transformed. When in doubt make it a note.
my %to_convert = (

	#tags
	'dbxref'   => 'db_xref',
	'parent'   => 'note',
	'name'     => 'label',
	'Name'     => 'label',
	'old-name' => 'obsolete_name',
	'id'       => 'note',
	'nat-host' => 'host',
	'genome'   => 'note'
	, #This sounds important. Perhaps should look into what it's supposed to be.
	  #keys
	'region' => 'source'
);
my $seq = Bio::SeqIO->new(
	-file   => $options->{fasta},
	-format => 'fasta'
);
my $seq_out = Bio::Seq->new( -seq => $seq->next_seq->seq );

my $gff = new Bio::Tools::GFF( -file => $options->{gff3}, -gff_version => 3 );
while ( my $feature = $gff->next_feature ) {
	my $primary_tag = $feature->primary_tag;
	if ( !defined( $valid_keys{$primary_tag} )
		|| $to_convert{$primary_tag} )
	{
		if ( $to_convert{$primary_tag} ) {
			$primary_tag = $to_convert{$primary_tag};
		}
		elsif ( $feature->has_tag('gbkey')
			&& ( $primary_tag && not $to_convert{$primary_tag} ) )
		{
			$primary_tag =
			  join( '', $feature->get_tag_values('gbkey') );

			#}elsif($valid_keys{$primary_tag}){
			#	$primary_tag = $1;
		}
		else {
			print STDERR
			  "WARN: Invalid GBK primary tag [$primary_tag]\n";
		}
	}
	if ( $feature->has_tag('gbkey') ) {
		$feature->remove_tag('gbkey');
	}
	$feature->set_attributes( -primary => $primary_tag );
	foreach my $tag_id ( $feature->get_all_tags ) {
		if ( !defined( $valid_tags{$tag_id} ) ) {
			if ( $to_convert{ lc($tag_id) } ) {
				$feature->add_tag_value(
					$to_convert{ lc($tag_id) },
					$feature->get_tag_values($tag_id)
				);

		#}elsif($valid_tags =~ /\b($tag_id)\b/i){
		#	$feature->add_tag_value($1,$feature->get_tag_values($tag_id));
			}
			else {
				print STDERR "WARN: Invalid key id [$tag_id]\n";
			}
			$feature->remove_tag($tag_id);
		}
	}
	$seq_out->add_SeqFeature($feature);
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'genbank',
	libCPT => $libCPT,
);
$crr_output->CRR(data => $seq_out);

=head1 DESCRIPTION

Given a GFF3 set of annotations and a fasta file, this tool will create a single genbank file out of them

=cut
