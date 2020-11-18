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
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
		[ 'atg', 'ATG is a valid start' ],
		[ 'ctg', 'CTG is a valid start' ],
		[ 'ttg', 'TTG is a valid start' ],
		[ 'gtg', 'GTG is a valid start' ],
		[
			'min_gene_length',
			'Minimum gene length.',
			{
				required => 1,
				validate => 'Int',
				default  => 30,
			}
		],
		['translate', 'Translate from DNA to proteins', { validate => 'Flag' }],
	],
	'outputs' => [
		[
			'fasta',
			'Found ORFs',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'orfs',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
	],
	'defaults' => [
		'appid'   => 'NaiveOrfFinder',
		'appname' => 'Naïve Orf Finder',
		'appvers' => '1.94',
		'appdesc' =>
'finds orf finders via a very basic start/stop analysis. No further processing is done. c.f. Emboss Sixpack',
	],
	'tests' => [
	],
);


use CPT::Bio;
my $bio = CPT::Bio->new();
my $seq_obj = ${ $bio->requestCopy( 'file' => $options->{file} ) };

use CPT::Bio::ORF;
my $orf = CPT::Bio::ORF->new(
	min_gene_length => $options->{min_gene_length},
	sc_atg => $options->{atg},
	sc_ctg => $options->{ctg},
	sc_ttg => $options->{ttg},
	sc_gtg => $options->{gtg},
);

my @seqs = $orf->run($seq_obj->seq());
#print Dumper map{ $_->translate() }@seqs;
if($options->{translate}){
	@seqs = map{ $_->translate() } @seqs;
}


use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'fasta',
	libCPT => $libCPT,
);
$crr_output->CRR(data => \@seqs);

=head1 DESCRIPTION

Basically this is EMBOSS six-pack, with the option of selecting your start codons

=cut
