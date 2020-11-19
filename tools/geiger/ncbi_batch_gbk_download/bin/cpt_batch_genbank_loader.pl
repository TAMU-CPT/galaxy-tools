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

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'accession|a',
			'Accession Number(s) to download. May be separated by space/comma/newline, or specified multiple times',
			{ multiple => 1, required => 0, validate => 'String' }
		],
		[
			'accession_file',
			'Accession list file',
			{ multiple => 0, required => 0, validate => 'File/Input' }
		],
		[
			'split',
			'If multiple files are downloaded, should they be split',
			{ validate => 'Flag' }
		],
	],
	'outputs' => [
		[
			'genbank',
			'Output Genbank Files',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'downloaded',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank'
			}
		],
	],
	'defaults' => [
		'appid'   => 'Batch_GBK_Downloader',
		'appname' => 'Batch Genbank Downloader',
		'appdesc' => 'Downloads genbank files from Genbank\'s servers',
		'appvers' => '1.94',
	]
);
use Bio::DB::GenBank;
my $gb = Bio::DB::GenBank->new(
	-retrievaltype => 'tempfile',
	-email         => 'cpt@tamu.edu'
);

my @accessions = map { split/[,\s]+/}@{ $options->{accession} };
if($options->{accession_file}){
    open(my $af_fh, '<', $options->{accession_file});
    foreach my $line(<$af_fh>){
        $line =~ s/[^A-Za-z0-9_-]//g;
        push(@accessions, $line);
    }
    close($af_fh);
}

my $seqio      = $gb->get_Stream_by_acc( \@accessions );


use CPT::OutputFiles;
my $gbk_output = CPT::OutputFiles->new(
	name => 'genbank',
	GGO => $ggo,
);

if(!$options->{split}){
	# Single output, could use varCRR to same effect.
	$gbk_output->CRR(data => $seqio);
}
else {
	while(my $seqobj = $seqio->next_seq){
		# Wow this is simple.
		my $name = $seqobj->accession_number();
		$name =~ s/_//g;
		$gbk_output->varCRR(data => $seqobj, filename => $name);
	}
}


=head1 NAME

Batch Genbank Downloader

=head1 DESCRIPTION

Allows for downloading of multiple genbank files from NCBI by accession numbers. This script produces a single genbank file with all of your records concatenated, which is allowed per specification.

Please note that it is unlikely CPT scripts in genbank will behave nicely with this sort of input (most are written around reading standalone genbank files). As of now (2014-03) I am working on fixing this from both ends: a script to split genbank files, and better handing of genbank database.

=cut
