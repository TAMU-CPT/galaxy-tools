#!/usr/bin/perl

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
			{ multiple => 1, validate => 'String' }
		],
		[
			'accession_file',
			'File containing a list of newline separated accession numbers',
			{ multiple => 1, validate => 'File/Input' }
		],
		[
			'split',
			'If multiple files are downloaded, should they be split',
			{ validate => 'Flag' }
		],
		[
			'email',
			'User email address (for NCBI accounting purposes)',
			{ validate => 'String', required => 1, default => '__user_email__', _galaxy_specific => 1, hidden => 1 }
		]
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

if(!$options->{galaxy} && $options->{email} eq '__user_email__'){
	die "Sorry about this, due to limitations of our library, you've not been informed about a required parameter, --email. Please specify --email with your email to let NCBI know who is making the request";
}

my @accessions = ();

if(defined($options->{accession})){
	push(@accessions, map { split/[,\s]+/}@{ $options->{accession} });
}
if(defined($options->{accession_file})){
	foreach my $af(@{$options->{accession_file}}){
		open(my $fh, '<', $af);
		while(<$fh>){
			chomp;
			push(@accessions, $_);
		}
		close($fh);
	}
}


if(scalar @accessions == 0){
	die 'You must specify at least one accession number';
}

use CPT::OutputFiles;
my $gbk_output = CPT::OutputFiles->new(
	name => 'genbank',
	GGO => $ggo,
);

# Batch them into N sequence jobs
my $iterations = (scalar @accessions) /500;
for(my $i = 0; $i < $iterations; $i++){
	print STDERR "On $i of $iterations\n";
	use Bio::DB::GenBank;
	my $gb = Bio::DB::GenBank->new(
		-retrievaltype => 'tempfile',
		-email         => 'cpt@tamu.edu'
	);
	my @acc_subset = @accessions[500*$i..(500*$i+500)];
	my $seqio      = $gb->get_Stream_by_acc( \@acc_subset );

	if(!$options->{split}){
		$gbk_output->varCRR(data => $seqio, filename => "batch_$i");
	}
	else {
		while(my $seqobj = $seqio->next_seq){
			# Wow this is simple.
			my $name = $seqobj->accession_number();
			$name =~ s/_//g;
			$gbk_output->varCRR(data => $seqobj, filename => $name);
		}
	}
}


=head1 NAME

Batch Genbank Downloader

=head1 DESCRIPTION

Allows for downloading of multiple genbank files from NCBI by accession numbers. This script produces a single genbank file with all of your records concatenated, which is allowed per specification.

Please note that it is unlikely CPT scripts in genbank will behave nicely with this sort of input (most are written around reading standalone genbank files). As of now (2014-03) I am working on fixing this from both ends: a script to split genbank files, and better handing of genbank database.

=cut
