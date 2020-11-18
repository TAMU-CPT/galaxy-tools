#!/usr/bin/perl
use strict;
use warnings;

use CPT;
use Data::Dumper;

my $phage_in_middle         = qr/^(?<host>.*)\s*phage (?<phage>.*)$/;
my $bacteriophage_in_middle = qr/^(?<host>.*)\s*bacteriophage (?<phage>.*)$/;
my $starts_with_phage =
  qr/^(bacterio|vibrio|Bacterio|Vibrio)?[Pp]hage (?<phage>.*)$/;
my $new_style_names = qr/(?<phage>v[A-Z]_[A-Z][a-z]{2}[A-Z]_.*)/;

use CPT::GalaxyGetOpt;
use Bio::SeqIO;
my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[
			'file', 'Input Genbank file',
			{
				validate => 'File/Input',
				required => 1,
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
	],
	'outputs' => [
		[
			'output',
			'Output Genome File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'out',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PhageRenamer',
		'appname' => 'Phage Renamer',
		'appdesc' => 'renames files according to specific rules',
		'appvers' => '1.96.0',
	],
);
use CPT::OutputFiles;
my $gbk_output = CPT::OutputFiles->new(
	name => 'output',
	GGO => $ggo,
);

my $seq_in = Bio::SeqIO->new(
	-file   => $options->{file},
	-format => 'Genbank',
);
my $genome = $seq_in->next_seq;

my $acc         = $genome->accession_number;
my $species     = $genome->species();
my $name        = $species->node_name();
my $parsed_name = name_parser($name);
my $output_filename;
if ( defined $parsed_name ) {
	$output_filename = $parsed_name;
}
else {
	print STDERR "Could not detect from $name\n";
	$output_filename = $acc;
}
$output_filename =~ s/\s/_/g;

$gbk_output->varCRR(data => $genome, filename => $output_filename);

sub name_parser {
	my ($parse_against) = @_;
	my ( $host, $phage );
	if ( $parse_against =~ $phage_in_middle ) {
		$host  = $+{host};
		$phage = $+{phage};
	}
	elsif ( $parse_against =~ $bacteriophage_in_middle ) {
		$host  = $+{host};
		$phage = $+{phage};
	}
	elsif ( $parse_against =~ $starts_with_phage ) {
		$phage = $+{phage};
	}
	elsif ( $parse_against =~ $new_style_names ) {
		$phage = $+{phage};
	}
	if ( defined $host ) {
		return "$phage [$host]";
	}
	else {
		return $phage;
	}
	return
}
