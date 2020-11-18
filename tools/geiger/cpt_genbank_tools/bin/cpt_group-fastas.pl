#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Fasta input file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			"center_around", "Center
			around which residue? This will occur for the first
			time this residue is seen in a given sequence",
			{ validate => 'String', required => 1 }
		],
		[
			"number_before",
			"How many resiudes before the target should be
		considered?",
			{
				validate => 'Int',
				required => 1,
				default  => 4
			}
		],
	],
	'outputs' => [
		[
			'report',
			'Grouped Sequence Report',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'report',
				data_format    => 'text/html',
				default_format => 'HTML',
			}
		],
	],
	'defaults' => [
		'appid'   => 'GroupSeqs',
		'appname' => 'Group Sequences',
		'appvers' => '1.94',
		'appdesc' => 'group sequences according to user queries',
	],
	'tests' => [
	],
);
use Bio::SeqIO;
use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqio = $bio->getSeqIO($options->{file});

my %data;
my $n            = $options->{number_before};
my $c            = $options->{center_around};
my $reg          = qr/(.{$n})$c/;
my $max_location = -1;
while ( my $seq = $seqio->next_seq ) {
	my ( $loc, $before );
	if ( $seq->seq() =~ $reg ) {
		$loc    = $-[0] + $n;
		$before = $1;
		if ( $loc > $max_location ) {
			$max_location = $loc;
		}
	}
	else {
		$before = "";
		$loc    = 0;
	}

	$data{ $seq->display_id() } = {
		seq => $seq->seq(),
		loc => $loc,
		bef => $before,
	};
}

#use CPT::Report;
#use CPT::Report::HTML;
#my $now_string = localtime;
#my $r          = CPT::Report::HTML->new(
#author => 'galaxy',
#date   => $now_string,
#title  => $options->{"appname"} . " Report"
#);

my $report     = "<!DOCTYPE html>\n<html><head></head><body>";
my $is_new_key = "";
foreach ( sort { $data{$b}{bef} cmp $data{$a}{bef} } keys(%data) ) {
	if ( $data{$_}{bef} ne $is_new_key ) {
		$report .= sprintf( '<h1>%s</h1>', $data{$_}{bef} );
		$is_new_key = $data{$_}{bef};
	}
	my $addition = $max_location - $data{$_}{loc};
	$report .= sprintf( '<p>%s</p>', ">$_\n" . $data{$_}{seq} . "\n" );
}
$report .= "</body></html>";

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'report',
	GGO => $ggo,
);
$crr_output->CRR(data => $report);

=head1 DESCRIPTION

Given a set of fasta sequences, this tool will group them acording to the first time it encounters a set of N residues which it should C<center_around>. All the sequences with that cluster of amino acids are then reported

=cut
