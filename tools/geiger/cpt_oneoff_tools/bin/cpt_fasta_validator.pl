#!/usr/bin/perl
#
#	   Code written by Eric Rasche
#			   mailto:rasche.eric@yandex.ru
#			   tel:404.692.2048
#			   http://eric.rasche.co.uk
#	   for
#			   Center for Phage Technology
#

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Digest::Crc32;
use Bio::Tools::SeqStats;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file|f' => 'Input file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
	],
	'outputs' => [
		[
			'report',
			'Validation Report',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'report',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'FastaValidator',
		'appname' => 'FASTA Validator',
		'appvers' => '1.94',
		'appdesc' =>
'checks FASTA files for formatting, duplicate sequences, and duplicate IDs',
	],
	'tests' => [
	],
);

my $crc = new Digest::Crc32();

# Load the fasta
open( my $fasta, '<', $options->{file} );

my ( $current_id, $current_annot, $line_number, $current_seq,
	$line_no_started_on );
$line_number = 0;
my %data;
my %crc32_seq_lookup;
my @errors;
while (<$fasta>) {
	$line_number++;
	chomp $_;
	if ( $_ =~ /^>([^ ]+)(.*)/ ) {
		if ( not defined($current_id) ) {
			$current_id    = $1;
			$current_annot = $2;
			$current_annot =~ s/^\s+|\s+$//g;
		}
		else {

			# Store last sequence
			# done in else, so as to not be called very first time
			store(
				id    => $current_id,
				seq   => $current_seq,
				note  => $current_annot,
				first => $line_no_started_on,
			);
		}
		$line_no_started_on = $line_number;

		# Start the new one
		$current_id    = $1;
		$current_annot = $2;
		$current_annot =~ s/^\s+|\s+$//g;
		$current_seq = "";
	}
	else {
		$current_seq .= $_;
	}
}
close($fasta);
store(
	id    => $current_id,
	seq   => $current_seq,
	note  => $current_annot,
	first => $line_no_started_on,
);

sub store {
	my (%d) = @_;
	if ( $data{ $d{id} } ) {
		push( @errors, [ 'duplicate_id', $d{id} ] );
	}
	my $hash = $crc->strcrc32( $d{seq} );
	if ( $crc32_seq_lookup{$hash} ) {
		push( @errors, [ 'duplicate_seq', $hash ] );
	}

	# Store everything with that hash
	push( @{ $crc32_seq_lookup{$hash} }, $d{id} );

	# store the sequence itself
	push(
		@{ $data{ $d{id} } },
		{
			id    => $d{id},
			seq   => $d{seq},
			note  => $d{note},
			first => $d{first},
			hash  => $hash,
		}
	);
}
use Data::Dumper;

# Generate report
my $report = "";
for my $error (@errors) {
	my @err = @{$error};
	if ( $err[0] eq 'duplicate_id' ) {
		my @hits    = @{ $data{ $err[1] } };
		my $no_hits = scalar @hits;
		$report .= "ERROR: Duplicate ID found [$no_hits]\n";
		for my $hit (@hits) {
			my %z = %{$hit};
			$report .=
			  sprintf( "\t%s on line %s\n", $z{id}, $z{first} );
		}
	}
	elsif ( $err[0] eq 'duplicate_seq' ) {

		#my @hits = @{$data{$err[1]}};
		my @hits    = @{ $crc32_seq_lookup{ $err[1] } };
		my $no_hits = scalar @hits;
		$report .= "ERROR: Duplicate sequence found [$no_hits]\n";
		for my $hit (@hits) {
			my @correct_hits;
			foreach my $subhit ( @{ $data{$hit} } ) {
				if ( ${$subhit}{hash} eq $err[1] ) {
					push( @correct_hits, $subhit );
					$report .= sprintf(
"\t%s on line %s with sequence \"%s...\"\n",
						${$subhit}{id},
						${$subhit}{first},
						substr(
							${$subhit}{seq}, 0, 30
						)
					);
				}
			}
		}

	}
}
if ( length($report) > 1 ) {
	print STDERR $report;
}
my @rows = ();
foreach my $key ( sort keys(%data) ) {
	foreach my $subhit ( @{ $data{$key} } ) {
		my %hit = %{$subhit};
		push(
			@rows,
			[
				$hit{id},   $hit{note}, $hit{first},
				$hit{hash}, $hit{seq}
			]
		);
	}
}
my %result = (
	'FASTA' => {
		'header' => [
			'ID', 'Note', 'Line starts on', 'CRC32 Digest',
			'Sequence'
		],
		'data' => \@rows,
	},
);

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'report',
        GGO => $ggo,
);
$output->CRR(data => \%result);

=head1 DESCRIPTION

This tool attempts to "validate" fasta files, i.e., identify redundant sequences, bad fasta IDs, duplicate fasta IDs, etc.

=cut
