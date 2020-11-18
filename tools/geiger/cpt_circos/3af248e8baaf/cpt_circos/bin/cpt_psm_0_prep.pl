#!/usr/bin/env perl
use strict;
use warnings;
use Storable;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Data::Dumper;
use File::Temp qw/tempfile tempdir/;


use CPT;
my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[ 'file', 'Input file', { validate => 'File/Input', multiple => 1, required => 1} ],
	],
	'outputs' => [
		[
			'cpt_psm_object',
			'Output PSM Object',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'cpt_psm',
				data_format    => 'Dummy',
				default_format => 'Dummy'
			}
		],
	],
	'defaults' => [
		'appid'   => 'PSM.Prep',
		'appname' => 'PSM Prep',
		'appdesc' => 'prepares data for the PSM Plotter',
		'appvers' => '1.94.2',
	],
	'tests' => [],
);



use CPT::Bio;
my $bio = CPT::Bio->new();

my @genbank_files = @{$options->{file}};

my %data = (
	file_list => [],
);

my $GLOBAL_UNLINK_VAR = 1;
my $tempdir = tempdir('cpt.psm2.XXXXXXX',CLEANUP => $GLOBAL_UNLINK_VAR);

foreach my $file(@genbank_files){
	my $seqio_object = Bio::SeqIO->new(-file => $file,-format=>'genbank');
	while(my $seqobj = $seqio_object->next_seq){
		my ( $fh, $path ) = tempfile('cds_export.XXXXXXXXX', UNLINK => $GLOBAL_UNLINK_VAR, DIR => $tempdir, SUFFIX => '.fa');

		my @gi_array;
		foreach my $feat ( $seqobj->get_SeqFeatures ) {
			if($feat->primary_tag eq 'CDS'){
				my $header = $bio->_getIdentifier($feat);
				my $seq = $bio->translate(
					$bio->intelligent_get_seq($feat));

				# Proteins come with translated stop codon
				$seq =~ s/\*//g;
				$seq =~ s/\+//g;
				$seq =~ s/#//g;
				
				push @gi_array, $header;
				print $fh ">$header\n$seq\n";
			}
		}
		$data{gbk}{$seqobj->display_id()}{'gi'} = \@gi_array;
		$data{gbk}{$seqobj->display_id()}{'fasta_location'} = $path;
		$data{gbk}{$seqobj->display_id()}{'gbk_location'} = $file;
		push(@{$data{file_list}}, $seqobj->display_id());
		close $fh;
	}
}



use IPC::Run3;

# Concatenate Fasta Files
my @fasta_files;
foreach(@{$data{file_list}}){
	push(@fasta_files, $data{gbk}{$_}{fasta_location});
}
my @command = ('cat', @fasta_files);
my ($merged_fa_fh, $merged_fa_path) = tempfile('merged.XXXXXXXXX', UNLINK => 1, DIR => $tempdir, SUFFIX => '.fa');
my ($in, $out, $err);
run3 \@command, \$in, \$out, \$err;
if($err){
	print STDERR $err;
}
print $merged_fa_fh $out;
close($merged_fa_fh);


# Create Blast Database
my ($blastdb_fh, $blastdb_path) = tempfile('blastdb.XXXXXXXXX', UNLINK => 1, DIR => $tempdir);
@command = ('makeblastdb', 
	'-dbtype', 'prot',
	'-in', $merged_fa_path,
	'-out', $blastdb_path,
);
my ($makeblast_in,$makeblast_out,$makeblast_err);
run3 \@command, \$makeblast_in, \$makeblast_out, \$makeblast_err;

# Blast files
foreach(@{$data{file_list}}){
	#push(@fasta_files, $data{gbk}{$_}{fasta_location});
	my @blast_cmd = ('blastp',
		'-query', $data{gbk}{$_}{fasta_location},
		'-out', $data{gbk}{$_}{fasta_location} . ".xml",
		'-outfmt', '5',
		'-db', $blastdb_path,
	);
	my ($blast_in,$blast_out,$blast_err);
	run3 \@blast_cmd, \$blast_in, \$blast_out, \$blast_err;
}

my $value;
my @data_tsv;

foreach(@{$data{file_list}}){
	#push(@fasta_files, $data{gbk}{$_}{fasta_location});
	my $file = $data{gbk}{$_}{fasta_location};

	my $in = new Bio::SearchIO(
		-format   => 'blastxml',
		-tempfile => 1,
		-file     => "$file.xml",
	);
	while( my $result = $in->next_result ) {
		while( my $hit = $result->next_hit ) {
			while( my $hsp = $hit->next_hsp ) {
				my $Identity = $hsp->percent_identity/100 * $hsp->length('query');
				my $IterationQueryLength = $result->query_length();
				my $HitLength = $hit->length();
				my $dice = (200*$Identity)/($HitLength + $IterationQueryLength);

				my $c = $result->query_description();#genome_a_header
				my $d = $hit->name();#genome_b_header
				# Skip self-self links
				next if($c eq $d);

				push (@data_tsv,
					[
						$c,
						$d,
						$hsp->evalue(),
						$dice,
					]
				);

			}#close hsp
		}#close hit
	}#close result
}#close for genomes

$data{hit_table} = \@data_tsv;

use CPT::OutputFiles;
my $psmout = CPT::OutputFiles->new(
	name => 'cpt_psm_object',
	libCPT => $libCPT,
);
my $psm_loc = $psmout->CRR(data => "none", extension => 'psm');
store \%data,$psm_loc;







=head1 NAME

PSM Prep

=head1 DESCRIPTION

This tool takes in 2 or more GenBank files, blasts, and prepares data structures for use in the companion tool: PSM Plotter. Select as many (multi)-gbk files as you I<might> want to plot. Once this tool is done, you can select any subset of those to plot then.

=cut
