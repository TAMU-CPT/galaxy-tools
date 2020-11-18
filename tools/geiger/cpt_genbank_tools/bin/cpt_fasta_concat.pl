#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;
use Bio::SeqIO;
use CPT::Bio;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
				multiple => 1,
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
		[
			'reverse', 'Sequence should
            be reversed',
			{
				required => 1,
				validate => 'Option',
				options  => { 'yes' => 'yes', 'no' => 'no' },
				multiple => 1
			}
		],
	],
	'outputs' => [
		[
			'merged',
			'merged fasta sequences',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'merged',
				data_format    => 'genomic/raw',
				default_format => 'Fasta'
			}
		],
		[
			'rejects',
			'non-merged fasta sequences',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'rejects',
				data_format    => 'genomic/raw',
				default_format => 'Fasta'
			}
		],
	],
	'defaults' => [
		'appid'   => 'fasta_concat',
		'appname' => 'Fasta concat',
		'appvers' => '1.94',
		'appdesc' =>
'concatenates FASTA files. Requires two (or more) inputs fasta files with identical numbers of sequences, and pairs them off going down the list.',
	],
	'tests' => [
		#{
			#test_name => "Normal+Reversed",
			#command_line =>
#"--file t/test-files/proteins2.fa --reverse no --file t/test-files/proteins2.fa --reverse yes",
			#outputs => {
				#"merged.fa" =>
				  #"a536f9e3d35463f3e7d525267becaa08",
			#}
		#},
		#{
			#test_name => "Normal+Normal",
			#command_line =>
#"--file t/test-files/proteins2.fa --reverse no --file t/test-files/proteins2.fa --reverse no",
			#outputs => {
				#"merged.fa" =>
				  #"bea65b31b50d78476016d2bc6d79df52",
			#}
		#},
		#{
			#test_name => "Reversed+Normal",
			#command_line =>
#"--file t/test-files/proteins2.fa --reverse yes --file t/test-files/proteins2.fa --reverse no",
			#outputs => {
				#"merged.fa" =>
				  #"7322f2e44cd6d816957115e4a25c0712",
			#}
		#},
		#{
			#test_name => "Reversed+Reversed",
			#command_line =>
#"--file t/test-files/proteins2.fa --reverse yes --file t/test-files/proteins2.fa --reverse yes",
			#outputs => {
				#"merged.fa" =>
				  #"f484e9959c105318fc5afd7ebfb78ab3",
			#}
		#},
	],
);

my $bio = CPT::Bio->new();
my @reversals;
foreach my $reverse ( @{ $options->{reverse} } ) {
	push( @reversals, ( $reverse eq 'yes' ? 1 : 0 ) );
}

my %feature_set;
my $file_idx = 0;
# For each file in our file set
my @orig_key_ordering;
foreach my $file ( @{ $options->{file} } ) {
	my $seqio = $bio->getSeqIO($file);
	# Read in all the fasta sequences
	while ( my $seqobj = $seqio->next_seq() ) {
		my $seq = $seqobj->seq();
		my $display_id = $seqobj->display_id();
		# Reverse if need be
		if ( $reversals[$file_idx] ) {
			$seq = reverse($seq);
		}
		# Remove stop codons
		$seq =~ s/\*//g;
		$seq =~ s/\+//g;
		$seq =~ s/#//g;
		# Remove some specific strings (BAD BAD BAD)
		$display_id =~ s/_outside_\d+$//g;
		$display_id =~ s/_second//g;

		#unless ( $feature_set{ $display_id } ) {
			#$feature_set{ $display_id } = '>' . $display_id . "\n";
		#}
		push(@orig_key_ordering, $display_id);
		push(@{$feature_set{ $display_id }},$seq);
	}
	$file_idx++;
}

# Sets of FASTA features
my @good;
my @rejected;

foreach(@orig_key_ordering) {
	my @seqs = @{$feature_set{$_}};
	my $seqobj = Bio::Seq->new(
		-display_id => $_,
		-seq => join('', @seqs),
	);
	if(scalar(@seqs) > 1){
		push(@good, $seqobj);
	}else{
		push(@rejected, $seqobj);
	}
}

use CPT::OutputFiles;
my $crroutput = CPT::OutputFiles->new(
        name => 'merged',
        GGO => $ggo,
);
$crroutput->CRR(data => \@good);

my $rejected_output = CPT::OutputFiles->new(
        name => 'rejects',
        GGO => $ggo,
);
$rejected_output->CRR(data => \@rejected);

=head1 DESCRIPTION

This tool will attempt to concatenate multiple fasta files by finding sequences with B<IDENTICAL> fasta IDs in each file. Sequences with a partner will be placed in the "merged" file, and those without partners will be placed in the rejects file

=cut
