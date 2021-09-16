#!/usr/bin/env perl
use strict;
use warnings;
use Storable;
use CPT::GalaxyGetOpt;
use CPT::Bio::NW_MSA;
use Data::Dumper;
use CPT::Circos::Conf;
use POSIX;


my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file', 'PSM2 Data File', { validate => 'File/Input', required => 1 } ],
		[],
		['Cutoffs'],
		['evalue' , 'Evalue cutoff' , { validate => 'Float' , default => 1e-4 } ] ,
		['dice'   , 'Dice cutoff'   , { validate => 'Float' , default => 50 } ]   ,
		[],
		['Alignment Options'],
		['mismatch'    , 'Mismatch Score' , { validate => 'Float' , default => -1} ] ,
		['gap_penalty' , 'Gap Penalty'    , { validate => 'Float' , default => '0.0' } ] ,
		['match'       , 'Match Score'    , { validate => 'Float' , default => 5 } ]  ,
	],
	'outputs' => [
		[
			'diff_table',
			'Output Comparison Table',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'genome_comp',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			},
        ],
        [
			'blastclust',
			'Output Blastclust Table',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'blastclust',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PSM.Comp',
		'appname' => 'PSM Comparison Table',
		'appdesc' => 'aligns and lists data from PSM Prep',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);


my %data_file = %{retrieve($options->{file})};

print STDERR "Aliging genomes\n";
my $msa = CPT::Bio::NW_MSA->new(
	gap_penalty => $options->{'gap_penalty'},
	match_score => $options->{'match'},
	mismatch_score => $options->{'mismatch'},
	bidi => 1,
);

my @hits = @{$data_file{hit_table}};
my @clusters;

foreach my $hit(@hits){
	my ($from, $to, $evalue, $dice) = @{$hit};
	if($evalue < $options->{evalue} && $dice > $options->{dice}){
		if($options->{verbose}){
			print "$from $to\n";
		}

        my $foundmatch = 0;
        foreach my $cluster(@clusters){
            if($from ~~ @{$cluster} || $to ~~ @{$cluster}){
                $foundmatch = 1;
                if(!($from ~~ @{$cluster})){
                    push(@{$cluster}, $from);
                }
                if(!($to ~~ @{$cluster})){
                    push(@{$cluster}, $to);
                }
            }
        }
        if($foundmatch == 0){
            push(@clusters, ["".($#clusters+2), $from, $to]);
        }
		$msa->add_relationship($from, $to);
	}
}
my @fixed_clusters;
foreach my $cluster (@clusters) {
    my ($idx, @values) = @{$cluster};
    push(@fixed_clusters, [$idx, join(',', @values)]);
}

my @user_ordering = keys($data_file{gbk});

foreach my $genome(@user_ordering){
	print STDERR "\tAligning $genome\n";
	my $gi_list_ref = $data_file{gbk}{$genome}{gi};#"GI" list
	if(! defined $gi_list_ref){
		warn "Could not find $genome genome in the data file. Please be sure you have correctly specified the name of a genome from a genbank file. (See the LOCUS line for the name).";
	}else{
		$msa->align_list($gi_list_ref);
	}
}

my @aligned_results = $msa->merged_array();
# Remove CRC64 hashes from sequences
foreach my $row(@aligned_results){
	$row = [map { s/_[A-F0-9]{16}$//; $_ } @{$row}];
    #my $key = ${$row}[0];
    #foreach my $cluster(@clusters){

    #}
}

my %table = (
	'Sheet1' => {
		header => \@user_ordering,
		data =>  \@aligned_results,
	}
);


use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'diff_table',
	GGO => $ggo,
);
$crr_output->CRR(data => \%table);

my %table2 = (
    'Sheet1' => {
        header => ['Cluster ID', 'Contents'],
        data => \@fixed_clusters,
    }
);
my $crr_output2 = CPT::OutputFiles->new(
	name => 'blastclust',
	GGO => $ggo,
);
$crr_output2->CRR(data => \%table2);

=head1 DESCRIPTION

Following the execution of the PSM Prep tool, this tool simply aligns the genomes and generates a table comparison the positions of all proteins. It can be very useful to figure out which genes are missing in which genomes.

=head2 IMPORTANT PARAMETERS

=over 4

=item C<mismatch>, C<gap_penalty>, C<match>

These parameters control the Needleman-Wunsch Multiple Sequence Alignment library's scoring scheme. Mismatch scores are generally negative and discourage unrelated proteins from being plotted in a line together. Match scores encourage related proteins to line up. Gap penalty is set at zero as we generally prefer gaps to mismatches in this tool; phage genomes are small and gaps are "cheap" to use, whereas mismatches can sometimes give an incorrect impression of relatedness. That said, how your plots look is completely up to you and we encourage experimentation!

=back

=cut
