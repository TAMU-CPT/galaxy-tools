#!/usr/bin/env perl
use strict;
use warnings;
use Text::Levenshtein qw/distance/;
use Carp;
use CPT::GalaxyGetOpt;

# PODNAME: phantasm_compare_cids.pl

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file', 'Genome CID list',
			{
				validate => 'File/Input',
				required => 1,
			}
		],
	],
	'outputs' => [
		[
			'comp_cid',
			'Scored CID pairs',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'comp_cids',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PHAnTASM.Comparison.CIDS',
		'appname' => 'CID Scoring',
		'appdesc' => 'produces from-to-score table of CIDS for input CID table',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

my %map;

# Load in data
open(my $cids_fh, '<', $options->{file});
while(<$cids_fh>){
	chomp;
	# Ignore commented lines
	next if ($_ =~ /^#/);
	# Remove quotes if there are any
	s/"//g;
	# Split by commas and tabs
	my ($seq, $cid) = split(/,|\t/,$_);
	$map{$seq} = revcomrot($cid);
}
close($cids_fh);

# Compare N vs N
my @ks = sort(keys(%map));
my $i = 0;
my @table_main;
foreach my $from(@ks){
	foreach my $to(@ks){
		if($i % 100 == 0){
			print STDERR (100 * $i/(scalar(@ks) * scalar(@ks))),"\n";
		}
		$i++;
		my $c = cidscore($from,$to);
		if($c == 0){$c = '0.0';}
		push(@table_main, [$from, $to, $c]);
	}
}

# Save Data
my %data = (
        'Sheet1' => {
		header => ['From', 'To', 'Score'],
		data => \@table_main,
        }
);
use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
        name => 'comp_cid',
        GGO => $ggo,
);
$csv_output->CRR(data => \%data);




sub cidscore {
	my ($from, $to) = @_;
	my @f = @{$map{$from}};
	my @t = @{$map{$to}};
	# Some high number
	my $minscore = 1_000;
	foreach(@t){
		my $p = scorepair($f[0],$_);
		if(abs($p) < $minscore){
			$minscore = $p;
		}
	}
	return $minscore;
}

sub scorepair {
	my ($a, $b) = @_;
	return distance($a,$b);
}

sub revcomrot {
	my ($cid) = @_;
	# The very first array element is actually our original! easy peasy
	return [rot($cid), rot(revcom($cid))];
}
sub revcom {
	my ($cid) = @_;
	#reverse
	my @b = reverse(unpack("(A2)*", $cid));
	my $result = join('',@b);
	# Complement
	$result =~ y/+-/-+/;
	return $result;
}
sub rot {
	my ($cid) = @_;
	my @result = ();
	for(my $i=0;$i<length($cid)/2;$i++){
		push(@result, substr($cid,$i,length($cid)) . substr($cid,0,$i));
	}
	return @result;
}





=head1 NAME

PHAnTASM CID Scoring Tool

=head1 DESCRIPTION

Given a list of CIDS

    #Genome	CID
    Lambda	c+m+r-r-m+m+
    Rogue	r-m+m+c+m+r-

This tool will score all vs. all pairings. The scoring process involves taking a query CID (let's use Lambda's in this example) and a subject CID (Rogue's). The subject CID is reversed and complemented

    r-m+m+c+m+r-
    # reversed + comp'd:
    r+m-c-m-m-r+

Then with both of these strings, they're rotated through all possible rotations. The trivial example of the string C<abcd> would produce

    abcd
    bcda
    cdab
    dabc

With a long list of reversed+complemented and rotations of these possibilities, we score our subject and queries with a simple calculation of Levenshtein distance. The minimum score "wins", and the others are discarded. This process is done to account for the possibility the genome was opened in a different location than our query, but still has an identical layout of genes.

=cut
