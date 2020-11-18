use Bio::Chado::Schema;
use Bio::Perl;
use strict;
use warnings;
use Bio::Tools::CodonTable;

# PODNAME: extract_by_query.pl

my $codon_table = Bio::Tools::CodonTable->new(-id=>11);

use CPT;
my $libCPT = CPT->new();
use CPT::GalaxyGetOpt;
my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[
			'database',
			'Database Name',
			{
				required => 1,
				validate => 'String'
			}
		],
		['query','Search string for searching through all tags in the database. Use "*" as a wildcard character, e.g., "*ISP*" will match anything with the characters ISP in it', { required => 1, validate => 'String'}],
		['translate','Translate to amino acid sequence (uses t/n table #11)'],
	],
	'outputs' => [
		[
			'sequence',
			'Extracted Features',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'query',
				data_format    => 'genomic/raw',
				default_format => 'Fasta'
			}
		],
	],
	'defaults' => [
		'appid'   => 'CHAPLIN_extract_fasta',
		'appname' => 'Extract Fasta from Query',
		'appdesc' => 'allows exporting a query against chado as fasta.',
	]
);

$options->{query} =~ s/\*/%/g;


my $dsn = "dbi:Pg:dbname=" . $options->{database} . ";host=cpt.tamu.edu;port=5432;sslmode=require";
my $user = "charm_admin";
my $password = "oNFkI0KyoGygRp8Zf7jOVIrR1VmsOWak";
my $chado = Bio::Chado::Schema->connect( $dsn, $user, $password );

my $results = $chado->resultset('Sequence::Featureprop')->search(
	{ value => { like => $options->{query} } },
	{ join => ['cvterm', 'feature'] }
);

my @data;
my @feature_ids;
my @fasta_sequences;
while(my $row = $results->next){
	# If it's a primary sequence, skip
	#if(defined $row->feature->seqlen && $row->feature->seqlen > 0){
		#print STDERR "Hit a primary sequence\t";
		#print STDERR $row->feature->feature_id;
		#print STDERR "\n";
		#next;
	#}

	# Grab the first referenced featureloc
	my @srcfeat = $row->feature->featureloc_features;
	#
	foreach(@srcfeat){
		my $srcfeat_id = $_->srcfeature_id;
		my $subquery = $chado->resultset('Sequence::Feature')->search(
			{ feature_id => $srcfeat_id }
		);
		my ($left, $right, $strand) = ($_->fmin, $_->fmax, $_->strand);
		while(my $subhit = $subquery->next){
			my $seq;
			if($strand == -1){
				$seq = substr($subhit->seq(),$left,($right-$left));
				$seq = revcom($seq)->seq;
			}else{
				$seq = substr($subhit->seq(),$left,($right-$left));
			}
			if($options->{translate}){
				$seq = $codon_table->translate($seq);
			}
			$seq =~ s/\*$//g;
			push(@data,
				{
					organsim      => $row->feature->organism->common_name,
					parent_length => $subhit->seqlen,
					tag           => $row->feature->uniquename,
					left          => $left,
					right         => $right,
					strand        => $strand,
					length        => length($seq),
					seq           => $seq,
					query_hit     => $row->value,
				}
			);

			my $tag =  $row->feature->organism->common_name . "_" . $options->{query};
			$tag =~ s/[^A-Za-z0-9_:.-]*//g;
			push(@fasta_sequences,
				sprintf(">%s [id=%s;length=%s;strand=%s;left=%s;right=%s;query=%s]\n%s", $tag, $row->feature->uniquename,length($seq),$strand,$left,$right,$row->value,$seq));
		}
	}
}
if($options->{verbose}){
	use Data::Format::Pretty::Console qw(format_pretty);
	print format_pretty(\@data);
}
#open(my $output,'>', 'out.fa');
#print $output join("\n",@fasta_sequences);
#close($output);
use CPT::OutputFiles;
my $data_out = CPT::OutputFiles->new(
	name   => 'sequence',
	GGO => $ggo,
);
$data_out->CRR(data => join("\n", @fasta_sequences));
