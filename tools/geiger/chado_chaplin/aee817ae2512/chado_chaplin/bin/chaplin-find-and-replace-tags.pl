use Bio::Chado::Schema;
use strict;
use warnings;
# PODNAME: search-by-feature-tag.pl
use CPT;
my $libCPT = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[
			'database',
			'Database Name',
			{
				required => 1,
				validate => 'String'
			}
		],
		[
			'query',
			'SQL Query String used to select rows for doing replacement on (* is wildcard character)',
			{
				required => 1,
				validate => 'String'
			}
		],
		[
			'pcre_regex_from',
			'Regular expression used to match the target string inside the querried data. E.g., a value of "abc" would constitute part of the complete regex s/abc/???/g',
			{
				required => 1,
				validate => 'String'
			}
		],
		[
			'pcre_regex_to',
			'Portion of the regular expression representing the text replacement. E.g., a value of "def" would constitute part of the complete regex s/???/def/g',
			{
				required => 1,
				validate => 'String'
			}
		],
		[
			'noop',
			'Dry-run/No-op: no database modifications are made',
			{
				validate => 'Flag',
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Find and Replace Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'changes',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'CHAPLIN_search_by_feature_tag',
		'appname' => 'Search by Feature Tag',
		'appdesc' => 'lists features in a database resulting from a given search',
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

my %response = (
	'Sheet1' => {
		header => ['Organism ID', 'Feature Name', 'Tag', 'Value'],
		data => [],
	}
);

my @data;
push(@data, ['Organism', 'Feature Name', 'Tag', 'Original Value', 'New Value']);
my @feature_ids;
my $from = $options->{pcre_regex_from};
my $to = $options->{pcre_regex_to};
while(my $row = $results->next){
	my $val = $row->value;
	$val =~ s/$from/$to/gmx;
	push(@data,
		[
			$row->feature->organism_id,
			$row->feature->name,
			$row->cvterm->name,
			$row->value,
			$val,
		]
	);
	if(! defined $options->{noop}){
		$row->value($val);
		$row->update();
	}
}

$response{Sheet1}{data} = \@data;

use CPT::OutputFiles;
my $data_out = CPT::OutputFiles->new(
	name   => 'results',
	libCPT => $libCPT,
);
$data_out->CRR(data => \%response);
