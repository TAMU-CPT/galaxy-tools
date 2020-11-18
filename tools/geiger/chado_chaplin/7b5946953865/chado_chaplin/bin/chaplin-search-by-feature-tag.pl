use Bio::Chado::Schema;
use strict;
use warnings;
use Data::Format::Pretty::Console qw(format_pretty);
# PODNAME: search-by-feature-tag.pl
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
		[
			'query',
			'Query String (* is wildcard character)',
			{
				required => 1,
				validate => 'String'
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Search Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'search_results',
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

my @data;
my @headers = ('Organism', 'Feature Name', 'Tag', 'Value');
my @feature_ids;
while(my $row = $results->next){
	push(@data,[$row->feature->organism->common_name, $row->feature->name, $row->cvterm->name, $row->value]);
}

if($options->{verbose}){
	print format_pretty(\@data);
}


my %results = (
	'Sheet1' => {
		header => \@headers,
		data => \@data,
	}
);
use CPT::OutputFiles;
my $data_out = CPT::OutputFiles->new(
	name   => 'results',
	GGO => $ggo,
);
$data_out->CRR(data => \%results);
