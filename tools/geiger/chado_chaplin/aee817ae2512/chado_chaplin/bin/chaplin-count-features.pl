use Bio::Chado::Schema;
use strict;
use warnings;

# PODNAME: count_features.pl

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
	],
	'outputs' => [
		[
			'results',
			'Feature Counts',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'feature_counts',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'CHAPLIN_count_features',
		'appname' => 'Count Features',
		'appdesc' => 'counts numbers of each feature type in the entire database',
	]
);



my $dsn = "dbi:Pg:dbname=" . $options->{database} . ";host=cpt.tamu.edu;port=5432;sslmode=require";
my $user = "charm_admin";
my $password = "oNFkI0KyoGygRp8Zf7jOVIrR1VmsOWak";
my $chado = Bio::Chado::Schema->connect( $dsn, $user, $password );


my $results = $chado->resultset('Sequence::Feature')->search();

my %cv_terms;
while(my $row = $results->next){
	$cv_terms{$row->type_id}++;
}

my %ppresults;
foreach my $k (keys(%cv_terms)){
	$results = $chado->resultset('Cv::Cvterm')->search({
		cvterm_id => $k
	});
	while(my $row = $results->next){
		$ppresults{$row->name} = $cv_terms{$k};
	}
}

if($options->{verbose}){
	use Data::Format::Pretty::Console qw(format_pretty);
	print format_pretty(\%ppresults);
}

my @data;
foreach(sort {$ppresults{$b} <=> $ppresults{$a}} keys %ppresults){
	push(@data,[$_,$ppresults{$_}]);
}
my %crr_data = (
	'Sheet1' => {
		header => ['tag','count'],
		data => \@data,
	}
);
use CPT::OutputFiles;
my $data_out = CPT::OutputFiles->new(
	name   => 'results',
	libCPT => $libCPT,
);
$data_out->CRR(data => \%crr_data);
