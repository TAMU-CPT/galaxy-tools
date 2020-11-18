use Bio::Chado::Schema;
#use Data::Format::Pretty::Console qw(format_pretty);
use strict;
use warnings;

# PODNAME: list_orgs.pl

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
	],
	'outputs' => [
		[
			'results',
			'Organism Listing',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'orgs',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'CHAPLIN_list_organisms',
		'appname' => 'List Organisms',
		'appdesc' => 'lists organisms in a chado database',
	]
);



my $dsn = "dbi:Pg:dbname=" . $options->{database} . ";host=cpt.tamu.edu;port=5432;sslmode=require";
my $user = "charm_admin";
my $password = "oNFkI0KyoGygRp8Zf7jOVIrR1VmsOWak";
my $chado = Bio::Chado::Schema->connect( $dsn, $user, $password );


my $results = $chado->resultset('Organism::Organism')->search();

my @data;

while(my $row = $results->next){
	push(@data, [$row->id, $row->genus, $row->species, $row->common_name]);
}

my %crr_data = (
	'Sheet1' => {
		header => ['ID', 'Genus', 'Species', 'Common Name'],
		data => \@data,
	}
);

use CPT::OutputFiles;
my $data_out = CPT::OutputFiles->new(
	name   => 'results',
	GGO => $ggo,
);
$data_out->CRR(data => \%crr_data);
