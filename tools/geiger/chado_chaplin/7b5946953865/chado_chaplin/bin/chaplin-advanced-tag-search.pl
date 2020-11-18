use Bio::Chado::Schema;
use strict;
use warnings;

# PODNAME: chaplin-find-orgs-without-tag.pl
use Data::Dumper;
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
			'Tag Value to be queried in the database. Information about the feature will be printed if it exists in that organism, or a notice will be printed if it was not found. * is a multiple character wildcard, and ? is a single character.',
			{
				required => 1,
				validate => 'String',
				multiple => 1,
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
				default        => 'orgs',
				data_format    => 'text/html',
				default_format => 'HTML'
			}
		],
	],
	'defaults' => [
		'appid'   => 'CHAPLIN_advanced_search',
		'appname' => 'Advanced Tag Search',
		'appdesc' => 'searches for organisms with and without a specified query string',
	]
);



my $dsn = "dbi:Pg:dbname=" . $options->{database} . ";host=cpt.tamu.edu;port=5432;sslmode=require";
my $user = "charm_admin";
my $password = "oNFkI0KyoGygRp8Zf7jOVIrR1VmsOWak";
my $chado = Bio::Chado::Schema->connect( $dsn, $user, $password );

my %color_scheme = (
	"CDS" => {
		"color" =>  "#1c86ee",
		"border" =>  1,
		"plot" =>  1
	}
);

my @queries_orig = @{$options->{query}};
my @queries;
foreach(@queries_orig){
	my $query = "$_";
	# Escape all %s
	$query =~ s/%/%%/;
	# Add in our custom %s
	$query =~ s/\*/%/g;

	# Escape all _s
	$query =~ s/_/__/;
	# Add in our custom _s
	$query =~ s/\?/_/g;
	push(@queries, $query);
}


my $results = $chado->resultset('Organism::Organism')->search(undef, { 'order_by' => {'-asc', 'species', '-asc', 'genus', '-asc', 'common_name' } });

use CPT::Report;
use CPT::Report::HTML;
my $report = CPT::Report::HTML->new();
$report->title("Query");



my %cvterms_to_lookup;
my @haves;
my @havenots;

# Search through orgs
while(my $row = $results->next){
	#$report->h1(sprintf("<i>%s. %s</i> (%s)", substr($row->genus,0,1), $row->species, $row->common_name));
	my %org_results = (
		features => query_organism($row->id),
		genus =>  $row->genus,
		species =>  $row->species,
		common_name =>  $row->common_name,
	);

	if(scalar @{$org_results{features}} > 0){
		push(@haves, \%org_results);
	}else{
		push(@havenots, \%org_results);
	}
}

my %tn = (
	primary_tag => "Primary Tag",
	start       => "Start",
	end         => "End",
	strand      => "Strand",
	timelastmodified => "Last Modified",
);

my %cv = load_cvterms();
$report->h1("Organisms with features matching " . join(", ", @queries_orig));
foreach my $organism(@haves){
	my %d = %{$organism};
	$report->h3(sprintf("%s %s", $d{genus},$d{species}));
	# Map of features
	my @fp = @{$d{features}};
	$report->a(generateMapOfFeatures(@fp));

	$report->a('<div class="row">');
	foreach(@fp){
		my %fi = %{$_};


		$report->a('<div class="col-md-6">');
		#$report->finalize_table();
		$report->a("<pre>" . genbankFeatureString(%fi) . "</pre>");
		$report->a('</div>');
	}
	$report->a('</div>');
}
$report->h1("Organisms without matching features");
$report->list_start('bullet');
foreach my $organism(@havenots){
	my %d = %{$organism};
	$report->list_element(sprintf("%s %s", $d{genus}, $d{species}));
}
$report->list_end();

sub genbankFeatureString{
	my (%fi) = @_;
	my $resp;

	if($fi{strand} == 1){
		$resp .= sprintf("     CDS             %s..%s\n", $fi{start},
			$fi{end})
	}else{
		$resp .= sprintf("     CDS             complement(%s..%s)\n",
			$fi{start}, $fi{end})
	}

	foreach(@{$fi{featprops}}){
		my ($id,$val,$rank) = @{$_};
		$resp .= sprintf("                     /%s=\"%s\"\n", $cv{$id}, $val);
	}
	return $resp;
}



use CPT::Plot::Base;
use Bio::SeqFeature::Generic;
sub generateMapOfFeatures {
	my (@features) = @_;
	my $min = 1_000_000_000;
	my $max = 0;
	foreach(@features){
		my %f = %{$_};
		if($f{start} < $min){
			$min = $f{start}
		}
		if($f{start} > $max){
			$max = $f{start}
		}
		if($f{end} < $min){
			$min = $f{end}
		}
		if($f{end} > $max){
			$max = $f{end}
		}
	}
	$min -= .1 * ($max - $min);
	$max += .1 * ($max - $min);
	$min = int($min);
	$max = int($max);

	my @real_features;
	foreach(@features){
		my $feat = Bio::SeqFeature::Generic->new(
			-start  => ${$_}{start}  , -end         => ${$_}{end} ,
			-strand => ${$_}{strand} , -primary_tag => 'CDS'
		);
		push(@real_features, $feat);
	}

	if(!defined $min || !defined $max){
		return 'Could not build SVG';
	}

	my $svg_control = CPT::Plot::Base->new(
		'_ft_count' => 0,
		'_internal_maxrowlength' => 0,
		'double_line_for_overlap' => 1,
		'genome_length' => 4000,
		'ils' => 000,
		'justified' => 'justify',
		'label' => undef,
		'label_callouts' => undef,
		'label_from' => 'numeric',
		'label_numbering_count' => 1,
		'label_numeric_features' => undef,
		'label_pos' => 'above',
		'label_query' => undef,
		'label_shrink_mode' => 'cutoff',
		'label_text_source' => 'locus_tag',
		'line_count' => 1,
		'opacity' => '0.7',
		'rows' => 1,
		'separate_strands' => 1,
		'split_factor' => '1.02',
		'view' => 'alt_artemis',
		'x_offset' => 30,
		'y_offset' => 70,
		'zoom' => 50,
		'color_scheme' => \%color_scheme,
		'features' => \@real_features,
		'start' => $min,
		'end' => $max,
	);
	$svg_control->init();
	$svg_control->partitionLines();
	$svg_control->createSVG();
	return $svg_control->getSVG->xmlify();
}


# Search within an org.
sub query_organism {
	my ($oid) = @_;
	my $features = $chado->resultset('Sequence::Featureprop')->search(
		{
			'feature.organism_id' => $oid,
			value => { -like => \@queries },
		},
		{ join => ['cvterm', 'feature'] , },
	);
	my $i = 0;
	my @hits;
	while(my $row = $features->next()){
		if($i == 0 && $options->{verbose}){
			print STDERR "Hit in $oid\n";
		}
		$i++;
		my $feature_ref = feature_info($row->feature_id());
		if(defined $feature_ref){
			push(@hits, $feature_ref);
		}
	}
	return \@hits;
}
sub feature_info {
	# Returns 1 feature.
	my ($fid) = @_;
	my $feature_info = $chado->resultset('Sequence::Feature')->search(
		{
			'feature_id' => $fid,
		},
		#{ join => ['featureloc_features'] , },
	)->first;

	my @featprops = $feature_info->featureprops;
	my $loc = $feature_info->featureloc_features->first;
	# Lazily evaluate
	foreach(@featprops){
		$cvterms_to_lookup{$_->type_id}++;
	}
	$cvterms_to_lookup{$feature_info->type_id}++;
	if(defined $loc){
		my %feature_info = (
			primary_tag => $feature_info->type_id,
			start       => $loc->fmin+ ($loc->strand == 1 ? 1:0),
			end         => $loc->fmax+ ($loc->strand == 1 ? 1:0),
			strand      => $loc->strand,
			featprops   => [map { [ $_->type_id, $_->value, $_->rank]} @featprops],
			timelastmodified => $feature_info->timelastmodified,
		);
		return \%feature_info;
	}
	return;
}
sub load_cvterms {
	# Lookup the cvterms
	my @cvterm_list = map { { cvterm_id => $_ } } keys(%cvterms_to_lookup);
	my $cvterms = $chado->resultset('Cv::Cvterm')->search(
		\@cvterm_list
	);
	my %cv;
	while(my $term = $cvterms->next){
		$cv{$term->cvterm_id} = $term->name;
	}
	return %cv;
}


use CPT::OutputFiles;
my $data_out = CPT::OutputFiles->new(
	name   => 'results',
	GGO => $ggo,
);
$data_out->CRR(data => $report->get_content);
