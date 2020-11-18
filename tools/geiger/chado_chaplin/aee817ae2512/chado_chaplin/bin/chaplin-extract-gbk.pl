use Bio::Chado::Schema;
use Data::Format::Pretty::Console qw(format_pretty);
use strict;
use warnings;


# PODNAME: extract-gbk.pl


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
			'cn',
			'Organism Common Name',
			{
				required => 1,
				validate => 'String'
			}
		],
	],
	'outputs' => [
		[
			'gbk',
			'Extracted Genomic Data',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'extracted',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank'
			}
		],
	],
	'defaults' => [
		'appid'   => 'CHAPLIN_extract_gbk',
		'appname' => 'Extract Genbank',
		'appdesc' => 'converts chado database entry to genbank file',
	]
);



my $dsn = "dbi:Pg:dbname=" . $options->{database} . ";host=cpt.tamu.edu;port=5432;sslmode=require";
my $user = "charm_admin";
my $password = "oNFkI0KyoGygRp8Zf7jOVIrR1VmsOWak";
my $chado = Bio::Chado::Schema->connect( $dsn, $user, $password );

my $results = $chado->resultset('Organism::Organism')->search(
	{ common_name => { -like => $options->{cn} } },
);

my @organism_ids;
my @found_orgs;
while(my $row = $results->next){
	push(@organism_ids, $row->id);
	push(@found_orgs,$row->common_name);
}

if(scalar @organism_ids > 1 || scalar @organism_ids == 0){
	die 'Please refine your search query so as to only select one organism. Found: [' . join(",", @found_orgs) . ']';
}


my $features = $chado->resultset('Sequence::Feature')->search(
	{ organism_id => $organism_ids[0] },
);



use Bio::SeqFeature::Generic;

# Must pre-create then hold features
my @features;
# store info for the eventual seqobj
my ($seq,$id);
# We need to resolve cvterms, however we cannot do this until after we've fetched features. It's a bit of a pain.
my %cvterms_to_lookup;
while(my $feat = $features->next){
	if($feat->seqlen){
		$seq = $feat->residues;
		$id = $feat->uniquename;
	}
	$cvterms_to_lookup{$feat->type_id}++;
	my @featprops = $feat->featureprops;
	my @featlocs = $feat->featureloc_features;
	foreach my $loc(@featlocs){
		foreach(@featprops){
			$cvterms_to_lookup{$_->type_id}++;
		}
		push(@features, {
			primary_tag => $feat->type_id,
			start       => $loc->fmin+ ($loc->strand == 1 ? 1:0),
			end         => $loc->fmax+ ($loc->strand == 1 ? 1:0),
			strand      => $loc->strand,
			featprops   => [map { [ $_->type_id, $_->value, $_->rank]} @featprops]
		});
	}
}


# Lookup the cvterms
my @cvterm_list = map { { cvterm_id => $_ } } keys(%cvterms_to_lookup);
my $cvterms = $chado->resultset('Cv::Cvterm')->search(
	\@cvterm_list
);
my %cv;
while(my $term = $cvterms->next){
	$cv{$term->cvterm_id} = $term->name;
}

# Load the qualifier mappings
use File::Spec;
#my $dir = File::ShareDir::dist_dir('CPT-CHAPLIN');
my $dir = ".";
my $qual_map_loc = File::Spec->catfile($dir,'qualifier_mapping');
open(my $qual_map,'<',$qual_map_loc);
my %qual_trans;
while(<$qual_map>){
	if($_ !~ /^#/){
		chomp $_;
		my ($a,$b) = split(/\s+/,$_);
		if($b){
			$qual_trans{$a} = $b;
		}
	}
}

# Create the sequence
use Bio::Seq;
my $seq_obj = Bio::Seq->new(
	-seq        => $seq,
	-display_id => $id,
);

# Create and add features
foreach(@features){
	my %f = %{$_};
	my %keys;
	# Translate the feature props into qualifiers
	foreach my $ref(@{$f{featprops}}){
		my ($id, $val0) = @{$ref};
		my $key = $cv{$id};
		# If there's a mapping we ought to use, use that
		if($qual_trans{$key}){
			$key = $qual_trans{$key};
		}
		unless($keys{$key}){
			$keys{$key} = [];
		}
		push($keys{$key},$val0);
	}
	# Create our new feature
	my $new_feat = new Bio::SeqFeature::Generic(
		-start       => $f{start},
		-end         => $f{end},
		-strand      => $f{strand},
		-primary_tag => $cv{$f{primary_tag}},
		-tag         => \%keys,
	);
	$seq_obj->add_SeqFeature($new_feat)
}

use CPT::OutputFiles;
my $data_out = CPT::OutputFiles->new(
	name   => 'gbk',
	libCPT => $libCPT,
);
$data_out->CRR(data => $seq_obj);
