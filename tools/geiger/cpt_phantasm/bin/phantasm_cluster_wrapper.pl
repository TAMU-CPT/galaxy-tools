#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:   404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#

use strict;
use warnings;

# PODNAME: phantasm_cluster_wrapper
use CPT::GalaxyGetOpt;
my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
	],
	'outputs' => [
		[
			'cluster_dendrogram',
			'hclust dendrogram of clusters',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'dendrogram',
				data_format    => 'Dummy',
				default_format => 'Dummy',
			}
		],
		[
			'newick_tree',
			'Newick Tree from Hierarchical Clustering',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'tree.newick',
				data_format    => 'Dummy',
				default_format => 'Dummy',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PHAnTASM.ClusterWrapper',
		'appname' => 'Hierarchical Clustering',
		'appdesc' => 'runs R stat module\'s hclust algorithm',
		'appvers' => '1.94',
	],
);

use File::ShareDir;
use File::Spec;
my $dir     = File::ShareDir::dist_dir('CPT-PHAnTASM');
my $rscript = File::Spec->catfile( $dir, 'phantasm_cluster.R' );

use CPT::OutputFiles;
my $png_output = CPT::OutputFiles->new(
	name => 'cluster_dendrogram',
	GGO => $ggo,
);
my ($png_loc) = $png_output->subCRR(
	filename => 'dendrogram', extension => 'png', data => '', 
	data_format => 'Dummy', format_as => 'Dummy',
);
my $tree_output = CPT::OutputFiles->new(
	name => 'newick_tree',
	GGO => $ggo,
);
my ($tree_loc) = $tree_output->subCRR(
	filename => 'genome_cluster', extension => 'newick', data => '', 
	data_format => 'Dummy', format_as => 'Dummy',
);

use IPC::Run3;
my @cmd = ('Rscript', $rscript, $options->{file},$png_loc, $tree_loc);
my ($in, $out, $err);
run3 \@cmd, \$in, \$out, \$err;
if($err){
	print "Error Output: $err";
}

=head1 NAME

PHAnTASM Genomic Hierarchical Clustering

=head1 DESCRIPTION

This tool accepts the output of the PHAnTASM Comparison Map generator in order to build and plot a dendrogram from the results. This tool is a thin wrapper around the C<stat> package in R's hclust method.

=cut
