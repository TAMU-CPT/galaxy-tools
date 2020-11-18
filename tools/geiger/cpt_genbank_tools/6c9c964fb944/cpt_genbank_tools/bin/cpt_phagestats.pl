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

use CPT;
use Data::Dumper;

my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[ 'file|f', 'Input file',
			{
				validate => 'File/Input',
				#file_format => ['Genbank'],
				required => 1,
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Genbank Info and Statistics',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'info',
				data_format    => 'text/html',
				default_format => 'HTML',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PhageInfo',
		'appname' => 'Genbank Genome Information',
		'appdesc' => 'calculates various statistics and information about a genbank file',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

use CPT::Report::HTML;
my $report = CPT::Report::HTML->new();

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio_object = $bio->getSeqIO($options->{file});
# For all genomes in the GBK file
while(my $seq_object = $seqio_object->next_seq){
	my @features = $seq_object->get_SeqFeatures();
	my @genes;
	my $bases = 0;
	my @density_map;
	use Statistics::Descriptive;
	my $stat = Statistics::Descriptive::Full->new();
	my %g;
	foreach(@features){
		if($_->primary_tag eq 'CDS'){
			my %tags = map {$_ => 1} $_->get_all_tags();
			if(!defined($tags{'pseudo'}) && !defined($tags{'pseudogene'})){
				push(@genes, $_);
				my $len = $_->end() - $_->start();
				$bases += $len;
				$stat->add_data($len);
				for(my $i = int($_->start/1000); $i <= int($_->end/1000); $i++){
					$density_map[$i] += 1;
				}
				foreach my $b(split //,$_->seq->seq){
					$g{$b}++;
				}
			}
		}
	}
	my $density = 0;
	foreach(my $i=0;$i<scalar@density_map;$i++){
		if(defined $density_map[$i]){
			$density += $density_map[$i];
		}
	}
	$density /= scalar @density_map;
	my %go;
	foreach my $b(split //, $seq_object->seq){
		$go{$b}++;
	}

	my $g_sum = 0;
	foreach(keys(%g)){
		$g_sum+=$g{$_};
	}
	my $go_sum = 0;
	foreach(keys(%go)){
		$go_sum+=$go{$_};
	}








	$report->h1("Overview of " . $seq_object->display_id());
	$report->list_start('bullet');
	$report->list_element('Number of bases: ' . $seq_object->length());
	$report->list_element('Number of features: ' . scalar(@features));
	$report->list_end();

	$report->h3("Genes (CDS features withouth a /pseudo or /pseudogene qualifier)");
	$report->list_start('bullet');
	$report->list_element('count: ' . scalar(@genes));
	$report->list_element('bases: ' . $bases);
	$report->list_element(sprintf('density %0.3f genes per kb', $density));
	$report->list_element('average length: ' . int($stat->mean));
	$report->list_element('gene sequence composition:');
		$report->list_start('bullet');
		foreach(sort(keys(%g))){
			$report->list_element(sprintf("%s content: %s (%0.3f%%)", $_ , $g{$_} , 100*$g{$_}/$g_sum ));
		}
		$report->list_end();
	$report->list_element(sprintf('GC Percentage: %0.3f', 100*($g{'G'} + $g{'C'})/$g_sum ));
	$report->list_end();

	$report->h3('Overall sequence composition');
	$report->list_start('bullet');
	foreach(sort(keys(%g))){
		$report->list_element(sprintf("%s content: %s (%0.3f%%)", $_ , $go{$_} , 100*$go{$_}/$go_sum ));
	}
	$report->list_element(sprintf('GC Percentage: %0.3f', 100*($go{'G'} + $go{'C'})/$go_sum ));
	$report->list_end();

}

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'results',
        libCPT => $libCPT,
);
$output->CRR(data => $report->get_content());

=head1 DESCRIPTION

Genomic Info, very similar to the "Genome Overview" in artemis.

=cut

