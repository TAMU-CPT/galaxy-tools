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

use CPT::GalaxyGetOpt;
use Data::Dumper;
use CPT::Bio;
use CPT::Circos::Conf;
my $bio = CPT::Bio->new();

my @argv_copy;
foreach(@ARGV){push(@argv_copy, "$_");}

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file|f', 'Input file', { validate => 'File/Input',
				file_format => ['genbank', 'embl', 'txt'],
			} ],
		[],
		['Track Configuration'],
		['track_key', 'Key to select from genbank data', { validate => 'Genomic/Tag', required => 1, multiple => 1 } ],
		['track_feature_filter_invert', 'Should the qualifier search be inverted?', { validate => 'Option', options => { 'yes', 'Yes', 'no', 'No' } , multiple => 1 } ],
		['track_feature_filter_hastag', 'Select a tag which should be present in that qualifier (e.g., signal/tmhelix/pseudo)', { validate => 'String' , multiple => 1 } ],
		['track_feature_filter_textquery', 'Specify text which MUST be in that tag', { validate => 'String' , multiple => 1 } ],
		['track_feature_filter_strand', 'Which strand should the feature appear on?', { validate => 'Option', options => { 'f', 'Forward', 'r', 'Reverse', 'a', 'Any' } , multiple => 1 } ],
		[],
		['enable_gc_skew_plot', 'Enable/Disable calculation of GC Skew Plot', { validate => 'Flag' } ],
		['gc_skew_plot_window_size', 'Window size for calculation of GC Skew', { validate => 'Int', min => 1000, default => 10000} ],
		['gc_skew_plot_step_size', 'Step size for calculation of GC Skew', { validate => 'Int', min => 200, default => 200 } ],
	],
	'outputs' => [
		[
			'output_circos_confs',
			'Circos Configuration Files',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'out',
				data_format    => 'archive',
				default_format => 'zip',
			}
		],
	],
	'defaults' => [
		'appid'   => 'CircularDNAPlotter',
		'appname' => 'Circos based DNAPlotter',
		'appdesc' => 'plots genomes similar to Artemis\'s DNAPlotter',
		'appvers' => '1.94.1',
	],
);

#perl cpt_dnaplotter.pl \
	#-f ../t/test-files/moon.gbk \
	#--track_key CDS --track_feature_filter_invert yes --track_feature_filter_hastag pseudo --track_feature_filter_strand f \
	#--track_key CDS --track_feature_filter_invert yes --track_feature_filter_hastag pseudo  --track_feature_filter_strand r \
	#--track_key CDS --track_feature_filter_hastag pseudo  --track_feature_filter_strand a \
	#--track_key tRNA  --track_feature_filter_strand a \
	#--track_key CDS  --track_feature_filter_hastag signal --track_feature_filter_strand a \
	#--track_key CDS  --track_feature_filter_hastag tmhelix --track_feature_filter_strand a


my @reorg_args = ();

my $cum_gc_ske_mean = 0;



my %latest = ();
for(my $i = 0; $i < scalar(@argv_copy); $i++){
	my $c = $argv_copy[$i];
	# We have entered a new one block
	if($c eq '--track_key'){
		# If we have loaded data
		if(scalar(keys(%latest)) > 0){
			my %copy;
			foreach(keys(%latest)){
				$copy{$_} = "" . $latest{$_};
			}
			push(@reorg_args, \%copy);
		}

		# Clean out latest to prep for new data
		foreach(keys(%latest)){
			delete $latest{$_};
		}
	}

	if($c =~ /^--track_(.*)/){
		$latest{$1} = $argv_copy[$i+1];
		# Artificially bump so we don't bother looking at the answer to
		# this question. We can "get away" with this because none of
		# the options are flags. However, I've disabled it in the event
		# that flags ARE introduced
		#$i++;
	}
}
push(@reorg_args, \%latest);
#$VAR1 = [
          #{
            #'feature_filter_invert' => 'yes',
            #'feature_plot_color' => '005500',
            #'feature_filter_strand' => 'f',
            #'feature_filter_hastag' => 'pseudo',
            #'key' => 'CDS'
          #},
          #{
            #'feature_filter_strand' => 'a',
            #'key' => 'RBS'
          #}
        #];
my @files = ();

my $number_of_tracks = 0;
sub register_track {
	#my $r0 = ( 90 - (10 * $number_of_tracks - 1)/2) / 100;
	#my $r1 = ( 90 - (10 * $number_of_tracks - 9)/2) / 100;
	my $r0 = ( 90 - (10 * $number_of_tracks - 1)/1) / 100;
	my $r1 = ( 90 - (10 * $number_of_tracks - 9)/1) / 100;
	$number_of_tracks++;
	return ($r0, $r1);
}

sub circosconf {
	my $cc = CPT::Circos::Conf->new();
	$cc->include('etc/colors_fonts_patterns.conf');
	# Features to plot along the genome
	$cc->include('ideogram.conf');
	# markings indicating position along genome
	$cc->include('ticks.conf');
	# Genome data
	$cc->set('karyotype', 'karyotype.conf');
	# Default image params are fine
	$cc->start_block('image');
	$cc->include('etc/image.conf');
	$cc->end_block();
	# ???
	$cc->set('chromosome_units', '1000');
	$cc->set('chromosome_display_default', 'yes');
	#$cc->include('highlights.conf');
	$cc->include('plots.conf');

	$cc->include('etc/housekeeping.conf');
	my $result = $cc->finalize();
	$cc = CPT::Circos::Conf->new();
	return $result;
}
sub ideogramconf{
	my $cc = CPT::Circos::Conf->new();
	$cc->start_block('ideogram');
	$cc->start_block('spacing');
	$cc->set('default','0u');
	$cc->set('break','0u');
	$cc->end_block();

	$cc->set('thickness', '20p');
	$cc->set('stroke_thickness', '2');
	$cc->set('stroke_color', 'black');
	$cc->set('fill','no');
	$cc->set('fill_color','black');
	$cc->set('radius','0.85r');
	$cc->set('show_label','yes');
	$cc->set('label_font','default');
	$cc->set('label_radius','dims(ideogram,radius) + 0.05');
	$cc->set('label_size','36');
	$cc->set('label_parallel','yes');
	$cc->set('label_case','upper');

	$cc->set('band_stroke_thickness','2');
	$cc->set('show_bands','yes');
	$cc->set('fill_bands','yes');
	$cc->end_block();

	return $cc->finalize();
}
sub generate_feature_table {
	my ($filename, %filter) = @_;
	print "Filtering on features\n";
	print Dumper \%filter;
	my $seqio_object = Bio::SeqIO->new(-file => $options->{file}, -format=>'genbank');
	# Only handing first sequence.
	my $seq_object = $seqio_object->next_seq;
	my $parent = $seq_object->display_id();
	# Feature data
	my @features;
	foreach my $feat($seq_object->get_SeqFeatures()){
		if($feat->primary_tag() eq $filter{key}){
			# If they've said "hastag" AND we do indeed have that tag AND we haven't inverted this filter.
			if(
				($filter{feature_filter_hastag} && $feat->has_tag($filter{feature_filter_hastag}) && !$filter{feature_filter_invert})
				||
				($filter{feature_filter_hastag} && $filter{feature_filter_invert} && !$feat->has_tag($filter{feature_filter_hastag}))
				||
				(! $filter{feature_filter_hastag})
			){
				if(
					!$filter{feature_filter_strand}
					||
					($feat->strand() == 1 && ( $filter{feature_filter_strand} eq 'f' ||  $filter{feature_filter_strand} eq 'a' ))
					||
					($feat->strand() == -1 && ( $filter{feature_filter_strand} eq 'r' ||  $filter{feature_filter_strand} eq 'a' ))
					||
					($feat->strand() == 0 && ( $filter{feature_filter_strand} eq 'a' ))
				){
					push(@features, join(' ', $parent, $feat->start, $feat->end));
				}
			}
		}
	}
	print "Found " . scalar @features . " features \n";
	push(@files, [ 'data/' . $filename, join("\n", @features) ] );
}
sub plotsconf{
	my $cc = CPT::Circos::Conf->new();

	$cc->start_block('plots');


	my $idx = 0;
	foreach(@reorg_args){
		my %filter = %{$_};
          #{
            #'feature_filter_invert' => 'yes',
            #'feature_plot_color' => '005500',
            #'feature_filter_strand' => 'f',
            #'feature_filter_hastag' => 'pseudo',
            #'key' => 'CDS'
          #},
		$idx++;
		my $filename = sprintf('org.features.%s.txt', $idx);
		generate_feature_table($filename, %filter);

		my ($r0,$r1) = register_track();
		$cc->start_block('plot');
		$cc->set('type','tile');
		$cc->set('file',$filename);
		$cc->set('orientation', 'center');
		$cc->set('thickness', '30');
		$cc->set('r1', $r1 . 'r');# '0.78r');
		$cc->set('r0', $r0 . 'r');# '0.72r');
		$cc->set('layers','3');
		$cc->set('layers_overflow','collapse');
		$cc->set('color','paired-6-qual-' . $idx);
		$cc->end_block();
	}

	if($options->{enable_gc_skew_plot}){
		my ($r0,$r1) = register_track();
		$cc->start_block('plot');
		$cc->set('type','histogram');
		$cc->set('file','gc.txt');
		$cc->set('r1',$r1 . 'r');
		$cc->set('r0',$r0 . 'r');
		$cc->set('fill_color','purple');
		$cc->set('orientation','in');
		$cc->start_block('rules');
		$cc->start_block('rule');
		$cc->set('condition','var(value) < 0');
		$cc->set('fill_color', 'green');
		$cc->end_block();
		$cc->end_block();
		$cc->end_block();
	}

	#$cc->start_block('plot');
	#$cc->set('type','histogram');
	#$cc->set('file','gc_cumulative.txt');
	#$cc->set('r1','0.6r');
	#$cc->set('r0','0.55r');
	#$cc->set('fill_color','purple');
	#$cc->set('orientation','out');
	#$cc->start_block('rules');
	#$cc->start_block('rule');
	#$cc->set('condition','var(value) < ' . $cum_gc_ske_mean);
	#$cc->set('fill_color', 'green');
	#$cc->end_block();
	#$cc->end_block();
	#$cc->end_block();

	$cc->end_block();
	return $cc->finalize();
}
sub ticksconf{
	my $cc = CPT::Circos::Conf->new();

	$cc->set('show_ticks','yes');
	$cc->set('show_tick_labels','yes');
	$cc->start_block('ticks');
	$cc->set('radius','1r');
	$cc->set('color','black');
	$cc->set('thickness','2p');
	$cc->set('multiplier','1e-3');
	$cc->set('format','%d');

	$cc->start_block('tick');
	$cc->set('spacing','1000u');
	$cc->set('size','10p');
	$cc->end_block();

	$cc->start_block('tick');
	$cc->set('spacing','10000u');
	$cc->set('size','15p');
	$cc->set('show_label','yes');
	$cc->set('label_size','20p');
	$cc->set('label_offset','10p');
	$cc->set('format','%d');
	$cc->end_block();
	$cc->end_block();
	return $cc->finalize();
}
sub karyotype {
	my @karyotype_data = ();
	my $seqio_object = Bio::SeqIO->new(-file => $options->{file}, -format=>'genbank');
	# Only handing first sequence.
	my $seq_object = $seqio_object->next_seq;
	# Main 'chromosome' data
	push(@karyotype_data, join(' ', 'chr', '-',$seq_object->display_id(),$seq_object->display_id(),0, $seq_object->length(),'black'));

	return join("\n", @karyotype_data);
}

sub gcgraph_cumulative {
	my @gcdata = ();
	my $seqio_object = Bio::SeqIO->new(-file => $options->{file}, -format=>'genbank');
	# Only handing first sequence.
	my $seq_object = $seqio_object->next_seq;

	my $parent = $seq_object->display_id();

	my $seq = $seq_object->seq();
	my $sep = int($options->{gc_skew_plot_window_size}/2);
	my $stepsep = int($options->{gc_skew_plot_step_size}/2);
	my $cumulative_gc_skew = 0;
	my @cumgc_vals;

	my $count = 0;
	foreach(my $i = $sep; $i < $seq_object->length() - $sep - $options->{gc_skew_plot_step_size}; $i += $options->{gc_skew_plot_step_size}){
		$count++;
		$cumulative_gc_skew += _calculate_gc_skew_for_seq(substr($seq,$i-$sep,$options->{gc_skew_plot_window_size})),
		push(@cumgc_vals, $cumulative_gc_skew);
		push(@gcdata, join(" ",
			$parent,
			$i - $stepsep,
			$i + $stepsep,
			$cumulative_gc_skew
		));
	}

	my $sum = 0;
	foreach(@cumgc_vals){$sum += $_;}
	$cum_gc_ske_mean = $sum / $count;

	return join("\n", @gcdata);
	# Main 'chromosome' data
}
sub gcgraph {
	my @gcdata = ();
	my $seqio_object = Bio::SeqIO->new(-file => $options->{file}, -format=>'genbank');
	# Only handing first sequence.
	my $seq_object = $seqio_object->next_seq;

	my $parent = $seq_object->display_id();

	my $seq = $seq_object->seq();
	my $sep = int($options->{gc_skew_plot_window_size}/2);
	my $stepsep = int($options->{gc_skew_plot_step_size}/2);
	foreach(my $i = $sep; $i < $seq_object->length() - $sep - $options->{gc_skew_plot_step_size}; $i += $options->{gc_skew_plot_step_size}){
		push(@gcdata, join(" ",
			$parent,
			$i - $stepsep,
			$i + $stepsep,
			_calculate_gc_skew_for_seq(substr($seq,$i-$sep,$options->{gc_skew_plot_window_size})),
		));
	}
	return join("\n", @gcdata);
	# Main 'chromosome' data
}
sub _calculate_gc_skew_for_seq {
	my ($seq) = @_;
	$seq = uc($seq);
	my %counts;
	foreach(split //,$seq){
		$counts{$_}++;
	}
	if($counts{G} + $counts{C} > 0){
		return ($counts{G} - $counts{C}) / ($counts{G} + $counts{C});
	}
	return 0;
}

push(@files, [ 'etc/karyotype.conf', karyotype() ]);
push(@files, [ 'etc/circos.conf', circosconf() ]);
push(@files, [ 'etc/ideogram.conf', ideogramconf() ]);
if($options->{enable_gc_skew_plot}){
	push(@files, [ 'data/gc.txt', gcgraph() ]);
	push(@files, [ 'data/gc_cumulative.txt', gcgraph_cumulative() ]);
}

push(@files, [ 'etc/plots.conf', plotsconf() ]);
push(@files, [ 'etc/ticks.conf', ticksconf() ]);


use Archive::Any::Create;
my $archive = Archive::Any::Create->new();
$archive->container('conf');
foreach(@files){
	my ($location, $content) = @{$_};
	$archive->add_file($location, $content);
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'output_circos_confs',
	GGO => $ggo,
);
$crr_output->CRR(data => $archive);

=head1 NAME

DNAPlotter

=head1 DESCRIPTION

Much like artemis's DNAPlotter, this tool plots genomes in a ciruclar fashion, and can plot gc-deviation tracks as well. The options are somewhat reduced compared to artemis, so if you need something that isn't available in this version please file a bug report.

=head1 USE

Each track has several parameters:

=over 4

=item track_key

This selects a set of features from a genbank file, e.g., CDSs or tRNAs

=item track_feature_filter_invert

This parameter will invert (negate) the results of whatever query parameters you specify after it.

=item track_feature_filter_hastag

Require that a feature has a specific tag....

=item track_feature_filter_textquery

...with this specific text in it

=item track_feature_filter_strand

Which strand should the feature be on (not inverted)

=back

Additionally, users are able to enable/disable GC skew plots. it's suggested that these are generally left alone, as they can quickly increase runtime.

=cut
