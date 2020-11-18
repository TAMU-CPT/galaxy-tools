#!/usr/bin/env perl
use strict;
use warnings;
use Storable;
use CPT;
use CPT::Bio::NW_MSA;
use Data::Dumper;
use CPT::Circos::Conf;
use POSIX;


my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[ 'file', 'PSM2 Data File', { validate => 'File/Input', required => 1 } ],
		[ 'user_ordering', 'List of genome IDs used in the analysis, can be comma/space/newline separated.', { validate => 'String', required => 1, multiple => 1,}],
		[],
		['Plot Options'],
		['percent_filled'   , 'Percentage of a whole block that an individual gene is'                  , { validate => 'Float', default=>'0.8', min => '0.1', max => '1.0' }],
		['ig_dist'          , 'Maximum length of links between genome comparisons'                      , { validate => 'Int', default => 100 }],
		['stroke_thickness' , 'Thickness of inter-genome links'                                         , { validate => 'Int', default => '2', min => 1, max => 10 } ],
		['every_nth'        , 'Plot every Nth gene a modified version of the main color for that genome', { validate => 'Int', default => '20'}],
		[],
		['Cutoffs'],
		['evalue' , 'Evalue cutoff' , { validate => 'Float' , default => 1e-4 } ] ,
		['dice'   , 'Dice cutoff'   , { validate => 'Float' , default => 50 } ]   ,
		[],
		['Alignment Options'],
		['mismatch'    , 'Mismatch Score' , { validate => 'Float' , default => -0.1 } ] ,
		['gap_penalty' , 'Gap Penalty'    , { validate => 'Float' , default => '0.0' } ] ,
		['match'       , 'Match Score'    , { validate => 'Float' , default => 5 } ]  ,
	],
	'outputs' => [
		[
			'output_circos_confs',
			'Output Circos Conf Object',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'psm3',
				data_format    => 'archive',
				default_format => 'tar.gz',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PSM.Plot',
		'appname' => 'PSM Plotter',
		'appdesc' => 'plots data from PSM Prep',
		'appvers' => '1.94.2',
	],
	'tests' => [
	],
);

my $percent_filled = $options->{percent_filled};
my $width          = 1000*$percent_filled;
my $spacing        = 1000-$width;

my $circos_base_location = '/opt/circos-0.66';



my $offset = ($width+$spacing)/2;
my $full_increment   = $width+$spacing;
my $halfwidth        = $width/2;

#my %option_map = (
	#'offset'           => ($width+$spacing)/2,
	#'full_increment'   => $width+$spacing,
	#'halfwidth'        => $width/2,
	#'gap_penalty'      => $options->{gap_penalty},
	#'match'            => $options->{match},
	#'heatmap'          => 1,
	#'heatmap_low'      => hex("0xCCCCCC"),
	#'heatmap_high'     => hex("0x333333"),
	#'heatmap_bucket'   => 8,
	#'every_nth'        => $options->{every_nth},
	#'user_ordering'    => $options->{user_ordering},
	#'dice'             => $options->{dice},
	#'evalue'           => $options->{evalue},
#);
#Color/name correspondance, to be used in writing the circos-0.63-pre1 files

my @user_ordering;
foreach(@{$options->{user_ordering}}){
	push(@user_ordering, split(/[,\s]+/, $_));
}
my @aligned_results;
my %fh_relationship;
my %precomputed_colour_hash;
my %protein_position_quicklookup;
my %response = ();
my %data_file = %{retrieve($options->{file})};

my %uo_idx;
for(my $i=0;$i<scalar @user_ordering;$i++){
	$uo_idx{$user_ordering[$i]} = $i;
}

align();
my %compmap_proteins = %{compmap_proteins()};
# a => "link text"
my %linkages = %{linkages()};
# a => b => "link text"

sub align{
	print STDERR "Aliging genomes\n";
	my $msa = CPT::Bio::NW_MSA->new(
		gap_penalty => $options->{'gap_penalty'},
		match_score => $options->{'match'},
		mismatch_score => $options->{'mismatch'},
		bidi => 1,
	);

	my @hits = @{$data_file{hit_table}};

	foreach my $hit(@hits){
		my ($from, $to, $evalue, $dice) = @{$hit};
		if($evalue < $options->{evalue} && $dice > $options->{dice}){
			$msa->add_relationship($from, $to);
		}
	}

	foreach my $genome(@user_ordering){
		my $gi_list_ref = $data_file{gbk}{$genome}{gi};#"GI" list
		$msa->align_list($gi_list_ref);
	}
	@aligned_results = $msa->merged_array();
}
sub compmap_proteins{
	my @Narr = ();#Keep count of how many items we've had in a single
	#column, so the modulus whe we're colouring them in will work properly,
	#rather than being on every Nth radial colum in the plot
	#
	my $_max = scalar @aligned_results;
	my %protein_files;

	for(my $i = 0; $i < scalar @aligned_results; $i++){
		# Get the current row from the PSM result object
		my @current_row = @{$aligned_results[$i]};
		for(my $j = 0; $j < scalar @current_row; $j++){
			if($current_row[$j] ne "-"){
				$protein_position_quicklookup{$current_row[$j]} = [($_max - $i -1),$j];
				$Narr[$j]++;
				my $color_str = '';
				if($Narr[$j] % $options->{every_nth} == 0){
					$color_str = 'fill_color=accent-8-qual-inv-' . ($j+1) . ',';
				}
				my $str = join(' ',
						'compmap ',
						(($_max - $i) * $full_increment - $halfwidth ),
						(($_max - $i) * $full_increment + $halfwidth ),
						"${color_str}f=". $current_row[$j]
					);
				$protein_files{$user_ordering[$j]} .= $str . "\n";
			}
		}
	}
	return \%protein_files;
}
sub linkages{
	my @hits = @{$data_file{hit_table}};
	my %links;
	foreach my $hit (@hits){
		my ($from, $to, $evalue, $dice) = @{$hit};
		if($evalue < $options->{evalue} && $dice > $options->{dice}){
			if(defined $protein_position_quicklookup{$from} && defined $protein_position_quicklookup{$to}){
				my ($theta0,$radius0) = @{$protein_position_quicklookup{$from}};
				my ($theta1,$radius1) = @{$protein_position_quicklookup{$to}};
				# If this is a self-self link, disable plotting because we don't care.
				# If ig_dist is disabled or distance is between them is less than our minimum
				if($radius0 != $radius1
					&& ($options->{'ig_dist'} == "-1" || abs($theta0-$theta1) <= $options->{'ig_dist'})
				){
					# Create the dataset
					my @row_data = ('compmap',
					);
					# We work under the assumption that all hits
					# are bi-directional, so we swap them to be
					# smallest first no matter what.
					if($radius1 < $radius0){
						my $tmp = $radius1;
						$radius1 = $radius0;
						$radius0 = $tmp;
						# We also want to add in reverse order
						push(@row_data,
							(($theta0+1)*$full_increment),
							(($theta1+1)*$full_increment),
						);
					}else{
						push(@row_data,
							(($theta1+1)*$full_increment),
							(($theta0+1)*$full_increment),
						);
					}

					# If it's a link with the same theta
					# value, then we'll go ahead and hide
					# behind the track to make it a little
					# prettier.
					my $zstr;
					if($theta0 == $theta1){
						$zstr="z=0";
					}else{
						$zstr="z=100";
					}
					
					# Create the additional row data
					push(@row_data,
						join(',', "dice=$dice", "color=" . colorstr($dice))
					);
					$links{$radius0}{$radius1} .= join(' ', @row_data) . "\n";
				}
			}
		}
	}
	return \%links;
}
sub colorstr {
	my ($dice) = @_;
	if($dice > 90) {
		return 'black';
	}else{
		return 'greys-9-seq-' . floor($dice / 10);
	}
}
sub circosconf {
	my $cc = CPT::Circos::Conf->new();
	$cc->include('etc/colors_fonts_patterns.conf');
	$cc->start_block('colors');
	$cc->set('accent-8-qual-inv-1', '42, 135, 42');
	$cc->set('accent-8-qual-inv-2', '111, 83, 150');
	$cc->set('accent-8-qual-inv-3', '182, 112, 46');
	$cc->set('accent-8-qual-inv-4', '178, 178, 53');
	$cc->set('accent-8-qual-inv-5', '13, 63, 128');
	$cc->set('accent-8-qual-inv-6', '183, 0, 96');
	$cc->set('accent-8-qual-inv-7', '116, 47, 0');
	$cc->set('accent-8-qual-inv-8', '45, 45, 45');
	$cc->end_block();
	# markings indicating position along genome
	$cc->include('example/etc/ideogram.conf');
	#$cc->include('rules.conf');
	# Genome data
	$cc->set('karyotype', 'karyotype.conf');
	# Default image params are fine
	$cc->start_block('image');
	$cc->include('etc/image.conf');
	$cc->end_block();
	#$cc->include('highlights.conf');
	$cc->include('plots.conf');
	#$cc->include('rules.conf');

	$cc->include('etc/housekeeping.conf');
	my $result = $cc->finalize();
	$cc = CPT::Circos::Conf->new();
	return $result;
}
sub karyotype {
	my @karyotype_data = (
		"chr - compmap compmap 0 ".((scalar @aligned_results)*1000+500)." white"
	);
	return join("\n", @karyotype_data);
}

my @files = ();

my $number_of_tracks = 0;
sub register_track {
	my ($r0,$r1) = calculate_individual_track($number_of_tracks);
	$number_of_tracks++;
	return ($r0, $r1);
}
sub calculate_individual_track {
	my ($idx) = @_;
	my $r0 = ( 90 - (10 * $idx - 1)/1) / 100;
	my $r1 = ( 90 - (10 * $idx - 9)/1) / 100;
	return ($r0, $r1);
}
sub genome_data {
	my $cc = CPT::Circos::Conf->new();
	# Map string back to position in array.
	#$cc->set('z',10);
	# Loop across all our protein data sets
	$cc->start_block('plots');
	foreach my $genome(@user_ordering){
		# Add protein file
		my $filename = sprintf('org.features.%s.txt', $genome);
		push(@files, [ 'data/'.$filename, $compmap_proteins{$genome}]);
		# Create associated tracks

		my ($r0,$r1) = register_track();
		$cc->start_block('plot');
		$cc->set('type','highlight');
		$cc->set('file', $filename);
		$cc->set('r0', $r0 .'r');
		$cc->set('r1', $r1 .'r');
		$cc->set('z', '50');
		$cc->set('fill_color','accent-8-qual-' . ($uo_idx{$genome} + 1));
		$cc->set('stroke_thickness', '1');
		$cc->set('stroke_color', 'black');
		$cc->end_block();
	}


	foreach my $from(@user_ordering){
		foreach my $to(@user_ordering){
			next if($from eq $to || $uo_idx{$from} > $uo_idx{$to});
		
			if($linkages{$uo_idx{$from}}{$uo_idx{$to}}){
				my $filename = sprintf('links.%s.%s.txt', $from, $to);
				push(@files, [ 'data/'.$filename, $linkages{$uo_idx{$from}}{$uo_idx{$to}}]);
				#push(@files, [ 'data/'.$filename, 'blaaaaaaah']);

				my ($r0a, $r0b) = calculate_individual_track($uo_idx{$to});
				my ($r1a, $r1b) = calculate_individual_track($uo_idx{$from});
		
				# If they're in this ordering, they will be pointing at
				# the "outsides" of each genome/protein track, so we
				# swap with the internal ones.
				if($r1b > $r0a){
					$r0a = $r1a;
					$r1b = $r0b;
				}

				$cc->start_block('plot');
				$cc->set('type','connector');
				$cc->set('thickness', $options->{stroke_thickness});
				$cc->set('file', $filename);
				if($r1b<$r0a){
					$cc->set('r0', $r1b .'r');
					$cc->set('r1', $r0a .'r');
				}else{
					$cc->set('r0', $r0a .'r');
					$cc->set('r1', $r1b .'r');
				}
				$cc->set('connector_dims', '0,0.3,0.4,0.3,0');
				$cc->set('color','black');
				$cc->end_block();
			}
		}
	}


	$cc->end_block();
	return $cc->finalize();
}
sub rulesconf {
	my ($self) = @_;
	my $cc = CPT::Circos::Conf->new();
	$cc->start_block('rules');
	for(my $i = 0; $i < 10; $i++){
		$cc->start_block('rule');
		$cc->set('importance', 10 - $i);
		$cc->set('condition', 'var(dice) > ' . (10*$i));
		if($i == 9){
			$cc->set('color', 'black');
		}else{
			$cc->set('color', 'gray' . (10 * ($i+1)));
		}
		#$cc->set('z',10-$i);
		$cc->end_block();
	}
	$cc->end_block();
}

push(@files, [ 'etc/karyotype.conf', karyotype() ]);
push(@files, [ 'etc/circos.conf', circosconf() ]);
#push(@files, [ 'etc/rules.conf', rulesconf() ]);
push(@files, [ 'etc/plots.conf', genome_data() ]);



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
	libCPT => $libCPT,
);
$crr_output->CRR(data => $archive);


=head1 NAME

PSM Plotter

=head1 DESCRIPTION

Following the execution of the PSM Prep tool, this tool plots a subset of those genomes as ciruclar tracks with protein-protein relationships plotted as links between the boxes representing proteins.

=head2 IMPORTANT PARAMETERS

=over 4

=item C<evalue>, C<dice>

Adjusting these parameters will affect which links are plotted. Links are heatmapped into bins based on dice score as that is the easiest measure to work with, and scales nicely from 0 to 100. For example, a link with a dice score of 20-29 would be plotted as 20% black (grey20), whereas a dice score of 90+ would be plotted as solid black

=item C<mismatch>, C<gap_penalty>, C<match>

These parameters control the Needleman-Wunsch Multiple Sequence Alignment library's scoring scheme. Mismatch scores are generally negative and discourage unrelated proteins from being plotted in a line together. Match scores encourage related proteins to line up. Gap penalty is set at zero as we generally prefer gaps to mismatches in this tool; phage genomes are small and gaps are "cheap" to use, whereas mismatches can sometimes give an incorrect impression of relatedness. That said, how your plots look is completely up to you and we encourage experimentation!

=item C<every_nth>

Every Nth gene in a genome will be plotted a slightly different color.

=back

=head2 Why Can't I Control Colors?

    Brewer colors compose Brewer palettes which have been manually defined by
    Cynthia Brewer for their perceptual properties.

    http://circos.ca/tutorials/lessons/configuration/colors/

Color palette choice is very important to creating an attractive and easy to read graphic. In the words of Dr. Krzywinski, L<Color palettes matter|http://mkweb.bcgsc.ca/jclub/biovis/brewer/colorpalettes.pdf>. Humans selecting from an RGB/HSV color palette tend to make poor choices, so we've removed the option in lieu of using the very attractive L<Brewer Palettes|http://colorbrewer2.org/>. Specifically, I've selected the 8 class qualtitative "Accent" color set, which has produced some very nice maps. If you would like the option of selecting amongst the other 8-class qualitative color sets, please file a bug report and let me know.


=cut
