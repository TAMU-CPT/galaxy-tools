#!/usr/bin/perl
use strict;
use CPT::GalaxyGetOpt;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[
			'blast', 'Input Blast 25 Column file',
			{
				validate => 'File/Input',
			}
		],
	],
	'outputs' => [
		[
			'blast_report',
			'Blast Report',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'blast_report',
				data_format    => 'text/html',
				default_format => 'HTML',
			}
		],
	],
	'defaults' => [
		'appid'   => 'BlastViz',
		'appname' => 'Blast Viz',
		'appdesc' => 'plots blast alignments in an easy to understand manner',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);


my %data;
my %lengths;
my %max_bitscore;

open(my $blast_input, '<', $options->{blast});
while (<$blast_input>) {
	next if /^\#/;  # ignore comments
	chomp;
	my ($qseqid, $sseqid, $pident, $length, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore, $sallseqid, $score, $nident, $positive, $gaps, $ppos, $qframe, $sframe, $qseq, $sseq, $qlen, $slen, $salltitles) = split /\t+/;
	push( @{$data{$qseqid}{$sseqid}}, [$qstart, $qend, int($bitscore)]);
	$lengths{$qseqid} = $qlen;
	$lengths{$sseqid} = $slen;

	
	$max_bitscore{$qseqid}{$sseqid} = $bitscore if not defined $max_bitscore{$qseqid};
	if($bitscore > $max_bitscore{$qseqid}{$sseqid}){
		$max_bitscore{$qseqid}{$sseqid} = $bitscore;
	}
}
close($blast_input);

use File::Copy;
use File::Temp qw/ tempdir tempfile /;
use File::Spec::Functions qw/catfile catdir/;
my @files_to_move;
my $dir = tempdir('cpt.viz.blast.XXXXXXX',CLEANUP => 1 );
my $blast_results = '';

foreach my $genome(keys %data){
	print STDERR "Processing $genome\n";
	$blast_results .= '<h2>' . $genome . "</h2>\n";
	my $panel = Bio::Graphics::Panel->new(
					      -length    => $lengths{$genome},
					      -width     => 1000, # 1kpx wide
					      -pad_left  => 10,
					      -pad_right => 10,
					     );
	my $full_length = Bio::SeqFeature::Generic->new(
							-start => 1,
							-end   => $lengths{$genome},
						);
	$panel->add_track($full_length,
			  -glyph   => 'arrow',
			  -tick    => 2,
			  -fgcolor => 'black',
			  -double  => 1,
			 );
	my %hit_ids;
	for my $subject(keys $data{$genome}){
		my $track = $panel->add_track(
					      -glyph     => 'graded_segments',
					      -label     => 1,
					      -bgcolor   => 'blue',
					      -key       => $subject,
					      -min_score => 0,
					      -max_score => $max_bitscore{$genome}{$subject},
					     );
		foreach my $row(@{$data{$genome}{$subject}}){
			my ($qstart, $qend, $bitscore) = @{$row};
			if($qend-$qstart > 5000){
				my $feature = Bio::SeqFeature::Generic->new(
								      -display_name => $subject,
								      -score        => $bitscore,
								      -start        => $qstart,
								      -end          => $qend
								     );
				$track->add_feature($feature);
				$hit_ids{$subject}++;
			}
		}
	}
	my ( $fh, $path ) = tempfile(
		'XXXXXXXXX',
		UNLINK => 1,
		DIR    => $dir,
		SUFFIX => '.png'
	);
	open(my $panel_out, '>', $path);
	print $panel_out $panel->png;
	close($panel_out);
	my $safe_fn = $genome;
	$safe_fn =~ s/[^A-Za-z0-9]//g;
	push(@files_to_move, [$path, "$safe_fn"]);
	$blast_results .= "<img src=\"$safe_fn.png\" />\n";

	$blast_results .= "<ul>";
	foreach(sort{$hit_ids{$b} <=> $hit_ids{$a}} keys(%hit_ids)){
		if($_ =~ /gi\|([0-9]+)\|/){
			$blast_results .= sprintf('<li><a target="_blank" href="http://www.ncbi.nlm.nih.gov/nuccore/%s">%s (%s)</a></li>', $1, $_, $hit_ids{$_});
		}
	}
	$blast_results .= "</ul>";
}

my $html = sprintf( '
<!DOCTYPE html>
<html>
<body>
<h1>Blast Results</h1>
%s
</body>
</html>
',
	$blast_results,
);

use CPT::OutputFiles;
my $blast_output = CPT::OutputFiles->new(
	name => 'blast_report',
	GGO => $ggo,
);
$blast_output->CRR(data => $html);

#We prepared an array of files to move with their original location and a new
#filename. Extensions are assumed to be PNG in this case.
foreach my $move(@files_to_move){
	my ($original, $new) = @{$move};
	my @produced_files = $blast_output->subCRR(data_format => "Dummy", format_as => "Dummy", filename=>$new, extension => 'png');
	move($original,$produced_files[0]);
}
