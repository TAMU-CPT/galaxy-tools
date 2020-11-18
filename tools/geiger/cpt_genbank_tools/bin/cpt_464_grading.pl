#!/usr/bin/perl
use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file', 'Input Student Genbank file',
			{
				validate => 'File/Input',
				#file_format => ['Genbank'],
				required => 1,
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
	],
	'outputs' => [
		[
			'corrections',
			'List of Good/Bad Annotations',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'corrections',
				data_format    => 'text/html',
				default_format => 'HTML',
			}
		],
	],
	'defaults' => [
		'appid'   => '464.GradingAndCorrections',
		'appname' => '464 Grading and Correction Utility',
		'appdesc' => 'makes numerous corrections to the genome and scores the student\'s work',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

my $only_ipr_regex = qr/^Interpro:IPR[0-9]+$/;
my $only_ipr_regex_bad = qr/^IPR[0-9]+$/;
my $possible_mult_ipr = qr/[^:]+IPR[0-9]+$/;
my $possible_mult_ipr2 = qr/.+Interpro:IPR[0-9]+$/;
my $tmhelix_regex = qr/^[0-9]+TMD(s?) \((?<nums>[0-9, -]+)\) N[ -](in|out), C[ -](in|out)/;
my $tmhelix_regex_maybe = qr/[0-9]+.*TMD.*\(.*\).*N.*(in|out).*C.*(in|out)/;
my $signal_regex = qr/signal.*cleavage.*/;
my $signal_regex_good = qr/signal peptidase .* cleavage site [0-9]+-[0-9]+/;




use CPT::Report;
use CPT::Report::HTML;
my $report = CPT::Report::HTML->new();

use Bio::SeqIO;
use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqio_object = $bio->getSeqIO($options->{file});
# For all genomes in the GBK file
while(my $seq_object = $seqio_object->next_seq){
	# For RBS Scoring
	my @rbs_map;
	foreach my $rbs($seq_object->get_SeqFeatures()){
		if($rbs->primary_tag eq 'RBS'){
			push(@rbs_map, $rbs);
		}
	}

	# Actual Scoring Process
	my %scores;
	foreach my $feat($seq_object->get_SeqFeatures()){
		if($feat->primary_tag eq 'CDS'){
			# InterPro
			my ($good_ref, $bad_ref) = iprscore($feat);
			push(@{$scores{InterPro}{good}}, @{$good_ref});
			push(@{$scores{InterPro}{bad}}, @{$bad_ref});

			# RBS
			my @rbs = locate_upstream_rbs($feat,\@rbs_map);
			if(@rbs && scalar @rbs == 1 && defined $rbs[0]){
				push(@{$scores{RBS}{good}}, sprintf("RBSs found, %s nt away from %s [%s]",  inter_feature_spread($feat, $rbs[0]), $bio->_getIdentifier($feat), $rbs[0]->seq->seq));
			}else{
				my $i = 0;
				foreach(@rbs){
					if(defined $_ && $_){
						$i++;
					}
				}
				push(@{$scores{RBS}{bad}}, sprintf("%s RBSs found for %s", $i, $bio->_getIdentifier($feat)));
			}

			# TMD
			($good_ref, $bad_ref) = tmd_score($feat);
			push(@{$scores{TMDs}{good}}, @{$good_ref});
			push(@{$scores{TMDs}{bad}}, @{$bad_ref});

			# Signal (SignalP, LipoP)
			($good_ref, $bad_ref) = signal_score($feat);
			push(@{$scores{SignalSequences}{good}}, @{$good_ref});
			push(@{$scores{SignalSequences}{bad}}, @{$bad_ref});
		}

		if($feat->primary_tag eq 'terminator'){
			my ($good_ref, $bad_ref) = terminators($feat);
			push(@{$scores{Terminators}{good}}, @{$good_ref});
			push(@{$scores{Terminators}{bad}}, @{$bad_ref});
		}
		if($feat->primary_tag eq 'tRNA'){
			my ($good_ref, $bad_ref) = trnas($feat);
			push(@{$scores{tRNAs}{good}}, @{$good_ref});
			push(@{$scores{tRNAs}{bad}}, @{$bad_ref});
		}
	}

	$report->h1("Report for " . $seq_object->display_id());
	foreach(sort(keys(%scores))){
		$report->h2($_);
		$report->h4("Good Annotations");
		my @good = @{$scores{$_}{good}};

		$report->list_start('bullet');
		foreach my $item(@good){
			$report->list_element($item);
		}
		$report->list_end();

		$report->h4("Bad Annotations");
		my @bad = @{$scores{$_}{bad}};

		$report->list_start('bullet');
		foreach my $item(@bad){
			$report->list_element($item);
		}
		$report->list_end();
	}
}

sub terminators {
	my ($f) = @_;
	my @good;
	my @bad;
	my @tags = $f->get_tag_values();
	if(scalar @tags > 0){
		push(@bad, sprintf("Terminator: %s had %s tags associated with it (%s)", $bio->_getIdentifier($f), scalar @tags, join(',',@tags)));
	}else{
		push(@good, sprintf("Terminator: %s was good", $bio->_getIdentifier($f)));
	}
	return (\@good,\@bad);
}

sub trnas {
	my ($f) = @_;
	my @good;
	my @bad;
	foreach my $tag($f->get_all_tags()){
		if($tag ne 'evidence' && $tag ne 'product'){
			push(@bad, sprintf("tRNA %s has bad tag associated with it %s", $bio->_getIdentifier($f), $tag));
		}
	}


	if(scalar(@bad) == 0){
		push(@good, sprintf("tRNA %s was good", $bio->_getIdentifier($f)));
	}

	return (\@good,\@bad);
	#tRNA            complement(101183..101259)
	#/evidence="tRNAscan-SE"
	#/product="tRNA-Pro"

}


sub signal_score {
	my ($f) = @_;
	my @good;
	my @bad;
	foreach my $tag($f->get_all_tags()){
		foreach my $value($f->get_tag_values($tag)){
			# If we have something that looks like an IPR number
			if($value =~ $signal_regex){
				if($tag eq 'signal' && $value =~ $signal_regex_good){
					push(@good, sprintf('%s="%s"', $tag, $value));
				}else{
					push(@bad, sprintf('%s="%s"',$tag, $value));
				}
				return (\@good, \@bad);
			}
		}
	}
	# If we'd seen it otherwise, we would've returned by now
	if($f->has_tag('signal')){
		foreach my $value($f->get_tag_values('signal')){
			push(@bad, sprintf('%s="%s"','signal', $value));
		}

	}
	return (\@good, \@bad);
}

sub tmd_score {
	my ($f) = @_;
	my @good;
	my @bad;
	foreach my $tag($f->get_all_tags()){
		foreach my $value($f->get_tag_values($tag)){
			# If we have something that looks like an IPR number
			if($value =~ $tmhelix_regex){
				if($tag eq 'tmhelix'){
					push(@good, sprintf('%s="%s"', $tag, $value));
				}else{
					push(@bad, sprintf('%s="%s"',$tag, $value));
				}
				return (\@good, \@bad);
			}
			if($value =~ $tmhelix_regex_maybe){
				push(@bad, sprintf('%s="%s"',$tag, $value));
				return (\@good, \@bad);
			}
		}
	}
	# If we'd seen it otherwise, we would've returned by now
	if($f->has_tag('tmhelix')){
		foreach my $value($f->get_tag_values('tmhelix')){
			push(@bad, sprintf('%s="%s"','tmhelix', $value));
		}

	}
	return (\@good, \@bad);

	#/tmhelix="2TMD (7-29 34-56), N-in, C-in"
}

sub iprscore {
	my ($f) = @_;
	my @good;
	my @bad;
	foreach my $tag($f->get_all_tags()){
		foreach my $value($f->get_tag_values($tag)){
			# If we have something that looks like an IPR number
			if($value =~ $only_ipr_regex){
				if($tag eq 'dbxref'){
					push(@good, sprintf('db_xref="%s"',$value));
				}else{
					push(@bad, sprintf('%s="%s"',$tag, $value));
				}
			}
			if($value =~ $only_ipr_regex_bad){
				push(@bad, sprintf('%s="%s"',$tag, $value));
			}
			# If we have something that looks like multiple dbxrefs on a single line
			if($value =~ $possible_mult_ipr){
				push(@bad, sprintf('%s="%s"',$tag, $value));
				my @tags = split(/[\s,]+/, $value);
			}
			if($value =~ $possible_mult_ipr2){
				push(@bad, sprintf('%s="%s"',$tag, $value));
			}
		}
	}
	return (\@good, \@bad);
}

sub locate_upstream_rbs {
	my ($feat, $rbs_ref) = @_;
	my @rbs_list = @{$rbs_ref};
	my @close;
	foreach my $rbs(@rbs_list){
		if(is_close($feat, $rbs)){
			push(@close, $rbs);
		}
	}
	if(scalar @close > 1){
		warn "Found MORE THAN ONE possible RBS for " . $bio->_getIdentifier($feat) . ": " . join("\n", map { "\t" . $bio->_getIdentifier($_) } @close ) . "\n";
	}
	return $close[0];
}

sub is_close{
	my ($feat_a, $feat_b) = @_;
	if($feat_a->strand() != $feat_b->strand()){
		return 0;
	}
	if(inter_feature_spread($feat_a, $feat_b) < 15){
		return 1;
	}
	return 0;
}

sub inter_feature_spread {
	my ($feat_a, $feat_b) = @_;
	if($feat_a->strand() != $feat_b->strand()){
		warn "Cannot compare distance for features on different strands";
		return -1;
	}
	if($feat_a->strand() == 1){
		return abs($feat_a->start() - $feat_b->end()) - 1;
	}
	if($feat_a->strand() == -1){
		return abs($feat_a->end() - $feat_b->start()) - 1;
	}
	return abs($feat_a->start() - $feat_b->end());
}







use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
	name => 'corrections',
	GGO => $ggo,
);
$csv_output->CRR(data => $report->get_content());

=head2 README

This is very much a prototype. If you see "incorrect" calls from this tool, PLEASE report them to me so they can be patched!

=head2 Functionality

This tool runs numerous checks over the correctness of annotations. These checks are mostly related to formatting, rather than an actual assessment of the annotation quality, as that tool is still in development.

=head3 InterPro

C<db_xref> fields which look like they contain an IPR number are checked for formatting. Specifically this should be C<Interpro:IPR######>. C<db_xref>s which look like C<IPR######> will be failed.

=head3 RBSs

If an RBS is annotated within 15 bases upstream of a gene, that's marked as good. If there are zero, or more than one RBS annotated, it's marked as bad.

=head3 TMDs

This check is for formatting of TMDs. TMDs should be annotated like the following:

    4TMDs (12-34, 56-72) N-in, C-out
    1TMD (12-34) N out, C in

the dashes were made optional, as is the 's' following 'TMD'

=head3 SignalP/LipoP

If a note contains C<signal.*cleavage>, then it was checked to make sure it was of the form:

    signal peptidase .* cleavage site ##-##

=head3 Terminators

All information must be removed from terminators.

=head3 tRNAs

Tags other than C<evidence> and C<product> must be removed from tRNAs.

=cut
