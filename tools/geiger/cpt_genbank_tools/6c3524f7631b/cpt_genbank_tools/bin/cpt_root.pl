#!/usr/bin/env perl
#
#       Code written by Eric Rasche 
#               mailto:rasche.eric@yandex.ru
#               tel:404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#
$|=1;
use strict;
use warnings;

use CPT;
use Bio::Tools::Run::StandAloneBlastPlus;
use DateTime::Format::Strptime;
use Bio::DB::EUtilities;
use Bio::DB::EUtilities;
use File::Copy;
use Bio::SeqIO;
use CPT::OutputFiles;


my $libCPT = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[
			'protein_f' => 'Input file containing a single protein sequence',
			{
				validate    => 'File/Input',
			}
		],
		[
			'protein_t' => 'Pasted single protein sequence',
			{
				validate    => 'String',
			}
		],
		[ 'name' => 'Name of the protein' ,  { required => 1, validate => 'String' }],
	],
	'outputs' => [
		[
			'results',
			'Root Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'root',
				data_format    => 'text/html',
				default_format => 'HTML'
			}
		],
	],
	'defaults' => [
		'appid'   => 'Root',
		'appname' => 'Root',
		'appdesc' => 'ferrets out original annotation of a given protein',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

if(!defined($options->{protein_f}) && !defined($options->{protein_t})){
	die 'Please specify protein in either file or text format';
}
use File::Temp qw/tempfile/;
my $template = 'cpt.root.XXXXXXX';

my %ids;
my @acc_sorted;
my $strp = DateTime::Format::Strptime->new(
    pattern   => '%d-%b-%Y',
    locale    => 'en_US',
    time_zone => 'America/Chicago',
);
my ($protein_query_fh, $protein_query) = tempfile($template);
if($options->{protein_t}){
	open(my $fh, '>', $protein_query);
	my $ps = $options->{protein_t};
	if($ps !~ /^>/){
		$ps = ">seq\n$ps";
	}
	print $fh $ps;
	close($fh);
}else{
	copy($options->{protein_f}, $protein_query);
}


	print STDERR "Searching\n";
my ($hist, @ids) = search($options->{name});
my $count = scalar @ids;
print STDERR "Fetch\n";
#my $protein_results = fetch($hist);
my $protein_results = 'cpt.root.LRXcXNl';
print STDERR "Parse\n";
my $merged_fasta = parse($protein_results,$count);
print STDERR "Blast\n";
blast($merged_fasta);
print STDERR "Reporting\n";
my $pandoc_report = report();

pandoc($pandoc_report);

sub blast{
	my ($merged_fasta) = @_;
	my $dir = File::Temp->newdir(TEMPLATE=>$template);
	my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
		-DB_DIR => $dir->dirname,
		-db_data => $merged_fasta,
		-create => 0,
	);
	$fac->make_db;
	#my ($fh = File::Temp->new($template, UNLINK => 1, SUFFIX => '.xml');
	my ($blast_results_fh, $blast_results) = tempfile($template, SUFFIX => '.xml');

	my $result = $fac->blastp(
		-query     => $protein_query,
		-outfile   => $blast_results,
		-outformat => 5
	);
	my @hits = $result->hits;
	foreach my $hit(@hits){
		while( my $hsp = $hit->next_hsp ) {
			my $Identity             = $hsp->percent_identity * $hsp->length('query');
			my $dice                 = (2*$Identity)/($hit->length + $result->query_length());
			$ids{$hit->accession}{blast} = [
			        $hsp->evalue,
			        $dice,
			        $hsp->seq_str('query'),
				$hsp->seq_str('match'),
				$hsp->seq_str('hit'),
			];
		}
	}
}
sub parse{
	my ($protein_results, $count) = @_;
	my $c = 0;
	my $seqio_object = Bio::SeqIO->new(-file =>$protein_results);
	my ($merged_fasta, $merged_fasta_filename) = tempfile($template);

	while(my $seq_object = $seqio_object->next_seq){
		my $anno_collection = $seq_object->annotation;
		my $acc = $seq_object->accession;
		$ids{$acc} = ();
		for my $key ( $anno_collection->get_all_annotation_keys ) {
			my @annotations = $anno_collection->get_Annotations($key);
			my $ref_count=0;
			for my $value ( @annotations ) {
				if($value->tagname eq 'date_changed'){
					#print "tagname : ", $value->tagname, "\n";
					#my $dt = $strp->parse_datetime($value->display_text);
					#print "  annotation value: ", $value->display_text, "\t";
					#print $dt->strftime("%Y%m%d %H:%M %Z")."\n";
					$ids{$acc}{'submission_date'} = $strp->parse_datetime($value->display_text);
				}
				if($value->tagname eq 'reference'){
					$ref_count++;
					if($ref_count == 1){
						$ids{$acc}{'first_reference'} = $value->hash_tree;
					}
				}
			}
		}
		my @genbank_text;
		for my $feat_object ($seq_object->get_SeqFeatures) {
			my $feature_text = sprintf("    %-20s",$feat_object->primary_tag);
			$feature_text .= sprintf("%d..%d\n",$feat_object->start,$feat_object->end);
			for my $tag ($feat_object->get_all_tags()){
				for my $value ($feat_object->get_tag_values($tag)){
					$feature_text .= "    " . " " x 20 . "/$tag=\"$value\"\n";
				}
			}
			push(@genbank_text,$feature_text);
			if ($feat_object->primary_tag eq "CDS") {
				print $merged_fasta '>' . $seq_object->accession,"\n";
				my $seq = $feat_object->seq->seq;
				$seq =~ s/(.{80})/$1\n/g;
				print $merged_fasta $seq,"\n";
			}
		}
		$ids{$acc}{'genbank_text'} = join("\n",@genbank_text);
	}
	close($merged_fasta);
	return $merged_fasta_filename;
}
sub search{
	my ($query) = @_;
	my $esearch = Bio::DB::EUtilities->new(
		-eutil      => 'esearch',
		-db         => 'protein',
		-term       => $query,
		-retmax     => 10000,
		-email      => 'cpt@tamu.edu',
		-usehistory => 'y',
	);
	my @ids = $esearch->get_ids;
	return ($esearch->next_History , @ids);
}
sub fetch{
	my $hist = shift;
	my $factory = Bio::DB::EUtilities->new(
	        -eutil   => 'efetch',
	        -db      => 'protein',
	        -rettype => 'gb',
	        -email   => 'cpt@tamu.edu',
		-history => $hist,
		-usehistory => 'y'
	);
	my ($fh, $filename) = tempfile( $template, UNLINK => 0 );
	print STDERR "STORING TO $filename\n";
	$factory->get_Response(-file=>$filename);
	return $filename;
}
sub pandoc{
	my ($pandoc_report_location) = @_;
	my $html_output = CPT::OutputFiles->new(
		name => 'results',
		libCPT => $libCPT,
	);
	my ($filename) = $html_output->CRR(data =>'');

	system("pandoc -t html5 --standalone --css='https://cpt.tamu.edu/md.css' $pandoc_report_location > $filename");
}
sub report{
	my $result;
	$result .= "# Summary Information\n";
	my %date_data;
	foreach(keys %ids){
		$date_data{$ids{$_}{submission_date}->strftime("%Y")}++
	}
	my @years = sort(keys %date_data);
	$result .="We found " . scalar(keys %ids) . " record(s), the first being in **" . (sort(keys %date_data))[0] . "** \n\nHistogram of Dates\n\n";
	
	foreach($years[0]..$years[-1]){
		if($date_data{$_}){
			$result .=  "    $_\t" . '*' x int(log(int($date_data{$_}))) . "\n";
		}
		else{
			$result .=  "    $_\t\n";
		}
	}
	$result .=  "\n\n";

	$result .=  "# Full Report\n";
	@acc_sorted = sort { DateTime->compare($ids{$a}{'submission_date'},$ids{$b}{'submission_date'}) } keys %ids;
	foreach my $acc (@acc_sorted){
		$result .=  "## $acc (" . $ids{$acc}{submission_date}->strftime("%Y-%m-%d") . ")\n";
		my $refref = $ids{$acc}{'first_reference'};
		if($refref && ref($refref) eq 'HASH'){
			$result .= "### Genbank Record\n\n";
			$result .= $ids{$acc}{'genbank_text'} . "\n\n";
			$result .=  "### Literature Reference\n\n";
			my ($authors,$journal,$title,$pubmed) = map { $ids{$acc}{'first_reference'}{$_} } qw(authors location title pubmed);
			unless($authors){$authors=" "}
			if($pubmed){
			$result .=  <<MD;
Title

:   [$title](http://www.ncbi.nlm.nih.gov/pubmed/$pubmed)

MD
			}else{
			$result .=  <<MD;
Title

:   $title

MD

			}

			$result .=  <<MD;
Authors

:   $authors

Journal

:   $journal

MD
		}
		my $blastref = $ids{$acc}{'blast'};
		if($blastref){
			my @blast = @{$ids{$acc}{'blast'}};
			$result .=  "\n### Blast Data\n";
			$result .=  <<MD;
Evalue

:   $blast[0]

Dice

:   $blast[1]

Match Sequence

:       $blast[2]
        $blast[3]
        $blast[4]

MD
		}
	}

	my ($pandoc_report_fh, $pandoc_report) = tempfile($template);
	print $pandoc_report_fh $result;
	return $pandoc_report;
}

=head1 DESCRIPTION

This tool attempts to find the first instance of a protein being annotated as "X", where X is your search query. If you were interested in the first time someone annotated a "Large Terminase" as a "Large Terminase", you would query on that. It is recommended you add an actual terminase sequence so as to identify which annotations are "true" annotations to your protein of interest and not misannotation by other researchers.

=cut
