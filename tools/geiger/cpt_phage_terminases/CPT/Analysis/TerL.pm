package CPT::Analysis::TerL;

# ABSTRACT: Guess phage packaging strategy based on homology to terminases (TerL) of phages with known packaging strategies

use strict;
use warnings;
use Data::Dumper;
use autodie;
use Moose;
use Bio::SearchIO;
use File::ShareDir;
use File::Spec;


has 'hmmer_evalue_cutoff' => (is => 'rw', isa => 'Int');
has 'blast_evalue_cutoff' => (is => 'rw', isa => 'Int');
has 'blast_dice_cutoff' => (is => 'rw', isa => 'Int');

has 'search_hmmer' => (is => 'rw', isa => 'Bool');
has 'search_blast' => (is => 'rw', isa => 'Bool');

has 'data_dir' => (is => 'rw', isa => 'Any');

has 'awk_string' => (is => 'ro', isa => 'Str', default => 'BEGIN{print "#row,query id,subject id,evalue,dice" } {qid=$1;sid=$2;percent_identical=$3;query_length=$4;subject_length=$5;evalue=$6;dice=(2*percent_identical*subject_length/ ( subject_length + query_length ) );printf("%d,%s,%s,%s,%s\n",FNR, qid, sid, dice, evalue);}');

has 'input_file' => (is => 'rw', isa => 'Str');

my @column_max = (20,10,20,10,20,60);
my %hits;

sub run{
	my ($self, %data) = @_;
	#$self->prepare(%data);
	my $dir = File::ShareDir::dist_dir('CPT-Analysis-TerL');
	$self->data_dir($dir);
	$self->input_file($data{input_file});


	if($self->search_hmmer()){
	       $self->hmmer();
	}
	if($self->search_blast()){
		$self->blast($self);
	}
	
	return $self->guess();
}

sub hmmer{
	my ($self) = @_;

	use Term::ProgressBar 2.00;
	my $progress = Term::ProgressBar->new({
		name => 'HMMER',
		count => 100,
		ETA => 'linear',
	});
	$progress->max_update_rate(1);
	
	my $hmmer_db_dir = File::Spec->catdir($self->data_dir(),'db','hmmer');
	my @hmmer_dbs = glob(File::Spec->catfile($hmmer_db_dir,"*.hmm"));

	my $i=0;
	foreach(@hmmer_dbs){
		$progress->update( 100* ($i++) / scalar @hmmer_dbs );
		my $db_name = substr($_,rindex($_,'/')+1,-4);
		my $output_filename = sprintf('%s.%s.out', substr($self->input_file(),0,rindex($self->input_file(),'.')), $db_name);
		my $query = sprintf('hmmsearch %s %s > %s', $_, $self->input_file(), $output_filename);
		system($query);

		my $in = Bio::SearchIO->new(-format => 'hmmer',
			-file   => $output_filename);
		while( my $result = $in->next_result ) {
			while( my $hit = $result->next_hit ) {
				while( my $hsp = $hit->next_hsp ) {
					my ($from, $to) = ($result->query_name,$hit->name);
					unless($hits{$to}{$from}){
						$hits{$to}{$from} = [];
					}
					push($hits{$to}{$from},
						{
							'type' => 'hmmer',
							'data' => {
								'evalue' => ($hsp->evalue eq '0'? '0.0' : $hsp->evalue ),
							}
						}
					);
				}
			}
		}
		unlink($output_filename);
	}
	$progress->update(100);
}

sub blast{
	my ($self) = @_;
	
	use Term::ProgressBar 2.00;
	my $progress = Term::ProgressBar->new({
		name => 'Blast',
		count => 100,
		ETA => 'linear',
	});
	$progress->max_update_rate(1);

	my $blast_db_dir = File::Spec->catdir($self->data_dir(),'db','blast');
	my @blast_dbs = map { substr($_,0,-4) }glob(File::Spec->catfile($blast_db_dir,"*.phr"));

	my $i=0;
	#my %hits;

	foreach my $blast_db(@blast_dbs){
		$progress->update(100 * ($i++) / scalar(@blast_dbs) );
		my $output_str = substr($blast_db,rindex($blast_db,'/')+1) . '.csv';
		my $query = sprintf('psiblast -query %s -db %s -evalue %s -outfmt "6 qseqid sseqid pident qlen slen evalue" | awk \'%s\' > %s',
				$self->input_file(), $blast_db, '1e-5', $self->awk_string(),$output_str);
		system($query);
		open(my $tmpfh, '<', $output_str);
		while(<$tmpfh>){
			chomp $_;
			if($_ !~ /^#/){
				my @line = split(/,/,$_);
				unless($hits{$line[1]}{$line[2]}){
					$hits{$line[1]}{$line[2]} = [];
				}
				push($hits{$line[1]}{$line[2]},
					{
						'type' => 'psiblast',
						'data' => {
							'evalue' => $line[4],
							'dice' => $line[3],
						}
					}
				);
			}
		}
		close($tmpfh);
		unlink($output_str);
		
	}
	$progress->update(100);
}

sub guess{
	my ($self) = @_;

	open(my $groupings_fh,'<',File::Spec->catfile($self->data_dir(),'groupings.tsv'));
	
	# Load groupings.tsv into memory
	my %data;
	while(<$groupings_fh>){
		if($_ !~ /^#/){
			chomp $_;
			my ($major, $minor, $hit) = split(/\t/,$_);
			unless($data{$major}{$minor}){
				$data{$major}{$minor} = []
			}
			push($data{$major}{$minor}, $hit);
		}
	}

	# Create a reverse lookup table
	my %rdata;
	foreach my $i(keys %data){
		foreach my $j ( keys %{ $data{$i}}){
			foreach my $k(@{ $data{$i}{$j}}){
				if(defined($k)){
					$rdata{$k} = [$i,$j];
				}else{
					print "$i $j\n";
				}
			}
		}
	}


	# Table printing stuff
	my @header = ('Major Category', 'Major hits', 'Minor Category', 'Minor hits', 'Analysis', 'Evidence Type', 'Evidence');
	my %output;

	# Loop across the input keys
	foreach my $input_key(keys %hits){
		my %guesses;
		my %guess_evidence;
		# And across all of the hits that the query hit to
		foreach my $against ( keys %{$hits{$input_key}}){
			# We look at the evidence
			my @evidence = @{$hits{$input_key}{$against}};

			my ($type_major, $type_minor);
			if($rdata{$against}){
				($type_major, $type_minor) = @{$rdata{$against}};
			}else{
				($type_major, $type_minor) = (substr($against,0,rindex($against,'_')) , substr($against,rindex($against,'_')+1));
			}
			# Prepare hashes.
			unless($guesses{$type_major}){
				$guesses{$type_major} = ();
			}
			unless($guesses{$type_major}{$type_minor}){
				$guesses{$type_major}{$type_minor} = ();
			}
			# Loop across the evidence
			foreach my $piece_of_evidence(@evidence){
				# Here is an example piece of evidence
				#    'GK_Gilmour_Gene43' => {
				#        'SP18' => [
				#             {
				#               'data' => {
				#                           'evalue' => '2e-08',
				#                           'dice' => '23.8627'
				#                         },
				#               'type' => 'psiblast'
				#             }
				#           ],
				#
				if($type_major !~ /subject i/){
					my %piece = %{$piece_of_evidence};
					if($self->validate_evidence($piece_of_evidence) > 0 ){
						$guess_evidence{$type_major}++;
						$guess_evidence{"$type_major$type_minor"}++;

						# If it's not defined, set up sub arrays.
						unless($guesses{$type_major}{$type_minor}{$piece{'type'}}){
							if($piece{'type'} eq 'psiblast'){
								$guesses{$type_major}{$type_minor}{$piece{'type'}} = {'evalue' => [], 'dice'=>[] };
							}else{
								$guesses{$type_major}{$type_minor}{$piece{'type'}} = {'evalue' => [] };
							}
						}
						# If it's zero, correct to zero.
						if($piece{'data'}{'evalue'} eq '0.0'){
							$piece{'data'}{'evalue'} = '0.0';
						}
						# Add our evalue
						push($guesses{$type_major}{$type_minor}{$piece{'type'}}{'evalue'}, $piece{'data'}{'evalue'});
						# And if psiblast, add dice
						if($piece{'type'} eq 'psiblast'){
							push($guesses{$type_major}{$type_minor}{$piece{'type'}}{'dice'}, $piece{'data'}{'dice'});
						}
					}
				}
			}
		}

		my @output_sheet;

		foreach my $major(keys %guesses){
			if($guess_evidence{$major}){
				foreach my $minor(keys %{$guesses{$major}})
				{
					if($guess_evidence{"$major$minor"}){
						foreach my $evidence_category(keys %{$guesses{$major}{$minor}}){
							if($evidence_category ne 'evidence'){
								# things like evalue, dice
								my %hits = %{$guesses{$major}{$minor}{$evidence_category}};
								foreach my $subtype(keys %hits){ # should be evalue, dice
									push(@output_sheet,
										[$major, $guess_evidence{$major}, $minor, $guess_evidence{"$major$minor"},
											$evidence_category, $subtype, join(',', @{$hits{$subtype}})]
									);
								}
							}
						}
					}
				}
			}
		}


		if(!scalar @output_sheet){
			@output_sheet = (
				['No evidence above threshold'],
			);
		}
		$output{$input_key} = {
			header => \@header,
			data => \@output_sheet,
		};
	}
	return \%output;
}

sub validate_evidence{
	my ($self, $piece_of_evidence) = @_;
	my %piece = %{$piece_of_evidence};
	#my ($self, $type, $subtype, $value) = @_;
	#             {
	#               'data' => {
	#                           'evalue' => '2e-08',
	#                           'dice' => '23.8627'
	#                         },
	#               'type' => 'psiblast'
	#             }
	if($piece{type} eq 'hmmer'){
		my $value = $piece{data}{evalue};
		if($value eq '0.0' ||  $value eq '0'){
			return 1;
		}elsif(! defined($value) ||  $value eq ''){
			return 0;
		}else{
			return (log($value) < $self->hmmer_evalue_cutoff()); # -64
		}
	}
	elsif($piece{type} eq 'psiblast'){
		my $evalue = $piece{data}{evalue};
		my $dice = $piece{data}{dice};

		my $evalue_return = 0;
		if($evalue eq '0.0' || $evalue eq '0'){
			$evalue_return = 1;
		}else{
			$evalue_return = (log($evalue) < $self->blast_evalue_cutoff()); #-140
		}

		return $evalue_return && ($dice > $self->blast_dice_cutoff()); # 30
	}
}

1;
