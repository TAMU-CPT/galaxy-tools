package CPT::Analysis::PPEPP;
use strict;
use warnings;
use Moose;
use IPC::Run3;
use File::ShareDir;
use File::Spec;

# ABSTRACT: Runs phage promoter prediction based on models from canonical phages

has 'dblocation'     => (is => 'rw', isa => 'Any');
has 'lookahead_dist' => (is => 'rw', isa => 'Num');
has 'evalue'         => (is => 'rw', isa => 'Num');
has 'inc'            => (is => 'rw', isa => 'Num');
has 'f1'             => (is => 'rw', isa => 'Num');
has 'f2'             => (is => 'rw', isa => 'Num');
has 'f3'             => (is => 'rw', isa => 'Num');
has 'domz'           => (is => 'rw', isa => 'Num');

my %hash = (
	"186" => {
		"186_P12" => "TAGCGAxxxxxxxxxxxxxxxxxxGTAGCC",	
		"186_PJ" => "TAGCCCxxxxxxxxxxxxxxxxxxGACAAT",	
		"186_PV" => "TGGTCTxxxxxxxxxxxxxxxxxxGGAAAA",	
	},
	"P1" => {
		"P1_P3cre" => "TTGATCxxxxxxxxxxxxxxxxxxxTAAAAT",
		"P1_P2cre" => "TTGATCxxxxxxxxxxxxxxxxxxxCATTAT",
		"P1_P1cre" => "TGTACAxxxxxxxxxxxxxxxxxxxCATTAT",
		"P1_P2mod_(PM2)" => "ATGCATxxxxxxxxxxxxxxxxxxxTAAACT",
		"P1_PkilA_(P53)" => "TGGATCxxxxxxxxxxxxxxxxxxxTCTAAT",
		"P1_PaskilA_(P53as)" => "TCGACAxxxxxxxxxxxxxxxxxxxTATATT",
		"P1_PparB" => "AAGATTxxxxxxxxxxxxxxxxxxxGACCGT",
		"P1_Ppar" => "TGTACAxxxxxxxxxxxxxxxxxxxTATTAT",
		"P1_PrepA" => "TTGCTAxxxxxxxxxxxxxxxxxxxTAGGAT",
		"P1_Plpa_(Pr94)" => "TTGCTAxxxxxxxxxxxxxxxxxxxTATTAT",
		"P1_Pc1_(P99a)" => "TCTATTxxxxxxxxxxxxxxxxxxxTATAAT",
		"P1_P1coi_(P99d)" => "TTGACAxxxxxxxxxxxxxxxxxxxAGTTTT",
		"P1_P4cre_(P99e)" => "TTGATTxxxxxxxxxxxxxxxxxxxTAAACT",
		"P1_P2ssb_(Pr21)" => "TTGTTAxxxxxxxxxxxxxxxxxxxTACATT",
		"P1_Pdmt" => "TTGTTTxxxxxxxxxxxxxxxxxxxTACCAT",
		"P1_P3pmgT" => "TTGATTxxxxxxxxxxxxxxxxxxxGAAAAT",
		"P1_P2c4_(P51a)" => "TTGTATxxxxxxxxxxxxxxxxxxxTATCTT",
		"P1_P1c4_(P51b)" => "TTGACCxxxxxxxxxxxxxxxxxxxTATAGT",
		"P1_P1ban" => "TTGCTCxxxxxxxxxxxxxxxxxxxTAATAT",
	},
	"P2" => {
		"P2_PF" => "TAGCCTxxxxxxxxxxxxxxxxxxxGAAAAT",
		"P2_PO" => "TGGCGGxxxxxxxxxxxxxxxxxxxGGAAAC",
		"P2_PP" => "TAGCGAxxxxxxxxxxxxxxxxxxxGTAGCC",
		"P2_PV" => "TAGCATxxxxxxxxxxxxxxxxxxxGCAATC",
	},
	"P4" => {
		"P4_PLL" => "TGTTGCxxxxxxxxxxxxxxxxxxxACAAAA",
		"P4_PSid" => "TCGTGTxxxxxxxxxxxxxxxxxxxCACAAT",
	},
	"T3" => {
		"T3_A0" => "TTGACTxxxxxxxxxxxxxxxxxxxTATCCT",
		"T3_A1" => "TTGACTxxxxxxxxxxxxxxxxxxxTATTAT",
		"T3_A2" => "TTGACAxxxxxxxxxxxxxxxxxxxTAAGAT",
		"T3_A3" => "TTGACAxxxxxxxxxxxxxxxxxxxTACGAT",
		"T3_B" => "ATGACCxxxxxxxxxxxxxxxxxxxTGTACT",
		"T3_C" => "TTGTCAxxxxxxxxxxxxxxxxxxxTGTAAC",
	},
	"T5" => {
		"T5_(Chinese)_4663-4737" => "TTGACAxxxxxxxxxxxxxxxxxxxTTTAAG",
		"T5_(Chinese)_4795-4880" => "TTGACAxxxxxxxxxxxxxxxxxxxTAAAAT",
		"T5_(Chinese)_5847-5936" => "CTAATGxxxxxxxxxxxxxxxxxxxTAGAGA",
		"T5_(Orsay)_4664-4738" => "TTGACAxxxxxxxxxxxxxxxxxxxTAAGAT",
		"T5_(Orsay)_4810-4867" => "TTGACAxxxxxxxxxxxxxxxxxxxTAAAAT",
		"T5_(Orsay)_5942-6010" => "CTAATGxxxxxxxxxxxxxxxxxxxTAGAGA",
		"T5_(Orsay)_3353-3416" => "TTAATCxxxxxxxxxxxxxxxxxxxTTATTA",
		"T5_(Russian)_4799-4873" => "TTGACAxxxxxxxxxxxxxxxxxxxTAAAAT",
		"T5_(Russian)_5938-6012" => "CTAATGxxxxxxxxxxxxxxxxxxxTAGAGA",
		"T5_(Russian)_4338-4664" => "TTAATAxxxxxxxxxxxxxxxxxxxTGCAAA",
	},
	"T7" => {
		"T7_A1_429-527" => "CTTAAAxxxxxxxxxxxxxxxxxxxTTACAG",
		"T7_A2_545-658" => "TTGACAxxxxxxxxxxxxxxxxxxxTAAGAT",
		"T7_A3_669-779" => "TTGACAxxxxxxxxxxxxxxxxxxxTACGAT",
	}
);
# Temp files
my $tmp_fasta_file = File::Temp->new(
	TEMPLATE => 'ppepp.fasta.tempXXXXX',
	DIR      => '/tmp/',
	UNLINK   => 1,
	SUFFIX   => '.fasta'
);
my $tmp_hmmer_file = File::Temp->new(
	TEMPLATE => 'ppepp.hmmer.tempXXXXX',
	DIR      => '/tmp/',
	UNLINK   => 1,
	SUFFIX   => '.hmmer'
);

use Bio::Seq;
use Bio::SeqIO;

sub getDataForFASTA {
	my ( $self, $gbk, $fasta ) = @_;
	my %fasta_seqs;
	my $seq = Bio::SeqIO->new(
		-file   => $gbk,
		-format => 'genbank'
	);
	my $seq_obj = $seq->next_seq;
	my $count   = 0;
	for my $feat_object ( $seq_obj->get_SeqFeatures ) {
		if ( $feat_object->primary_tag eq "CDS" ) {
			my $upstream = lc( $seq_obj->subseq( $feat_object->start - $self->lookahead_dist()+1, $feat_object->start - 1 ) );
			my $name = "";
			if ( $feat_object->has_tag("gene") ) {
				my @values = $feat_object->get_tag_values('gene');
				$name = $values[0];
			}
			elsif ( $feat_object->has_tag('locus') ) {
				my @values = $feat_object->get_tag_values('locus');
				$name = $values[0];
			}
			else {
				$name = "gp$count";
			}

			my @tmp = (
				$name,    #name
				$feat_object->start,
				$feat_object->end,
				( ( $feat_object->strand == "1" ) ? "1" : "0" ),    #strand
				lc( $seq_obj->subseq( $feat_object->start - $self->lookahead_dist()+1 - 1, $feat_object->start - 1 ) ),
			);
			$fasta_seqs{"$tmp[0]_$tmp[1]_$tmp[2]_$tmp[3]"} = $tmp[4];
			print $fasta ">$tmp[0]_$tmp[1]_$tmp[2]_$tmp[3]\n$tmp[4]\n";
			$count++;
		}
	}
	return %fasta_seqs;
}

sub analyze_sequence{
	my ($self, $file) = @_;
	$self->fix_dblocation();
	my %fasta_seqs = $self->getDataForFASTA($file, $tmp_fasta_file);

	my ($in,$out,$err);
	my @cmd = ('hmmsearch', '--F1',$self->f1(),'--F2', $self->f2(),'--F3',
		$self->f3(),'-E', $self->evalue(), '--incE', $self->inc(),
		'--domZ', $self->domz(), '--notextw', '--tblout',
		$tmp_hmmer_file,$self->dblocation(), $tmp_fasta_file);
	run3 \@cmd, \$in,\$out,\$err;
	if($err){
		print "Tried to execute: " . join(" ", @cmd);
		die "Error running hmmsearch: $err"
	}

	open(my $hmmer_tbl,'<',$tmp_hmmer_file);
	my @output;
	while(<$hmmer_tbl>){
		#print $_;
		chomp $_;
		if(substr($_,0,1) ne '#'){
			push(@output, $_);
		}
	}
	close($hmmer_tbl);

	my $regex = qr/^(?<target>[^ ]*)\s+-\s+(?<query>[^ ]*)\s+-\s+(?<evalue>[0-9.]*)\s+/;
	my %hits_by_model;
	foreach my $row(@output){
	#	# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
	#	#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
	#	NinC_41081_41953_1   -          186                  -                1.7    2.3   0.0         2    2.1   0.0   1.6   1   1   0   1   1   1   0 -
		if($row =~ $regex){
			#printf STDERR "Hit from %s to model %s\n", $+{target}, $+{query};
			my $new_seq = Bio::Seq->new( -display_id =>  $+{target}, -seq => $fasta_seqs{$+{target}});
			push(@{$hits_by_model{$+{query}}},[$new_seq, $row]);
			#print join(" | ", $+{target}, $+{query}, $+{evalue}) , "\n";
		}
	}


	my @data;
	foreach my $model(keys(%hits_by_model)){
		#print "# Hits to $model\n\n";
		my @model_hits = @{$hits_by_model{$model}};
		foreach my $sequence_data(@model_hits){
			my ($sequence, $string) = @{$sequence_data};
			push(@data,parse_cols2($sequence->display_id(), $string));
		}
	}
	return @data;
}

sub parse_cols2{
	my ($id,$string) = @_;
	my @cur_row = split('\s+',$string);
	return \@cur_row;
}

sub fix_dblocation{
	my ($self) = @_;
	if(! defined($self->dblocation()) || !$self->dblocation()){
		my $dir = File::ShareDir::dist_dir('CPT-Analysis-PPEPP');
		my $promoterdb = File::Spec->catfile(
			File::Spec->catdir($dir,'pdbs'),
			'promoterdb'
		);
		$self->dblocation($promoterdb);
	}
}

no Moose;
1;
