#!/usr/bin/perl

use strict;
use warnings;
use autodie;
use CPT::GalaxyGetOpt;
use CPT::External::TMHMM;
use CPT::External::LipoP;
use IPC::Run3;
use File::Temp qw/ tempfile tempdir /;
use File::Copy;
use CGI qw/:standard/;

my $ggo = CPT::GalaxyGetOpt->new();

=head1 NAME

GetSpanin - attempt to identify spanins in a whole genome sequence

=head1 DESCRIPTION

According to Dr. Young:

 > The new FindSpanin program should do this:
 > check all such orfs  for lipobox matching a table of possible Lipoboxes,
 > which can be supplied by Rohit/Manoj and should be settable.
 >
 > Also, if all else fails, it simply reports all orfs that have Cys residues
 > and meet the two criteria.  This is to cover our asses if a new lipobox is
 > present.  Has happened already.
 >
 > Then off it goes.
 >
 > Output would be hits in BOTH lipo search...which should include all possible
 > lipoproteins, including osp and usp candidates, and all possible isp candidates.
 >
 > IF there are isp-osp adjacencies or overlappencies or embeddencies, then the
 > program should cry hallelujah and say Gotcha.


=cut

#
#The new FindSpanin program should do this:
#If off, then check all such orfs for lipobox matching a table of possible Lipoboxes, which can be supplied by Rohit/Manoj and should be settable.
#Also, if all else fails, it simply reports all orfs that have Cys residues and meet the two criteria.  This is to cover our asses if a new lipobox is present.  Has happened already.
#
#Then off it goes.
#
#Output would be hits in BOTH lipo search...which should include all possible lipoproteins, including osp and usp candidates.
#And all possible isp candidates.
#IF there are isp-osp adjacencies or overlappencies or embeddencies, then the program should cry hallelujah and say Gotcha.
#These are my thoughts.

# DONE

#Find every possible start codon (irrespective of shine-dalgarno).   I would take out the choice here, just specify ATG, GTG and TTG for everybody.
#
#Offer ± or both strand search, default ±.
#Offer a parameter about largest possible Lipobox signal length...default 30 (distance from start to first AA of lipobox).  "osp_signal"
#Offer a parameter about longest N terminal cytoplasmic domain for the i-spanin...default 10 (distance to first AA of tmd of ispanin)  "isp-nterm"
#Offer a parameter about smallest ORF length for o-spanin (in codons, make it default 30 )..."osp_min"
#Offer a parameter about smallest ORF length for i-spanin (in codons, make it default 60)..."isp_min"
#Offer a toggle for LipoP.  Default off.  (If on, then submit all orfs that meet osp_signal and osp-min criteria are submitted to LipoP.  But LipoP can be too narrow minded.)
my $options = $ggo->getOptions(
		'options' => [

			[ 'file', 'Input genome sequence as fasta file', { required => 1, validate => 'File/Input', file_format => ['fasta'] } ],
			[ 'strand', 'Which strand to search', { required => 1, validate => 'Option', options => {'s', 'Sense', 'a', 'Antisense', 'b','Both'} } ],
			[],
			['Filtering'],
			[ 'isp_min'    , 'smallest ORF length for i-spanin (in codons)'           , { required => 1 , validate => 'Int' , default => 60 } ] ,
			[ 'isp_nterm_mindist'  , 'Min distance to first AA of tmd of ispanin (in codons)'     , { required => 1 , validate => 'Int' , default => 10 } ] ,
			[ 'isp_nterm_maxdist'  , 'Max distance to first AA of tmd of ispanin (in codons)'     , { required => 1 , validate => 'Int' , default => 30 } ] ,
			[ 'osp_min'    , 'smallest ORF length for o-spanin (in codons)'           , { required => 1 , validate => 'Int' , default => 30 } ] ,
			[ 'osp_signal_mindist' , 'Min distance from start to first AA of lipobox (in codons)' , { required => 1 , validate => 'Int' , default => 10 } ] ,
			[ 'osp_signal_maxdist' , 'Max distance from start to first AA of lipobox (in codons)' , { required => 1 , validate => 'Int' , default => 30 } ] ,
			[],
			['Extra'],
			#[ 'start_codon', 'Which start codons to look for. ATG => ATG, TTG, GTG', { required => 1, validate => 'Option', options => {'ATG', 'ATG/TTG/GTG', 'ACTG', 'ATG/TTG/CTG/GTG'}, default => 'ATG' } ],
			[ 'cys_failover', 'failover to anything with a cys residue', { validate => 'Flag' } ],
			[ 'use_lipop', 'Use LipoP for additional information', { validate => 'Flag' } ],
			[],
			['Results'],
			[ 'max_isp_osp_distance', 'Maximum distance between end of i-spanin, and start of o-spanin', { required => 1, validate => 'Int', default => 10 } ], ],

	'outputs' => [
		[
			'spanin_report',
			'Output HTML File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'report',
				data_format    => 'text/html',
				default_format => 'HTML',
			}
		],
		[
			'putative_isps',
			'Putative ISPs ',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'putative_isps',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
		[
			'putative_osps',
			'Putative OSPs ',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'putative_osps',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
	],
	'defaults' => [
		'appid'   => 'FindSpanin',
		'appname' => 'FindSpanin',
		'appdesc' => 'does what it says on the tin',
	],
	'tests' => [
	],
);

my @putative_osps;
my @putative_isps;

use CPT::OutputFiles;

use CPT::Bio;
my $bio = CPT::Bio->new();
my $seq_obj = ${ $bio->requestCopy( 'file' => $options->{file} ) };

my %table;
use CPT::Bio::ORF;
use Digest::MD5 qw(md5_hex);

my $posps = CPT::Bio::ORF->new(
	'min_gene_length' => $options->{osp_min},
	sc_atg => 1,
	sc_ctg => 0,
	sc_ttg => 1,
	sc_gtg => 1,
);
@putative_osps = $posps->run($seq_obj->seq());
# Copy data to P.ISPs
foreach(@putative_osps){
	push(@putative_isps,$_);
}

# Filter

open( my $fh, '>>', 'hashtable.tsv' );
@putative_isps = filter_isps_for_tmds( @putative_isps );
close($fh);
@putative_osps = filter_osps_for_lipoboxes( @putative_osps );

# Store
save_fasta( 'putative_osps', \@putative_osps );
save_fasta( 'putative_isps', \@putative_isps );
check_candidates( \@putative_isps, \@putative_osps );

sub check_candidates {
	my ( $isp_ref, $osp_ref ) = @_;
	my @isps = @{$isp_ref};
	my @osps = @{$osp_ref};
	my %istarts;
	my %ostarts;
	my %iends;
	my %oends;
	my @pairs;

	#len,start,end,strand,seq
	foreach ( @isps ) {
		#'desc' => '[15925-16112; 63 aa long]',
		if($_->desc =~ /\[(\d+)-(\d+);/){
			my $start = $1;
			push( @{ $istarts{$start} }, $_ );
		}
	}
	foreach ( @isps ) {
		if($_->desc =~ /\[(\d+)-(\d+);/){
			my $end = $2;
			push( @{ $iends{$end} }, $_ );
		}
	}
	foreach ( @osps ) {
		if($_->desc =~ /\[(\d+)-(\d+);/){
			my $start = $1;
			push( @{ $ostarts{$start} }, $_ );
		}
	}
	foreach ( @osps ) {
		if($_->desc =~ /\[(\d+)-(\d+);/){
			my $end = $2;
			push( @{ $oends{$end} }, $_ );
		}
	}

	# Prep HTML output
	my $response = '<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="generator" content="pandoc">
  <title></title>
  <!--[if lt IE 9]>
    <script src="http://html5shim.googlecode.com/svn/trunk/html5.js"></script>
  <![endif]-->
</head>
<body>

<style type="text/css">
a,abbr,acronym,address,applet,article,aside,audio,b,big,blockquote,body,canvas,caption,center,cite,code,dd,del,details,dfn,div,dl,dt,em,embed,fieldset,figcaption,figure,footer,form,h1,h2,h3,h4,h5,h6,header,hgroup,html,i,iframe,img,ins,kbd,label,legend,li,mark,menu,nav,object,ol,output,p,pre,q,ruby,s,samp,section,small,span,strike,strong,sub,summary,sup,table,tbody,td,tfoot,th,thead,time,tr,tt,u,ul,var,video{margin:0;padding:0;border:0}body{font-family:Helvetica,arial,freesans,clean,sans-serif;font-size:14px;line-height:1.6;color:#333;background-color:#fff;padding:20px;max-width:960px;margin:0 auto}body>:first-child{margin-top:0!important}body>:last-child{margin-bottom:0!important}blockquote,dl,ol,p,pre,table,ul{margin:15px 0}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#000}h2{font-size:24px;border-bottom:1px solid #ccc;color:#000}h3{font-size:18px}h4{font-size:16px}h5{font-size:14px}h6{color:#777;font-size:14px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6,body>h1:first-child,body>h1:first-child+h2,body>h2:first-child,body>h3:first-child,body>h4:first-child,body>h5:first-child,body>h6:first-child{margin-top:0;padding-top:0}h1+p,h2+p,h3+p,h4+p,h5+p,h6+p{margin-top:10px}a{color:#4183C4;text-decoration:none}a:hover{text-decoration:underline}ol,ul{padding-left:30px}ol li ol:first-of-type,ol li ul:first-of-type,ol li>:first-child,ul li ol:first-of-type,ul li ul:first-of-type,ul li>:first-child{margin-top:0}ol ol,ol ul,ul ol,ul ul{margin-bottom:0}dl{padding:0}dl dt{font-size:14px;font-weight:700;font-style:italic;padding:0;margin:15px 0 5px}dl dt:first-child{padding:0}dl dt>:first-child{margin-top:0}dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}dl dd>:first-child{margin-top:0}dl dd>:last-child{margin-bottom:0}code,pre,tt{font-size:12px;font-family:Consolas,"Liberation Mono",Courier,monospace}code,tt{margin:0;padding:0;white-space:nowrap;border:1px solid #eaeaea;background-color:#f8f8f8;border-radius:3px}pre>code{margin:0;padding:0;white-space:pre;border:0;background:0 0}pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:0}blockquote{border-left:4px solid #DDD;padding:0 15px;color:#777}blockquote>:first-child{margin-top:0}blockquote>:last-child{margin-bottom:0}hr{clear:both;margin:15px 0;height:0;overflow:hidden;border:0;background:0 0;border-bottom:4px solid #ddd;padding:0}table th{font-weight:700}table td,table th{border:1px solid #ccc;padding:6px 13px}table tr{border-top:1px solid #ccc;background-color:#fff}table tr:nth-child(2n){background-color:#f8f8f8}.markdown-body img{max-width:100%}
input {
	display:none;
}
span#content {
	display:none;
}
input#show:checked ~ span#content {
	display:block;
}
input#hide:checked ~ span#content {
	display:none;
}
</style>
	';

	$response .= h1("Results");

	# actual checking
	foreach my $pmi ( @isps ) {
		my @hits = find_osps( $pmi, \@isps, \@osps, \%ostarts );
		my ($len, $start, $end, $strand, $seq) = convert_to_old_format($pmi);

		my $should_display = 0;
		my $subresponse = '';
		if ( scalar @hits > 0 ) {
			$subresponse = h2("Putative ISP: $pmi\n");
			my @simulated_fasta_thingy;
			push @simulated_fasta_thingy, format_header( start => $start, end => $end, strand => $strand );
			push @simulated_fasta_thingy, $bio->translate($seq);
			$subresponse .= p("ISPs semi-aligned to OSPs");
			foreach my $osp_hit_obj (@hits) {
				my ( $osp_hit, $type ) = @{$osp_hit_obj};
				my ( $olen, $ostart, $oend, $ostrand, $oseq ) = convert_to_old_format($osp_hit);
				# Remove those in the same reading frame
				if(($start-$ostart) % 3 == 0){
					next;
				}
				$should_display = 1;
				#push @simulated_fasta_thingy , sprintf("&gt;%s_%s_%s_%s", $ostart,$oend,$ostrand,$type);
				push @simulated_fasta_thingy, format_header( start => $ostart, end => $oend, strand => $ostrand, type => $type );
				my $diff;
				if ( $strand eq 'F' ) {
					$diff = int( ( $ostart - $start ) / 3 );
				}
				else {
					$diff = int( ( $oend - $end ) / 3 );
				}
				push @simulated_fasta_thingy, " " x $diff . $bio->translate($oseq);
			}
			$subresponse .= p("Information on ISPs and possible OSPs");
			foreach my $osp_hit_obj (@hits) {
				my ( $osp_hit, $type ) = @{$osp_hit_obj};
				my ( $olen, $ostart, $oend, $ostrand, $oseq ) = convert_to_old_format($osp_hit);
			}
			$subresponse .= pre( join( "\n", @simulated_fasta_thingy ) );
		}
		if ($should_display){
			$response .= $subresponse;
		}
	}

	$response .= "</body></html>";
	my $output = CPT::OutputFiles->new(
		name => 'spanin_report',
		GGO => $ggo,
	);
	$output->CRR(data => $response);
}

sub format_header {
	my (%o) = @_;
	if ( $o{strand} eq 'F' ) {
		if ( $o{type} ) {
			return sprintf( '&gt;%s_%s_%s_%s', $o{start}, $o{end}, $o{strand}, $o{type} );
		}
		else {
			return sprintf( '&gt;%s_%s_%s', $o{start}, $o{end}, $o{strand} );
		}
	}
	else {
		if ( $o{type} ) {
			return sprintf( '&gt;%s_%s_%s_%s', $o{end}, $o{start}, $o{strand}, $o{type} );
		}
		else {
			return sprintf( '&gt;%s_%s_%s', $o{end}, $o{start}, $o{strand} );
		}
	}
}

sub bounds {
	my ( $a, $b, $c ) = @_;
	return $a ge $b && $a le $c;
}

sub convert_to_old_format {
	my ($feature ) = @_;
	my ($len, $start, $end, $strand, $seq);
	#'desc' => '[15925-16112; 63 aa long]',
	if($feature->desc =~ /\[(\d+)-(\d+);/){
		$start = $1;
		$end = $2;
		$len = abs($start - $end);
	}
	# display_id' => 'orf00878_R
	if($feature->display_id =~ /(.*)_([RF])/){
		$strand = $2;
	}
	$seq = $feature->seq();
	return $len, $start, $end, $strand, $seq;
}

sub find_osps {
	my ( $pmi, $isps_ref, $osps_ref, $ostarts_ref ) = @_;
	my @isps    = @{$isps_ref};
	my @osps    = @{$osps_ref};
	my %ostarts = %{$ostarts_ref};
	my ( $len, $start, $end, $strand, $seq ) = convert_to_old_format( $pmi );

	my @possible_hits;

	if ( $strand eq 'F' ) {
		for ( my $j = $start + 1 ; $j < $end + $options->{max_isp_osp_distance} ; $j++ ) {
			if ( defined( $ostarts{$j} ) ) {
				my @possible_matching_osps = @{ $ostarts{$j} };
				foreach my $pmo (@possible_matching_osps) {
					my @a = convert_to_old_format($pmo);
					if ( $a[3] eq $strand ) {
						my $type;
						if ( bounds( $a[1], $start, $end ) && bounds( $a[2], $start, $end ) ) {
							$type = "embedded";
						}
						elsif ( bounds( $a[1], $start, $end ) && !bounds( $a[2], $start, $end ) ) {
							$type = "overlapped";
						}
						elsif ( $a[1] > $end ) {
							$type = "adjacent";
						}
						if ($type) {
							push( @possible_hits, [ $pmo, $type ] );
						}
					}
				}
			}
		}
	}
	else {
		for ( my $j = $end - 1 ; $j > $start - $options->{max_isp_osp_distance} ; $j-- ) {
			if ( defined( $ostarts{$j} ) ) {

				#printf STDERR "\t%s: %s\n", $j, $ostarts{$j};
				my @possible_matching_osps = @{ $ostarts{$j} };
				foreach my $pmo (@possible_matching_osps) {
					my @a = convert_to_old_format($pmo);
					if ( $a[3] eq $strand ) {
						my $type;
						if ( bounds( $a[2], $start, $end ) && bounds( $a[1], $start, $end ) ) {
							$type = "embedded";
						}
						elsif ( bounds( $a[2], $start, $end ) && !bounds( $a[1], $start, $end ) ) {
							$type = "overlapped";
						}
						elsif ( $a[2] < $start ) {
							$type = "adjacent";
						}
						if ($type) {
							push( @possible_hits, [ $pmo, $type ] );
						}
					}
				}
			}
		}
	}
	return @possible_hits;
}

sub save_fasta {
	my ($output_name, $data_ref ) = @_;
	my $output = CPT::OutputFiles->new(
		name => $output_name,
		GGO => $ggo,
	);
	$output->CRR(data => $data_ref);
}

sub filter_osps_for_lipoboxes {
	my (@sequences) = @_;
	my @good_sequences;
	printf STDERR "[OSP] Analysing %s sequences\n", scalar(@sequences);

	my @lipoboxes = (
		qr/[ILMFTV][^REKD][GAS]C/    # (L-V-or other large hydrophobic)-X (not charged) – (G or A or S) – C

	);
	for ( my $i = 0 ; $i < scalar @sequences; $i++ ) {
		printf STDERR "\b\b\b\b%3.3d%%", int( $i * 100 / scalar(@sequences) );
		my $individual_seq = $sequences[$i];
		my $protein_sequence = $individual_seq->translate->seq;

		my @lipobox_locations;
		#my @data = @{ $local_pusps{ $keyset[$i] } };

		## LipoP
		if ( $options->{use_lipop} ) {
			my $lipop = CPT::External::LipoP->new();
			$lipop->analyze(
				sprintf(
					">%s\n%s\n",
					$individual_seq->display_id,
					$protein_sequence,
				)
			);
			my $cleavage = $lipop->cleavage();
			if ( defined $cleavage ) {
				my ( $start, $end ) = @{$cleavage};
				push @lipobox_locations, $start;
			}
			$lipop->cleanup();
		}
		## Built in Lipobox list
		foreach my $regex (@lipoboxes) {
			if ( uc( $protein_sequence ) =~ $regex ) {
				# http://stackoverflow.com/questions/87380/how-can-i-find-the-location-of-a-regex-match-in-perl
				push( @lipobox_locations, @- );
				#$should_keep = 1;
			}
		}

		## Failover
		my $should_keep = 0;
		if ( $options->{cys_failover} ) {
			if ( uc( $protein_sequence ) =~ /C/ ) {
				push( @lipobox_locations, @- );

				#$should_keep = 1;
			}
		}

		foreach (@lipobox_locations) {
			if ( $_ < $options->{osp_signal_maxdist} && $_ > $options->{osp_signal_mindist} ) {
				$should_keep = 1;
			}
		}
		if ($should_keep) {
			push(@good_sequences, $individual_seq);
		}

	}
	printf STDERR "%s => %s\n", scalar( @sequences ), scalar( @good_sequences );
	return @good_sequences;
}

sub lookup_prep {
	open( my $fhr, '<', 'hashtable.tsv' );
	while (<$fhr>) {
		chomp $_;
		my ( $md5, $hits, $first ) = split( /\t/, $_ );
		if ( $hits != 0 ) {
			$table{$md5} = [ $hits, $first ];
		}
		else {
			$table{$md5} = [0];
		}
	}
	close($fhr);
	printf STDERR "Loaded %s cached entries\n", scalar( keys(%table) );
}

sub lookup_count {
	my ($md5) = @_;
	if ( $table{$md5} ) {
		my @tmp = @{ $table{$md5} };
		return $tmp[0];
	}
	return;
}

sub lookup_pos {
	my ($md5) = @_;
	if ( $table{$md5} ) {
		my @tmp = @{ $table{$md5} };
		return $tmp[1];
	}
	return;
}

sub filter_isps_for_tmds {
	my (@sequences) = @_;
	my @good_sequences;
	printf STDERR "[ISP] Analysing %s sequences\n", scalar(@good_sequences);

	lookup_prep();
	print STDERR "Prepped lookup table\n";
	for ( my $i = 0 ; $i < scalar @sequences; $i++ ) {
		printf STDERR "\b\b\b\b%3.3d%%", int( $i * 100 / scalar(@sequences) );
		my $individual_seq = $sequences[$i];
		my $protein_sequence = $individual_seq->translate->seq;

		my $hash        = md5_hex( $protein_sequence );
		my $lookup_tmds = lookup_count($hash);                #`perl lookup.pl $hash`;

		my $tmds;
		my $tmhmm;
		# If it's defined in our database of pre-run TMHMM results
		if ( defined $lookup_tmds ) {
			# Get the start location
			my $first_start = lookup_pos($hash);
			if ( defined($first_start)
				&& $first_start < $options->{isp_nterm_maxdist}
				&& $first_start > $options->{isp_nterm_mindist}) {
				push(@good_sequences, $individual_seq);
			}
		}
		else {
			print STDERR "> New sequence, TMHMM\n";
			$tmhmm = CPT::External::TMHMM->new();
			$tmhmm->analyze(
				sprintf(
					">%s\n%s\n",
					$individual_seq->display_id,
					$protein_sequence,
				)
			);
			my $tmds = $tmhmm->num_predicted();
			if ( $tmds > 0 ) {
				my @locations   = @{ $tmhmm->predicted_locations() };
				my $should_keep = 0;
				my @starts;
				foreach (@locations) {
					my ( $start, $end ) = @{$_};
					push( @starts, $start );
				}
				my @sorted_starts = sort { $a <=> $b } (@starts);
				my $first_start = shift @sorted_starts;
				printf $fh "%s\t%s\t%s\n", $hash, $tmds, $first_start;
				if ($first_start < $options->{isp_nterm_maxdist}
					&& $first_start > $options->{isp_nterm_mindist}) {
					push(@good_sequences, $individual_seq);
				}
			}
			else {
				printf $fh "%s\t%s\t%s\n", $hash, $tmds, '';
			}

			$tmhmm->cleanup();
		}

	}
	printf STDERR "%s => %s\n", scalar( @sequences ), scalar( @good_sequences );
	return @good_sequences;
}
