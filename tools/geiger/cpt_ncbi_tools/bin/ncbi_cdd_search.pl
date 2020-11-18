#!/usr/bin/env perl

use strict;
use warnings;
use LWP::UserAgent;
use Getopt::Std;
my $bwrpsb = "http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi";

use Data::Dumper;
use CPT::GalaxyGetOpt;
my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Input (protein) fasta file',
			{
				validate => 'File/Input',
				file_format => ['fasta'],
			}
		],
	],
	'outputs' => [
		[
			'output',
			'Raw CDD Output',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'cdd',
				data_format    => 'text/plain',
				default_format => 'TXT',
			}
		],
	],
	'defaults' => [
		'appid'   => 'NCBICDDSearch',
		'appname' => 'NCBI CDD Search',
		'appdesc' => 'Searches the NCBI CDD database',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

#my @queries = <STDIN>;
my @queries;
open(my $fh, '<', $options->{file});
while(<$fh>){
	push(@queries, $_);
}
close($fh);

my $output;

###############################################################################
# submitting the search
###############################################################################
my $rid;
{
	my $browser = LWP::UserAgent->new;
	my $response = $browser->post(
		$bwrpsb,
		[
			'useid1' => "true",
			'maxhit' => "1000",
			'filter' => "true",
			'db'     => "cdd",
			'evalue' => "0.001",
			'cddefl' => "true",
			'qdefl'  => "true",
			'dmode'  => "all",
			'clonly' => "false",
			'tdata'  => "hits",
			( map {; queries => $_ } @queries )
		],
	);
	die "Error: ", $response->status_line
		unless $response->is_success;

	if($response->content =~ /^#cdsid\s+([a-zA-Z0-9-]+)/m) {
		$rid =$1;
		$output .= "Search with Request-ID $rid started.\n";
	} else {
		die "Submitting the search failed,\n can't make sense of response: $response->content\n";
	}
}
###############################################################################
# checking for completion, wait 5 seconds between checks
###############################################################################

$|++;
my $done = 0;
my $status = -1;
while ($done == 0) {
	sleep(5);
	my $browser = LWP::UserAgent->new;
	my $response = $browser->post(
		$bwrpsb,
		[
			'tdata' => "hits",
			'cdsid' => $rid
		],
	);
	die "Error: ", $response->status_line
		unless $response->is_success;

	if ($response->content =~ /^#status\s+([\d])/m) {
		$status = $1;
		if ($status == 0) {
			$done = 1;
			$output .= "Search has been completed, retrieving results ..\n";
		} elsif ($status == 3) {
		} elsif ($status == 1) {
			die "Invalid request ID\n";
		} elsif ($status == 2) {
			die "Invalid input - missing query information or search ID\n";
		} elsif ($status == 4) {
			die "Queue Manager Service error\n";
		} elsif ($status == 5) {
			die "Data corrupted or no longer available\n";
		}
	} else {
		die "Checking search status failed,\ncan't make sense of response: $response->content\n";
	}

}

###############################################################################
# retrieve and display results
###############################################################################
{
	my $browser = LWP::UserAgent->new;
	my $response = $browser->post(
		$bwrpsb,
		[
			'tdata'  => "hits",
			'cddefl' => "true",
			'qdefl'  => "true",
			'dmode'  => "all",
			'clonly' => "false",
			'cdsid'  => $rid
		],
	);
	die "Error: ", $response->status_line
		unless $response->is_success;

	$output .= $response->content . "\n";
}


use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'output',
	GGO => $ggo,
);
$crr_output->CRR(data => $output);
