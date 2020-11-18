#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Carp;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		['ldap_url' , 'LDAP Server URL' , {required => 1 , validate => 'String'    , default => 'ldaps://00-ldap-biobio.tamu.edu'}] ,
		['ldap_base', 'LDAP Server Base', {required => 1 , validate => 'String'    , default => 'ou=People,dc=tamu,dc=edu'}] ,
		['mode', 'Selection mode', { required => 1, validate => 'Option', options => {'cn', 'Search for CNs, i.e., people', 'ou', 'Only list OUs (groups/organizationalUnits)'}, default => 'cn' } ],
	],
	'outputs' => [
		[
			'sequin_options',
			'Converted LDAP tree for use in sequin as options in a tree-select',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'sequin',
				data_format    => 'Dummy',
				default_format => 'Dummy'
			}
		],
	],
	'defaults' => [
		'appid'   => 'galaxyhelper.sequin_options',
		'appname' => 'Galaxy Helper: Sequin Options',
		'appdesc' => 'generates .loc file for cpt_sequin.pl',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);


use File::Temp;
my $temp = File::Temp::tempdir(TEMPLATE=>'cpt_sequin.XXXXXXXXX',DIR=>'/tmp/',CLEANUP=>0);
my $submission_title = $options->{genbank_submission_title};
#############################################################
# Grab LDAP Data
#############################################################
use Net::LDAPS;
my $ldap = Net::LDAPS->new($options->{ldap_url}) or die "$@";
my $mesg = $ldap->bind; # an anonymous bind
my %ldap_data = (
);

# Get people and OUs
$mesg = $ldap->search(   # perform a search
	base   => $options->{ldap_base},
	filter => "(|(objectClass=organizationalUnit)(objectClass=inetOrgPerson))"
);


use Data::Diver qw( Dive DiveRef DiveVal DiveError );
for ( my $i = 0 ; $i < $mesg->count ; $i++ ) {
	my $entry = $mesg->entry($i);
	my %ldap = ();
	
	my $dn = $entry->dn();
	my @dn_parts = split /,/, $dn;
	my @arr = reverse(@dn_parts[1..scalar(@dn_parts)-1]);
	
	# If the location isn't yet defined, create it.
	if(not defined DiveVal(\%ldap_data, @arr)){
		${DiveRef( \%ldap_data, @arr)} = {};
	}

	# If it's a person (cn=), then set to 1.
	if($dn_parts[0] =~ /^cn=/ && $options->{mode} eq 'cn'){
		${DiveRef( \%ldap_data, @arr, $dn_parts[0])} = 1;
	}

}

use CPT::OutputFiles;
my $crroutput = CPT::OutputFiles->new(
        name => 'sequin_options',
        GGO => $ggo,
);
my @output_files = $crroutput->CRR(extension=> "loc", data => " ");
my $output_file = $output_files[0];

use IO::File;
open(my $output, '>', $output_file);
print $output '#' . join("\t", ('id', 'Name', 'dn')) . "\n";
recurse($output, split(/,/,$options->{ldap_base}));
sub recurse{
	my ($output, @keys) = @_;
	my $ref = DiveRef(\%ldap_data, reverse(@keys));
	my $id = join('', @keys);
	$id =~ s/[^A-Za-z0-9]//g;
	my $tail = scalar@keys - 1 - scalar(split(/,/, $options->{ldap_base}));
	my $path = join('/', map { (split '=', $_)[1] } @keys[1..$tail]);
	my $display = $keys[0];
	if($keys[0] =~ '^cn'){
		$display = substr($keys[0],3);
	}
	printf $output "%s\t%s [%s]\t%s\n", $id, $display, $path,  join(',', @keys);
	if(ref $ref ne 'SCALAR'){
		# However if they're an OU (i.e, hash b/c cn's are set to 1 (scalar))
		my @k = keys %{${$ref}};
		# Recurse through all of those.
		foreach my $key(@k){
			recurse($output,$key,@keys);
		}
	}
}
close($output);

=head1 NAME

Galaxy Helper - .loc dump of users in LDAP

=head1 DESCRIPTION

This script lives to connect to an LDAP database and then dump a list of users as a .loc file for use in galaxy. This file can then be used to allow users to select an LDAP entry from a dropdown in galaxy.

Once you've generated the file, you'll need to place it in C<$GALAXY_ROOT/tool-data/sequin.loc>

Tools can access this data by declaring parameters like the following:

    <param help="select a user from the LDAP database" name="dn" label="User DN" type="drill_down" multiple="False" display="radio" hierarchy="recurse" from_file="sequin.xml" optional="False"/>

=cut
