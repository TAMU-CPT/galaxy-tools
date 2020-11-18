#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Carp;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		['ldap_url'                   , 'LDAP Server URL'                      , {required => 1 , validate => 'String'    , default => 'ldaps://00-ldap-biobio.tamu.edu'}] ,
		['ldap_base'                  , 'LDAP Server Base'                     , {required => 1 , validate => 'String'    , default => 'ou=People,dc=tamu,dc=edu'}] ,
		['genbank_submission_title'   , 'Genbank Submission Title'             , {required => 1 , validate => 'String'}]                               ,
		['genbank_record_author'      , 'Genbank record authors (LDAP DN)'     , {required => 1 , validate => 'String'    , multiple => 1  }]                              ,
		['genbank_record_contact'     , 'Genbank record contact (LDAP DN)'     , {required => 1 , validate => 'String'    , multiple => 1  }]                              ,
		['genbank_record_affiliation' , 'Genbank record affiliation (LDAP DN)' , {required => 1 , validate => 'String' }] ,
		['paper_submission_title'     , 'Paper Submission Title'               , {required => 1 , validate => 'String' }] ,
		['paper_record_author'        , 'Paper record authors (LDAP DN)'       , {required => 1 , validate => 'String'    , multiple => 1  }]                              ,
		['paper_record_contact'       , 'Paper record contact (LDAP DN)'       , {required => 1 , validate => 'String'    , multiple => 1  }]                              ,
		['paper_record_affiliation'   , 'Paper record affiliation (LDAP DN)'   , {required => 1 , validate => 'String' }] ,
		['genbank_file'               , 'Genbank File'                         , {required => 1 , validate => 'String' }] ,
	],
	'outputs' => [
		[
			'sequin',
			'Sequin Output',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'sequin',
				data_format    => 'archive',
				default_format => 'zip'
			}
		],
		[
			'discrep',
			'NCBI tbl2asn Discrepencies',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'discrep',
				data_format    => 'text/plain',
				default_format => 'TXT'
			}
		],
	],
	'defaults' => [
		'appid'   => 'sequin',
		'appname' => 'Sequin',
		'appdesc' => 'generates files which can be sent to NCBI for genome submission',
		'appvers' => '1.94',
	],
	# Test cases would go here....if I had one.
	# This code is harder to test, because tbl2asn includes date-specific
	# information in one of this files it generates. As a result, the test
	# case would go "out of date" every day and be completely meaningless.
	#
	# If I could write a test case, the options would look like:
	#
	# perl sequin.pl
	#	--genbank_submission_title "Karma"
	#	--genbank_record_author "cn=Eric Rasche,ou=staff,ou=People,dc=tamu,dc=edu"
	#	--genbank_record_contact "cn=Eric Rasche,ou=staff,ou=People,dc=tamu,dc=edu"
	#	--genbank_record_affiliation "ou=CPT,ou=CPT,ou=People,dc=tamu,dc=edu"
	#	--paper_submission_title "Some Paper Title"
	#	--paper_record_author "cn=Eric Rasche,ou=staff,ou=People,dc=tamu,dc=edu"
	#	--paper_record_contact "cn=Eric Rasche,ou=staff,ou=People,dc=tamu,dc=edu"
	#	--paper_record_affiliation "ou=CPT,ou=CPT,ou=People,dc=tamu,dc=edu"
	#	--genbank_file ../t/test-files/Karma.gbk
	#
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
    people => {},
    ou => {},
);

# Get People
$mesg = $ldap->search(   # perform a search
    base   => $options->{ldap_base},
    filter => "objectClass=inetOrgPerson"
);


for ( my $i = 0 ; $i < $mesg->count ; $i++ ) {
    my $entry = $mesg->entry($i);
    my %ldap = ();
    foreach my $i(qw(mail givenName sn generationQualifier initials)){
        $ldap{$i} = $entry->get_value($i);
        if($ldap{$i}){
            $ldap{$i} =~ s/\$/\n/g;
        }
    }
    $ldap_data{people}{$entry->dn()} = \%ldap;
}

# Get orgs/groups
$mesg = $ldap->search(   # perform a search
    base   => $options->{ldap_base},
    filter => "objectClass=organizationalUnit"
);

for ( my $i = 0 ; $i < $mesg->count ; $i++ ) {
    my $entry = $mesg->entry($i);
    my %ldap = ();
    foreach my $i(qw(ou l postalAddress postalCode st street telephoneNumber)){
        $ldap{$i} = $entry->get_value($i);
        if($ldap{$i}){
            $ldap{$i} =~ s/\$/\n/g;
        }
    }
    $ldap_data{ou}{$entry->dn()} = \%ldap;
}
#############################################################
# Create the FASTA File
#############################################################
my $genbank_displayid;
open(my $fsa,'>',"$temp/$submission_title.fsa");
use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio_object = $bio->getSeqIO($options->{genbank_file});

while(my $seq_object = $seqio_object->next_seq){
	$genbank_displayid = $seq_object->display_id();
	my $sequence = $seq_object->seq();
	$sequence =~ s/(.{80})/$1\n/g;
	printf $fsa ">%s [organism=Phage %s][gcode=11]\n", $options->{genbank_submission_title}, $options->{genbank_submission_title};
	print $fsa $sequence;
	print $fsa "\n";
}
close($fsa);
#############################################################
# Create the TBL file
#############################################################
open(my $tbl,'>',"$temp/$submission_title.tbl");

printf $tbl ">Feature %s\n",$options->{genbank_submission_title};
my $seqio_object_gbk = $bio->getSeqIO($options->{genbank_file});
while(my $seq_object = $seqio_object_gbk->next_seq){
	for my $feat_object ($seq_object->get_SeqFeatures) {
		my $pk = $feat_object->primary_tag;
		if($pk ne 'source'){
			my $start = $feat_object->start;
			my $end = $feat_object->end;
			if($feat_object->strand == -1){
				printf $tbl "%s\t%s\t" , $end, $start;
			}else{
				printf $tbl "%s\t%s\t" , $start, $end;
			}
			print $tbl $pk,"\n";
			for my $tag ($feat_object->get_all_tags) {
				for my $value ($feat_object->get_tag_values($tag)) {
                    printf $tbl "\t\t\t%s\t%s\n", $tag, $value;
				}
			}
		}
	}
}

close($tbl);


#############################################################
# Create the SBT File
#############################################################
open(my $sbt,'>',"$temp/$submission_title.sbt");


my $contact_string = create_contact_string(4, $options->{genbank_record_contact});
my $affil_string  = create_affil_string(4, $options->{genbank_record_affiliation});
my $affil_string2 = create_affil_string(5, $options->{paper_record_affiliation});
my $record_author_str = create_contact_string(4,$options->{genbank_record_author});
my $paper_author_str  = create_contact_string(6,$options->{paper_record_author});

my $publication_status = "unpublished";
my $publication_title = $options->{paper_submission_title};

# Now we need to generate the template file
my $TBL_0 = <<TBL;
Submit-block ::= {
  contact {
    contact {
      name name {
$contact_string
      },
      affil std {
$affil_string
      }
    }
  },
  cit {
    authors {
      names std {
$record_author_str
      },
      affil std {
$affil_string
      }
    }
  },
  subtype new
}
TBL

print $sbt $TBL_0;

my $TBL_2 = <<TBL;
Seqdesc ::= pub {
  pub {
    gen {
      cit "$publication_status",
      authors {
        names std {
          {
$paper_author_str
          }
        },
        affil std {
$affil_string2
        }
      },
      title "$publication_title"
    }
  }
}
TBL
print $sbt $TBL_2;


close($sbt);
#############################################################
# Run tbl2asn
#############################################################
# Depends on the debian package ncbi-tools-bin
my @command = ("tbl2asn",'-p', $temp, '-M', 'n','-Z',  $temp . '/discrep' , '-r', $temp);
use IPC::Run3;
my ($in,$out,$err);
run3 \@command, \$in, \$out, \$err;
# Hiding b/c it always prints "[tbl2asn] Validating YYY"
if($err){
	#print STDERR $err;
}
use File::Copy;
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'discrep',
	GGO => $ggo,
);
my $discrep_location = $crr_output->CRR(data => "none", extension => "txt");
copy("$temp/discrep", $discrep_location);

#############################################################
# Post process
#############################################################
use Archive::Any::Create;
my $archive = Archive::Any::Create->new();
# Put in a named folder sequin_Karma
my $new_name = sprintf('SQN_%s%s',
	$options->{genbank_submission_title} ,
	substr($file,length($temp))
);
$archive->container($new_name);
my @files = glob("$temp/*");
foreach my $file(@files){
	#my ($location, $content) = @{$_};
	my $content;
	open(my $fh, '<', $file);
	while(<$fh>){ $content .= $_ }
	close($fh);
	$archive->add_file($content);
}

use CPT::OutputFiles;
my $crr_output2 = CPT::OutputFiles->new(
	name => 'sequin',
	GGO => $ggo,
);
$crr_output2->CRR(data => $archive);
rmdir($temp);

#############################################################
# Helper Functions
#############################################################
sub create_contact_string{
	my ($offset,$author_dn_ref) = @_;
	my @dns = @{$author_dn_ref};
	my @output;
	foreach my $dn(@dns){
		my %data = ();
		$data{last} = $ldap_data{people}{$dn}{'sn'};
		$data{first} = $ldap_data{people}{$dn}{'givenName'};

		# Initials SUCK.
		my $unparsed_suffix = $ldap_data{people}{$dn}{'initials'};
		my @letters;
		if($unparsed_suffix){
		    $unparsed_suffix =~ s/[^A-Za-z]//g; # Remove everything that isn't a capital letter
		    @letters = split(//, $unparsed_suffix); # split into array
		}
		unshift @letters, substr($data{first},0,1); #Prepend first letter of first name
		$data{initials} = join('.', map{uc($_)}@letters);
		$data{suffix} = $ldap_data{people}{$dn}{'generationQualifier'};
		foreach my $key(keys %data){
		    unless($data{$key}){
			$data{$key} = "";
		    }
		}
		my $str =
					'  ' x $offset . "{\n" .
					'  ' x ($offset+1) . "name name {\n" .
					create_string(($offset+2), \%data, [qw(last first initials suffix)]) . "\n" .
					'  ' x ($offset+1) . "}\n" .
					'  ' x $offset . "}";
		push(@output,$str);
	}
	return	join( ",\n" , @output	);
}
sub create_affil_string{
	my ($offset,$affil_dn) = @_;
	if($ldap_data{ou}{$affil_dn}){
		my %ldap_hit = %{$ldap_data{ou}{$affil_dn}};
		my %affil_info = (
			affil   => $ldap_hit{ou},
			div     => '',
			city    => $ldap_hit{l},
			'sub'   => $ldap_hit{st},
			country => 'USA',
			street  => $ldap_hit{street},
			email   => 'cpt@tamu.edu', #TODO, get this into the LDAP database
			fax     => '',
			phone   => $ldap_hit{postalCode},
			zip     => $ldap_hit{telephoneNumber},
		);
		return '  ' x $offset . "{\n" .
					'  ' x ($offset+1) . "name name {\n" .
					create_string(($offset+2),\%affil_info, [keys(%affil_info)]) . "\n" .
					'  ' x ($offset+1) . "}\n" .
					'  ' x $offset . "}";
    }else{
        die 'You selected a person (cn= entry) as the affiliation. This is incorrect. Please select an organisation (ou= entry)';
    }
}
sub create_string{
	my ($length,$hash_ref,$arr_ref) = @_;
	my %hash = %{$hash_ref};
	my @arr = @{$arr_ref};
	#return join(",\n",map { if($hash{$_}){'  ' x $length . $_ . ' "' . $hash{$_} . '"'}else{()}  } @arr);
	return join(",\n",map { '  ' x $length . $_ . ' "' . ($hash{$_} ) . '"' } @arr);
}


=head1 WHY YOU CARE

Have you ever had to submit more than one GenBank file? Have you ever had to submit dozens? Do you worry that you might make mistakes? This tool avoids B<ALL OF THAT> by letting you select the name from our directory, and have the person's information filled in correctly and automatically. This tool, on average, saves about 20-30 minutes of work due simply to how easy it is to fill out a couple text boxes, rather than the complex and I<error prone> workflow normal to Sequin users.

=head1 DESCRIPTION

Generate Sequin/tbl2asn output for submission to genbank. Files are created as an archive, which can be downloaded and edited further. If the organisation is NOT the CPT, B<PLEASE> submit your organisation's contact information and address and we will add it to the directory.

=cut
