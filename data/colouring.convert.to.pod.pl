use strict;
use warnings;
use CPT;
my $libcpt = CPT->new();
my %hash = %{$libcpt->JSONYAMLopts(file=> 'colouring.yaml')};
foreach(sort(keys(%hash))){
	printf "=item %s\n\n", $hash{$_}{title};
	printf "    Contains text: \"" . join('", "', @{$hash{$_}{members}}) ."\"\n\n";
	if(defined($hash{$_}{custom})){
		printf "    Special rules: \n";
		foreach my $special(keys $hash{$_}{custom}){
			printf "        - \"%s\" but not \"%s\"\n", join('" or "', @{$hash{$_}{custom}{$special}{is}}), join('" or "', @{$hash{$_}{custom}{$special}{isnot}});
		}
	}
	my ($r,$g,$b)= split / /, $hash{$_}{color};

	printf "Color: L<rgb(%s,%s,%s)|https://duckduckgo.com/?q=rgb+%s+%s+%s>\n\n", $r,$g,$b,$r,$g,$b;
}
