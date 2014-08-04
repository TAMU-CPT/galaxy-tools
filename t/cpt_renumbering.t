#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 8;
use IPC::Run3 qw(run3);

my ( @base, @cmd, $in, $out, $err );

@base = ('perl', 'bin/cpt_renumbering.pl');
my %result_files = (
  "Default" => {
    command_line => "--file test-data/inputs/single.gbk --string_prefix gp ",
    outputs => {
      "change_table" => ["renumbered.Sheet1.csv", "test-data/outputs/renumbered.csv"],
      "genbank" => ["renumbered.gbk", "test-data/outputs/renumbered.gbk"],
    },
  },
  "cluster by tag 'locus_tag'" => {
    command_line => "--clustered_by_tag locus_tag --file test-data/inputs/single.gbk --string_prefix gp ",
    outputs => {
      "change_table" => ["renumbered.Sheet1.csv", "test-data/outputs/renumbered.cluster.csv"],
      "genbank" => ["renumbered.gbk", "test-data/outputs/renumbered.cluster.gbk"],
    },
  },
);

foreach ( keys(%result_files) ) {
  # run with the command line
  my @cmd1 = ( @base, split( / /, $result_files{$_}{command_line} ) );
  run3 \@cmd1, \$in, \$out, \$err;
  if($err){ print STDERR "Exec STDERR: $err"; }
  if($out){ print STDERR "Exec STDOUT $out"; }
  # and now compare files
  foreach my $file_cmp ( keys( %{$result_files{$_}{outputs}} ) ) {
    my ($gen, $static) = @{$result_files{$_}{outputs}{$file_cmp}};
    my @diff = ( "diff", $gen, $static );
    my ($in_g, $out_g, $err_g);
    run3 \@diff, \$in_g, \$out_g, \$err_g;
    if($err_g) { print STDERR "err_g $err_g\n"; }
    if($out_g) { print STDOUT "out_g $out_g\n"; }
    chomp $out_g;
    is( -e $gen, 1, "[$_] Output file must exist"); 
    is( length($out_g), 0, "[$_] Checking validity of output '$file_cmp'" );
    unlink $gen;
  }
}