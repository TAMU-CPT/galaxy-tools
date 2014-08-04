#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 6;
use IPC::Run3 qw(run3);

my ( @base, @cmd, $in, $out, $err );

@base = ('perl', 'bin/cpt_gc.pl');
my %result_files = (
  "Default GBK" => {
    command_line => "--file test-data/inputs/multi.gbk ",
    outputs => {
      "results" => ["gc.Sheet1.csv", "test-data/outputs/gc.default.csv"],
    },
  },
  "Default GBK CDSs" => {
    command_line => "--file test-data/inputs/multi.gbk --tag CDS ",
    outputs => {
      "results" => ["gc.Sheet1.csv", "test-data/outputs/gc.cds.csv"],
    },
  },
  "Default .fa" => {
    command_line => "--file test-data/inputs/multi.cds.fa ",
    outputs => {
      "results" => ["gc.Sheet1.csv", "test-data/outputs/gc.fa.csv"],
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