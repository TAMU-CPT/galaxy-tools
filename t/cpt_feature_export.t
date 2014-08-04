#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 10;
use IPC::Run3 qw(run3);

my ( @base, @cmd, $in, $out, $err );

@base = ('perl', 'bin/cpt_feature_export.pl');
my %result_files = (
  "Default DNA export" => {
    command_line => "--file test-data/inputs/single.gbk --tag CDS ",
    outputs => {
      "fasta" => ["export.fa", "test-data/outputs/export.default.fa"],
    },
  },
  "translated protein export" => {
    command_line => "--file test-data/inputs/single.gbk --tag gene --translate  ",
    outputs => {
      "fasta" => ["export.fa", "test-data/outputs/export.translate.fa"],
    },
  },
  "split features" => {
    command_line => "--file test-data/inputs/split-feature.gbk --tag CDS ",
    outputs => {
      "fasta" => ["export.fa", "test-data/outputs/export.split.fa"],
    },
  },
  "split features, translated" => {
    command_line => "--file test-data/inputs/split-feature.gbk --tag CDS --translate  ",
    outputs => {
      "fasta" => ["export.fa", "test-data/outputs/export.splittrans.fa"],
    },
  },
  "Upstream/Downstream" => {
    command_line => "--file test-data/inputs/split-feature.gbk --n_bases_downstream 50 --n_bases_upstream 10 --tag CDS --translate  ",
    outputs => {
      "fasta" => ["export.fa", "test-data/outputs/export.up_down.fa"],
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