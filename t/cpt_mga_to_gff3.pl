#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 1;
use IPC::Run3 qw(run3);

my ( @base, @cmd, $in, $out, $err );

@base = ('perl', '-MDevel::Cover', 'bin/cpt_mga_to_gff3.pl');
my %result_files = (
  "Default" => {
    command_line => "--file test-data/inputs/mga.txt ",
    outputs => {
      "results" => ["mga.gff3", "test-data/outputs/mga.gff3"],
    },
  },
);

foreach ( keys(%result_files) ) {
  # run with the command line
  my @cmd1 = ( @base, split( / /, $result_files{$_}{command_line} ) );
  run3 \@cmd1, \$in, \$out, \$err;
  if($err){ print STDERR "Exec STDERR: $err"; }
  # and now compare files
  foreach my $file_cmp ( keys( %{$result_files{$_}{outputs}} ) ) {
    my ($gen, $static) = @{$result_files{$_}{outputs}{$file_cmp}};
    my @md5_generated = ( "md5sum", $gen );
    my @md5_static = ("md5sum", $static);
    my ($in_g, $out_g, $err_g);
    my ($in_s, $out_s, $err_s);
    run3 \@md5_generated, \$in_g, \$out_g, \$err_g;
    run3 \@md5_static, \$in_s, \$out_s, \$err_s;
    if($err_g) { print STDERR "err_g $err_g\n"; }
    if($err_s) { print STDERR "err_s $err_s\n"; }
    chomp $out_g;
    chomp $out_s;
    $out_g =~ s/ .*//;
    $out_s =~ s/ .*//;
    is( $out_g, $out_s, "[$_] Checking validity of output '$file_cmp'" );
    unlink $gen;
  }
}