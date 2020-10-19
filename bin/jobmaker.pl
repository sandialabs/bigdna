#! /usr/bin/perl
use strict; use warnings; use List::Util 'shuffle';

my $genomespath = '/projects/islands/hybrid';
die "Usage: perl $0 <RUN.settings file>\n" unless $ARGV[0] and -f $ARGV[0];
my ($setsfile, $jobsfile) = ($ARGV[0], $ARGV[0]);  $jobsfile =~ s/\.settings$/.jobs/; die "$setsfile, $jobsfile\n" if $setsfile eq $jobsfile;

mkdir 'runs';
mkdir 'sums';
mkdir 'finalReports';
open IN, $setsfile;
my $header = <IN>; $header =~ s/^#//; chomp $header;
my ($i, %cats, @jobs, @sums) = (0);
for (split /\t/, $header) {$cats{$_} = $i; $i ++}
while (<IN>) {
#Run QUICK_P3 OVERLAP_MIN_SIZE OVERLAP_MAX_SIZE OPTIMIZE SOLUTION TNT_USE PRIMER_NUM_RETURN Dataset End
 chomp; my @f = split "\t"; next if /^#/;
 my $sets = "TYPE=SETTING\n";
 for (qw/OPTIMIZE OVERLAP_MAX_SIZE OVERLAP_MIN_SIZE PRIMER_NUM_RETURN SOLUTION TNT_USE QUICK_P3/) {
  # if QUICK_P3=basic: hiL22,loT57,hiT63,hiH47,hiM5,loG20
  $sets .= "$_=$f[$cats{$_}]\n";
 } 
 my ($file, @isles) = ($f[$cats{Dataset}]); #, '3082.gff');
 mkdir "runs/$f[$cats{Run}]";
 for (`cat sequences/$file.gff`) {  # 3082 with 3-10 frags
  my @g = split "\t";
  /fasta=([^;]+)/; my $fasta = $1;
  my $annot = $fasta; $annot =~ s/genome.fa/isleannots.gff/;
  /ID=([^;]+)/; my $id = $1;
  $id =~ s/[\/\+\|]/_/g;  # Weird character
  my $out = $sets;
  $out .= "=\nTYPE=SEGMENT\nNAME=$id\nFASTA=$fasta\nENTRY=$g[0]\nL=$g[3]\nR=$g[4]\nORIENT=$g[6]\n";
  if ($f[$cats{Dataset}] eq 'eighths') {
   $out .= "LINEAR=1\n";
  }
  if ($f[$cats{End}] !~ /[^\d]/) {$out .= "L_TOLERANCE=$f[$cats{End}]\nR_TOLERANCE=$f[$cats{End}]\n"}
  elsif ($f[$cats{End}] eq 'gene') {$out .= "L_GENE=1\nR_GENE=1\nANNOT=$annot\n"}
  elsif ($f[$cats{End}] eq 'delta') {$out .= "DELTA_INT=1\nANNOT=$annot\n"}
  $out .= "=\n";
  mkdir "runs/$f[$cats{Run}]/$id";
  open OUT, ">runs/$f[$cats{Run}]/$id/config"; print OUT $out; close OUT;
  push @jobs, "perl bin/bigdna.pl runs/$f[$cats{Run}]/$id/config &> runs/$f[$cats{Run}]/$id/err\n";
 }
 push @sums, "tail -q -n 1 runs/$f[$cats{Run}]/\*/bigdna.log > sums/$f[$cats{Run}]\n";
}
#for (`cat jobs/indels.jobs`) {push @jobs, $_}
open JOBS, ">$jobsfile"; print JOBS join('', shuffle @jobs); close JOBS;
$jobsfile =~ s/jobs$/tails/;
open SUM, ">$jobsfile"; print SUM join('', @sums); close SUM;
