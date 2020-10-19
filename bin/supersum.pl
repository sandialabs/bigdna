use strict; use warnings;
die "Usage: perl $0 SETS.settings > SETS.supersum\n" unless $ARGV[0] and -f $ARGV[0];
my ($setsfile, %ccref, %nocc) = ($ARGV[0]);  # 'cc' means the subset of 200 used for exhaustive and finding cognates among 3082
for (`cat sequences/200.gff`) {next unless /ID=([^;]+)/; $ccref{$1} ++}
for (qw/bas21PF bas61PF hiH21PF hiH41PF hiH41UF hiH61PF hiL21PF loT61PF hiL61PF hiLloT21PF hiLloT41PF hiLloT41UF hiLloT61PF hiM21PF hiM41PF hiM41UF hiM61PF hiT21PF hiT41PF hiT41UF hiT61PF loG21PF loG41PF loG41UF loG61PF loT21PF loT41PF loT41UF loT61PF/) 
{$nocc{$_} ++}
open SETS, $setsfile;
my $header = <SETS>; $header =~ s/^#//; chomp $header;

my ($i, %settingss, @setorder) = (0);
#for (split /\t/, $header) {$sets{$_} = $i; push @setsorder, $_; $i ++}
for (split /\t/, $header) {push @setorder, $_}
#my @outsetsorder = qw/Run OPTIMIZE OVERLAP_MAX_SIZE OVERLAP_MIN_SIZE PRIMER_NUM_RETURN SOLUTION TNT_USE TREATMENTS Dataset/;
my @catorder = qw/flen frag rend done back tryfrst frst trylast last tryint int pen uni tnt_test tnt_fail tnt_nohit tnt_multi time solved reject/;
my @anal = qw//;
print join("\t", @setorder, @catorder, @anal), "\n";  # Output header line
while (<SETS>) {
 #Run TREATMENT OVERLAP_MIN_SIZE OVERLAP_MAX_SIZE OPTIMIZE SOLUTION TNT_USE PRIMER_NUM_RETURN Dataset
 next if /^#/;
 my ($cc, $line, $linecc);
 chomp; my @f = split "\t";
 my (%sets);
 for (@setorder) {
  $sets{$_} = shift @f;
  $line .= "$sets{$_}\t";
  if ($_ eq 'Dataset' and $sets{$_} eq '3082') {$linecc .= "200\t"; $cc ++}
  elsif ($_ eq 'Run') {$linecc .= "$sets{$_}-200\t"}
  else {$linecc .= "$sets{$_}\t"}
 }
 my ($run, $gff, $fragsdone, $ccfragsdone, %tots, %cctots, %d) = ($sets{Run}, $sets{Dataset}, 0);
 for (`cat sums/$run`) {
  chomp; my @g = split "\t";
  warn "$run $g[0]\n" unless scalar(@g) == 21;
  my $isle = shift @g;
  my %cats;
  for (@catorder) {$cats{$_} = shift @g}
  warn "$run $isle no frag\n" unless $cats{"frag"};
  if ($cats{done}) {
   $fragsdone += $cats{frag};
   $ccfragsdone += $cats{frag} if $cc and $ccref{$isle};
   for (qw/pen uni/) {$cats{$_} *= $cats{frag}}
  }
  for (@catorder) {
     $tots{$_} += $cats{$_};
   $cctots{$_} += $cats{$_} if $cc and $ccref{$isle};
  }    
 }
 for (qw/pen uni/) {
    $tots{$_} /= $fragsdone;
  $cctots{$_} /= $ccfragsdone if $cc;
 }
 for (@catorder) {$line .= "$tots{$_}\t"}
 print "$line\n";
 if ($cc and not $nocc{$run}) {
  for (@catorder) {$linecc .= "$cctots{$_}\t"}
  print "$linecc\n";
 }
 #$i ++; last if $i == 5;
}
