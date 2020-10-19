use strict; use warnings;
die "Usage: perl $0 RunSum > out\n" unless $ARGV[0] and -f $ARGV[0];
my ($run, %ccref) = ($ARGV[0]);  # 'cc' means the subset of 200 used for exhaustive and finding cognates among 3082
for (`cat sequences/200.gff`) {next unless /ID=([^;]+)/; $ccref{$1} ++}

my ($i) = (0);
my @catorder = qw/flen frag rend done back tryfrst frst trylast last tryint int pen uni tnt_test tnt_fail tnt_nohit tnt_multi time solved reject/;
print join("\t", @catorder), "\n";  # Output header line
my ($cc) = (0);
my (%fragsdone, %tots, %d);
for (`cat $run`) {
 chomp; my @g = split "\t";
 my $fragCt = $g[2];
 warn "$run $g[0]\n" unless scalar(@g) == 21;
 my $isle = shift @g;
 if ($cc and not $ccref{$isle}) {next}
 $tots{$fragCt}{n} ++;
 my %cats;
 for (@catorder) {$cats{$_} = shift @g}
 warn "$run $isle no frag\n" unless $cats{"frag"};
 if ($cats{done}) {
  $fragsdone{$fragCt} += $cats{frag};
  for (qw/pen uni/) {$cats{$_} *= $cats{frag}}
 }
 for (@catorder) {
  $tots{$fragCt}{$_} += $cats{$_};
 }
}
for my $fragCt (sort {$a <=> $b} keys %tots) {
 for (qw/pen uni/) {
  $tots{$fragCt}{$_} /= $fragsdone{$fragCt};
 }
 my $line = $fragCt;
 for ('n', @catorder) {$line .= "\t$tots{$fragCt}{$_}"}
 print "$line\n";
}
