#! usr/bin/perl
use strict; use warnings; use File::Spec; use Time::HiRes qw(time); use POSIX qw(strftime);

my ($verbose, $noprimers, $map, $configfile, %config, %segs, %coord, $setting, $recur, $p3begin, $ret, $optP);
my (@solved, @frags, $fCt, $fLen, $seg, %seen, %ass, %p3s, $tntcall, %tnts, %tntbads, %stats, %windows);
my $t0 = time;
my $invoke = "Called \"$0 @ARGV\" on " . localtime;
my %paths = (bin => File::Spec->rel2abs($0)); $paths{bin} =~ s/\/([^\/]+)$//; my $scriptname = $1;
Options(); # see bottom of script; help and other messages, Getopt
chdir $paths{out};
open LOG, ">bigdna.log";
Log("$invoke\nConfig file: $paths{in}/$configfile\n") if $verbose;
ReadDefault();
my $set = $config{SETTING};
ReadConfig();
LoadTreats();
SetCheck();
CollectSeqs();
for (qw/REGION SEGMENT/) {$config{$_} = ()}  # Free some memory, all reqd info now in %segs
Rebuilds();
Splice();

#SUBROUTINES
sub Recurse {  # Full backdown for 1st-solution triggered by return value 1 at Solution()
 $recur ++;
 for my $i (1 .. $fCt) {
  next if $frags[$i];                           # Find first untreated fragment 
  my ($r, $s) = ($recur, scalar(@solved));
  Log("Begin recursion R$r F$i\n") if $verbose;
  my $pcrs = Primer3($i);                       # List of candidate PCRs for fragment
  for my $pcr (@{$pcrs}) {                      # Try them all
   next if $seen{$$pcr{R}[0]};                  # Skip right end already found to fail
   $frags[$i] = $pcr;                           # Test valid PCR
   if ($i == $fCt and $ret = Solution()) {      # SOLUTION! Primer3 filled last PCR (Test with Tnt)
    if ($ret) {
     Log("Start $ret drops from R$r\n") if $verbose;
     pop @frags;
     return $ret-1;
    }
   }
   $ret = Recurse();                            # RECURSE forward
   if ($ret) {                                  # Continuing to drop back
    pop @frags;
    return $ret-1;
   }
  }
  $frags[$i] = '';                              # All pcrs have failed; empty this fragment
  $s = @solved-$s;
  Log("End R$r, $s new-solved, lastR=$frags[$i-1]{R}[0]\n") if $verbose;
  unless ($s) {                                 # No new solutions came from this recursion; ie it failed
   $seen{$frags[$i-1]{R}[0]} ++;                # Mark failed right end to prevent retesting
   Log("Mark $frags[$i-1]{R}[0] failed\n") if $verbose;
  }
  return 0;                                     # Backtrack one level
 }
}

sub Stats {
 my @cats = qw/flen frag rend done back tryfrst frst trylast last tryint int pen uni tnt_test tnt_fail tnt_nohit tnt_multi time solved reject/;
 @stats{qw/tnt_test tnt_fail/} = (scalar(keys %tnts), scalar(keys %tntbads));
 for (@cats) {$stats{$_} = 0 unless $stats{$_}}
 $stats{back} = 1 if $stats{done} and $stats{tryint} + $stats{trylast} > $stats{int} + $stats{last};
 Log("$stats{tnt_test} tntblasted pcrs; $stats{tnt_fail} rejected\n") if $verbose and $$set{TNT_USE} ne 'off';
 Log(join ("\t", 'Segment', @cats) . "\n");
 my $log = $seg; for (@cats) {$log .= "\t$stats{$_}"} $log .= "\n"; Log($log);
}

sub Rebuilds {  # Rebuild each segment of assembly (including single-PCR segments)
 for (sort {$segs{$a}{order} <=> $segs{$b}{order}} keys %segs) {
  ($seg, $recur) = ($_, 0);
  my $s = $segs{$seg};
  %seen = (); %p3s = (); %tnts = (); %tntbads = (); @solved = (); 
  %windows = (Lin => 1, Lout => 1, Rin => $$s{len}, Rout => $$s{len});  # Tolerance for rebuild termini, default is forcing to extreme termini
  my (@L, @R);
  if ($$s{DELTA_INT} && $$s{annots}) {  # Find int gene closest to an end
   for (qw/L_GENE R_GENE L_TOLERANCE R_TOLERANCE/) {delete $$s{$_}}
   for (@{$$s{annots}}) {
    if ($$_[8] =~ /annot=[PRSY]/) {  # Int gene
     my ($diff, $side) = ($$s{len}-$$_[3]+1, 'R');
     ($diff, $side) = ($$_[4], 'L') if $diff > $$_[4];
     @windows{qw/diff Iend IL IR/} = ($diff, $side, $$_[3]-1, $$_[4]+1) unless $windows{Iend} and $windows{diff} < $diff;
    }
    push @L, $$_[3];
    unshift @R, $$_[4];
    if ($windows{Iend} and $windows{Iend} eq 'L' and not $windows{InL} and $$_[3] >= $windows{IR}) {$windows{InL} = $$_[3]}
   }
   @R = sort {$b <=> $a} @R;
   if    ($windows{Iend} and $windows{Iend} eq 'L') {
    @windows{qw/Lout Lin/} = @windows{qw/IR InL/};
    for (@R) {next if $_ > $$s{len} -99; @windows{qw/Rout Rin/} = ($$s{len} -99, $_); last}
   }
   elsif ($windows{Iend} and $windows{Iend} eq 'R') {
    for (@R) {next if $_ > $windows{IL}; @windows{qw/Rout Rin/} = ($windows{IL}, $_); last}
    for (@L) {next if $_ < 100; @windows{qw/Lout Lin/} = (100, $_); last}
   }
  }
  if ($$s{L_GENE} && $$s{annots}) {$windows{Lin} = $$s{annots}[0][3]}
  elsif ($$s{L_TOLERANCE}) {$windows{Lin} = $$s{L_TOLERANCE}}
  if ($$s{R_GENE} && $$s{annots}) {  # Find rightmost gene right end (R_GENE) and shortest terminal integrase gene deletion (DELTA_INT)
   for (@{$$s{annots}}) {push @R, $$_[4]};
   @R = sort {$b <=> $a} @R;
   $windows{Rin} = $R[0];
  } elsif ($$s{R_TOLERANCE}) {$windows{Rin} = $$s{len}-$$s{R_TOLERANCE}+1}
  #for (sort keys %windows) {print "$_ $windows{$_}\n"}
  $fCt = ($windows{Rout}-$windows{Lout}+1)/$$set{PCR_MAX_SIZE}; $fCt ++ unless $fCt == int $fCt; $fCt = int $fCt;
  $fLen = int(($windows{Rout}-$windows{Lout}+1)/$fCt);
  %stats = (len => $$s{len}, frag => $fCt, flen => $fLen);
  my $max = $fLen + 2* $$set{REBUILD_WINDOW};
  my $opt = $$set{opt};
  $p3begin =~ s/(PRIMER_PRODUCT_SIZE_RANGE=\d+-)\S*/$1$max/;
  $$set{PRIMER_PRODUCT_SIZE_RANGE} = '200-' . ($fLen + 2* $$set{REBUILD_WINDOW}); 
  Log("$segs{$seg}{len}-bp Segment $seg rebuilding with $fCt PCRs of ~$fLen: end windows L:$windows{Lout}-$windows{Lin}, R:$windows{Rin}-$windows{Rout}\n");
  if ($fCt > 2) {unless (RightTest()) {$stats{rend} = 1; Log(Timer() . "FATAL: RightTest failed for Segment $seg\n"); Stats(); exit}}
  @frags = (); $frags[0]{R}[0] = 0;  # $frags[0] is never used (but "preceding" right end should be defined for %seen)
  Recurse();  # Start recursion
  unless (@solved) {Log(Timer() . "FATAL: $segs{$seg}{len}-bp Segment $seg $fCt-PCR rebuilding failed\n"); Stats(); exit}
  if ($$set{SOLUTION} eq 'exhaustive') {
   if ($opt eq 'pen') {$opt = 'uni'} else {$opt = 'pen'}  # Find and report optimal by undesired criterion first
   @solved = sort {$$a{$opt} <=> $$b{$opt}} @solved;  # Post-exhaustive non-opt not tnt-checked!
   my %sco; @sco{qw/pen uni/} = ($solved[0]{pen}/$fCt, $solved[0]{uni}/$fCt);
   $stats{$opt} = $sco{$opt};
   Log("$opt-optimizing: penalty=$sco{pen}, nonuniformity=$sco{uni}:$solved[0]{sum}\n");
   $opt = $$set{opt};
   @solved = sort {$$a{$opt} <=> $$b{$opt}} @solved;  # Leave solved[0] as the best by OPTIMIZE
   if ($$set{TNT_USE} eq 'post-exhaustive') {
    Log("Tnt launched, post-exhaustive\n") if $verbose; 
    while (@solved) {
     my ($sol, @tnt) = ($solved[0]{sol});
     for (@{$sol}) {my @s = @{$_}; my $sum = "$s[0]-$s[3],$s[1],$s[4]"; push @tnt, [$s[0], $s[2], $s[3], $s[5], $sum];}
     if (Tnt(\@tnt)) {  # Misprime, requires backtracking
      $stats{reject} ++;
      Log("Rejected for mispriming: $solved[0]{sum}\n") if $verbose;
      shift @solved;
      next;
     }
     last;
    }
   }
  }
  unless (@solved) {Log(Timer() . "FATAL: $segs{$seg}{len}-bp Segment $seg $fCt-PCR rebuilding failed after post-exhaustive tntBlast\n"); Stats(); exit}
  my %sco; @sco{qw/pen uni/} = ($solved[0]{pen}/$fCt, $solved[0]{uni}/$fCt);
  for (qw/pen uni/) {$stats{$_} = $sco{$_} unless $stats{$_}}
  Log("Final $opt-optimizing: penalty=$sco{pen}, nonuniformity=$sco{uni}:$solved[0]{sum}\n");
  Log("$stats{reject} rejects of $stats{solved} solutions\n") if $verbose and $stats{reject};
  @stats{qw/time done/} = (sprintf("%.6f", time-$t0), 1);
  Stats();
  $segs{$seg}{pcr} = $solved[0]{sum};
 }
}

sub Solution {  # Return value is number of backtracks to perform immediately
 my ($sum, @uni, @pen, %t);
 $stats{solved} ++;
 for my $i (1..$#frags) {
  $t{pen} += $frags[$i]{pen};
  $t{uni} += $frags[$i]{uni};
  push @{$t{sol}}, [@{$frags[$i]{L}}[0,1,2], @{$frags[$i]{R}}[0,1,2]];
  $sum .= " F$i:$frags[$i]{L}[0]-$frags[$i]{R}[0] $frags[$i]{L}[1],$frags[$i]{R}[1];";
 }
 if ($$set{TNT_USE} eq 'per-solution') {
  Log("Tnt per-solution launched for " . scalar(@frags) . " fragments, per-solution\n") if $verbose;
  for my $i (1..$#frags) {push @{$t{tnt}}, [@{$frags[$i]{L}}[0,2], @{$frags[$i]{R}}[0,2], $frags[$i]{sum}]}  # [0] Lcoord [1] Lseq [2] Rcoord [3] Rseq [4] pcrSummary
  $ret = Tnt($t{tnt});
  if ($ret) {  # Misprime, requires backtracking
   Log("Rejected for mispriming: $sum\n") if $verbose;
   $stats{reject}++;
   return $ret;
  }
 }
 Log("$seg SOLVED: Penalty=$t{pen}; Nonuni=$t{uni};$sum\n") if $verbose;
 push @solved, {pen => $t{pen}, uni => $t{uni}, sol => $t{sol}, sum => $sum};
 if ($$set{SOLUTION} eq 'exhaustive') {return 0} else {return $#frags}  # Full backdown for 1st
}

sub Tnt {  # Tests all PCRs of a rebuild solution OR (per-recursion) all PCRs from a Primer3 run
 return if $$set{TNT_USE} eq 'off';  # XXX?
 my ($call, $pcrs, %expected, %found, $expect, @in, $hit) = ("$tntcall -d $segs{$seg}{FASTA} -l " . (2*$fLen) . " > /dev/null", @_);
 my $backtracks = 0;  # How many recursion levels to drop to correct mispriming PCR
 for (my $i=$#{$pcrs}; $i >= 0; $i--) { #  Contents of each PCR are: [0] Lend, [1] Lseq, [2] Rend, [3] Rseq, [4] ExpectLabel
  $expect = ${$pcrs}[$i][4];  # Each $expect should be a label like this: 'Lend-Rend,lenL,lenR'
  $backtracks = $#{$pcrs} -$i if $tntbads{$expect};  # Backtracks to fix leftmost bad PCR
  next if $tnts{$expect};
  $expect =~ /(,.*)/; 
  %{$expected{$expect}} = (ct => 0, lens => $1);
  push @in, join(' ', $expect, @{$$pcrs[$i]}[1,3]) . "\n";  # @in are tnt inputs; empty if all previously tested
 }
 return 1 if $backtracks and $$set{TNT_USE} eq 'post-exhaustive';  # Post-exhaustive: return 0=OK, 1=bad
 unless (@in) {  # All Pcrs have previously been tested
  return $backtracks;  # Per-solution: backtrack or not; Per-recursion: return value unimportant but saves a Tnt run
 }
 open INPUT, ">tntIn"; print INPUT join('', @in); close INPUT;
 unlink "tntOut";
 system $call;
 unless (-f "tntOut") {Log("WARNING: " . scalar(@in) . " tntblast inputs, no output for $call\n"); return}
 for (`cat tntOut`) {
  if (/^name = (\S+)/) {($expect, $hit) = ($1, '')}
  elsif (/^amplicon range = (\d+) \.\. (\d+)/) {  # Tntblast output is zero-based
   my ($L, $R) = ($1+1 -$segs{$seg}{FA_L}+1, $2+1 -$segs{$seg}{FA_L}+1);
   ($L, $R) = ($segs{$seg}{FA_R} -($2+1) + 1, $segs{$seg}{FA_R} -$1) if $segs{$seg}{FA_ORIENT} eq '-';
   for ($L, $R) {$_ = $_ % $segs{$seg}{entrylen}}
   $hit = "$L-$R" . $expected{$expect}{lens};
  }
  elsif (/^>(\S+)/) {  # $1 is the ENTRY hit
   next if $tnts{$expect};
   unless ($hit) {Log("Entry but no hit for $expect\n") if $verbose; next}
   my $entryhit = $1;
   if ($entryhit eq $segs{$seg}{ENTRY} and $expected{$hit}) {
    Log("expected $expect hits other expected $hit\n") if $verbose and $hit ne $expect;
    $expected{$hit}{ct} ++;
   } else {$tnts{$expect} =1; $tntbads{$expect} =1; Log("$expect misprimes at $hit on $entryhit\n") if $verbose}
  }
 }
 unlink "tntOut", 'tntIn';
 for (keys %expected) {
  next if $tnts{$_};
  $tnts{$_} = 1;
  if    ($expected{$_}{ct} == 0) {$stats{tnt_nohit} ++; Log("No expected PCR $_ for segment $seg found by tntblast\n") if $verbose}
  elsif ($expected{$_}{ct} > 1)  {$stats{tnt_multi} ++; $tntbads{$_} = 1; Log("Expected PCR $_ for segment $seg found multiply by tntblast\n") if $verbose}
 }
 $backtracks = 0;
 for (my $i=$#{$pcrs}; $i >= 0; $i--) { #
  $backtracks = $#{$pcrs} -$i if $tntbads{$$pcrs[$i][4]};  # Backtracks to fix leftmost bad PCR
 }
 return $backtracks;
}

sub SetCheck {
 for (sort keys %{$set}) {$p3begin .= "$_=$$set{$_}\n" if /^PRIMER_/}
 $p3begin .= "PRIMER_PRODUCT_SIZE_RANGE=200-x\n";
 return if $$set{TNT_USE} eq 'off';
 if ($$set{TNT_USE} eq 'post-exhaustive' and $$set{SOLUTION} eq 'first') {
  $$set{TNT_USE} = 'pre-solution';
  Log("TNT_USE post-exhaustive not allowed with SOLUTION=first, resetting to 'pre-solution'\n");
 }
 my %tnt_args = (qw/TNT_MIN_TM -e TNT_CLAMP --primer-clamp TNT_MAX_GAP --max-gap TNT_MAX_MISMATCH --max-mismatch/);
 $tntcall = 'tntblast -i tntIn -o tntOut';
 for (keys %tnt_args) {if (defined $$set{$_}) {$tntcall .= " $tnt_args{$_} $$set{$_}"}}
 Log("Tntblast basic call (database and max length customized for each segment): $tntcall\n") if $verbose;
}

sub Splice {
 my (@primers, $last);
 for (sort {$segs{$a}{order} <=> $segs{$b}{order}} keys %segs) {
  $seg = $_;
  my $e = $segs{$seg}{ENTRY};
  while ($segs{$seg}{pcr} =~ s/ F(\d+):(\d+)-(\d+) (\d+),(\d+);//) {
   my ($frag, $L, $R, $sL, $sR) = ($1, $2, $3, $4, $5);
   my @seqs = (substr($segs{$seg}{seq}, $2-1, $4), Revcomp(substr($segs{$seg}{seq}, $3-$5, $5)));
   if ($last and $1 == 1) {($primers[-1][0], $seqs[0]) = Frankenstein (@{$primers[-1]}[0,5], $seqs[0], $L, $sL)}
   push @primers, [$seqs[0], $seg, $frag, 'F', $L, $sL, "$2-" . ($2+$4-1)];  # [0] primerSeq [1] seg [2] frag [3] dir [4] 5'end [5] len [6] span
   push @primers, [$seqs[1], $seg, $frag, 'R', $R, $sR, "$3-" . ($3-$5+1)]; 
  }
  $last = $primers[-1][0] if @primers;
  unless ($$set{LINEAR}) {  # Circularize
   ($primers[-1][0], $primers[0][0]) = Frankenstein(@{$primers[-1]}[0,5], @{$primers[0]}[0,4,5]);
  }
 }
 open OUT, ">primers.txt"; for my $p (@primers) {print OUT "$$p[0]\t$segs{$$p[1]}{order}_$$p[2]$$p[3]\t$$p[1]:$$p[6]\n"} close OUT;
}

sub Frankenstein {
 my ($seqL, $lenL, $seqR, $startR, $lenR) = @_;
 my @franks = (lc(Revcomp($seqR)) . $seqL, lc(Revcomp($seqL)) . $seqR);
 my $diff = $$set{OVERLAP_MIN_SIZE} - $lenL - $lenR;  # Length from L segment to make up to OVERLAP_MIN
 if ($diff > 0) {$franks[0] = lc(Revcomp(substr($segs{$seg}{seq}, $startR+$lenR-1, $diff))) . $franks[0]}
 return (@franks);
}

sub Log {print LOG $_[0]} #warn $_[0] if $verbose

sub RightTest {  # Feasibility test for rightmost PCR during rebuilding; avoids slow discovery that the last PCR is impossible; no tntBlast testing
 my ($len, $in, @pcrs, @tntin) = ($windows{Rout});
 my @ok = ($len - $fLen - int($$set{REBUILD_WINDOW}/2) +2);
 $ok[1] = $ok[0] + $$set{REBUILD_WINDOW} -1;
 $ok[0] = $windows{Lout} if $ok[0] < $windows{Lout};
 push @ok, @windows{qw/Rin Rout/};
 if ($ok[2] == $ok[3]) {$in .= "SEQUENCE_FORCE_RIGHT_START=$ok[2]\n"; $ok[2] -= $$set{PRIMER_MAX_SIZE} - $optP}
 $in .= "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=" . join(',', $ok[0], $ok[1]-$ok[0]+$optP, $ok[2]-$optP+1, $len-$ok[2]+$optP) . "\n";
 Log($in) if $verbose;
 open P3, ">p3in"; print P3 "$p3begin${in}SEQUENCE_TEMPLATE=$segs{$seg}{seq}\n=\n"; close P3;
 system "primer3_core p3in > p3out";
 unless (-f 'p3out') {Log("WARNING: No Primer3 output for seg $seg RightTest\n"); unlink 'p3in'; return 0;}
 for (`cat p3out`) {
  push (@{$pcrs[$2]{$1}}, $3    ) if /^PRIMER_(R|L)[^_]+_(\d+)_SEQUENCE=(\S+)/;
  unshift (@{$pcrs[$2]{$1}}, $3, $4) if /^PRIMER_(R|L)[^_]+_(\d+)=(\d+),(\d+)/;
 }
 unlink 'p3in', 'p3out';
 for (my $j=$#pcrs; $j >= 0; $j--) {splice(@pcrs, $j, 1) if $pcrs[$j]{L}[0] > $ok[1] or $pcrs[$j]{R}[0] < $ok[2]}  # Primer3 5' control
 Log("RightTest on segment $seg: " . scalar(@pcrs) . " pairs\n") if $verbose;
 return scalar(@pcrs);
}

sub Primer3 {
 my ($i) = @_;
 my ($target, $len, $in, @pcrs, @tntin, @ok) = ($i * $fLen, $segs{$seg}{len});  # @ok holds ranges for acceptable 5' primer ends
 if ($i == 1) {      # Special L end instructions for leftmost
  @ok = @windows{qw/Lout Lin/};
  if ($ok[0] == $ok[1]) {$ok[1] += $$set{PRIMER_MAX_SIZE} - $optP; $in .= "SEQUENCE_FORCE_LEFT_START=$ok[0]\n"}
 } else {            # All other fragments of rebuilding
  @ok = (${$frags[$i-1]}{R}[0] - $$set{OVERLAP_MAX_SIZE} + 1, ${$frags[$i-1]}{R}[0] - $$set{OVERLAP_MIN_SIZE} + 1);
 }
 if ($i == $fCt) {  # Special R end instructions for rightmost
  push @ok, @windows{qw/Rin Rout/};
  if ($ok[2] == $ok[3]) {$ok[2] -= $$set{PRIMER_MAX_SIZE} - $optP; $in .= "SEQUENCE_FORCE_RIGHT_START=$ok[3]\n"}
  $ok[3] = $len;
 } else {            # All other fragments of rebuilding
  push @ok, $target - int($$set{REBUILD_WINDOW}/2), $target + int($$set{REBUILD_WINDOW}/2 - 1);
 }
 my $p3 = "$ok[0],$ok[3]";
 return $p3s{$p3} if $p3s{$p3};
 $in .= "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=" . join(',', $ok[0], $ok[1]-$ok[0]+$optP, $ok[2]-$optP+1, $ok[3]-$ok[2]+$optP) . "\n";
 Log($in) if $verbose;
 open P3, ">p3in"; print P3 "$p3begin${in}SEQUENCE_TEMPLATE=$segs{$seg}{seq}\n=\n"; close P3;
 system "primer3_core p3in > p3out";
 unless (-f 'p3out') {P3stats($i, 0, "WARNING: No Primer3 output for seg $seg $p3\n"); $p3s{$p3} = \@pcrs; return \@pcrs;}
 for (`cat p3out`) {
  $pcrs[$1]{pen} = $2 if /^PRIMER_PAIR_(\d+)_PENALTY=(\S+)/;
  push (@{$pcrs[$2]{$1}}, $3    ) if /^PRIMER_(R|L)[^_]+_(\d+)_SEQUENCE=(\S+)/;
  unshift (@{$pcrs[$2]{$1}}, $3, $4) if /^PRIMER_(R|L)[^_]+_(\d+)=(\d+),(\d+)/;
 }
 for my $j (0..$#pcrs) {
  $pcrs[$j]{uni} = abs($target - $pcrs[$j]{R}[0]);
  $pcrs[$j]{sum} = "$pcrs[$j]{L}[0]-$pcrs[$j]{R}[0],$pcrs[$j]{L}[1],$pcrs[$j]{R}[1]";
 }
 if ($$set{OPTIMIZE} eq 'uniform' and $i != $fCt) {  # Sort by closest to center of internal right window
  @pcrs = sort {$$a{uni} <=> $$b{uni}} @pcrs;
 }
 for (my $j=$#pcrs; $j >= 0; $j--) {splice(@pcrs, $j, 1) if $pcrs[$j]{L}[0] > $ok[1] or $pcrs[$j]{R}[0] < $ok[2] or $seen{$pcrs[$j]{R}[0]} or $tnts{$pcrs[$j]{sum}}}
 if ($$set{TNT_USE} eq 'per-recursion') {
  Log("Tnt launched for " . scalar(@pcrs) . " pcrs, per-recursion\n") if $verbose;
  for (@pcrs) {push @tntin, [@{$$_{L}}[0,2], @{$$_{R}}[0,2], $$_{sum}]}
  Tnt(\@tntin) if @tntin;
  for (my $j=$#pcrs; $j >= 0; $j--) {splice(@pcrs, $j, 1) if $tntbads{$pcrs[$j]{sum}}}
 }
 my $ct = scalar(@pcrs); P3stats($i, $ct, "Primer3 search on segment $seg: fragment $i left $ok[0]-$ok[1], right $ok[2]-$ok[3]: $ct pairs\n");
 unlink 'p3in', 'p3out' unless $seg eq 'Ape1.38.V';  # CHANGE IN FINAL XXX
 $p3s{$p3} = \@pcrs;
 return \@pcrs;  # Keys: pen uni L(pos,len) R(pos,len)
}

sub P3stats {
 my ($i, $ct, $log) = @_;
 if ($i == 1) {$stats{tryfrst} ++; $stats{frst} ++ if $ct}
 elsif ($i == $fCt) {$stats{trylast} ++; $stats{last} ++ if $ct}
 else {$stats{tryint} ++; $stats{int} ++ if $ct}
 Log($log) if $verbose;
}

sub ReadConfig {  # L=0 > 1; R=0 > full-length
 my ($order, $type, $name) = (1);
 for (`cat $configfile`) {
  next if /^#/;
  my ($key, $value) = ('', '');
  chomp;
  if (/^(TYPE)=(\S+)/) {  # Top record lines
   ($key, $value) = ($1, $2);
   $type = $value;
   Log("Undefined type $type\n"), exit unless $type =~ /^SETTING|REGION|SEGMENT$/;
   $setting .= "=\n" unless $type eq 'SETTING';
  } elsif (/^(NAME)=(.+)/) {
   ($key, $value, $name) = ($1, $2, $2);
   Log("No type specified in configfile $configfile line $_\n"), exit unless $type;
   Log("Type $type name $value duplicated\n"), exit if $config{$type}{$value};
   if ($type eq 'SEGMENT') {$segs{$value}{order} = $order; $order ++;}
   @{$config{$type}{$value}}{qw/FASTA ENTRY ANNOT L R ORIENT/} = ('', '', '', 1, 0, '+');  # R=0 changes to full length once known
  } elsif (/^\s*=\s*$/) {  # Record end
   ($type, $name) = ('', ''); next
  } elsif (/^([^=#]+)=(.+)/) {  # Tag, stopping early if SETTING
   ($key, $value) = ($1, $2);
   Log("No type specified in configfile $configfile line $_\n"), exit unless $type;
   if ($type eq 'SETTING') {$config{$type}{$key} = $value; next}
   Log("No name specified in configfile $configfile line $_\n"), exit unless $name;
   $config{$type}{$name}{$key} = $value;
   for (qw/FASTA ANNOT/) {  # Convert relative path (from config file dir) into absolute path; check file existence
    next unless $key eq $_;
    $value = File::Spec->rel2abs($value);
    unless (-f $value) {Log("No $_ $value\n"); exit}
   }
   if ($value =~ /[^\d]/ and ($key eq 'L' or $key eq 'R')) {Log("L and R values must be zero or positive integer\n"); exit}
   if ($key eq 'ORIENT' and $value !~ /^[\+\-]$/) {Log("ORIENT must be + or -\n"); exit}
   $config{$type}{$name}{L} = 1 unless $config{$type}{$name}{L};
  }
  $setting .= "$key=$value\n" unless $type eq 'SETTING';
 }
 for (qw/PRIMER_OPT_SIZE PCR_MAX_SIZE REBUILD_WINDOW OVERLAP_MAX_SIZE OVERLAP_MIN_SIZE OPTIMIZE/) 
 {Log("FATAL: $_ not set\n"), exit unless defined $config{SETTING}{$_}}
 $config{SETTING}{opt} = $config{SETTING}{OPTIMIZE}; $config{SETTING}{opt} =~ s/^(...).*/$1/;  # Short form often useful 
 $optP = $$set{PRIMER_OPT_SIZE};
 Log("No TYPE=SEGMENT line(s) in configfile $configfile\n"), exit unless %segs;
}

sub LoadTreats {
 my %treats = qw/hiL PRIMER_MAX_SIZE loT PRIMER_MIN_TM hiT PRIMER_MAX_TM hiH PRIMER_MAX_HAIRPIN_TH hiM PRIMER_MAX_POLY_X loG PRIMER_MIN_GC/;
 if ($$set{QUICK_P3}) {
  for my $t (split ',', $$set{QUICK_P3}) {
   $t =~ s/(\d+)$//;
   my $value = $1;
   next unless $treats{$t};  # Unrecognized abbreviation
   $$set{$treats{$t}} = $value;
   Log("Key Primer3 setting $treats{$t} ($t) = $value\n") if $verbose;
  }
 }
 if (%{$set}) {
  $setting .= "=\nTYPE=SETTING\n";
  for (sort keys %{$set}) {$setting .= "$_=$$set{$_}\n"}
 }
 Log($setting . "=\n") if $verbose;
}

sub CollectSeqs {  # Progresses to REGION if necessary to find FASTA
 for my $s (sort {$segs{$a}{order} <=> $segs{$b}{order}} keys %segs) {
  $seg = $config{SEGMENT}{$s};
  if ($$seg{REFERENCE}) {  # In case via REGION, first get REGION then SEGMENT
   my $name = $$seg{REFERENCE};
   my $reg = $config{REGION}{$name};
   unless ($$reg{seq}) {  # Get REGION's seq etc.
    my @params = @{$reg}{qw/FASTA ENTRY L R ORIENT/};
    @{$reg}{qw/seq ENTRY entrylen len source/} = GetSeq(@params);
    $$reg{annots} = GetAnnot($$reg{ANNOT}, @params) if $$reg{ANNOT};
   }
   $$seg{FASTA} = $$reg{FASTA};
   my @params =	($name, '', @{$seg}{qw/L R ORIENT/});
   @{$seg}{qw/seq ENTRY entrylen len source/} = GetSeq(@params, $$reg{seq}, $$reg{entrylen});
   $$seg{R} = $$seg{len} unless $$seg{R};
   if ($$seg{ORIENT} eq $$reg{ORIENT}) {$$seg{FA_ORIENT} = '+'} else {$$seg{FA_ORIENT} = '-'}
   if ($$reg{ORIENT} eq '+') {@{$seg}{qw/FA_L FA_R/} = ($$reg{L}+$$seg{L}-1, $$reg{L}+$$seg{R}-1)}
   else {@{$seg}{qw/FA_L FA_R/} = ($$reg{R}-$$seg{R}+1, $$reg{R}-$$seg{L}+1)}
   $$seg{annots} = GetAnnot('', '', @params, $$reg{annots}) if @{$$reg{annots}};
   $$seg{source} .= "=$$reg{source}";
  } else {  # Fasta specified directly for segment
   Log("No FASTA nor REFERENCE for segment $s\n"), exit unless $$seg{FASTA};
   my @params = @{$seg}{qw/FASTA ENTRY L R ORIENT/};
   @{$seg}{qw/seq ENTRY entrylen len source/} = GetSeq(@params);
   $$seg{R} = $$seg{len} unless $$seg{R};
   @{$seg}{qw/FA_L FA_R FA_ORIENT/} = @{$seg}{qw/L R ORIENT/};
   $$seg{annots} = GetAnnot($$seg{ANNOT}, @params) if $$seg{ANNOT};
  }
  for (keys %{$seg}) {$segs{$s}{$_} = $$seg{$_}}
 }
}

sub GetSeq {
 my ($fa, $entry, $L, $R, $ori, $full, $entrylen, $collect, $seq) = @_;
 unless ($full) {
  for (`cat $fa`) {
   if (/^>(\S+)/) {
    last if $collect;
    $entry = $1 unless $entry;
    $collect ++ if $entry eq $1;
   } elsif ($collect) {chomp; $full .= $_}
  }
 }
 Log("No seq for entry $entry in fasta $fa\n"), exit unless $full;
 my $len = length($full);
 $entrylen = $len unless $entrylen;
 $R = $len unless $R;
 while ($R > $len) {$R -= $len}
 Log("R coord $R must be from 1 to length(=$len) of entry $entry in fasta $fa\n"), exit if $R > $len or $R < 1;
 if ($L <= $R) {
  Log("R coord $R must be greater than L coord $L < length $len\n"), exit unless $R > $L;
  $seq = uc substr($full, $L-1, $R-$L+1);
 } else {  # Wrap around circle
  $seq = uc(substr($full, $L-1) . substr($full, 0, $R));
 }
 if ($ori eq '-') {$seq = Revcomp($seq)}
 return ($seq, $entry, $entrylen, length($seq), "$fa/$entry:${L}_${R}_$ori");
}

sub GetAnnot {
 my ($file, $fa, $entry, $L, $R, $ori, $ann, @annots) = @_;
 @annots = @{$ann} if $ann;
 if ($file) {
  for (`cat $file`) {
   my @f = split "\t";
   if ($f[0] eq $entry and $f[3] > $L and $f[4] < $R) {push @annots, [@f]}
  }
 }
 for my $ann (@annots) {
  if ($ori eq '-') {
   ($$ann[3], $$ann[4]) = ($R-$$ann[4]+1, $R-$$ann[3]+1);
   $$ann[6] =~ tr/\+\-/-+/;
  } else {
   ($$ann[3], $$ann[4]) = ($$ann[3]-$L+1, $$ann[4]-$L+1);
  }
 }
 @annots = reverse @annots if $ori eq '-';
 return \@annots;
}

sub Revcomp {my $seq = $_[0]; $seq =~ tr/ACGT/TGCA/; $seq = reverse $seq; return $seq}
sub ReadDefault {
 $paths{lib} = $paths{bin}; $paths{lib} =~ s/[^\/]*$/lib/;
 for (`cat $paths{lib}/default.txt`) {chomp; $config{SETTING}{$1} = $2 if /^([^#]\S*)=(.+)/}
 $config{SETTING}{PRIMER_FIRST_BASE_INDEX} = 1;  # Insist on 1-based indexing
 my $path = `which primer3_core`; chomp $path;
 Log("FATAL: No primer3_core in paths: $path\n"), exit unless $path =~ /\//;
 $path =~ s/primer3_core$/primer3_config\//;
 $config{SETTING}{PRIMER_THERMODYNAMIC_PARAMETERS_PATH} = $path;
}

sub Timer{my $t = sprintf("%.6f", time-$t0); return "[t=$t] "}

sub Options {
my $version = '1.0 (Nov 2020)';
#   '    |    '    |    '    |    '    |    '    |    '    |    '    |    '    |
my $help = <<END;
$scriptname version $version
Usage: perl $scriptname [options] CONFIG_FILE
 -map       Create a map of the fragments named -n.png or, if -n is not given, -fConstruct.png
 -out       <output directory> (default uses same directory as CONFIG_FILE)
 -noprimers Omit output primers.txt file
 -verbose   Print log messages to screen also
  Additional options: -help, -version, -authors, -license
  See README.txt for further details and CONFIG_FILE instructions.

END
 die $help if @ARGV == 0;
 use Getopt::Long;
 my $options_okay = GetOptions(
  'help'     => sub {print $help; exit},
  'version'  => sub {print "$scriptname version $version\n"; exit},
  'authors'  => sub {print "AUTHORS: Ivan Vuong, Catherine Mageeney, Kelly Williams (kpwilli\@sandia.gov)\n"; exit},
  'license'  => sub {print `cat $paths{bin}/../LICENSE`; exit},
  'verbose'  => \$verbose,
  'out=s'    => \$paths{out},
  'noprimers'=> \$noprimers,
  'map'      => \$map,
 );
 die $help if !$options_okay;
 $configfile = $ARGV[-1];
 die "No configuration file $configfile\n" unless -f $configfile;
 $paths{in} = File::Spec->rel2abs($configfile); $paths{in} =~ s/\/[^\/]*$//;
 $configfile =~ s/.*\///;
 $paths{out} = $paths{in} unless $paths{out};
 mkdir $paths{out};
 die "Output directory $paths{out} unavailable\n" unless -e $paths{out} and -d $paths{out};
}
