use strict; use warnings;

mkdir "runs";
mkdir "runs/deletes";
mkdir "runs/inserts";
open JOBS, ">jobs/indels.jobs";
for (`cat sequences/3082.gff`) {
 my @f = split "\t";
 next unless m|ID=([^;]+).*fasta=../../../([^;]+)|;
 my ($acc, $L, $R, $id, $fasta, $annfile, $len) = (@f[0,3,4], $1, $2, $2, $f[4]-$f[3]+1);
 $annfile =~ s/genome.fa/isleannots.gff/;
 #die "($acc, $L, $R, $id, $fasta, $annfile)\n";
 my @anns;
 for (`cat $annfile`) {
  my @g = split "\t";
  push @anns, $_ unless $g[3]>$f[4] or $g[4]<$f[3];
 }
 my @g = split "\t", $anns[int($#anns/2)];
 #warn "$id\t$L-", $g[3]-1, ",", $g[4]+1, "-$R\n";
 #die scalar(@anns), " $g[3] $g[4]\n"
 my $segL = join("\n", 
 "TYPE=SETTING",
 "OVERLAP_MIN_SIZE=40",
 "OVERLAP_MAX_SIZE=80",
 "PRIMER_NUM_RETURN=7",
 "TNT_USE=per-recursion",
 "SOLUTION=first",
 "QUICK_P3=hiL33",
 "=",
 "TYPE=SEGMENT",
 "NAME=$id/L",
 "FASTA=../../../$fasta",
 "ENTRY=$acc",
 "L=$L",
 "R=" . ($g[3]-1),
 "R_GENE=1",
 "ANNOT=../../../$annfile",
 "=");
 my $ins = join("\n",
 "TYPE=SEGMENT",
 "NAME=gfp",
 "FASTA=../../../sequences/gfp.fa",
 "L_TOLERANCE=17",
 "R_TOLERANCE=14",
 "=");
 my $segR = join("\n",
 "TYPE=SEGMENT",
 "NAME=$id/R",
 "FASTA=../../../$fasta",
 "ENTRY=$acc",
 "L=" . ($g[4]+1),
 "R=$R",
 "L_GENE=1",
 "ANNOT=../../../$annfile",
 "=");
 mkdir "runs/deletes/$id";
 open OUT, ">runs/deletes/$id/config"; print OUT "$segL\n$segR"; close OUT;
 mkdir "runs/inserts/$id";
 open OUT, ">runs/inserts/$id/config"; print OUT "$segL\n$ins\n$segR"; close OUT; 
 print JOBS "perl bin/bigdna.pl runs/deletes/$id/config &> runs/deletes/$id/err\nperl bin/bigdna.pl runs/inserts/$id/config &> runs/inserts/$id/err\n";
 #die "runs/inserts/$id/config\n";
}
open TAILS, ">jobs/indels.tails";
print TAILS "grep -Pv '^(F|S|[0-9])' runs/inserts/*bigdna.log > sums/inserts\ngrep -Pv '^(F|S|[0-9])' runs/deletes/*bigdna.log > sums/deletes\n";
