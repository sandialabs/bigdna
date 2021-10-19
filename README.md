# bigdna
BigDNA suggests primers for the long PCRs (SEGMENTs) needed for a big DNA assembly
It takes a configuration file as input and outputs the files primers.txt (listing suggested oligos for assembly) and bigdna.log.
bigDNA relies heavily on Primer3; tntBlast is also recommended unless you turn its use off with the TNT_OFF configuration tag.
Run "perl [path/]bigdna.pl -h" for additional options.
Sample call: $ perl bigdna.pl configfile

Installation on Unix / OS X
1. Download and unpack the source from https://github.com/sandialabs/bigdna
2. Install dependencies and place them in your PATH<br>
	Primer3 https://github.com/primer3-org<br>
	tntBLAST (the single CPU build is sufficient) https://public.lanl.gov/jgans/tntblast/tntblast_doc.html<br>
3. Test (INSTALL is the path to your bigdna folder)<br>
	Examine the folder INSTALL/testing/simple and its file INSTALL/testing/simple/config<br>
	Run: $ perl INSTALL/bin/bigdna.pl INSTALL/testing/simple/config<br>
	Check INSTALL/testing/simple for bigdna.log summarizing the run and primers.txt with suggested primer sequences<br>

Configuration file specifications:
The config file is based on the simple Boulder-IO multi-record format (one key=value pair tag per line, except for lone-equals-sign lines that separate records) with further modifications:
You'll find valid examples of config files in the "testing" directory of this package.
Each record requires a TYPE=value first line, where allowed values are SETTING, SEGMENT and REGION. At least one SEGMENT record is required that points to a FASTA file (see "REGION records" below). The FASTA value should be an **absolute** path/file to a fastA file (INSTALL is a shortcut to the home bigdna directory).
We recommend that each FASTA file reflect the entire sequence of the DNA prep template for PCR of SEGMENT, e.g., a whole bacterial genome.
In sample config files FASTA entries use INSTALL to denote the path to the bigdna folder, but any valid file path is OK

Tags allowed for the TYPE=SETTING record are:<br>
	PCR_MAX_SIZE=integer,<br>
	FRAGMENT_MAX_NUM=integer,<br>
	OVERLAP_MAX_SIZE=integer,<br>
	OVERLAP_MIN_SIZE=integer,<br>
	SOLUTION='first' or 'exhaustive',<br>
	OPTIMIZE='penalty', 'uniform', 'maxpenalty' or 'maxuniform' (optimization metric, default = penalty),<br>
	RETEST=1 (turns off the block to previously seen 5' ends during rebuilding),<br>
	LINEAR=1 (default is to circularize the final product),<br>
	TNT_USE='off', 'per-recursion', 'per-solution', or (for exhaustive mode) 'post-exhaustive',<br>
	THERMODYNAMICS=comma-separated list from hiLx,loTx,hiHx,hiMx,hiTx,loGx<br>
		where x is value for PRIMER_ tags (respectively) MAX_SIZE, MIN_TM, MAX_HAIRPIN_TH, MAX_POLY_X, MAX_TM or MIN_GC<br>
		or they start with 'PRIMER_' (valid Primer3 tags) or 'TNT_' (to configure tntBlast if on).<br>
These tags override any same-key tags of the global lib/defaults.txt file.

For multi-SEGMENT assembly, each SEGMENT record must appear in the order (and orientation) for the desired assembly.<br>
	REGION records are intended to facilitate reuse when deletions or insertions use multiple SEGMENTS, or to allow a more convenient sequence coordinate system.<br>
	REGION and SEGMENT records  must have a second line, NAME=value, where the value must be unique for that type and could be as simple as a serial number.<br>
	For REGION records, a FASTA tag is required (value is a path/file to a fastA file, either absolute or relative to the directory of the config file).<br>
	Each SEGMENT record must have either a FASTA tag or a single REFERENCE=region-name that leads to a FASTA.

Other optional keys in REGION and SEGMENT records are:<br>
 ENTRY: strongly suggested for pointing to a particular entry in the FASTA, otherwise the first FASTA entry will be used.<br>
 ANNOT: path/file to a gff3 formatted file (either absolute or relative to the directory of the config file); should only accompany a matching FASTA.<br>
 L: 1-based  left endpoint of the SEGMENT or REGION; L=0 or omitted L=value is converted to L=1.<br>
 R: 1-based right endpoint of the SEGMENT or REGION; R=0 or omitted R=value is converted to the rightmost postion of the FASTA(ENTRY) or REFERENCE.<br>
 ORIENT: orientation relative to the FASTA or reference REGION; only + or - values allowed; taken as + if omitted.<br>

Additional assembly option keys for SEGMENT records (omitting both types means the terminus of the PCR product will be forced to the extreme L or R end):
 L_TOLERANCE, R_TOLERANCE: the terminus of the PCR product is allowed within a window, whose size is fixed by the given value.<br>
 L_GENE, R_GENE:           the terminus of the PCR product is allowed within a window, whose size is determined by the distal-most gene end; requires access to a ANNOT file.<br>
 DELTA_INT: deletes the closest integrase gene to an end of the SEGMENT and 100 bp at the other end (dedicated to a particular phage engineering problem)<br>

Sample Configuration file<br>
TYPE=SETTING<br>
PRIMER_MAX_SIZE=33<br>
=<br>
TYPE=SEGMENT<br>
NAME=Eco837.22.Z<br>
FASTA=INSTALL/sequences/ecoli.fa<br>
ENTRY=NC_000913.2<br>
L=2753966<br>
R=2775998<br>
ORIENT=+<br>
L_TOLERANCE=200<br>
R_TOLERANCE=200<br>
=<br>

Output (primers.txt) from above<br>
agagactcccgctgtaacctCGCGAAGTCCGAAGAGAACT	1_1F	Eco837.22.Z:109-128<br>
CATGTGTCGACGCAACGATC	1_1R	Eco837.22.Z:6687-6668<br>
CCGCAACAGTTGGTGACTTG	1_2F	Eco837.22.Z:6640-6659<br>
CTGCGCTCATCGTTCGAAAG	1_2R	Eco837.22.Z:15783-15764<br>
CTCCAGAGGTAGGCCACGTA	1_3F	Eco837.22.Z:15744-15763<br>
agttctcttcggacttcgcgAGGTTACAGCGGGAGTCTCT	1_3R	Eco837.22.Z:21913-21894<br>

Each test X in the testing folder can be run as follows: $ perl bin/bigdna.pl testing/X/config

To reproduce all data generated for our publication,
from the bigdna directory:<br>
$ perl bin/indel.pl  # Populates 'runs' directory with insertion/deletion config files<br>
$ perl bin/jobmaker.pl jobs/runs.settings  # Populates 'runs' directory with main job set, jobs/tails files<br>
[EDIT jobs/batch: BATCH=runs, set --jobs to number of cpu available]<br>
$ bash jobs/batch   # Runs 191082 jobs<br>
$ perl bin/jobmaker.pl jobs/eighths.settings  # Populates 'runs' directory with 8 eighth-Mycoplasma jobs<br>
[EDIT jobs/batch to BATCH=eighths, to --jobs 8]<br>
$ bash jobs/batch  # Runs 8 eighths jobs<br>
The jobs/batch script requires installation of GNU parallel (https://www.gnu.org/software/parallel/)<br>
If you don't plan to repeat this validation, you can save space by deleting the 'sequences/genomes' directory
