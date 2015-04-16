#!/usr/bin/env perl
# EPN, Fri Mar 13 14:17:55 2015
# dnaorg_fetch_dna_wrapper.pl: fetch sequences from GenBank
# 
# This script can be used to fetch CDS sequences of proteins
# that link to a specific Gene database symbol, or to fetch
# sequences for which the accessions are listed in a file.
# In the latter case the accessions can either be protein,
# in which case the CDSs will be fetched, or DNA (requiring
# the -dna option).
# 
# usage: perl dnaorg_fetch_dna_wrapper.pl <symbol>
#   OR 
# usage: perl dnaorg_fetch_dna_wrapper.pl -d <symbol> -alist <file with list of protein accessions>
#   OR
# usage: perl dnaorg_fetch_dna_wrapper.pl -d <symbol> -dna -alist <file with list of nucleotide accessions>
# 
# In the first case (if -alist is not used) <symbol> should be from the Gene database.
# Otherwise (if -alist is used) it can be anything and will serve only as the name
# of the output directroy.
#
# Creates a directory called <symbol> and populates it with several output
# files, including: 
#  <symbol>.log: list of all output files and brief descriptions
#  <symbol>.cmd: list of all commands run by this script
#  <symbol>.sum: copy of all stdout from this script
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

my ($seconds, $microseconds) = gettimeofday();
my $start_secs = ($seconds + ($microseconds / 1000000.));
my $executable = $0;

# initialize variables that can be changed with cmdline options
my $do_force   = 0;     # set to '1' with -f, overwrite output dir if it exists
my $be_verbose = 0;     # set to '1' with -v, output commands as they're run to stdout
my $out_dir    = undef; # set to a value <s> with -d <s>
my $do_old     = 0;     # set to '1' with -old
# different 'modes', if all are false we run in default mode
my $do_num_proteins_mode = 0;     # set to '1' if -np option used.
my $do_acclist_mode      = 0;     # set to '1' if -alist option used.
my $do_uniprot_xref      = 0;     # set to '1' if -up option used.
my $do_dna               = 0;     # set to '1' if -dna option is used
my $acclist_file         = undef; # if -alist is used, this is defined as $symbol

&GetOptions( "f"       => \$do_force, 
             "v"       => \$be_verbose,
             "d=s"     => \$out_dir,
             "np"      => \$do_num_proteins_mode,
             "alist"   => \$do_acclist_mode,
             "up"      => \$do_uniprot_xref,
             "dna"     => \$do_dna,
             "old"     => \$do_old);

my $usage;
$usage  = "dnaorg_fetch_dna_wrapper.pl [OPTIONS] <symbol (from GENE database, e.g. 'infB'>\n";
$usage .= "dnaorg_fetch_dna_wrapper.pl [OPTIONS]      -d <symbol> -alist <file with list of protein accessions>\n";
$usage .= "dnaorg_fetch_dna_wrapper.pl [OPTIONS] -dna -d <symbol> -alist <file with list of nucleotide accessions>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-f         : force; if dir <symbol> exists, overwrite it\n";
$usage .= "\t\t-v         : be verbose; output commands to stdout as they're run\n";
$usage .= "\t\t-d <s>     : define output directory as <s>, not <symbol>\n";
$usage .= "\t\t-np        : determine number of matching protein accessions then exit\n";
$usage .= "\t\t-alist     : <symbol> is really a list of accessions, requires -d option too\n";
$usage .= "\t\t-up        : additional run experimental code for fetching non-CDS UniProt CDS via xrefs\n";
$usage .= "\t\t-dna       : -alist contains nucleotide accessions, not protein (requires -alist and -d)\n";
$usage .= "\t\t-old       : use extract_fasta_multi_exon instead of esl-fetch-cds.pl\n";
$usage .= "\n";

if(scalar(@ARGV) != 1) { die $usage; }
my ($symbol) = (@ARGV);

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if($do_force) { 
  $opts_used_short .= "-f ";
  $opts_used_long  .= "# option:  forcing overwrite of $symbol dir/file [-f]\n"; 
}
if($be_verbose) {
  $opts_used_short .= "-v ";
  $opts_used_long  .= "# option:  verbose mode [-v]\n"; 
}
if(defined $out_dir) { 
  $opts_used_short .= "-d $out_dir ";
  $opts_used_long  .= "# option:  output directory specified as $out_dir [-d]\n"; 
}
if($do_num_proteins_mode) { 
  $opts_used_short .= "-np ";
  $opts_used_long  .= "# option:  determining number of matching protein accessions, then exiting\n"; 
}
if($do_acclist_mode) { 
  $opts_used_short .= "-alist ";
  $opts_used_long  .= "# option:  $symbol is a list of accessions, not a symbol\n"; 
}
if($do_uniprot_xref) { 
  $opts_used_short .= "-up ";
  $opts_used_long  .= "# option:  trying to fetch non-CDS UniProt using xrefs\n";
}
if($do_dna) { 
  $opts_used_short .= "-dna ";
  $opts_used_long  .= "# option:  -alist contains nucleotide accessions, not protein\n";
}
if($do_old) { 
  $opts_used_short .= "-old ";
  $opts_used_long  .= "# option:  using old strategy (extract_fasta_multi_exon) [-old]\n"; 
}

# check for incompatible option combinations:
if($do_acclist_mode) { 
  if(! defined $out_dir) { 
    die "ERROR the -alist option requires the -d option";
  }
  if($do_num_proteins_mode) { 
    die "ERROR the -alist option is incompatible with -np";
  }
}
if($do_dna) { 
  if(! $do_acclist_mode) { 
    die "ERROR the -dna option requires the -alist option";
  }
  if($do_uniprot_xref) { 
    die "ERROR the -dna option is incompatible with the -up option";
  }
}

my $exec_dir  = "/panfs/pan1/dnaorg/programs";
my $idstat    = "/netopt/genbank/subtool/bin/idstat";
my $fetch_cds = "$exec_dir/esl-fetch-cds.pl";
my $idfetch   = "/netopt/ncbi_tools64/bin/idfetch";
my $query;              # an argument to a '-query' option in a edirect tool cmdline
my $cmd;                # a command to run with system() in RunCommand()
my $desc;               # description of a step
my $nlines;             # number of lines in a file
my $nlost;              # number of lost elements (often accessions) in a step
my $ncreated;           # number of created elements (often accessions) in a step, should always be 0
my $descwidth = 80;     # width of description column
my $nsecs;              # number of seconds a command took
my $errmsg;             # error message

###############
# Preliminaries
###############
# check if our output dir $symbol exists
if(! defined $out_dir) { 
  $out_dir = $symbol . "/";
}
else { 
  if($out_dir !~ m/\/$/) { $out_dir .= "/"; } # add '/'
}
if(-d $out_dir) { 
  if($do_force) { RunCommand("rm -rf $out_dir", $be_verbose, undef); }
  else          { die "ERROR directory named $out_dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $out_dir) { 
  if($do_force) { RunCommand("rm $out_dir", $be_verbose, undef); }
  else          { die "ERROR a file named $out_dir already exists. Remove it, or use -f to overwrite it."; }
}

if($do_acclist_mode) { 
  $acclist_file = $symbol;
  $symbol       = $out_dir; # this must be defined, we died above it if wasn't
  $symbol       =~ s/\///;
  if(! -e $acclist_file) { die "ERROR no file $acclist_file exists"; }
  if(! -s $acclist_file) { die "ERROR file $acclist_file is empty"; }
}


# create the dir
RunCommand("mkdir $out_dir", $be_verbose, undef);

# open the log and summary files:
my $log_file = $out_dir . $symbol . ".log";
my $sum_file = $out_dir . $symbol . ".sum";
my $cmd_file = $out_dir . $symbol . ".cmd";
my $log_FH; # file handle for log file output
my $sum_FH; # file handle for summary file output
my $cmd_FH; # file handle for command file output
open($log_FH, ">" . $log_file) || die "ERROR unable to open $log_file for writing";
open($sum_FH, ">" . $sum_file) || die "ERROR unable to open $sum_file for writing";
open($cmd_FH, ">" . $cmd_file) || die "ERROR unable to open $cmd_file for writing";
OutputFileInfo($log_file, "log file containing list and description of all output files", "", $log_FH); 
OutputFileInfo($sum_file, "summary file containing all standard output", "", $log_FH); 
OutputFileInfo($cmd_file, "command file all executed commands", "", $log_FH); 

# output banner
my $script_name = "dnaorg_fetch_dna_wrapper.pl";
my $script_desc = "Fetch DNA sequences from GenBank";
PrintToStdoutAndFile("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n", $sum_FH);
PrintToStdoutAndFile("# $script_name: $script_desc\n", $sum_FH);
PrintToStdoutAndFile("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n", $sum_FH);
PrintToStdoutAndFile("# command: $executable $opts_used_short$symbol\n", $sum_FH);
PrintToStdoutAndFile(sprintf("# date:    %s\n", scalar localtime()), $sum_FH);
if($opts_used_long ne "") { 
  PrintToStdoutAndFile("$opts_used_long", $sum_FH);
}
if(! $do_acclist_mode) { 
  PrintToStdoutAndFile("# symbol:  $symbol\n", $sum_FH);
}
PrintToStdoutAndFile("#\n", $sum_FH);


#####################################################################
# Stage 1: Fetch the CDS sequences (or nucleotide sequences if -dna).
#####################################################################
if($do_num_proteins_mode) { 
  PrintToStdoutAndFile("# Special mode (-np): determining number of matching protein accessions, then exiting.\n", $sum_FH);
}
elsif($do_dna) { 
  PrintToStdoutAndFile("# Stage 1: preparing and fetching nucleotide sequences\n", $sum_FH);
}
else { 
  PrintToStdoutAndFile("# Stage 1: preparing and fetching CDS sequences\n", $sum_FH);
}
PrintToStdoutAndFile("#\n", $sum_FH);
# create column headers
my $desc_w       = length($symbol) + 63;
my $desc_dashes = "#";
for(my $i = 0; $i < $desc_w-1; $i++) { $desc_dashes .= "-"; }
PrintToStdoutAndFile(sprintf("#%-*s  %10s  %10s  %10s  %10s  %s\n", $desc_w-1, " description", "\# output", "\# lost", "\# created", "seconds", "output-file-name"), $sum_FH);
PrintToStdoutAndFile(sprintf("%s  %10s  %10s  %10s  %10s  %s\n", $desc_dashes, "----------", "----------", "----------", "----------", "-----------------"), $sum_FH);


###################################################################################
# Step 1.1: 
# if $do_alist: 
#     if $do_dna: Fetch all nucleotide accessions listed in the accession file
#     else:       Fetch all protein    accessions listed in the accession file
# else:           Fetch all protein    accessions that match a query to our symbol 
###################################################################################
# TODO: experiment with different queries to pull in aliases/synonyms,
# possibly with command line options. Also experiment with iterating 
# this step to pull in more accessions.
my $allacc_file     = $out_dir . $symbol . ".all.acc";
my $allacconly_file = $out_dir . $symbol . ".all.acconly";

$query       = "\"$symbol [GENE]\"";
if($do_acclist_mode) { 
  $cmd         = "cat $acclist_file | sort > $allacc_file";
  $desc        = "Sorted accessions from file $acclist_file";
}
else { # default (normal) case
  $cmd         = "esearch -db protein -query $query | efetch -format acc | sort > $allacc_file";
  $desc        = "Protein accessions fetched from protein database";
}
$nsecs       = RunCommand($cmd, $be_verbose, $cmd_FH);
$nlines      = GetNumLinesInFile($allacc_file);
OutputFileInfo($allacc_file, $desc, $cmd, $log_FH);

$cmd    = "cat $allacc_file | sed 's/\\.[0-9]*//' > $allacconly_file";
$nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);

$nlost       = 0;
$ncreated    = $nlines;
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($allacc_file), $nlost, $ncreated, $nsecs, $allacc_file), $sum_FH);
if($ncreated == 0) { # this is an okay result
  PrintToStdoutAndFile("#\n", $sum_FH);
  PrintToStdoutAndFile("# No accessions fetched. Exiting.\n", $sum_FH);
  Conclude($start_secs, $sum_FH, $log_FH);
  exit 0;
}
# End of section TODO refers to

##############################################################
# Step 1.2: Remove suppressed accessions (according to idstat) 
##############################################################
my $acc_file         = $out_dir . $symbol . ".acc";
my $acc_file_lost    = $acc_file . ".lost";
my $acc_file_created = $acc_file . ".created";
if($do_dna) { $desc  = "Non-suppressed nucleotide accessions fetched from nuccore database"; }
else        { $desc  = "Non-suppressed protein accessions fetched from protein database"; }
my $idstat_file = $allacconly_file . ".idstat";
$cmd = "$idstat -i PUBSEQ_OS -A $allacconly_file > $idstat_file";
$nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);

# parse the idstat input to determine which accessions are not suppressed
removeGivenIdstat($allacc_file, $idstat_file, $acc_file);

# which accessions were lost because they're dead or withdrawn or suppressed?
$cmd = "awk '{ print \$1 }' $acc_file | comm -2 -3 $allacc_file - > $acc_file_lost";
RunCommand($cmd, $be_verbose, $cmd_FH);

# which accessions were created relative to the allacc file? (should be none)
$cmd = "awk '{ print \$1 }' $acc_file | comm -2 -3 - $allacc_file > $acc_file_created";
RunCommand($cmd, $be_verbose, $cmd_FH);

# enforce nothing was lost or created
($nlost, $ncreated, $errmsg) = CheckLostAndCreated($acc_file_lost, -1, $acc_file_created, 0); 
# in above subroutine call: -1 says it's okay if we lost some sequences, '0' says it's not okay if we created some sequences

OutputFileInfo($acc_file_lost,    "Protein accessions ($nlost) that are dead|withdrawn|suppressed according to idstat.", $cmd, $log_FH);
OutputFileInfo($acc_file_created, "Protein accessions ($ncreated) added by idstat processing (in error).", $cmd, $log_FH);

# output information on this step to stdout and sum file
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($acc_file), $nlost, $ncreated, $nsecs, $acc_file), $sum_FH);
if($errmsg ne "") { die $errmsg; }

if($do_num_proteins_mode) { 
  PrintToStdoutAndFile("#\n", $sum_FH);
  PrintToStdoutAndFile("Number-of-proteins: $ncreated\n", $sum_FH);
  Conclude($start_secs, $sum_FH, $log_FH);
  exit 0;
}
if(GetNumLinesInFile($acc_file) == 0) { # this is an okay result
  PrintToStdoutAndFile("#\n", $sum_FH);
  PrintToStdoutAndFile("# No (non-suppressed) accessions fetched. Exiting.\n", $sum_FH);
  Conclude($start_secs, $sum_FH, $log_FH);
  exit 0;
}

###################################################################################
# Step 1.3: Determine the CDS sequences/coordinates that correspond to our proteins 
# from step 1.1 and fetch them.
###################################################################################
my $efa_file         = $out_dir . $symbol . ".efa";
my $efa_file_lost    = $efa_file . ".lost";
my $efa_file_created = $efa_file . ".created";
if(! $do_dna) { 
  $cmd = "cat $acc_file | epost -db protein -format acc | efetch -format gpc | xtract -insd CDS coded_by codon_start | grep . | sort ";
  if($do_old) { 
    $cmd .= " | $exec_dir/coded_by2extract_fasta_multi_exon.pl > $efa_file";
  }
  else { 
    $cmd .= " > $efa_file";
  }
  $desc = "Protein accessions that have CDS annotation";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($efa_file, $desc, $cmd, $log_FH);

# which accessions were lost relative to the previous step? (we can lose some, these should have no CDS)
  if($do_old) { 
    $cmd = "awk '{ print \$NF }' $efa_file | comm -2 -3 $acc_file - > $efa_file_lost";
  }
  else {
    $cmd = "awk '{ print \$1 }' $efa_file | comm -2 -3 $acc_file - > $efa_file_lost";
  }
  RunCommand($cmd, $be_verbose, $cmd_FH);

# which accessions were created relative to the previous step? (should be none)
  if($do_old) { 
    $cmd = "awk '{ print \$NF }' $efa_file | comm -2 -3 - $acc_file > $efa_file_created";
  }
  else { 
    $cmd = "awk '{ print \$1 }' $efa_file | comm -2 -3 - $acc_file > $efa_file_created";
  }
  RunCommand($cmd, $be_verbose, $cmd_FH);

# enforce nothing was lost or created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($efa_file_lost, -1, $efa_file_created, 0); 
# in above subroutine call: -1 says it's okay if we lost some sequences, '0' says it's not okay if we created some sequences

  OutputFileInfo($efa_file_created, "New protein accessions ($ncreated) returned by CDS retrieval (in error).", $cmd, $log_FH);
  OutputFileInfo($efa_file_lost, "Protein accessions ($nlost) without CDS annotation.", $cmd, $log_FH);

# output information on this step to stdout and sum file
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($efa_file), $nlost, $ncreated, $nsecs, $efa_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }
} # end of 'if(! $do_dna)'

#############################################################
#EPN, Wed Mar 25 14:22:09 2015
# Experimental code for fetching UniProt CDS sequences, currently only run if -up
if($do_uniprot_xref) { # note: if this is true, then $do_dna must be false, we enforce it above
################################
# Try to fetch any uniprot sequences that we can 
  ($seconds, $microseconds) = gettimeofday();
  my $tmp_start_secs = ($seconds + ($microseconds / 1000000.));
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  my $up_root  = $out_dir . $symbol . ".up.out";
  my $tmpfile0 = $up_root . ".0";
  my $tmpfile1 = $up_root . ".1";
  my $tmpfile2 = $up_root . ".2";
  my $tmpfile3 = $up_root . ".3";
  my $tmpfile4 = $up_root . ".4";
  # get the protein lengths
  $cmd = "cat $efa_file_lost | grep -v ^WP\_ | grep -v ^YP\_ | epost -db protein -format acc | efetch -format gpc | xtract -insd INSDSeq_length | grep . | sort > $tmpfile0";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  # parse the lengths
  open(IN, $tmpfile0) || die "ERROR unable to open $tmpfile0";
  my %prot_len_H = (); # lengths of each protein, value: protein accession.version, value: length in AA
  while(my $line = <IN>) { 
    chomp $line;
    my ($prot, $len) = split(/\s+/, $line);
    if($len !~ m/^\d+$/) { die "ERROR unrecognized protein length information: $line\n"; }
    $prot_len_H{$prot} = $len;
    # print "$prot $len\n";
  }
  close(IN);
  
  # get the xrefs lines
  $cmd = "cat $efa_file_lost | grep -v ^WP\_ | grep -v ^YP\_ | epost -db protein -format acc | efetch -format gpc | xtract -insd INSDSeq_source-db | grep xrefs | sort > $tmpfile1";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  # parse the information we gathered on xrefs for the possible aliens in $efa_file_lost_maybe_alien_xrefs
  # Based on Alejandro's /panfs/pan1.be-md.ncbi.nlm.nih.gov/dnaorg/2015.03.19/get_xrefs_row/extract_xref_seqs.pl
  # we store the xrefs in a hash that maps the putative nucleotide accessions to their mother protein
  # accessions
  open(IN,  "<$tmpfile1") || die "ERROR unable to open $tmpfile1 for reading";
  open(OUT, ">$tmpfile2") || die "ERROR unable to open $tmpfile2 for writing";
  my $line;
  my %xref2prot_HA = (); # hash of arrays, key: xref accession, value: array of protein accessions that have this xref accession

  while(defined($line = <IN>)) { 
  # O02771.1	UniProtKB: locus SMN_CANFA, accession O02771; class: standard. created: Jul 15, 1998. sequence updated: Jul 1, 1997. annotation updated: Mar 4, 2015. xrefs: U50746.1, AAB58318.1, NP_001003226.1 xrefs (non-sequence databases): UniGene:Cfa.191, ProteinModelPortal:O02771, SMR:O02771, STRING:9615.ENSCAFP00000031438, PaxDb:O02771, PRIDE:O02771, GeneID:403896, KEGG:cfa:403896, CTD:39844, eggNOG:NOG296671, HOGENOM:HOG000232199, HOVERGEN:HBG000211, InParanoid:O02771, KO:K13129, Reactome:REACT_234653, NextBio:20817388, Proteomes:UP000002254, GO:0015030, GO:0005737, GO:0005829, GO:0097504, GO:0005654, GO:0005634, GO:0032797, GO:0034719, GO:0030018, GO:0003723, GO:0007409, GO:0007019, GO:0000387, InterPro:IPR010304, InterPro:IPR002999, Pfam:PF06003, SMART:SM00333, PROSITE:PS50304
    chomp $line;
    my $paccn = $line;
    $paccn =~ s/\s+.*$//; # remove everything after first space (giving O02771.1)
    my $xrefs = $line;
    my ($xref_info) = ($line=~m/xrefs:(.*?)xrefs/);
    my @id_A = split /\,/, $xref_info;
    foreach my $id (@id_A) { 
      my ($no_space_id) = ($id=~m/([A-Z0-9_\.]+)/);
      # we need to ensure we can map $no_space_id to exactly one protein accession
      if(! exists $xref2prot_HA{$no_space_id}) { 
        @{$xref2prot_HA{$no_space_id}} = (); # initialize
      }
      push(@{$xref2prot_HA{$no_space_id}}, $paccn);
      print OUT $paccn . " " . $no_space_id . "\n";
    }
  }
  close(IN);
  close(OUT);  
  
  # $tmpfile2 has a list of possible nucleotide accessions corresponding to the uniprot accessions
  # note this step uses protein accessions
  $cmd = "awk '{ print \$1 }' $tmpfile2 | epost -db protein -format acc | efetch -format gpc | xtract -insd Protein gene gene_synonym | grep . | sort > $tmpfile3";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  
  # read $tmpfile3 to get the gene names for each protein accession:
  open(IN, $tmpfile3) || die "ERROR unable to open $tmpfile3";
  my %prot2gene_HA = (); # value is a protein accession.version, value is an array of 'gene' and 'gene_synonym's for that protein
  while(my $line = <IN>) { 
    chomp $line;
    my (@elA) = split(/\s+/, $line);
    if(scalar(@elA) >= 2) { 
      my $prot = shift @elA; # first token is protein accession
      if(exists ($prot2gene_HA{$prot})) { die "ERROR read protein $prot twice in $tmpfile3"; }
      $prot2gene_HA{$prot} = [@elA]; # remaining tokens are 'gene's or 'gene_synonym's
    }
    else { # no 'gene's or 'gene synonym's for this protein, currently we do nothing, just skip it
      ;
    }
  }
  close(IN);
  
  # now for each potential nucleotide accession for protein X in $tmpfile2, try to find a CDS for the 
  # gene that protein X codes for
  $cmd = "awk '{ print \$2 }' $tmpfile2 | epost -db nuccore -format acc | efetch -format gpc | xtract -pattern INSDSeq -ACCN INSDSeq_accession-version -group INSDFeature -match INSDFeature_key:CDS -pfx \"\\n\" -element \"\&ACCN\" -block INSDFeature -match INSDFeature_key:CDS -element INSDFeature_location -subset INSDQualifier -match INSDQualifier_name:gene -element INSDQualifier_value | grep . | cat - > $tmpfile4";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  
  # now parse $tmpfile4
  # example lines:
  # CU329670.1	complement(<1..5662)	tlh1	
  # CU329670.1	complement(5726..6331)	
  # CU329670.1	complement(join(822332..822429,822489..822849))	smn1	
  my $up_efa_file = $up_root . ".uniprot-matches.efa";
  # EPN, Wed Apr  8 08:53:28 2015: decided to output all matches to $up_efa_file, we'll test them all with esl-test-cds-translate-vs-fetch.pl 
  # we may want to undo this in the future though
  #my $tmpfile6 = $up_root . ".notby3.6";
  #my $tmpfile7 = $up_root . ".lendiffer.7";
  open(IN,   "<$tmpfile4") || die "ERROR unable to open $tmpfile4 for reading";
  open(OUT5, ">$up_efa_file") || die "ERROR unable to open $up_efa_file for writing";
  #open(OUT6, ">$tmpfile6") || die "ERROR unable to open $tmpfile6 for writing";
  #open(OUT7, ">$tmpfile7") || die "ERROR unable to open $tmpfile7 for writing";
  my %coords_for_prot_HA = (); # coord strings for we've found for each protein, key: protein accession, value: array of matching CDS coord strings
  while($line = <IN>) { 
    chomp $line;
    my ($nt, $coords, $gene) = split(/\s+/, $line);
    # first we have to doctor the coords string
    # examples: 
    # nt:U50746.1   coords:4..867 has to become U50746.1:4..867
    # nt:U43883.1	  join(U43876.1:608..688,U43877.1:104..175,U43878.1:118..237,U43879.1:84..284,U43880.1:69..221,U43881.1:103..198,U43882.1:53..163,209..259)
    # has to become join(U43876.1:608..688,U43877.1:104..175,U43878.1:118..237,U43879.1:84..284,U43880.1:69..221,U43881.1:103..198,U43882.1:53..163,U43883.1:209..259)
    my $head = "";
    my $tail = "";
    if($coords =~ /(^.+\()/) { $head = $1; } # save everything up to final '(', if anything
    if($coords =~ /(\).*$)/) { $tail = $1; } # save everything starting at first ')', if anything
    $coords =~ s/^.+\(//; # remove head
    $coords =~ s/\).*$//; # remove tail
    my $new_coords = $head;
    foreach my $tok (split(',', $coords)) { 
      if($tok !~ m/^\S+\:/) { # no accession.version in this token, add it
        $tok = $nt . ":" . $tok; 
      }
      $new_coords .= $tok . ",";
    }
    # we added one too many ',', remove the final one
    $new_coords =~ s/\,$//;
    $new_coords .= $tail;
    $coords = $new_coords;
    if(defined $gene) { # this indicates we found INSDQualifier_name 'gene' for this NT accession
      foreach my $prot (@{$xref2prot_HA{$nt}}) { 
        if(exists $prot2gene_HA{$prot}) {
          for(my $g = 0; $g < scalar(@{$prot2gene_HA{$prot}}); $g++) { 
            if($gene eq $prot2gene_HA{$prot}[$g]) { # match
              my $already_printed = 0;
              # determine CDS length from coords string
              my $cds_len  = LengthFromCoords($coords);
              if(! exists $prot_len_H{$prot}) { 
                die "ERROR protein $prot somehow created in error";
              }
              my $prot_len = $prot_len_H{$prot};
              if(! exists $coords_for_prot_HA{$prot}) { 
                # initialize
                @{$coords_for_prot_HA{$prot}} = (); 
              }
              else { # search for a match in the existing coords strings for this protein
                for(my $z = 0; $z < scalar(@{$coords_for_prot_HA{$prot}}); $z++) { 
                  if($coords_for_prot_HA{$prot}[$z] eq $coords) { 
                    $already_printed = 1;
                  }
                }
              }
              if(! $already_printed) { 
                # currently, output all matches to OUT5
                print OUT5 ("$prot\t$coords\n");
                # THIS BLOCK ALTERNATIVELY OUTPUTS TO 3 DIFFERENT FILES
                ## check that the lengths make sense
                #if($cds_len % 3 == 0) { # should be...
                #  if((sprintf("%d", ($cds_len / 3)) eq $prot_len) || # matches exactly
                #    (sprintf("%d", ($cds_len / 3)) eq ($prot_len+1))) { # CDS presumably includes stop codon
                #    print OUT5 ("$prot\t$coords\t$prot_len\t$cds_len\n");
                #  }
                #  else { 
                #    print OUT7 ("$prot\t$prot_len\t$nt\t$cds_len\t$coords\n");
                #  }
                #}
                #else { # cds_len a factor of 3...
                #  print OUT6 ("$prot\t$coords\t$prot_len\t$cds_len\n");
                #}
                # END OF BLOCK
                push(@{$coords_for_prot_HA{$prot}}, $coords);
              }
            }
          }
        }
      }
    }
  }
  close(OUT);
  close(IN);
  ($seconds, $microseconds) = gettimeofday();
  my $tmp_end_secs = ($seconds + ($microseconds / 1000000.));
  
  # determine how many proteins we output >= 1 CDS for in $up_efa_file, and how many we didn't
  my $nxref = GetNumLinesInFile($tmpfile1);
  $nlost = $nxref - scalar(keys %coords_for_prot_HA);

  $desc = "Following xrefs to find CDS for accessions with direct CDS links";
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10s  %10.1f  %s\n", $desc_w, $desc, $nxref, $nlost, "?", $tmp_end_secs - $tmp_start_secs, $up_efa_file), $sum_FH);

} # end of 'if($do_uniprot_xref)'
# End of experimental code block for fetching UniProt CDS sequences, currently NOT run
################################

####################################################################################
# Step 1.4: Convert the extract_fasta_multi_exon input file to idfetch input format.
# We need to do this because idfetch can't read from a pipe, if it could we 
# wouldn't need to create this file.
####################################################################################
my $idfetch_file = $out_dir . $symbol . ".idfetch";
if($do_old) { 
  $cmd  = "awk '{ print \$1 }' $efa_file | sort | uniq > $idfetch_file";
  $desc = "Unique CDS source sequences (idfetch input file)";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($idfetch_file, $desc, $cmd, $log_FH); 
# don't output a summary of this step 
}

########################################################
# Step 1.5: Extract the sequences using esl-fetch-cds.pl
########################################################
my $fa_file         = $out_dir . $symbol . ".fa";
my $fa_file_names   = $fa_file . ".names";
my $fa_file_lost    = $fa_file . ".lost";
my $fa_file_created = $fa_file . ".created";
my $cds_or_nt       = ($do_dna) ? "nucleotide" : "CDS";
$desc  = "CDS sequences in FASTA format";

# Preliminary step: we need a version of $acc_file that has version information removed
my $acconly_file = $out_dir . $symbol . ".acconly"; # same as $acc_file but with version info removed
$cmd  = "cat $acc_file | sed 's/\\.[0-9]*//' > $acconly_file";
RunCommand($cmd, $be_verbose, $cmd_FH);

if($do_old) { 
  $cmd = "$idfetch -t 5 -c 1 -G $idfetch_file | $exec_dir/id_fasta.pl | $exec_dir/extract_fasta_multi_exon $efa_file > $fa_file";
}
elsif($do_dna) { 
  $cmd  = "$idfetch -t 5 -c 1 -G $acc_file | $exec_dir/id_fasta.pl > $fa_file";
  $desc = "nucleotide sequences in FASTA format";
}
else { # normal case
  $cmd = "$fetch_cds -odir $out_dir $efa_file > $fa_file";
}
$nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($fa_file, $desc, $cmd, $log_FH);

# get list of protein accessions we have CDS sequences (or just nucleotide sequences if $do_dna) for in fasta file
# we'll use this to make sure we fetched all of them
if($do_old) { 
  $cmd = "grep ^\\> $fa_file | sed 's/^>//' | awk '{ print \$1 }' | sed 's/^.*\://' | sort | sed 's/\\.[0-9]*//' > $fa_file_names";
}
else {
  $cmd = "grep ^\\> $fa_file | sed 's/^>//' | awk '{ print \$1 }' | sed 's/\:.*\$//' | sort | sed 's/\\.[0-9]*//' > $fa_file_names";
}
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($fa_file_names, "Protein accessions for each $cds_or_nt sequence from FASTA file from $fa_file.", $cmd, $log_FH);

# which accessions were lost relative to the previous step? (should be none)
if($do_old) { 
  $cmd = "awk '{ print \$NF }' $efa_file | sed 's/\\.[0-9]*//' | comm -2 -3 - $fa_file_names > $fa_file_lost";
}
else { 
  if($do_dna) { $cmd = "comm -2 -3 $acconly_file $fa_file_names > $fa_file_lost"; }
  else        { $cmd = "awk '{ print \$1 }' $efa_file | sed 's/\\.[0-9]*//' | comm -2 -3 - $fa_file_names > $fa_file_lost"; }
}
RunCommand($cmd, $be_verbose, $log_FH);
OutputFileInfo($fa_file_lost, "$cds_or_nt sequences not fetched", $cmd, $log_FH);

# which accessions were created relative to the previous step? (should be none)
if($do_old) { 
  $cmd = "awk '{ print \$NF }' $efa_file | sed 's/\\.[0-9]*//' | comm -2 -3 $fa_file_names - > $fa_file_created";
}
else { 
  if($do_dna) { $cmd = "comm -2 -3 $fa_file_names $acconly_file > $fa_file_created"; }
  else        { $cmd = "awk '{ print \$1 }' $efa_file | sed 's/\\.[0-9]*//' | comm -2 -3 $fa_file_names - > $fa_file_created"; }
}
RunCommand($cmd, $be_verbose, $log_FH);
OutputFileInfo($fa_file_created, "New $cds_or_nt sequences created when fetching sequences (in error)", $cmd, $log_FH);

# enforce nothing was lost or created
($nlost, $ncreated, $errmsg) = CheckLostAndCreated($fa_file_lost, 0, $fa_file_created, 0); # '0' are max allowed lines in each of these files
if($errmsg ne "") { die $errmsg; }

# output information on this step to stdout and sum file
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($fa_file_names), $nlost, $ncreated, $nsecs, $fa_file), $sum_FH);


###################################
# Stage 2: Validate Stage 1 results
###################################
PrintToStdoutAndFile("#\n", $sum_FH);
PrintToStdoutAndFile("# Stage 2: validating stage 1 results\n", $sum_FH);
PrintToStdoutAndFile("#\n", $sum_FH);
PrintToStdoutAndFile(sprintf("#%-*s  %10s  %10s  %10s  %10s  %s\n", $desc_w-1, " description", "\# output", "\# lost", "\# created", "seconds", "output-file-name"), $sum_FH);
PrintToStdoutAndFile(sprintf("%s  %10s  %10s  %10s  %10s  %s\n", $desc_dashes, "----------", "----------", "----------", "----------", "-----------------"), $sum_FH);

# For step 2.4 we need a de-versioned $efa_file_lost
my $efa_file_lost_acconly = $efa_file_lost . ".acconly"; # same as $efa_file_lost but with version info removed
if(! $do_dna) { 
  $cmd  = "cat $efa_file_lost | sed 's/\\.[0-9]*//' > $efa_file_lost_acconly";
  RunCommand($cmd, $be_verbose, $cmd_FH);
}

########################################################################
# Step 2.1: Validate all accessions exist in protein or nuccore database
########################################################################
my $exists_file         = $out_dir . $symbol . ".exists";
my $exists_file_lost    = $exists_file . ".lost";
my $exists_file_created = $exists_file . ".created";
my $db             = ($do_dna) ? "nuccore" : "protein";
my $protein_or_dna = ($do_dna) ? "DNA" : "protein";
$desc = "Validating that all accessions have a corresponding $protein_or_dna";
$cmd  = "cat $acc_file | epost -db $db -format acc | efetch -format acc | sed 's/\\.[0-9]*//' | sort > $exists_file";
$nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($exists_file, $desc, $cmd, $log_FH);

# $exists_file format, multiple lines of:
# <accession>

# which accessions don't have a protein/nucleotide record? (should be none)
$cmd = "comm -2 -3 $acconly_file $exists_file > $exists_file_lost";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($exists_file_lost, sprintf("Accessions without %ss.", $protein_or_dna), $cmd, $log_FH);

# which accessions are created by this step? (should be none)
$cmd = "comm -2 -3 $exists_file $acconly_file > $exists_file_created";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($exists_file_created, "Accessions unexpectedly created by edirect command.", $cmd, $log_FH);

# enforce nothing was lost or created
($nlost, $ncreated, $errmsg) = CheckLostAndCreated($exists_file_lost, 0, $exists_file_created, 0); 
# in above subroutine call: '0's say it's not okay if we lost or created some accessions

# output information on this step to stdout and sum file
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($exists_file), $nlost, $ncreated, $nsecs, $exists_file), $sum_FH);
if($errmsg ne "") { die $errmsg; }


################################################
# Step 2.2: Validate all accessions have a locus
################################################
my $lexists_file         = $out_dir . $symbol . ".lexists";
my $lexists_file_lost    = $lexists_file . ".lost";
my $lexists_file_created = $lexists_file . ".created";

$desc = "Validating that all accessions have a corresponding locus";
$cmd  = "cat $acc_file | epost -db $db -format acc | efetch -format gpc | xtract -insd INSDSeq_locus | grep . | sort | sed 's/\\.[0-9]*//' > $lexists_file";
$nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($lexists_file, $desc, $cmd, $log_FH);

# $lexists_file format, multiple lines of:
# <accession> <locus>

# which accessions don't have a locus? (should be none)
$cmd   = "awk '{ print \$1 }' $lexists_file | comm -2 -3 $acconly_file - > $lexists_file_lost";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($lexists_file_lost, "Accessions without a locus.", $cmd, $log_FH);

# which accessions are created by this step? (should be none)
$cmd   = "awk '{ print \$1 }' $lexists_file | comm -2 -3 - $acconly_file > $lexists_file_created";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($lexists_file_created, "Accessions unexpectedly created by edirect command.", $cmd, $log_FH);

# enforce nothing was lost or created
($nlost, $ncreated, $errmsg) = CheckLostAndCreated($lexists_file_lost, 0, $lexists_file_created, 0); 
# in above subroutine call: '0's say it's not okay if we lost or created some accessions

# output information on this step to stdout and sum file
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($lexists_file), $nlost, $ncreated, $nsecs, $lexists_file), $sum_FH);
if($errmsg ne "") { die $errmsg; }


#############################################################################
# Step 2.3: Validate we've fetched sequences for all proteins that have a CDS
#############################################################################
my $has_cds_file         = $out_dir . $symbol . ".has_cds";
my $has_cds_file_lost    = $has_cds_file . ".lost";
my $has_cds_file_created = $has_cds_file . ".created";
if(! $do_dna) { 
  $desc = "Validating that we fetched sequences for all proteins with a CDS";
  $cmd  = "cat $exists_file | epost -db protein -format acc | efetch -format gpc | xtract -pattern INSDSeq -ACCN INSDSeq_accession-version -group INSDSeq_feature-table -match INSDFeature_key:CDS -element \"&ACCN\" | sort | sed 's/\\.[0-9]*//' > $has_cds_file";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($has_cds_file, $desc, $cmd, $log_FH);

  # $has_cds_file, multiple lines of:
  # <accession>

  # which accessions did we fetch that don't seem to have CDS? (should be none)
  $cmd = "comm -2 -3 $fa_file_names $has_cds_file > $has_cds_file_lost";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($has_cds_file_lost, "accessions we fetched which don't seem to have a CDS.", $cmd, $log_FH);

  # which accessions have CDS that we didn't fetch? (should be none)
  $cmd = "comm -2 -3 $has_cds_file $fa_file_names > $has_cds_file_created";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($has_cds_file_created, "accessions with a CDS that we failed to fetch.", $cmd, $log_FH);

  # enforce nothing was lost or created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($has_cds_file_lost, 0, $has_cds_file_created, 0);
  # in above subroutine call: '0's say it's not okay if we lost or created some accessions

  # output information on this step to stdout and sum file
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($has_cds_file), $nlost, $ncreated, $nsecs, $has_cds_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }
} # end of 'if($do_dna)'

#########################################################################################
# Step 2.4: Validate we haven't fetched any sequences for proteins that do not have a CDS
#########################################################################################
my $no_cds_file         = $out_dir . $symbol . ".no_cds";
my $no_cds_file_lost    = $no_cds_file . ".lost";
my $no_cds_file_created = $no_cds_file . ".created";
if(! $do_dna) { 
  $desc = "Validating that no seqs were fetched for any proteins w/o a CDS";
  # same command as step 2.3 with 'avoid' replacing 'match'
  $cmd  = "cat $exists_file | epost -db protein -format acc | efetch -format gpc | xtract -pattern INSDSeq -ACCN INSDSeq_accession-version -group INSDSeq_feature-table -avoid INSDFeature_key:CDS -element \"&ACCN\" | sort | sed 's/\\.[0-9]*//' > $no_cds_file";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_file, $desc, $cmd, $log_FH);

  # $no_cds_file, multiple lines of:
  # <accession>

  # which accessions did we not fetch that seem to have CDS? (should be none)
  $cmd = "comm -2 -3 $efa_file_lost_acconly $no_cds_file > $no_cds_file_lost";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_file_lost, "accessions we did not fetch that seem to have a CDS.", $cmd, $log_FH);

  # which accessions don't have CDSs that we did fetch? (should be none)
  $cmd = "comm -2 -3 $no_cds_file $efa_file_lost_acconly > $no_cds_file_created";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_file_created, "accessions without CDS that we did fetch.", $cmd, $log_FH);

  # enforce nothing was lost or created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($no_cds_file_lost, 0, $no_cds_file_created, 0);
  # in above subroutine call: '0's say it's not okay if we lost or created some accessions

  # output information on this step to stdout and sum file
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($no_cds_file), $nlost, $ncreated, $nsecs, $no_cds_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }
} # end of 'if(! $do_dna)'


######################################################################
# Step 2.5: Validate we can explain all proteins that don't have a CDS
######################################################################
# Look at all accessions we weren't about to fetch a CDS from, all should be either:
# (A) start with 'WP_', these are 100% identical at AA level to some other protein and have no CDS
# (B) start with 'YP_' and are dead according to idstat (return an ERROR with idstat)
# (C) are UniProt accessions (TODO: confirm we can get a CDS from these)
my $no_cds_type_a_file         = $out_dir . $symbol . ".no_cds_type_a";           # type A accessions
my $no_cds_type_b_file         = $out_dir . $symbol . ".no_cds_type_b";           # type B accessions
my $no_cds_type_c_file         = $out_dir . $symbol . ".no_cds_type_c";           # type C accessions
my $no_cds_type_a_b_c_file     = $out_dir . $symbol . ".no_cds_type_a_b_c";       # sorted concatenation of type A, B, and C accessions
my $no_cds_type_file_lost      = $out_dir . $symbol . ".no_cds_type.lost";        # accessions that are not type A, B, or C
my $no_cds_type_file_created   = $out_dir . $symbol . ".no_cds_type.created";     # accessions somehow created in this step

if(! $do_dna) { 
  ($seconds, $microseconds) = gettimeofday(); # so we can time this step
  my $step_start_secs = ($seconds + ($microseconds / 1000000.));
  
  # determine type A accessions
  $cmd   = "grep ^WP\_ $no_cds_file | cat - > $no_cds_type_a_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_a_file, "accessions without CDS that begin with \"WP\_\" (type A)", $cmd, $log_FH);
  
  # determine type B accessions
  $cmd   = "grep ^YP\_ $no_cds_file | cat - > $no_cds_type_b_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_b_file, "accessions without CDS that begin with \"YP\_\" (type B)", $cmd, $log_FH);
  
  # determine type C accessions, note that these all begin with a letter followed by a number
  # Reference: http://www.uniprot.org/help/accession_numbers
  # on that page, the regex that matches all UniProt accessions: [OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
  # here we break it down into 3 separate grep calls:
  # Note also that these have no versions
  # ^[OPQ][0-9][A-Z0-9]{3}[0-9]$
  # ^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$
  # ^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){2}$
  $cmd   = "grep '\^[OPQ][0-9][A-Z0-9]\\{3\\}[0-9]\$' $no_cds_file | cat - > $no_cds_type_c_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  $cmd   = "grep '\^[A-NR-Z][0-9][A-Z][A-Z0-9]\\{2\\}[0-9]\$' $no_cds_file | cat - >> $no_cds_type_c_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  $cmd   = "grep '\^[A-NR-Z][0-9]([A-Z][A-Z0-9]\\{2\\}[0-9])\\{2\\}\$' $no_cds_file | cat - >> $no_cds_type_c_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  
  # type A, type B, and type C files should comprise all the accessions, we can't explain any that are missing

  # combine the type A, type B and type C files
  $cmd = "cat $no_cds_type_a_file $no_cds_type_b_file $no_cds_type_c_file | sort > $no_cds_type_a_b_c_file";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_a_b_c_file, "accessions which are missing CDS and can be classified as type A, B or C.", $cmd, $log_FH);
  
  # which accessions are not type A, type B or type C?
  $cmd = "comm -2 -3 $no_cds_file $no_cds_type_a_b_c_file > $no_cds_type_file_lost";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_file_lost, "accessions which are missing CDS for reasons we don't understand.", $cmd, $log_FH);
  
  # which accessions were somehow created during this step? (should be none)
  $cmd = "comm -2 -3 $no_cds_type_a_b_c_file $no_cds_file > $no_cds_type_file_created";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_file_lost, "accessions which were unexpectedly created when investigating missing CDSs.", $cmd, $log_FH);

  # enforce nothing was lost or created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($no_cds_type_file_lost, 0, $no_cds_type_file_created, 0);
  # in above subroutine call: '0's say it's not okay if we lost or created some accessions

  # output information on this step to stdout and sum file
  ($seconds, $microseconds) = gettimeofday(); # so we can time this step
  my $step_end_secs = ($seconds + ($microseconds / 1000000.));

  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, "Validating that proteins w/o CDS are explainable (type A, B or C)", GetNumLinesInFile($no_cds_type_a_b_c_file), $nlost, $ncreated, $step_end_secs - $step_start_secs, $no_cds_type_a_b_c_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }
}

##########
# Conclude
##########
Conclude($start_secs, $sum_FH, $log_FH);

exit 0;

#############
# SUBROUTINES
#############
#
# Subroutine: RunCommand()
# Args:       $cmd:            command to run, with a "system" command;
#             $be_verbose:     '1' to output command to stdout before we run it, '0' not to
#             $FH:             file handle to output commands to
#
# Returns:    amount of time the command took, in seconds
# Dies:       if $cmd fails

sub RunCommand {
  my $sub_name = "RunCommand()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cmd, $be_verbose, $FH) = @_;

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }
  if(defined $FH) { 
    print $FH $cmd . "\n";
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  system($cmd);
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if($? != 0) { die "ERROR command failed:\n$cmd\n"; }

  return ($stop_time - $start_time);
}

# Subroutine: GetNumLinesInFile()
# Purpose:    Counts the number of non-blank lines 
#             in file $file.
# Args:       $file: file to get number of lines of

sub GetNumLinesInFile {
  my $sub_name = "GetNumLinesInFile()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($file) = (@_);

  if(! -e $file) { die "ERROR file $file does not exist in CheckNumLinesInFile.\n"; }

  my $nlines = `grep -c . $file`;
  chomp $nlines;

  return $nlines;
}

# Subroutine: CheckLostAndCreated()
# Purpose:    Determine number of lines in 'lost' and 'created' files and check
#             if it exceeds our limits. 
# Args:       $lost_file:    lost file, each line represents a lost accession/id
#             $maxlost:      max number of allowed lost lines, -1 for no max
#             $created_file: created file, each line represents a newly created accession/id
#             $maxcreated:   max number of allowed created lines, -1 for no max
# 
# Returns: Three values:
#          $nlost:    number of lines in $lostfile
#          $ncreated: number of lines in $createdfile
#          $errmsg:   "" if no error, else filled with reason caller should die
#                     will be != "" if we find > maximum number of lines in 
#                     $lost_file and $created_file
#     
# Dies:    if $lostfile or $createdfile do not exist
sub CheckLostAndCreated { 
  my $sub_name = "OutputLostAndCreated()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($lost_file, $maxlost, $created_file, $maxcreated) = (@_);

  if(! -e $lost_file)    { die "ERROR in $sub_name, lost file $lost_file does not exist"; }
  if(! -e $created_file) { die "ERROR in $sub_name, created file $created_file does not exist"; }

  $nlost    = GetNumLinesInFile($lost_file);
  $ncreated = GetNumLinesInFile($created_file);
  
  my $errmsg = "";
  if(($maxlost    != -1) && ($nlost    > $maxlost))    { $errmsg = "FAILED: #lost ($nlost) > $maxlost in file $lost_file"; }
  if(($maxcreated != -1) && ($ncreated > $maxcreated)) { $errmsg = "FAILED: #created ($ncreated) > $maxcreated in file $created_file"; }

  return ($nlost, $ncreated, $errmsg);
}

# Subroutine: PrintToStdoutAndFile()
# Purpose:    Output a string to stdout and to a file handle.
# Args:       $string:    string to print to stdout and $FH
#             $FH:        file handle to print to, unless undef
#
# Returns: void
sub PrintToStdoutAndFile { 

  my $sub_name = "PrintToStdoutAndFile()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($string, $FH) = (@_);

  print $string;
  if(defined $FH) { print $FH $string; }

  return;
}

# Subroutine: OutputFileInfo()
# Purpose:    Add output file name to an array and description
#             to a hash, and output it to stdout and FH (if defined)
# Args:       $name:    name of output file
#             $desc:    description of file
#             $cmd:     command used to generate file $name
#             $FH:      file handle to output to
# Returns: void
sub OutputFileInfo {
  my $sub_name  = "OutputFileInfo()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($name, $desc, $cmd, $FH) = (@_);

  if(! defined $FH) { die "ERROR file handle not defined in $sub_name"; }

  print $FH ("--------------------------------------------------\n");
  print $FH ("filename:        $name\n");
  print $FH ("description:     $desc\n");
  if($cmd ne "") { 
    print  $FH ("command:         $cmd\n");
    printf $FH ("number-of-lines: %d\n", GetNumLinesInFile($name));
  }
  else { 
    print $FH ("command:         N/A\n");
  }
  print $FH ("--------------------------------------------------\n");

  return;
}

# Subroutine: LengthFromCoords()
# Purpose:    Determine the length of a CDS given a coordinate 
#             string $coords.
# Args:       $coords:  the coords string
# Returns:    the length (in nt) of the CDS described by $coords
sub LengthFromCoords {
  my $sub_name  = "LengthFromCoords()"; 
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords) = (@_);
  #print("IN LengthFromCoords(): $coords\n");

  my $orig_coords = $coords; # only nec so we can output it if there's an error
  # example $coords values:
  # U77714.1:1..867
  # complement(U77714.1:1..867)
  # join(U43876.1:608..688,U43877.1:104..175)
  # complement(join(U43876.1:608..688,U43877.1:104..175))
  # 
  # the 'join(' and 'complement' are useless
  $coords =~ s/^complement\(//;
  $coords =~ s/^join\(//;
  $coords =~ s/\)$//;
  $coords =~ s/\)$//;
  my @exon_A = split(",", $coords);

  if(scalar(@exon_A) == 0) { # single exon CDS
    push(@exon_A, $coords);
  }
  my $len = 0;
  foreach my $exon (@exon_A) { 
    if($exon =~ /^\S+\.?\d*\:\<?(\d+)\.\.\>?(\d+)$/) { 
      my ($start, $end) = ($1, $2);
      #print("\t$start $end\n");
      if($end < $start) { 
        die "ERROR found end coordinate ($end) less than start coord ($start) in coords string $orig_coords";
      }
      $len += ($end-$start+1);
    }
    else { 
      die "ERROR unable to parse exon $exon in coords string $orig_coords when trying to get length in LengthFromCoords()";
    }
  }
  # printf("LengthFromCoords(): $orig_coords $len\n");
  return $len;
}

# Subroutine: removeGivenIdstat()

# Purpose: Given an output file from 'idstat' ($idstat_file) and a
#             file that lists all accessions that 'idstat' was run on
#             (*in order*) ($allacc_file) output a new file with
#             some accessions removed ($out_file). Remove all accessions
#             that meet any of the following criteria:
#             - have 'Dead' in the 'State' column.
#             - have any value except 'No' in the 'Wdrw' (withdrawn) column
#             - have 'Yes' or 'Perm' in the 'Supp' (suppressed) column
#
# Args:       $idstat_file: the idstat output file
#             $allacc_file: the accession file containing all accession+versions 
#                           that idstat was run on the accession-onlys of to produce
#                           idstat_file
#             $out_file:    file to output non-suppressed accession+versions to
# Returns:    number of seconds required to do this parsing
# Dies:       if we have a problem parsing $idstat_file
sub removeGivenIdstat {
  my $sub_name  = "removeSuppressedGivenIdstat()";
  my $nargs_exp = 3;
  
  my ($allacc_file, $idstat_file, $out_file) = (@_);
  my ($seconds, $microseconds) = gettimeofday();
  my $start_secs = ($seconds + ($microseconds / 1000000.));

  my @accn_A  = (); # array of accessions, in order read from idstate fule
  my @state_A = (); # array of 'State' column values same order as accn_A
  my @wdrw_A  = (); # array of 'Wdrw' column values same order as accn_A
  my @supp_A  = (); # array of 'Supp' column values same order as accn_A

  # first store 'suppressed' value for all accessions
  open(IDSTAT, $idstat_file) || die "ERROR unable to open $idstat_file"; 
  while(my $line = <IDSTAT>) { 
    #example:
    # Accession        GI         SatKey    SatName   State Owner                 Wdrw Supp Date Loaded           Username
    # ---------------- ---------- --------- --------- ----- --------------------- ---- ---- --------------------- --------
    # AAH45158         28175331   122024391 NCBI      Alive NCBI-MGC              No   No   9/23/2014 3:10:49PM   kans

    # chew up everything until we hit a line that starts with '---'
    if($line =~ m/^\-\-\-/) { # next line is special, its the most recent entry for an accession
      my $spec_line = <IDSTAT>;
      # AAH45158         28175331   122024391 NCBI      Alive NCBI-MGC              No   No   9/23/2014 3:10:49PM   kans
      my @elA = split(/\s+/, $spec_line);
      my ($accn, $state, $withdrawn, $suppressed) = ($elA[0], $elA[4], $elA[6], $elA[7]);
      push(@accn_A,  $accn);
      push(@state_A, $state);
      push(@wdrw_A,  $withdrawn);
      push(@supp_A,  $suppressed);
    }
  }
  close(IDSTAT);

  open(ACC,       $allacc_file) || die "ERROR unable to open $allacc_file for reading";
  open(OUT, ">" . $out_file)    || die "ERROR unable to open $out_file for writing";
  my $ctr = 0;
  while(my $line = <ACC>) { 
    if($line =~ m/\w/) { # skip blank lines
      if((! defined $accn_A[$ctr]) || (! defined $supp_A[$ctr])) { 
        die "ERROR ran out of parsed idstat accessions prematurely!"; 
      }
      chomp $line;
      my $accn_version = $line;
      my $accn2match = $accn_A[$ctr];
      if($accn_version !~ m/$accn2match\.?\d*/) { 
        die "ERROR, problem parsing $idstat_file and $allacc_file in $sub_name, trying to match $accn2match to $accn_version"; 
      }

      my $do_remove = 0; # set to '1' below if removed
      # parse 'state' column
      if($state_A[$ctr]   eq "Dead") { $do_remove = 1; } # dead, remove it

      # parse 'wdrw' column
      if($wdrw_A[$ctr]    ne "No")   { $do_remove = 1; }

      # parse 'supp' column
      if($supp_A[$ctr]    eq "Yes")  { $do_remove = 1; } # suppressed, remove it
      elsif($supp_A[$ctr] eq "Perm") { $do_remove = 1; } # permanently suppressed, remove it
      elsif($supp_A[$ctr] eq "No")   { ; } # not suppressed
      else                           { die "ERROR unable to parse Supp column for $accn2match, read $supp_A[$ctr]"; }

      if(! $do_remove) { 
        print OUT $accn_version . "\n";
      }

      $ctr++;
    }
  }
  close ACC;
  close OUT;

  ($seconds, $microseconds) = gettimeofday(); # so we can time this step
  my $end_secs = ($seconds + ($microseconds / 1000000.));

  return ($end_secs - $start_secs);
}
  

# Subroutine: Conclude()
# Purpose:    Print out conclusion text and close file handles in preparation for exit.
# Args:       $start_secs: number of seconds since epoch and start of this script running
#             $sum_FH:     file handle for summary file
#             $log_FH:     file handle for log file
# Returns:    void
sub Conclude {
  my $sub_name  = "Conclude()";
  my $nargs_exp = 3;
  
  my ($start_seconds, $sum_FH, $log_FH) = (@_);
  ($seconds, $microseconds) = gettimeofday();
  my $end_secs = ($seconds + ($microseconds / 1000000.));

  PrintToStdoutAndFile("#\n", $sum_FH);
  PrintToStdoutAndFile("# Output files created by this script with brief descriptions listed in log file:  $log_file\n",    $sum_FH);
  PrintToStdoutAndFile("# This output printed to stdout written to summary file:                           $sum_file\n",    $sum_FH);
  PrintToStdoutAndFile("# All commands executed by this script listed in cmd file:                         $cmd_file\n",    $sum_FH);
  PrintToStdoutAndFile("# All output files created in directory:                                           \.\/$out_dir\n", $sum_FH);
  PrintToStdoutAndFile(sprintf("# Total seconds elapsed:                                                           %.1f\n", $end_secs-$start_secs), $sum_FH);
  PrintToStdoutAndFile("#\n",     $sum_FH);
  PrintToStdoutAndFile("#[ok]\n", $sum_FH);
  
  close $log_FH;
  close $sum_FH;
  return;
}
