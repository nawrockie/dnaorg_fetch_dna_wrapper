# EPN, Thu Apr 16 09:23:11 2015
# pwd: /panfs/pan1/dnaorg/programs/15_0415_dnaorg_fetch_dna_wrapper
# 
# Example of running dnaorg_fetch_dna_wrapper.pl to fetch
# nucleotide sequences. The program can be used to:
# mode A) fetch CDS sequences of proteins that link to a gene symbol
# mode B) fetch nucleotide sequences of noncoding genes that link to a gene symbol
# mode C) fetch CDS sequences of protein accessions listed in a file
# mode D) fetch nucleotide sequences of nucleotide accessions listed in a file
#
# This is a wrapper Perl script which calls other programs,
# including those housed/described in the following directories:
#
# /panfs/pan1/dnaorg/programs/15_0219_edirect/
# /panfs/pan1/dnaorg/programs/15_0324_esl_fetch_cds/
# 
##################
# Previous version
##################
# dnaorg_fetch_cds_wrapper.pl: see /panfs/pan1/dnaorg/programs/15_0310_dnaorg_fetch_cds_wrapper/
#
#######################
# More information
#######################
#
# For more information, see the 00NOTES.* files in the above listed
# 'dnaorg/programs/' directories.
#
# Also, see /home/nawrocke/notebook/15_0415_dnaorg_fetch_dna_wrapper/00LOG.txt
# and /home/nawrocke/notebook/15_0310_dnaorg_fetch_cds_wrapper/00LOG.txt
# for notes on development and testing of this program and its predecessor.
# 
#######################
# Prerequisites
#######################
# 
# Directories that include the BioEasel perl modules must be part of your
# $PERL5LIB environment variable in order for esl-fetch-cds.pl to work.
# 
# For bash shell users
source /panfs/pan1/dnaorg/programs/setup-bio-easel.bash.sh
# For C shell or C shell compatible users
source /panfs/pan1/dnaorg/programs/setup-bio-easel.csh.sh
#
#######################
# Usage and options
#######################
# The default output is informative about how to use the script:
# $ perl dnaorg_fetch_dna_wrapper.pl 
##
##dnaorg_fetch_dna_wrapper.pl:
##
## This script fetches DNA sequences from GenBank. It can be
## run in three different modes:
##
## Mode 1: fetch CDS sequences for protein database records that link
##         to a symbol for a protein coding gene in the Gene database,
##         or if no matches exist in protein, repeat search in nuccore
##         and fetch all nucleotide records that link to the noncoding
##         symbol in the Gene database.
##
##   usage: 'perl dnaorg_fetch_dna_wrapper.pl [OPTIONS] <symbol>'
##
##
## Mode 2: fetch CDS sequences for a list of protein accessions
##
##   usage: 'perl dnaorg_fetch_dna_wrapper.pl [OPTIONS] -d <name_for_outdir> -plist <file_with_list_of_protein_accessions>'
##
##
## Mode 3: fetch nucleotide sequences for a list of nuccore accessions
##
##   usage: 'perl dnaorg_fetch_dna_wrapper.pl [OPTIONS] -d <name_for_outdir> -ntlist <file_with_list_of_nucleotide_accessions>'
##
## NOTE: NCBI Gene symbols can include whitespace (' ') but this script
##       does not allow them in the <symbol> command line argument. To
##       specify a ' ' in the <symbol> replace it with a '~' on the
##       command line. (As of 04/28/15 no '~' exist in NCBI Gene symbols
##
## BASIC OPTIONS:
##  -f      : force; if dir <symbol> exists, overwrite it
##  -v      : be verbose; output commands to stdout as they're run
##  -d <s>  : define output directory as <s>, not <symbol>
##  -nt     : search for matches to symbol in ONLY the nuccore db, not the protein db
##  -notnt  : if zero matches for symbol in the protein db, DO NOT retry using the nuccore db
##  -nosyn  : only include records for which the symbol is the primary symbol, not the synonym
##
## OPTIONS THAT ENABLE ALTERNATE MODES:
##  -plist  : <symbol> is really a list of protein accessions,    requires -d option too
##  -ntlist : <symbol> is really a list of nucleotide accessions, requires -d option too
##  -num    : determine number of matching protein/nucleotide accessions, then exit
##
## OPTIONS THAT AFFECT FAILURE/WARNING OF PRE-DETERMINED SYMBOLS:
##  -ffile <f> : fail          if a symbol listed in <f> is used as input symbol [default: /panfs/pan1/dnaorg/share/all.multi-synonym.list]
##  -sfile <f> : require symbol be primary symbol for symbols listed in file <f> [default: /panfs/pan1/dnaorg/share/all.primary-and-synonym.list]
##  -noffile   : do not fail for symbols listed in a file
##  -nosfile   : do not skip synonyms for symbols listed in a file
##
## EXPERIMENTAL/ADVANCED OPTIONS:
##  -up     : additional run experimental code for fetching non-CDS UniProt CDS via xrefs
##  -old    : use extract_fasta_multi_exon instead of esl-fetch-cds.pl
##
##
## This script will create a directory called <symbol> (for modes 1 and 2)
## or called <name_for_outdir> (for modes 3 and 4) and populate it with
## several output files, including:
##
##  <symbol>.log: list of all output files and brief descriptions
##  <symbol>.cmd: list of all commands run by this script
##  <symbol>.sum: copy of all stdout from this script
##
## For modes 2 and 3, the files will start with <name_for_outdir>
## instead of <symbol>.
##
#############################
# Example commands and output
#############################
# 
# Example 1: fetch CDS sequences of proteins that link to a gene symbol:
# 
# COMMAND:
dnaorg_fetch_dna_wrapper.pl smn1
#
# OUTPUT (each line has been prefixed with a '#'):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## dnaorg_fetch_dna_wrapper.pl: Fetch DNA sequences from GenBank
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## command: ./dnaorg_fetch_dna_wrapper.pl -f smn1
## date:    Mon Apr 20 13:47:26 2015
## option:  forcing overwrite of smn1 dir/file [-f]
## symbol:  smn1
##
## Stage 1: preparing and fetching sequences
##
## description                                                               # output      # lost   # created     seconds  output-file-name
##-----------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Protein_accessions_fetched_from_protein_database                                 155           0         155         4.1  smn1/smn1.all.acc
#Non-suppressed_protein_accessions_fetched_from_protein_database                  155           0           0        10.5  smn1/smn1.acc
#Protein_accessions_that_have_CDS_annotation                                      145          10           0         9.3  smn1/smn1.efa
#CDS_sequences_in_FASTA_format                                                    145           0           0        51.1  smn1/smn1.fa
##
## Stage 2: validating stage 1 results
##
## description                                                               # output      # lost   # created     seconds  output-file-name
##-----------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Validating_that_all_accessions_have_a_corresponding_protein_record               155           0           0        54.9  smn1/smn1.exists
#Validating_that_all_accessions_have_a_corresponding_locus                        155           0           0         6.0  smn1/smn1.lexists
#Validating_that_we_fetched_sequences_for_all_proteins_with_a_CDS                 145           0           0         9.3  smn1/smn1.has_cds
#Validating_that_no_seqs_were_fetched_for_any_proteins_without_a_CDS               10           0           0        10.7  smn1/smn1.no_cds
#Validating_that_proteins_without_CDS_are_explainable_(type_A,_B_or_C)             10           0           0         0.1  smn1/smn1.no_cds_type_a_b_c
##
## Output files created by this script with brief descriptions listed in log file:  smn1/smn1.log
## This output printed to stdout written to summary file:                           smn1/smn1.sum
## All commands executed by this script listed in cmd file:                         smn1/smn1.cmd
## All output files created in directory:                                           ./smn1/
## Total seconds elapsed:                                                           105.3
##
## Type of records: protein
##[ok]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Explanation of columns:
# '# output':  number of accessions that are returned from each step. 
# '# lost'  :  number of accessions from previous step not returned in current step (stage 1)
#              OR number of accessions not validated (stage 2)
# '# created': number of new accessions created at this step. This should be 0 for
#              ALL but the first step (script will die in error if it's not)
# 'seconds':   seconds to complete this step
# 'output-file-name': output file created by this step.
# 
# -------------------------------
# Example 2: fetch nucleotide sequences of DNA sequences that link to a noncoding gene symbol:
# COMMAND: 
#           
dnaorg_fetch_dna_wrapper.pl SNORA71D
#
# Note: alternatively the -nt option could be used to specify that only the 
# nuccore database be searched.
#
# OUTPUT (each line has been prefixed with a '#'):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## dnaorg_fetch_dna_wrapper.pl: Fetch DNA sequences from GenBank
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## command: ./dnaorg_fetch_dna_wrapper.pl -f SNORA71D
## date:    Mon Apr 20 13:51:10 2015
## option:  forcing overwrite of SNORA71D dir/file [-f]
## symbol:  SNORA71D
##
## Stage 1: preparing and fetching sequences
##
## description                                                                   # output      # lost   # created     seconds  output-file-name
##---------------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Protein_accessions_fetched_from_protein_database                                       0           0           0         0.9  SNORA71D/SNORA71D.all.acc
#Nucleotide_accessions_fetched_from_nuccore_database                                    5           0           5         2.1  SNORA71D/SNORA71D.all.acc
#Non-suppressed_nucleotide_accessions_fetched_from_nuccore_database                     5           0           0         1.2  SNORA71D/SNORA71D.acc
#Nucleotide_sequences_in_FASTA_format                                                   5           0           0       135.3  SNORA71D/SNORA71D.fa
##
## Stage 2: validating stage 1 results
##
## description                                                                   # output      # lost   # created     seconds  output-file-name
##---------------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Validating_that_all_accessions_have_a_corresponding_DNA_record                         5           0           0       137.3  SNORA71D/SNORA71D.exists
#Validating_that_all_accessions_have_a_corresponding_locus                              5           0           0         3.0  SNORA71D/SNORA71D.lexists
##
## Output files created by this script with brief descriptions listed in log file:  SNORA71D/SNORA71D.log
## This output printed to stdout written to summary file:                           SNORA71D/SNORA71D.sum
## All commands executed by this script listed in cmd file:                         SNORA71D/SNORA71D.cmd
## All output files created in directory:                                           ./SNORA71D/
## Total seconds elapsed:                                                           147.4
##
## Type of records: nucleotide
##[ok]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Note that there are only two files created for the validation stage.
# Currently, there are fewer validation steps performed for nucleotide
# sequences. The other protein validation steps do not apply to 
# nucleotide sequences.
#
# -------------------------------
# Example 3: fetch a list of protein accessions using the --plist option:
#
# COMMAND:
dnaorg_fetch_dna_wrapper.pl -d samp5 -plist samp5.list
#
# OUTPUT (each line prefixed with a '#'):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## dnaorg_fetch_dna_wrapper.pl: Fetch DNA sequences from GenBank
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## command: ./dnaorg_fetch_dna_wrapper.pl -d samp5 -plist samp5
## date:    Mon Apr 20 14:04:09 2015
## option:  output directory specified as samp5 [-d]
## option:  samp5.list is a list of protein accessions, not a symbol [-plist]
##
## Stage 1: preparing and fetching sequences
##
## description                                                                # output      # lost   # created     seconds  output-file-name
##------------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Sorted_accessions_from_file_samp5.list                                              5           0           5         0.0  samp5/samp5.all.acc
#Non-suppressed_protein_accessions_fetched_from_protein_database                     5           0           0         1.3  samp5/samp5.acc
#Protein_accessions_that_have_CDS_annotation                                         5           0           0         3.2  samp5/samp5.efa
#CDS_sequences_in_FASTA_format                                                       5           0           0         6.5  samp5/samp5.fa
##
## Stage 2: validating stage 1 results
##
## description                                                                # output      # lost   # created     seconds  output-file-name
##------------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Validating_that_all_accessions_have_a_corresponding_protein_record                  5           0           0         8.5  samp5/samp5.exists
#Validating_that_all_accessions_have_a_corresponding_locus                           5           0           0         3.4  samp5/samp5.lexists
#Validating_that_we_fetched_sequences_for_all_proteins_with_a_CDS                    5           0           0         3.2  samp5/samp5.has_cds
#Validating_that_no_seqs_were_fetched_for_any_proteins_without_a_CDS                 0           0           0         3.2  samp5/samp5.no_cds
#Validating_that_proteins_without_CDS_are_explainable_(type_A,_B_or_C)               0           0           0         0.1  samp5/samp5.no_cds_type_a_b_c
##
## Output files created by this script with brief descriptions listed in log file:  samp5/samp5.log
## This output printed to stdout written to summary file:                           samp5/samp5.sum
## All commands executed by this script listed in cmd file:                         samp5/samp5.cmd
## All output files created in directory:                                           ./samp5/
## Total seconds elapsed:                                                           23.1
##
## Type of records: protein
##[ok]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This example could be repeated with a file listing nucleotide
# accessions instead of protein accessions with the use of the -ntfile
# option instead of -pfile.
#
# -------------------------------
# Example 4: fetch only the number of sequences that match, not 
#            the sequences themselves:
#
# COMMAND:
dnaorg_fetch_dna_wrapper.pl -num SNORA71D
#
# OUTPUT (each line prefixed with a '#'):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## dnaorg_fetch_dna_wrapper.pl: Fetch DNA sequences from GenBank
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## command: ./dnaorg_fetch_dna_wrapper.pl -f -num SNORA71D
## date:    Mon Apr 20 14:14:25 2015
## option:  forcing overwrite of SNORA71D dir/file [-f]
## option:  determining number of matching protein accessions, then exiting [-num]
## symbol:  SNORA71D
##
## Special mode (-num): determining number of matching protein accessions, then exiting.
##
## description                                                                   # output      # lost   # created     seconds  output-file-name
##---------------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Protein_accessions_fetched_from_protein_database                                       0           0           0         0.7  SNORA71D/SNORA71D.all.acc
#Nucleotide_accessions_fetched_from_nuccore_database                                    5           0           5         2.0  SNORA71D/SNORA71D.all.acc
#Non-suppressed_nucleotide_accessions_fetched_from_nuccore_database                     5           0           0         0.3  SNORA71D/SNORA71D.acc
##
#Number_of_nucleotide_records: 5
##
## Output files created by this script with brief descriptions listed in log file:  SNORA71D/SNORA71D.log
## This output printed to stdout written to summary file:                           SNORA71D/SNORA71D.sum
## All commands executed by this script listed in cmd file:                         SNORA71D/SNORA71D.cmd
## All output files created in directory:                                           ./SNORA71D/
## Total seconds elapsed:                                                           3.0
##
## Type of records: nucleotide
##[ok]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##############
# Output files
##############
# 
# About 30 files will be created in a newly created directory
# called '<symbol>'. These will all be prefixed with '<symbol>'.
# The file <symbol>.log contains a description of all of these files
# along with the command that generated it.
# For example:
#
#--------------------------------------------------
#filename:        smn1/smn1.acc
#description:     Protein accessions fetched from protein database
#command:         esearch -db protein -query "smn1 [GENE]" | efetch -format acc | sort > smn1/smn1.acc
#number-of-lines: 156
#--------------------------------------------------
# 
# Other important files include:
#
# <symbol>.sum: the standard output from the program (what is printed above as OUTPUT)
# <symbol>.cmd: all the commands run by the script, one per line
# <symbol>.acc: accessions that survive step 1.1
# <symbol>.fa:  all CDS sequences in fasta format
# 
# For complete list with descriptions see <symbol>.log.
# 
#
######################
# Programs called by dnaorg_fetch_cds_wrapper.pl
# - These programs are called by default. An additional 2 programs
#   are called and listed below if the -old option is used.
# - This list excludes unix commands, e.g. 'grep', 'comm', 'sort'.
# - idfetch section copied from /panfs/pan1/dnaorg/programs/15_0303_extract_fasta2_multiple_exons/00NOTES.sh
######################
#
# esl-fetch-cds.pl
# type:     Perl script that uses the Inline C-equipped 
#           Bio-Easel perl modules (meaning it uses C too).
# author:   EPN
# location: /panfs/pan1/dnaorg/programs/15_0324_esl_fetch_cds/esl-fetch-cds.pl (script)
#           /panfs/pan1/dnaorg/programs/Bio-Easel/blib/lib, /panfs/pan1/dnaorg/programs/Bio-Easel/blib/arch (perl modules)
#
# Fetches CDS sequences given an input file that specifies the coordinates
# (the format of a 'coded_by' INSDQualifier field of a 'CDS' INSDFeature_key
# in an NCBI 'gpc' format. An example of fetching such information using edirect
# tools:
#
# 'esearch -db protein -query AAH45158.1 | efetch -format gpc | xtract -insd CDS coded_by' 
# AAH45158.1	BC045158.1:4..870
# 
# For more info: 'esl-fetch-cds.pl'
# 
# NOTE: this program is called internally by esl-fetch-cds.pl, not by
#       the wrapper dnaorg_fetch_cds.pl.
#
######################
#
# idstat    
# type:     ?
# author:   ?
# location: /netopt/genbank/subtool/bin/idstat
# 
# Reports statistics/status of specified genbank accessions/identifiers. 
# Used here to make sure that all proteins without a CDS return 
# an error from idstat indicating those proteins don't actually
# exist in GenBank (NOTE: I'm a little unsure of what it actually means
# when an accession returns an error, how could it not exist if edirect
# fetched it?).
# 
# Example:
# > /netopt/genbank/subtool/bin/idstat -i PUBSEQ_OS -a O02771.1
#   /netopt/genbank/subtool/bin/idstat.pl : ERROR : Accession/OSLT >O02771.1< was not found in ID : skipping this sequence
#
######################
# Programs called by dnaorg_fetch_cds_wrapper.pl if the -old option
# is used.
# - id_fasta.pl and extract_fasta_multi_exon sections copied from /panfs/pan1/dnaorg/programs/15_0303_extract_fasta2_multiple_exons/00NOTES.sh
#
######################
# id_fasta.pl 
# type:     Perl script
# author:   ?
# location: local, in this directory
# 
# Renames fasta sequences so the IDs are Genbank IDs,
# modifying only the name (first token attached to the '>')
# and regurgitating everything else. Original name is 
# included as the first token in the sequence description,
# followed by the original description.
# 
# Usage: id_fasta.pl <fasta file, from 'idfetch'>
#
######################
#
# extract_fasta_multi_exon
# type:     C program
# location: local, in this directory
# author:   Richa Agarwala [extract_fasta]
#           Eric Nawrocki [modified original by Richa]
#
# Extracts specified intervals/strands from DNA sequences in a
# provided fasta file. Those intervals can be multi-exonic (this is
# the main difference between extract_fasta_multi_exon and
# extract_fasta).
#
# Usage:
# > ./extract_fasta_multi_exon
# extract_fasta_multi_exon <interval_list> [<fa_file>]
# Format of <interval_list>, each line must look like this:
#
# <accession/id> <num-pieces (n)> <start_1> <end_1> <start_2> <end_2> ... <start_n> <end_n> <strand> <optional-extra-string-to-add-to-name>
#
# For all n <start> must be <= <end> and 
# <start_j> must be > <end_i> for all i < j.
#
# Compilation:
#  > gcc -o extract_fasta_multi_exon extract_fasta_multi_exon.c
#
# Source of 'extract_fasta_multi_exon':
# /home/nawrocke/notebook/15_0303_dnaorg_extract_fasta_multiple_exons/
# (notes in the 00LOG.txt file in the above dir).
#
# Source of original 'extract_fasta':
# ~agarwala/programs/extract_fasta.c
#
############################
# Input files
############################
#
# NONE
#
#############################################
# Last updated: EPN, Mon Apr 20 14:23:40 2015
#############################################
#
# This directory is under git control. Use
# git commands to see revision history.
#
############################################
