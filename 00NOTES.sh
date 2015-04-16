# EPN, Thu Apr 16 09:23:11 2015
# pwd: /panfs/pan1/dnaorg/programs/15_0415_dnaorg_fetch_dna_wrapper
# 
# Example of running dnaorg_fetch_dna_wrapper.pl to fetch
# nucleotide sequences. The program can be used to:
# mode A) fetch CDS sequences of proteins that link to a gene symbol
# mode B) fetch CDS sequences of protein accessions listed in a file
# mode C) fetch DNA sequences of nucleotide accessions listed in a file
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
# Also, see /home/nawrocke/notebook/15_0310_dnaorg_fetch_cds_wrapper/00LOG.txt
# for notes on development and testing of this program.
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
# $ perl dnaorg_fetch_dna_wrapper.pl 
# dnaorg_fetch_cds_wrapper.pl [OPTIONS] <symbol (from GENE database, e.g. 'infB'>
# dnaorg_fetch_cds_wrapper.pl [OPTIONS]      -d <symbol> -alist <file with list of protein accessions>
# dnaorg_fetch_cds_wrapper.pl [OPTIONS] -dna -d <symbol> -alist <file with list of nucleotide accessions>
#	OPTIONS:
#		-f         : force; if dir <symbol> exists, overwrite it
#		-v         : be verbose; output commands to stdout as they're run
#		-d <s>     : define output directory as <s>, not <symbol>
#		-np        : determine number of matching protein accessions then exit
#		-alist     : <symbol> is really a list of accessions, requires -d option too
#		-up        : additional run experimental code for fetching non-CDS UniProt CDS via xrefs
#		-dna       : -alist contains nucleotide accessions, not protein (requires -alist and -d)
#		-old       : use extract_fasta_multi_exon instead of esl-fetch-cds.pl
#
#############################
# Example commands and output
#############################
# 
#-------------------------------------------------------------------
# mode A) fetch CDS sequences of proteins that link to a gene symbol
#-------------------------------------------------------------------
dnaorg_fetch_dna_wrapper.pl smn1
#
# This will create a directory called 'smn1' and populate it with
# about 30 files. Many of these will be empty, and are only created to
# check for errors or unexpected results. The output of the example
# command and descriptions of several of the important files are
# below. 
# 
# What the script actually does:
# Stage 1: Fetches CDS sequences for the given <symbol>.
#   Step 1.1: Fetches all protein accessions that match a query to our symbol 
#   Step 1.2: Determines the CDS sequences/coordinates that correspond to the proteins from 1.1
#   Step 1.3: Converts the extract_fasta_multi_exon input file to idfetch input format.
#   Step 1.4: Extract the sequences using esl-fetch-cds.pl (or extract_fasta_multi_exon if -old)
#
# Stage 2: Validation of Stage 1 results
#   Step 2.1: Validate all accessions have proteins
#   Step 2.2: Validate all accessions have a locus
#   Step 2.3: Validate we've fetched sequences for all proteins that have a CDS
#   Step 2.4: Validate we haven't fetched any sequences for proteins that do not have a CDS
#   Step 2.5: Validate we can explain all proteins that don't have a CDS
# 
# Throughout these stages, many files are compared to make sure their
# differences are expected. If they're unexpected, the script dies in 
# error.
# 
# OUTPUT (each line has been prefixed with a '#'):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## dnaorg_fetch_cds_wrapper.pl: Fetch CDS sequences for a given protein gene symbol
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## command: ./dnaorg_fetch_dna_wrapper.pl smn1
## date:    Thu Apr 16 09:31:48 2015
## symbol:  smn1
##
## Stage 1: preparing and fetching CDS sequences
##
## description                                                          # output      # lost   # created     seconds  output-file-name
##------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Protein accessions fetched from protein database                            155           0         155         4.5  smn1/smn1.all.acc
#Non-suppressed protein accessions fetched from protein database             155           0           0        16.0  smn1/smn1.acc
#Protein accessions that have CDS annotation                                 145          10           0         8.0  smn1/smn1.efa
#CDS sequences in FASTA format                                               145           0           0        54.6  smn1/smn1.fa
##
## Stage 2: validating stage 1 results
##
## description                                                          # output      # lost   # created     seconds  output-file-name
##------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Validating that all accessions have a corresponding protein                 155           0           0        59.0  smn1/smn1.exists
#Validating that all accessions have a corresponding locus                   155           0           0         8.1  smn1/smn1.lexists
#Validating that we fetched sequences for all proteins with a CDS            145           0           0         7.8  smn1/smn1.has_cds
#Validating that no seqs were fetched for any proteins w/o a CDS              10           0           0         7.1  smn1/smn1.no_cds
#Validating that proteins w/o CDS are explainable (type A, B or C)            10           0           0         0.1  smn1/smn1.no_cds_type_a_b_c
##
## Output files created by this script with brief descriptions listed in log file:  smn1/smn1.log
## This output printed to stdout written to summary file:                           smn1/smn1.sum
## All commands executed by this script listed in cmd file:                         smn1/smn1.cmd
## All output files created in directory:                                           ./smn1/
## Total seconds elapsed:                                                           111.1
##
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
#-------------------------------------------------------------------
# mode B) HERE HERE HERE
#-------------------------------------------------------------------
#
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
######################################
# Special (alternative) running modes:
######################################
#
# -np option: determine number of matching protein accessions then exit.
# 
# Example:
#
# $ dnaorg_fetch_cds_wrapper.pl -np -f BCAR4
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## dnaorg_fetch_cds_wrapper.pl: Fetch CDS sequences for a given protein gene symbol
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## command: ./dnaorg_fetch_cds_wrapper.pl -f -np BCAR4
## date:    Mon Mar 30 10:28:02 2015
## option:  forcing overwrite of BCAR4 dir/file [-f]
## option:  determining number of matching protein accessions, then exiting
## symbol:  BCAR4
##
## Special mode (-np): determining number of matching protein accessions, then exiting.
##
## description                                                           # output      # lost   # created     seconds  output-file-name
##-------------------------------------------------------------------  ----------  ----------  ----------  ----------  -----------------
#Protein accessions fetched from protein database                               5           0           5         2.3  BCAR4/BCAR4.acc
##
#Number-of-proteins: 5
##
## Output files created by this script with brief descriptions listed in log file:  BCAR4/BCAR4.log
## This output printed to stdout written to summary file:                           BCAR4/BCAR4.sum
## All commands executed by this script listed in cmd file:                         BCAR4/BCAR4.cmd
## All output files created in directory:                                           ./BCAR4/
## Total seconds elapsed:                                                           2.3
##
##[ok]
#
# For single token output that gives number of protein accessions and
# cleans up all created files:
#
# $ dnaorg_fetch_cds_wrapper.pl -d tmp -np -f BCAR4 | grep ^N | awk '{ print $2 }'; rm -rf tmp 
# 5
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
# Last updated: EPN, Mon Mar 30 10:30:09 2015
#############################################
# Log of changes
#############################################
# EPN, Mon Apr 13 10:04:16 2015
#
# - Modification to Apr 10 update, dead
#   and withdrawn accessions are also 
#   removed in addition to 'suppressed' ones.
#
#  xref: /panfs/pan1/dnaorg/programs/15_0310_dnaorg_fetch_cds_wrapper/00LOG.txt
#
# previous version: ./bkups/15_0413-1-before-update/
# updated  version: ./bkups/15_0413-2-after-update/
#############################################
# EPN, Fri Apr 10 16:05:52 2015
#
# - Added a new step that removes suppressed 
#   (according to idstat) accessions early
#   in the script (immediately after fetching
#   the accessions).
#
#  xref: /panfs/pan1/dnaorg/programs/15_0310_dnaorg_fetch_cds_wrapper/00LOG.txt
#
# previous version: ./bkups/15_0410-1-before-update/
# updated  version: ./bkups/15_0410-2-after-update/
##############################################
# EPN, Tue Apr  7 15:54:41 2015
#
# - Added the -alist option for inputting a list of accessions to
#   fetch CDS for, instead of a gene symbol
# - If zero accessions are fetched, script now exits without error,
#   previously it exited in error in this situation.
# - codon_start now fetched and added to fetched sequence names
# - incomplete cds information '<' and '>' before first and final 
#   coordinate of CDS now preserved in names in output sequence file
#
# previous version: ./bkups/15_0407-1-before-update/
# updated  version: ./bkups/15_0407-2-after-update/
