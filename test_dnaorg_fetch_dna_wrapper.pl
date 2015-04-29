use strict;
use warnings;

#my $usage = "perl test_dnaorg_fetch_dna_wrapper.pl <executable-dir> <testfile-dir> <'1' for verbose mode, '0' for quiet mode>\n";
my $usage = "perl test_dnaorg_fetch_dna_wrapper.pl <executable-dir> <'1' for verbose mode, '0' for quiet mode>\n";
if(scalar(@ARGV) != 2) { die $usage; }

#my ($execdir, $testdir, $be_verbose) = @ARGV;
my ($execdir, $be_verbose) = @ARGV;

if($be_verbose ne "1" && $be_verbose ne "0") { die $usage; }

$execdir =~ s/\/$//; # remove trailing '/' if one exists
#$testdir =~ s/\/$//; # remove trailing '/' if one exists

# definitions of variables
my $fail_ntests       =  2; # number of 'fail' tests,    the tests that are     expected to return a failure (non-zero exit status)
my $nofail_ntests     =  5; # number of 'no fail' tests, the tests that are not expected to return a failure (non-zero exit status)
my @fail_descA        = (); # brief descriptions of each fail test, for informative output
my @nofail_descA      = (); # brief descriptions of each no fail test, for informative output
my @fail_argA         = (); # name of single command line argument for each fail test
my @nofail_argA       = (); # name of single command line argument for each nofail test
my @fail_dirA         = (); # name of output directories for each fail test
my @nofail_dirA       = (); # name of output directories for each no fail test
my @fail_inputA       = (); # full input (options plus command line arguments) for each fail test
my @nofail_inputA     = (); # full input (options plus command line arguments) for each no fail test
#my @nofail_outputA    = (); # output files for each no fail test
my @fail_nerrlinesA   = (); # number of expected stderr output lines for fail tests
my @fail_stdout_okayA = (); # '1' if it's okay for some stdout to be printed for each fail test
my @nofail_noutlinesA = (); # number of expected stdout output lines for no fail tests
my ($f, $n);                # counter over fail and no fail tests, respectively
my $cmd;                    # a command to run with system()

# number of expected lines to stderr if script fails
my $fail_exp_nerrlines = 64;

# number of expected lines for successful execution of script
my $nofail_exp_noutlines_protein    = 36;
my $nofail_exp_noutlines_nucleotide = 33;
my $nofail_exp_noutlines_num        = 28;
my $nofail_extra_lines_skip_synonyms = 5;

# check executable files exist
my $wrapper  = $execdir . "/dnaorg_fetch_dna_wrapper.pl";
if(! -e $wrapper) { die "ERROR $wrapper does not exist in $execdir"; }

# Description of tests, included here for reference
# tests that are expected to fail (testing the error checking code in the program)
$fail_descA[0]  = "wrong number cmdline args";
$fail_descA[1]  = "symbol listed in fail file";
# make sure we have descriptions for all fail tests
for($f = 0; $f < $fail_ntests; $f++) { 
  if((! defined $fail_descA[$f]) || ($fail_descA[$f] eq "")) { 
    die "ERROR (in test script) fail test " . ($f+1) . " does not have a description";
  }
}

# tests that are expected to not fail
$nofail_descA[0]  = "protein symbol mode (primary symbol and alias)";
$nofail_descA[1]  = "nucleotide symbol mode";
$nofail_descA[2]  = "list mode";
$nofail_descA[3]  = "num mode";
$nofail_descA[4]  = "-nosyn option, nucleotide symbol mode";
# make sure we have descriptions for all no fail tests
for($n = 0; $n < $nofail_ntests; $n++) { 
  if((! defined $nofail_descA[$n]) || ($nofail_descA[$n] eq "")) { 
    die "ERROR (in test script) nofail test " . ($n+1) . " does not have a description";
  }
}

# set input command line arguments and options:
$fail_argA[0]   = "";
$fail_argA[1]   = "SECX";

$fail_dirA[0]   = "fail0";
$fail_dirA[1]   = "fail1";

$fail_inputA[0] = "-f -d $fail_dirA[0] $fail_argA[0]";
$fail_inputA[1] = "-f -d $fail_dirA[1] $fail_argA[1]";

# set input command line arguments and options:
$nofail_argA[0] = "USP17L5";
$nofail_argA[1] = "SNORA71D";
$nofail_argA[2] = "samp5.list";
$nofail_argA[3] = "SNORA71D";
$nofail_argA[4] = "SNORA71D";

$nofail_dirA[0] = "nofail0";
$nofail_dirA[1] = "nofail1";
$nofail_dirA[2] = "nofail2";
$nofail_dirA[3] = "nofail3";
$nofail_dirA[4] = "nofail4";

$nofail_inputA[0] = "-f -d $nofail_dirA[0] $nofail_argA[0]";
$nofail_inputA[1] = "-f -d $nofail_dirA[1] $nofail_argA[1]";
$nofail_inputA[2] = "-f -plist -d $nofail_dirA[2] $nofail_argA[2]";
$nofail_inputA[3] = "-f -num -d $nofail_dirA[3] $nofail_argA[3]";
$nofail_inputA[4] = "-f -nosyn -d $nofail_dirA[4] $nofail_argA[4]";

## check that required inputs/outputs exist
#for($f = 0; $f < $fail_ntests; $f++) { 
#  $fail_inputA[$f]  = $testdir . "/fail." . ($f+1) . ".in"; # off-by-one
#  if(! -e $fail_inputA[$f]) { die "ERROR " . $fail_inputA[$f] . "(input file) does not exist"; }
#  $fail_descA[$f] .= " [input file: " . $fail_inputA[$f] . "]";
#}
#for($n = 0; $n < $nofail_ntests; $n++) { 
#  $nofail_inputA[$n]   = $testdir . "/nofail." . ($n+1) . ".in";  # off-by-one
#  $nofail_outputA[$n]  = $testdir . "/nofail." . ($n+1) . ".out"; # off-by-one
#  if(! -e $nofail_inputA[$n])  { die "ERROR " . $nofail_inputA[$n] . " (input file) does not exist"; } 
#  if(! -s $nofail_outputA[$n]) { die "ERROR " . $nofail_outputA[$n] . " (output file) does not exist"; }
#  $nofail_descA[$n] .= " [input file: " . $nofail_inputA[$n] . ", output_file: " . $nofail_outputA[$n] . "]";
#}

# set default number of stderr lines for fail tests,
# and whether or not it's okay for each fail test to
# create some stdout, or not
for($f = 0; $f < $fail_ntests; $f++) { 
  $fail_stdout_okayA[$f] = 0; # not okay, by default
}
$fail_nerrlinesA[0] = $fail_exp_nerrlines;
$fail_nerrlinesA[1] = 2;
# set exceptions: 
#$fail_stdout_okayA[10] = 0; # this test should print some to stdout before failing

$nofail_noutlinesA[0]  = $nofail_exp_noutlines_protein + $nofail_extra_lines_skip_synonyms;
$nofail_noutlinesA[1]  = $nofail_exp_noutlines_nucleotide;
$nofail_noutlinesA[2]  = $nofail_exp_noutlines_protein;
$nofail_noutlinesA[3]  = $nofail_exp_noutlines_num;
$nofail_noutlinesA[4]  = $nofail_exp_noutlines_nucleotide + 2; # -nosyn gives 2 extra lines of output

# Run fail tests (these should create stderr output)
for($f = 0; $f < $fail_ntests; $f++) { 
  my $stderr = "fail." . ($f+1) . ".stderr";
  my $stdout = "fail." . ($f+1) . ".stdout";

  # run command and ensure it has expected non-zero exit status (die if it does not)
  $cmd = "$wrapper $fail_inputA[$f] 2> $stderr > $stdout";
  RunCommand($cmd, 1, $fail_descA[$f], $be_verbose); # '1' says: failure is expected

  # check that command had expected number of lines in stderr and 
  # either did or did not have stdout output
  CheckNumLinesInFile($stderr, $fail_nerrlinesA[$f], $fail_descA[$f]);
  if(! $fail_stdout_okayA[$f]) { 
    CheckNumLinesInFile($stdout, 0, $fail_descA[$f]); # 0 says that $stdout should either not exist or exist but with 0 lines
  }

  # clean up
  for my $file ($stderr, $stdout) { 
    if(-e $file) { unlink $file; }
  }
  $cmd = "rm -rf $fail_dirA[$f]";
  RunCommand($cmd, 0, "removing $fail_dirA[$f]", 0);
}

# Run no fail tests (these should not fail and should give output)
for($n = 0; $n < $nofail_ntests; $n++) { 
  my $stderr = "nofail." . ($n+1) . ".stderr"; # off-by-one, to match input file name
  my $stdout = "nofail." . ($n+1) . ".stdout"; # off-by-one, to match input file name

  # run command and ensure it has expected non-zero exit status (die if it does not)
  $cmd = "$wrapper $nofail_inputA[$n] 2> $stderr > $stdout";
  RunCommand($cmd, 0, $nofail_descA[$n], $be_verbose); # '0' says: failure is NOT expected

  # check that command had expected stderr and stdout output
  CheckNumLinesInFile($stderr, 0,                      $nofail_descA[$n]); # 0 says that $stderr should either not exist or exist but with 0 lines 
  CheckNumLinesInFile($stdout, $nofail_noutlinesA[$n], $nofail_descA[$n]);

  # make sure output file is what we expect
  #CompareFilesWithDiff($stdout, $nofail_outputA[$n], $nofail_descA[$n]);

  for my $file ($stderr, $stdout) { 
    if(-e $file) { unlink $file; }
  }
  $cmd = "rm -rf $nofail_dirA[$n]";
  RunCommand($cmd, 0, "removing $nofail_dirA[$n]", 0);
}
print "PASS [$fail_ntests tests failed as expected and $nofail_ntests tests succeeded as expected]\n";

exit 0;  


# Subroutine: RunCommand()
# Args:       $cmd:            command to run, with a "system" command;
#             $expect_failure: '1' if $cmd should fail (return non-zero exit status)
#             $desc:           description to print if command return expected exit status
#             $be_verbose:     '1' to output command before we run it, '0' not to
# Dies:       if $cmd fails and $expect_failure == 0
#             OR
#             if $cmd does not fail and $expect_failure == 1

sub RunCommand {
  my $sub_name = "RunCommand()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cmd, $expect_failure, $desc, $be_verbose) = @_;

  if($be_verbose) { printf("Test: $desc\n"); }   
  #print Running cmd: $cmd (expect failure: $expect_failure)\n", $desc);
  system($cmd);
  if($expect_failure) { 
    if($? == 0) { die "ERROR command with desc $desc did not fail.\nActual command:\n$cmd\n"; }
  }
  else { 
    if($? != 0) { die "ERROR command with desc $desc failed.\nActual command:\n$cmd\n"; }
  }
  return;
}

# Subroutine: CheckNumLinesInFile()
# Args:       $file:       command to run, with a "system" command;
#             $nlines_exp: number of lines expected
#                          if '0',  file can either not exist or exist but be empty. 
#                          if '-1', file should not exist
#             $desc:       description to print if command return expected exit status
# Dies:       if file does not have expected number of lines

sub CheckNumLinesInFile {
  my $sub_name = "CheckNumLinesInFile()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($file, $nlines_exp, $desc) = @_;
  if($nlines_exp < -1) { die "ERROR $sub_name, negative number of expected lines"; }

  my $die_msg = "This means the following test failed: $desc\n";

  if($nlines_exp == -1) { 
    if(-e $file) { die "ERROR file $file exists (it should not exist).\n$die_msg"; }
  }
  elsif($nlines_exp == 0) { 
    if(-s $file) { die "ERROR file $file exists and is non-empty (it should either not exist or be empty).\n$die_msg"; }
  }
  else { # $nlines_exp > 0
    if(! -e $file) { die "ERROR file $file does not exist but it should.\n$die_msg"; }
    my $nlines = `wc -l $file | awk '{ print \$1 }'`;
    chomp $nlines;
    if($? != 0) { die "ERROR problem getting number of lines in file $file.\n$die_msg"; }
    if($nlines != $nlines_exp) { die "ERROR, wrong number of lines ($nlines != $nlines_exp) in file $file.\n$die_msg"; }
  }
  return;
}

# Subroutine: CompareFilesWithDiff()
# Args:       $file1: first file
#             $file2: second file
#             $desc:  description to print if files differ
# Dies:       if 'diff' says $file1 and $file2 are different

sub CompareFilesWithDiff {
  my $sub_name = "CompareFilesWithDiff()";
  my $nargs_exp = 3;

  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($file1, $file2, $desc) = @_;

  my $die_msg = "This means the following test failed: $desc\n";

  my $diff_output = `diff $file1 $file2`;

  if($? != 0 || $diff_output =~ m/\w/) { 
    die "ERROR diff says $file1 and $file2 differ.\n$die_msg" . $diff_output; 
  }
  return;
}
  

