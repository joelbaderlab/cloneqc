#!/usr/bin/perl -w
use strict;
use lib "/usr/local/bioperl-1.4";
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Root::Exception;
use Error qw(:try);
use Table;

# get and install clustalw (easy) ftp.ebi.ac.uk/pub/software/unix/clustalw
# install io_lib from staden package v. 1.8.12 or earlier
# copy os.h and config.h to the /usr/local/include/io_lib directory
# in os.h, change quotes around config.h to <..> or vice versa 
# instructions for this hack are in bioperl-ext
# install bioperl and bioperl-ext
# bioperl has a built-in wrapper for clustalw http://doc.bioperl.org/bioperl-run/Bio/Tools/Run/Alignment/Clustalw.html

# direction = forward / reverse read
# revcomflag = 0 (no reverse complement) / 1 (reverse complement)

my $DATADIR = "";
my $FIRSTFILE = "";
my $WINDOW = 10;    # window for quality scores, take min of scores +/- WINDOW


my @CLPARAMS = ('ktuple' => 4, 'OUTFILE' => 'clustalw.tmp');

# format = abi or seq (seq means seq and qual files)
# datadir = directory where files reside, and also where output goes
# firstfile = first file to process, optional
sub getParams {
  my ($format, $datadir, $firstfile) = ("", "", "");
  my $nparams = scalar(@ARGV);
  if ( ($nparams < 2) || ($nparams > 3) ) {
    die "Usage: $0 abi|seq <datadirectory> [firstfile]\n";
  }
  $format = $ARGV[0];
  $datadir = $ARGV[1];
  ($datadir =~ /\/$/) || ($datadir .= "/");
  if ($nparams >= 3) {
    $firstfile = $ARGV[2];
  }
  return($format, $datadir, $firstfile);
}

sub getSeqHash {
  my ($file) = @_;
  my $io = Bio::SeqIO->new(-file => "< $file",
			   -format => "fasta");
  my %ret = ();
  while (my $seqobj = $io->next_seq() ) {
    my $key = $seqobj->primary_id();
    (!exists($ret{$key})) || die "Error: reused desired id $key\n";
    $ret{$key} = $seqobj;
  }
  print scalar(keys %ret) . " desired sequences read from $file\n";
  my $ret = \%ret;
  return($ret);
}

sub slidingWindow {
  my ($old, $window) = @_;
  my @new = ();
  my $n = scalar(@{$old});
  for (my $i = 0; $i < $n; ++$i) {
    $new[ $i ] = 100;
    for (my $j = $i - $window; $j <= $i + $window; ++$j) {
      if (($j < 0) || ($j >= $n)) { next; }
      my $val = $old->[$j];
      if ($val < $new[ $i ]) { $new[ $i ] = $val; }
    }
  }
  return(\@new);
}

sub parseStadenExp {
  my ($str) = @_;
  my @qual = ();
  my @line = split(/\n/, $str);
  while (my $rec = shift(@line)) {
    if ($rec =~ /^AV/) {
      unshift(@line, $rec);
      last;
    }
  }
  while (my $rec = shift(@line)) {
    if ($rec !~ /^AV/) {
      unshift(@line, $rec);
      last;
    }
    my @tok = split(/\s+/, $rec);
    my $avstr = shift(@tok);
    push(@qual, @tok);
  }
  while (my $rec = shift(@line)) {
    if ($rec =~ /^SQ/) {
      last;
    }
  }
  my $seq = "";
  while (my $rec = shift(@line)) {
    if ($rec =~ /^\/\//) { last; }
    $rec =~ s/\s+//g;
    $seq .= $rec;
  }

  my $newqual = slidingWindow(\@qual, $WINDOW);

  return($seq, $newqual);
}

# read an abi file and convert it to Bio::Seq object
# carry along quality scores
sub abi2bioseq {
    my ($file) = @_;
    my $cmd = "convert_trace abi exp < $file";
    print "$cmd\n";
    my $str = `$cmd`;
    # $str =~ s/\s+//g;
    my ($seq, $qual) = parseStadenExp($str);
    (length($seq) == scalar(@$qual)) || die "abi2bioseq error in $file length mismatch\n";
    # print join(" ", "seq", length($seq), ":" , $seq) . "\n";
    # print join(" ", "qual", scalar(@$qual), ":", @$qual) . "\n";
    my @base = split("/", $file);
    my $base = pop(@base);
    if (length($seq) == 0) {
      print "Warning: empty sequence\n";
      $seq = "NNNNN";
      $qual = [ 0, 0, 0, 0, 0 ];
    }
    my $bioseq = Bio::Seq->new(-display_id => $base, -seq => $seq, 
			       -alphabet => 'dna' );
    return($bioseq, $qual);
}

# clone -> read (F/R) -> seq , qual
sub getCloneFilesByName {
  my %ret = ();
  foreach my $clone qw(plate1_A02 plate1_A03) {
    foreach my $dir qw(F R) {
      my $file = $DATADIR . $clone . "_${dir}.ab1";
      if (! -e $file) {
	print "Warning: skipping $file\n";
	next;
      }
      $ret{$clone}{$dir}{file} = $file;
    }
  }
  my $nclone = scalar(keys %ret);
  print "$nclone clones\n";
  my $ret = \%ret;
  return($ret);
}

sub getCloneFilesByDir {
  my ($dirname, $firstfile) = @_;
  my %ret = ();
  opendir(DIR, $dirname) || die "could not open dir $dirname: $!";
  print "opened $dirname ...\n";
  my @files = readdir(DIR);
  closedir(DIR);
  foreach my $file (@files) {
    print "file: $file\n";
    ($file =~ /^(.*)_([FR])\.ab1$/) || next;
    print "*** matched $file ***\n";
    if ($file lt $firstfile) {
      print "... waiting for $firstfile\n";
      next;
    }
    my $fullpath = $dirname . $file;
    my $clone = $1;
    my $read = $2;
    $ret{$clone}{$read}{abifile} = $fullpath;
  }
  my $nclone = scalar(keys %ret);
  print "$nclone clones\n";
  my $ret = \%ret;
  return($ret);
}

sub getCloneSeqs {
  my ($clonehash) = @_;
  foreach my $clone (sort keys %$clonehash) {
    foreach my $read (sort keys %{$clonehash->{$clone}}) {
      my $file = $clonehash->{$clone}{$read}{abifile};
      if (! -e $file) {
	print "Warning: skipping $file\n";
	next;
      }
      my ($seqobj, $qual) = abi2bioseq($file);
      @{$clonehash->{$clone}{$read}}{qw(seqobj qual)} = ($seqobj, $qual);
    }
  }
  return($clonehash);
}


sub getCloneSeqQual {
  my ($dirname, $firstfile) = @_;
  my %ret = ();
  opendir(DIR, $dirname) || die "could not open dir $dirname: $!";
  print "opened $dirname ...\n";
  my @files = readdir(DIR);
  closedir(DIR);
  foreach my $file (@files) {
    print "file: $file\n";
    ($file =~ /^(.*)_([FR])\.trimmed\.(seq)$/) 
      || ($file =~ /^(.*)_([FR])\.trimmed\.(qual)$/) || next;
    print "*** matched $file ***\n";
    if ($file lt $firstfile) {
      print "... waiting for $firstfile\n";
      next;
    }
    my $fullpath = $dirname . $file;
    my $clone = $1;
    my $read = $2;
    my $type = $3;
    $type .= "file";
    $ret{$clone}{$read}{$type} = $fullpath;
  }

  # now process the clones that have a seq and a qual
  foreach my $clone (sort keys %ret) {
    foreach my $read (sort keys %{$ret{$clone}}) {
      if (! (exists($ret{$clone}{$read}{seqfile}) &&
	     exists($ret{$clone}{$read}{qualfile})) ) {
	print "dropping $clone $read missing seq or qual\n";
	delete($ret{$clone}{$read});
      }
      my $io = Bio::SeqIO->new(-file => "< $ret{$clone}{$read}{seqfile}",
			       -format => "fasta");
      my $seqobj = $io->next_seq();
      my @qual = ();
      my $qualfile = $ret{$clone}{$read}{qualfile};
      open(IN, "< $qualfile");
      my $header = <IN>;
      ($header =~ /^>/) || die "bad header for $qualfile\n";
      while (<IN>) {
	chomp;
	my @tok = split;
	push(@qual, @tok);
      }
      close(IN);
      (scalar(@qual) == length($seqobj->seq())) ||
	die "mismatch in length for $clone $read\n";
      my $newqual = slidingWindow(\@qual, $WINDOW);
      @{$ret{$clone}{$read}}{qw(seqobj qual)} = ($seqobj, $newqual);
    }
  }
  my $nclone = scalar(keys %ret);
  print "$nclone clones\n";
  my $ret = \%ret;
  return($ret);
}


sub seqobj2str {
  my ($seqobj) = @_;
  my $str = ">" . $seqobj->display_id() . "\n" . $seqobj->seq() . "\n";
  return($str);
}

sub printDesiredInfo {
  my ($tbl) = @_;
  print "\n*** Desired Sequences ***\n";
  foreach my $id (sort keys %$tbl) {
    print seqobj2str($tbl->{$id});
  }
}

sub printCloneInfo {
  my ($tbl) = @_;
  print "\n*** Clones ***\n";
  foreach my $id (sort keys %$tbl) {
    print "\n";
    foreach my $read (sort keys %{$tbl->{$id}}) {
      foreach my $file (qw(abifile seqfile qualfile)) {
	exists($tbl->{$id}{$read}{$file}) || next;
	print "$id $read $file $tbl->{$id}{$read}{$file}\n";
      }
    }
  }
}

sub getAln {
  my ($seqs) = @_;
  my $base = "tmp_clustalw";
  my $fasta = $DATADIR . $base . ".fasta";
  my $out = $DATADIR . $base . ".out";
  my $aln = $DATADIR . $base . ".aln";
  my $tmp = Bio::SeqIO->new(-file => "> $fasta" , -format => "fasta");
  foreach my $seq (@$seqs) {
    $tmp->write_seq($seq);
  }
  my $cmd = "clustalw $fasta > $out";
  # print "... $cmd\n";
  system("$cmd");
  my $alnio = Bio::AlignIO->new(-file => "< $aln", -format => "clustalw");
  my $alnobj = $alnio->next_aln();
  my $pctid = $alnobj->percentage_identity();
  return($alnobj, $pctid);
}

sub getBestMatch {
  my ($desiredHash, $seqobj, $guess) = @_;
  my $revobj = $seqobj->revcom();
  my $bestpctid = -1;
  my $bestkey = "";
  my $bestrc = "";
  my @pair = ( $seqobj , $revobj );
  my @keys = sort keys %$desiredHash;
  if (defined($guess)) { unshift(@keys, $guess); }
  foreach my $key (@keys) {
    if (!exists($desiredHash->{$key})) { next; }
    foreach my $revcom (0, 1) {
      my ($alnobj, $pctid) = getAln([$desiredHash->{$key}, $pair[$revcom]]);
      if ($pctid > $bestpctid) {
	$bestpctid = $pctid;
	$bestkey = $key;
	$bestrc = $revcom;
      }
    }
    if ($bestpctid > 95) { last; }
  }
  return($bestkey, $bestrc, $bestpctid);
}


sub findMutations {
  my ($alnobj, $seqs, $quals) = @_;

  # results table, column qc
  # PK = column
  # columns
  # start, end = start and end locations in the desired sequence
  # seq_1, seq_2, seq_3 = sequence data
  # qual_1, qual_2, qual_3 = quality scores
  # extra = extra sequence data and quality scores
  # qcval = final qc status for this mutation
  # names = desired sequence and reads
  my %colqc = ();

  # sequences in the alignment are 1 = desired, 2 = read 1, 3 = read 2, ...
  # sequences in seqs are 0 = desired, 1 = read 1, 2 = read 2, ...
  # quality scores are 0 = read 1, 1 = read 2, ...
  # renumber so that all the numbering follows the alignment
  my $nseq = $alnobj->no_sequences();
  my @alnseq = ();
  my @origseq = ();
  my @qual = ();
  for (my $i = 1; $i <= $nseq; ++$i) {
    $alnseq[$i] = $alnobj->get_seq_by_pos($i);
    $origseq[$i] = $seqs->[$i - 1];
    if ($i == 1) {
      $qual[1] = [ (100) x $origseq[1]->length ];
    } else {
      $qual[$i] = $quals->[$i - 2];
    }
  }
  print "aligned sequences:\n";
  my @id = ();
  for (my $i = 1; $i <= $nseq; ++$i) {
    $id[$i] = $alnseq[$i]->id();
    print join(" ", $i, $id[$i]) . "\n";
  }
  my $length = $alnobj->length();
  my $pctid = $alnobj->percentage_identity();
  print "alignment length $length pctid $pctid\n";

  # match line: * = match, space = mismatch
  # my $matchline = $alnobj->match_line();
  # print "match line: <$matchline>\n";

  # consensus : [ACGT] = match, ? = mismatch
  my $consensus = $alnobj->consensus_string(100);
  # print "consensus : <$consensus>\n";

  # name and length of the desired sequence
  my $name = $id[1];
  my $len1 = $origseq[1]->length();
  print "desired sequence: $name ($len1 nt)\n";
  # print $alnseq[1]->seq() . "\n";
  # print $origseq[1]->seq() . "\n";

  # start and end columns of the desired sequence
  my $col1 = $alnobj->column_from_residue_number($name, 1);
  my $col2 = $alnobj->column_from_residue_number($name, $len1);
  print "col1 $col1 col2 $col2\n";

  my @mutnlist = ( );

  for (my $col = $col1; $col <= $col2; ++$col) {
    my $consch = substr($consensus, $col - 1, 1);
    if ($consch ne "?") { next; }
    my $mesg = "$col $consch";
    my @alnch = ( );
    my @score = ( );
    my @loc = ( );
    for (my $i = 1; $i <= $nseq; ++$i) {
      my $loc = $alnseq[$i]->location_from_column($col);
      my $str = "loc-undef";
      my $ch = "-";
      my $score = 0;
      if (defined($loc)) {
	my ($start, $end) = ($loc->start, $loc->end);
	$str = join("-", $start, $end, $loc->location_type);
	for (my $j = $start; $j <= $end; ++$j) {
	  my $val = exists($qual[$i]->[$j - 1]) ? $qual[$i]->[$j - 1] : 0;
	  $score += $val;
	}
	$score /= ($end - $start + 1);
      }
      $alnch[$i] = $alnseq[$i]->subseq($col, $col);
      $score[$i] = $score;
      $loc[$i] = $loc;
      $mesg .= "\t$str";      
    }

    $colqc{$col} = { };
    if (defined($loc[1])) {
      @{$colqc{$col}}{qw(start end)} = ($loc[1]->start, $loc[1]->end);
    } else {
      @{$colqc{$col}}{qw(start end)} = ("undef", "undef");
    }
    for (my $i = 1; ($i <= $nseq) && ($i <= 3); ++$i) {
      $colqc{$col}{"seq" . $i} = $alnch[$i];
      $colqc{$col}{"qual" . $i} = $score[$i];
      $colqc{$col}{"name" . $i} = $id[$i];
    }
    $colqc{$col}{seq_extra} = $colqc{$col}{qual_extra} = $colqc{$col}{name_extra} = ".";
    if ($nseq > 3) {
      $colqc{$col}{seq_extra} = join(",", @alnch[4..$nseq]);
      $colqc{$col}{qual_extra} = join(",", @score[4..$nseq]);
      $colqc{$col}{name_extra} = join(",", @id[4..$nseq]);
    }

    # now to evaluate
    my $evalstr = "";
    my $qcstr = "CHECK";
    my (@match, @mismatch) = ( );
    my %mutn = ();
    for (my $i = 2; $i <= $nseq; ++$i) {
      if ($alnch[$i] eq $alnch[1]) {
	push(@match, $i);
      } else {
	push(@mismatch, $i);
	++$mutn{$alnch[$i]};
      }
    }
    my $pos1 = $colqc{$col}{start};
    my $ch1 = $alnch[1];
    my $newch = join("|", sort keys %mutn);
    # more information for a check
    # give each character together with quality score
    my $checkch = "[" .
      join( "",
	    map { $alnch[$_] . "(" . $score[$_] . ")" }
	    (2 .. $nseq) ) . "]" ;
    # if there is at least one match, check for a basecalling error
    my $GOODSCORE = 40;
    my $BADSCORE = 30;
    if (@match > 0) {
      my $basecallerror = 0;
      foreach my $i (@match) {
	if ($score[$i] >= $GOODSCORE) {
	  ++$basecallerror;
	}
      }
      foreach my $i (@mismatch) {
	if ($score[$i] >= $BADSCORE) {
	  $basecallerror = 0;
	}
      }
      if ($basecallerror) {
	$evalstr = "BASECALLERROR";
	$qcstr = "PASS";
      } else {
	$evalstr = "CHECK BASECALL";
	$qcstr = "CHECK";
	push(@mutnlist, ["chk", $pos1, $ch1, $checkch] );
      }
    } elsif (scalar(keys %mutn) == 1) { 
      # all the sequences agree, a mutation
      $qcstr = "FAIL";
      if ($alnch[1] eq "-") {
	$evalstr = "INSERTION $newch";
	push(@mutnlist, ["ins", $pos1, $ch1, $newch]);
      } elsif ($newch eq "-") {
	$evalstr = "DELETION $alnch[1]";
	push(@mutnlist, ["del", $pos1, $ch1, $newch]);
      } else {
	$evalstr = "SUBSTITUTION $alnch[1] -> $newch";
	push(@mutnlist, ["sub", $pos1, $ch1, $newch]);
      }
    } else {
      $evalstr = "CHECK $alnch[1] -> [$newch]";
      push(@mutnlist, ["chk", $pos1, $ch1, $checkch]);
      $qcstr = "CHECK";
    }
    $colqc{$col}{evalstr} = $evalstr;
    $colqc{$col}{qcstr} = $qcstr;


    $mesg .=  "\t" .
      " " . join("", @alnch[1..$nseq]) .
      " " . join(",", @score[1..$nseq]) . "\n";
    print $mesg;
    my $mutnliststr = (@mutnlist == 0) ? "NoMutns" :
      join(" ", @{$mutnlist[$#mutnlist]});
    print "$mutnliststr\n";
  }
  my $ret = Table->as(\%colqc);

  # join insertion / deletion runs
  my $lasttype = "";
  my $lastpos = -10;
  my $lastch = "";
  my $lastnew = "";
  my @newmutnlist = ();
  foreach my $rec (@mutnlist) {
    my ($type, $pos, $ch, $newch) = @$rec;
    # if same type and same pos, add the information
    if ( ($type eq $lasttype) && ($pos <= $lastpos + 1) ) {
      $lastpos = $pos;
      $lastch .= $ch;
      $lastnew .= $newch;
    } else {
      # append to the list
      if ($lasttype ne "") {
	push(@newmutnlist, [ $lasttype, $lastpos, $lastch, $lastnew ]);
      }
      # update the memory
      $lasttype = $type;
      $lastpos = $pos;
      $lastch = $ch;
      $lastnew = $newch;
    }
  }
  if ($lasttype ne "") {
    push(@newmutnlist, [ $lasttype, $lastpos, $lastch, $lastnew ]);
  }
  
  # mutantion counts
  my %mutncnt = ();
  my @mutntype = qw(ins del sub chk);
  foreach my $key (@mutntype) { $mutncnt{$key} = 0; }
  foreach my $rec (@newmutnlist) {
    ++$mutncnt{ $rec->[0] };
  }
  my $mutnstr = (@newmutnlist == 0) ? "" :
    join(" ", map { join(":", @$_) } @newmutnlist);
  
  return($ret, \%mutncnt, $mutnstr); 
}

my $GUESS = "";  # global holding the last bbid
sub validateClone {
  my ($desiredHash, $reads) = @_;
  # get the best match for each read
  my %qc = ();
  my %match = ();
  my %cntrc = ();
  my $match_id = undef;
  my $EMPTY = "EMPTY";
  my $emptythresh = 50;
  my @seqobjs = ();
  my @quals = ();
  my @tracenames = ();
  foreach my $rd (sort keys %$reads) {
    print "... read $rd\n";
    my $obj = $reads->{$rd}{seqobj};
    my $seqlen = length($obj->seq);
    my ($match, $revcom, $pctid) = ($EMPTY, 0, 0);
    if ($seqlen >= $emptythresh) {
      ($match, $revcom, $pctid) = getBestMatch($desiredHash, $obj, $GUESS);
    }
    print "read $rd bestseq $match revcom $revcom pctid $pctid\n";
    $match_id = $match;
    $GUESS = $match;
    $match{$match}++;
    $cntrc{$revcom}++;
    my $newobj = $obj;
    if ($revcom) {
      $newobj = $obj->revcom();
      my $newid = $newobj->display_id() . "_revcom";
      $newobj->display_id($newid);
    }
    push(@seqobjs, $newobj );
    push(@quals, ($revcom == 0) ? $reads->{$rd}{qual} :
	 [ reverse(@{$reads->{$rd}{qual}}) ]);
    push(@tracenames, $newobj->display_id());
  }
  # check that all the reads have the same match
  my $matchqc = "PASS";
  if ($match_id eq $EMPTY) {
    $matchqc = "FAIL:EMPTY";
  } elsif (scalar(keys %match) != 1) {
    $matchqc = "FAIL:MULTIPLE";
  }
  $qc{matchqc} = $matchqc;
  $qc{bb_id} = join(" ", sort keys %match);
  # check that there is at least one rc and one non-rc read
  $cntrc{0} += 0;
  $cntrc{1} += 0;
  $qc{revcomqc} = ( ($cntrc{0} > 0) && ($cntrc{1} > 0) ) ? "PASS" : "FAIL";

  for (my $i = 0; ($i < @tracenames) && ($i < 2); ++$i) {
    $qc{"read" . ($i + 1)} = $tracenames[$i];
  }
  $qc{"read_extra"} = (@tracenames > 2) ?
    join(";", @tracenames[2..$#tracenames]) : ".";

  if (($qc{matchqc} ne "PASS") ||
      ($qc{revcomqc} ne "PASS")) { return(\%qc); }

  # align everything
  unshift(@seqobjs, $desiredHash->{$match_id});
  my ($alnobj, $pctid) = getAln([@seqobjs]);
  $qc{pctid} = $pctid;
  $qc{length} = $desiredHash->{$match_id}->length();

  my ($colqc, $mutncnt, $mutnstr) = ( {} , "NA", "");
  my $error = 0;
  try {
    ($colqc, $mutncnt, $mutnstr) = 
      findMutations($alnobj, [ @seqobjs ], [ @quals ]);
  } otherwise {
    print "alignment error, failing this clone\n";
    print "error text:\n-----\n" . join("\n", @_) . "-----\n";
    $qc{mutnqc} = "FAIL:ALIGNCRASH";
    $qc{mutncnt} = "NA";
    $qc{overallqc} = "FAIL";
    $error = 1;
  } finally { };  # need the trailing semicolon!!!
  
  if ($error) { return(\%qc); }
  my %cnt = ();
  foreach my $col (keys %$colqc) {
    ++$cnt{ $colqc->{$col}{qcstr} };
  }
  my $mutnqc = "PASS";
  if (exists($cnt{FAIL})) {
    $mutnqc = "FAIL";
  } elsif (exists($cnt{CHECK})) {
    $mutnqc = "CHECK";
  }
  $qc{mutnqc} = $mutnqc;
  $qc{mutnstr} = $mutnstr;

  my $ntot = 0;
  foreach my $type (keys %$mutncnt) {
    $qc{"n_" . $type} = $mutncnt->{$type};
    $ntot += $mutncnt->{$type};
  }
  $qc{n_tot} = $ntot;

  ++$cnt{ $qc{matchqc} };
  ++$cnt{ $qc{revcomqc} };

  my $qcstr = "PASS";
  if (exists($cnt{FAIL})) {
    $qcstr = "FAIL";
  } elsif (exists($cnt{CHECK})) {
    $qcstr = "CHECK";
  }
  $qc{overallqc} = $qcstr;
  
  return(\%qc);

}
    

# validate a list of clones against a list of desired sequences
# clonecolqc
# PK = (cloneid, column in alignment)
sub validateAll {

  my ($desiredhash, $clonehash) = @_;

  my %cloneqc = ();

  my $nclone = scalar(keys %$clonehash);
  my $ndesired = scalar(keys %$desiredhash);

  print "validating $nclone clones against $ndesired desired sequences\n";

  my @FIELD = 	qw(bb_id 
		   length
		   overallqc
		   mutnqc
		   revcomqc
		   matchqc
		   pctid
		   read1
		   read2
		   read_extra
		   n_ins
		   n_del
		   n_sub
		   n_chk
		   n_tot
		   mutnstr
		   );
  
  my $logfile = $DATADIR . "log.xls";
  open(LOG, "> $logfile");
  print LOG join("\t", "clone", @FIELD) . "\n";
  close(LOG);
  foreach my $clone (sort keys %$clonehash) {
    print "\nworking on clone $clone ...\n";
    my ($qc) = validateClone($desiredhash, $clonehash->{$clone});
    foreach my $f (@FIELD) {
      exists($qc->{$f}) || ($qc->{$f} = "NA");
    }
    $cloneqc{$clone} = $qc;

    open(LOG, ">> $logfile");
    print LOG join("\t", $clone, map { $qc->{$_} } @FIELD) . "\n";
    close(LOG);
  }
  my $cloneqc = Table->as(\%cloneqc);
  return($cloneqc);
}


###
# main
###


# get the data directory and the first file to process
(my $format, $DATADIR, $FIRSTFILE) = getParams();
print "$0 $format $DATADIR $FIRSTFILE\n";

# get the list of desired sequences
#my $desiredfile = "3L.3_23.BAG2007f-BBs.fasta";
my $desiredfile = "3R1.fasta";
#my $desiredfile = "tmp.fasta";
my $desiredhash = getSeqHash($desiredfile);

# get clone files
my $clonehash = {};
if ($format eq "abi") {
  $clonehash = getCloneFilesByDir($DATADIR, $FIRSTFILE);
  getCloneSeqs($clonehash);
} elsif ($format eq "seq") {
  $clonehash = getCloneSeqQual($DATADIR, $FIRSTFILE);
}

printCloneInfo($clonehash);

# printDesiredInfo($desiredhash);

my ($cloneqc) = validateAll($desiredhash, $clonehash);

$cloneqc->write($DATADIR . "cloneqc.xls", 1,
		[qw(bb_id 
		    length
		    overallqc
		    mutnqc
		    revcomqc
		    matchqc
		    pctid
		    read1
		    read2
		    read_extra
		    n_ins
		    n_del
		    n_sub
		    n_chk
		    n_tot
		    mutnstr
		   )]
	       );
