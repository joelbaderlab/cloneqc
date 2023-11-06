#!/usr/bin/perl -w
package Table;

use strict;
use Exporter;

our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw( ) ] );
our @EXPORT_OK = (@{$EXPORT_TAGS{all}});
our @EXPORT = qw();

# undef would be preferable for missing data
# but it's a hassle to wrap prints and comparisons in defined()'s
# hence ... the NA string
my $NA = "NA";

# bless a hash as a table
sub as {
    my ($class, $obj) = @_;
    my $reftype = ref($obj);
    print "blessing a $reftype as a $class\n";
    my $self = bless $obj, $class;
    return($self);
}

sub new {
    my ($class, $file, $nkey) = @_;
    (defined($nkey) && ($nkey > 0)) || die "Error: need nkey > 0 in $class new\n";
    defined($file) || die "Error: file name undefined\n";
    open(IN, "< $file") || die "could not read $file\n";
    my $header = <IN>;
    while ($header =~ /^\#/) {
      $header = <IN>;
    }
    # chomp($header);
    $header =~ s/\r?\n//;
    my @col = split /\t/, $header;
    my $nextra = scalar(@col) - $nkey;
    ($nextra >= 0) || die "too few fields to fill PK: $_\n";
    my @index = (0..($nkey-1));
    my @field = @col[$nkey..$#col]; 
    print "reading $file as a $class\n";
    print "primary key: " . join(" ", @col[@index]) . "\n";
    print "data fields: " . join(" ", @field) . "\n";

    
    my %ret = ();
    my $table = \%ret;
    my %seen = ();
    while (<IN>) {
	# chomp;
	$_ =~ s/\r?\n//;
	my @tok = split /\t/;
	# (@tok == @col) || print "WARNING: bad token count: $_\n";
	if (@tok < @col) {
	    my $cnt = scalar(@col) - scalar(@tok);
	    for (my $i = 0; $i < $cnt; ++$i) {
		push(@tok, "");
	    }
	}
	my @key = @tok[@index];
	my $pk = join("\t", @key);
	!exists($seen{$pk}) || die "repeated PK: $pk\n";
	++$seen{$pk};
	my $ref = \%ret;
	foreach my $k (@key) {
	    exists($ref->{$k}) || ($ref->{$k} = {});
	    $ref = $ref->{$k};
	}
	if ($nextra > 0) {
	    @{$ref}{@field} = @tok[$nkey..$#col];
	} else {
	    ++$ref;
	}
    }
    close(IN);
    print scalar(keys %seen) . " PK's in $file\n";
    my $self = bless \%ret, $class;
    return($self);
}    
    

sub writeRows {
  my ($table, $fh, $col, $nkey, $list, $sort) = @_;
  $nkey >= 0 || die "bad nkey $nkey\n";
  if ($nkey == 0) {
    # print $fh join("\t", @{$list}) . "\n";
    defined($table) || die "table not defined for list\n";
    defined($col) || die "column not defined\n";
    foreach my $c (@{$col}) {
      exists( $table->{$c} ) || die "$c does not exist\n";
      if (!defined($table->{$c} )) {
	print "WARNING: $c exists but not defined\n";
	print "keys: " . join(" ", @{$list}) . "\n";
	$table->{$c} = "";
      }
    }
    my $valuestr = (@{$col} > 0) ? join("\t", @{$table}{@{$col}}) :
      $table;
    print $fh join("\t", @{$list}, $valuestr) . "\n";
  } else {
    my @keys = keys %{$table};
    my $numeric = 1;
    foreach my $k (@keys) {
      if ($k !~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/) {
	$numeric = 0;
	last;
      }
    }
    if ($numeric) {
      @keys = sort {$a <=> $b} @keys;
    } else {
      @keys = sort @keys;
    }

    #if ($sort eq "alpha") {
    # @keys = sort @keys;
    #} elsif ($sort eq "numeric") {
    #@keys = sort {$a <=> $b} @keys;
    #} else {
    #  die "unknown sort mode: $sort (should be 'alpha' or 'numeric')\n";
    #}

    foreach my $key (@keys) {
      push(@{$list}, $key);
      writeRows($table->{$key}, $fh, $col, $nkey - 1, $list, $sort);
      pop(@{$list});
    }
  }
}

sub getFields {
    my ($self, $nkey) = @_;
    my $ptr = $self;
    for (my $i = 0; $i < $nkey; ++$i) {
	ref($ptr) || die "not a ref, nkey = $nkey, i = $i\n";
	my ($key) = keys %{$ptr};
	$ptr = $ptr->{$key};
    }
    my @field = ();
    if (ref($ptr)) {
	@field = sort keys %{$ptr};
    }
    my $ret = [ @field ];
    return($ret);
}

sub write {
    my ($self, $file, $nkey, $order, $sort) = @_;
    defined($nkey) || die "need to specify number of keys\n";
    defined($sort) || ($sort = "alpha");
    open(OUT, "> $file") || die "could not open $file\n";
    print "writing table to $file\n";
    # get the fields
    my $col = defined($order) ? $order : getFields($self, $nkey);
    my @col = (@{$col} > 0) ? @{$col} : "value";
    my @key = map { "key" . $_ } (1..$nkey);
    print OUT join("\t", @key, @col). "\n";
    writeRows($self, *OUT, $col, $nkey, [ ], $sort);
    close(OUT);
}

sub writeRaggedRows {
    my ($table, $fh, $nkey, $list) = @_;
    $nkey >= 0 || die "bad nkey $nkey\n";
    if ($nkey == 0) {
	my $str = join(";", map { $_ . "=" . $table->{$_} }
		       sort keys %{$table});
	print $fh join("\t", @{$list}, $str) . "\n";
    } else {
	foreach my $key (sort keys %{$table}) {
	    push(@{$list}, $key);
	    writeRaggedRows($table->{$key}, $fh, $nkey - 1, $list);
	    pop(@{$list});
	}
    }
}

sub writeRagged {
    my ($self, $file, $nkey) = @_;
    open(OUT, "> $file") || die "could not open $file\n";
    print "writing table to $file\n";
    # get the fields
    my @key = map { "key" . $_ } (1..$nkey);
    print OUT join("\t", @key, "properties"). "\n";
    writeRaggedRows($self, *OUT, $nkey, [ ]);
    close(OUT);
}

sub sumRows {
    my ($self, $ret, $nkey) = @_;
    if ($nkey == 0) {
	foreach my $f (keys %{$self}) {
	    my $x = $self->{$f};
	    $ret->{n}{$f} += 1;
	    $ret->{avg}{$f} += $x;
	    $ret->{stdev}{$f} += $x * $x;
	}
    } else {
	foreach my $k (keys %{$self}) {
	    sumRows($self->{$k}, $ret, $nkey - 1);
	}
    }
}
	    

sub collapseRows {
    my ($self, $ret, $oldkey, $newkey) = @_;
    if ($newkey > 0) {
	foreach my $k (keys %{$self}) {
	    $ret->{$k} = {};
	    collapseRows($self->{$k}, $ret->{$k}, $oldkey - 1, $newkey - 1);
	}
    } else {
	$ret->{n} = {};
	$ret->{avg} = {};
	$ret->{stdev} = {};
	sumRows($self, $ret, $oldkey);
	foreach my $f (keys %{$ret->{n}}) {
	    my ($sum0, $sum1, $sum2) =
		($ret->{n}{$f}, $ret->{avg}{$f}, $ret->{stdev}{$f});
	    my $avg = $sum1 / $sum0;
	    my $stdev = 0;
	    if ($sum0 > 1) {
		$stdev = $sum2 - ($sum1 * $sum1 / $sum0);
		$stdev = sqrt($stdev / $sum0);
	    }
	    $ret->{avg}{$f} = $avg;
	    $ret->{stdev}{$f} = $stdev;
	}
    }
}

# collapse the number of keys from old to new
sub collapse {
    my ($self, $oldkey, $newkey) = @_;
    my $ret = {};
    bless $ret, ref($self);
    collapseRows($self, $ret, $oldkey, $newkey);
    return($ret);
}

sub leftOuterJoinRows {
    my ($master, $slave, $nkey, $col) = @_;
    if ($nkey == 0) {
	my @col = ();
	if (defined($col)) {
	    @col = @{$col};
	} else {
	    @col = keys %{$slave};
	}
	my @val = ();
	foreach my $c (@col) {
	    !exists($master->{$c}) || die "overwriting column $c in master\n";
	    push(@val, exists($slave->{$c}) ? $slave->{$c} : "NA");
	}
	@{$master}{@col} = @val;
    } else {
	foreach my $k (keys %{$master}) {
	    if (exists($slave->{$k})) {
		leftOuterJoinRows($master->{$k}, $slave->{$k}, $nkey - 1, $col);
	    } else {
		fill($master->{$k}, $nkey - 1, $col, $NA);
	    }
	}
    }
}

sub leftOuterJoin {
    my ($self, $slave, $nkey, $col) = @_;
    my $newcol = defined($col) ? $col : getFields($slave, $nkey);
    leftOuterJoinRows($self, $slave, $nkey, $newcol);
}

sub subsetRows {
    my ($self, $subset, $nkey) = @_;
    if ($nkey == 1) {
	foreach my $k (keys %{$self}) {
	    exists($subset->{$k}) || delete($self->{$k});
	}
    } else {
	foreach my $k (keys %{$self}) {
	    my $delete = 1;
	    if (exists($subset->{$k})) {
		subsetRows($self->{$k}, $subset->{$k}, $nkey - 1);
		$delete = scalar(keys %{$self->{$k}}) > 0;
	    }
	    if ($delete) { delete($self->{$k}); }
	}
    }
}

sub subset {
    my ($self, $subset, $nkey) = @_;
    subsetRows($self, $subset, $nkey);
}

sub fill {
    my ($self, $nkey, $col, $val) = @_;
    if ($nkey == 0) {
	foreach my $c (@{$col}) {
	    exists($self->{$c}) || ($self->{$c} = $val);
	}
    } else {
	foreach my $k (keys %{$self}) {
	    fill($self->{$k}, $nkey-1, $col, $val);
	}
    }
}

sub tableTest {
    my $table = tableRead("times.txt", [ qw(a b times)], 2);
    my $plus = tableRead("plus.txt", [ qw(a b plus minus) ], 2);
    innerjoin($table, $plus, 2, [ qw(plus minus) ]);
    tableWrite($table, "junk.txt", [ qw(a b plus times minus) ], 2);
    tableWrite($plus, "junk1.txt", [ qw(a b minus plus) ], 2);
}

1;
