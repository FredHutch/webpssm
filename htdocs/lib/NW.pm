# needleman-wunsch for V3 PSSM pre-processing
# MAJ 12/10/03

package NW;
use Exporter;
@ISA = ('Exporter');
@EXPORT = ('nw','loadpam', 'mkseq', 'pralign');

# 
# slide numbers refer to powerpoint presentation
# describing the needleman wunsch algorithm
# "stupid poster.ppt"
# a very brief description is on slide 7
#


# codes for incremental path vectors
# 
# slide 8, right panel shows the path vector possibilities;
# more than one vector can occupy a cell, so the path vectors
# are coded in binary, such that ORing them will allow representation
# of multiple vectors

$II = 4;
$IZ = 2;
$ZI = 1;

#to prevent addition of gaps to sequence $s, set $first_seq_fixed:
# (not required in PSSM use)
$first_seq_fixed = 0;

#default aa align
# set NW parameters; these could be tweaked
$unit = 1; # number of aa's/element to align
$gxpen = 5; # gap extension penalty ; guess from BioEdit
$gpen = 20; # gap penalty ; guess from BioEdit-

# this function scores a single site pairwise comparison
sub wpam { # use PAM matrix to compare aa pairs
    my ($x, $y) = @_;
    return -$pam{$x}{$y}; # pam is similarity, want distance
}

# wh, wv scores gaps in first, second sequence 
sub wh {
    return $gxpen if shift( @_) =~ /-/; # reduced penalty to traverse preexisting gap
    return $gpen;
}

sub wv {
    return $gxpen if shift( @_) =~ /-/; # reduced penalty to traverse preexisting gap
    return $gpen;
}


sub nw { #align a pair of seqs by Needleman-Wunsch
    $s = shift @_;
    $t = shift @_;
    
    @s = unpack( "a$unit" x (length($s)/$unit), $s);
    @t = unpack( "a$unit" x (length($t)/$unit), $t);
    unshift @s, "0";
    unshift @t, "0";
    # init path matrix edges
    $dist[0][0] = 0;
    for ($a=0,$i=1;$i<=$#t;$i++) {
	$a += wh($t[$i]);
	$dist[0][$i] = $a;
	$pthm[0][$i] = $ZI;
    }
    for ($a=0,$i=1;$i<=$#s;$i++) {
	$a += wv($s[$i]);
	$dist[$i][0] = $a;
	$pthm[$i][0] = $IZ;
    }
    # calc internal path matrix values
    for ($i=1; $i<=$#s; $i++) {
	for ($j=1; $j<=$#t; $j++) {
	    # for cell ($i, $j), calculate the vector(s) that
	    # lead to the minimum score increment
	    # put inc scores in keys to hash, whose values are the 
	    # variable names of the coresponding vector codes
	    my %x=();
	    $x{$dist[$i-1][$j] + wv($s[$i])} .= "$IZ " unless ($first_seq_fixed && ($i<=$j));
	    $x{$dist[$i-1][$j-1] + wpam($s[$i], $t[$j])} .= "$II ";
	    $x{$dist[$i][$j-1] + wh($t[$j])} .= "$ZI " unless ($first_seq_fixed);
	    # get min. incremental score, and set cell in distance mx
	    $a = [sort {$a<=>$b} keys %x]->[0]; # capture minimum
	    $dist[$i][$j] = $a;
	    # now put corresponding vector code in path matrix
	    $pthm[$i][$j] = 0;
	    foreach ( @{[split(/ /, $x{$a})]} ) {
		$pthm[$i][$j] |= $_;
	    }
	}
    }
    # traverse recursively all isoscoring, minimum paths, following 
    # possibilities indicated by incremental path vectors in
    # path matrix
    # graphical example in slide 9
    @pathset = ();
    fpath($#s, $#t, ());
    # pathset is now array of refs to two-element arrays containing the 
    # the gap patterns of sequences $s and $t, for each of the minimum
    # scoring alignments

    1;
    return ($dist[$#s][$#t], [@pathset]); # return min score and \@pathset
}

sub fpath { # traverse paths in alignment path matrix
# sets @pathset elements to refs to arrays containing the 
# gap template for the first and second sequences.
    my ($i, $j, @path) = @_;
    
    if ($i+$j == 0) { # finished
	my($x,$y);
	$x = $y = "";
	for (my $l=0; $l < scalar @path; $l+=2) { # transpose path
	    $x .= $path[$l];
	    $y .= $path[$l+1];
	}
	push @pathset, [($x, $y)];
	return;
    } else {
	my @leaf = ();
	push @leaf, [(1,1)] if $pthm[$i][$j] & $II;
	push @leaf, [(1,0)] if $pthm[$i][$j] & $IZ;
	push @leaf, [(0,1)] if $pthm[$i][$j] & $ZI;

	foreach (@leaf) {
	    my @v = @$_;
	    fpath( $i-$v[0], $j-$v[1], @v, @path);
	}
    }
}

sub avg {
    my @a = @_;
    return 0 unless scalar @a;
    my $a = 0;
    foreach (@a) {
	$a += $_;
    }
    return $a / (scalar @a);
}

sub island {
    my $i = shift @_;
    my $j = shift @_;
    return ($pthm[$i-1][$j] & $II) && ($pthm[$i-2][$j-1] & $IZ);
}

sub loadpam {
    $fn = shift @_; # PAM filename
    open (PAM, $fn) or die "Can't find PAM matrix $fn\n";
    %pam = ();
    for ($line = <PAM>; $line =~ /\#/; $line = <PAM>) {
	 1;
     }
    while ($line !~ /Q/) { # find amino acid name row
	$line = <PAM>;
	chop;
    }
    $line =~ s/^\s+//;
    @aa = split(/\s+/, $line);
    foreach $aa (@aa) {
	$line = <PAM>; chop;
	@scr = split(/\s+/, $line);
	die "AA rows out of order\n" if ($scr[0] != $_);
	shift @scr;
	foreach $aa2 (@aa) {
	    $pam{$aa}->{$aa2} = shift @scr;
	}
    }
}

sub mkseq { # return a gap-inserted seq, using string of bases and binary gap template, last arg opt. num bases
# call: mkseq( $seq, $template )
    my $seq = shift @_;
    my @t = split ("", shift(@_));
    my $unit = shift @_;
    my $i=0, $rtn="";
    $unit = 1 if !$unit;
    foreach (@t) {
	if ($_) {
	    $rtn .= substr($seq, $i, $unit); # next element if '1'
	    $i+=$unit;
	} else {
	    $rtn .= "-" x $unit; # gap if '0'
	}
    }
    return $rtn;
}

sub pralign {
# print a single alignment
# call with a pointer to alignment array and seq list
    my @parms = @_;
    my @align = @{$parms[0]};
    my @seqlist = @{$parms[1]};
    for (my $i = 0; $i < $#align; $i++) {
	print [split(/ /, $seqlist{$seqlist[$i]})]->[0], "\t";
	print mkseq($seqlist[$i],$align[$i]),"\n";
    }
print "score: $align[-1]\n";
}
