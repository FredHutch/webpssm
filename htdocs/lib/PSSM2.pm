package PSSM;

use NW;
use Path;
use FileHandle;
use strict 'refs';
use strict 'vars';

BEGIN { push(@INC,'.'); }

sub new {
	my($class,$inputhash) = @_;
	#my $path = "/data/htdocs/pssm";
	my $path = $Path::documentroot;
	# which PAM?
	$inputhash->{pamfnstem} = "HIVPAM10" if(!$inputhash->{pamfnstem});
	# consensi based on Resch et al. 2001, Virology 288:51-62 data
	$inputhash->{x4cons} = "CTRPNNNTRKSIHIGPGRAFYTTGRIIGDIRQAHC" 
		if(!$inputhash->{x4cons});
	$inputhash->{r5cons} = "CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC" 
		if(!$inputhash->{r5cons});

	my(@adist);
	open DIST, "$path/$inputhash->{pamfnstem}\.dist" 
		|| die "Alignment score distribution file problem: $!\n";
	while (<DIST>) {
		# print STDERR "PAMFNSTEM: $_\n";
    	push(@adist, $_);
	}
	close(DIST);

	@adist = sort {$a <=> $b} @adist;
	loadpam("$path/$inputhash->{pamfnstem}\.pam");

	my $self = {
			pamfnstem	=> $inputhash->{pamfnstem},
			x4cons		=> $inputhash->{x4cons},
			r5cons		=> $inputhash->{r5cons},
			adist		=> \@adist,
	};

	return bless ($self, $class);
}

sub quit {
	my($self) = @_;
	undef($self);
	return();
}

sub setmatrix { 
	# get matrix and assoc. empirical distributions from specified file
	# sets %pssm, @x4dist, @r5dist, $x4cut, $r5cut
	# call: setmatrix( $matrixfilenamestem )

	my($self,$fnstem) = @_;
	
	my $path = $Path::documentroot."/matrix";
	my($matrix_aref) = parse_file($self,"$path/$fnstem.matrix");
	my(%pssm,$x4cut,$r5cut,$optcut,@x4dist,@r5dist,$i);
	foreach(@{$matrix_aref}) {
		# print STDERR "SETMATRIX: $_\n";
		my(@a);
		## skip first line; assume it contains site numbers
		if(!$i) { ++$i; next; } 
		++$i;

		@a = split /\s+/;
		# print STDERR "SETMATRIX: @a\n";

		if (/xcut/) { # get the x4 cutoff value
			$x4cut = $a[1]; 
		} elsif (/rcut/) { # get the r5 cutoff value
			$r5cut = $a[1];
		} elsif (/optcut/) {
			$optcut = $a[1];
		} else {
			# assume first string on line is one-letter 
			my($a) = shift(@a); 
			# AA designator
			$pssm{$a} = \@a; # rest are site scores for that aa
		}
	}
	undef($i);

	if (-f "$path/$fnstem.x4dist") {
		my($x4dist_aref) = parse_file($self,"$path/$fnstem.x4dist");

		foreach(@{$x4dist_aref}) {
			push(@x4dist,$_);
		}
		@x4dist = sort {$a <=> $b} @x4dist;

		$self->{x4dist} = \@x4dist;
	}

	if (-f "$path/$fnstem.r5dist") {
		my($r5dist_aref) = parse_file($self,"$path/$fnstem.r5dist");
		foreach(@{$r5dist_aref}) {
			push(@r5dist, $_);
		}
		@r5dist = sort {$a <=> $b} @r5dist;

		$self->{r5dist} = \@r5dist;
	}

	$self->{pssm} = \%pssm;
	$self->{r5cut} = $r5cut;
	$self->{x4cut} = $x4cut;
	$self->{optcut} = $optcut;

	return();
}


# prepseq
# align V3 seq to consensus V3s, to produce a PSSM-scorable sequence
# using Needleman-Wunsch and PAM matrix
#

sub prepseq { 
	# call: prepseq( \@seqs ), @seqs are the V3s to be prepped
	# return: \@(scorable seqs), \@(alignment scores), 
	# \@(percentile of score)

	my($self,$seqs) = @_;
	my(@ret, $ret, @scr, @pct);

	# do alignment to R5 or X4 consensus based on whether "G[KR][IV]"
	# motifs present: yes, use X4 cons; no, R5 cons.
	foreach (@{$seqs}) {
		my($cons) = (/G[KR][IV]/) ? $self->{x4cons} : $self->{r5cons};
		# get rid of pre-inserted gaps
		s/-//g;
		my($score, $templates) = nw($cons, $_);
		push @scr, $score;
		push(@pct, pctile($self, $score, $self->{adist}));
		# create scorable sequence(s)
		$ret = "";
		foreach my $tpl (@{$templates}) {
			my($i,$j,$s);
			for ($i=$j=0,$s=""; $i<=length($$tpl[0]); $i++) {
				if (substr($$tpl[0], $i, 1)==1) {
					if (substr($$tpl[1],$i,1)==1) {
						$s .=  substr($_,$j++,1);
					} else { # gap 
						$s .= "-";
					}
				} else { # delete bases in target
					$j++;
				}
			}
			# separate isoscoring scorable seqs by spaces
			$ret .= "$s "; # to be split later
		}
		chop($ret); #lose the last space
		push(@ret,$ret);
	}
	return (\@ret,\@scr,\@pct);
}

sub dopssm { 
	# calculate scores and predictions for an array of V3 seqs
	# call: dopssm( \@seqs )
	# return ([@(score)],[@(predictions)],[@(x4 pctiles)],[@(r5 pctiles)],
	#         [@(aas at 11 & 25)],
	#         [@(no. pos-charges)],[@(net pos-charges)]);
	
	my($self,$seqs) = @_;
	return () if (!$self->{pssm});
	my(@scr, $scr, @net, @pos, @x4pct, @r5pct, @gen, $gen, @pred);

	my($len) = length($seqs->[0]);

	foreach (@{$seqs}) {
		# print STDERR "PSSM dopssm: '$_'\n";
		my(@a) = split ''; ## score
		
		#foreach(@a) {
		#	print STDERR "'$_' ";
		# }
		# print STDERR "\n";
		
		my($i);
		for ($scr=$i=0;$i<=$len;$i++) {
			$scr += $self->{pssm}->{$a[$i]}->[$i];
			# print STDERR "PSSM: $self->{pssm}->{$a[$i]}->[$i]\n";
			# print STDERR "SCR: $scr\n";
		}
		undef($i);
		
		push(@scr,$scr);
	
		# charge info
		@a = /[RKH]/g; # length of @a is no. positive charges
		my(@b) = /[DE]/g; # length of @b is no. neg. charges
		push(@pos, scalar @a);
		push(@net, (scalar @a - scalar @b));
	
		# genotype
		$gen = substr($_,10,1) . substr($_,24,1);
		push(@gen,$gen);
	
		# percentiles
		push(@x4pct, pctile($self, $scr, $self->{x4dist}));
		push(@r5pct, pctile($self, $scr, $self->{r5dist}));
	
		# prediction (1 if predicted X4, 0 if R5)

		#optcut matrices predict w/ single value, don't use canonical sites
		if (defined $self->{optcut}) {
			push @pred, ($scr >= $self->{optcut}) ? 1 : 0;
		} 
		#other matrices use two cutoff values and canonical sites to predict
		else {
			push(@pred, (($scr > $self->{x4cut}) 
					|| (($scr > $self->{r5cut}) 
					&& ($gen =~ /[RK]/))) ? 1 : 0);
		}
	}
	return (\@scr,\@pred,\@x4pct,\@r5pct,\@gen,\@pos,\@net);
}


sub pctile {
	# calc estimated percentile of number in a given empirical distribution
	# (distribution assumed sorted ascending)
	# call: pctile( $self, $x, \@dist )
    
    my($self,$x,$dist) = @_;
    my($i);
    return unless ref $dist;
    for ($i=0; $i<=scalar(@{$dist}); $i++) {
		last if ($x < $dist->[$i]);
    }
    $i = scalar(@{$dist}) if ($i > scalar(@{$dist}));
    
    return ($i/scalar(@{$dist}));
}


sub parse_file {
	my($self,$file) = @_;
	my(@file,$file_line);
	my($file_fh) = new FileHandle "$file";
	if (!$file_fh) {
		die "couldn't open $file";
	}
	while(my($c) = $file_fh->getc) {
		if($c !~ /\n/) {
			$file_line .= $c;
			if($file_fh->eof) { 
				$file_fh->close;
				push(@file,$file_line);
				$file_line = undef;
				last;
			}
		} elsif($c =~ /\n/) {
			push(@file,$file_line);
			$file_line = undef;
			if($file_fh->eof) { 
				$file_fh->close; 
				$file_line = undef;
				last; 
			}
		}
	}
	return(\@file);
}

1;
