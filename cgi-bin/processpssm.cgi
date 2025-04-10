#!/usr/bin/perl

use strict;
use DBI;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use lib "$ENV{'DOCUMENT_ROOT'}/lib";
use Path;
use File::Basename;

my $q = new CGI;
print $q->header;

my $jobid = $q->param('jobid') || time().int(rand(90)+10);
my $email = $q->param('email');
my $seqs = $q->param('seqs');
my $seqFile = $q->param('seqFile');
my $datatype = $q->param('datatype');
my $matrix = $q->param("matrix") || "x4r5";
my $dopssm = $q->param('dopssm');
my $pamfnstem = $q->param('pamfnstem') || "HIVPAM10";
my $x4cons = $q->param('x4cons') || "CTRPNNNTRKSIHIGPGRAFYTTGRIIGDIRQAHC";			
my $r5cons = $q->param('r5cons') || "CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC";
my $remote_addr = '';

if ($dopssm) {
	$remote_addr = $ENV{'REMOTE_ADDR'};		
	if ($email =~ /^\s*$/) {
		$email = '';
	}else {
		$email =~ s/^\s+//;
		$email =~ s/\s+$//;
	}
}

print <<END_HTML;

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
        "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
	<title>WebPSSM Result</title>
	<link href="/stylesheets/pssm.css" media="screen" rel="Stylesheet" type="text/css" />
	<link rel="stylesheet" href="/stylesheets/spin.css">
</head>
<body>
	<div id="wrap">

    <div id="header">
	    <div class="spacer">&nbsp;</div>    
		<span class="logo">WebPSSM Result</span>   
    </div>
    
    <div id="nav">
	</div>
	
	<div class="spacer">&nbsp;</div>
	
END_HTML

my $uploadBase = $Path::uploadbase;
my $uploadDir = $uploadBase."/$jobid";
my $csvfile = $uploadDir."/$jobid.pssm.txt";
my $htmlfile = $uploadDir."/$jobid.pssm.html";
my $errfile = $uploadDir."/$jobid.err.txt";
my $toggle = $uploadDir.'/toggle';

my $errMsg = '';
if ($dopssm) {	
	# upload files
	mkdir $uploadDir;
	chmod 0755, $uploadDir;
	my $uploadFile = $uploadDir."/$jobid.seq.fasta";
	my ($cleanLinsRef, $seqFlag);	
	if ($seqs) {	# paste sequences
		$seqs =~ s/\r\n/\n/g;
		$seqs =~ s/\r/\n/g;
		my @lines = split /\n/, $seqs;
		($cleanLinsRef, $seqFlag, $errMsg) = CleanSeqs (\@lines, $datatype);
		if ($seqFlag eq 'error') {
			open ERR, ">$errfile" or die "couldn't open $errfile: $!\n";
			print ERR "$errMsg\n";
			close ERR;
		}else {
			WriteFile ($uploadFile, $cleanLinsRef);
		}		
	}elsif($seqFile) {	# upload file
		my $fh = $q->upload("seqFile");
		my $FileLines = GetFileLines ($fh);
		($cleanLinsRef, $seqFlag, $errMsg) = CleanSeqs ($FileLines, $datatype);
		if ($seqFlag eq 'error') {
			open ERR, ">$errfile" or die "couldn't open $errfile: $!\n";
			print ERR "$errMsg\n";
			close ERR;
		}else {
			WriteFile ($uploadFile, $cleanLinsRef);
		}		
	} 
	unless ($seqFlag eq 'error') {
		my @paramsArray;
		push @paramsArray, $pamfnstem, $x4cons, $r5cons, $matrix, $jobid, $email, $uploadFile, $uploadDir, $remote_addr, $datatype;
		
		# fork.
		my $pid = fork();
		die "Failed to fork: $!" unless defined $pid;

		if ($pid == 0) {
			# Execute the background process.
			close(STDIN);
			close(STDOUT);
			close(STDERR);
			open(STDIN,  "</dev/null");
			open(STDOUT, ">/dev/null");
			open(STDERR, ">/dev/null");
			exec("perl", "webpssm4.pl", @paramsArray);
			exit(0);
		}
	}	
}

if (-s $errfile) {
	open ERR, $errfile or die "couldn't open $errfile: $!\n";
	while (my $line = <ERR>) {
		chomp $line;
		print $line;
	}
	close ERR;
}else {
	if(!-e $toggle) {
		print "<div>";
		print "<h3>Your job #$jobid is being processed  </h3>";
			print "<div class=\"spinner\">";
				#print "<div class=\"circle\"></div>";
			print "</div>";
			print "</div>";
		if($email) {
			print "<p>Result will be sent to <b>$email</b> when the job finishes.";
			print "<p>You can close browser if you want.";
		}else {
			print "<p>Please wait here to watch the progress of your job.</p>";
			print "<p>This page will update itself automatically until job is done.</p>";
		}
		print "<script>";
			print "function autoRefresh() {";
				print "location.href = \"processpssm.cgi?jobid=$jobid&email=$email\";";
			print "}";
			print "setInterval('autoRefresh()', 10000);";
		print "</script>";
	}

	if (-e $toggle and -s $csvfile) {
		print_result ($htmlfile);
	}
}

PrintFooter();

#####################################################################
sub GetFileLines {
	my $fh = shift;
	my $line = "";
	my @buffer = <$fh>;
 	foreach my $element (@buffer) {
 		$line .= $element;
 	}
 	if ($line =~ /\r\n/) {
		$line =~ s/\r//g;
	}elsif ($line =~ /\r/) {
		$line =~ s/\r/\n/g;
	}
	my @fileLines = split /\n/, $line;
	return \@fileLines;
}


######################################################################################
sub WriteFile {
	my ($uploadFile, $fileLines) = @_;
	open OUT, ">$uploadFile" or die "couldn't open $uploadFile: $!\n";
	foreach my $line (@$fileLines) {
		print OUT $line,"\n";		
	}
}


######################################################################################
sub CleanSeqs {
	my $lines = shift;
	my $datatype = shift;
	my (%nameStatus, $seq, @cleanLines);
	my $flag = my $seqCount = 0;
	my $seqFlag = 'success';
	my $errMsg = '';
	foreach my $line (@$lines) {
		next if $line =~ /^\s*$/;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		if (!$flag) {
			if ($line =~ /^>(.*)$/) {
				my $name = $1;
				$nameStatus{$name} = 1;
				$name =~ s/\W/_/g;
				$line = '>'.$name;
				push @cleanLines, $line;
				$flag = 1;
				$seqCount++;
			}else {
				print "<p>Error: The input sequences are not in fasta format.</p>";
				PrintFooter();
			}
		}else {
			if ($line =~ /^>(.*)$/) {
				my $name = $1;
				$seqCount++;
				if ($seqCount > 10000) {
					$errMsg = "<p>Error: The maximum number of input sequences is 10000. Please divide your input into smaller inputs and submit them again.</p>";
					$seqFlag = 'error';
					last;
				}
				if (!$nameStatus{$name}) {
					$nameStatus{$name} = 1;
				}else {
					$errMsg = "<p>Error: The input sequences have the same name of $name. Sequence name must be unique.</p>";
					$seqFlag = 'error';
					last;
				}				
				my $seqNogaps = $seq;
				$seqNogaps =~ s/\-//g;
				if (!$seqNogaps) {	# in case of no sequence after name
					pop @cleanLines;	# remove the name
				}else {
					if ($datatype eq 'aa' && length $seqNogaps > 50) {
						$errMsg = "<p>Error: The input sequences are probably not V3 loop amino acid sequences. The length of some sequences is longer than 50. 
						The length of V3 loop is around 35 amino acids. Please check your sequences. You data will not been processed.</p>";
						$seqFlag = 'error';
						last;
					}elsif ($datatype eq 'aa' && $seqNogaps =~ /[^A-Za-z]/) {
						$errMsg = "<p>Error: The input sequences contain unrecognized amino acid characters (e.g., #, *, ?). Please check your sequences. You data will not been processed.</p>";
						$seqFlag = 'error';
						last;
					}elsif ($datatype eq 'nt' && (length $seqNogaps > 150 || length $seqNogaps < 50)) {
						$errMsg = "<p>Error: The input sequences are probably not V3 loop nuclotide sequences. The length of some sequences is longer than 150 or shorter than 50. 
						The length of V3 loop is around 105 nt. Please check your sequences. You data will not been processed.</p>";
						$seqFlag = 'error';
						last;
					}elsif ($datatype eq 'nt' && $seqNogaps =~ /[ILFPQEX]/i) {
						$errMsg = "<p>Error: The input sequences are probably not nuclotide sequences. 
						Sequences contain unrecognized characters for nucleotide sequences that may be caused by choosing wrong data type of your input sequences.</p>";
						#PrintFooter();
						$seqFlag = 'error';
						last;
					}else {
						push @cleanLines, $seq;
					}			
				}
				$name =~ s/\W/_/g;
				$line = '>'.$name;
				push @cleanLines, $line;
				$seq = '';
			}else {
				$seq .= $line;
			}
		}		
	}

	if ($seqFlag eq 'success') {
		# last sequence
		my $seqNogaps = $seq;
		$seqNogaps =~ s/\-//g;
		if (!$seqNogaps) {	# in case of no sequence after name
			pop @cleanLines;	# remove the name
		}else {
			if ($datatype eq 'aa' && length $seqNogaps > 50) {
				$errMsg = "<p>Error: The input sequences are probably not V3 loop amino acid sequences. The length of some sequences is longer than 50. 
				The length of V3 loop is around 35 amino acids. Please check your sequences. You data will not been processed.</p>";
				$seqFlag = 'error';
			}elsif ($datatype eq 'aa' && $seqNogaps =~ /[^A-Za-z]/) {
				$errMsg = "<p>Error: The input sequences contain unrecognized amino acid characters (e.g., #, *, ?). Please check your sequences. You data will not been processed.</p>";
				$seqFlag = 'error';
			}elsif ($datatype eq 'nt' && (length $seqNogaps > 150 || length $seqNogaps < 50)) {
				$errMsg = "<p>Error: The input sequences are probably not V3 loop nuclotide sequences. The length of some sequences is longer than 150 or shorter than 50. 
				The length of V3 loop is around 105 nt. Please check your sequences. You data will not been processed.</p>";
				$seqFlag = 'error';
			}elsif ($datatype eq 'nt' && $seqNogaps =~ /[ILFPQEX]/i) {
				$errMsg = "<p>Error: The input sequences are probably not nuclotide sequences. 
				Sequences contain unrecognized characters for nucleotide sequences that may be caused by choosing wrong data type of your input sequences.</p>";
				$seqFlag = 'error';
			}else {
				push @cleanLines, $seq;
			}
		}
		$seq = '';
	}

	return (\@cleanLines, $seqFlag, $errMsg);
}


######################################################################################
sub print_result {
	my ($htmlfile) = @_;
	open HTML, $htmlfile or die "couldn't open $htmlfile: $!\n";
	while (my $line = <HTML>) {
		print $line;
	}
	close HTML;
	print "<a href=/index.html>Return</a><P>\n";
	print "<P>";
	print qq!
<table border=0 cellpadding=7 bgcolor=#CCCCCC width=400 class=table>
<tr><td align=left valign=top>
<div class=th align=center style=color:black>Explanation of results</div>
<div class=data>
<a name=SCORE><a href=\#TOP>score</a>:<BR>
Sequence PSSM score
<P>
<a name=PRED><a href=\#TOP>pred</a>:<BR>
1: CXCR4 use predicted; 0: CCR5 use only predicted
<P>
<a name=X4><a href=\#TOP>x4.pct</a>:<BR>
Percentile of score within score distribution of known X4 viruses.
<P>
<a name=R5><a href=\#TOP>r5.pct</a>:<BR>
Percentile of score within score distribution of known R5 viruses.
<P>
<a name=GENO><a href=\#TOP>geno</a>:<BR>
Amino acid residues at V3 sites 11 and 25; basic residues (R or K) at
either or both of these sites is predictive of CXCR4 use.
<P>
<a name=POS><a href=\#TOP>pos.chg</a>:<BR>
Total number of positively charged (R/K/H) amino acid residues. Higher
positive charge has been correlated with likelihood of CXCR4 use.
<P>
<a name=NET><a href=\#TOP>net.chg</a>:<BR>
Number of positively charged (R/K/H) amino acid residues, minus number
of negatively charged (D/E) residues.
<P>
<a name=PCT><a href=\#TOP>percentile</a>:<BR>
If this number is over .95, this sequence is most likely not a V3.
</div>
</td></tr></table>
</div>
!;
	print "<br><a href=/index.html>Return</a></center></div></div><br>\n";
}


######################################################################################
sub PrintFooter {
	print "<div id='footer'>
				<p class='copyright'>&copy; 2025 Fred Hutch Cancer Center. All rights reserved.</p>
			</div>
			</div>
		</body>
		<html>";
	exit 0;
}
