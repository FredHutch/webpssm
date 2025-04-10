#!/usr/bin/perl

use strict;

use Bio::Seq;
use Bio::Tools::IUPAC;
use lib "$ENV{'DOCUMENT_ROOT'}/lib";
use Path;
use PSSM2;
use FileHandle;
use Data::Dumper;
use Email::Sender::Simple qw(sendmail);
use Email::Sender::Transport::SMTP qw();
use Email::Simple;
use Email::Simple::Creator; # For creating the email

# overwrite the existing file or not. Default is to overwrite
# chanage the value to 0 if you do not want to overwrite an existing file.
my $g_overwrite=1;

# if you want to restrict upload to files with certain extentions, change
# the value of $g_restrict_by_ext=1 and ALSO modify the @g_allowed_ext if you
# want to add other allowable extensions.
my $g_restrict_by_ext=1;
# case insensitive, so file with Jpeg JPEG GIF gif etc will be allowed
my @g_allowed_ext=("fas,fsa");


my $g_debug = 0;
my $g_title = "File upload";
my $g_upload_path = 'foo';

BEGIN { push(@INC,'.'); }

my $pamfnstem   = shift;
my $x4cons      = shift;
my $r5cons      = shift;
my $matrix      = shift;
my $jobid       = shift;
my $email       = shift;
my $uploadFile  = shift;
my $uploadDir   = shift;
my $addr        = shift;
my $datatype    = shift;
my(%params, $pssmObj, $parsefile, $nt_aa_names, $aaNameSeq, $ntNames);
_set_params($pamfnstem,$x4cons,$r5cons,$matrix,\%params);

my $htmlfile = $uploadDir."/$jobid.pssm.html";
my $csvfile = $uploadDir."/$jobid.pssm.txt";
my $statsfile = $Path::statsbase."/webpssm.log";
my $logfile = $uploadDir."/$jobid.log";
open LOG, ">", $logfile or die "couldn't open $logfile\n";

print LOG "pamfnstem: $pamfnstem\n";
print LOG "x4cons: $x4cons\n";
print LOG "r5cons: $r5cons\n";
print LOG "matrix: $matrix\n";
print LOG "jobid: $jobid\n";
print LOG "email: $email\n";
print LOG "uploadFile: $uploadFile\n";
print LOG "uploadDir: $uploadDir\n";
print LOG "addr: $addr\n";
print LOG "datatype: $datatype\n";
print LOG "logFile: $logfile\n";
print LOG "params: ".Dumper(\%params)."\n";
close LOG;

if ($datatype eq 'nt') {	# nucleotide sequence - need to translate aa sequence first
	my $aaFile = $uploadDir."/$jobid.aa.fasta";
	($ntNames, $nt_aa_names, $aaNameSeq) = translate($uploadFile, $aaFile);
	$uploadFile = $aaFile;
	if (keys %{$aaNameSeq} > 10000) {
		my $errFile = $uploadDir."/$jobid.err.txt";
		open ERR, ">$errFile" or die "couldn't open $errFile: $!\n";
		print ERR "<p>Error: The translated amino acid sequences exceed the the maximum number of amino acid sequences of 10000. Please check your input nucleotide sequences and divide them into smaller inputs.</p>\n";
		exit;
	}
}
$pssmObj = PSSM->new(\%params);
my(@pssm_results) = do_pssm($jobid,\%params,$uploadFile,$pssmObj);
print_html($htmlfile, $jobid,\%params,$pssmObj,$datatype,$ntNames,$nt_aa_names,$aaNameSeq,@pssm_results);
print_csv($csvfile, $jobid,\%params,$datatype,$ntNames,$nt_aa_names,$aaNameSeq,@pssm_results);
print_stats($statsfile, $jobid, $uploadFile, $email, $addr);

#create a file to indicate the status of the finished job.
open TOGGLE, ">", "$uploadDir/toggle" or die "couldn't create file toggle\n";
close TOGGLE;
	
if ($email) {
	my $emailbody = "<p>Your job #$jobid has finished on our server. Please click 
	<a href=https://webpssm.fredhutch.org/cgi-bin/processpssm.cgi?jobid=$jobid>
	here</a> to get result. The result will be kept for 5 days after this message was sent.</p>
	<p>If you have any questions please email to mullspt\@uw.edu. Thanks.</p>";	

	# Create the email
	my $cemail = Email::Simple->create(
		header => [
			#To => '"Recipient Name" <recipient@fredhutch.org>',
			#From => '"Sender Name" <sender@fredhutch.org>',
			To => $email,
			From => 'webpssm@fredhutch.org',
			Subject => "Your WebPSSM #$jobid Result",
		],
		body => $emailbody,
	);
	$cemail->header_set( 'Content-Type' => 'Text/html' );
	$cemail->header_set( 'Reply-To' => 'mullspt@uw.edu' );
	
	# Configure the SMTP transport
	my $transport = Email::Sender::Transport::SMTP->new({
		host => 'mx.fhcrc.org', # Your SMTP server address
		port => 25, # Common ports are 25, 465, or 587
		ssl => 0, # Set to 1 if SSL is required
		# sasl_username => 'your_username', # Your SMTP username
		#sasl_password => 'your_password', # Your SMTP password
	});
	
	# Send the email
	eval {
		sendmail($cemail, { transport => $transport });
		print "Email sent successfully!\n";
	};
	if ($@) {
		die "Failed to send email: $@\n";
	}
}


sub parse_file {
	my($file) = @_;
	open(FH, $file);
	undef $/;
	my $lines = <FH>;
	#print STDERR "\n parse_file: " . $lines . "\n";
	close FH;
	return $lines;
}


sub translate {
	my ($uploadFile, $aaFile) = @_;
	my (@ntNames, %ntNameSeq, %nt_aa_names, %aaNameSeq, @aaNames);
	my $name = '';
	open IN, $uploadFile or die "couldn't open $uploadFile: $!\n";
	while (my $line = <IN>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>(.*)$/) {
			$name = $1;
			push @ntNames, $name;
		}else {
			$line =~ s/\-//g;
			$ntNameSeq{$name} .= uc $line;
		}
	}
	close IN;
	my $flag = 0;
	foreach my $ntName (@ntNames) {
		my $ntSeq = $ntNameSeq{$ntName};
		my $aaIdx = 0;
		my %aaSeqStatus = ();
		my $ambiseq = new Bio::Seq (-seq => $ntSeq, -alphabet => 'dna');
		my $stream = new Bio::Tools::IUPAC (-seq => $ambiseq);
		while (my $uniqueseq = $stream->next_seq()) {
			my $aaObj = $uniqueseq->translate();
			my $aaSeq = $aaObj->seq();
			if (!$aaSeqStatus{$aaSeq}) {
				$aaSeqStatus{$aaSeq} = 1;
				$aaIdx++;				
				my $aaName = $ntName.'_'.$aaIdx;
				push @aaNames, $aaName;
				push @{$nt_aa_names{$ntName}}, $aaName; 
				$aaNameSeq{$aaName} = $aaSeq;
				if ($aaIdx > 10000) {
					$flag = 1;
					last;
				}
			}			
		}
		last if $flag;
	}
	unless ($flag) {
		open OUT, ">$aaFile" or die "couldn't open $aaFile: $!\n";
		foreach my $aaName (@aaNames) {
			my $aaSeq = $aaNameSeq{$aaName};
			unless ($aaSeq =~ /\*/) {
				print OUT ">$aaName\n";
				print OUT "$aaSeq\n";
			}
		}
		close OUT;
	}	
	return (\@ntNames, \%nt_aa_names, \%aaNameSeq);
}


sub do_pssm {
	my($jobid, $params,$uploadFile,$pssmObj) = @_;
	my(@a,$i,$inc,@newSeqids);
	
	# print STDERR "\n do_pssm: \n" . $raw_seqs . "\n";
	
	my($seqidAref, $seqsAref) = fasta_parser($jobid, $pssmObj,$uploadFile);
	
	# foreach my $j (0 .. @{$seqidAref}-1) {
	# 	print STDERR " do_pssm: '$seqidAref->[$j]' is '$seqsAref->[$j]'\n";
	# }
	
	## $pssmObj->setmatrix("x4r5"); # could also be "sinsi"
	$pssmObj->setmatrix($params->{matrix});
	my($prseqs,$alscr,$alpct) = $pssmObj->prepseq($seqsAref);
	
	#foreach(@{$prseqs}) {
	#	print STDERR "do_pssm prseqs: $_\n";
	#}
	
	for (@a=(),$i=0; $i<=@{$prseqs}-1; $i++) {
		# has alternative alignments
    	push (@a, $i) if $prseqs->[$i] =~ / /;
	}

	foreach (@a) {
	    $_ += $inc; # put index in register based on previous splice
		# get the alternative aligns
    	my(@alt) = split(/ /, $prseqs->[$_]);

		## Add a notation to the pssmObj about our alternative aligns 
		$pssmObj->{alts}->{$seqidAref->[$_]} = [ @alt ];

    	# replicate info for this sequence appropriately
    	splice(@{$prseqs}, $_, 1, @alt);
    	splice(@{$seqidAref}, $_, 1, ($seqidAref->[$_]) x @alt);
    	splice(@{$alscr}, $_, 1, ($alscr->[$_]) x @alt);
    	splice(@{$alpct}, $_, 1, ($alpct->[$_]) x @alt);
    	$inc += @alt-1;
	} 
	
	my($scr,$pred,$x4pct,$r5pct,$gen,$pos,$net) = $pssmObj->dopssm($prseqs);
	return($prseqs,$seqidAref,$seqsAref,$scr,$pred,$x4pct,$r5pct,$gen,$pos,
		$net,$alpct);
}


sub fasta_parser {
	my($jobid, $pssmObj,$uploadFile) = @_;	
	my(@seqids,@seqs,$lastline,$lcount,$arcount,$switch);
	my(@fasta);
	my($fasta_aref) = $pssmObj->parse_file($uploadFile);
	foreach(@{$fasta_aref}) {
		push(@fasta,$_);
	}
	@fasta = sort {$a <=> $b} @fasta;

	foreach(@fasta) {
		# 0print STDERR "FASTA: $_\n";
		++$lcount;
		s/\s+$//;
		if (!$_) {
			# blank line, just skip it!
			print STDERR "fasta_parser: Blank line.\n";
			next;
		} elsif ($_ =~ /^\>/ && $lastline !~ /^\>/) {
			# print STDERR "fasta_parser: got line with seqid: $_\n";
			if ($switch) { 
				++$arcount; 
			} else { 
				$switch = 1; 
				$arcount = 0; 
			}
			# a new sequence id
			$lastline = $_;
			s/^\>//;
			my($id,@misc) = split(/\s+/);
			$id =~ s/[\s,]+//g;
			## Do we need a uniq id for each seq?
			#$id = _gen_uniqid(\@seqids) if(!$id);
			my $idx = $arcount + 1;
			$id = 'seq_'.$idx if !$id;
			
			# print STDERR  "id: " . $id . "\n";
			
			#$lastseqid = $id;
			$seqids[$arcount] = $id;
		} elsif ($_ !~ /^\>/ && $_ =~ /^[A-Za-z\-]+$/ && $lastline =~ /^\>/) {
			# first part of the actual sequence
			# print STDERR "fasta_parser: got seq: $_\n";
			s/\s+//g;
			$seqs[$arcount] = uc($_);
			$lastline = $_;
		} elsif ($_ !~ /^\>/ && $lastline !~ /^\>/ && $_ =~ /^[A-Za-z\-]+$/) {
			# more lines of a sequence
			# print STDERR "fasta_parser: got multiline seq: $_\n";
			s/\s+//g;
			$seqs[$arcount] .= uc($_);
			$lastline = $_;
		} else {
			print STDERR "Input not in FASTA format";
			return("ERROR: Input not in FASTA format (line $lcount)"
					. ": \"$_\"");
		}
	}
	return(\@seqids,\@seqs);
}

sub _gen_uniqid {
	my($seqs) = @_;
	my($id,$i,$j);
	until($i) {
		$id = int(rand(1000));
		$id = "seq_" . $id;
		foreach(@{$seqs}) {
			if($_ eq $id) { ++$j; }
		}
		++$i if(!$j);
	}
	return($id);
}

sub _set_params {
	my($pamfnstem,$x4cons,$r5cons,$matrix,$params) = @_;
	$params->{pamfnstem} = $pamfnstem;
	$params->{x4cons} = $x4cons;
	$params->{r5cons} = $r5cons;
	$params->{matrix} = $matrix;
	return();
}


sub print_csv {
	my($csvfile, $jobid,$params,$datatype,$ntNames,$nt_aa_names,$aaNameSeq,$prseqs,$seqidAref,$seqsAref,$scr,$pred,$x4pct,$r5pct,$gen,$pos,$net,$alpct)
		= @_;
#	my($url) = "data/results/$jobid.pssm.txt";	
#	my($file) = $uploadDir."/$jobid.pssm.txt";
	
	open(CSV, ">$csvfile") || print STDERR "Cannot open $csvfile\: $!\n";
	print CSV "WebPSSM ($params->{matrix}) Results for id $jobid\n";
	if ($datatype eq 'aa') {
		print CSV "name\tscore\tpred\tx4.pct\tr5.pct\tgeno\tpos.chg\tnet.chg\tpercentile\n";	
		my($i);
		for ($i=0; $i<=@{$prseqs}-1; $i++) {
			$scr->[$i] = sprintf("%.2f", $scr->[$i]);
			$x4pct->[$i] = sprintf("%.2f", $x4pct->[$i]);
			$r5pct->[$i] = sprintf("%.2f", $r5pct->[$i]);
			$alpct->[$i] = sprintf("%.2f", $alpct->[$i]);
			print CSV "$seqidAref->[$i]\t";	# sequence name
			print CSV "$scr->[$i]\t";		# score	
			print CSV "$pred->[$i]\t";		# predicted
			print CSV "$x4pct->[$i]\t";		# x4 percentile
			print CSV "$r5pct->[$i]\t";		# r5 percentile
			print CSV "$gen->[$i]\t";		# genotype
			print CSV "$pos->[$i]\t";		# pos.chg
			print CSV "$net->[$i]\t";		# net.chg
			print CSV "$alpct->[$i]\n";		# percentile
		}
	}else {	#nucleotide sequences
		print CSV "name\ttranslated amino acid sequence\tscore\tpred\tx4.pct\tr5.pct\tgeno\tpos.chg\tnet.chg\tpercentile\n";
		my $i = 0;
		foreach my $ntName (@$ntNames) {
			foreach my $aaName (@{$nt_aa_names->{$ntName}}) {
				my $aaSeq = $aaNameSeq->{$aaName};
				if ($aaSeq =~ /\*/) {
					print CSV "$ntName\t$aaSeq\n";
				}else {
					while ($aaName eq $seqidAref->[$i]) {
						$scr->[$i] = sprintf("%.2f", $scr->[$i]);
						$x4pct->[$i] = sprintf("%.2f", $x4pct->[$i]);
						$r5pct->[$i] = sprintf("%.2f", $r5pct->[$i]);
						$alpct->[$i] = sprintf("%.2f", $alpct->[$i]);
						print CSV "$ntName\t";			# sequence name
						print CSV "$aaSeq\t";			# translated aa seqeuence
						print CSV "$scr->[$i]\t";		# score	
						print CSV "$pred->[$i]\t";		# predicted
						print CSV "$x4pct->[$i]\t";		# x4 percentile
						print CSV "$r5pct->[$i]\t";		# r5 percentile
						print CSV "$gen->[$i]\t";		# genotype
						print CSV "$pos->[$i]\t";		# pos.chg
						print CSV "$net->[$i]\t";		# net.chg
						print CSV "$alpct->[$i]\n";		# percentile
						++$i;
					}
				}
			}
		}
		unless ($i == scalar @$prseqs) {
			die "something not right\n";
		}
	}	
	close(CSV);
#	return($url);
}

sub print_stats {
	my ($statsfile, $jobid, $seqFile, $email, $addr) = @_;
	my $seqCount = 0;
	my $date = `date`;
	chomp $date;
	open IN, $seqFile or die "couldn't open $seqFile: $!\n";
	while (<IN>) {
		$seqCount++ if />/;
	}
	close IN;
	open STATS, ">>$statsfile" or die "couldn't open $statsfile: $!\n";
	print STATS "$jobid\n";
	print STATS $jobid,"\t",$seqCount,"\t",$date,"\t",$addr,"\t",$email,"\n";
	close STATS;
}

sub print_html {
	my($htmlfile, $jobid,$params,$pssmObj,$datatype,$ntNames,$nt_aa_names,$aaNameSeq,$prseqs,$seqidAref,$seqsAref,$scr,$pred,$x4pct,$r5pct,$gen,
		$pos,$net,$alpct) = @_;
	my($time) = scalar localtime(time);
	my($gmtime) = scalar gmtime();
#	my($url) = "data/results/$jobid.pssm.html";	
#	my $file = $uploadDir."/$jobid.pssm.html";
#	my $csvFile = $uploadDir."/$jobid.pssm.txt";
	
	open(HTML, ">$htmlfile") || print STDERR "Cannot open $htmlfile\: $!\n";
	print HTML "<a name=TOP>\n";
	print HTML "<div align='center'>\n";
	print HTML "<div class='title'>WebPSSM results for id $jobid ($gmtime GMT)</div>\n";
	print HTML "<div class='data'>Using the $params->{'matrix'} matrix.</div><br>\n";
	print HTML "<div class='nav'>Download these results as a tab-delimited file\n"; 
	print HTML "<a href=download.cgi?id=$jobid>here</a>.<div><br>\n";
	print HTML "<div class='nav'><a href=/index.html>Return</a></div><br>\n";
	print HTML "<table cellspacing=0 cellpadding=7 class=table>\n";
	print HTML "<tr bgcolor=#666666>\n";
	print HTML "<th class='th'>name</td>\n";
	if ($datatype eq 'nt') {
		print HTML "<th class='th'>translated amino acid sequence</td>\n";	
	}
	print HTML "<th><a class='th' href=\#SCORE>score</a></td>\n";
	print HTML "<th><a class='th' href=\#PRED>pred</a></td>\n";
	print HTML "<th><a class='th' href=\#X4>x4.pct</a></td>\n";
	print HTML "<th><a class='th' href=\#R5>r5.pct</a></td>\n";
	print HTML "<th><a class='th' href=\#GENO>geno</a></td>\n";
	print HTML "<th><a class='th' href=\#POS>pos.chg</a></td>\n";
	print HTML "<th><a class='th' href=\#NET>net.chg</a></td>\n";
	print HTML "<th><a class='th' href=\#PCT>percentile</a></td>\n";
	print HTML "</tr>\n";
	my($i,$switch);
	if ($datatype eq 'aa') {
		for ($i=0; $i<=@{$prseqs}-1; $i++) {
			$scr->[$i] = sprintf("%.2f", $scr->[$i]);
			$x4pct->[$i] = sprintf("%.2f", $x4pct->[$i]);
			$r5pct->[$i] = sprintf("%.2f", $r5pct->[$i]);
			if(!$switch) {
				print HTML "<tr>\n";
				$switch = 1;
			} else {
				print HTML "<tr bgcolor=#CCCCCC>\n";
				$switch = 0;
			}
			if($pssmObj->{alts}->{$seqidAref->[$i]})  {
				print HTML "<td><a class='data' href=\"\#$seqidAref->[$i]\">";
				print HTML "$seqidAref->[$i]</a>\#</font></td>\n";
			} else {
				print HTML "<td class='data'>$seqidAref->[$i]</td>\n";
			}
			print HTML "<td align=right class='data'>$scr->[$i]</td>\n";
			print HTML "<td align=right class='data'>$pred->[$i]</td>\n";
			print HTML "<td align=right class='data'>$x4pct->[$i]</td>\n";
			print HTML "<td align=right class='data'>$r5pct->[$i]</td>\n";
			print HTML "<td align=right class='data'>$gen->[$i]</td>\n";
			print HTML "<td align=right class='data'>$pos->[$i]</td>\n";
			print HTML "<td align=right class='data'>$net->[$i]</td>\n";
			my($pct) = sprintf("%.2f", $alpct->[$i]);
			print HTML "<td align=right class='data'>$pct</td>\n";
			print HTML "</tr>\n";
		}
	}else {
		foreach my $ntName (@$ntNames) {
			foreach my $aaName (@{$nt_aa_names->{$ntName}}) {
				
				my $aaSeq = $aaNameSeq->{$aaName};
				if ($aaSeq =~ /\*/) {
					if(!$switch) {
						print HTML "<tr>\n";
						$switch = 1;
					} else {
						print HTML "<tr bgcolor=#CCCCCC>\n";
						$switch = 0;
					}
					print HTML "<td class='data'>$ntName</td>\n";
					print HTML "<td class='data'>$aaSeq</td>\n";
					for (my $i = 0; $i < 8; $i++) {
						print HTML "<td class='data'></td>\n";
					}
					print HTML "</tr>\n";
				}else {
					while ($aaName eq $seqidAref->[$i]) {
						$scr->[$i] = sprintf("%.2f", $scr->[$i]);
						$x4pct->[$i] = sprintf("%.2f", $x4pct->[$i]);
						$r5pct->[$i] = sprintf("%.2f", $r5pct->[$i]);
						if(!$switch) {
							print HTML "<tr>\n";
							$switch = 1;
						} else {
							print HTML "<tr bgcolor=#CCCCCC>\n";
							$switch = 0;
						}
						print HTML "<td class='data'>$ntName</td>\n";
						if($pssmObj->{alts}->{$seqidAref->[$i]})  {
							print HTML "<td><a class='data' href=\"\#$aaSeq\">";
							print HTML "$aaSeq</a>\#</font></td>\n";
						} else {
							print HTML "<td class='data'>$aaSeq</td>\n";
						}
						print HTML "<td align=right class='data'>$scr->[$i]</td>\n";
						print HTML "<td align=right class='data'>$pred->[$i]</td>\n";
						print HTML "<td align=right class='data'>$x4pct->[$i]</td>\n";
						print HTML "<td align=right class='data'>$r5pct->[$i]</td>\n";
						print HTML "<td align=right class='data'>$gen->[$i]</td>\n";
						print HTML "<td align=right class='data'>$pos->[$i]</td>\n";
						print HTML "<td align=right class='data'>$net->[$i]</td>\n";
						my($pct) = sprintf("%.2f", $alpct->[$i]);
						print HTML "<td align=right class='data'>$pct</td>\n";
						print HTML "</tr>\n";
						++$i;
					}
				}
			}
		}
		unless ($i == scalar @$prseqs) {
			die "something not right\n";
		}
	}
	print HTML "</table><P>\n";

	if(keys(%{$pssmObj->{alts}})) {
		print HTML "<table border=0><tr><td align=left class=data>\n";
		print HTML "<strong>Multiple alignments (#) were generated for the";
		print HTML " following sequences:</strong>";
		print HTML "<ul>\n";
		foreach my $nseqid (keys %{$pssmObj->{alts}}) {
			if ($datatype eq 'nt') {
				print HTML "<a name=$aaNameSeq->{$nseqid}>";
				print HTML "<li> $aaNameSeq->{$nseqid}</a>\n";
			}else {
				print HTML "<a name=$nseqid>";
				print HTML "<li> $nseqid</a>\n";
			}
			print HTML "<ul>\n";
			foreach my $align (@{$pssmObj->{alts}->{$nseqid}}) {
				print HTML "<li> $align\n";
			}
			print HTML "</ul>\n";
		}
		print HTML "</ul>\n";
		print HTML "</td></tr></table><BR>\n";
	}
	close HTML;
}
