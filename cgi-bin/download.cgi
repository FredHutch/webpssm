#!/usr/bin/perl -w

use CGI ':standard';
use CGI::Carp qw(fatalsToBrowser); 
use lib "$ENV{'DOCUMENT_ROOT'}/lib";
use Path;

my $cgi = new CGI;
my $id = $cgi->param('id');
my $uploadbase = $Path::uploadbase;
my $fileName = $id.'.pssm.txt';

if ($id !~ /^\d+$/) {
	print "Content-type: text/html\n\n";
	print "Invalid input, process terminated.<br>";
	exit;
}
my $files_location = "$uploadbase/$id";

my @fileholder;

if ($id eq '') { 
	print "Content-type: text/html\n\n"; 
	print "You must specify a file to download."; 
} else {
	open(DLFILE, "<", "$files_location/$fileName") || Error('open', 'file'); 
	@fileholder = <DLFILE>; 
	close (DLFILE) || Error ('close', 'file'); 	
	print "Content-Type:application/x-download\n"; 
	print "Content-Disposition:attachment;filename=$fileName\n\n";
	print @fileholder;
}

sub Error {
	print "Content-type: text/html\n\n";
	print "The server can't $_[0] the $_[1]: $! \n";
	exit;
}
