#!/usr/bin/perl

use strict;
use File::Path;

my $dir = '/usr/local/apache2/htdocs/outputs';
opendir DIR, $dir;
while (my $subdir = readdir DIR) {
	next if $subdir =~ /^\./;
	$subdir = $dir.'/'.$subdir;
	my $rmflag = 1;
	if (-d $subdir) {
		my $togglefile = $subdir."/toggle";
		if(-e $togglefile) { # job finished
			if (-M $togglefile < 5) {
				$rmflag = 0;
			}
		}else { # job is running or failed for some reasons
			opendir SUBDIR, $subdir;
			while (my $file = readdir SUBDIR) {
				my $fullname = $subdir."/".$file;
				if((-f $fullname) && (-M $fullname < 5)) {
					$rmflag = 0;
					last;
				}
			}
			closedir SUBDIR;
		}
		if ($rmflag) {
			rmtree($subdir);
		}
	}	
}
closedir DIR;
