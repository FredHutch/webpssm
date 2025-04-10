package Path;

use strict;
use warnings;
use Carp qw[croak carp];
use Data::Dumper;


=head1 NAME

Common -- package for common parameters used in WebPSSM

=head1 SYNOPSIS


=head1 METHODS


=cut

our $documentroot = $ENV{'DOCUMENT_ROOT'};
our $uploadbase = "$documentroot/outputs";
our $statsbase = "$documentroot/stats";


1; #TRUE!!
 
