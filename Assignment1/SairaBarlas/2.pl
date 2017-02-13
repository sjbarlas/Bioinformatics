use strict;
use warnings;
use autodie;

open my $fh, '<', 'file.fasta';           # Open filehandle in read mode

while ( my $line = <$fh> ) {              # Loop over line by line

    print $line                           # Print line if it matches pattern
      if $line =~ /h..hc.c/;              # '.' in a regular expression matches
                                          # (almost) anything
}

close $fh; 