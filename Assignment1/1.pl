!#/usr/bin/perl -w
use strict;
use warnings;
while (my $line = <DATA>) {
    chomp $line;
    print "processing $line\n";
    my ($id, $rna_sq) = split(/\s+/, $line);

    while ($rna_sq =~ /atg/g) {
        # $start and $stop are 0-based indexes
        my $start = pos($rna_sq) - 3; # back up to include the start sequence

        # discard if no stop sequence can be found
        last unless $rna_sq =~ /tga|taa|tag/g;

        my $stop    = pos($rna_sq);
        my $genelen = $stop - $start;

        my $gene    = substr($rna_sq, $start, $genelen);
        print "\t" . join(' ', $id, $start+1, $stop, $gene, $genelen) . "\n";
    }
}
