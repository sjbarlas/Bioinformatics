use strict;
use warnings;

# Search for the longest open reading frame for this DNA.
print "\nHere is the largest ORF, from 5' to 3':\n" ;
local $_ = $RNA_seq ;
while ( / AUG /g ) {
    my $start = pos () - 2 ;
    if ( / UGA|UAA|UAG /g ) {
        my $stop = pos ;
        $gene = substr ( $_ , $start - 1 , $stop - $start + 1 ), $/ ;
        print "$gene" ;
    }
}

# The next set of commands translates the ORF found above for an amino acid seq.
print "\nThe largest reading Frame is:\t\t\t" . $protein { "gene" } . "\n" ;
sub translate {
    my ( $gene , $reading_frame ) = @_ ;
    my %protein = ();
    for ( $i = $reading_frame ; $i < length ( $gene ); $i += 3 ) {
        $codon = substr ( $gene , $i , 3 );
        $amino_acid = translate_codon( $codon );
        $protein { $amino_acid }++;
        $protein { "gene" } .= $amino_acid ;
    }
    return %protein ;
}

sub translate_codon {
if ( $_ [ 0 ] =~ / GC[AGCU] /i )             { return 'A';} # Alanine;
if ( $_ [ 0 ] =~ / UGC|UGU /i )              { return 'C';} # Cysteine
if ( $_ [ 0 ] =~ / GAC|GAU /i )              { return 'D';} # Aspartic Acid;
if ( $_ [ 0 ] =~ / GAA|GAG /i )              { return 'Q';} # Glutamine;
if ( $_ [ 0 ] =~ / UUC|UUU /i )              { return 'F';} # Phenylalanine;
if ( $_ [ 0 ] =~ / GG[AGCU] /i )             { return 'G';} # Glycine;
if ( $_ [ 0 ] =~ / CAC|CAU /i )              { return 'His';} # Histine (start codon);
if ( $_ [ 0 ] =~ / AU[AUC] /i )              { return 'I';} # Isoleucine;
if ( $_ [ 0 ] =~ / AAA|AAG /i )              { return 'K';} # Lysine;
if ( $_ [ 0 ] =~ / UUA|UUG|CU[AGCU] /i )     { return 'Leu';} # Leucine;
if ( $_ [ 0 ] =~ / AUG /i )                  { return 'M';} # Methionine;
if ( $_ [ 0 ] =~ / AAC|AAU /i )              { return 'N';} # Asparagine;
if ( $_ [ 0 ] =~ / CC[AGCU] /i )             { return 'P';} # Proline;
if ( $_ [ 0 ] =~ / CAA|CAG /i )              { return 'G';} # Glutamine;
if ( $_ [ 0 ] =~ / AGA|AGG|CG[AGCU] /i )     { return 'R';} # Arginine;
if ( $_ [ 0 ] =~ / AGC|AGU|UC[AGCU] /i )     { return 'S';} # Serine;
if ( $_ [ 0 ] =~ / AC[AGCU] /i )             { return 'T';} # Threonine;
if ( $_ [ 0 ] =~ / GU[AGCU] /i )             { return 'V';} # Valine;
if ( $_ [ 0 ] =~ / UGG /i )                  { return 'W';} # Tryptophan;
if ( $_ [ 0 ] =~ / UAC|UAU /i )              { return 'Y';} # Tyrosine;
if ( $_ [ 0 ] =~ / UAA|UGA|UAG /i )          { return "***" ;} # Stop Codons;
}
Am I missing something?
