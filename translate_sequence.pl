#!/usr/bin/perl -w


#   The following perl script will:
# 	ask the user for a file and extract out the DNA sequences, if any, contain in it
#   Translate a DNA sequence into its amino acid sequence for all six reading frames
#   Display the amino acid sequence to the screen in an appropriate format
#   print the amino acid sequence to a word document, chosen by the user, in the form of a string  

# require to declare variable explicitly 
use strict;
use warnings;




	# Initialize variables 
	# to store the DNA sequence in an array; to store the DNA sequence as a string 
	#the reverse compliment of this sequence and the translated protein sequence
my @file_data = (  );
my $dna = '';
my $revcom = '';
my $protein = '';


# Enter the file name; pass the file name to a function which returns the DNA sequence as an array  
	print "enter the name of the sequence file : ";
	my $FileName = <>;
	chomp $FileName;
	@file_data = get_file_data($FileName);


# Extract the sequence data from the contents of the file "sample.dna"
	$dna = extract_sequence_from_fasta_data(@file_data);


# enter the name of the file to output the results 
	print "enter the name of the file to output results: ";
	my $FileName_Out = $FileName."translation.doc";



# Translate the DNA to protein in six reading frames
#   and print the protein in lines 60 characters long
	print "\n -------Reading Frame 1--------\n\n";
	$protein = translate_frame($dna, 1);
	print_sequence($protein, 60);
	print_sequence_file($protein, 60, "1 5' to 3'", $FileName_Out);


	print "\n -------Reading Frame 2--------\n\n";
	$protein = translate_frame($dna, 2);
	print_sequence($protein, 60);
	print_sequence_file($protein, 60, "2 5' to 3'", $FileName_Out);


	print "\n -------Reading Frame 3--------\n\n";
	$protein = translate_frame($dna, 3);
	print_sequence($protein, 60);
	print_sequence_file($protein, 60, "3 5' to 3'", $FileName_Out);




# get the reverse complement of the primary DNA sequence
	$revcom = revcom($dna);



# print the amino acid sequence for the frames of the reverse compliment
	
	print "\n -------Reading Frame 4--------\n\n";
	$protein = translate_frame($revcom, 1);
	print_sequence($protein, 60);
	print_sequence_file($protein, 60, "4 3' to 5'", $FileName_Out);


	print "\n -------Reading Frame 5--------\n\n";
	$protein = translate_frame($revcom, 2);
	print_sequence($protein, 60);
	print_sequence_file($protein, 60, "5 3' to 5'", $FileName_Out);


	print "\n -------Reading Frame 6--------\n\n";
	$protein = translate_frame($revcom, 3);
	print_sequence($protein, 60);	
	print_sequence_file($protein, 60, "6 3' to 5'", $FileName_Out);

exit;






#  ****************************  functions of the perl script********************************

# From Chapter 8

#
# codon2aa
#
# 





#********************************* codon2aa ******************************************

	# A subroutine to translate a DNA, 3-character sequence, a codon to an amino acid using hash table lookup
	# in this function the codons are extracted from the sequence using the substring function
	# in could just as well be done in other ways such as the  unpack  function discussed in class 
	# (adapted from From Chapter 8 mastering perl in bioinformatics)	
	 

# 	pass a dna sequence 
#   
#   
#	
#   returns the translated amino acid sequence as a string 

#********************************************************************************************


sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '#',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine  using a space a it is easier to see in a quick inspection
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{

            print STDERR "Bad codon \"$codon\"!!\n";
            exit;
    }
}


 


#********************************* extract_sequence_from_fasta_data ******************************************

	# A subroutine to extract FASTA sequence data from an array
	# (adapted from From Chapter 8 mastering perl in bioinformatics)	
	 

# 	pass array containing the contents of the file 
#   
#   
#	
#   returns only the DNA sequences

#********************************************************************************************

sub extract_sequence_from_fasta_data {

    my(@fasta_file_data) = @_;

    use strict;
    use warnings;

    # Declare and initialize variables
    my $sequence = '';

    foreach my $line (@fasta_file_data) {

        # discard blank line
        if ($line =~ /^\s*$/) {
            next;

        # discard comment line
        } elsif($line =~ /^\s*#/) {
            next;

        # discard fasta header line
        } elsif($line =~ /^>/) {
            next;

        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }

    # remove non-sequence data (in this case, whitespace) from $sequence string
    $sequence =~ s/\s//g;

    return $sequence;
}






#********************************* get_file_data ******************************************

	# A Subroutine to Read FASTA Files and to get data from a file
	# each line will contain the number of characters passed to the function ($length)
	# (adapted from From Chapter 8 mastering perl in bioinformatics)	
	 

# 	pass the file name of the FASTA file
#   
#   
#	
#   returns all the data in the file as an array

#********************************************************************************************

sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}





#********************************* print_sequence ******************************************

	# A subroutine to format and print sequence data
	# each line will contain the number of characters passed to the function ($length)
	# (adapted from From Chapter 8 mastering perl in bioinformatics)	
	 

# 	pass amino acid sequence
#   pass the length of the sequence
#   
#	
#   returns nothing

#********************************************************************************************

sub print_sequence {

    my($sequence, $length) = @_;

    use strict;
    use warnings;

    # Print sequence in lines of $length
   for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ) {
        print substr($sequence, $pos, $length), "\n\n";
   
	 
   
   
    }
}




#********************************* print_sequence_file ******************************************

	# A subroutine to output the translated sequence to a file 
	 

# 	pass amino acid sequence
#   pass the frame number 
#   pass the name of the output file 
#	
#   returns nothing

#********************************************************************************************

sub print_sequence_file {

    my($sequence, $length, $number, $FileName) = @_;

    use strict;
    use warnings;
	
	#enter the name of the file passed to the sub routine and open it for writing
	
	unless( open(MYFILE, ">>", $FileName) ) {
        print STDERR "Cannot open file \"filename\"\n\n";
        exit;
    }
	
	
	# append the Amino acid sequence to the file paases to the function ($filename) 
	
	print MYFILE "\n\n----------------------- Reading Frame ******** $number -----------------------------\n\n";
	
	print MYFILE "$sequence\n\n";
	
	print MYFILE "\n\n------------------- end of frame translation ------------------------------\n\n";
	close MYFILE 
}









#*********************************  dna2peptide ******************************************

	# A subroutine to translated (convert) a DNA sequence into its equivalent amino acid sequence
	# (adapted from From Chapter 8 mastering perl in bioinformatics) 

# 	pass the DNA sequence
#	
#   returns the translated protein sequence as a string

#********************************************************************************************

sub dna2peptide {

    my($dna) = @_;

    use strict;
    use warnings;
  
    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return $protein;
}


#*********************************  translate_frame ******************************************

	# A subroutine to translate a frame of DNA
	# (adapted from From Chapter 8 mastering perl in bioinformatics) 

# 	pass the DNA sequence
#	the starting position of the reading frame
#	the ending position of the reading frame 
#
#   returns the translated sequence as a amino acid sequence as a string

#********************************************************************************************



sub translate_frame {

    my($seq, $start, $end) = @_;

    my $protein;

    # To make the subroutine easier to use, you won't need to specify
    #  the end point-it will just go to the end of the sequence
    #  by default.
    unless($end) {
        $end = length($seq);
    }

    # Finally, calculate and return the translation
        return dna2peptide ( substr ( $seq, $start - 1, $end -$start + 1) );
}


 

#*********************************  revcom ******************************************

	# A subroutine to compute the reverse complement of DNA sequence
	# (adapted from From Chapter 8 mastering perl in bioinformatics) 

# 	pass the DNA sequence
#	
#   returns the reverse compliment of the DNA sequence as a string

#********************************************************************************************


sub revcom {

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}
