#!/usr/bin/perl

my $DNA= $ARGV[0];

@dna = $DNA;
print "@dna[0]\t\t-> DNA\n";


for ($i = 0 ; $i < length(@dna[0]); $i++ )  { 
    @c[$i] = substr( @dna[0], $i, 1);
    #print "@c[$i]";
    
    if(@c[$i] eq 'A') { 
        $counterpart=$counterpart."T";
        $rna=$rna."A";
    }
    
    elsif(@c[$i] eq 'T') { 
        $counterpart=$counterpart."A";
        $rna=$rna."U";
    }
        
    elsif(@c[$i] eq 'C') { 
        $counterpart=$counterpart."G";
        $rna=$rna."C";
        
    }
    
    elsif(@c[$i] eq 'G') { 
        $counterpart=$counterpart."C";
        $rna=$rna."G";
    }
}

print "$counterpart\t\t-> Counterpart of DNA\n";
print "$rna\t\t-> RNA\n";

my $j = 0;
for ($i = 0 ; $i < length ($rna); $i = $i+3) {
    @codon[$j] = substr( $rna, $i, 3);
    print "@codon[$j] ";
    $j++;
}

print "-> Codon\n";  

my(%genetic_code_table) = (
   
    'UUU' => 'F',    
    'UUC' => 'F',
    'UUA' => 'L',
    'UUG' => 'L',
    'CUU' => 'L',
    'CUC' => 'L',
    'CUA' => 'L',
    'CUG' => 'L',
    'AUU' => 'I',
    'AUC' => 'I',
    'AUA' => 'I',
    'AUG' => 'M',
    'GUU' => 'V',
    'GUC' => 'V',
    'GUA' => 'V',
    'GUG' => 'V',
    'UCU' => 'S',
    'UCC' => 'S',
    'UCA' => 'S',
    'UCG' => 'S',
    'CCU' => 'P',
    'CCC' => 'P',
    'CCA' => 'P',
    'CCG' => 'P',
    'ACU' => 'T',
    'ACC' => 'T',
    'ACA' => 'T',
    'ACG' => 'T',
    'GCU' => 'A',
    'GCC' => 'A',
    'GCA' => 'A',
    'GCG' => 'A',
    'UAU' => 'Y',
    'UAC' => 'Y',
    'UAA' => '.',
    'UAG' => '.',
    'CAU' => 'H',
    'CAC' => 'H',
    'CAA' => 'Q',
    'CAG' => 'Q',
    'AAU' => 'N',
    'AAC' => 'N',
    'AAA' => 'K',
    'AAG' => 'K',
    'GAU' => 'D',
    'GAC' => 'D',
    'GAA' => 'E',
    'GAG' => 'E',
    'UGU' => 'C',
    'UGC' => 'C',
    'UGA' => '.',
    'UGG' => 'W',
    'CGU' => 'R',
    'CGC' => 'R',
    'CGA' => 'R',
    'CGG' => 'R',
    'AGU' => 'S',
    'AGC' => 'S',
    'AGA' => 'R',
    'AGG' => 'R',
    'GGU' => 'G',
    'GGC' => 'G',
    'GGA' => 'G',
    'GGG' => 'G',
);

for (my $j = 0 ; $j < length ($rna)/3; $j++) {
     
    if( exists ($genetic_code_table{@codon[$j]}) ) {
        @result[$j] = $genetic_code_table{@codon[$j]};
        print "$genetic_code_table{@codon[$j]}   ";
    }
}

open(RESULT, "+>result.txt") or die "Couldn't open file result.txt, $!";

print RESULT @result;

close(RESULT) || die "Couldn't close file properly";
print "\n\n";
