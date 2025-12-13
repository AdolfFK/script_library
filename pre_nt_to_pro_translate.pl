#!/bin/perl 
##########################################################
#
# This program translates nucleic acid into protein sequence
# Mail wangjyafk@126.com
# Version 1.0
##########################################################
use strict; 
use warnings;
use File::Basename;

my $ARGVlength=@ARGV;
my @usage = qq(
===================== pre_nt_to_pro_translate.pl =====================
                            Version 1.0

Usage:   perl pre_nt_to_pro_translate.pl \
         .../mutation.statistics file \  
         .../output list 
         
)or die $!;
#Example:perl pre_nt_to_pro_translate.pl rsva-f.blast.pre.mafft.uniq.fasta rsva-f
###check option


if (($ARGVlength == 0) or ($ARGV[0] eq "-h")) 
  {
    print join("\n",@usage)."\n";
    exit;
  }

### Parameter assignment 
my $arrayfile = "$ARGV[0]";
my $mutationout = "$ARGV[1]";
if ( !(-e $mutationout) ) {
  #print "NO";
  mkdir( $mutationout ) or die "Can not creat shelloutdir directory, $!";
}

###############
## main program
###############
open(IN1,"$arrayfile") or die("can not open the file!");
open(OUT1,">$mutationout/$arrayfile.1.fasta");
open(OUT2,">$mutationout/$arrayfile.2.fasta");
open(OUT3,">$mutationout/$arrayfile.3.fasta");
open(OUT4,">$mutationout/$arrayfile.4.fasta");

my %hash;
my $k;
my $seq;
while (<IN1>) {
	chomp;
	$_ =~ s/[\r\n]+\Z//g;
	if ($_ =~ /^>/) {
		$k = $_;
	}else{
		$seq = $_;
		$hash{$k} .= $seq;
	}
}
close IN1;

foreach my $id (keys %hash){
	my $seq = $hash{$id};
	my $seq1 = $seq;
	my @temp1 = split //,$seq1;
	my $protein1 = proteinseq(@temp1);
	my $count1 = () = $protein1 =~ /-/g;
	my @matches1 = $protein1 =~ /(M[^-]+)/g;
	my ($max1,$len1) = getmaxstr(@matches1);

	my $seq2 = substr($seq1,1);
	my @temp2 = split //,$seq2;
	my $protein2 = proteinseq(@temp2);
	my $count2 = () = $protein2 =~ /-/g;
	my @matches2 = $protein2 =~ /(M[^-]+)/g;
	my ($max2,$len2) = getmaxstr(@matches2);

	my $seq3 = substr($seq1,2);
	my @temp3 = split //,$seq3;
	my $protein3 = proteinseq(@temp3);
	my $count3 = () = $protein3 =~ /-/g;
	my @matches3 = $protein3 =~ /(M[^-]+)/g;
	my ($max3,$len3) = getmaxstr(@matches3);

	my $seq4 = reverse $seq;
	my @temp4 = split //,$seq4;
	my $protein4 = proteinseq(@temp4);
	my $count4 = () = $protein4 =~ /-/g;
	my @matches4 = $protein4 =~ /(M[^-]+)/g;
	my ($max4,$len4) = getmaxstr(@matches4);

	my $seq5 = substr($seq4,1);
	my @temp5 = split //,$seq5;
	my $protein5 = proteinseq(@temp5);
	my $count5 = () = $protein5 =~ /-/g;
	my @matches5 = $protein5 =~ /(M[^-]+)/g;
	my ($max5,$len5) = getmaxstr(@matches5);

	my $seq6 = substr($seq5,2);
	my @temp6 = split //,$seq6;
	my $protein6 = proteinseq(@temp6);
	my $count6 = () = $protein6 =~ /-/g;
	my @matches6 = $protein6 =~ /(M[^-]+)/g;
	my ($max6,$len6) = getmaxstr(@matches6);

	my $mincount = $count1;
	my $minmax = $max1;
	my $minlength = $len1;
	my @arrayall = ($count1, $count2, $count3, $count4, $count5, $count6);
	my @arraymax = ($max1, $max2, $max3, $max4, $max5, $max6);
	my @arraylen = ($len1, $len2, $len3, $len4, $len5, $len6);


	for (my $i = 0; $i <= $#arrayall; $i++) {
		if ($arrayall[$i] <= $mincount) {
			$mincount = $arrayall[$i];
			$minmax = $arraymax[$i];
			$minlength = $arraylen[$i];
		}
		print OUT2 "$arrayall[$i]\t";
		print OUT3 "$arraymax[$i]\t";
		print OUT4 "$arraylen[$i]\t";
	}
	print OUT2 "$id\n";
	print OUT3 "$id\n";
	print OUT4 "$id\n";
	
	print OUT1 "$id\n$minmax\n";
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;

sub getmaxstr{
	my @matches = @_;
	my $max_length = 0;
	my $max_substring = "";

	foreach my $match (@matches) {
	    my $length = length($match);
	    if ($length > $max_length) {
	        $max_length = $length;
	        $max_substring = $match;
	    }
	}
	return ($max_substring,$max_length);
}


sub proteinseq {
	my @tmp = @_;
	my $protein = '';
	# print $#tmp,"\t";
	my $seq = join(" ", @tmp);
	$seq =~ s/ //g;
	# print $seq,"\t";
	for(my $i=0; $i<(length($seq)-2);$i+=3)
	{
		my $codon=substr($seq,$i,3);
		$protein.=codon2aa($codon);
	}	

	return $protein;
}


#***************************************************************************************************
# codon2aa   
#   
# A subroutine to translate a DNA 3-character codon to an amino acid Version 3, using hash lookup
#   
#***************************************************************************************************  
sub codon2aa   
{   
    my($codon) = @_;   
   
    $codon = uc $codon;#uc=uppercase;lc=lowercase  
                   #That is, case conversion, uc means converting all lowercase to uppercase  
               	   #'lc' converts all uppercase to lowercase 
    
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
    'TAA' => '-',    # Stop   
    'TAG' => '-',    # Stop   
    'TGC' => 'C',    # Cysteine   
    'TGT' => 'C',    # Cysteine   
    'TGA' => '-',    # Stop   
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
    'ATG' => 'M',    # Methionine   
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
   
    if(exists $genetic_code{$codon})   
    {   
        return $genetic_code{$codon};   
    }  
    else  
    {   
   		return "X";
        # print STDERR "Bad codon \"$codon\"!!\n";   
        # exit;   
    }   
} 
