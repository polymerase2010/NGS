#!/usr/bin/perl

use warnings; 
use strict;
open(IR,$ARGV[0])||die;
open(NEW,">$ARGV[0].prot")||die;
my $dna;
while(<IR>){
	chomp;
	if(!/^chr/){print NEW "$_\n";next;}
	my @line=split /\t/;
	if($line[14] eq 'synonymous' && $line[12]=~/(\d+)/){
		$dna=substr($line[47],6+($1%3),3);
		my $protein='';
		my $codon;
		for(my $i=0; $i<(length($dna)-2);$i+=3) {
			$codon=substr($dna,$i,3);
			$protein.=&codon2aa($codon);
		}
		my $p=int(($1-1)/3)+1;
		print $1%3,"\t$1\t$line[47]\t$dna\t$protein\n";
		$line[13]="p.$protein$p";
	}
	print NEW join("\t",@line),"\n";
}
close IR;
close NEW;
#*****************************************************************************************#
# codon2aa #
# A subroutine to translate a DNA 3-character codon to an amino acid
# Version 3, using hash lookup
sub codon2aa {
	 my($codon) = @_; $codon = uc $codon;
#uc=uppercase;lc=lowercase
#涔熷氨鏄ぇ灏忓啓杞崲锛寀c琛ㄧず灏嗘墍鏈夌殑灏忓啓 杞崲涓哄ぇ鍐?
#lc灏嗘墍鏈夌殑澶у啓杞崲涓哄皬鍐?
	my(%genetic_code) = ('TCA' => 'S', # Serine 
			'TCC' => 'S', # Serine 
			'TCG' => 'S', # Serine 
			'TCT' => 'S', # Serine 
			'TTC' => 'F', # Phenylalanine 
			'TTT' => 'F', # Phenylalanine 
			'TTA' => 'L', # Leucine 
			'TTG' => 'L', # Leucine 
			'TAC' => 'Y', # Tyrosine 
			'TAT' => 'Y', # Tyrosine 
			'TAA' => '_', # Stop
			'TAG' => '_', # Stop 
			'TGC' => 'C', # Cysteine 
			'TGT' => 'C', # Cysteine 
			'TGA' => '_', # Stop 
			'TGG' => 'W', # Tryptophan 
			'CTA' => 'L', # Leucine 
			'CTC' => 'L', # Leucine 
			'CTG' => 'L', # Leucine 
			'CTT' => 'L', # Leucine 
			'CCA' => 'P', # Proline 
			'CCC' => 'P', # Proline 
			'CCG' => 'P', # Proline 
			'CCT' => 'P', # Proline 
			'CAC' => 'H', # Histidine 
			'CAT' => 'H', # Histidine 
			'CAA' => 'Q', # Glutamine 
			'CAG' => 'Q', # Glutamine 
			'CGA' => 'R', # Arginine 
			'CGC' => 'R', # Arginine 
			'CGG' => 'R', # Arginine 
			'CGT' => 'R', # Arginine 
			'ATA' => 'I', # Isoleucine 
			'ATC' => 'I', # Isoleucine 
			'ATT' => 'I', # Isoleucine 
			'ATG' => 'M', # Methionine 
			'ACA' => 'T', # Threonine 
			'ACC' => 'T', # Threonine 
			'ACG' => 'T', # Threonine 
			'ACT' => 'T', # Threonine 
			'AAC' => 'N', # Asparagine 
			'AAT' => 'N', # Asparagine 
			'AAA' => 'K', # Lysine 
			'AAG' => 'K', # Lysine 
			'AGC' => 'S', # Serine 
			'AGT' => 'S', # Serine 
			'AGA' => 'R', # Arginine 
			'AGG' => 'R', # Arginine 
			'GTA' => 'V', # Valine 
			'GTC' => 'V', # Valine 
			'GTG' => 'V', # Valine 
			'GTT' => 'V', # Valine 
			'GCA' => 'A', # Alanine 
			'GCC' => 'A', # Alanine 
			'GCG' => 'A', # Alanine 
			'GCT' => 'A', # Alanine 
			'GAC' => 'D', # Aspartic Acid 
			'GAT' => 'D', # Aspartic Acid 
			'GAA' => 'E', # Glutamic Acid 
			'GAG' => 'E', # Glutamic Acid 
			'GGA' => 'G', # Glycine 
			'GGC' => 'G', # Glycine 
			'GGG' => 'G', # Glycine 
			'GGT' => 'G', # Glycine 
			);
	if(exists $genetic_code{$codon}) { 
		return $genetic_code{$codon}; 
	} else { 
		print STDERR "Bad codon \"$codon\"!!\n"; exit; 
	}
}
 #***************************************************************************************** #
