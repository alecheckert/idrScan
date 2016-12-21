'''
lcCalculator.py -- determine the local bias of individual protein sequences
for subsets of the amino acids
'''
import pandas as pd
import numpy as np
import argparse
import math

sample_sequence = 'GAFWGLVFGLGVGLLRMILEFSYPAPACGEVDRRPAVLKDFHYLYFAILLCGLTAIVIVIVSLCTTPIPEEQ'
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'Z']

def makeMasker(target_aa):
    '''
    Make a sequence mask object, which is helpful in making binary
    masks for protein sequences.

    INPUT
        target_aa : list of char

    RETURNS
        dict, protein sequence masker

    '''
    masker = {}
    for aa in amino_acids:
        if aa in target_aa:
            masker[aa] = 1
        else:
            masker[aa] = 0
    return masker

def aaCount(sequence, masker):
	return sum([masker[i] for i in sequence])

def highestCount(sequence, masker, window_size=30):
	'''
	For a protein sequence, find the subsequence with the highest bias
	toward a subset of the amino acids. Then returnt the count of amino
	acids in that subsequence.
	
	INPUT
		sequence : string, protein sequence
		masker : dict, protein sequence masker object
		window_size : int, length of subsequence

	RETURNS
		int, frequency of the set of amino acids in the subsequence
			with the highest frequency

	'''
	x = [masker[i] for i in sequence]
	current = (None, None, None, 0)
	if len(sequence)<window_size:
		return (0, len(sequence), sequence, sum(x))
	for start in range(len(sequence)-window_size):
		count = sum(x[start:start+window_size])
		if count > current[3]:
			current = (start, start+window_size, sequence[start:start+window_size], count) 
	return current

def lcCalculator(sequence_file, target_aa, window_size=30, protein_column='protein', write_in_place=True, target_aa_name=None):
	'''
	Look at a set of protein sequences, find the subregions in each that
	are most biased toward a subset of the amino acids, then for each
	protein record the count of those amino acids in those subregions.

	INPUT
		sequence_file : string, name of Pandas CSV containing protein
			sequences
		target_aa : list of char, subset of amino acids to consider
		window_size : int, length of subsequences
		protein_column : string, name of column in *sequence_file* 
			containing protein sequences
		write_in_place : bool, whether or not to write to the file
		target_aa_name : string, name of the subset of amino acids. This
			becomes the name of the new column in the output dataframe.

	RETURNS
		Pandas DataFrame, contents of file containing the new counts
		
	'''
	f = pd.read_csv(sequence_file)
	masker = makeMasker(target_aa)
	if not protein_column in f.columns:
		print "Could not find the column ``%s'' in %s." % (protein_column, sequence_file)
		exit(1)
	aa_string = ''
	if target_aa_name:
		aa_string = target_aa_name
	else:
		for aa in target_aa:
			aa_string += '%s_' % aa
	peak_column = '%speak' % aa_string
	for i in f.index:
		if (type(f.ix[i,protein_column]!=type('')) and type(f.ix[i,protein_column])!=str) or (len(f.ix[i,protein_column])==0):
			print type(f.ix[i,protein_column])
			f.ix[i,peak_column] = 0
		else:
			start, end, subseq, count = highestCount(f.ix[i,protein_column], masker=masker, window_size=window_size)
			print subseq, count ##just for fun
			f.ix[i,peak_column] = count
	if write_in_place:
		f.to_csv(sequence_file, index=False)
	return f

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Look at a set of protein sequences. For each, find the subsequence that is most biased toward a subset of the amino acids, then record the count of amino acids in that subsequence.')
	parser.add_argument('infile', type=str, help='name of Pandas CSV containing protein sequences')
	parser.add_argument('-a' '--aminoacid', type=str, action='append', dest='target_amino_acids', help='amino acids to look for')
	parser.add_argument('-c', '--column', type=str, help='name of column in the infile containing protein sequences. Default "protein".', default='protein')
	parser.add_argument('-w', '--window', type=int, help='length of subsequences to consider. Default 30.', default=30)
	parser.add_argument('-p', '--write_in_place', action='store_true', help='write results to the input file', default=False)
	parser.add_argument('-o', '--output_column', type=str, help='name of the column in the output dataframe to write to. Default None, meaning that output column is set to the target amino acid names.', default=None)

	args = parser.parse_args()
	lcCalculator(args.infile, target_aa=args.target_amino_acids, window_size=args.window, protein_column=args.column, write_in_place=args.write_in_place, target_aa_name=args.output_column)
