'''
localEntropy.py -- calculate the mean local sequence entropy of 
various protein sequences
'''
import pandas as pd
import numpy as np
import argparse
import math

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def entropy(sequence):
	'''
	Calculate the entropy of a protein sequence.

	INPUT
		sequence : string, protein sequence

	RETURNS
		float, entropy

	'''
	count_dict = {i: 0 for i in amino_acids}
	for aa in sequence:
		try:
			count_dict[aa] += 1
		except KeyError:
			continue
	count_list = [count_dict[i] for i in count_dict]
	count_sum = sum(count_list)
	count_list = [float(i)/count_sum for i in count_list if i!=0]
	return abs(sum([i*math.log(i) for i in count_list]))

def localEntropy(sequence, window=20, consider_short_sequences=False):
	'''
	Calculate the mean sequence entropy within a defined window size
	across all subsequences of a protein sequence.

	INPUT
		sequence : string, protein sequence
		window : int, length of subsequences to consider
		consider_short_sequences : bool, whether or not to calculate
			the entropy of sequences that are shorter than *window*.
			If *False*, returns *None* when it encounters such sequences.

	RETURNS
		float, mean local entropy

	'''
	if len(sequence)<window and consider_short_sequences:
		return entropy(sequence)
	elif len(sequence)<window and not consider_short_sequences:
		return False
	else:
		entropies = []
		for i in range(len(sequence)-window+1):
			subsequence = sequence[i:(i+window)]
			entropies.append(entropy(subsequence))
		return float(sum(entropies))/len(entropies)

def runLocalEntropy(filename, window=20, protein_column='protein', consider_short_sequences=False):
	'''
	Run the local entropy algorithm on a set of protein sequences in a 
	Pandas-style CSV file.

	INPUT
		filename : string, name of file containing protein sequences. Must
			be in Pandas CSV format.
		window : int, length of subsequences for local entropy calculation
		protein_column : string, name of column in *filename* containing
			protein sequences
		consider_short_sequences : bool, whether or not to return entropy
			scores for proteins that are shorter than *window*

	RETURNS
		Pandas DataFrame object, updated contents of *filename* containing
			the ``local_entropy'' column.

	'''
	f = pd.read_csv(filename)
	for i in f.index:
		try:
			f.ix[i,'local_entropy'] = localEntropy(f.ix[i,'protein'])
		except TypeError:
			f.ix[i,'local_entropy'] = False
		#print f.ix[i,'description'], f.ix[i,'protein'], f.ix[i,'local_entropy']
	f.to_csv(filename,index=False)
	return f

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='calculate mean local entropy over subsequences of proteins')
	parser.add_argument('infile', type=str, help='CSV containing protein sequences')
	parser.add_argument('-w', '--window', type=int, help='length of subsequences for the entropy calculation. Default 20.', default=20)
	parser.add_argument('-c', '--column', type=str, help='name of protein column in *filename* containing the protein sequences. Default "protein".', default='protein')
	parser.add_argument('-s', '--short', action='store_true', help='whether or not to calculate entropy for sequences shorter than *window*', default=False)
	args = parser.parse_args()
	print "Calculating local entropy..."
	runLocalEntropy(args.infile, window=args.window, protein_column=args.column, consider_short_sequences=args.short)
	print "Finished"
