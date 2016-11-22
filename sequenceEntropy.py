'''
sequenceEntropy.py -- calculate the entropy of protein sequences
'''
import pandas as pd
import math
import argparse
import os

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def entropy(sequence, discard_non_amino_acid=False):
	'''
	Calculate the entropy of a peptide sequence.

	INPUT
		sequence : string, a peptide sequence. Note that the program does not consider non-amino
			acid characters present in the sequence.

	RETURNS
		float : sequence entropy

	'''
	if type(sequence)!=type(''):
		return float('NaN')
	else:
		freqs = [float(sum([1 for i in sequence if i==aa]))/len(sequence) for aa in amino_acids]
		return -1 * sum([i*math.log(i) for i in freqs if i>0])

def sequenceEntropy(filename, peptide_column='protein', write_to_file=False, outname=None):
	'''
	Given a CSV encoding a set of protein sequences, adds an ``entropy'' column
	with the sequence entropy.

	INPUT 
		filename : string, name of Pandas-style CSV containing the protein information
		peptide_column : name of the column in *filename* that contains the protein sequence
		write_to_file : bool, whether to write the information to an output file
		outname : string, name of file to write to

	RETURNS
		Pandas DataFrame object, with the added 'protein' column

	'''
	df = pd.read_csv(filename)
	if peptide_column not in df.columns:
		print "Did not find the column %s in dataframe %s" % (peptide_column, filename)
		exit(1)
	entropies = pd.Series([])
	for i in df.index:
		entropies.ix[i] = entropy(df.ix[i,peptide_column])
	df['entropy']=entropies
	if write_to_file:
		if not outname:
			outname='%s_entropy.csv' % filename.replace('.csv','')
		df.to_csv(outname,index=False)
	return df

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='calculate the entropy of peptide sequences contained in a CSV')
	parser.add_argument('infile', type=str, help='name of Pandas-style CSV with peptide sequences, or name of directory containing CSV files')
	parser.add_argument('-c', '--column', type=str, help='name of peptide-containing column in the CSV', default='protein')
	args = parser.parse_args()
	if os.path.isdir(args.infile):
		fs = [i for i in os.listdir(args.infile) if '.csv' in i]
		for f in fs:
			try:
				sequenceEntropy('%s/%s' % (args.infile, f), peptide_column=args.column, write_to_file=True, outname='%s/%s' % (args.infile, f))
				print 'successfully calculated entropy for %s' % f
			except pandas.io.common.CParserError:
				continue
	else:
		sequenceEntropy(args.infile, peptide_column=args.column, write_to_file=True, outname=args.infile)
	print "Finished"
