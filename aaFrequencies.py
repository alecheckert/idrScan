'''
aaFrequencies.py -- calculate the frequency of each amino acid in a
string contained in a Pandas DataFrame
'''
import pandas as pd
import math
import argparse

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def letterFrequencies(word, name=None):
	'''
	Return a dictionary with the frequencies of each letter in a string.
	
	INPUT
		word : string, containing the aa's to analyze

	RETURNS
		pandas Series object containing frequencies of each letter

	'''
	if type(word)!=type(''):
		return pd.Series(name=name)
	letterdict = {}
	for letter in word:
		letterdict[letter] = letterdict.get(letter,0)+1
	for letter in letterdict:
		letterdict[letter] = float(letterdict[letter])/len(word)
	return pd.Series(letterdict, name=name)

def aaFrequencies(infile, target_column='protein'):
	'''
	Given a CSV with protein information (exons or transcripts), adds
	columns with the frequencies of each amino acid in the protein.

	INPUT
		infile : string, name of target CSV
		target_column : string, column in *infile* with the protein sequence

	RETURNS
		Pandas DataFrame with additional frequency information

	'''
	f = pd.read_csv(infile)
	freqs_df = pd.DataFrame()
	for i in f.index:
		freqs = letterFrequencies(f.ix[i,target_column], name=i)
		freqs_df = freqs_df.append(freqs)
		print "Calculated aa frequencies for %r" % i
	f = pd.concat([f,freqs_df], axis=1)
	f.to_csv(infile, index=False)
	return f

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='calculate the frequency of each letter in a string column')
	parser.add_argument('infile', type=str, help='CSV containing strings to be analyzed')
	parser.add_argument('-c', '--column', type=str, help="column in the CSV that contains the string to be analyzed. Default ``protein''", default='protein')
	args = parser.parse_args()
	aaFrequencies(args.infile, target_column=args.column)
