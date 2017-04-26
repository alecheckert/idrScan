'''
pairAlign.py -- search two DataFrames for proteins homologous to each other by
performing exhaustive pairwise local alignment
'''
import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo
import argparse

def topScore(seq_1, seq_2, opening_cost=-12, extension_cost=-2):
	'''
	Perform a local alignment and return the top scoring local alignment.
	Uses the BLOSUM62 matrix and an affine gap penalty.

	Note on the 
	
	INPUT
		seq_1, seq_2 : strings, the sequences to be aligned
		opening_cost : int, the cost to open a gap
		extension_cost : int, the cost to extend a gap

	RETURNS
		float, the score for the top-scoring alignment

	'''
	top_score = pairwise2.align.localds(seq_1, seq_2, MatrixInfo.blosum62, opening_cost, extension_cost, score_only=True)
	return top_score

def topScoreGlobal(seq_1, seq_2, opening_cost=-12, extension_cost=-2):
	'''
	Perform a global alignment.

	INPUT
		seq_1, seq_2: strings, the sequences to be aligned
	
	RETURNS
		float, the score for the global alignment

	'''
	top_score = pairwise2.align.globalds(seq_1, seq_2, MatrixInfo.blosum62, opening_cost, extension_cost, score_only=True)
	return top_score

def pairAlign(file_1, file_2, protein_col='protein', name_col='gene_id', outfile=None):
	'''
	For two Pandas DataFrame objects containing protein sequences, perform
	an exhaustive pairwise local alignment between all pairs of sequences.
	Return an ordered list with the scores of the top-scoring alignment for
	each pair.

	INPUT
		file_1 : CSV with Pandas DataFrame containing protein sequences
		file_2 : CSV with Pandas DataFrame containing protein sequences
		protein_col : string, name of the column with the protein sequences
			in df_1 and df_2
		name_col : string, column with the names of the protein sequences
		outfile : string, the file to write the results to

	RETURNS
		Pandas DataFrame with columns [seq_in_df_1, seq_in_df_2, top_score]

	'''
	df_1 = pd.read_csv(file_1)
	df_2 = pd.read_csv(file_2)
	results = []
	u = 0
	for i in df_1.index:
		for j in df_2.index:
			seq_i = df_1.ix[i,protein_col]
			seq_j = df_2.ix[j,protein_col]
			score = topScore(seq_i,seq_j)
			results.append((df_1.ix[i,name_col], df_2.ix[j,name_col], score))
			u+=1
			if u%50==0:
				print "%d sequences aligned..." % u
	results = np.asarray(results)
	return_frame = pd.DataFrame(results, columns=['mate_1', 'mate_2', 'score'])
	if outfile:
		return_frame.to_csv(outfile,index=False)
	return return_frame

def pairAlignGlobal(file_1, file_2, protein_col='protein', name_col='gene_id', outfile=None):
	'''
	Similar to pairAlign(), but performs a global alignment instead of a local alignment.

	INPUT
		file_1, file_2 : strings, CSVs containing the dataframes with protein sequences
		protein_col : string, the name of the column containing the protein sequences
		name_col : string, name of the column containing the sequence identities
		outfile : string, the file to write the results to

	RETURNS
		Pandas DataFrame with columns [seq_name_1, seq_name_2, global_alignment_score]

	'''
	df_1 = pd.read_csv(file_1)
	df_2 = pd.read_csv(file_2)
	results = []
	u = 0
	for i in df_1.index:
		for j in df_2.index:
			seq_i = df_1.ix[i,protein_col]
			seq_j = df_2.ix[j,protein_col]
			score = topScoreGlobal(seq_i, seq_j)
			results.append((df_1.ix[i,name_col], df_2.ix[j,name_col], score))
			u+=1
			if u%50==0:
				print '%d sequences aligned...' % u
	results = np.asarray(results)
	return_frame = pd.DataFrame(results, columns=['mate_1', 'mate_2', 'global_score'])
	if outfile:
		return_frame.to_csv(outfile,index=False)
	return return_frame

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='conduct local alignment between all pairs of sequences in two dataframes and return the top alignment scores')
	parser.add_argument('file_1', type=str, help='CSV containing the first set of protein sequences')
	parser.add_argument('file_2', type=str, help='CSV containing the second set of protein sequences')
	parser.add_argument('-o', '--outfile', type=str, help='name of the output file. Default is a synthesis of the input filenames', default=None)
	parser.add_argument('-c', '--column', type=str, help="name of the column in the CSVs containing the proteins sequences. Default ``protein''", default='protein')
	parser.add_argument('-n', '--name', type=str, help="name of the column in the CSVs containing the names of the sequences. Default ``gene_id''", default='gene_id')
	parser.add_argument('-g', '--global_align', action='store_true', help='perform a global alignment rather than a local alignment')

	args = parser.parse_args()
	if args.outfile==None:
		args.outfile = 'pairwise_alignment_%s_%s.csv' % (args.file_1.replace('.csv',''), args.file_2.replace('.csv', ''))
	print "Opening files..."
	try:
		if args.global_align:
			print "Performing global alignments..."
			pairAlignGlobal(args.file_1, args.file_2, protein_col=args.column, name_col=args.name, outfile=args.outfile)
		else:
			print "Performing local alignments..."
			pairAlign(args.file_1, args.file_2, protein_col=args.column, name_col=args.name, outfile=args.outfile)
		print "Finished"
	except (TypeError, pd.io.common.CParserError) as e_all:
		print "Could not read the input files."
		exit(1)
