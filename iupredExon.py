'''
iupredExon.py -- calculating the average structural disorder in each 
exon of a multi-exon transcript CSV
'''
import pandas as pd
import argparse
import numpy as np
import os
import math

def sequenceIupred(sequence, kind='short'):
	'''
	Given a sequence, calculate the IUPRED score of each residue.

	INPUT
		sequence : string, list of amino acids in the sequence
		kind : string, one of the three IUPRED options

	RETURNS
		list of lists, IUPRED score object

	'''
	temp_fasta = '_temp.fasta'
	temp_iupred = '_temp.iupred'

	o = open(temp_fasta,'w')
	o.write('>temp\n%s' % sequence)
	o.close()
	os.system("/Users/alecheckert/bin/iupred/iupred %s %s > %s" % (temp_fasta, kind, temp_iupred))
	scores = readIupred(temp_iupred)
	return scores

def readIupred(iupred_file):
	'''
	Reads output of the type given by IUPRED. *iupred_file* is a 
	file containing output exactly as given by the *iupred* command.

	'''	
	g = open(iupred_file, 'r')
	glines = g.read().split('\n')
	glines = [i for i in glines if len(i)>0]
	glines = [i for i in glines if i[0]!='#']
	glines = [i.split(' ') for i in glines]
	glines = [[tryConvert(j) for j in i if len(j)>0] for i in glines if len(i)>0]
	return glines

def tryConvert(arg):
	'''
	Converts *arg* to int or float if possible.

	'''
	try:
		if arg.isdigit():
			arg = int(arg)
			return arg
		elif all([i.isalnum() for i in arg.split('.')]):
			arg = float(arg)
			return arg
		else:
			return arg
	except (AttributeError, TypeError, ValueError) as e3:
		return arg

def transcriptSequence(transcript_file, peptide_column='protein', rank_column='rank'):
	'''
	Given a transcript file containing the protein sequences of 
	individual exons, returns the sequence of the entire protein.

	INPUT
		transcript_file : string, the name of the file containing the
			exon sequences, as CSV
		peptide_column : string, name of the column in the CSV
			corresponding to the peptide sequences
		rank_column : string, name of the column in the CSV corresponding
			to the exon rank

	RETURNS
		string, sequence of full protein

	'''
	f = pd.read_csv(transcript_file)
	f = f.sort_values(by=rank_column)
	f = f.set_index(rank_column, drop=False)
	protein = ''
	for i in f.index:
		if type(f.ix[i,peptide_column])==type(''):
			protein += f.ix[i,peptide_column]
	return protein

def iupredTranscript(transcript_file, peptide_column='protein', rank_column='rank', kind='short'):
	'''
	Given a transcript file containing the peptide sequences of individual
	exons, returns the IUPRED scores for the residues in the full
	transcript sequence.

	INPUT
		transcript_file : string, the name of the transcript file (CSV)
		peptide_column, rank_column : strings, the names of the
			columns in the CSV containing the peptide sequences and
			exon rank, respectively.
		kind : string, IUPRED calculation type

	RETURNS
		list of lists, IUPRED data object

	'''
	sequence = transcriptSequence(transcript_file, peptide_column=peptide_column, rank_column=rank_column)
	scores = sequenceIupred(sequence, kind=kind)
	score_df = pd.DataFrame(scores, columns=['residue', 'aa', 'residue_iupred'])
	return score_df

def pyiupred(transcript_file, peptide_column='protein', rank_column='rank', kind='short'):
	'''
	Takes a transcript file that contains exon peptide sequences, and adds a new
	column that is the average IUPRED score for the residues in that exon.

	INPUT
		transcript_file : string, name of the file containing the exon protein sequences
			in CSV format
		peptide_column : string, name of column in *transcript_file* that contains the
			exon sequences
		rank_column : string, name of column in *transcript_file* that contains the
			exon ranks
		kind : string, one of the three IUPRED types

	RETURNS
		Pandas DataFrame, the transcript information with the new transcript_iupred
			column

	'''
	scores = iupredTranscript(transcript_file, peptide_column=peptide_column, rank_column=rank_column, kind=kind)
	df = pd.read_csv(transcript_file)
	df = df.sort_values(by=rank_column)
	df = df.set_index(rank_column, drop=False)
	df.ix[1,'protein_start']=0
	if type(df.ix[1,'protein'])==type(''):
		df.ix[1,'protein_end']=len(df.ix[1,'protein'])
	else:
		df.ix[1,'protein_end']=0
	if len(df)==1:
		pass
	else:
		for i in range(2, len(df.index)+1):
			df.ix[i, 'protein_start'] = df.ix[i-1, 'protein_end']
			if type(df.ix[i,'protein'])==type(''):
				df.ix[i, 'protein_end'] = df.ix[i, 'protein_start'] + len(df.ix[i,'protein'])
			else:
				df.ix[i,'protein_end']=df.ix[i,'protein_start']
	for i in df.index:
		exon_scores = scores.ix[df.ix[i,'protein_start']:df.ix[i,'protein_end']-1]['residue_iupred'].values
		if len(exon_scores)==0:
			mean_iupred = float('NaN')
			stdev_iupred = float('NaN')
		else:
			mean_iupred = float(sum(exon_scores))/len(exon_scores)
			stdev_iupred = math.sqrt(float(sum([(mean_iupred-j)**2 for j in exon_scores]))/len(exon_scores))
		df.ix[i,'exon_iupred'] = mean_iupred
		df.ix[i,'exon_iupred_stdev'] = stdev_iupred
	df.to_csv(transcript_file, index=False)
	return df

def pyiupredFiles(transcript_files, peptide_column='protein', rank_column='rank', kind='kind'):
	'''
	Executes pyiupred() once for each of the files in *transcript_files*. See the
	pyiupred() function for documentation.

	'''
	for h in transcript_files:
		print "Calculating exon IUPRED scores for file %s..." % h
		try:
			pyiupred(h, peptide_column=peptide_column, rank_column=rank_column, kind=kind)
		except (KeyError, ValueError) as e2:
			print "Could not calculate IUPRED scores for file %s" % h
			continue
		

if __name__=='__main__':
	#x = pyiupred('ENST00000637978.csv', 'temp_outfile.csv')
	parser = argparse.ArgumentParser(description='Calculate IUPRED intrinsic disorder scores on individual exons in a transcript file')
	parser.add_argument('infile', type=str, help='name of transcript file or directory containing transcript files. Transcript files must be Pandas-style CSV with the sequences of the individual exons.')
	parser.add_argument('-k', '--kind', type=str, help='kind of IUPRED energy function to use. Default short.', default='short')
	parser.add_argument('--peptide_column', type=str, help='name of column containing exon peptide sequences in the transcript CSV', default='protein')
	parser.add_argument('--rank_column', type=str, help='name of column containing exon ranks in the transcript CSV', default='rank')

	args = parser.parse_args()
	if os.path.isdir(args.infile):
		transcript_files = ['%s/%s' % (args.infile, k) for k in os.listdir(args.infile)]
		pyiupredFiles(transcript_files, peptide_column=args.peptide_column, rank_column=args.rank_column, kind=args.kind)
		if os.path.isfile('_temp.fasta'):
			os.system('rm _temp.fasta')
		if os.path.isfile('_temp.iupred'):
			os.system('rm _temp.iupred')
	else:
		pyiupred(args.infile, peptide_column=args.peptide_column, rank_column=args.rank_column, kind=args.kind)
	if os.path.isfile('_temp.fasta'):
		os.system('rm _temp.fasta')
	if os.path.isfile('_temp.iupred'):
		os.system('rm _temp.iupred')
	print "Finished"
