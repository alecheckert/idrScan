'''
findDisorderedRegions.py -- find discrete regions with high IUPred scores in
a CSV of protein sequences.

Requires the ``iupred'' command on the shell path. For Unix-like systems:

	1. Download ``iupred'' from https://iupred2a.elte.hu/download
	2. Expand and place in suitable location (e.g. ~/bin)
	3. Compile ``iupred'' (preferably with gcc)
	4. Append path to iupred folder to environment variables
	   PATH and IUPred_PATH

'''
__author__ = 'Alec Heckert'
import sys
import pandas as pd
import argparse
import os
import math
from io import StringIO
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage.filters import gaussian_filter
sns.set_style('white')

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def wrapup(outname):
	plt.tight_layout()
	plt.savefig(outname, dpi=400)
	plt.close()
	os.system('open %s' % outname)

def sequenceIupred(sequence, kind='short'):
	'''
	Given a sequence, calculate the IUPred score of each residue.
	Assumes ``iupred'' command is in shell path.

	INPUT
		sequence	:	str, single letter amino acid sequence
		kind		:	str, `short', `long', or `glob'

	RETURNS
		pandas.DataFrame with columns
			`pos'	:	int, index in sequence
			`AA'	:	char, residue
			`score'	:	float, the IUPred score for that residue

	'''
	temp_fasta, temp_iupred = '_temp.fasta', '_temp.iupred'
	with open(temp_fasta, 'w') as outf:
		outf.write('>temp\n%s' % sequence)
	os.system("iupred %s %s > %s" % (temp_fasta, kind, temp_iupred))
	scores = readIupred(temp_iupred)
	os.system('rm %s %s' % (temp_fasta, temp_iupred))
	return scores

def readIupred(iupred_file):
	'''
	Reads output of ``iupred'' command.

	INPUT
		iupred_file		:	str, file containing IUPred output
	
	RETURNS
		pandas.DataFrame with columns
			`pos'	:	int, index in sequence
			`AA'	:	char, residue
			`score'	:	float, the IUPred score for that residue

	'''
	with open(iupred_file, 'r') as g:
		fstring = StringIO('\n'.join(['pos AA score'] + [i for i in g.read().split('\n') if (len(i)>0) and (i[0]!='#')]))
	return pd.read_csv(fstring, delim_whitespace=True, index_col=0)

def iupredTrace(protein_file, kind='short', protein_column='protein', id_column='name'):
	'''
	Performs ``iupred'' on each sequence in a CSV, saving output
	to individual trace files.

	INPUT
		protein_file	:	str, CSV name with protein sequences
		kind			:	str, `short', `long', or `glob'
		protein_column	:	str, name of the column containing protein sequences
		id_column		:	str, name of the column containing IDs to be
							used in file name

	'''
	if not os.path.isfile(protein_file):
		print("Could not find file %s" % protein_file)
		exit(1)
	else:
		f = pd.read_csv(protein_file)
		if protein_column not in f.columns:
			msg_fmt = "Successfully found file %s, but could not find a ``%s'' column in that file."
			msg = msg_fmt % (protein_file, protein_column)
			exit(1)
		u = 0
		for i in f.index:
			trace_df = sequenceIupred(f.ix[i,protein_column], kind=kind)
			if id_column in f.columns:
				trace_df.to_csv('%s_iupred_trace.csv' % f.ix[i, id_column], index=False)
			else:
				trace_df.to_csv('sequence_%d_iupred_trace.csv' % u, index=False)
			u += 1

def findIupredRegions(
		sequence, 
		kind='short', 
		identifier=None, 
		threshold=0.5, 
		plot=False, 
		filter_sigma=None
	):
	trace_df = sequenceIupred(sequence, kind=kind)
	trace = trace_df.score.tolist()
	if filter_sigma:
		trace = gaussian_filter(trace, filter_sigma)
	regions = []
	i = 0
	while i < len(trace):
		if trace[i] > threshold:
			start = i
			i += 1
			while i < len(trace) and trace[i] > threshold:
				i += 1
			regions.append((start, i))
		i += 1
	if plot:
		plt.plot(trace_df.index, trace, color='m', linestyle='-', linewidth=3, alpha=0.7)
		for start, stop in regions:
			plt.axvspan(start, stop, alpha=0.5, color='r')
		plt.xticks(plt.gca().get_xticks())
		plt.ylim((0, 1))
		plt.yticks([0.1*i for i in range(11)])
		plt.ylabel('IUPred score')
		plt.xlabel('Residue')
		plt.xlim((0, plt.gca().get_xlim()[1]))
		if identifier:
			wrapup('%s_iupred_regions.png' % identifier)
		else:
			wrapup('temp.png')
	return regions

def sequenceComposition(
		regions_file, 
		sequence_file, 
		protein_column='protein'
	):
	r_df = pd.read_csv(regions_file)
	s_df = pd.read_csv(sequence_file)
	for i in r_df.index:
		identifier = r_df.ix[i, 'protein_index']
		sequence = s_df.ix[identifier, protein_column]
		subsequence = sequence[r_df.ix[i, 'region_start'] : r_df.ix[i, 'region_end']]
		for aa in amino_acids:
			r_df.ix[i, '%s_freq' % aa] = float(sum([1 for j in subsequence if j==aa]))/len(subsequence)
		r_df.ix[i, 'complexity'] = complexity(subsequence)
		r_df.ix[i, 'protein_complexity'] = complexity(sequence)
	r_df.to_csv(regions_file, index=False)

def complexity(sequence):
	return - sum([(float(sequence.count(aa))/len(sequence)) * math.log(float(sequence.count(aa)) / len(sequence), 2) for aa in amino_acids if aa in sequence])
		

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Find disordered regions based on IUPred score in a CSV of sequences'
	)
	parser.add_argument(
		'infile', 
		type=str, 
		help='CSV containing full protein sequences'
	)
	parser.add_argument(
		'-t',
		'--threshold',
		type=float,
		help='iupred score threshold to call disordered regions. Default 0.5',
		default=0.5
	)
	parser.add_argument(
		'-s',
		'--sigma',
		type=float,
		help='Kernel width to use for smoothing, if desired. Default None.',
		default=None
	)
	parser.add_argument(
		'-p',
		'--protein_column',
		type=str,
		default='protein',
		help='Name of column containing the protein sequences'
	)
	parser.add_argument(
		'-i',
		'--id_column',
		type=str,
		default='name',
		help='Name of column containing sequence IDs to be used in filenames'
	)
	parser.add_argument(
		'-k',
		'--kind',
		type=str,
		default='short',
		help='Kind of IUPred score to calculate. Default is short. Possibilities short, long, glob'
	)
	parser.add_argument(
		'-l', 
		'--plot', 
		action='store_true', 
		default=False, 
		help='Plot individual protein traces. Default False'
	)
	parser.add_argument(
		'-c',
		'--seq_composition',
		action = 'store_true',
		help = 'Also calculate the sequence composition of each disordered region. Default False',
		default = False
	)
		
		

	args = parser.parse_args()
	f = pd.read_csv(args.infile)
	regions_filename = '%s_disordered_regions.csv' % args.infile.replace('.csv', '')
	regions_file = open(regions_filename, 'w')
	regions_file.write('protein_id,protein_index,protein_length,region_start,region_end,region_length\n')
	if args.protein_column not in f.columns:
		print('Could not find column %s' % args.protein_column)
		exit(1)
	for _idx in f.index:
		if args.id_column in f.columns:
			identifier = f.ix[_idx, args.id_column]
		else:
			identifier = str(_idx)
		sequence = f.ix[_idx, args.protein_column]
		regions = findIupredRegions(sequence, kind=args.kind, filter_sigma=args.sigma, threshold=args.threshold, identifier=identifier, plot=args.plot)
		for a, b in regions:
			regions_file.write('%s,%d,%d,%d,%d,%d\n' % (identifier, _idx, len(sequence), a, b, b-a))
		if _idx % 10 == 0:
			sys.stdout.write('Finished with %d/%d sequences\r' % (_idx, len(f)))
			sys.stdout.flush()
	regions_file.close()
	if args.seq_composition:
		print('Calculating sequence composition of IDRs...')
		sequenceComposition(regions_filename, args.infile, protein_column=args.protein_column)
