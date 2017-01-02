'''
iupredProtein.py -- calculate IUPRED scores for entire proteins (non-exons)
'''
import pandas as pd
import argparse
import os
import math

def sequenceIupred(sequence, kind='short'):
    '''
    Given a sequence, calculate the IUPRED score of each residue.

    INPUT
        sequence : string, amino acids in the sequence
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
    os.system('rm %s' % temp_fasta)
    os.system('rm %s' % temp_iupred)
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

def iupredProtein(protein_file, kind='short', protein_column='protein'):
	'''
	Given a file with the entire protein sequences for a list of transcripts,
	calculate the mean IUPred score for each protein.

	INPUT
		protein_file : string, name of file containing protein sequences under
			the column ``protein''

	RETURNS
		pandas DataFrame, containing the extra iupred_protein column

	'''
	if not os.path.isfile(protein_file):
		print "Could not find file %s" % protein_file
		exit(1)
	else:
		f = pd.read_csv(protein_file)
		if protein_column not in f.columns:
			print "Successfully found file %s, but could not find a ``%s'' column in that file." % (protein_file, protein_column)
			exit(1)
		for i in f.index:
			iupred_protein_list = sequenceIupred(f.ix[i,protein_column], kind=kind)
			scores = [j[2] for j in iupred_protein_list]
			mean_iupred = float(sum(scores))/len(scores)
			stdev_iupred = math.sqrt(float(sum([(mean_iupred-k)**2 for k in scores]))/len(scores))
			f.ix[i,'iupred_mean']=mean_iupred
			f.ix[i,'iupred_stdev']=stdev_iupred
			if 'description' in f.columns:
				print "Successfully calculated IUPred scores for %s" % f.ix[i,'description']
		f.to_csv(protein_file, index=False)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='calculate IUPred scores for whole protein sequences')
	parser.add_argument('infile', type=str, help='CSV containing the full protein sequences')
	parser.add_argument('-p', '--protein_column', type=str, help="name of the column containing the protein sequences in the infile. Default is ``protein''", default='protein')
	parser.add_argument('-k', '--kind', type=str, help="kind of IUPred score to calculate. default is ``short''", default='short')
			
	args = parser.parse_args()
	try:
		iupredProtein(args.infile, kind=args.kind, protein_column=args.protein_column)
	except (TypeError, pd.io.common.CParserError) as e_all:
		print "Could not read the input files."
		exit(1)
