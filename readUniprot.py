'''
readUniprot.py -- extract information from FASTA files and include information
from Uniprot TSV export files
'''
import pandas as pd
import argparse

def readFasta(fasta_file, name_column='name', sequence_column='sequence'):
	'''
	Read a FASTA and return as a Pandas DataFrame.
	
	INPUT
		fasta_file : string, the name of the FASTA file
		name_column : string, name of the FASTA header column in the
			output file
		sequence_column : string, name of the sequence column in the
			output file

	RETURNS
		Pandas DataFrame with columns [name, sequence]

	'''
	f = open(fasta_file,'r')
	flines = f.read().replace(',','&')
	f.close()
	n_seq = len([i for i in flines if i=='>'])
	flines = [fline for fline in flines.split('\n') if len(fline)>0]
	df = pd.DataFrame(columns=[name_column,sequence_column], index=range(n_seq))
	k = 0
	for i in range(len(flines)):
		if flines[i][0]=='>':
			sequence = ''
			j = i+1
			while flines[j][0]!='>' and j!=len(flines)-1:
				sequence += flines[j]
				j+=1
			df.ix[k,name_column] = flines[i].replace('>','')
			df.ix[k,sequence_column] = sequence
			k+=1
	return df

def readUniprot(fasta_file, metadata_file, delim='\t', write_to_file=False, outfilename='readUniprotDefault.csv'):
	'''
	Read a FASTA and corresponding metadata downloaded from Uniprot and
	associate the information in the metadata file with the sequences.
	
	INPUT
		fasta_file : string, FASTA containing the sequences
		metadata_file : string, file containing metadata. The ``Entry'' column
			of this TSV and the header of the FASTA file is used to associate
			sequences with metadata.
		delim : char, the delimiter used in the metadata file
		write_to_file : bool, whether to write the results to an output file
	
	RETURNS
		Pandas DataFrame with columns [fasta_header, sequence, [metadata]]

	'''
	fasta_f = readFasta(fasta_file, name_column='fasta_header')
	meta_f = pd.read_csv(metadata_file, sep=delim)
	for i in meta_f.index:
		if len(fasta_f[fasta_f['fasta_header'].str.contains(meta_f.ix[i,'Entry'])])!=0:
			meta_f.ix[i,'fasta_header'] = fasta_f.ix[i,'fasta_header']
			meta_f.ix[i,'sequence'] = fasta_f.ix[i,'sequence']
	return meta_f
			

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='read FASTA files and UNIPROT files and write to CSV')
	parser.add_argument('fasta_file', type=str, help='name of FASTA file')
	parser.add_argument('-u', '--uniprot', type=str, help='name of the metadata file')
	parser.add_argument('-o', '--outfile', type=str, help='name of the file to write the results to. Default <file_prefix>.csv')
	parser.add_argument('-d', '--delim', type=str, help='kind of delimiter used in the metadata file. Default tab.', default='\t')

	args = parser.parse_args()
	if not args.uniprot:
		print "Reading FASTA..."
		try:
			f = readFasta(args.fasta_file)
		except IOError:
			print 'Could not find the FASTA file %s' % args.fasta_file
	else:
		print "Synthesizing FASTA and metadata files..."
		try:
			f = readUniprot(args.fasta_file, args.uniprot, delim=args.delim, outfilename=args.outfile)
		except IOError:
			print 'Could not find the FASTA file %s or the UNIPROT file %s' % (args.fasta_file, args.uniprot)
			exit(1)
	print "Writing data..."
	if args.outfile:
		f.to_csv(args.outfile,index=False)
	else:
		outfile = '%s.csv' % args.fasta_file.split('.')[0]
		f.to_csv(outfile,index=False)

	
