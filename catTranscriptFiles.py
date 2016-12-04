'''
catTranscriptFiles.py -- for a directory of transcript files, find the full
protein for each transcript and write the full proteins to a new CSV
'''
import pandas as pd
import os
import argparse

def catTranscriptFiles(filenames, outfile):
	'''
	Find the full protein sequences for a set of transcript files and write them
	to a new CSV.

	INPUT
		filenames : list of string, names of transcript files or directories containing
			transcript files
		outfile : string, name of file to write full protein sequences to

	RETURNS
		Pandas DataFrame, containing the outfile information

	'''
	target_files = []
	for filename in filenames:
		if os.path.isdir(filename):
			subfilenames = os.listdir(filename)
			for subfilename in subfilenames:
				if '.csv' in subfilename:
					target_files.append('%s/%s' % (filename,subfilename))
		elif os.path.isfile(filename):
			target_files.append(filename)
		else:
			continue
	if len(target_files)==0:
		print "Failed to find CSVs in the specified files or directories"
		exit(1)
	if len(target_files)==1:
		print "Only one CSV specified"
		new_row = getFullProtein(target_files[0])
		result = pd.DataFrame(columns=new_row.index)
		result = result.append(new_row)
		result.to_csv(outfile, index=False)
		print "Successfully concatenated file %s" % target_files[0]
		return result
	else:
		new_row = getFullProtein(target_files[0])
		print new_row
		result = pd.DataFrame(columns=new_row.index)
		result = result.append(new_row)
		print "Successfully concatenated file %s" % target_files[0]
		for f in target_files[1:]:
			new_row = getFullProtein(f)
			print new_row #debugging
			result = result.append(new_row)
			print "Successfully concatenated file %s" % f
		result.to_csv(outfile,index=False)
		return result
	for f in target_files:
		try:
			dataframes.append(pd.read_csv(f))
		except pd.io.common.CParserError:
			print "File %s does not appear to be CSV" % f
			continue
	

def getFullProtein(transcript_file, column_drops=['constitutive', 'startPhase', 'endPhase', 'protein_start', 'protein_end']):
	'''
	Concatenate the exons from one transcript file into a single protein, and return
	a DataFrame containing the full protein and various information about the protein.

	INPUT
		transcript_file : string, name of the target file
		column_drops : list of string, columns to be selectively excluded from the
			final dataframe

	RETURNS
		Pandas DataFrame, containing full protein information (one row)

	'''
	try:
		f = pd.read_csv(transcript_file)
		f = f.sort_values(by='rank')
		f = f.set_index('rank', drop=False)
	except (TypeError, KeyError, ValueError, pd.io.common.CParserError) as e4:
		print "Could not read file %f" % f
		return False
	full_protein = ''
	for i in f.index:
		try:
			full_protein += f.ix[i,'protein']
		except TypeError:
			continue
	#do some violence to the dataframe
	f.ix[1,'protein']=full_protein
	if len(f.index)==1:
		return f
	else:
		columns = [j for j in f.columns if j!='protein']
		for column in columns:
			if not all([f.ix[i,column]==f.ix[1,column] for i in range(2, len(f.index))]):
				f = f.drop(column, 1)
			elif column in column_drops:
				f = f.drop(column, 1)
			else:
				continue
		return f.loc[1]

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='for a set of transcript files or directories containing transcript files, make a new file with the full protein sequences of each transcript')
	parser.add_argument('-i', '--infile', type=str, nargs='+', help='names of transcript files or directories containing transcript files to concatenate')
	parser.add_argument('-o', '--outfile', type=str, help='file to write to', required=True)

	args = parser.parse_args()
	try:
		catTranscriptFiles(args.infile, outfile=args.outfile)
	except (TypeError, KeyError, ValueError) as e3:
		print "TypeError encountered. Make sure that all inputs are CSVs."
		exit(1)
