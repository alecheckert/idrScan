'''
iupredProtein.py -- calculate IUPRED scores for entire proteins (non-exons)
'''
import pandas as pd
import argparse
import os
import math

from StringIO import StringIO


def sequenceIupred(sequence, kind='short'):
    '''
    Given a sequence, calculate the IUPRED score of each residue.
    Assumes that the *iupred* command is in the shell path.

    INPUT
        sequence : string, amino acids in the sequence
        kind : string, one of the three IUPRED options

    RETURNS
        A DataFrame with columns "pos" "AA" and "score"
        where pos is the index into the amino acid sequence,
        AA is the amino acid single letter code and score is the
        IUPred disorder score you requested.
    '''
    temp_fasta = '_temp.fasta'
    temp_iupred = '_temp.iupred'

    with open(temp_fasta, 'w') as outf:
        outf.write('>temp\n%s' % sequence)

    os.system("iupred %s %s > %s" % (temp_fasta, kind, temp_iupred))

    scores = readIupred(temp_iupred)
    os.system('rm %s' % temp_fasta)
    os.system('rm %s' % temp_iupred)
    return scores


def readIupred(iupred_file):
    '''
    Reads output of the type given by IUPRED. *iupred_file* is a
    file containing output exactly as given by the *iupred* command.

    RETURNS
        A DataFrame with columns "pos" "AA" and "score"
        where pos is the index into the amino acid sequence,
        AA is the amino acid single letter code and score is the
        IUPred disorder score you requested.
    '''
    with open(iupred_file, 'r') as g:
        glines = g.read().split('\n')
        glines = [i for i in glines if len(i) > 0]
        glines = [i for i in glines if i[0] != '#']

    glines = ['pos AA score'] + glines
    fstring = '\n'.join(glines)
    fakefile = StringIO(fstring)
    return pd.read_csv(fakefile, delim_whitespace=True, index_col=0)


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
            msg_fmt = "Successfully found file %s, but could not find a ``%s'' column in that file."
            msg = msg_fmt % (protein_file, protein_column)
            exit(1)

        for i in f.index:
            iupred_protein_df = sequenceIupred(
                f.ix[i, protein_column], iupred_dir, kind=kind)
            scores = iupred_protein_df.score
            f.ix[i, 'iupred_mean'] = scores.mean()
            f.ix[i, 'iupred_stdev'] = scores.std()

            if 'description' in f.columns:
                desc = f.ix[i, 'description']
                print "Successfully calculated IUPred scores for %s" % desc
        f.to_csv(protein_file, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='calculate IUPred scores for whole protein sequences')
    parser.add_argument('infile', 
                        type=str,
                        help='CSV containing the full protein sequences')
    parser.add_argument('-p',
                        '--protein_column', type=str, default='protein',
                        help="Name of the column containing the protein sequences in the infile.")
    parser.add_argument('-k', '--kind', type=str, default='short',
                        help="kind of IUPred score to calculate. default is ``short''")

    args = parser.parse_args()
    try:
        iupredProtein(args.infile, kind=args.kind,
                      protein_column=args.protein_column)
    except:
        print "Could not read the input files."
        exit(1)
