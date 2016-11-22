======================================================================================================
idrScan: a Python utility for calculating and comparing intrinsic disorder scores for individual exons.
======================================================================================================

Motivation.
~~~~~~~~~~~

Protein disorder prediction algorithms often include information from an entire protein sequence to calculate the order and disorder of individual residues. This class of algorithms includes the popular predictor IUPRED. However, such an approach presents a challenge for protein disorder prediction individual exons, which are isolated from the rest of the protein sequence. ``idrScan`` is a set of tools operating within the Pandas framework to assign disorder scores to individual exons within the context of the complete protein sequence. It is designed to work with the ``translateExon`` module, accepting that module's files as inputs.

All of the ``idrScan`` functions are designed to operate on CSV files containing exon information and to be run from the command line.

Component Programs.
~~~~~~~~~~~~~~~~~~~

``sequenceEntropy.py``: calculate the simple information (base 2) entropy of individual exons. Entropy is an effective way of detecting homopolymeric sequences. Short sequences (less than 50 amino acids) have entropy sharply skewed by sampling effects.

``iupredExon.py``: calculate the mean IUPRED score and its standard deviation for individual exons in a transcript file. IUPRED is a predictor that assigns each amino acid an energetic contribution to the stability of the protein; this contribution is a function of the frequencies of all other amino acids in the sequence. As such, IUPRED scores will be different when the algorithm is run on an individual exon or on an entire protein sequence. 
