# BLAST--MUSCLE_ALIGNMENT--DOMAIN_SEARCH--PHYLOGENETIC TREE

The script main.py performs a BLASTp research, making use of Biopython, between a fasta file (query) containing one 
or several sequences, and several genbank format files, which  must be all together kept in a directory. 
This genbanks will be used by the program to creat a Database.  The user can choose the coverage and identity values to filter the results. 
The e-value has been fixed in <= 0.00001. The output of BLASTp is kept in a folder called "Results" where you will find diferent result_files: 
	- results_blastp : this file contains all the information obtain by performing Bio.Blast in 6 format
	- asked_results_blastp : this file extracts the columns (qseid sseqid pident qcovs evalue) 
	- finalresult_blastp : this file adds the subject sequences to the 'asked_results_blastp' file 
With the resulting sequences, this script creates for each query a multifasta file, containing the query
sequence and subjects found for each one. This files, in "Results_blastp/", will be used then to do a Muscle alignment that
will allow the user to get a Neighbour-joining phylogenetic tree for each query to know the distance to the other proteins.
This alignment file called "alignments.fa" is in "Muscle_alignment/".  The program also provides for each query, a file containing the protein domains found by parsing the database "Prosite.dat" that must be uploaded in advance by the user. 
This results can be found in "Data_prosite" in different files with .txt extension. Previously, the script also creates a file with all the proteins and some related information extract from Prosite.


#Instalation

Not required

#USAGE

	Usage: main.py	[-h] 

- [h]: help* : this is opens a function with a help message about the program. 

	Usage: main.py [fasta file] [directory/] [identity cut-off] [coverage cut-off]

- [fasta file] : A FASTA format file containing one or more protein sequences to perform BLASTp. 
- [directory/] : A directory with diferent GenBank format files to create a Database for BLASTp.
- [identity cut-off] : Percent identity must be a number between 0  and 100
- [coverage cut-off] : Query Coverage per suject must be a number between 0 and 100
