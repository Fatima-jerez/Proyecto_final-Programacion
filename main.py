#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os

import blast as bs
import  muscle_align as ms
import dominios as dm

def ayuda():
        #esta función muestra un mensaje de ayuda más simple
        print("Usage: main.py [ fasta file ] [ directory/ ] [ identity cut-off ] [ coverage cut-off ] \n")
        print(" or main.py [-h to get futher information")

def ayuda2():
	#versión de ayuda más desarrollada
	help = (
	"""
	WELCOME TO HELP OPTION
	-------------------------------------------------------------------------------------------------------------
        Usage: main.py [-h ]
	
	- [ h ] : print in terminal the help option of the program
        -------------------------------------------------------------------------------------------------------------
	Usage: main.py [fasta file] [directory/] [identity cut-off] [coverage cut-off]

	- [fasta file] : A FASTA format file containing one or more protein sequences to perform BLASTp. 
	- [directory/] : A directory with diferent GenBank format files to create a Database for BLASTp.
	- [identity cut-off] : Percent identity must be a number between 0  and 100
	- [coverage cut-off] : Query Coverage per suject must be a number between 0 and 100
	-------------------------------------------------------------------------------------------------------------
	"""
	)
	print(help)
	sys.exit()

def input_control():
	"""
	Esta función nos permite comprobar que el número de argumentos
	que se ha introducido sea el correcto, que se encuentren en el formato
	que se requiere para ejecutar los módulos y distinguir cuando el usuario
	invoca la opción de ayuda
	"""
	arguments = sys.argv
	for argument in arguments:
		if argument == "-h":
			ayuda2()
			sys.exit()
	input=sys.argv[1]
	input_folder=sys.argv[2]
	if len(sys.argv) == 3:
		coverage=50
		perc=30
	elif len(arguments) == 5:
		perc=sys.argv[3]
		coverage=sys.argv[4]
		try:
			coverage=float(coverage)
			perc=float(perc)
		except:
			print('"Coverage" and "Identity" must be numbers between 0 and 100')
			sys.exit()
		if perc>100 or perc<0:
			print("ERROR: main.py -h to see the help option")
			sys.exit()
		if coverage>100 or coverage<0:
			print('ERROR: main.py -h to see the help option')
			sys.exit()
	else:
		print('ERROR: Incorrect number of arguments')
		ayuda()
		sys.exit()
	return(input, input_folder, perc, coverage)


def main():
	try:
		input, input_folder, perc, coverage=input_control()
		print('')
		print("This script performs a BLASp research with the query introduced by the user and a Database with genome sequences"+"\n"+" With the resulting proteins the program will extract in different files, the domains of this proteins found in 'Prosite.dat'"+"\n"+"and will perform an alignment and phylogenetic tree to see the relation between them.")
		print("----------------------------------------------------------------------------------------------------------------------------------------")
		print("Working on your request")
		print("")
		print("Almost got it")
		print("")
		print('---BLASp research in database---')
		print('---Phylogenetic tree created with protein results---')
		print('---Protein domains in "Data_prosite/"---')
		
		#IMPORT FROM BLAST-SCRIPT
		bs.is_fasta(input)
		query_dic=bs.dicc(input)
		convertir=bs.convert(input_folder)
		database=bs.DB()
		bs.dicc(input)
		bs.blast(input)
		bs.fasta_muscle(perc, coverage, query_dic)
		#IMPORT FROM MUSCLE
		ms.muscle_align()
		#IMPORT FROM DOMINIOS	
		dm.doc()
		dictionary=dm.dicc_domains()
		dm.parsear(dictionary)
		print('')
		print('--------------------------------------------------------------------------------')
		print(" YOUR RESULTS ARE READY ")
		print("--------------------------------------------------------------------------------")
		print(" --> --> Please read de README.md file for better understanding the program work")
	except:
		sys.exit()

main()


