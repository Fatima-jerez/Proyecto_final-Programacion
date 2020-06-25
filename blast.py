#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from subprocess import call, Popen, PIPE
from re import sub

def is_fasta(input):
        """
        comprueba que el archivo que hemos introducido es fasta y en caso contrario aborta el script
	"""        
	result = False
	for line in input:
		if line.startswith(">"):
			result = True
                        break
			if result == False:
				sys.exit("\nError: Query file with incorrect format. Write -h to see the help option")

#Creación del archivo multifasta a partir de los Genbank almacenados en la carpeta que introduce el usuario 
def convert(input_folder):
	"""
	Con esta función vamos a convertir los archivos en formato GenBank
	almacenados en la carpeta que introduzca el usuario como segundo argumento
	a un único archivo multifasta
	"""
#Comprobamos que el archivo que queremos crear no existe previamente, en cuyo caso lo eliminamos
	if os.path.isfile("fasta_seq.fa") == True:
		os.remove("fasta_seq.fa")
	else:
		pass
#Utilizamos la opción -Append- para realizar el open, que nos permite añadir información 
#si ya se existe la carpeta o en caso de que no exista, la crea
	folder=open("fasta_seq.fa", "a")
	for file in os.listdir(input_folder):
		if os.path.isfile(input_folder+file) == True:
#Comprobamos que sean archivos los documentos que están almacenados
#en folder y que hemos recogido en una lista mediante listdir
			with open(input_folder+file, "r") as input_handle:
				for record in SeqIO.parse(input_handle, "genbank"):
					for feature in record.features:
						if feature.type == 'CDS':
							try:
								protein = feature.qualifiers['translation'][0]
							except:
								protein = "No protein"
							if protein != "No protein":
								folder.write(">"+feature.qualifiers['locus_tag'][0]+"\n"+protein+"\n")
	folder.close()

def DB():
	"""
	Esta función nos permite crear la base de datos con la que se realizará el 
	Blastp a partir del archivo multifasta
	"""
#Creamos una carpeta en la que almacenar los archivos de la base de datos
	directory = 'Database/'
	dir = os.path.dirname(directory)
	if not os.path.exists(dir):
		try:
			os.mkdir('Database/')	
		except OSError:
			print("Creation of the directory failed")

	multifasta="fasta_seq.fa"
	db='Database/DataBase'
	DBcreate_args=['makeblastdb', '-dbtype', 'prot', '-parse_seqids',  '-in', multifasta, '-out', db]
	try:
		call(DBcreate_args)
	except Exception as error:
		print ('BLAST database not created')
		print (error)
		sys.exit()
	return db
	

#almacenamos en un diccionario las secuencias del query que hemos obtenido
def dicc(input):
	"""
	Esta función nos permite guardar en un diccionario los nombres de las query,
	siendo estos las "key" y las secuencias correspondientes el "value"
	"""
	seq_query=dict()
        qseq=str()
	new_file=open(input, "r")
	count=0
        for line in new_file:
		if line.startswith(">"):
			if count != 0:
				seq_query[q_id]=qseq
				qseq=str()
			#Utilizamos strip() sin parámetros que nos permite eliminar de los string los espacios en blanco
			q_id=line.strip()
			q_id=q_id[1:]
			count=1
		elif line.startswith("\n") == False:
			qseq=qseq+line.strip()
		else:
			continue


	seq_query[q_id]=qseq
	new_file.close()
	return seq_query

#Una vez que hemos creado la base de datos a partir del multifasta
#la utilizamos para realizar el blast con las secuencias query
#introducidas por el usuario como primer argumento y con la identidad
#y coverage como tercer y cuarto argumento, respectivamente. 

def blast(input):
	database="Database/DataBase"
	dir='Results/'
        if os.path.isdir(dir):
                pass
        else:
                os.mkdir(dir)

	blast_args = NcbiblastpCommandline(query=input, out="Results/result_blastp", db=database, outfmt="\'6 qseqid sseqid pident qcovs qlen slen length bitscore evalue\'")
#Este formato nos aporta la siguiente información
#qseqid sseqid pident qcovs qlen slen length bitscore evalue
	try:
		print(blast_args)
		stdout, stderr = blast_args()
	except:
		print("Error: Blast couldn't be performed")
		sys.exit()
	result2=open("Results/asked_results_blastp", "w")
	f=open("Results/result_blastp", "r")
	seq=open("Results/sseq", "w")
	for line in f:
		fields=line.split("\t")
		#Fijo el valor del e-value que corresponde a la columna 10 del archivo
		if float(fields[8]) <= 0.00001:
			#Selecciono el qseqid sseqid pident qcovs evalue
			result2.write(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[8])
			results=open("Results/asked_results_blastp", "r")
			for fig in results:
				colum=fig.split("\t")
				seq.write(colum[1] + "\n")
#Guardamos en el arvhio "sseq" los nombres de los subject que nos han salido del blast
#para buscar a continuación las secuencias en la base de datos que hemos creado
#y poder añadirla al archivo "asked_results_blastp"
	f.close()
	seq.close()
	result2.close()
	results.close()
	blasp = ['blastdbcmd', '-entry_batch', 'Results/sseq', '-db', database, '-outfmt', '%s', '-out', 'Results/subjects_seq']
	
	try:
		call(blasp)
	except:
		print("Error: Blast couldn't be performed")
		sys.exit()
	f.close()
	seq.close()
	
	alin=open("Results/subjects_seq", "r")
	seqls=list()
	for item in alin:
		seqls.append(item)
	alin.close()

	output=open("Results/asked_results_blastp", "r")
	blastp=open("Results/finalresult_blastp", "w")
	cont = 0

	#Unimos el subjectID con la secuencia correspondiente en el archivo "finalresult_blastp"
	for str in output:
		str=str.strip()
		blastp.write(str + "\t" + seqls[cont])
		cont = cont + 1

	output.close()
	blastp.close()

	subjects_blast=open("Results/subjects_blast.fa", "w")
	inputfile=open("Results/finalresult_blastp", "r")
	for column in inputfile:
		columns=column.split("\t")
		subjects_blast.write(">" + columns[1] + "\n" + columns[5])
	subjects_blast.close()
	inputfile.close()

def fasta_muscle(perc, coverage, query_dic):
	"""
	Esta función nos permite crear los archivos fasta para realizar los 
	alineamientos y el árbol con muscle. La forma que tendrá será los
	key que habíamos almacenado en el diccionario y su secuencia
	correspondiente y los subject, a partir del archivo de Blast
	"""
	dir='Results_blastp/'
	if os.path.isdir(dir):
		pass
	else:
		os.mkdir(dir)
		#Creamos la carpeta Results_blastp  dentro de la cual habrá un archivo para con el nombre de cada query
		#que contendrá la secuencia query y sus subjects
	for key in query_dic:
		f_query="Results_blastp/"+key+".fa"
		f_sequence=open(f_query, "w")
		f_sequence.write(">" + key + "\n" + query_dic[key] + "\n")
		entrada=open('Results/finalresult_blastp', "r")
		for line in entrada:
			columns=line.split("\t")
			if columns[0] == key:
		#Nuestra primera columna del blast corresponde a los nombres de las queries del diccionario
		#nuestra segunda columna del blast es el coverage y la tercera el porcentaje de identidad
				if float(columns[3])>=float(coverage) and float(columns[2])>=float(perc):
			#Escribimos en el mismo archivo donde habías escrito anteriormennte
			#las queries, las secuencias resultantes del blast
					f_sequence.write(">" + columns[1] + "\n" + columns[5])

		f_sequence.close()
		entrada.close()

#input=sys.argv[1]
#input_folder=sys.argv[2]
#coverage=sys.argv[3]
#perc=sys.argv[4]

#is_fasta(input)
#seq_query=dicc(input)
#convert(input_folder)
#DB()
#dicc(input)
#blast(input)
#seq_query=dicc(input)
#fasta_muscle(seq_query, perc, coverage)

