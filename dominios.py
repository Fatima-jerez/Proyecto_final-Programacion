#!/usr/bin/python
# -*- coding:utf-8 -*-

import sys
import os
import re
from Bio.ExPASy import Prosite
from Bio import SeqIO


def doc():
	"""
	Esta función nos permite parsear el archivo Prosite.dat
	para extraer la el nombre, accesion, patrón y descripción
	de las diferentes proteínas
	"""
	input_file='prosite.dat'
	path='Data_prosite/'
	if os.path.isdir(path) == True:
		pass
	else:
		os.mkdir(path)
	information=open(path+"Data", "w")
	handle = open(input_file, "r")
	records = Prosite.parse(handle)
	for record in records:
		information.write(record.name + "\t")
		information.write(record.accession + "\t")
		information.write(record.pattern + "\t")
		information.write(record.description + "\n")
	information.close()
	handle.close()


def dicc_domains():
	"""
	Esta función permite crear un diccionario en el que almacenar
	como  key los accesion del archivo "Data_prosite/Data" y como valor
	los dominios correspondientes a cada key. Y a partir de este diccionario, 
	se buscan en las secuencias subject resultantes del Blastp los dominios
	que presentan cada una de las proteínas, que se almacenaran para cada
	una de ellas en un archivo .txt
	"""
	diccionario=dict()
	archivo=open('Data_prosite/Data', "r")
	for column in archivo:
		item=column.split("\t")
		name=item[0]
		accesion=item[1]
		pattern=item[2].strip()
		pattern=pattern.replace("(", "{")
		pattern=pattern.replace(")", "}")
		pattern=pattern.replace("x", ".")
		pattern=pattern.replace("-", "")
		pattern=pattern.strip()
		description=item[3]
		#Como algunos accesion no tienen patrón, no los incluimos en el diccionario
		if pattern == "":
			pass
		else:
			diccionario[accesion]=[pattern, name, description]
	archivo.close()
	return diccionario

def parsear(dictionary):
	#Parseamos ahora los archivo fasta resultantes de Blastp con las secuencias de proteínas de las que queremos sus dominios
	dir="Results_blastp/"
	for inputfile in os.listdir(dir):
		input_file=open(dir + inputfile, "r")
		output="Data_prosite/"+inputfile[:-3]+"_domains.txt"
		out_file=open(output, "w")
		for lines in input_file:
			if lines.startswith(">"):
				subID=lines
				subID=subID[1:]
				out_file.write(subID + "\n")
			elif lines.startswith("\n") == False:
				seq = lines.strip()
				count=0
				for key in dictionary:
					domains=re.compile(dictionary[key][0])
					findings=domains.search(seq)
					if findings == None:
						pass
					else:
						count += 1
						result=findings.group()
						out_file.write("Nº DOMAIN" + " " + str(count) + " " + "DOMAIN \n" + "DOMAIN_NAME " + dictionary[key][1] + "\n" + "ACCESSION " + key + "\n" + "PATTERN " + str(result) + "\n" + "DESCRIPTION " + dictionary[key][2] + "\n")
			else:
				pass

	input_file.close()
	out_file.close()
