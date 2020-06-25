#!/usr/bin/python
# -*- coding:utf-8 -*-

import sys
import os
from subprocess import call

def muscle_align():
        """
        Esta función permite a partir de los archivos fasta obtenidos
        del blast para cada una de los subject, poder realizar un alineamiento
        con muscle.
        """
        DIR="Results_blastp/"
        dir='Muscle/'
        if not os.path.exists(dir):
                try:
                        os.mkdir(dir)
                except:
                        pass 

		for file in os.listdir(DIR):
			input_file=DIR + file
			output_file=dir+file[:-3]+"_m.align"
			muscle_cmd = ['muscle', '-in', input_file, '-out', output_file]
			
			try:
				call(muscle_cmd)
			except:
				print("Alignment not performed")
				sys.exit()
			output='Muscle/'+file[:-3]+'_tree.nw'
			tree = ['muscle', '-maketree', '-in', output_file, '-out', output, '-cluster', '-neighborjoining'] 
			try:
				call(tree)
			except:
				print("No phylogenitc tree")
				sys.exit()
