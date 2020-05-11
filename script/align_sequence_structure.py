#!/usr/bin/env python

#script to build sequence, structural and combined trees
#input fasta file with sequences and directory containing structures


#necessary packages
import subprocess
import sys
import re
import math
import numpy
import glob
from Bio import AlignIO

#import structural functions
from id_matrix import struc_matrix



#function to combine sequence and structure distance 
#matrices into a supermatrix
def supermatrix(seq_dist, struc_dist):
	super = list();

	print("process files...")

	sq_dist = seq_dist

	st_dist = struc_dist

	print("build matrix...")
	for i in range(len(st_dist)):

		#use only sequence if distance is less than 0.4
		if(sq_dist[i] < 0.4):
			mean = sq_dist
		else:
			mean = [numpy.mean(j) for j in zip(st_dist[i],sq_dist[i])]

		super.append(mean)
	return(super)

#function to process the structural alignment file
#produces distance matrix
def process_struc(struc_file):
	f = open(struc_file);

	ident = list();

	for line in f:

		if("#" in line):
			next;
		elif(re.search(r'^\d+\n', line)):
			next
		else:
			vals = re.findall(r'\d+\.\d+', line)
			ident.append(vals);
	
	ident = [list(map(float,i)) for i in ident]

	f.close();
	return(ident)

#function to normalise distance matrix
#uses min max normalisation
def normalise(dist_matrix):
	ndist_matrix = list()

	for dist_list in dist_matrix:

		nlist = list()
		max_val =  max(dist_list)
		min_val = min(dist_list)

		for dist in dist_list:
			ndist = (dist - min_val)/(max_val - min_val)
			nlist.append(ndist)

		ndist_matrix.append(nlist)

	return(ndist_matrix)

#function to process sequence alignemnt file in clustalW format
#returns identity matrix
def process_seq(seq_file):
	f = open(seq_file);

	ident = list();

	for line in f:

		if("#" in line):
			next;
		else:
			vals = re.findall(r'\d+\.\d+', line)
			ident.append(vals);
	
	ident.pop(0)
	ident = [list(map(float,i)) for i in ident]
	
	f.close()

	return(ident)

#funtion to process sequence distance file produced by tree puzzle
#returns distance matrix
def process_treepuz_dist(dist_file):
	f = open(dist_file);

	dist = list();

	columns = 0

	col = None

	for line in f:
		
		if(re.search(r'^([A-Z]|[a-z]|[0-9]){4}', line)):
			columns = columns + 1

			if (col != None):
				dist.append(col)

			col = list()

			names = re.findall(r'^([A-Z]|[a-z]|[0-9]){4}',line)
			vals = re.findall(r'\d+\.\d+', line)
			
			for  val in vals:
				col.append(float(val))

		elif(columns > 0):
			vals = re.findall(r'\d+\.\d+', line)
			for  val in vals:
				col.append(float(val))

	dist.append(col)

	f.close();
	return(normalise(dist))

#function to align strutures and produce distance matrix
#	p = True prints filenames being processed
#returns distance matrix
def struc_matrix(files, p = False):

	size = len(files)

	dist_matrix = [[0]*size for i in range(size)]

	count_1 = 0

	for file_1 in files:
		count_2 = 0
		if (p == True):
			print("file 1:")
			print(file_1)

		name_1 = clean_file_name(file_1)

		for file_2 in files:
			if (p == True):
				print ("file 2")
				print (file_2)
			name_2 = clean_file_name(file_2)

			if (file_1 != file_2):
				command = "./TMalign " + file_1 +" " + file_2 + " | tee alignment >> alignments/" + name_1 + "_" + name_2
				print(command)

				subprocess.call(command, shell=True)

				file = "alignment"

				scores = (get_score(file))

				dist_matrix[count_1][count_2] = round(1 - float(scores[0]),4)

			else:
				dist_matrix[count_1][count_2] = 0

			count_2 = count_2 + 1
		count_1 = count_1  + 1

	return (dist_matrix)


#function converts identity matrix to distance matrix
#for conversion of sequence identity matrix
def id_to_dist(ident):
	dist = list();
	for row in ident:
		nrow = [(1-(i/100)) for i in row]
		dist.append(nrow)

	return(dist)

#funtion to retrieve protein ids from alignment files
#l input indicates long name or not
#	l=True when labels have not been changed after 
#	downloading from pdb
#returns list of labels
def get_names(file, l=False):
	f = open(file);

	names = list();

	for line in f:

		if("#" in line):
			next;
		elif(re.search(r'^\d+\n', line)):
			next;
		elif(re.search(r'^\n', line)):
			next;
		elif(l==True):
			name = line.split("_");
			name = name[0];
			name = name.split(": ");
			name = name[1];
			name = name.lower();
			names.append(name);
		else:
			name = line.split(": ")
			name = name[1];
			name = name.split(" ")
			name = name[0]
			name = name.lower();
			names.append(name);

	f.close();
	return(names)

#function to find labels for proteins in tree puzzle output file
#returns list of labels
def get_treepuz_names(file):
	f = open(file)

	names = list()

	for line in f:
		if(re.search(r'^([A-Z]|[a-z]|[0-9]){4}', line)):
			name = line.split(" ")
			name = name[0]
			name = name.lower()

			names.append(name)

	f.close()
	return(names)

#get TM score from structural alignment
def get_score(file):
	alignment = open(file)

	aligned = ()
	scores = []


	for line in alignment:
		if "TM-score=" in line:
			split = line.split(" (")
			split = split[0].split("= ")
			scores.append(split[1])

	return scores

#function to write distance matrix as mega file
#	super is the distance matrix to be written
#	names is the protein ids
#	filename is file to be written
def write_file(super,names,filename):
	f = open(filename, "w+")

	f.write("#mega\n!TITLE distance matrix; \n!Format DataType=distance;\n"
		+ "!Description \n\tproteins;\n\n");

	for name in range(len(names)):
		f.write("#" + names[name] + "\n")

	f.write("\n")

	for r_num in range(len(super)):
		row = super[r_num];

		for i in range(r_num):
			f.write(str(row[i]) + "\t");

		f.write("\n")


#function to strip pdb file names and produce just pdb ids
#returns list of pdb ids from filenames
def clean_file_names(pdb_files):
	names = list()

	for file in pdb_files:
		file = file.split("/")[1]
		name = file.split(".pdb")[0]

		names.append(name)

	return names

#function to add directory and pdb extension to pdb ids
#returns list of files
def pdb_ext(pdb_files, direc):
	files = list()

	for file in pdb_files:
		nfile = direc + "/" + file + ".pdb"

		files.append(nfile)

	return files

#function to convert fasta alignment to phylip alignemnt
#writes new file "file_out" based on "file_in"
def fasta_to_phylip(file_in, file_out):

	input_handle = open(file_in, "rU")
	output_handle = open(file_out, "w")

	alignments = AlignIO.parse(input_handle, "fasta")
	AlignIO.write(alignments, output_handle, "phylip")

	output_handle.close()
	input_handle.close()

	return("done")



####################### GET FILES TO AILIGN #######################

fasta_file = sys.argv[1]

pdb_files = (glob.glob(sys.argv[2] + "/*.pdb"))


############################ SEQUENCE ALIGN ############################

email = "beni.cawood@gmail.com"

#run clusalW alignment
command = "python clustalo.py --email " + email + " --outfmt fa " + fasta_file
print("Running Clustal..."
subprocess.call(command, shell=True)

#get ouput distance file
files = glob.glob("clustalo*.pim")
seq_file = files[0]

#when using only clustalW and not trimming uncomment this section
#names = get_names(seq_file)
#seq = process_seq(seq_file)
#seq_dist = id_to_dist(seq)


#trim with gblocks
file = glob.glob("clustalo*.fasta")
print ("Running Gblocks...")
command = "./Gblocks " + file[0] + " -t=p -b4=5 -b5=a -b2=6"
subprocess.call(command, shell=True)


#run tree puzzle
gfile = glob.glob("clustalo*.fasta-gb")
fasta_to_phylip(gfile[0], "tree_in.phy")
print("Running Tree puzzle...")
command = "puzzle tree_in.phy"
subprocess.call(command, shell=True)


#get ouputted distance file from tree puzzle
dfile = glob.glob("tree_in.phy.dist")
seq_dist = process_treepuz_dist(dfile[0])
names = get_treepuz_names("tree_in.phy.dist")



######################### STRUCTURAL ALIGNMENT ######################### 

#use sequence ids to correclty order structural alignment
pdb_files = pdb_ext(names, sys.argv[2])
#perform structural alignemnts 
struc_dist = struc_matrix(pdb_files,p=False)


########################## WRITE OUTPUT FILES ########################## 

#write combined distance matrix
write_file(supermatrix(seq_dist,struc_dist),names, "supermatrix.meg")

#write structural matrix
write_file(struc_dist, names, "distmatrix.meg")

#write sequence matrix
write_file(seq_dist, names,"seqmatrix.meg")


############################# BUILD TREES #############################

print("build tree...")
#build combined tree
subprocess.call("megacc -a infer_ME_distances.mao -d supermatrix.meg", shell=True)
#build structural tree
subprocess.call("megacc -a infer_ME_distances.mao -d distmatrix.meg", shell=True)
#build sequence tree
subprocess.call("megacc -a infer_ME_distances.mao -d seqmatrix.meg", shell=True)
print("done!")


#print trees in newark fomat in terminal
files = glob.glob("M10CC_Out/*.nwk")
tree_file = files[0];

print("tree in newark format: \n")
with open(tree_file, 'r') as tree:
    print(tree.read())




