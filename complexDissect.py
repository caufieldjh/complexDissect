#!/usr/bin/python
#complexDissect.py
'''
A small utility for comparing sets of protein complexes.

Secondary goal: predict conservation of complexes across species using
protein orthology.

INPUT: 
Two files. Assumed to have names like "complexes*.txt" 
but will take input for names
Each needs to have the following format.
ProteinID	ComplexMembership
ID1	A
ID2	A
ID3	B
ID4	B
ID4	C

One of these files is for the "experimental" set.
The second file is the model set.
Identifiers should be Uniprot IDs but don't have to be as long as both
files use the same type of unique identifier (i.e., not gene names).
IDs will be converted to Uniprot IDs.
Proteins may belong to more than one complex, as shown above, so
interactors may not be unique but lines will be.

"Short" format files can also be parsed. These have the following format:
ComplexID	Proteins
A	ID1	ID2
B	ID3	ID4
C	ID4

Each line contains a complex name/identifier, followed by each of the
components of the complex, separated by tabs.
So if the data looks like this (e.g., Hu et al 2009, table S7):

CplxID	Complex_members
1	b1234	b1235	b1236

then it can be converted into "long" format.

OUTPUT:
A file indicating, for each "experimental" complex, conservation
in the model, both within a complex and across the full set.
For example, proteins A and B may both be present in the same complex
in the model set but protein C, though in that complex in the
experimental set, may be in a completely different complex.

Uses ecoli.txt ID conversion file provided by Uniprot/Swiss-Prot.

'''
import glob, os, re, requests, sys, urllib2
from collections import Counter
from datetime import date

##Options

##Classes
class ProteinComplex():
	name = "defaultname"

##Functions
def getEcoliIDs():
	baseURL = "http://www.uniprot.org/docs/"
	convfilename = "ecoli.txt"
	convfilepath = baseURL + convfilename
	if os.path.isfile(convfilename): 
		print("Found E. coli ID conversion file on disk: %s" % convfilename)
		id_file = open(convfilename)
	else:
		print("Downloading E. coli ID mapping file from %s" % convfilepath)
		response = urllib2.urlopen(convfilepath)
		id_file = open(os.path.basename(convfilename), "w+b") #Start local file
		chunk = 10240
		while 1:
			data = (response.read(chunk)) #Read 10 Kb at a time
			id_file.write(data)
			if not data:
				print("\n%s file download complete." % convfilename)
				break
			sys.stdout.write(".")
			
	id_file.seek(0)
	capturing = 0	#Don't want the header so we iterate past it
	first_id = "b0001" #For E. coli, this should be the first identifier
	ecoli_ids = {}
	for line in id_file:
		if first_id in line:
			capturing = 1
		if capturing == 1 and not line.strip():	#That is, if we hit an empty line
			break
		if capturing == 1:
			line_content = (line.rstrip()).split()
			#print(line_content)
			bcode = line_content[0][:5]
			if line_content[1][:2] == "JW":
				jwcode = line_content[1]
				uniprotac = line_content[3]
			else:
				jwcode = "NA"
				uniprotac = line_content[2]
			gene_ids = [bcode, jwcode]
			ecoli_ids[uniprotac] = gene_ids
		
	print("Mapped E. coli IDs.")
	id_file.close()
	return ecoli_ids
	
def convertToLong(filename):
	#Converts short complex files (one complex per line)
	#to long form (one complex member per line)
	complex_members = {} #Complex names are keys, members are values
	with open(filename) as shortfile:
		shortfile.readline()	#skip the header
		for line in shortfile:
			line_content = (line.rstrip()).split()
			complex_members[line_content[0]] = line_content[1:]
			
	long_file_name = filename[0:-4] + "_long.txt"
	with open(long_file_name, "w+b") as longfile:
		longfile.write("ProteinID\tComplexMembership\n")
		for complex_name in complex_members:
			for complex_member in complex_members[complex_name]:
				longfile.write(complex_member + "\t" + complex_name + "\n")
	return(long_file_name)
	
def compareSets(filename1, name1, filename2, name2, mode):
	#file1 is the "experimental" file
	#file2 is the model file
	#name variables are used to identify complexes, as their names may
	#be very similar or even just integers.
	#Mode determines if we need to convert interactor IDs.
	
	exp_complexes = {} #Complex names are keys, members are values
	model_complexes = {} #Complex names are keys, members are values
	exp_conservation = {}	#Complex names are keys, conservation values are values
	
	#Load the files and store as dicts
	with open(filename1) as file1:
		file1.readline()	#skip the header
		for line in file1:
			line_content = (line.rstrip()).split()
			complex_name = name1 + "_" + line_content[1]
			if complex_name in exp_complexes:
				exp_complexes[complex_name].append(line_content[0])
			else:
				exp_complexes[complex_name] = [line_content[0]]
	with open(filename2) as file2:
		file2.readline()	#skip the header
		for line in file2:
			line_content = (line.rstrip()).split()
			complex_name = name2 + "_" + line_content[1]
			if complex_name in model_complexes:
				model_complexes[complex_name].append(line_content[0])
			else:
				model_complexes[complex_name] = [line_content[0]]
	
	#Now compare the exp. complexes to the model
	for complex_name in exp_complexes:	#For each complex in the experimental set
		exp_conservation[complex_name] = [0,0]	#No conservation by default
		any_conserved = 0 #Coverage of this complex across the whole model set
		for complex_member in exp_complexes[complex_name]:	#For each complex member
			members_counted = []	#List of members already found, so we don't double-count
			for model_complex in model_complexes:	#Search across all model complexes
				complex_con = 0	#The maximum coverage of any model complex
				con_members_conserved = 0
				if complex_member in model_complexes[model_complex]:	#If the complex member is in the model set
					con_members_conserved = con_members_conserved +1	#Increment its conservation by one
					if complex_member not in members_counted:
						any_conserved = any_conserved +1
					members_counted.append(complex_member)
					complex_con = float(con_members_conserved)/len(model_complexes[model_complex])
					if complex_con > exp_conservation[complex_name][0]:	#Set maximum complex conservation to the largest value
						exp_conservation[complex_name][0] = complex_con
		any_coverage = float(any_conserved)/len(exp_complexes[complex_name])
		exp_conservation[complex_name][1] = any_coverage
	
	#Output for each exp. complex:
	#maximum conservation in a single complex
	#conservation across the complete model set, in any complex		
	compared_file_name = "complexes_compared_" + name1 + "_vs_" + name2 + ".txt"
	with open(compared_file_name, "w+b") as compared_file:
		compared_file.write(name1 + "_Complex\tMaxComplexCon\tSetCon\n")
		for complex_name in exp_conservation:
			compared_file.write("%s\t%5.4f\t%5.4f\n" % (complex_name,
								exp_conservation[complex_name][0],
								exp_conservation[complex_name][1]))
	
	return compared_file_name
	
def compareSpecies(filename1, name1):
	print("Work in progress.")
	compared_file_name = "temp.txt"
	return compared_file_name	

##Main
print("Ready to compare protein complex sets.")
mode = "default"
is_this_ecoli = raw_input("Will you need to convert E. coli locus IDs? Y/N\n")
if is_this_ecoli.lower() == "y":
	getEcoliIDs()
	mode = "ecoli"
else:
	print("OK.")
	
convert_choice = raw_input("Do you need to convert a file to long format? Y/N\n")
if convert_choice.lower() == "y":
	have_short_file = False
	while not have_short_file:
		short_file_name = raw_input("Name of the file to convert?\n")
		if not os.path.isfile(short_file_name):
			print("Couldn't find that file. Try again.")
		else:
			have_short_file = True
	long_file_name = convertToLong(short_file_name)
	print("Long format file is available: %s" % long_file_name)

complex_file_list = glob.glob('complexes*.txt')	
if len(complex_file_list) >2:
	sys.exit("Found more than two complex files. Check for duplicates.")
if len(complex_file_list) == 2:
	filename1 = complex_file_list[0]
	filename2 = complex_file_list[1]
	print("Found two protein complex files:\n%s\n%s" % (filename1, filename2))
	name1 = raw_input("Please provide a short name for the first (experimental) set.\n")
	name2 = raw_input("Please provide a short name for the second (model) set.\n")
	compared_file_name = compareSets(filename1, name1, filename2, name2, mode)
if len(complex_file_list) <2:
	print("Less than two protein complex files found. Provide their names, please.")
	have_file_1 = False
	have_file_2 = False
	while not have_file_1:
		filename1 = raw_input("Name of the experimental complex file?\n")
		if not os.path.isfile(filename1):
			print("Couldn't find that file. Try again.")
		else:
			name1 = raw_input("Please provide a short name for this set.\n")
			have_file_1 = True
	while not have_file_2:
		filename2 = raw_input("Name of the model complex file?\n")
		if not os.path.isfile(filename2):
			print("Couldn't find that file. Try again.")
		elif filename2 == filename1:
			print("Can't compare set to itself. Try again.")
		else:
			name2 = raw_input("Please provide a short name for this set.\n")
			have_file_2 = True
	
model_comparison_choice = raw_input("Do you need to compare the complex set to a model set? Y/N\n")
if model_comparison_choice.lower() == "y":
	model_comparison = True
else:
	model_comparison = False

taxon_comparison_choice = raw_input("Do you need to compare the complex set across species? Y/N\n")
if taxon_comparison_choice.lower() == "y":
	taxon_comparison = True
else:
	taxon_comparison = False

if model_comparison == True:
	compared_file_name = compareSets(filename1, name1, filename2, name2, mode)
	print("Model comparison complete. See %s." % compared_file_name)

if taxon_comparison == True:
	taxcompare_file_name = compareSpecies(filename1, name1)
	print("Broad taxonomic comparison complete. See %s." % taxcompare_file_name)



