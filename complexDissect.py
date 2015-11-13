#!/usr/bin/python
#complexDissect.py
'''
A small utility for comparing sets of protein complexes.

Secondary goal: predict conservation of complexes across species using
protein orthology.

Tertiary goal: allows exploration of specific protein complexes to 
retreive details like best match, size, and functional details.

Requires several Gb of disk space to store mapping files.

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

IN PROGRESS: Enabling broad cross-genome complex conservation prediction.
0. Map components to Uniprot IDs if not done so already
	For E. coli, use ecoli.txt.
1. Map Uniprot IDs of components to eggNOG IDs.
	Uses eggnog4.protein_id_conversion.tsv and members files.
	(Code re-use happens here.)
	Map UPID to eggNOG protein ID, then to NOG ID.
2. Determine all taxids component is present in.
	This information is built into the eggNOG protein ID, so it
	can actually be obtained in the previous step.
	Produce a binary matrix of components (as OGs) vs. taxids.
3. Combine components into their respective complexes.
	Produce a matrix of complexes vs. taxids, where values are 
	each (components present/components in model complex).
4. Use output to produce tree and heatmap. This happens outside the
	software. Model on Caufield et al. 2015 PLoS Comp Bio.

'''
import glob, gzip, os, re, requests, sys, urllib2
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
	#filename1 is the "experimental" file
	#filename2 is the model file
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
	compared_file_name = "compared_complexes_" + name1 + "_vs_" + name2 + ".txt"
	with open(compared_file_name, "w+b") as compared_file:
		compared_file.write(name1 + "_Complex\tMaxComplexCon\tSetCon\n")
		for complex_name in exp_conservation:
			compared_file.write("%s\t%5.4f\t%5.4f\n" % (complex_name,
								exp_conservation[complex_name][0],
								exp_conservation[complex_name][1]))
	
	return compared_file_name
	
def compareSpecies(filename1, name1, id_conversion):
	#filename1 is the "experimental" file
	#name1 is the short name of the set
	#id_conversion is a dictionary with UniprotAC's as keys
	#and values of two types. For E. coli those types are bcode and
	#JW-code, in that order.
	
	#This method requires the eggNOG map files and downloads them if needed
	
	def get_eggnog_maps(): 
		#Download and unzip the eggNOG ID conversion file 
		#Filters file to just Uniprot IDs; the resulting file is the map file.
		#One Uniprot ID may correspond to multiple OGs - e.g. COG1234,COG3810,COG9313. 
		#these cases are considered OGs in their own right as this may indicate a pattern of conserved sequences on its own 
		baseURL = "http://eggnogdb.embl.de/download/latest/"
		convfilename = "eggnog4.protein_id_conversion.tsv.gz "	#File contains ALL database identifiers and corresponding proteins
		
		convfilepath = baseURL + convfilename
		outfilepath = convfilename[0:-3]
		have_mapping_file = 0
		if os.path.isfile(convfilename): 
			print("Found compressed ID conversion file on disk: %s" % convfilename)
			decompress_convfile = 1
			have_mapping_file = 1
		if os.path.isfile(outfilepath): 
			print("Found ID conversion file on disk: %s" % outfilepath)
			decompress_convfile = 0
			have_mapping_file = 1
		if have_mapping_file == 0:
			print("Downloading latest ID mapping file - this file is ~400 Mb compressed so this may take some time.")
			print("Downloading from %s" % convfilepath)
			response = urllib2.urlopen(convfilepath)
			compressed_file = open(os.path.basename(convfilename), "w+b") #Start local compressed file
			chunk = 1048576
			while 1:
				data = (response.read(chunk)) #Read one Mb at a time
				compressed_file.write(data)
				if not data:
					print("\n%s file download complete." % convfilename)
					compressed_file.close()
					break
				sys.stdout.write(".")
			decompress_convfile = 1
			
		if decompress_convfile == 1:
			print("Decompressing map file. Lines written, in millions:")
			#Done in chunks since it's a large file
			with gzip.open(convfilename) as infile: #Open that compressed file, read and write to uncompressed file
				outfile = open(outfilepath, "w+b")
				linecount = 0
				for line in infile:
					outfile.write(line)
					linecount = linecount +1
					if linecount % 100000 == 0:
							sys.stdout.write(".")
					if linecount % 1000000 == 0:
							sys.stdout.write(str(linecount/1000000))
				infile.close()
			newconvfilename = outfilepath
			outfile.close()
		
		#Download and decompress member NOG files (2 of them)
		nogURL = baseURL + "data/NOG/"
		nogfilename = "NOG.members.tsv.gz"
		bactnogURL = baseURL + "data/bactNOG/"
		bactnogfilename = "bactNOG.members.tsv.gz" 
		all_nog_locations = [[nogURL, nogfilename], [bactnogURL, bactnogfilename]]
		
		for location in all_nog_locations:
			baseURL = location[0]
			memberfilename = location[1]
			memberfilepath = baseURL + memberfilename
			outfilepath = memberfilename[0:-3]
			if os.path.isfile(memberfilename): 
				print("\nFound compressed NOG membership file on disk: %s" % memberfilename)
				decompress_memberfile = 1
			if os.path.isfile(outfilepath): 
				print("\nFound NOG membership file on disk: %s" % outfilepath)
				decompress_memberfile = 0
			else:
				print("\nDownloading NOG membership file - this may take some time.")
				print("Downloading from %s" % memberfilepath)
				response = urllib2.urlopen(memberfilepath)
				compressed_file = open(os.path.basename(memberfilename), "w+b") #Start local compressed file
				chunk = 1048576
				while 1:
					data = (response.read(chunk)) #Read one Mb at a time
					compressed_file.write(data)
					if not data:
						print("\n%s file download complete." % memberfilename)
						compressed_file.close()
						break
					sys.stdout.write(".")
				decompress_memberfile = 1
				
			if decompress_memberfile == 1:
				print("Decompressing NOG membership file %s" % memberfilename)
				#Done in chunks since it's a large file
				with gzip.open(memberfilename) as infile: #Open that compressed file, read and write to uncompressed file
					outfile = open(outfilepath, "w+b")
					linecount = 0
					for line in infile:
						outfile.write(line)
						linecount = linecount +1
						if linecount % 100000 == 0:
							sys.stdout.write(".")
						if linecount % 1000000 == 0:
							sys.stdout.write(str(linecount/1000000))
					infile.close()
				outfile.close()
				
		#Clean up by removing compressed files
		
		all_compressed_files = [convfilename, nogfilename, bactnogfilename]
		for filename in all_compressed_files:
			if os.path.isfile(filename):
				print("\nRemoving compressed file: %s" % filename)
				os.remove(filename)
		
		#Load and filter the ID conversion file as dictionary
		print("Parsing ID conversion file. Lines read, in millions:")
		with open(convfilename[0:-3]) as infile:
			id_dict = {}	#Dictionary of eggNOG protein IDs with database IDs as keys
			#Gets filtered down to relevant database IDs (i.e., Uniprot IDs)
			linecount = 0
			for line in infile:
				linecount = linecount +1
				line_raw = ((line.rstrip()).split("\t"))	#Protein IDs are split for some reason; merge them
				one_id_set = [line_raw[0] + "." + line_raw[1], line_raw[2], line_raw[3]]
				if "UniProt_AC" in one_id_set[2]:
					id_dict[one_id_set[1]] = one_id_set[0]
				if linecount % 100000 == 0:
					sys.stdout.write(".")
				if linecount % 1000000 == 0:
					sys.stdout.write(str(linecount/1000000))
			infile.close()
	
		#Use filtered ID conversion input to map to NOG members
		print("\nReading NOG membership files.")
		all_nog_filenames = [nogfilename[0:-3], bactnogfilename[0:-3]]
		nog_members = {}	#Dictionary of NOG ids with protein IDs as keys (need to split entries for each)
		nog_count = 0
		og_taxid_map = {} #OG ids are keys, all taxids they have members in are values
		for filename in all_nog_filenames:
			temp_nog_members = {}	#We will have duplicates within each set but don't want to lose the information.
			print("Reading from %s" % filename)
			with open(filename) as infile:
				for line in infile:
					nog_count = nog_count +1
					line_raw = ((line.rstrip()).split("\t"))
					nog_id = line_raw[1]
					line_members = line_raw[5].split(",")
					og_taxid_map[nog_id] = []
					for protein_id in line_members:			#The same protein could be in more than one OG at the same level
						taxid = (protein_id.split("."))[0]
						if taxid not in og_taxid_map[nog_id]:
							og_taxid_map[nog_id].append(taxid)
						if protein_id in temp_nog_members:
							temp_nog_members[protein_id] = temp_nog_members[protein_id] + "," + nog_id
						else:
							temp_nog_members[protein_id] = nog_id
				infile.close()
			nog_members.update(temp_nog_members)
		
		upids_length = str(len(id_dict))
		nogs_length = str(nog_count)
		proteins_length = str(len(nog_members))
		
		print("Mapping %s Uniprot IDs to %s NOGs through %s eggNOG protein IDs:" % (upids_length, nogs_length, proteins_length))
		upid_to_NOG = {}	#Conversion dictionary. Values are OGs, keys are UPIDs.
		mapped_count = 0	#upids mapped to nogs.
		for upid in id_dict:
			if id_dict[upid] in nog_members:
				upid_to_NOG[upid] = nog_members[id_dict[upid]]
				mapped_count = mapped_count +1
				if mapped_count % 100000 == 0:
					sys.stdout.write(".")
				if mapped_count % 1000000 == 0:
					sys.stdout.write(str(mapped_count/1000000))
			
		#Use this mapping to build map file, named "uniprot_og_maps_*.txt"
		print("Writing map file.")
		nowstring = (date.today()).isoformat()
		mapfilename = "uniprot_og_maps_" + nowstring + ".txt"
		mapfile = open(mapfilename, "w+b")
		for mapping in upid_to_NOG:
			mapfile.write(mapping + "\t" + upid_to_NOG[mapping] + "\n")	#Each line is a uniprot ID and an OG id
		mapfile.close()
		
		#Use the taxids to build a map file for OGs to each taxid they have a member in, named "og_to_taxid_*.txt"
		print("Writing taxonomy matrix file.")
		taxonfilename = "og_to_taxid_" + nowstring + ".txt"
		taxonfile = open(taxonfilename, "w+b")
		for mapping in og_taxid_map:
			taxonfile.write(mapping + "\t" + " ".join(og_taxid_map[mapping]) + "\n")	#Each line is an OG then a list of taxids 
			#Taxids separated by space
		taxonfile.close()
	
	print("***Species comparison: Work in progress.***")
	component_con_file_name = name1 + "_component_conservation.txt"
	cplx_con_file_name = name1 + "_complex_conservation.txt"
	
	exp_complexes = {} #Complex names are keys, members are values
	uniprot_to_og = {} #ALL available Uniprot IDs are keys, members are eggNOG OG IDs
		#Can only include UPIDs mapped to eggNOG IDs
	og_to_taxid = {} #All the taxids this OG has members in, as per the member files
	unified_components = []		#All protein complex components, as UPIDs
	component_ogs = {}	#Component UPIDs are keys, members are corresponding OGs
	unmapped_components = []	#Any UPID not found in the eggNOG ID conversion file
	og_list = []	#Just the unique OGs in use
	all_taxids = [] #All taxids in use. Don't care about parents or children here.
	all_cplx_taxids = [] #All taxids relevant to the complex components.
	ready_to_map = False
	
	while not ready_to_map:
		map_file_list = glob.glob("uniprot_og_maps_*.txt")
		taxon_file_list = glob.glob("og_to_taxid_*.txt")
		if len(map_file_list) >0 and len(taxon_file_list) >0:
			if len(map_file_list) == 1:
				print("Found map file %s on disk." % map_file_list[0])
				print("Loading file...")
				with open(map_file_list[0]) as id_map_file:
					linecount = 0
					for line in id_map_file:
						linecount = linecount +1
						if linecount % 10000 == 0:
							sys.stdout.write(".")
						line_content = (line.rstrip()).split("\t")
						uniprot_to_og[line_content[0]] = line_content[1]
			else:
				sys.exit("Found more than one map file on disk. Please check for duplicates.")
			if len(taxon_file_list) == 1:
				print("\nFound OG vs taxon file %s on disk." % taxon_file_list[0])
				print("Loading file...")
				with open(taxon_file_list[0]) as taxon_file:
					linecount = 0
					for line in taxon_file:
						linecount = linecount +1
						if linecount % 10000 == 0:
							sys.stdout.write(".")
						line_content = (line.rstrip()).split("\t")
						og_to_taxid[line_content[0]] = line_content[1].split(" ")
						for taxid in line_content[1].split(" "):
							if taxid not in all_taxids:
								all_taxids.append(taxid)
			else:
				sys.exit("\nFound more than one OG vs taxon file on disk. Please check for duplicates.")
			ready_to_map = True
		else:
			print("A protein map or a taxon file is missing. Rebuilding them...")
			get_eggnog_maps()
	
	with open(filename1) as file1:
		file1.readline()	#skip the header
		for line in file1:
			line_content = (line.rstrip()).split()
			complex_name = name1 + "_" + line_content[1]
			if complex_name in exp_complexes:
				exp_complexes[complex_name].append(line_content[0])
			else:
				exp_complexes[complex_name] = [line_content[0]]
	
	for name in exp_complexes:
		for component in exp_complexes[name]:
			for identifier in id_conversion:
				expanded_ids = id_conversion[identifier]
				if component in expanded_ids:
					#print(component + "\t" + identifier)
					unified_components.append(identifier)
	
	print("\nSearching for conservation of %s unique proteins." % len(unified_components))
	
	for component in unified_components:
		if component in uniprot_to_og:
			new_og = uniprot_to_og[component]
			component_ogs[component] = new_og
			if new_og not in og_list:
				og_list.append(new_og)
		else:
			component_ogs[component] = component
			unmapped_components.append(component)
	
	for og in og_list:
		if og in og_to_taxid:
			for taxid in og_to_taxid[og]:
				if taxid not in all_cplx_taxids:
					all_cplx_taxids.append(taxid)
					
	print("Mapped complex components to %s OGs." % len(og_list))
	print("%s complex components did not map to OGs." % len(unmapped_components))
	print("Orthologs of these components are found across %s taxids." % len(all_cplx_taxids)) 
	
	with open(component_con_file_name, "w+b") as compared_file:
		compared_file.write("\t" + "\t".join(all_cplx_taxids) + "\n")
		for og in og_list:
			compared_file.write(og + "\n")
		
	with open(cplx_con_file_name, "w+b") as compared_file:
		compared_file.write("\t" + "\t".join(all_cplx_taxids) + "\n")
		for complex_name in exp_complexes:
			compared_file.write(complex_name + "\n")

	compared_file_names = [component_con_file_name, cplx_con_file_name]
	return compared_file_names	

##Main
print("Ready to compare protein complex sets.")
mode = "default"
is_this_ecoli = raw_input("Will you need to convert E. coli locus IDs? Y/N\n")
if is_this_ecoli.lower() == "y":
	id_conversion = getEcoliIDs()
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
	
model_comparison_choice = raw_input("Compare the experimental complex set to the model set? Y/N\n")
if model_comparison_choice.lower() == "y":
	name2 = raw_input("Please provide a short name for the second (model) set.\n")
	compared_file_name = compareSets(filename1, name1, filename2, name2, mode)
	print("Model comparison complete. See %s." % compared_file_name)
else:
	model_comparison = False

taxon_comparison_choice = raw_input("Compare the experimental complex set across species? Y/N\n")
if taxon_comparison_choice.lower() == "y":
	taxcompare_file_names = compareSpecies(filename1, name1, id_conversion)
	print("Broad taxonomic comparison complete. \n" + 
		"See %s for component conservation and %s for complex conservation."
		% (taxcompare_file_names[0], taxcompare_file_names[1]))
else:
	taxon_comparison = False

