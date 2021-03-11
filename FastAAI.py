#!/usr/bin/env python

"""
########################################################################
# Authors:	   Carlos Ruiz, Kenji Gerhardt
# Intitution:   Georgia Institute of Technology
# Version:	  1.0.0
# Date:		 March 11, 2021

# Description: Calculates the average amino acid identity using k-mers
from single copy genes. It is a faster and more accurate version of the 
regular AAI (Blast or Diamond) and the hAAI implemented in MiGA.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import subprocess, argparse, multiprocessing, datetime, shutil
import textwrap, pickle, gzip
import numpy as np
from tempfile import TemporaryDirectory
from random import randint
from pathlib import Path
from sys import argv
from sys import exit
import os
from itertools import islice
from functools import partial
import copy
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import CountVectorizer



################################################################################
"""---1.0 Define Functions---"""
# --- Run prodigal ---
# ------------------------------------------------------
def run_prodigal(input_file):
	"""
	Runs prodigal, compares translation tables and stores faa files

	Arguments:
	   input_file -- Path to genome FastA file
	
	Returns:
		output -- Path to amino acid fasta result
	"""
	# Predict proteins with translation tables 4 and 11
	file_path = Path(input_file)
	filename = file_path.name
	
	base_name = os.path.splitext(os.path.basename(str(filename)))[0]
	
	folder = Path(output_root_directory + "/predicted_proteins/")
	protein_output = folder / (base_name + '.faa')
	output_11 = folder / (base_name + '.faa.11')
	temp_output = folder / (base_name + '.temp')
	
	subprocess.call(["prodigal", "-i", str(file_path), "-a", str(output_11), 
					"-p", "meta", "-q", "-o", str(temp_output)])
					
	output_4 = folder / (base_name + '.faa.4')
	temp_output = folder / (base_name + '.temp')
	subprocess.call(["prodigal", "-i", str(file_path), "-a", str(output_4), 
					"-p", "meta", "-g", "4", "-q", "-o", str(temp_output)])

	# Compare translation tables
	length_4 = 0
	length_11 = 0
	with open(output_4, 'r') as table_4:
		for line in table_4:
			if line.startswith(">"):
				continue
			else:
				length_4 += len(line.strip())

	with open(output_11, 'r') as table_11:
		for line in table_11:
			if line.startswith(">"):
				continue
			else:
				length_11 += len(line.strip())
	
	#We can get rid of the temp file immediately, we won't be using it again
	temp_output.unlink()
	
	#Select the winning translation table and remove the other immediately. Record the winner
	if (length_4 / length_11) >= 1.1:
		output_11.unlink()
		chosen_protein = open(output_4, 'r')
		table_11 = False
	else:
		output_4.unlink()
		chosen_protein = open(output_11, 'r')
		table_11 = True
		
	destination = open(protein_output, "w")
	
	for line in chosen_protein:
		if line.startswith(">"):
			destination.write("{}".format(line))
		else:
			line = line.replace('*', '')
			destination.write("{}".format(line))
	
	destination.close()
	chosen_protein.close()
	
	# Remove the winning intermediate file.
	if table_11:
		output_11.unlink()
	else:
		output_4.unlink()

	return str(protein_output)
# ------------------------------------------------------

# --- Run prodigal for viruses ---
# ------------------------------------------------------
def run_prodigal_virus(input_file):
	"""
	Runs prodigal, compares translation tables and stores faa files

	Arguments:
	   input_file -- Path to genome FastA file
	
	Returns:
		output -- Path to amino acid fasta result
	"""
	# Predict proteins with translation tables 4 and 11
	file_path = Path(input_file)
	filename = file_path.name
	
	base_name = os.path.splitext(os.path.basename(str(filename)))[0]
	
	folder = Path(output_root_directory + "/predicted_proteins/")
	#Here we have no selection between tables 4 and 11 so we just create an intermediate.
	intermediate_protein_output = folder / (base_name + '.intermediate.faa')
	#And here we have the intended output name after cleaning.
	final_protein_output = folder / (base_name + '.faa')
	temp_output = folder / (base_name + '.temp')
	subprocess.call(["prodigal", "-i", str(file_path), "-a", str(intermediate_protein_output), 
					"-p", "meta", "-q", "-o", str(temp_output)])
	
	# Remove intermediate files
	temp_output.unlink()
	
	chosen_protein = open(intermediate_protein_output, 'r')
	destination = open(final_protein_output, "w")
	
	for line in chosen_protein:
		if line.startswith(">"):
			destination.write("{}".format(line))
		else:
			line = line.replace('*', '')
			destination.write("{}".format(line))
	
	destination.close()
	chosen_protein.close()

	intermediate_protein_output.unlink()

	return str(protein_output)
# ------------------------------------------------------

# --- Run hmmsearch ---
# ------------------------------------------------------
def run_hmmsearch(input_file):
	"""
	Runs hmmsearch on the set of SCGs and select the
	best Archaea or Bacterial model
	
	Arguments:
		input_file -- Path to protein FastA file
	
	Returns:
		output -- Path to hmmsearch hits table
	"""
	file_path = Path(input_file)
	folder = Path(output_root_directory + "/hmms/")
	filename = file_path.name
	
	base_name = os.path.splitext(os.path.basename(str(filename)))[0]
	
	hmm_output = folder / (base_name + '.hmm')
	temp_output = folder / (base_name + '.temp')
	script_path = Path(__file__)
	script_dir = script_path.parent
	hmm_complete_model = script_dir / "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"
	subprocess.call(["hmmsearch", "--tblout", str(hmm_output), "-o", str(temp_output), "--cut_tc", "--cpu", "1",
					str(hmm_complete_model), str(file_path)])
	temp_output.unlink()
	return str(hmm_output)
# ------------------------------------------------------

# --- Filter HMM results for best matches ---
# ------------------------------------------------------
def hmm_filter(scg_hmm_file):
	"""
	Filters HMM results for best hits per protein
	
	Arguments:
		SCG_HMM_file {file path} -- Path to HMM results file
	
	Returns:
		outfile -- Path to filtered files
	"""
	hmm_path = Path(scg_hmm_file)
	filename = hmm_path.name
	
	base_name = os.path.splitext(os.path.basename(str(filename)))[0]
	
	folder = Path(output_root_directory + "/filtered_hmms/")
	outfile = folder / (base_name + '.hmm.filt')
	hmm_hit_dict = {}
	with open(scg_hmm_file, 'r') as hit_file:
		for line in hit_file:
			if line.startswith("#"):
				continue
			else:
				hit = line.strip().split()
				protein_name = hit[0]
				score = float(hit[8])
				if protein_name in hmm_hit_dict:
					if score > hmm_hit_dict[protein_name][0]:
						hmm_hit_dict[protein_name] = [score, line]
					elif score < hmm_hit_dict[protein_name][0]:
						continue
					else:
						if randint(2) > 0:
							hmm_hit_dict[protein_name] = [score, line]
				else:
					hmm_hit_dict[protein_name] = [score, line]
	with open(outfile, 'w') as output:
		for hits in hmm_hit_dict.values():
			output.write("{}".format(hits[1]))
	return str(outfile)
# ------------------------------------------------------

# --- Find Kmers from HMM results ---
# ------------------------------------------------------
def kmer_extract(input_files):
	"""
	Extract kmers from protein files that have hits
	in the HMM searches.
	
	Arguments:
		SCG_HMM_file {file path} -- Path to filtered HMM results.
	
	Returns:
		[genome_kmers] -- Dictionary of kmers per gene. 
	"""
	final_filename = input_files[0]
	protein_file = input_files[1]
	scg_hmm_file = input_files[2]
	
	#inputs now include the folder for the query.
	scg_hmm_file = scg_hmm_file.replace("/hmms/", "/filtered_hmms/")
	
	positive_matches = {}
	positive_proteins = []
	with open(scg_hmm_file, 'r') as hmm_input:
		for line in hmm_input:
			line = line.strip().split()
			protein_name = line[0]
			model_name = line[3]
			score = line[8]
			if model_name in positive_matches:
				if score > positive_matches[model_name][1]:
					positive_matches[model_name] = [protein_name, score]
				else:
					continue
			else:
				positive_matches[model_name] = [protein_name, score]
	for proteins in positive_matches.values():
		positive_proteins.append(proteins[0])
	scg_kmers = read_kmers_from_file(protein_file, positive_proteins, 4)
	for accession, protein in positive_matches.items():
		scg_kmers[accession] = scg_kmers.pop(protein[0])
	genome_kmers = {final_filename : scg_kmers}
	return genome_kmers
# ------------------------------------------------------

# --- Extract kmers from protein sequences ---
# ------------------------------------------------------
def read_kmers_from_file(filename, positive_hits, ksize):
	scg_kmers = {}
	store_sequence = False
	protein_name = ""
	protein_sequence = ""
	with open(filename) as fasta_in:
		for line in fasta_in:
			if line.startswith(">"):
				if store_sequence == True:
					kmers = build_kmers(protein_sequence, ksize)
					scg_kmers[protein_name] = kmers
				protein_sequence = ""
				store_sequence = False
				line = line.replace(">", "")
				protein_name = line.strip().split()[0]
				if protein_name in positive_hits:
					store_sequence = True
			else:
				if store_sequence == True:
					protein_sequence += line.strip()
				else:
					continue
			if store_sequence == True:
				kmers = build_kmers(protein_sequence, ksize)
				scg_kmers[protein_name] = kmers
	return scg_kmers
# ------------------------------------------------------

# --- Extract kmers from viral protein sequences ---
# ------------------------------------------------------
def read_viral_kmers_from_file(input_information):
	final_filename = input_information[0]
	protein_file = input_information[1]
	kmer_size = input_information[2]
	
	scg_kmers = set()
	protein_sequence = ""
	store_sequence = False
	number_of_proteins = 0
	with open(protein_file) as fasta_in:
		for line in fasta_in:
			if line.startswith(">"):
				number_of_proteins += 1
				if store_sequence == True:
					kmers = build_viral_kmers(protein_sequence, kmer_size)
					scg_kmers.update(kmers)
					protein_sequence = ""
				else:
					protein_sequence = ""
					store_sequence = True
			else:
				protein_sequence += line.strip()
			if store_sequence == True:
				kmers = build_viral_kmers(protein_sequence, kmer_size)
				scg_kmers.update(kmers)
	genome_kmers = {final_filename : [number_of_proteins, ','.join(list(scg_kmers))]}
	return genome_kmers
# ------------------------------------------------------

# --- Build Kmers ---
# ------------------------------------------------------
def build_kmers(sequence, ksize):
	kmers = []
	n_kmers = len(sequence) - ksize + 1

	for i in range(n_kmers):
		kmer = sequence[i:i + ksize]
		kmers.append(kmer)
	kmers_set = ','.join(set(kmers))
	return kmers_set
# ------------------------------------------------------

# --- Build Viral Kmers ---
# ------------------------------------------------------
def build_viral_kmers(sequence, ksize):
	kmers = []
	n_kmers = len(sequence) - ksize + 1

	for i in range(n_kmers):
		kmer = sequence[i:i + ksize]
		kmers.append(kmer)
	kmers_set = set(kmers)
	return kmers_set
# ------------------------------------------------------

# --- Create global dictionary with unique kmers and indices for each one ---
# ------------------------------------------------------

#Function for splitting a range into groups of consecutive n-sized chunks. Covers the whole range, even when a chunk is undersized relative to the others.
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield [i, i + n]

#Function for processing a chunk of a kmer dict's strings into indices.
def dissect_kmers(args):
	index_start, index_end = args[0], args[1]
	genomes = list(current_global_dict.keys())[index_start:index_end]
	kmers = {}
	for g in genomes:
		for marker_protein_id in current_global_dict[g]:
			kmer_list = current_global_dict[g][marker_protein_id].split(',')
			for kmer in kmer_list:
				try:
					kmers[kmer]
				except:
					kmers[kmer] = True
	
	return set(list(kmers.keys()))
	
def global_unique_kmers_parallel(kmer_dictionaries, threads):
	"""
	Extract every kmer in the whole dataset
	Create global dictionary with unique kmers and indices for each one
	
	Arguments:
		kmer_dict {dict} -- Dictionary with kmers for each marker protein per input file
	
	Returns:
		[global_kmer_index_dictionary] -- Dictionary with a unique index per kmer
	"""
	# Make this dictionary global regardless of quer == reference or not
	print("Indexing unique kmers... ", end = '', flush = True)
	
	global global_kmer_index_dictionary
	global_kmer_index_dictionary = {}
	
	global current_global_dict
	
	#initialize unique kmer set
	base_kmers = set()
	
	#Process query, query and ref dict as needed.
	for kd in kmer_dictionaries:
		#Make the current dict to parse global for the parallel procs
		current_global_dict = kd
		
		#Split the current dict into chunks, 1/thread
		total_genomes = len(current_global_dict)
		chunk_size = total_genomes // threads
		chunk_indices = chunks(range(0, total_genomes), chunk_size)
		
		#Parallel to the chunks.
		pool = multiprocessing.Pool(threads)
		collections = pool.map(dissect_kmers, chunk_indices)
		pool.close()
		pool.join()
		
		#merge kmer results into the parent set, no dupes.
		base_kmers = base_kmers.union(*collections)
			
		#Cleanup.
		del current_global_dict
	
	#Create a dictionary with unique indices
	counter = 0
	for kmer in base_kmers:
		try:
			global_kmer_index_dictionary[kmer]
		except:
			global_kmer_index_dictionary[kmer] = counter
			counter += 1	
	
	print("done!", len(global_kmer_index_dictionary), "unique kmers found.")
	#We made it global, so this is OK
	return None
	
#Old, serial version
def global_unique_kmers(kmer_dictionaries):
	"""
	Extract every kmer in the whole dataset
	Create global dictionary with unique kmers and indices for each one
	
	Arguments:
		kmer_dict {dict} -- Dictionary with kmers for each marker protein per input file
	
	Returns:
		[global_kmer_index_dictionary] -- Dictionary with a unique index per kmer
	"""
	# Make this dictionary global regardless of quer == reference or not
	print("Indexing unique kmers")
	global global_kmer_index_dictionary
	global_kmer_index_dictionary = {}
	
	counter = 0
	for kmer_dict in kmer_dictionaries:
		for marker_protein_id in kmer_dict.values():
			for kmer_list in marker_protein_id.values():
				kmer_list = kmer_list.split(',')
				for kmer in kmer_list:
					try:
						global_kmer_index_dictionary[kmer]
					except:
						global_kmer_index_dictionary[kmer] = counter
						counter += 1
# ------------------------------------------------------

# --- Create global viral dictionary with unique kmers and indices for each one ---
# ------------------------------------------------------
#Function for processing a chunk of a kmer dict's strings for viral genomes
def dissect_kmers_viral(args):
	index_start, index_end = args[0], args[1]
	genomes = list(current_global_dict.keys())[index_start:index_end]
	kmers = {}
	for g in genomes:
		#TODO - make sure this actually works.
		kmer_list = current_global_dict[g]
		for kmer in kmer_list[1].split(','):
			try:
				kmers[kmer]
			except:
				kmers[kmer] = True
					
	return set(list(kmers.keys()))

#Function for extracting the unique kmers from a set of genomes.
def global_unique_viral_kmers_parallel(kmer_dictionaries, threads):
	"""
	Extract every kmer in the whole dataset
	Create global dictionary with unique kmers and indices for each one
	
	Arguments:
		kmer_dict {dict} -- Dictionary with kmers for each marker protein per input file
	
	Returns:
		[global_kmer_index_dictionary] -- Dictionary with a unique index per kmer
	"""
	# Make this dictionary global regardless of quer == reference or not
	print("Indexing unique kmers")
	global global_kmer_index_dictionary
	global_kmer_index_dictionary = {}
	
	global current_global_dict
	
	base_kmers = set()
	
	for kd in kmer_dictionaries:
		current_global_dict = kd
		
		total_genomes = len(current_global_dict)
		chunk_size = total_genomes // threads 
		chunk_indices = chunks(range(0, total_genomes), chunk_size)
		
		#Parallel to the chunks...
		pool = multiprocessing.Pool(threads)
		collections = pool.map(dissect_kmers_viral, chunk_indices)
		pool.close()
		pool.join()
		
		
		base_kmers = base_kmers.union(*collections)

	counter = 0
	for kmer in base_kmers:
		try:
			global_kmer_index_dictionary[kmer]
		except:
			global_kmer_index_dictionary[kmer] = counter
			counter += 1
	
	del current_global_dict
	#We made it global, so this is OK
	return None

#old, serial version
def global_unique_viral_kmers(kmer_dictionaries):
	"""
	Extract every kmer in the whole dataset
	Create global dictionary with unique kmers and indices for each one
	
	Arguments:
		kmer_dict {dict} -- Dictionary with kmers for each marker protein per input file
	
	Returns:
		[global_kmer_index_dictionary] -- Dictionary with a unique index per kmer
	"""
	# Make this dictionary global regardless of quer == reference or not
	print("Indexing unique kmers")
	global global_kmer_index_dictionary
	global_kmer_index_dictionary = {}
	counter = 0
	for kmer_dict in kmer_dictionaries:
		for kmer_list in kmer_dict.values():
			for kmer in kmer_list[1].split(','):
				try:
					global_kmer_index_dictionary[kmer]
				except:
					global_kmer_index_dictionary[kmer] = counter
					counter += 1
# ------------------------------------------------------

# --- Convert kmers to indices ---
# ------------------------------------------------------
def convert_kmers_to_indices(kmer_dict):

#global_kmer_index_dictionary already exists out in the ether
	for genome in kmer_dict:
		for protein_marker in kmer_dict[genome]:
			kmer_index = []
			for kmer in kmer_dict[genome][protein_marker].split(','):
				kmer_index.append(global_kmer_index_dictionary[kmer])
			kmer_index = np.sort(np.unique(np.array(kmer_index, dtype=np.int32)))
			kmer_dict[genome][protein_marker] = kmer_index

	return kmer_dict
# ------------------------------------------------------

# --- Convert viral kmers to indices ---
# ------------------------------------------------------
#Same as above.
def convert_viral_kmers_to_indices(kmer_dict):
	for genome in kmer_dict:
		kmer_index = []
		for kmer in kmer_dict[genome][1].split(','):
			kmer_index.append(global_kmer_index_dictionary[kmer])
		kmer_index = np.sort(np.unique(np.array(kmer_index, dtype=np.int32)))
		kmer_dict[genome][1] = kmer_index

	return kmer_dict
# ------------------------------------------------------

# --- Transform kmer dictionaries to index dictionaries ---
# ------------------------------------------------------
def smart_arguments(kmer_dict, temporary_working_directory, single_dataset):
	#kmer_dict = convert_kmers_to_indices(kmer_dict)
	#Get skip indices
	smartargs = []
	genome_ids = list(kmer_dict.keys())
	for i in range(0, len(genome_ids)):
		if single_dataset == True:
			smartargs.append((temporary_working_directory, genome_ids[i], i))
		else:
			smartargs.append((temporary_working_directory, genome_ids[i]))
		
	return smartargs
# ------------------------------------------------------

# --- Transform viral kmer dictionaries to index dictionaries ---
# ------------------------------------------------------
def smart_arguments_viral(kmer_dict, temporary_working_directory, single_dataset):
	#kmer_dict = convert_viral_kmers_to_indices(kmer_dict)
	#Get skip indices
	smartargs = []
	genome_ids = list(kmer_dict.keys())
	for i in range(0, len(genome_ids)):
		if single_dataset == True:
			smartargs.append((temporary_working_directory, genome_ids[i], i))
		else:
			smartargs.append((temporary_working_directory, genome_ids[i]))
		
	return smartargs
# ------------------------------------------------------

# --- Parse kAAI when query == reference using slower memory conservative algorithm ---
# ------------------------------------------------------
def single_kaai_parser(arguments):
	"""
	Calculates the Jaccard distances using single protein markers shared by two genomes
	
	Arguments:
		arguments {tuple} -- Tuple with the temporary folder, the query id and the index of said query_id
	
	Returns:
		[Path to output] -- Path to output file
	"""
	temporary_folder = arguments[0]
	query_id = arguments[1]
	skip_first_n = arguments[2]
	
	temporary_folder = Path(str(temporary_folder.name))
	temporary_file = Path(query_id).name + '.faai.temp'
	temporary_output = temporary_folder / temporary_file
	
	query_scg_list = np.array(list(query_kmer_dictionary[query_id].keys()))

	with open(temporary_output, 'w') as out_file:
		#for target_genome, scg_ids in query_kmer_dictionary.items():
		for target_genome in list(query_kmer_dictionary.keys())[skip_first_n:]:
			# Get number and list of SCG detected in reference
			target_scg_list = np.array(list(query_kmer_dictionary[target_genome].keys()))
			shorter_genome = min(len(query_scg_list), len(target_scg_list))
			#If self, 1.0 similarity.
			if query_id == target_genome:
					out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
					1.0, 0.0, len(query_scg_list), len(target_scg_list), 100))
					continue
			
			jaccard_similarities = []
			# Get shared proteins (scgs)
			final_scg_list = np.intersect1d(query_scg_list, target_scg_list)
			# Extract a list of kmers for each SCG in the list
			query_kmer_list = list(map(query_kmer_dictionary[query_id].get, final_scg_list))
			reference_kmer_list = list(map(query_kmer_dictionary[target_genome].get, final_scg_list))
			# Calculate the jaccard index
			for accession in range(len(query_kmer_list)):
				union = len(np.union1d(query_kmer_list[accession], reference_kmer_list[accession]))
				intersection = len(query_kmer_list[accession]) + len(reference_kmer_list[accession]) - union
				jaccard_similarities.append(intersection / union)
			
			
			if len(jaccard_similarities) > 0:
				jaccard_similarities = np.array(jaccard_similarities, dtype=np.float_)
				try:
					mean = np.mean(jaccard_similarities)
					var = np.std(jaccard_similarities)
					
					if 0.0058 < mean < 0.843:
						aai_est = round(kaai_to_aai(mean), 2)
					if mean <= 0.0058:
						aai_est = "<30%"
					if mean >= 0.843:
						aai_est = ">90%"
					
					
					out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
							round(mean, 4), round(var, 4),
							len(jaccard_similarities), shorter_genome, aai_est))
				except:
					out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
							"NA", "NA", "NA", "NA", "NA"))
			else:
				out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
							"NA", "NA", "NA", "NA", "NA"))
												
	return temporary_output
# ------------------------------------------------------


# --- Parse kAAI when query == reference using faster, high mem algorithm---
# ------------------------------------------------------
def single_kaai_parser_fast(arguments):

	temporary_folder = arguments[0]
	query_id = arguments[1]
	skip_first_n = arguments[2]
	
	temporary_folder = Path(str(temporary_folder.name))
	temporary_file = Path(genome_index[query_id]).name + '.faai.temp'
	temporary_output = temporary_folder / temporary_file
	
	valid_genome_list = np.sort(np.array(list(kmer_dict.keys())[skip_first_n:], dtype = np.int32))
	
	jaccard_lists_per_genome = {}
	for gen in valid_genome_list:
		jaccard_lists_per_genome[gen] = []
	
	genome_intersection_counts = {}
	for gen in valid_genome_list:
		genome_intersection_counts[gen] = 0
		
	#We use this for keeping track of the number of times a genome appeared in protein tracking, which should be the number of intersecting protein families per genome.
	#number_of_prots_per_genome = genome_intersection_counts.copy()
	
	query_protein_ct = len(kmer_dict[query_id])
	
	#Go through the set of proteins and kmers for each protein in the query genome
	for protein in kmer_dict[query_id]:
		
		#We want to reinitialize the dict each time without recreating it from scratch
		non_ref = genome_intersection_counts.copy()
		
		#size of union = size of query protein kmers + size of target protein kmers - size of intersecting kmers
		query_union_length = len(kmer_dict[query_id][protein])
		
		#This loop handles one full protein, which is one item in our jaccard index. This handles the protein for all genomes, though.
		for kmer in kmer_dict[query_id][protein]:
			#This gets us down to genomes that are actually found
			for genome in intersect1d_searchsorted(inverted_dict[kmer][protein], valid_genome_list):
				#They have to be sharing the protein at this point, so we can use the single protein key.
				non_ref[genome] += 1
					
			
		#We still have the protein referenced here, so we can reuse it.
		for genome in non_ref:
			
			try: 
				#Non_ref has genomes that have a count of zero, which means the protein we're looking at never showed up in it.
				#That means trying to access that protein for that genome will fail, so if it does, we want to continue but not modify the 0 size intersection we pre-set.
				#If it succeeds, we get it to do some work we want done anyway.
				target_union_length = len(kmer_dict[genome][protein])
			except:
				continue
			
			#number_of_prots_per_genome[genome] += 1
			
			jaccard_lists_per_genome[genome].append((non_ref[genome] / (query_union_length + target_union_length - non_ref[genome])))
			
	#Open file here - access the name of the genome using the dict...
	
	handle = open(temporary_output, "w")
	
	#Write contents of jaccards per genome pairing.
	#We still need to get the number of kmers in the query vs the target for each protein, here...
	for idx in jaccard_lists_per_genome:
		
		#this is the number of possible protein matches with the query genome - between 0 and the shorter lsit of prots
		possible_matches = min(query_protein_ct, len(kmer_dict[idx]))
		
		#For each genome, this is the list of the jaccard index of each protein [from 0 to possible matches]
		if query_id == idx:
			#Case: self
			jaccard_mean = 1.0
			jaccard_stddev = 0.0
			aai_est = 100
			#This is the number of observed protein matches with the query genome
			number_intersecting_prots = possible_matches
			
		
		elif len(jaccard_lists_per_genome[idx]) > 0:
			#Case: not self; any intersections at all.
			jaccard_vec = np.array(jaccard_lists_per_genome[idx], dtype=np.float_)
			jaccard_mean = np.mean(jaccard_vec)
			jaccard_stddev = np.std(jaccard_vec)
		
			if 0.0058 < jaccard_mean < 0.843:
				aai_est = round(kaai_to_aai(jaccard_mean), 2)
			if jaccard_mean <= 0.0058:
				aai_est = "<30%"
			if jaccard_mean >= 0.843:
				aai_est = ">90%"
		
			#Pretty print
			jaccard_mean = round(jaccard_mean, 4)
			jaccard_stddev = round(jaccard_stddev, 4)
		
			#This is the number of observed protein matches with the query genome
			number_intersecting_prots = len(jaccard_vec)
		else:
			#Case: no intersections.
			#No point in throwing a warning or wasting the func call if there's nothing to do.
			jaccard_mean = "NA"
			jaccard_stddev = "NA"
			aai_est = "NA"
			#This is the number of observed protein matches with the query genome
			number_intersecting_prots = "NA"
			possible_matches = "NA"
		

		
		print(genome_index[query_id], genome_index[idx], jaccard_mean, jaccard_stddev, number_intersecting_prots, possible_matches, aai_est, sep = "\t", file = handle)
		
	handle.close()
	
	return temporary_output


# --- Parse viral kAAI when query == reference ---
# ------------------------------------------------------
def single_virus_kaai_parser(arguments):
	"""
	Calculates Jaccard distances on kmers from viral proteins
	
	Arguments:
		query_id {str} -- Id of the query genome
	
	Returns:
		[Path to output] -- Path to output file
	"""

	temporary_folder = arguments[0]
	query_id = arguments[1]
	skip_first_n = arguments[2]

	temporary_folder = Path(str(temporary_folder.name))
	temporary_file = Path(query_id).name + '.faai.temp'
	temporary_output = temporary_folder / temporary_file
	# Get query kmers
	proteins_query = query_kmer_dictionary[query_id][0]
	kmers_query = query_kmer_dictionary[query_id][1]

	# Start comparison with all genomes in the query dictionary
	with open(temporary_output, 'w') as out_file:
		for target_genome in list(query_kmer_dictionary.keys())[skip_first_n:]:
			# If self, 1.0 similarity
			if query_id == target_genome:
				out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
					1.0, proteins_query, proteins_query))
				continue

			jaccard_index = None
			proteins_reference = query_kmer_dictionary[target_genome][0]
			kmers_reference = query_kmer_dictionary[target_genome][1]
			# Calculate the Jaccard Index
			union = len(np.union1d(kmers_query, kmers_reference))
			intersection = len(kmers_query) + len(kmers_reference) - union
			jaccard_index = intersection/union
			out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
					jaccard_index, proteins_query, proteins_reference))
	return temporary_output
# ------------------------------------------------------

# --- Parse kAAI when query != reference ---
# ------------------------------------------------------
def double_kaai_parser(arguments):
	"""
	Calculates the Jaccard distances using single protein markers shared by two genomes
	
	Arguments:
		arguments {tuple} -- Tuple with the temporary folder, the query id and the index of said query_id
	
	Returns:
		[Path to output] -- Path to output file
	"""
	temporary_folder = arguments[0]
	query_id = arguments[1]
	
	temporary_folder = Path(str(temporary_folder.name))
	temporary_file = Path(query_id).name + '.faai.temp'
	temporary_output = temporary_folder / temporary_file
	
	query_scg_list = np.array(list(query_kmer_dictionary[query_id].keys()))

	with open(temporary_output, 'w') as out_file:
		for target_genome in list(reference_kmer_dictionary.keys()):
			# Get number and list of SCG detected in reference
			target_scg_list = np.array(list(reference_kmer_dictionary[target_genome].keys()))
			shorter_genome = min(len(query_scg_list), len(target_scg_list))
			#If self, 1.0 similarity.
			if query_id == target_genome:
					out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
					1.0, 0.0, len(query_scg_list), len(target_scg_list), 100))
					continue
			
			jaccard_similarities = []
			# Get shared proteins (scgs)
			final_scg_list = np.intersect1d(query_scg_list, target_scg_list)
			# Extract a list of kmers for each SCG in the list
			query_kmer_list = list(map(query_kmer_dictionary[query_id].get, final_scg_list))
			reference_kmer_list = list(map(reference_kmer_dictionary[target_genome].get, final_scg_list))
			# Calculate the jaccard index
			for accession in range(len(query_kmer_list)):
				union = len(np.union1d(query_kmer_list[accession], reference_kmer_list[accession]))
				intersection = len(query_kmer_list[accession]) + len(reference_kmer_list[accession]) - union
				jaccard_similarities.append(intersection / union)
			
			# Allow for numpy in-builts; they're a little faster.
			if len(jaccard_similarities) > 0:
				jaccard_similarities = np.array(jaccard_similarities, dtype=np.float_)
				try:
					mean = np.mean(jaccard_similarities)
					var = np.std(jaccard_similarities)
					
					if 0.0058 < mean < 0.843:
						aai_est = round(kaai_to_aai(mean), 2)
					if mean <= 0.0058:
						aai_est = "<30%"
					if mean >= 0.843:
						aai_est = ">90%"
						
					out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
							round(mean, 4), round(var, 4),
							len(jaccard_similarities), shorter_genome, aai_est))
				except:
					out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
							"NA", "NA", "NA", "NA", "NA"))
			else:
				out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
							"NA", "NA", "NA", "NA", "NA"))
	return temporary_output
# ------------------------------------------------------

# --- Parse kAAI when query == reference using faster, high mem algorithm---
# ------------------------------------------------------
def double_kaai_parser_fast(arguments):

	temporary_folder = arguments[0]
	query_id = arguments[1]
	#skip_first_n = arguments[2]
	
	temporary_folder = Path(str(temporary_folder.name))
	temporary_file = Path(genome_index[query_id]).name + '.faai.temp'
	temporary_output = temporary_folder / temporary_file
	
	#These should be globals ---
	#We always want all genomes in this case, so we ignore any skip_first_n
	valid_genome_list = np.sort(np.array(list(reference_kmer_dict.keys()), dtype = np.int32))
	
	#We could pre-make all of these outside this function and copy each time, but these won't take long to make anyway.
	jaccard_lists_per_genome = {}
	for gen in valid_genome_list:
		jaccard_lists_per_genome[gen] = []
	
	genome_intersection_counts = {}
	for gen in valid_genome_list:
		genome_intersection_counts[gen] = 0
		
	
	#We use this for keeping track of the number of times a genome appeared in protein tracking, which should be the number of intersecting protein families per genome.
	#number_of_prots_per_genome = genome_intersection_counts.copy()
	
	query_protein_ct = len(query_kmer_dict[query_id])
	
	#Go through the set of proteins and kmers for each protein in the query genome
	for protein in query_kmer_dict[query_id]:
		
		#We want to reinitialize the dict each time without recreating it from scratch
		non_ref = genome_intersection_counts.copy()
		
		#size of union = size of query protein kmers + size of target protein kmers - size of intersecting kmers
		query_union_length = len(query_kmer_dict[query_id][protein])
		
		#Inverted dict is already ref dict here, so no need to change...
		#This loop handles one full protein, which is one item in our jaccard index. This handles the protein for all genomes, though.
		for kmer in query_kmer_dict[query_id][protein]:
			#This gets us down to genomes that are actually found
			for genome in intersect1d_searchsorted(inverted_dict[kmer][protein], valid_genome_list):
				#They have to be sharing the protein at this point, so we can use the single protein key.
				non_ref[genome] += 1
					
			
		#We still have the protein referenced here, so we can reuse it.
		for genome in non_ref:
			
			try: 
				#Non_ref has genomes that have a count of zero, which means the protein we're looking at never showed up in it.
				#That means trying to access that protein for that genome will fail, so if it does, we want to continue but not modify the 0 size intersection we pre-set.
				#If it succeeds, we get it to do some work we want done anyway.
				target_union_length = len(reference_kmer_dict[genome][protein])
			except:
				continue
			
			#number_of_prots_per_genome[genome] += 1
			
			jaccard_lists_per_genome[genome].append((non_ref[genome] / (query_union_length + target_union_length - non_ref[genome])))
			
	#Open file here - access the name of the genome using the dict...
	
	handle = open(temporary_output, "w")
	
	#Write contents of jaccards per genome pairing.
	#We still need to get the number of kmers in the query vs the target for each protein, here...
	for idx in jaccard_lists_per_genome:
		
		#this is the number of possible protein matches with the query genome - between 0 and the shorter lsit of prots
		possible_matches = min(query_protein_ct, len(reference_kmer_dict[idx]))
		
		#For each genome, this is the list of the jaccard index of each protein [from 0 to possible matches]
		if query_id == idx:
			#Case: self
			jaccard_mean = 1.0
			jaccard_stddev = 0.0
			aai_est = 100
			#This is the number of observed protein matches with the query genome
			number_intersecting_prots = possible_matches
			
		
		elif len(jaccard_lists_per_genome[idx]) > 0:
			#Case: not self; any intersections at all.
			jaccard_vec = np.array(jaccard_lists_per_genome[idx], dtype=np.float_)
			jaccard_mean = np.mean(jaccard_vec)
			jaccard_stddev = np.std(jaccard_vec)
		
			if 0.0058 < jaccard_mean < 0.843:
				aai_est = round(kaai_to_aai(jaccard_mean), 2)
			if jaccard_mean <= 0.0058:
				aai_est = "<30%"
			if jaccard_mean >= 0.843:
				aai_est = ">90%"
				
			#Pretty print.
			jaccard_mean = round(jaccard_mean, 4)
			jaccard_stddev = round(jaccard_stddev, 4)
		
			#This is the number of observed protein matches with the query genome
			number_intersecting_prots = len(jaccard_vec)
		else:
			#Case: no intersections.
			#No point in throwing a warning or wasting the func call if there's nothing to do.
			jaccard_mean = "NA"
			jaccard_stddev = "NA"
			aai_est = "NA"
			#This is the number of observed protein matches with the query genome
			number_intersecting_prots = "NA"
			possible_matches = "NA"
		
		print(genome_index[query_id], genome_index[idx], jaccard_mean, jaccard_stddev, number_intersecting_prots, possible_matches, aai_est, sep = "\t", file = handle)
		
	handle.close()
	
	return temporary_output


# --- Parse viral kAAI when query != reference ---
# ------------------------------------------------------
def double_viral_kaai_parser(arguments):
	"""
	Calculates Jaccard distances on kmers from viral proteins
	
	Arguments:
		query_id {str} -- Id of the query genome
	
	Returns:
		[Path to output] -- Path to output file
	"""
	temporary_folder = arguments[0]
	query_id = arguments[1]

	temporary_folder = Path(str(temporary_folder.name))
	temporary_file = Path(query_id).name + '.faai.temp'
	temporary_output = temporary_folder / temporary_file
	# Get query kmers
	proteins_query = query_kmer_dictionary[query_id][0]
	kmers_query = query_kmer_dictionary[query_id][1]

	# Start comparison with all genomes in the query dictionary
	with open(temporary_output, 'w') as out_file:
		for target_genome in reference_kmer_dictionary.keys():
			# If self, 1.0 similarity
			if query_id == target_genome:
				out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
					1.0, proteins_query, proteins_query))
				continue

			jaccard_index = None
			proteins_reference = reference_kmer_dictionary[target_genome][0]
			kmers_reference = reference_kmer_dictionary[target_genome][1]
			# Calculate the Jaccard Index
			union = len(np.union1d(kmers_query, kmers_reference))
			intersection = len(kmers_query) + len(kmers_reference) - union
			jaccard_index = intersection/union
			out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
					jaccard_index, proteins_query, proteins_reference))
	return temporary_output
# ------------------------------------------------------

# --- Query == Reference initializer function ---
# ------------------------------------------------------
def single_dictionary_initializer(_dictionary):
	"""
	Make dictionary available for multiprocessing
	"""
	global query_kmer_dictionary
	query_kmer_dictionary = _dictionary
# ------------------------------------------------------

# --- Query != Reference initializer function ---
# ------------------------------------------------------
def two_dictionary_initializer(_query_dictionary, _reference_dictionary):
	"""
	Make dictionary available for multiprocessing
	"""
	global query_kmer_dictionary
	global reference_kmer_dictionary
	query_kmer_dictionary = _query_dictionary
	reference_kmer_dictionary = _reference_dictionary
# ------------------------------------------------------

# --- Merge kmer dictionaries ---
# ------------------------------------------------------
def merge_dicts(dictionaries):
	"""
	Given any number of dicts, shallow copy and merge into a new dict,
	precedence goes to key value pairs in latter dicts.
	"""
	result = {}
	for kmer_dictionary in dictionaries:
		result.update(kmer_dictionary)
	return result
# ------------------------------------------------------

# --- Estimate % AAI for 30 < similarity_score < 90 ---
# ------------------------------------------------------
def kaai_to_aai(kaai):
	
	# Transform the kAAI into estimated AAI values
	aai_hat = (-0.3087057 + 1.810741 * (np.exp(-(-0.2607023 * np.log(kaai))**(1/3.435))))*100
	return aai_hat
# ------------------------------------------------------

#Function for fast intersection. Assumes input arrays are unique and ascending sorted. Returns intersecting values.
def intersect1d_searchsorted(A,B):
	idx = np.searchsorted(B,A)
	idx[idx==len(B)] = 0
	return A[B[idx] == A]	


#--- Functions for the fast algorithm for kaai calculations ------
#Function for obtaining the string names of the genomes and proteins found in the kmer dict(s). Returns genome_index, protein_index dicts.
#No point in parallelizing. This is almost instant on thousands of genomes.
def obtain_genome_and_protein_family_indices(kmer_dicts):
	genome_index = {}
	protein_index = {}
	
	genome_counter = 0
	protein_counter = 0
	
	for dict in kmer_dicts:
		for genome in dict:
			try: 
				genome_index[genome]
			except:
				genome_index[genome] = genome_counter
				genome_counter += 1
			for protein in dict[genome]:
				try: 
					protein_index[protein]
				except:
					protein_index[protein] = protein_counter
					protein_counter += 1

	return genome_index, protein_index

#Function that swaps the string names of genomes and proteins for unique integer indices.	
def kmer_dict_index_swapper(kmer_dict, genome_index, protein_index):
	new_kmer_dict = {}

	for genome in kmer_dict:
		try: 
			new_kmer_dict[genome_index[genome]]
		except:
			new_kmer_dict[genome_index[genome]] = {}
		for protein in kmer_dict[genome]:
			try:
				#Check if the protein family's index is in the kmer dict.
				new_kmer_dict[genome_index[genome]][protein_index[protein]]
			except:
				#Assign the kmer string to the now index-only version of the kmer dict.
				new_kmer_dict[genome_index[genome]][protein_index[protein]] = kmer_dict[genome][protein]
				
	return new_kmer_dict

#Function for reversing the keys and values of the genome_index dict so that genome name can be accessed by integer index.
def gswap(g):
	reversed = {}
	for genome in g:
		reversed[g[genome]] = genome
		
	return(reversed)
	

#Get shared assets global.
def multi_initializer_single(_genome_idx, _kmer_dict, _inverted_dict):
	global genome_index
	global kmer_dict
	global inverted_dict
	genome_index = _genome_idx
	kmer_dict = _kmer_dict
	inverted_dict = _inverted_dict

#Get shared assets global - inverted is the inverted REFERENCE dict, not query. Query doesn't need inversion.
def multi_initializer_double(_genome_idx, _query_kmer_dict, _reference_kmer_dict, _inverted_dict):
	global genome_index
	global query_kmer_dict
	global reference_kmer_dict
	global inverted_dict
	genome_index = _genome_idx
	query_kmer_dict = _query_kmer_dict
	reference_kmer_dict = _reference_kmer_dict
	inverted_dict = _inverted_dict

	
'''
Function which inverts the kmer dictionary structure and returns a dict of results.

Structure of the kmer dict is a dict of dicts of arrays of form genome : protein_markers : [list of kmers per protein]
Structure of the tracker dictionary is a dict of unique kmers found across the set of all input proteins and arbitrary unique integer indices of form kmer : index

The return structure is a dict of dicts of arrays as the input kmer dict, but is of the form kmer : protein_marker : [genomes]
Each unique kmer is represented once, and keys to a dict of the set of protein markers in which that kmer is found.
Each protein marker keys to a list of genomes in which that protein marker is found.
'''
def unique_kmers_inverter(kmer_dict):

	invert_dict = {}
	
	#Easy prep work
	for kmer in global_kmer_index_dictionary:
		invert_dict[global_kmer_index_dictionary[kmer]] = {}
		
	#This can probably be parallelized.
	
	#Create the kmer:protein_fam:[genomes] dict; prepare the indexer
	for genome in kmer_dict:
		for pf in kmer_dict[genome]:
			#Because we pre-sort these as unique, each genome, pf should only go into the inverted dict zero or one times at each kmer.
			for kmer in kmer_dict[genome][pf]:
				try:
					invert_dict[kmer][pf].append(genome)
				except:
					invert_dict[kmer][pf] = [genome]
	
	#Switch to numpy array; sort the genomes for intersect1d_searchsorted speed increase
	for kmer in invert_dict:
		for pf in invert_dict[kmer]:
		#pure python list -> numpy array
			invert_dict[kmer][pf] = np.sort(np.array(invert_dict[kmer][pf], dtype=np.int32))
		
		
	return invert_dict
#----------------------------------------------------------------------------


#Function for acquiring input file paths and returning them as a list
#Used in parse_query_and_reference_inputs
def create_file_list(input_file_or_dir):
	#If the input is a directory, list the full path of each file in the dir.
	if os.path.isdir(input_file_or_dir):
		input_list = os.listdir(input_file_or_dir)
		for i in range(len(input_list)):
			input_list[i] = os.path.normpath(input_file_or_dir+"/"+input_list[i])
	#If the input is a file, the expectation is that the file has a path to the desired input files.
	else:
		input_list = []
		with open(input_file_or_dir, 'r') as input_fh:
			for line in input_fh:
				input_list.append(line.strip())
				
	return input_list


#Function that strips all extensions, no matter how many, from a path, then returns the last part of that path (after the final '/')
#Used in parse_query_and_reference_inputs
def get_input_base_name(filepath):
	
	#Get the last part of the path
	current_base_name = os.path.basename(filepath)
	previous_base_name = os.path.basename(filepath)
	
	#Remove extension if present. Leaves a naked string alone.
	current_base_name = os.path.splitext(current_base_name)[0]
	
	#Keep going until there's no more extensions to remove; never happens if there was no ext. in the first place
	while current_base_name != previous_base_name:
		previous_base_name = os.path.splitext(previous_base_name)[0]
		current_base_name = os.path.splitext(current_base_name)[0]
		
	return current_base_name
	

#Function for loading a database; 
#checks if the path is directly a database or if it's a file containing a list of paths
#Loads one or several DBs as a result
#used in check_and_load_database
def load_database(database_path_or_list_of_paths):
	kmer_dict = None
	kmer_dict_list = []
	
	#If the passed argument is a database, this will be set to true by the try block
	direct_reference = False
	
	#Attempt to directly open the path as a gzipped dict.
	try:		
		with gzip.open(database_path_or_list_of_paths, 'rb') as database_handle:
			temp_dict = pickle.load(database_handle)
		if isinstance(temp_dict,dict):
			kmer_dict_list.append(temp_dict)
		else:
			print("")
			exit("One of the database files appear to have the wrong format. Please provide a correctly formated databases.")
		kmer_dict = merge_dicts(kmer_dict_list)
		
		direct_reference = True
		print("done!")
	except:
		pass
	
	#or attempt to open the dict(s) from a list of paths in a file if the above fails
	if not direct_reference:
		with open(database_path_or_list_of_paths) as database_files:
			for db_location in database_files:
				if Path(db_location.strip()).is_file():
					with gzip.open(db_location.strip(), 'rb') as database_handle:
						temp_dict = pickle.load(database_handle)
						if isinstance(temp_dict,dict):
							kmer_dict_list.append(temp_dict)
						else:
							print("")
							exit("One of the database files appear to have the wrong format. Please provide a correctly formated databases.")
							
		kmer_dict = merge_dicts(kmer_dict_list)
		print("Database loaded!")
		
	return kmer_dict, kmer_dict_list


################################################################################
'''--- 2.0 Utility functions isloating flow of control for various steps of FastAAI ---'''

# --- Prepare tidy output directories for the user ---
# ------------------------------------------------------
def clean_outputs(output_directory, query_database, query_genomes, query_proteins, query_hmms, reference_database, reference_genomes, reference_proteins, reference_hmms):
	
	if not os.path.exists(output_directory):
		return None
		
	#handle users supplying '/' or no '/'
	output_directory = os.path.normpath(output_directory + "/FastAAI")
		
	#Base directory
	if not os.path.exists(output_directory):
		os.mkdir(output_directory)
		
	
	#Which directories need making?
	'''
	If any genomes are supplied, we need directories for predicted proteins, hmms, databases, and final results
	If any proteins are supplied, we need directories for hmms, databases, and final results
	If any hmms are supplied, we need directories for databases, and final results
	If only databases are supplied, we only need a directory for final results
	
	The logic cascades the directories needed from the lowest supplied input to the final results.
	'''
	proteins_and_onward = query_genomes != None or reference_genomes != None
	hmms_and_onward = (query_proteins != None or reference_proteins != None) or proteins_and_onward
	database_and_onward = (query_hmms != None or reference_hmms != None) or hmms_and_onward
	
	
	queries_started_before_hmms = (query_genomes != None) or (query_proteins != None and query_hmms is None)
	references_started_before_hmms = (reference_genomes != None) or (reference_proteins != None and reference_hmms is None)
	
	#We want to omit making the hmms folder if neither query nor reference needs it, which only happens if both query and reference start at or after HMM.
	started_from_hmms = ((query_hmms != None and query_proteins != None) or (reference_hmms != None and reference_proteins != None)) and not (queries_started_before_hmms or references_started_before_hmms)
	
	
	if proteins_and_onward:
		if not os.path.exists(os.path.normpath(output_directory + "/predicted_proteins")):
			os.mkdir(output_directory + "/predicted_proteins")

	if hmms_and_onward:
		if (not os.path.exists(os.path.normpath(output_directory + "/hmms")) and not started_from_hmms):
			os.mkdir(output_directory + "/hmms")
		if not os.path.exists(os.path.normpath(output_directory + "/filtered_hmms")):
			os.mkdir(output_directory + "/filtered_hmms")

	if database_and_onward:
		if not os.path.exists(os.path.normpath(output_directory + "/databases")):
			os.mkdir(output_directory + "/databases")
		
	#We always need this one
	if not os.path.exists(os.path.normpath(output_directory + "/final_result")):
		os.mkdir(output_directory + "/final_result")
	
	#Root path.
	return output_directory
# ------------------------------------------------------


#Function for determining if appropriate args have been passed to FastAAI; if options are valid, return the provided inputs as query and reference inputs and True if inputs are identical.
def check_inputs(query_genomes, query_proteins, query_hmms, query_database, virus, reference_genomes, reference_proteins, reference_hmms, reference_database):
	# Check if no query was provided
	if query_genomes == None and query_proteins == None and query_hmms == None and query_database == None:
		exit('Please prove a file with a list of queries, e.g., --qg, --qp, --qh, or --qd)')
	# Check query inputs
	query_input = None
	if query_hmms != None:
		if virus == True:
			exit("If you are comparing viruses, please start from the genome or protein files.")
		query_input = query_hmms
		if query_proteins != None:
			print("Starting from query hmmsearch results.")
			print("You also provided the list of protein files used for hmmsearch.")
		elif query_proteins == None:
			print("You chose to start from pre-computed hmmsearch results for your queries (--qh).")
			print("However, I also need the location of the query proteins used for hmmsearch.")
			exit("Please provide them with --qp.")
	elif query_proteins != None:
		query_input = query_proteins
		print("Starting from query proteins.")
	elif query_genomes != None:
		query_input = query_genomes
		print("Starting from query genomes.")
	elif query_database != None:
		query_input = query_database
		print("Starting from the pre-indexed query database.")
	# Check if no reference was provided
	if reference_genomes == None and reference_proteins == None and reference_hmms == None and reference_database == None:
		exit('Please prove a file with a list of references, e.g., --rg, --rp, --rh, or --rd)')
	# Check reference inputs
	reference_input = None
	if reference_hmms != None:
		if virus == True:
			exit("If you are comparing viruses, please start from the genome or protein files.")
		reference_input = reference_hmms
		if reference_proteins != None:
			print("Starting from reference hmmsearch results.")
			print("You also provided the list of protein files used for hmmsearch.")
		elif reference_proteins == None:
			print("You chose to start from pre-computed hmmsearch results for your references (--rh).")
			print("However, I also need the location of the query proteins used for hmmsearch.")
			exit("Please provide them with --rp.")
	elif reference_proteins != None:
		reference_input = reference_proteins
		print("Starting from reference proteins.")
	elif reference_genomes != None:
		reference_input = reference_genomes
		print("Starting from reference genomes.")
	elif reference_database != None:
		reference_input = reference_database
		print("Starting from the pre-indexed reference database.")
		
	return query_input, reference_input

	
#Function for checking presence, correctness of passed databases; load if present and correct and return results
def check_and_load_database(same_inputs, query_database, reference_database):
	#* Database Parsing is the same regardless of bacterial or viral genomes
	# If using pre-indexed databases, check if they are valid files.
	# ------------------------------------------------------
	# If any of the starting points is from database, then store the
	# kmer structures in the corresponding dictionaries.
	# Otherwise read the file list and get the filenames
	query_kmer_dict = None
	query_kmer_dict_list = []
	reference_kmer_dict = None
	reference_kmer_dict_list = []
		
	#direct load of a database file where the DB is passed as a path, not a list of paths.
	if query_database != None:
		print("Loading query database... ", end = '', flush = True)
		query_kmer_dict, query_kmer_dict_list = load_database(query_database)
	
	#The remainder of the code anticipates results of None for these args both if the refs provided are anything but a pre-indexed DB and if refs==query
	if reference_database != None and not same_inputs:
		print("Loading reference database... ", end = '', flush = True)
		reference_kmer_dict, reference_kmer_dict_list = load_database(reference_database)
			
	return query_kmer_dict, query_kmer_dict_list, reference_kmer_dict, reference_kmer_dict_list


#Function for collecting input sequences and producing file names for partial outputs
def parse_query_and_reference_inputs(output_root_directory, same_inputs, virus, query_database, query_input, query_genomes, query_proteins, query_hmms, reference_database, reference_input, reference_genomes, reference_proteins, reference_hmms):
	# Get files from the query and reference lists and then
	# create a dictionary with resulting filenames and a list with dictionary keys
	# The structure of the dictionary is:
	# original_query, proteins, hmms, filtered_hmms 
	# ------------------------------------------------------
	# First parse the query:
	query_list = []
	query_file_names = {}
	# For bacterial genomes
	if virus == False:
		if query_database != None:
			pass
		else:
			query_list = create_file_list(query_input)
			
			#Previously this was being done repeatedly for each genome.
			if query_hmms != None:
				query_protein_list = create_file_list(query_proteins)
			
			for index, query in enumerate(query_list):
				query_name = str(Path(query).name)
					
				#Query base name will be the file name without the last extension and without parent directories.
				query_base_name = get_input_base_name(query_name)
				
				predicted_protein_name = os.path.normpath(output_root_directory + "/predicted_proteins/" + query_base_name + '.faa')
				raw_hmm_name = os.path.normpath(output_root_directory + "/hmms/" + query_base_name + '.hmm')
				filtered_hmm_name = os.path.normpath(output_root_directory + "/filtered_hmms/" + query_base_name + '.hmm.filt')
					
				if query_hmms != None:
					query_file_names[query_name] = [None, query_protein_list[index], query, filtered_hmm_name]
				elif query_proteins != None:
					query_file_names[query_name] = [None, query, raw_hmm_name, filtered_hmm_name]
				elif query_genomes != None:
					query_file_names[query_name] = [query, predicted_protein_name, raw_hmm_name, filtered_hmm_name]
	# For viral genomes
	else:
		if query_database != None:
			pass
		else:
			query_list = create_file_list(query_input)
						
			for index, query in enumerate(query_list):
				query_name = str(Path(query).name)
				
				query_base_name = get_input_base_name(query_name)
				
				predicted_protein_name = os.path.normpath(output_root_directory + "/predicted_proteins/" + query_base_name + '.faa')
				
				if query_proteins != None:
					query_file_names[query_name] = [None, query]
				elif query_genomes != None:
					query_file_names[query_name] = [query, predicted_protein_name]

	# Then parse the references:
	reference_list = []
	reference_file_names = {}
	if same_inputs == True:
			pass
	else:
		# For bacterial genomes
		if virus == False:
			if reference_database != None:
				pass
			else:
				reference_list = create_file_list(reference_input)
				
				if reference_hmms != None:
					reference_list = create_file_list(reference_proteins)
				
				for index, reference in enumerate(reference_list):
					reference_name = str(Path(reference).name)
						
					reference_base_name = get_input_base_name(reference_name)
					
					predicted_protein_name = os.path.normpath(output_root_directory + "/predicted_proteins/" + reference_base_name + '.faa')
					raw_hmm_name = os.path.normpath(output_root_directory + "/hmms/" + reference_base_name + '.hmm')
					filtered_hmm_name = os.path.normpath(output_root_directory + "/filtered_hmms/" + reference_base_name + '.hmm.filt')
						
					if reference_hmms != None:					
						reference_file_names[reference_name] = [None, reference_protein_list[index], reference, filtered_hmm_name]
					elif reference_proteins != None:
						reference_file_names[reference_name] = [None, reference, raw_hmm_name, filtered_hmm_name]
					elif query_genomes != None:
						reference_file_names[reference_name] = [reference, predicted_protein_name, raw_hmm_name, filtered_hmm_name]
		# For viral genomes
		else:
			if reference_database != None:
				pass
			else:
				reference_list = create_file_list(reference_input)
				
				for index, reference in enumerate(reference_list):
					reference_name = str(Path(reference).name)
						
					reference_base_name = get_input_base_name(reference_name)
					predicted_protein_name = os.path.normpath(output_root_directory + "/predicted_proteins/" + reference_base_name + '.faa')
						
					if reference_proteins != None:
						reference_file_names[reference_name] = [None, reference]
					elif query_genomes != None:
						reference_file_names[reference_name] = [reference, predicted_protein_name]
	
	return 	query_list, query_file_names, reference_list, reference_file_names

	
#Function that creates dictionaries containing genomes, protein markers, and kmer lists for inputs. 
#Predicts proteins from input genomes with prodigal if needed. 
#Runs hmmer on proteins, predicted or input, if needed.
#Returns hmmer results if hmms are input.
def produce_kmer_dictionaries(_output_root_directory, threads, same_inputs, virus, query_kmer_dict, query_list, query_genomes, query_proteins, query_hmms, reference_kmer_dict, reference_list, reference_genomes, reference_proteins, reference_hmms, query_file_names, reference_file_names):
	# Pre-index queries
	
	#We need this repeatedly in many functs.
	global output_root_directory
	output_root_directory = _output_root_directory
	
	if query_kmer_dict == None:
		print("Beginning query processing.")
		# If using bacterial genomes
		if virus == False:
			if query_hmms != None:
				query_hmm_results = query_list
			elif query_proteins != None:
				query_protein_files = query_list
				print("Searching against HMM models... ", end = "", flush = True)
				try:
					pool = multiprocessing.Pool(threads)
					query_hmm_results = pool.map(run_hmmsearch, query_protein_files)
				finally:
					pool.close()
					pool.join()
					print("done!")
			elif query_genomes != None:
				print("Predicting proteins... ", end = "", flush = True)
				# Predict query proteins 
				try:
					pool = multiprocessing.Pool(threads)
					query_protein_files = pool.map(run_prodigal, query_list)
				finally:
					pool.close()
					pool.join()
					print("done!")
				print("Done!")
				print("Searching against HMM models... ", end = "", flush = True)
				# Run hmmsearch against proteins predicted
				try:
					pool = multiprocessing.Pool(threads)
					query_hmm_results = pool.map(run_hmmsearch, query_protein_files)
				finally:
					pool.close()
					pool.join()
					print("done!")
				print("Done!")
			print("Filtering query hmmsearch results... ", end = "", flush = True)
			# Filter query HMM search results
			try:
				pool = multiprocessing.Pool(threads)
				pool.map(hmm_filter, query_hmm_results)
			finally:
				pool.close()
				pool.join()
				print("done!")
			#
			#
			#
			#
			#
			#TODO these kmer_extracts should take the results of the filter funct, or get names earlier
			print("Extracting kmers from query proteins... ", end = "", flush = True)
			# Finding kmers for all queries
			query_information = []
			for name, values in query_file_names.items():
				query_information.append((name, values[1], values[3]))
			try:
				pool = multiprocessing.Pool(threads)
				kmer_results = pool.map(kmer_extract, query_information)
			finally:
				pool.close()
				pool.join()
				print("done!")
			query_kmer_dict = merge_dicts(kmer_results)
			del kmer_results
		# If using viral genomes
		else:
			if query_genomes != None:
				print("Predicting proteins... ", end = "", flush = True)
				# Predict query proteins 
				try:
					pool = multiprocessing.Pool(threads)
					query_protein_files = pool.map(run_prodigal_virus, query_list)
				finally:
					pool.close()
					pool.join()
					print("done!")
				print("Done!")
			elif query_proteins != None:
				query_protein_files = query_list
			print("Extracting kmers from query proteins... ", end = "", flush = True)
			query_information = []
			for name, values in query_file_names.items():
				query_information.append((name, values[1], 4))
			try:
				pool = multiprocessing.Pool(threads)
				kmer_results = pool.map(read_viral_kmers_from_file, query_information)
			finally:
				pool.close()
				pool.join()
				print("done!")
			query_kmer_dict = merge_dicts(kmer_results)
			del kmer_results
			
	# Pre-index references (if different from queries)
	if same_inputs == False and reference_kmer_dict == None:
		print("Beginning reference processing.")
		# If using bacterial genomes
		if virus == False:
			if reference_hmms != None:
				reference_hmm_results = reference_list
			elif reference_proteins != None:
				reference_protein_files = reference_list
				print("Searching against HMM models... ", end = "", flush = True)
				try:
					pool = multiprocessing.Pool(threads)
					reference_hmm_results = pool.map(run_hmmsearch, reference_protein_files)
				finally:
					pool.close()
					pool.join()
					print("done!")
			if reference_genomes != None:
				print("Predicting proteins... ", end = "", flush = True)
				# Predict reference proteins 
				try:
					pool = multiprocessing.Pool(threads)
					reference_protein_files = pool.map(run_prodigal, reference_list)
				finally:
					pool.close()
					pool.join()
					print("done!")
				print("Done!")
				print("Searching against HMM models... ", end = "", flush = True)
				# Run hmmsearch against proteins predicted
				try:
					pool = multiprocessing.Pool(threads)
					reference_hmm_results = pool.map(run_hmmsearch, reference_protein_files)
				finally:
					pool.close()
					pool.join()
					print("done!")
				print("Done!")
			print("Filtering reference hmmsearch results... ", end = "", flush = True)
			# Filter reference HMM search results
			try:
				pool = multiprocessing.Pool(threads)
				pool.map(hmm_filter, reference_hmm_results)
			finally:
				pool.close()
				pool.join()
				print("done!")
			print("Extracting kmers from reference proteins... ", end = "", flush = True) 
			# Finding kmers for all queries
			reference_information = []
			for name, values in reference_file_names.items():
				reference_information.append((name, values[1], values[3]))
			try:
				pool = multiprocessing.Pool(threads)
				kmer_results = pool.map(kmer_extract, reference_information)
			finally:
				pool.close()
				pool.join()
				print("done!")
			reference_kmer_dict = merge_dicts(kmer_results)
			del kmer_results
		# If using viral genomes
		else:
			if query_genomes != None:
				print("Predicting proteins... ", end = "", flush = True)
				# Predict query proteins 
				try:
					pool = multiprocessing.Pool(threads)
					query_protein_files = pool.map(run_prodigal, query_list)
				finally:
					pool.close()
					pool.join()
					print("done!")
			elif query_proteins != None:
				query_protein_files = query_list
			print("Extracting kmers from query proteins... ", end = "", flush = True)
			reference_information = []
			for name, values in reference_file_names.items():
				reference_information.append((name, values[1], 4))
			try:
				pool = multiprocessing.Pool(threads)
				kmer_results = pool.map(read_viral_kmers_from_file, reference_information)
			finally:
				pool.close()
				pool.join()
				print("done!")
			reference_kmer_dict = merge_dicts(kmer_results)
			del kmer_results
	
	#reference_kmer_dict will be empty if same inputs is true.
	return query_kmer_dict, reference_kmer_dict

	
#Function for writing compressed database outputs for quick start in subsequent runs
def save_databases(output_root_directory, same_inputs, query_database, query_input, reference_database, reference_input, query_kmer_dict, reference_kmer_dict):
	
	#If these are supplied as a directory, we want to make sure the name is kept regardless of the ending '/' that's likely with an autocomplete.
	if query_input.endswith("/"):
		query_input = query_input[:-1]
	if reference_input.endswith("/"):
		reference_input = reference_input[:-1]
	
	
	query_base_name = os.path.basename(os.path.splitext(query_input)[0])
	reference_base_name = os.path.basename(os.path.splitext(reference_input)[0])
	
	if same_inputs == True and query_database == None:
		print("Saving pre-indexed database... ", end = "", flush = True)
		query_database_name = output_root_directory + "/databases/" + query_base_name + '.db.gz'
		with gzip.open(query_database_name, 'wb') as database_handle:
			pickle.dump(query_kmer_dict, database_handle, protocol=4)
		print("done!")
	if same_inputs == False and query_database == None and reference_database == None:
		print("Saving pre-indexed databases... ", end = "", flush = True)
		query_database_name = output_root_directory + "/databases/" + query_base_name + '.db.gz'
		reference_database_name = output_root_directory + "/databases/" + reference_base_name + '.db.gz'
		with gzip.open(query_database_name, 'wb') as database_handle:
			pickle.dump(query_kmer_dict, database_handle, protocol=4)
		with gzip.open(reference_database_name, 'wb') as database_handle:
			pickle.dump(reference_kmer_dict, database_handle, protocol=4)
		print("done!")
	elif same_inputs == False and query_database == None:
		print("Saving pre-indexed query database... ", end = "", flush = True)
		query_database_name = output_root_directory + "/databases/" + query_base_name + '.db.gz'
		with gzip.open(query_database_name, 'wb') as database_handle:
			pickle.dump(query_kmer_dict, database_handle, protocol=4)
		print("done!")
	elif same_inputs == False and reference_database == None:
		print("Saving pre-indexed reference database... ", end = "", flush = True)
		reference_database_name = output_root_directory + "/databases/" + reference_base_name + '.db.gz'
		with gzip.open(reference_database_name, 'wb') as database_handle:
			pickle.dump(reference_kmer_dict, database_handle, protocol=4)
		print("done!")
				
	return None

#Prepare kmer dict for fast sort.
def fastsort_prep_single(kmer_dict):
	print("Using fast search algorithm.")
	
	print("Making genome and protein indexes... ", end = '', flush = True)
	genome_index, protein_index = obtain_genome_and_protein_family_indices([kmer_dict])
	print("done!")
	
	print("Converting genomes and proteins to indices... ", end = '', flush = True)
	kmer_dict = kmer_dict_index_swapper(kmer_dict, genome_index, protein_index)
	print("done!")
	
	print("Preparing kmers for fast search algorithm... ", end = '', flush = True)
	#Create a dict with kmer as the key, list of lists with list_1 = genomes, list_2 = protein families, where the L.O.L. has each kmer's groups on it
	inverted_dict = unique_kmers_inverter(kmer_dict)
	print("done!")
	
	genome_index = gswap(genome_index)
	
	return genome_index, kmer_dict, inverted_dict

	
def fastsort_prep_double(query_dict, reference_dict):
	print("Using fast search algorithm.")
	
	print("Making genome and protein indexes... ", end = '', flush = True)
	genome_index, protein_index = obtain_genome_and_protein_family_indices([query_dict, reference_dict])
	print("done!")
	
	print("Converting genomes and proteins to indices... ", end = '', flush = True)
	query_dict = kmer_dict_index_swapper(query_dict, genome_index, protein_index)
	reference_dict = kmer_dict_index_swapper(reference_dict, genome_index, protein_index)
	print("done!")
	
	print("Preparing reference kmers for fast search algorithm... ", end = '', flush = True)
	#Create a dict with kmer as the key, list of lists with list_1 = genomes, list_2 = protein families, where the L.O.L. has each kmer's groups on it
	#TODO parallize this. Split the kmer dict into [threads] chunks, have each process parse its chunk from a global dict, then merge the dicts.
	inverted_dict = unique_kmers_inverter(reference_dict)
	print("done!")
	
	genome_index = gswap(genome_index)
	
	return genome_index, query_dict, reference_dict, inverted_dict
	
	
#Function for calculating pairwise jaccard distances between the kmers of input genomes.
def calculate_jaccard_distances(threads, temporary_working_directory, same_inputs, virus, query_kmer_dict, reference_kmer_dict, fast):
	print("Preparing distance calculations.")

	if virus == False:
		if same_inputs == True:
			
			#Create global kmer index dictionary "global_kmer_index_dictionary" used in the pool functions.
			global_unique_kmers_parallel([query_kmer_dict], threads)
			
			#Separate the unique kmers and the index swapping.
			print("Preparing K-mers for searching... ", end = '', flush = True)
			query_kmer_dict = convert_kmers_to_indices(query_kmer_dict)
			print("done!")
			
			#The time gain from the faster alg isn't worth the prep time unless the set of genomes is pretty big.
			if fast:
				
				g, query_kmer_dict, inversion = fastsort_prep_single(query_kmer_dict)
				
				query_smart_args_tempdir = smart_arguments(query_kmer_dict, temporary_working_directory, single_dataset=True)
				
				print("Beginning FastAAI pairwise calculations... ", end = '', flush = True)
				
				try:
					pool = multiprocessing.Pool(threads, initializer = multi_initializer_single, initargs = (g, query_kmer_dict, inversion))
					Fraction_Results = pool.map(single_kaai_parser_fast, query_smart_args_tempdir)
				finally:
					pool.close()
					pool.join()
				
				print("done!")
				
			else:
			
				query_smart_args_tempdir = smart_arguments(query_kmer_dict, temporary_working_directory, single_dataset=True)
				
				print("Beginning FastAAI pairwise calculations... ", end = '', flush = True)
				
				try:
					pool = multiprocessing.Pool(threads, initializer = single_dictionary_initializer, initargs = (query_kmer_dict,))
					Fraction_Results = pool.map(single_kaai_parser, query_smart_args_tempdir)
				finally:
					pool.close()
					pool.join()
				
				print("done!")
			
		else:
			#Get the set of unique kmers and globalize the dict
			global_unique_kmers_parallel([query_kmer_dict, reference_kmer_dict], threads)
			
			#Switch the lists of kmers as comma-sep strings to integer numpy arrays
			print("Preparing query kmers for searching...")
			query_kmer_dict = convert_kmers_to_indices(query_kmer_dict)
			print("Queries processed.")
			#I have no idea why the ref dict doesn't want to play with the parallel.
			print("Preparing reference kmers for searching...")
			reference_kmer_dict = convert_kmers_to_indices(reference_kmer_dict)
			print("References processed.")
			
			#Skipped for testing. Normally length check is a good idea.
			if fast:
				
				g, query_kmer_dict, reference_kmer_dict, inversion = fastsort_prep_double(query_kmer_dict, reference_kmer_dict)
				
				#Get the smart arguments for parallelization.
				query_smart_args_tempdir = smart_arguments(query_kmer_dict, temporary_working_directory, single_dataset=False)
				_ref_smart_args_tempdir = smart_arguments(reference_kmer_dict, temporary_working_directory, single_dataset=False)
				
				print("Beginning FastAAI pairwise calculations... ", end = '', flush = True)
				
				try:
					pool = multiprocessing.Pool(threads, initializer = multi_initializer_double, initargs = (g, query_kmer_dict, reference_kmer_dict, inversion))
					Fraction_Results = pool.map(double_kaai_parser_fast, query_smart_args_tempdir)
				finally:
					pool.close()
					pool.join()
					
				print("done!")
				
			else:
				
				query_smart_args_tempdir = smart_arguments(query_kmer_dict, temporary_working_directory, single_dataset=False)
				_ref_smart_args_tempdir = smart_arguments(reference_kmer_dict, temporary_working_directory, single_dataset=False)
				
				print("Beginning FastAAI pairwise calculations... ", end = '', flush = True)
				
				try:
					pool = multiprocessing.Pool(threads, initializer = two_dictionary_initializer, initargs = (query_kmer_dict, reference_kmer_dict))
					Fraction_Results = pool.map(double_kaai_parser, query_smart_args_tempdir)
				finally:
					pool.close()
					pool.join()
					
				print("done!")
			
	else:
		if same_inputs == True:
			
			global_unique_viral_kmers_parallel([query_kmer_dict], threads)
			
			query_kmer_dict = convert_viral_kmers_to_indices(query_kmer_dict)
			
			query_smart_args_tempdir = smart_arguments_viral(query_kmer_dict, temporary_working_directory, single_dataset=True)
			
			print("Beginning FastAAI pairwise calculations... ", end = '', flush = True)
			try:
				pool = multiprocessing.Pool(threads, initializer = single_dictionary_initializer, initargs = (query_kmer_dict,))
				Fraction_Results = pool.map(single_virus_kaai_parser, query_smart_args_tempdir)
			finally:
				pool.close()
				pool.join()
				
			print("done!")
				
		else:
			
			global_unique_viral_kmers_parallel([query_kmer_dict, reference_kmer_dict], threads)
			
			query_kmer_dict = convert_viral_kmers_to_indices(query_kmer_dict)
			reference_kmer_dict = convert_viral_kmers_to_indices(reference_kmer_dict)
			
			query_smart_args_tempdir = smart_arguments_viral(query_kmer_dict, temporary_working_directory, single_dataset=False)
			_ref_smart_args_tempdir = smart_arguments_viral(reference_kmer_dict, temporary_working_directory, single_dataset=False)
			
			print("Beginning FastAAI pairwise calculations... ", end = '', flush = True)
			try:
				pool = multiprocessing.Pool(threads, initializer = two_dictionary_initializer, initargs = (query_kmer_dict, reference_kmer_dict))
				Fraction_Results = pool.map(double_viral_kaai_parser, query_smart_args_tempdir)
			finally:
				pool.close()
				pool.join()
				
			print("done!")

	return Fraction_Results

	
#Function for concatenating the jaccard distance outputs from calculate_jaccard_distances' per-genome temp outputs.
def merge_outputs(Fraction_Results, output, output_root_directory):
	print("Merging results... ", end = '', flush = True)
	with open(output_root_directory + "/final_result/" + output, 'w') as outfile:
		for file in Fraction_Results:
			with open(file) as Temp:
				shutil.copyfileobj(Temp, outfile)
			file.unlink()
	print("done!")
	return None

			
"""---3.0 Main Function---"""
#Parser modification for printing usage on no args.
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

#Separate options function for cleanliness in main.
def options():
	# Setup parser for arguments.
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''This script calculates the average amino acid identity using k-mers\n'''
						'''from single copy genes. It is a faster version of the regular AAI\n'''
						'''(Blast or Diamond) and the hAAI implemented in MiGA.\n\n'''
						
			'''Usage: ''' + argv[0] + ''' [OPTIONS] -q[g/p/h/d] -r[g/p/h/d]\n'''
			'''-q and -r specify query and reference, respectively.\n'''
			'''g/p/h/d specifies genomes, proteins, hmm models, and FastAAI pre-indexed databases, respectively.\n'''
			'''Example use: -qg [directory containing genomes in FASTA format] -rd [FastAAI database from a previous project]\n'''
			
			'''Optional Database Parameters: See ''' + argv[0] + ' -h\n'
			'''For more help, see''' + argv[0] + '--manual')
			
	mandatory_options = parser.add_argument_group('Mandatory i/o options. You must select an option for the queries and one for the references.')
	#---------------------------------------------------------------------------------------------------------------------------------------------
	mandatory_options.add_argument('--qg', dest='query_genomes', action='store', required=False,
									help='File with list of query genomes or a directory with only query genomes.')
	#---------------------------------------------------------------------------------------------------------------------------------------------								
	mandatory_options.add_argument('--qp', dest='query_proteins', action='store', required=False,
									help='File with list of query proteins or a directory with only query proteins.')
	#---------------------------------------------------------------------------------------------------------------------------------------------								
	mandatory_options.add_argument('--qh', dest='query_hmms', action='store', required=False,
									help=textwrap.dedent('''
									File with list of pre-computed query hmmsearch results.
									If you select this option you must also provide a file with 
									a list of protein files for the queries (with --qp).
									Directories containing these files are also acceptable.
									'''))
	#---------------------------------------------------------------------------------------------------------------------------------------------								
	mandatory_options.add_argument('--qd', dest='query_database', action='store', required=False,
									help='File with list of pre-indexed query databases or a path to a single pre-indexed query database.')
	#---------------------------------------------------------------------------------------------------------------------------------------------								
	mandatory_options.add_argument('--rg', dest='reference_genomes', action='store', required=False,
									help='File with list of reference genomes or a directory with only reference genomes.')
	#---------------------------------------------------------------------------------------------------------------------------------------------								
	mandatory_options.add_argument('--rp', dest='reference_proteins', action='store', required=False,
									help='File with list of reference proteins or a directory with only reference proteins.')
	#---------------------------------------------------------------------------------------------------------------------------------------------								
	mandatory_options.add_argument('--rh', dest='reference_hmms', action='store', required=False,
									help=textwrap.dedent('''
									File with list of pre-computed reference hmmsearch results.
									If you select this option you must also provide a file with 
									a list of protein files for the references (with --qp).
									Directories containing only these files are also acceptable.
									'''))
	#---------------------------------------------------------------------------------------------------------------------------------------------								
	mandatory_options.add_argument('--rd', dest='reference_database', action='store', required=False,
									help='File with list of pre-indexed reference databases or a path to a single pre-indexed reference database.')
	#---------------------------------------------------------------------------------------------------------------------------------------------								
	mandatory_options.add_argument('-o', '--output', dest='output', default = "kaai_comparisons.txt", action='store', required=False, help='Output file. By default kaai_comparisons.txt')
	#---------------------------------------------------------------------------------------------------------------------------------------------
	mandatory_options.add_argument('-d', '--directory', dest='directory', default = os.getcwd(), action='store', required=False, help='Place outputs in this directory. Default is the working directory.')
	#---------------------------------------------------------------------------------------------------------------------------------------------
	additional_input_options = parser.add_argument_group('Behavior modification options.')
	#---------------------------------------------------------------------------------------------------------------------------------------------										
	additional_input_options.add_argument('-i', '--index', dest='index_db', action='store_true', required=False, 
									help='Only index and store databases for later use. No AAI comparisons will be performed')
	#---------------------------------------------------------------------------------------------------------------------------------------------										
	misc_options = parser.add_argument_group('Miscellaneous options')
	#---------------------------------------------------------------------------------------------------------------------------------------------
	misc_options.add_argument('--virus', dest='virus', action='store_true', required=False,
									help='Toggle virus-virus comparisons. Use only with viral genomes or proteins.')
	#---------------------------------------------------------------------------------------------------------------------------------------------							
	misc_options.add_argument('-t', '--threads', dest='threads', action='store', default=1, type=int, required=False,
									help='Number of threads to use, by default 1')
	#---------------------------------------------------------------------------------------------------------------------------------------------							
	misc_options.add_argument('-f', '--fast', dest='fast', action='store_true', default=False, required=False,
									help='Use a faster method for calculating AAI. Consumes more memory. Only available for non-viral genomes.')
	#---------------------------------------------------------------------------------------------------------------------------------------------
	
	misc_options.add_argument('-man', '--manual', dest='manual', action='store_true', default=False, required=False,
									help='Print an extended usage manual concerning inputs and exit program.')
	
	args = parser.parse_args()
	
	return parser, args

#Returns a manual with extended usage descriptions.
def extended_usage_manual():
	#Newlines are arranged so that the output looks like the following text w/o newline chars
	manual = '''
		FastAAI Version 1.0
		
		FastAAI accepts a flexible range of inputs. 
		You can start both queries and references from:
		
		* Genomes in nucleotide FASTA format.
		* Proteins in animo acid FASTA format.
		* HMM models produced by a HMMER hmmsearch.
		* 1+ FastAAI pre-indexed database(s).
		
		FastAAI accepts these inputs in two ways:\n
		(1) A path to a directory containing only files of the input type.
		(2) A text file containing paths to files of the input type, one input file per line.
		
		If a directory is provided, all files in the directory will be used.
		
		Query inputs will be searched against reference inputs, and any
		combination of query-reference starting points is acceptable,
		i.e. (one of) query genomes/proteins/hmms/databases against
		(one of) reference genomes/proteins/hmms/databases will all work.
		FastAAI will handle any necessary intermediate steps, for any starting point.
		
		To specify an appropriate input, you simply combine the -q or -r flag
		with the file type you are selecting, e.g. -qg [queries] for query genomes,
		-rp [references] for reference proteins, and so on.
		
		To start from HMM files, both the HMMs and proteins are required, so
		you will need to provide both a -qh file and a -qp file to start from query hmms,
		and both a -rh file/directory and -rp file to start from reference hmms. In other cases,
		only one file per query/reference is required.
		
		If you specify the same file for both the query and the reference, FastAAI will
		create only a single database for that input, and will perform an all-vs-all
		comparison for genomes present in the input. This is the most efficient way to
		perform an all-vs-all comparison, and if used with the -i option (stop after making a FastAAI database)
		is the correct way to create a pre-indexed database for a single dataset.
		
		FastAAI will create an output directory named FastAAI, with child directories for
		proteins, hmms, filtered hmms, databases, and final results as needed. Using the
		-d option to specify an output directory will determine the location of the FastAAI directory,
		and giving a name to the results with the -o option will determine the name of the AAI comparisons
		file that will be found in [directory]/FastAAI/final_results/.
		
		Thank you for using FastAAI!
		'''
	
	return manual
	
#Main function calls options, prepares output directories, 
def main():

	#Get options block
	# ------------------------------------------------------
	help_message, args = options()
	
	#This is the case where a user supplies no options; prints usage
	if len(argv) == 1:
		exit(help_message.print_help())
		
	
	query_genomes = args.query_genomes
	reference_genomes = args.reference_genomes
	query_proteins = args.query_proteins
	reference_proteins = args.reference_proteins
	query_hmms = args.query_hmms
	reference_hmms = args.reference_hmms
	query_database = args.query_database
	reference_database = args.reference_database
	output = args.output
	output_directory = args.directory
	index_db = args.index_db
	threads = args.threads
	virus = args.virus
	fast_algorithm = args.fast
	
	print_manual = args.manual
	
	
	if print_manual:
		print(extended_usage_manual())
		exit()
		
	
	# ------------------------------------------------------

	
	# Check user input validity; select the appropriate non-empty inputs as query and reference; determine if query and reference files are identical.
	# ------------------------------------------------------
	query_input, reference_input = check_inputs(query_genomes, query_proteins, query_hmms, query_database, virus, reference_genomes, reference_proteins, reference_hmms, reference_database)
	
	#We don't need to initialize same inputs as T/F, nor do we need to check if a bool == T
	same_inputs = (query_input == reference_input)
	
	#correct messaging depending on the user's choices.
	if same_inputs and index_db:
		print('You specified the same query and reference files.')
		print('I will produce a single database as an output.')
	if same_inputs and not index_db:
		print('You specified the same query and reference files.')
		print('I will perform an all vs. all comparison :)')
	# ------------------------------------------------------
	
	
	# Create temporary working directory for outputs before final
	# ------------------------------------------------------
	temporary_working_directory = TemporaryDirectory()
	# ------------------------------------------------------
	
	#Wait until we do checks to make directories.
	output_root_directory = clean_outputs(output_directory, query_database, query_genomes, query_proteins, query_hmms, reference_database, reference_genomes, reference_proteins, reference_hmms)
	if output_root_directory == None:
		exit("The chosen output directory could not be found. This directory must already exist.")
	
	
	#Nothing meaningful happens before this point, so we shouldn't start the timer until here.
	time_format = "%d/%m/%Y %H:%M:%S"
	timer = datetime.datetime.now()
	printable_time = timer.strftime(time_format)
	print("FastAAI started on {}".format(printable_time))
		
	#Check to see if the user supplied a database of pre-indexed genomes and proteins from a previous kaai run; ensures database validity; loads database if present
	# ------------------------------------------------------
	query_kmer_dict, query_kmer_dict_list, reference_kmer_dict, reference_kmer_dict_list = check_and_load_database(same_inputs, query_database, reference_database)
	# ------------------------------------------------------
	
	
	#Collect sequence identifiers, file names as needed from inputs
	# ------------------------------------------------------
	query_list, query_file_names, reference_list, reference_file_names = parse_query_and_reference_inputs(output_root_directory, same_inputs, virus, query_database, query_input, query_genomes, query_proteins, query_hmms, reference_database, reference_input, reference_genomes, reference_proteins, reference_hmms)
	# ------------------------------------------------------

	
	# Pre-index and store databases
	# ------------------------------------------------------
	query_kmer_dict, reference_kmer_dict = produce_kmer_dictionaries(output_root_directory, threads, same_inputs, virus, query_kmer_dict, query_list, query_genomes, query_proteins, query_hmms, reference_kmer_dict, reference_list, reference_genomes, reference_proteins, reference_hmms, query_file_names, reference_file_names)
	# ------------------------------------------------------
	
	
	# Create or database(s) and compress it(them)
	# ------------------------------------------------------
	save_databases(output_root_directory, same_inputs, query_database, query_input, reference_database, reference_input, query_kmer_dict, reference_kmer_dict)
	# ------------------------------------------------------
	
	
	# Calculate Jaccard distances
	# ------------------------------------------------------
	if index_db == True:
		#Intentionally results in an early exit from FastAAI; this is just for producing the databases that can be later used for faster restarts
		print("Finished pre-indexing databases.")
		print("Next time you can run the program using only these files with --qd and(or) --rd.")
	else:
		#Calculates genome-genome jaccard similarity based on shared kmer sets
		Fraction_Results = calculate_jaccard_distances(threads, temporary_working_directory, same_inputs, virus, query_kmer_dict, reference_kmer_dict, fast_algorithm)
	# ------------------------------------------------------
	
	# Merge results into a single output
	# ------------------------------------------------------
		merge_outputs(Fraction_Results, output, output_root_directory)
	
	#Creating the DB and stopping there is also a correct finish, so we want this to print no matter how we got here.
	timer = datetime.datetime.now()
	printable_time = timer.strftime(time_format)
	print("FastAAI finished on {}".format(printable_time))


	
	
	
if __name__ == "__main__":
	main()
