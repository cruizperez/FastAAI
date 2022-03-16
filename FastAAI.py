#!/usr/bin/env python3

################################################################################
"""---0.0 Import Modules---"""
import subprocess 
import argparse 
import datetime 
import shutil
import textwrap
import multiprocessing
import pickle
import gzip
import tempfile
#Shouldn't play any role.
#from random import randint

#We could probably remove Path, too.
from pathlib import Path
#This as well
from functools import partial
import time
from collections import defaultdict
import sys
import os
from math import floor
import sqlite3
#numpy dependency
import numpy as np
import io
import random


#Takes a bytestring from the SQL database and converts it to a numpy array.
def convert_array(bytestring):
	return np.frombuffer(bytestring, dtype = np.int32)

def convert_float_array_16(bytestring):
	return np.frombuffer(bytestring, dtype = np.float16)

def convert_float_array_32(bytestring):
	return np.frombuffer(bytestring, dtype = np.float32)
	
def convert_float_array_64(bytestring):
	return np.frombuffer(bytestring, dtype = np.float64)

	
#Iterator for agnostic reader
class agnostic_reader_iterator:
	def __init__(self, reader):
		self.handle_ = reader.handle
		self.is_gz_ = reader.is_gz
		
	def __next__(self):
		if self.is_gz_:
			line = self.handle_.readline().decode()
		else:
			line = self.handle_.readline()
		
		#Ezpz EOF check
		if line:
			return line
		else:
			raise StopIteration

#File reader that doesn't care if you give it a gzipped file or not.
class agnostic_reader:
	def __init__(self, file):
		self.path = file
		
		with open(file, 'rb') as test_gz:
			#Gzip magic number
			is_gz = (test_gz.read(2) == b'\x1f\x8b')
		
		self.is_gz = is_gz
		
		if is_gz:
			self.handle = gzip.open(self.path)
		else:
			self.handle = open(self.path)
			
	def __iter__(self):
		return agnostic_reader_iterator(self)
		
	def close(self):
		self.handle.close()

#FastAAI database class. This is the final database
class fastaai_database:
	def __init__(self, path):
		#open SQL db and load in
		
		self.path = path
		self.exists = os.path.exists(path)
		
		self.child = None
		self.connection = None
		self.cursor = None
		
		self.child_connection = None
		self.child_cursor = None
		
		self.accessions = None
		#self.genomes = None
		
		#gak stands for 'genome_accession_kmer_counts'
		self.gak = None
		self.genome_index = None
		#Go from index to name
		self.reverse_genome_index = None
		self.protein_counts_by_genome = None
		
		#self.accession_set = None
		
		self.verbosity = False
		
	#Open an SQL connection
	def activate_connection(self, with_converter = True):
		# Converts np.array to TEXT when inserting
		##sqlite3.register_adapter(np.ndarray, adapt_array)

		#Converts byte string to numpy ndarray(int32) upon read from DB.
		if with_converter:
			sqlite3.register_converter("array", convert_array)
			self.connection = sqlite3.connect(self.path, detect_types=sqlite3.PARSE_DECLTYPES)
			
		else:
			#sqlite3.register_converter("array", convert_array)
			self.connection = sqlite3.connect(self.path)
			
		self.cursor = self.connection.cursor()
		self.exists = True
	
	#Close an SQL connection
	def close_connection(self):
		self.cursor.close()
		self.connection.close()
		#True cleanup - even a closed SQL connection obj cannot be passed to multiple processors, but a nonetype can.
		self.cursor = None
		self.connection = None
		
	def initialize_parent_database(self):
		if not self.exists:
			print("I need to be activated first!")
		else:
			#DB exists. Add metadata tables if needed.
			self.cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='genome_index' ''')
			if self.cursor.fetchone()[0]!=1 : 
				self.cursor.execute('''CREATE TABLE genome_index
				   (genome text, gen_id INTEGER PRIMARY KEY, protein_count INTEGER)''')
				self.connection.commit()
				
			self.cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='genome_acc_kmer_counts' ''')
			if self.cursor.fetchone()[0]!=1 : 
				self.cursor.execute('''CREATE TABLE genome_acc_kmer_counts
				   (genome INTEGER, accession INTEGER, count INTEGER)''')
				self.connection.commit()
				
	#Access an existing master database
	def activate_child_connection(self, child):
		#Don't try to connect unless it exists. This should never fail.
		if os.path.exists(child):
			self.child = child
			self.child_connection = sqlite3.connect(self.child, detect_types=sqlite3.PARSE_DECLTYPES)
			self.child_cursor = self.child_connection.cursor()
		else:
			print("Child database:", child, "not found!")
			
	#Close access to master DB
	def close_child_connection(self):
		if self.child_cursor is not None:
			self.child_cursor.close()
			self.child_connection.close()
			self.child_cursor = None
			self.child_connection = None
			self.child = None
				
	def add_child_to_parent(self, acc, child_db, remove = True, selected_kmers = None, genomes_too = False, just_genomes = False, update_gak = False):
		accession_index = generate_accessions_index()
		
		create_command = "CREATE TABLE IF NOT EXISTS " + acc + " (kmer INTEGER PRIMARY KEY, genomes array)"
		
		if not just_genomes:
			self.cursor.execute(create_command)
			self.connection.commit()
		
		if genomes_too or just_genomes:
			create_command = "CREATE TABLE IF NOT EXISTS " + acc + "_genomes (genome INTEGER PRIMARY KEY, kmers array)"
			self.cursor.execute(create_command)
			self.connection.commit()
		
		attach = "attach '"+child_db+"' as toMerge"
		
		if selected_kmers is not None:
			add = "INSERT OR REPLACE INTO " + acc + " SELECT * FROM toMerge." + acc + " WHERE kmer in ({kmers})".format(kmers = ','.join(['?']*len(selected_kmers)))
		else:
			add = "INSERT OR REPLACE INTO " + acc + " SELECT * FROM toMerge." + acc
		
		if genomes_too or just_genomes:
			add_genomes = "INSERT OR REPLACE INTO " + acc + "_genomes" + " SELECT * FROM toMerge." + acc+"_genomes"
			if update_gak:
				sql_acc_num = acc.replace("_", ".")
				sql_acc_num = accession_index[sql_acc_num]
				#Return num bytes, which is always 4*as many as there are entries, as the dtype is int32. See unique_kmers.
				gak_sql = 'INSERT OR REPLACE INTO genome_acc_kmer_counts SELECT genome, ' + str(sql_acc_num) + ', length(kmers)/4 FROM toMerge.' + acc + '_genomes'
		
		detach = "detach toMerge"
		
		self.cursor.execute(attach)
		self.connection.commit()
		
		if not just_genomes:
			if selected_kmers is not None:
				self.cursor.execute(add, selected_kmers)
			else:
				self.cursor.execute(add)
			
			self.connection.commit()	
		
		if genomes_too or just_genomes:
			self.cursor.execute(add_genomes)
			self.connection.commit()
			if update_gak:
				self.cursor.execute(gak_sql)
				self.connection.commit()
				
		self.cursor.execute(detach)
		self.connection.commit()
		
		if remove:
			os.remove(child_db)

	def add_genomes_first(self, accession, kmer_dict):
		kmer_lists = []
		for genome in kmer_dict:
			kmer_lists.append((genome, kmer_dict[genome].tobytes()))
			
		sql_friendly_accession = accession.replace(".", "_")
		
		#self.cursor.execute(" DROP TABLE IF EXISTS " + sql_friendly_accession + "_genomes")
		
		self.cursor.execute("CREATE TABLE IF NOT EXISTS " + sql_friendly_accession + "_genomes (genome INTEGER PRIMARY KEY, kmers array)")
		self.connection.commit()
		
		self.cursor.executemany("INSERT OR REPLACE INTO " + sql_friendly_accession + "_genomes VALUES (?, ?) ", kmer_lists)
			
		self.connection.commit()
		
		return sql_friendly_accession
		
			
	def load_genome_index(self):		
		self.genome_index = {}
		self.reverse_genome_index = {}
		self.protein_counts_by_genome = {}
		
		sql_command = ("SELECT genome, gen_id, protein_count FROM genome_index")
		
		#Break resist.
		gen = None
		id = None
		protein_count = None
		
		for result in self.cursor.execute(sql_command).fetchall():
			gen = result[0]
			id = result[1]
			protein_count = result[2]
			
			self.genome_index[gen] = id
			self.reverse_genome_index[id] = gen
			self.protein_counts_by_genome[id] = protein_count
				
		del gen
		del id
		del protein_count
		
	def load_accessions(self, permitted_genomes = None, permitted_accessions = None):		
		#self.protein_counts_by_genome = None
		
		self.gak = defaultdict(lambda: defaultdict())
		self.accessions = set()
				
				
		#It's possible to do both of these. Don't.
		if permitted_genomes is not None:
			sql_command = "SELECT * FROM genome_acc_kmer_counts WHERE genome IN ({genomes})".format(genomes=','.join(['?']*len(permitted_genomes)))
			#data type is very important to SQL
			sql_friendly = [int(permitted_genomes[i]) for i in range(0, len(permitted_genomes))]
			for result in self.cursor.execute(sql_command, sql_friendly).fetchall():
				genome, accession, kmer_ct = result[0], result[1], result[2]
				self.gak[genome][accession] = kmer_ct
				
		if permitted_accessions is not None:
			sql_command = "SELECT * FROM genome_acc_kmer_counts WHERE accession IN ({accessions})".format(accessions=','.join(['?']*len(permitted_accessions)))
			#data type is very important to SQL
			#sql_friendly = [int(permitted_accessions[i]) for i in range(0, len(permitted_genomes))]
			for result in self.cursor.execute(sql_command, permitted_accessions).fetchall():
				genome, accession, kmer_ct = result[0], result[1], result[2]
				self.gak[genome][accession] = kmer_ct
				
		#Normal case
		if permitted_accessions is None and permitted_genomes is None:
			sql_command = "SELECT * FROM genome_acc_kmer_counts"
			for result in self.cursor.execute(sql_command).fetchall():
				genome, accession, kmer_ct = result[0], result[1], result[2]
				self.gak[genome][accession] = kmer_ct

		#un-defaultdict
		self.gak = dict(self.gak)
		for genome in self.gak:
			self.gak[genome] = dict(self.gak[genome])
			self.accessions = self.accessions.union(self.gak[genome].keys())
			
		self.accessions = tuple(self.accessions)
		
	def just_accessions(self):
		converter = generate_accessions_index()
		acc_sql = "SELECT name FROM sqlite_master WHERE type='table'"
		tables = [item[0] for item in self.cursor.execute(acc_sql).fetchall()]
		
		genome_tables = []
		for table in tables:
			if table.endswith('_genomes'):
				genome_tables.append(table)
		
		for table in genome_tables:
			tables.pop(tables.index(table))

		tables.pop(tables.index('genome_acc_kmer_counts'))
		tables.pop(tables.index('genome_index'))
		
		#Back to indicies.
		tables = [converter[table.replace('_', '.')] for table in tables]
		
		self.accessions = tuple(tables)	
		
	def unload_genomes_and_accessions(self):
		self.gak = None
		self.genome_index = None
		#Go from index to name
		self.reverse_genome_index = None
		self.protein_counts_by_genome = None
		
#Child database class. This is only used during database builds and merges. Designed to take one single accession at a time and produce a correctly formatted table of kmers and accessions.
class child_database:
	def __init__(self, path, parent):
		#open SQL db and load in
		
		self.path = path
		self.exists = False
		
		self.parent = parent
		self.parent_exists = os.path.exists(parent)
		
		self.connection = None
		self.cursor = None
		
		self.parent_connection = None
		self.parent_cursor = None
		
		self.verbosity = False
		
	#Open an SQL connection
	def activate_child_connection(self):
		# Converts np.array to TEXT when inserting
		##sqlite3.register_adapter(np.ndarray, adapt_array)

		# Converts TEXT to np.array when selecting
		sqlite3.register_converter("array", convert_array)
		
		self.connection = sqlite3.connect(self.path, detect_types=sqlite3.PARSE_DECLTYPES)
		self.cursor = self.connection.cursor()
		self.exists = True
	
	#Close an SQL connection
	def close_child_connection(self):
		self.cursor.close()
		self.connection.close()
		#True cleanup - even a closed SQL connection obj cannot be passed to multiple processors, but a nonetype can.
		self.cursor = None
		self.connection = None
		
	def initialize_child_database(self):
		if not self.exists:
			print("I need to be activated first!")
		else:
			#DB exists. Add metadata tables.
			self.cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='genome_index' ''')
			if self.cursor.fetchone()[0]!=1 : 
				self.cursor.execute('''CREATE TABLE genome_index
				   (genome text, gen_id integer, protein_count integer)''')
				self.connection.commit()
				
			self.cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='genome_acc_kmer_counts' ''')
			if self.cursor.fetchone()[0]!=1 : 
				self.cursor.execute('''CREATE TABLE genome_acc_kmer_counts
				   (genome integer, accession integer, count integer)''')
				self.connection.commit()
				
	
	#Access an existing master database
	def activate_parent_connection(self):
		if os.path.exists(self.parent):
			self.parent_exists = True
			#sqlite3.register_adapter(np.ndarray, adapt_array)
			# Converts TEXT to np.array when selecting
			sqlite3.register_converter("array", convert_array)
			self.parent_connection = sqlite3.connect(self.parent, detect_types=sqlite3.PARSE_DECLTYPES)
			self.parent_cursor = self.parent_connection.cursor()
	
	#Close access to master DB
	def close_parent_connection(self):
		if self.parent_cursor is not None:
			self.parent_cursor.close()
			self.parent_connection.close()
			self.parent_cursor = None
			self.parent_connection = None
	
	def add_genomes_first(self, accession, kmer_lists):
	
		#kmer_lists = []
		#Shoot... gotta pass the args
		
		#for file in prepared_files:
		#	if accession in file.best_hits_kmers:
		#		kmer_lists.append((genome_index[file.basename], file.best_hits_kmers[accession].tobytes()))
	
		sql_friendly_accession = accession.replace(".", "_")
		
		self.cursor.execute(" DROP TABLE IF EXISTS " + sql_friendly_accession + "_genomes")
		
		self.cursor.execute("CREATE TABLE " + sql_friendly_accession + "_genomes (genome INTEGER PRIMARY KEY, kmers array)")
		self.connection.commit()
		
		self.cursor.executemany(" INSERT INTO " + sql_friendly_accession + "_genomes VALUES (?, ?) ", kmer_lists)
			
		self.connection.commit()
		
		return sql_friendly_accession
		
	
	def add_accession(self, accession, insert_kmers):
		sql_friendly_accession = accession.replace(".", "_")
		
		if self.parent_exists:
			parent_kmers = {}
			#Check to see if this acc. is already in parent DB
			table_exists = (self.parent_cursor.execute(" SELECT count(name) FROM sqlite_master WHERE type='table' AND name=(?)", (sql_friendly_accession,)).fetchone()[0] == 1)
			#If the accession is in the parent DB
			if table_exists:
				#Select the records where the kmers are in the new kmers to be added - we don't have to modify the ones that aren't.
				search_command = "SELECT * FROM "+ sql_friendly_accession + " WHERE kmer IN ({kmers})".format(kmers=','.join(['?']*len(insert_kmers)))
				
				#Convert the kmers in the current insert list to the correct type for sql to match them
				selection = tuple([int(key) for key in insert_kmers.keys()])
				
				for item in self.parent_cursor.execute(search_command, selection).fetchall():
					#Get the kmer for this parent
					k = item[0]
					#If the record would be modified in the parent, combine the to-add (which will replace the row) with the existing data. Otw. the record is unaffected and we can ignore it.
					if k in insert_kmers:
						insert_kmers[k] = np.union1d(insert_kmers[k], item[1])
			
			
			#Free up the space.
			del parent_kmers
		
		formatted_kmers = []
		
		#Translate the ndarray into its constituent byte data
		for kmer in insert_kmers:
			formatted_kmers.append((int(kmer), insert_kmers[kmer].tobytes(), ))
			
		del insert_kmers
			
		#Remove the child if it exists - it shouldn't ever exist because these child DBs should be deleted upon being added to the parent, but might if a run was stopped halfway.
		self.cursor.execute(" DROP TABLE IF EXISTS " + sql_friendly_accession)
		
		self.cursor.execute("CREATE TABLE " + sql_friendly_accession + " (kmer INTEGER PRIMARY KEY, genomes array)")
		self.connection.commit()
		
		self.cursor.executemany(" INSERT INTO " + sql_friendly_accession + " VALUES (?, ?) ", formatted_kmers)
			
		self.connection.commit()
			
		del formatted_kmers
		
		return sql_friendly_accession
	
	
#Holds partial results for calculating AAI.
class calculation_database:
	def __init__(self, path, precision):
		#open SQL db and load in
		
		self.path = path
		self.exists = False

		self.connection = None
		self.cursor = None
		
		self.genomes = None

		self.verbosity = False
		
		self.precision = precision
		
	#Open an SQL connection
	def activate_connection(self):
		# Converts np.array to TEXT when inserting
		##sqlite3.register_adapter(np.ndarray, adapt_array)

		# Converts TEXT to np.array when selecting
		if self.precision == "low":
			sqlite3.register_converter("array", convert_float_array_16)
		if self.precision == "med":
			sqlite3.register_converter("array", convert_float_array_32)
		if self.precision == "high":
			sqlite3.register_converter("array", convert_float_array_64)
		
		self.connection = sqlite3.connect(self.path, detect_types=sqlite3.PARSE_DECLTYPES)
		self.cursor = self.connection.cursor()
		self.exists = True
	
	#Close an SQL connection
	def close_connection(self):
		self.cursor.close()
		self.connection.close()
		#True cleanup - even a closed SQL connection obj cannot be passed to multiple processors, but a nonetype can.
		self.cursor = None
		self.connection = None
		
	def initialize_database(self):
		if not self.exists:
			print("I need to be activated first!")
		else:
			#DB exists. Add metadata tables.
			self.cursor.execute("DROP TABLE IF EXISTS jaccards")
			self.connection.commit()
			self.cursor.execute("CREATE TABLE jaccards (genome INTEGER PRIMARY KEY, jaccards array)")
			self.connection.commit()
	
'''
Class for handling all of the raw genome/protein/protein+HMM file inputs when building a database.

Takes a file or files and processes them from genome -> protein, protein -> hmm, prot+HMM -> kmerized protein best hits as numpy int arrays according to the kmer_index

'''
class input_file:
	def __init__(self, input_path, output, verbosity):
		#starting path for the file; irrelevant for protein and hmm, but otherwise useful for keeping track.
		self.path = input_path
		#Output directory starts with this
		self.output = os.path.normpath(os.path.basename(output) + "/")
		#For printing file updates, this is the input name
		self.name = os.path.basename(input_path)
		#original name is the key used for the genomes index later on.
		self.original_name = os.path.basename(input_path)
		#This is the name that can be used for building files with new extensions.
		if input_path.endswith(".gz"):
			#Remove .gz first to make names consistent.
			self.basename = os.path.splitext(os.path.basename(input_path[:-3]))[0]
		else:
			self.basename = os.path.splitext(os.path.basename(input_path))[0]
		#'genome' or 'protein' or 'protein and HMM' 
		self.status = None
		#These will keep track of paths for each stage of file for us.
		self.genome = None
		self.protein = None
		self.hmm = None
		
		self.best_hits = None
		self.best_hits_kmers = None
		
		self.protein_count = 0
		self.protein_kmer_count = {}
		
		self.trans_table = None
		self.start_time = None
		self.end_time = None
		self.err_log = ""
		#doesn't get updated otw.
		self.initial_state = "protein+HMM"
		
		self.verbose = verbosity
		
		#r_scripts_loc = os.path.dirname(sys.modules['metapop'].__file__) + "/metapop_r/"
		#"00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"
		self.hmm_path = None
		try:
			#Try to locate the data bundled as it would be with a pip/conda install.
			script_path = os.path.dirname(sys.modules['fastAAI_HMM_models'].__file__)
			hmm_complete_model = script_path + '/00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm'
			self.hmm_path = str(hmm_complete_model)
			#Check that the file exists or fail to the except.
			fh = open(self.hmm_path)
			fh.close()
		except:
			#Look in the same dir as the script; old method/MiGA friendly
			script_path = Path(__file__)
			script_dir = script_path.parent
			hmm_complete_model = script_dir / "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"
			self.hmm_path = str(hmm_complete_model)
	
	#Functions for externally setting status and file paths of particular types
	def set_genome(self, path):
		self.status = 'genome'
		self.genome = path
	
	def set_protein(self, path):
		self.status = 'protein'
		self.protein = path

	def set_hmm(self, path):
		if self.protein is None:
			print("Warning! I don't have a protein yet, so this HMM will be useless to me until I do!")
		self.status = 'protein and hmm'
		self.hmm = path

	#Runs prodigal, compares translation tables and stores faa files
	def genome_to_protein(self):
		if self.genome is None:
			print(self.name, "wasn't a declared as a genome! I can't make this into a protein!")
		else:
			folder = Path(self.output + "/predicted_proteins")
			protein_output = folder / (self.basename + '.faa')
			output_11 = folder / (self.basename + '.faa.11')
			output_4 = folder / (self.basename + '.faa.4')
			temp_output = folder / (self.basename + '.temp')
			
			intermediate = folder / (self.basename + '_genome_intermediate.fasta')
			
			#total_bases = 0
			
			genome_parser = agnostic_reader(self.genome)
			
			if genome_parser.is_gz:
				#File was a gzip; decompress it to an intermediate file and then run prodigal; delete after
				#print("unzipping input...")
				midpoint = open(intermediate, "w")
				#Count input bases and write an unzipped file for prodigal's sake.
				for line in genome_parser:
					#if not line.startswith(">"):
					#	total_bases += len(line.strip())
					midpoint.write(line)
					
				midpoint.close()
				
			else:
				#File is already unzipped, just point to it
				intermediate = self.genome
				#Count input bases
				#for line in genome_parser:
				#	if not line.startswith(">"):
				#		total_bases += len(line.strip())
			
			genome_parser.close()
			'''
			A chunk of code originally indended to match GTDBtk's table selection criteria.
			if total_bases > 100000:
				#training mode
				subprocess.call(["prodigal", "-i", str(intermediate), "-a", str(output_11), "-q", "-o", str(temp_output)])				
				subprocess.call(["prodigal", "-i", str(intermediate), "-a", str(output_4), "-g", "4", "-q", "-o", str(temp_output)])
			else:
				#Metagenome mode for very short genomes.
				subprocess.call(["prodigal", "-i", str(intermediate), "-p", "meta", "-a", str(output_11), "-q", "-o", str(temp_output)])				
				subprocess.call(["prodigal", "-i", str(intermediate), "-p", "meta", "-a", str(output_4), "-g", "4", "-q", "-o", str(temp_output)])
			'''
			
			subprocess.call(["prodigal", "-i", str(intermediate), "-a", str(output_11), "-q", "-o", str(temp_output)])				
			subprocess.call(["prodigal", "-i", str(intermediate), "-a", str(output_4), "-g", "4", "-q", "-o", str(temp_output)])
			
			#We can get rid of the temp file immediately, we won't be using it
			temp_output.unlink()
			if genome_parser.is_gz:
				#If the file was copied, delete. Otw. this would delete the input and we don't want that.
				intermediate.unlink()

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
			
			#Select the winning translation table and remove the other. Open the winner.
			if (length_4 / length_11) >= 1.1:
				output_11.unlink()
				self.trans_table = "4"
				chosen_protein = open(output_4, 'r')
				table_11 = False
			else:
				output_4.unlink()
				self.trans_table = "11"
				chosen_protein = open(output_11, 'r')
				table_11 = True
				
			destination = open(protein_output, "w")
			
			#Clean the winning output.
			for line in chosen_protein:
				if line.startswith(">"):
					destination.write("{}".format(line))
				else:
					line = line.replace('*', '')
					destination.write("{}".format(line))
			
			destination.close()
			chosen_protein.close()
			
			# Remove the winning intermediate file, since we have the cleaned output
			if table_11:
				output_11.unlink()
			else:
				output_4.unlink()
			
			self.set_protein(str(protein_output))

	#run hmmsearch on a protein	
	def protein_to_hmm(self):
		if self.protein is None:
			print(self.name, "wasn't a declared as a protein! I can't make this into an HMM!")
		else:
		
			folder = Path(self.output + "/hmms")
						
			hmm_output = folder / (self.basename + '.hmm')
			temp_output = folder / (self.basename + '.temp')
			
			intermediate = folder / (self.basename + '_protein_intermediate.faa')
			
			current_protein = ""
			current_seq = ""
			
			protein_parser = agnostic_reader(self.protein)
			
			#File was a gzip; decompress it to an intermediate file and then run prodigal; delete after
			#Keeps track of \n chars in the protein sequences.
			line_ct = 0
			midpoint = open(intermediate, "w")
			
			for line in protein_parser:
				if line.startswith(">"):
					if len(current_seq) > 0:
						if len(current_seq) < 100000:
							midpoint.write(current_protein)
							midpoint.write(current_seq)
						else:
							self.err_log += "Protein " + current_protein.strip().split()[0][1:] + " was observed to have >100K amino acids ( " + str(len(current_seq) - line_ct) + " AA found ). It was skipped. "
							#print("Protein", current_protein.strip()[1:], "was observed to have >100K amino acids (", len(current_seq) - line_ct, "AA found ).", file = sys.stderr)
							#print("HMMER cannot handle sequences that long, and the protein is almost certainly erroneous, anyway.", file = sys.stderr)
							#print("The protein will be skipped, and FastAAI will continue without it.", file = sys.stderr)
							
					current_protein = line
					current_seq = ""
					line_ct = 0
				else:
					line_ct += 1
					current_seq += line
			
			protein_parser.close()
			
			#Finally, last prot
			if len(current_seq) > 0:
				if len(current_seq) < 100000:
					midpoint.write(current_protein)
					midpoint.write(current_seq)
				else:
					self.err_log += "Protein " + current_protein.strip().split()[0][1:] + " was observed to have >100K amino acids ( " + str(len(current_seq) - line_ct) + " AA found ). It was skipped. "
					#print("Protein", current_protein.strip()[1:], "was observed to have >100K amino acids (", len(current_seq) - line_ct, "AA found ).", file = sys.stderr)
					#print("HMMER cannot handle sequences that long, and the protein is almost certainly erroneous, anyway.", file = sys.stderr)
					#print("The protein will be skipped, and FastAAI will continue without it.", file = sys.stderr)
			
			midpoint.close()
			
			#Should locate the DBs regardless of path.
			script_path = Path(__file__)
			script_dir = script_path.parent
			hmm_complete_model = script_dir / "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"
			
			subprocess.call(["hmmsearch", "--tblout", str(hmm_output), "-o", str(temp_output), "--cut_tc", "--cpu", "1",
							str(hmm_complete_model), str(intermediate)])
							
			temp_output.unlink()
			intermediate.unlink()
			
			self.set_hmm(str(hmm_output))
	
	def prot_and_hmm_to_besthits(self):
		prots = []
		accs = []
		scores = []
		f = agnostic_reader(self.hmm)
		for line in f:
			if line.startswith("#"):
				continue
			else:
				segs = line.strip().split()
				prots.append(segs[0])
				accs.append(segs[3])
				scores.append(segs[8])
			
		f.close()
		
		hmm_file = np.transpose(np.array([prots, accs, scores]))
		
		#hmm_file = np.loadtxt(hmm_file_name, comments = '#', usecols = (0, 3, 8), dtype=(str))
		#Sort the hmm file based on the score column in descending order.
		hmm_file = hmm_file[hmm_file[:,2].astype(float).argsort()[::-1]]
		
		#Identify the first row where each gene name appears, after sorting by score; 
		#in effect, return the highest scoring assignment per gene name
		#Sort the indices of the result to match the score-sorted table instead of alphabetical order of gene names
		hmm_file = hmm_file[np.sort(np.unique(hmm_file[:,0], return_index = True)[1])]
		
		#Filter the file again for the unique ACCESSION names, since we're only allowed one gene per accession, I guess?
		#Don't sort the indices, we don't care about the scores anymore.
		hmm_file = hmm_file[np.unique(hmm_file[:,1], return_index = True)[1]]
		
		self.best_hits = dict(zip(hmm_file[:,0], hmm_file[:,1]))
		
		self.best_hits_kmers = {}
		current_seq = ""
		current_prot = ""
		is_besthit = False
		
		prot = agnostic_reader(self.protein)
		
		for line in prot:
			
			if line.startswith(">"):
				if len(current_seq) > 0:
					kmer_set = unique_kmers(current_seq, 4)
					self.protein_kmer_count[current_prot] = kmer_set.shape[0]
					self.protein_count += 1
					self.best_hits_kmers[current_prot] = kmer_set
				#Select the best hit accession for this protein and just record that. We do not care about the names of the proteins.
				current_prot = line[1:].strip().split(" ")[0]
				if current_prot in self.best_hits:
					current_prot = self.best_hits[current_prot]
					is_besthit = True
				else:
					is_besthit = False
				current_seq = ""
			else:
				if is_besthit:
					current_seq += line.strip()
				
		prot.close()
	
		#Final iter. doesn't happen otw.
		if current_prot in self.best_hits:
			kmer_set = unique_kmers(current_seq, 4)
			#kmer_set = [kmer_index[k] for k in  kmer_set]
			self.protein_kmer_count[current_prot] = kmer_set.shape[0]
			self.protein_count += 1
			self.best_hits_kmers[current_prot] = kmer_set
			
		self.status = "finished preprocessing"
	
	def preprocess(self):
		#There's no advancement stage for protein and HMM
		if self.status == 'genome':
			start_time = curtime()
			#report = True
			if self.start_time is None:
				self.start_time = start_time
			
			if self.initial_state == "protein+HMM":
				self.initial_state = "genome"
			
			self.genome_to_protein()
			
			
		if self.status == 'protein':
			start_time = curtime()
			#report = True
			if self.start_time is None:
				self.start_time = start_time
				
			if self.initial_state == "protein+HMM":
				self.initial_state = "protein"
			
			self.protein_to_hmm()
			
		if self.status == 'protein and hmm':
			start_time = curtime()
			
			if self.start_time is None:
				self.start_time = start_time
				
			self.prot_and_hmm_to_besthits()
			
		#Add an end time if either genome -> protein -> HMM or protein -> HMM happened.
		if self.start_time is not None:
			end_time = curtime()
			self.end_time = end_time
		else:
			#Start was protein+HMM. There was no runtime, and intitial state is p+hmm
			#self.initial_state = "protein+HMM"
			self.start_time = "N/A"
			self.end_time = "N/A"
			
		#Protein not generated on this run.
		if self.trans_table is None:
			self.trans_table = "unknown"
	
	'''
	Viral functions
	'''
	#No translation table comparison for viruses. Slightly reduced logic.
	def viral_genome_to_protein(self):
		if self.genome is None:
			print(self.name, "wasn't a declared as a genome! I can't make this into a protein!")
		else:
			folder = Path(self.output + "/predicted_proteins")
			intermediate_protein_output = folder / (self.basename + '.intermediate.faa')
			final_protein_output = folder / (self.basename + '.faa')
			temp_output = folder / (self.basename + '.temp')
	
			subprocess.call(["prodigal", "-i", str(self.genome), "-a", str(intermediate_protein_output), "-p", "meta", "-q", "-o", str(temp_output)])
			
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
			
			self.protein = str(protein_output)
			self.status = 'protein'
	

'''
Preprocessing functions
				
Read directories, advance files to hmms as needed.				
'''
#Toy function for passing to a pool
def do_advance(input_file_object):
	input_file_object.preprocess()
	return input_file_object

def initialize_preproc(index):
	global kmer_index
	kmer_index = index
	
#Function which takes an input list
def advance_inputs(genomes = None, proteins = None, hmms = None, genomes_file = None, proteins_file = None, hmms_file = None, output = "FastAAI", threads = 1, verbose = False, db_name = ""):
	inputs = []
	
	hmm_broke = False
	
	if genomes_file is not None:
		fh = agnostic_reader(genomes_file)
		
		for line in fh:
			clean = line.strip()
			if not os.path.exists(clean):
				print("I can't find file", clean, "Are you sure this file exists and can be found from your current directory using the path you supplied in the input file?")
			else:
				current_file = input_file(clean, output, verbose)
				current_file.set_genome(clean)
				inputs.append(current_file)
				del current_file
			
		fh.close()
		
	if proteins_file is not None:
		fh = agnostic_reader(proteins_file)
		
		for line in fh:
			#GOTOGOTO
			print(line)
		
			clean = line.strip()
			if not os.path.exists(clean):
				print("I can't find file", clean, "Are you sure this file exists and can be found from your current directory using the path you supplied in the input file?")
			else:
				current_file = input_file(clean, output, verbose)
				current_file.set_protein(clean)
				inputs.append(current_file)
				del current_file
			
		fh.close()
		
	if hmms_file is not None:	
		fh = agnostic_reader(hmms_file)
		
		hmm_pairs = []
		
		for line in fh:
			clean = line.strip()
			if not os.path.exists(clean):
				print("I can't find file", clean, "Are you sure this file exists and can be found from your current directory using the path you supplied in the input file?")
			else:
				hmm_pairs.append(clean)
				
		fh.close()
		
		if len(hmm_pairs) != len(inputs):
			print("Protein and HMM file counts differ! There must be one HMM per protein, generated from its paired protein! These pairs must be in the same order in your input file!")
			hmm_broke = True
		else:
			for h, i in zip(hmm_pairs, inputs):
				i.set_hmm(h)
		
	if genomes is not None:
		set = os.listdir(genomes)
		#Sort is used to ensure lexicographic ordering.
		set.sort()
		set = [os.path.normpath(genomes + "/" + file) for file in set]
		
		for file in set:
			if not os.path.exists(file):
				print("I can't find", file, "Are you sure this file exists in the directory you supplied?")
			else:
				current_file = input_file(file, output, verbose)
				current_file.set_genome(file)
				inputs.append(current_file)
				del current_file
		
	if proteins is not None:
		set = os.listdir(proteins)
		set.sort()
		set = [os.path.normpath(proteins + "/" + file) for file in set]
		
		for file in set:
			if not os.path.exists(file):
				print("I can't find", file, "Are you sure this file exists in the directory you supplied?")
			else:
				current_file = input_file(file, output, verbose)
				current_file.set_protein(file)
				inputs.append(current_file)
				del current_file
			
	if hmms is not None:
		set = os.listdir(hmms)
		set.sort()
		set = [os.path.normpath(hmms + "/" + file) for file in set]
		
		hmm_pairs = []
		
		for file in set:
			if not os.path.exists(file):
				print("I can't find", file, "Are you sure this file exists in the directory you supplied?")
			else:
				hmm_pairs.append(file)
				
		if len(hmm_pairs) != len(inputs):
			print("Protein and HMM file counts differ! There must be one HMM per protein, generated from its paired protein! These must be in the same alphabetical order in their respective directories!")
			hmm_broke = True
		else:
			for h, i in zip(hmm_pairs, inputs):
				i.set_hmm(h)
	
	if hmm_broke:
		print("FastAAI can't proceed without matching HMM and protein pairs.")
		inputs = None
		return inputs
		
	total_counts = len(inputs)
	count = 0
	last_pct = 0
	
	if verbose:
		print("")
	#progress bar - possible dangerous use of the return to line start sequence.
		try:
			percentage = 0
			sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% (Genome ' + str(count) + " of " + str(total_counts) + ') at ' + curtime()+"\n")
			sys.stdout.flush()
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass
	
	results = []
	
	kmer_index_ = create_kmer_index()
	pool = multiprocessing.Pool(threads, initializer=initialize_preproc, initargs = (kmer_index_,))
	
	for res in pool.imap(do_advance, inputs):
		results.append(res)
		if verbose:
		#progress bar - possible dangerous use of the return to line start sequence.
			try:
				count += 1
				percentage = (count/total_counts)*100
				if int(percentage/2) > last_pct or partition == total_partitions:
					sys.stdout.write('\033[A')
					sys.stdout.flush()
					sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% (Genome ' + str(count) + " of " + str(total_counts) + ') at ' + curtime()+"\n")
					sys.stdout.flush()
					
				last_pct = int(percentage/2)
			except:
				#It's not really a big deal if the progress bar cannot be printed.
				pass
	
	pool.close()
	pool.join()
	
	inputs = results

	log_time = curtime()
	
	if os.path.exists(os.path.normpath(output + "/logs/" + os.path.splitext(os.path.basename(db_name))[0] + "_preprocessing_log.txt")):
		preproc_log = open(os.path.normpath(output + "/logs/" + os.path.splitext(os.path.basename(db_name))[0] + "_preprocessing_log.txt"), "a")
	else:
		preproc_log = open(os.path.normpath(output + "/logs/" + os.path.splitext(os.path.basename(db_name))[0] + "_preprocessing_log.txt"), "w")
		print("log_date", "genome_name", "started_as_a", "start_time", "end_time", "protein_translation_table", "errors", sep = "\t", file = preproc_log)
	for i in inputs:
		print(log_time, i.basename, i.initial_state, i.start_time, i.end_time, i.trans_table, i.err_log, sep = "\t", file = preproc_log)
	preproc_log.close()
	
	return inputs
	
'''
Utility functions
'''
def prepare_directories(output, status, build_or_query):
	preparation_successful = True
	
	if not os.path.exists(output):
		try:
			os.mkdir(output)
		except:
			print("")
			print("FastAAI tried to make output directory: '"+ output + "' but failed.")
			print("")
			print("Troubleshooting:")
			print("")
			print("    (1) Do you have permission to create directories in the location you specified?")
			print("    (2) Did you make sure that all directories other than", os.path.basename(output), "already exist?")
			print("")
			preparation_successful = False
	
	if preparation_successful:
		try:
			if status == 'genome':
				if not os.path.exists(os.path.normpath(output + "/" + "predicted_proteins")):
					os.mkdir(os.path.normpath(output + "/" + "predicted_proteins"))
				if not os.path.exists(os.path.normpath(output + "/" + "hmms")):
					os.mkdir(os.path.normpath(output + "/" + "hmms"))
			
			if status == 'protein':
				if not os.path.exists(os.path.normpath(output + "/" + "hmms")):
					os.mkdir(os.path.normpath(output + "/" + "hmms"))
					
			if not os.path.exists(os.path.normpath(output + "/" + "logs")):
				os.mkdir(os.path.normpath(output + "/" + "logs"))
			
			if build_or_query == "build":
				if not os.path.exists(os.path.normpath(output + "/" + "database")):
					os.mkdir(os.path.normpath(output + "/" + "database"))
			
			if build_or_query == "query":		
				if not os.path.exists(os.path.normpath(output + "/" + "results")):
					os.mkdir(os.path.normpath(output + "/" + "results"))

			
		except:
			print("FastAAI was able to create or find", output, "but couldn't make directories there.")
			print("")
			print("This shouldn't happen. Do you have permission to write to that directory?")
			
		
	return preparation_successful	

def check_out_input_files(genomes, proteins, hmms, gf, pf, hf):
	#Check only one method of supply was used per file type
	if (genomes is not None) and (gf is not None):
		print("Supply genomes either by directory or by file, not both.")
		return None
	if (proteins is not None) and (pf is not None):
		print("Supply proteins either by directory or by file, not both.")
		return None
	if (hmms is not None) and (hf is not None):
		print("Supply HMMs either by directory or by file, not both.")
		return None
	
	#check that not both proteins and genomes supplied in any combo.
	if ((genomes is not None) and (pf is not None))\
	or ((gf is not None) and (proteins is not None))\
	or ((genomes is not None) and (proteins is not None))\
	or ((gf is not None) and (pf is not None)):
		print("Supply either genomes or proteins, not both. You can supply proteins and HMMs, but not genomes and proteins.")
		return None
		
	#Check that if hmms are given, so are proteins
	if (hmms is not None) or (hf is not None):
		if (proteins is None) and (pf is None):
			print("If you supply HMMs, you also have to supply the proteins from which they were generated.")
			return None
		
	#Determine status
	if (genomes is not None) or (gf is not None):
		print("Starting from genomes")
		start = 'genome'
			
	else:
		if 	(hmms is not None) or (hf is not None):
			print("Starting from proteins and HMMs")
			start = 'protein and HMM'
			
		else:
			print("Starting from proteins")
			start = 'protein'
			
	return start


#Build DB from genomes
	
def unique_kmers(seq, ksize):
	n_kmers = len(seq) - ksize + 1
	kmers = []
	for i in range(n_kmers):
		kmers.append(kmer_index[seq[i:i + ksize]])
	#We care about the type because we're working with bytes later.
	return np.unique(kmers).astype(np.int32)

#Quickly creates a dict of all poss. tetramers in a fixed, alphabetical order.
#This can be used to index kmers so that the indices are identical (and thus interchangable) on separate runs of this program.	
def create_kmer_index():
	valid_chars = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', '*']
	#This meshgrid method will produce all unique tetramers from AAAA to **** in a consistent order.
	#Rightmost char to leftmost, A to * in the same order as valid_chars
	kmer_index_ = np.stack(np.meshgrid(valid_chars, valid_chars, valid_chars, valid_chars), -1).reshape(-1, 4)
	#Unless someone is passing more than 2.1 billion genomes, int32 will be enough.
	kmer_index_ = dict(zip([''.join(kmer_index_[i,]) for i in range(0, kmer_index_.shape[0])], np.arange(kmer_index_.shape[0], dtype = np.int32)))
	
	return kmer_index_

def split_seq(seq, num_grps):
	newseq = []
	splitsize = 1.0/num_grps*len(seq)
	for i in range(num_grps):
		newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
	return newseq
	
#gives the max and min index needed to split a list of (max_val) genomes into 
def split_indicies(max_val, num_grps):
	newseq = []
	splitsize = 1.0/num_grps*max_val
	for i in range(num_grps):
		newseq.append(((round(i*splitsize)), round((i+1)*splitsize)))
	return newseq
	
def list_to_index_dict(list):
	result = {}
	counter = 0
	for item in list:
		result[item] = counter
		counter += 1
	return result

def generate_accessions_index():
	list_of_poss_accs = list_to_index_dict(['PF01780.19', 'PF03948.14', 'PF17144.4', 'PF00830.19', 'PF00347.23', 'PF16906.5', 'PF13393.6',
	'PF02565.15', 'PF01991.18', 'PF01984.20', 'PF00861.22', 'PF13656.6', 'PF00368.18', 'PF01142.18', 'PF00312.22', 'PF02367.17',
	'PF01951.16', 'PF00749.21', 'PF01655.18', 'PF00318.20', 'PF01813.17', 'PF01649.18', 'PF01025.19', 'PF00380.19', 'PF01282.19',
	'PF01864.17', 'PF01783.23', 'PF01808.18', 'PF01982.16', 'PF01715.17', 'PF00213.18', 'PF00119.20', 'PF00573.22', 'PF01981.16',
	'PF00281.19', 'PF00584.20', 'PF00825.18', 'PF00406.22', 'PF00177.21', 'PF01192.22', 'PF05833.11', 'PF02699.15', 'PF01016.19',
	'PF01765.19', 'PF00453.18', 'PF01193.24', 'PF05221.17', 'PF00231.19', 'PF00416.22', 'PF02033.18', 'PF01668.18', 'PF00886.19',
	'PF00252.18', 'PF00572.18', 'PF00366.20', 'PF04104.14', 'PF04919.12', 'PF01912.18', 'PF00276.20', 'PF00203.21', 'PF00889.19',
	'PF02996.17', 'PF00121.18', 'PF01990.17', 'PF00344.20', 'PF00297.22', 'PF01196.19', 'PF01194.17', 'PF01725.16', 'PF00750.19',
	'PF00338.22', 'PF00238.19', 'PF01200.18', 'PF00162.19', 'PF00181.23', 'PF01866.17', 'PF00709.21', 'PF02006.16', 'PF00164.25',
	'PF00237.19', 'PF01139.17', 'PF01351.18', 'PF04010.13', 'PF06093.13', 'PF00828.19', 'PF02410.15', 'PF01176.19', 'PF02130.17',
	'PF01948.18', 'PF01195.19', 'PF01746.21', 'PF01667.17', 'PF03874.16', 'PF01090.19', 'PF01198.19', 'PF01250.17', 'PF17136.4',
	'PF06026.14', 'PF03652.15', 'PF04019.12', 'PF01201.22', 'PF00832.20', 'PF01264.21', 'PF03840.14', 'PF00831.23', 'PF00189.20',
	'PF02601.15', 'PF01496.19', 'PF00411.19', 'PF00334.19', 'PF00687.21', 'PF01157.18', 'PF01245.20', 'PF01994.16', 'PF01632.19',
	'PF00827.17', 'PF01015.18', 'PF00829.21', 'PF00410.19', 'PF00833.18', 'PF00935.19', 'PF01992.16'])
	
	return list_of_poss_accs

#Master function for building or adding to a DB with genomes.	
def add_inputs(output_path, parent_path, existing_index, threads, verbose, prep_args):
	
	genomes, proteins, hmms, gf, pf, hf, db_name = prep_args[0], prep_args[1], prep_args[2], prep_args[3], prep_args[4], prep_args[5], prep_args[6]
	
	print("")
	print("FastAAI is formatting your files to be saved to your database.")
	
	#Let's push this to the inputs section.
	inputs = advance_inputs(genomes = genomes, proteins = proteins, hmms = hmms, genomes_file = gf, proteins_file = pf, hmms_file = hf, output = output_path, threads = threads, verbose = verbose, db_name = db_name)
	
	if inputs is None:
		return False
		
	kmer_index = None
	
	#global genome_index
	genome_index = {}
	next_index = 0
	
	#Build upon the genome indexing of an existing DB
	if existing_index is not None:
		genome_index = existing_index
		#zero indexing makes this the next number to add.
		next_index = len(existing_index)
		
	final_db = fastaai_database(parent_path)
	final_db.activate_connection()
	final_db.initialize_parent_database()
		
	#This goes to the genome_index table
	protein_counts_to_add = []
	genome_acc_kmer_counts_to_add = []
		
	acc_index = generate_accessions_index()
		
	readied_kmers_by_acc = defaultdict(lambda: defaultdict(lambda: None))
	
	#unique_accessions = set()
	for file in inputs:
		
		genome = file.basename
		
		#Collect all of the accessions actually found. Will usually be 122 for reasonably sized datasets.
		#unique_accessions = unique_accessions.union(set(file.best_hits.values()))
		#Avoid adding duplicate genomes
		if genome not in genome_index:
			protein_counts_to_add.append((genome, next_index, file.protein_count))
			for prot in file.protein_kmer_count:
				genome_acc_kmer_counts_to_add.append((next_index, acc_index[prot], file.protein_kmer_count[prot]))
			genome_index[genome] = next_index
			next_index += 1
		
		this_index = genome_index[genome]
		for acc in file.best_hits_kmers:
			readied_kmers_by_acc[acc][this_index] = file.best_hits_kmers[acc]
		#Clean up space
		file.best_hits_kmers = None

	inputs = None
	
	#Default dicts can't be pickled.
	readied_kmers_by_acc = dict(readied_kmers_by_acc)
	
	genomes_per_acc = {}
	for acc in readied_kmers_by_acc:
		readied_kmers_by_acc[acc] = dict(readied_kmers_by_acc[acc])
		genomes_per_acc[acc] = list(readied_kmers_by_acc[acc].keys())
		final_db.add_genomes_first(acc, readied_kmers_by_acc[acc])
		readied_kmers_by_acc[acc] = None
		
	readied_kmers_by_acc = None
		
	add_genomes = "INSERT OR REPLACE INTO genome_index VALUES (?, ?, ?)"
	add_proteins = "INSERT OR REPLACE INTO genome_acc_kmer_counts VALUES (?, ?, ?)"
	
	final_db.cursor.executemany(add_genomes, protein_counts_to_add)
	final_db.cursor.executemany(add_proteins, genome_acc_kmer_counts_to_add)
	final_db.connection.commit()
	
	final_db.cursor.execute("CREATE INDEX IF NOT EXISTS kmer_acc ON genome_acc_kmer_counts (genome, accession);")
	final_db.connection.commit()
	
	protein_counts_to_add = None
	genome_acc_kmer_counts_to_add = None
			
	unique_accessions = list(genomes_per_acc.keys())
	child_args = []
	for i in range(0, len(unique_accessions)):
		accession = unique_accessions[i]
		name = "accession_" + unique_accessions[i] + "_partition_" + str(i)
		child_path = os.path.normpath(output_path+"/temp")
		child_args.append([accession, name, child_path, parent_path, genomes_per_acc[accession], genome_index])
			
	print("")
	print("Formatting data to add to database at", curtime())
	
	#Add partition, output, parent DB data.
	if not os.path.exists(os.path.normpath(output_path+"/temp")):
		try:
			os.mkdir(os.path.normpath(output_path+"/temp"))
		except:
			print("Output directory failed to create! Cannot continue.")
			return False
	
	if verbose:
		print("")
		count = 0
		total_counts = len(child_args)
		try:
			log_time = curtime()
			percentage = (count/total_counts)*100
			sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(count) + " of " + str(total_counts) + ' ) at ' + curtime() + "\n")
			sys.stdout.flush()
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass
	
	last_pct = 0
	
	quiverfull = []
	
	pool = multiprocessing.Pool(threads)
	
	for result in pool.imap_unordered(produce_children, child_args):
		acc = result[0]
		child = result[1]
		
		quiverfull.append([acc, child])
		
		if verbose:
			count += 1
			try:
				percentage = (count/total_counts)*100
				log_time = curtime()
				sys.stdout.write('\033[A')
				sys.stdout.flush()
				sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(count) + " of " + str(total_counts) + ' done at '+ curtime() + " )\n")
				sys.stdout.flush()
			except:
				#It's not really a big deal if the progress bar cannot be printed.
				pass
		
	pool.close()
	pool.join()
	
	print("")
	print("Adding data to final database.")
	
	if verbose:
		print("")
		
		count = 0
		total_counts = len(child_args)
		try:
			percentage = (count/total_counts)*100
			
			("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(count) + " of " + str(total_counts) + ' done at '+ curtime() + " )\n")
			sys.stdout.flush()		
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass
	
	last_pct = 0
	
	for result in quiverfull:
		acc = result[0]
		child = result[1]
		final_db.add_child_to_parent(acc, child)

		if verbose:
			count += 1
			try:
				percentage = (count/total_counts)*100
				log_time = curtime()
				sys.stdout.write('\033[A')
				sys.stdout.flush()
				sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(count) + " of " + str(total_counts) + ' done at '+ curtime() + " )\n")
				sys.stdout.flush()
			except:
				#It's not really a big deal if the progress bar cannot be printed.
				pass
	
	
	print("")
	#print("Cleaning up...")
	#final_db.connection.execute("VACUUM")
	
	final_db.close_connection()
	
	os.rmdir(os.path.normpath(output_path+"/temp"))
	
	return True

#genome_index is global already
def produce_children(args):
	acc = args[0]
	partition = args[1]
	output_base = args[2]
	parent_db = args[3]
	genomes_in_this_acc = args[4]
	genome_index = args[5]
	
	parental_database = fastaai_database(parent_db)
	
	sql_friendly_accession = acc.replace('.', '_')
	
	read_parent_sql = "SELECT * FROM " + sql_friendly_accession + "_genomes WHERE genome IN ({genomes})".format(genomes=','.join(['?']*len(genomes_in_this_acc)))
	
	parental_database.activate_connection()
	
	genomes_for_this_acc = dict(parental_database.cursor.execute(read_parent_sql, genomes_in_this_acc).fetchall())
	
	parental_database.close_connection()
		
	child_db = os.path.normpath(output_base + "/" + partition + ".db")
	
	this_child = child_database(child_db, parent_db)
	
	this_child.activate_child_connection()
	#this_child.initialize_child_database()
	this_child.activate_parent_connection()
	
	#Keys are genomes as indices, values are numpy arrays of kmers. This makes tuples.
	#this_child.add_genomes_first(acc, zip(genomes_for_this_acc.keys(), genomes_for_this_acc.values()))
	
	#Here's where we add the genomes as such to the children, too.
	readied_kmers = defaultdict(lambda: [])
	for genome in genomes_for_this_acc:
		for kmer in genomes_for_this_acc[genome]:
			readied_kmers[kmer].append(genome)
		#cleanup space
		genomes_for_this_acc[genome] = None
	
	del genomes_for_this_acc
	
	readied_kmers = dict(readied_kmers)
	for kmer in readied_kmers:
		readied_kmers[kmer] = np.array(readied_kmers[kmer], dtype = np.int32)
	
	sql_friendly_accession = this_child.add_accession(acc, readied_kmers)
		
	this_child.close_parent_connection()
	this_child.close_child_connection()
	
	del readied_kmers
	
	return [sql_friendly_accession, child_db]
	
#Build or add to a FastAAI DB
def build_db_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module allows you to create a FastAAI database from one or many genomes, proteins, or proteins and HMMs, or add these files to an existing one.
	
	Supply genomes OR proteins OR proteins AND HMMs as inputs.
	
	If you supply genomes, FastAAI will predict proteins from them, and HMMs will be created from those proteins
	If you supply only proteins, FastAAI will create HMM files from them, searching against FastAAI's internal database
	If you supply proteins AND HMMs, FastAAI will directly use them to build the database.\n
	You cannot supply both genomes and proteins
	''')

	parser.add_argument('-g', '--genomes',  dest = 'genomes', default = None, help =  'A directory containing genomes in FASTA format.')
	parser.add_argument('-p', '--proteins', dest = 'proteins', default = None, help = 'A directory containing protein amino acids in FASTA format.')
	parser.add_argument('-m', '--hmms',     dest = 'hmms', default = None, help =     'A directory containing the results of an HMM search on a set of proteins.')
	parser.add_argument('-d', '--database', dest = 'db_name', default = "FastAAI_database.sqlite.db", help =  'The name of the database you wish to create or add to. The database will be created if it doesn\'t already exist and placed in the output directory. FastAAI_database.sqlite.db by default.')
	
	parser.add_argument('-o', '--output',   dest = 'output', default = "FastAAI", help = 'The directory to place the database and any protein or HMM files FastAAI creates. By default, a directory named "FastAAI" will be created in the current working directory and results will be placed there.')
						
	parser.add_argument('--genome_file',    dest = 'gf', default = None, help =  'Alternative way to supply genomes. A file containing paths to your genome files, 1 per line.')
	parser.add_argument('--protein_file',   dest = 'pf', default = None, help = 'Alternative way to supply proteins. A file containing paths to your protein files, 1 per line.')
	parser.add_argument('--hmm_file',       dest = 'hf', default = None, help =     'Alternative way to supply HMMs. A file containing paths to your HMM files, 1 per line.')
		
	parser.add_argument('--threads',  dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')
	parser.add_argument('--verbose',        dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
		
	args, unknown = parser.parse_known_args()
	
	return parser, args
		
def build_db(genomes, proteins, hmms, db_name, output, threads, gf, pf, hf, verbose):
	
	start = check_out_input_files(genomes, proteins, hmms, gf, pf, hf)
	
	#If something failed, we stop.
	if start is None:
		return False
	
	good_to_go = prepare_directories(output, start, "build")
	
	if not good_to_go:
		return False
	
	#Check if the db contains path info. Incl. windows version.
	if "/" not in db_name and "\\" not in db_name:
		final_database = os.path.normpath(output + "/database/" + db_name)
	else:
		#If the person insists that the db has a path, let them.
		final_database = db_name
		
	#We'll skip trying this if the file already exists.
	existing_genome_IDs = None
	try:
		if os.path.exists(final_database):
			parent = fastaai_database(final_database)
			parent.activate_connection()
			
			existing_genome_IDs = {}
			sql_command = "SELECT genome, gen_id FROM genome_index"
			for result in parent.cursor.execute(sql_command).fetchall():
				genome = result[0]
				id = int(result[1])
				existing_genome_IDs[genome] = id
			
			parent.close_connection()
	except:
		print("You specified an existing file to be a database, but it does not appear to be a FastAAI database.")
		print("FastAAI will not be able to continue. Please give FastAAI a different database name and continue.")
		print("Exiting.")
		return False
	
	
	prep_args = [genomes, proteins, hmms, gf, pf, hf, db_name]
	
	#inputs, output_path, parent_path, existing_index, threads
	success = add_inputs(output, final_database, existing_genome_IDs, threads, verbose, prep_args)
	
	if success:
		print("Database build complete!")
		
	return success


#DB query functionality - unlimited version
def do_query_vs_target_aai_only(query_name, target_name, threads, output, precision, verbose):
	if not os.path.exists(os.path.normpath(output+"/temp")):
		os.mkdir(os.path.normpath(output+"/temp"))
		
	if precision == "low":
		jacc_precision = np.float16
	if precision == "med":
		jacc_precision = np.float32
	if precision == "high":
		jacc_precision = np.float64
	
	#Save the file paths.
	query = fastaai_database(query_name)
	target = fastaai_database(target_name)
	
	query.activate_connection()
	query.just_accessions()
	query_len = query.cursor.execute("SELECT Count(*) FROM genome_index").fetchall()[0][0]
	#query.close_connection()
	target.activate_connection()
	target.just_accessions()
	target_len = target.cursor.execute("SELECT Count(*) FROM genome_index").fetchall()[0][0]
	#target.close_connection()
	
	print("FastAAI will search", query_len, "query genomes against", target_len, "target genomes.")
	
	print("")
	print("FastAAI is preparing your AAI search... ", end = '', flush = True)
	
	accessions_in_common = list(set(query.accessions).intersection(target.accessions))
	
	query.accessions = None
	target.accessions = None
	
	query.close_connection()
	target.close_connection()
	
	load_args = [(query, target, acc) for acc in accessions_in_common]
	
	loads = []
	ordered_accs = []
	
	pool = multiprocessing.Pool(threads)
	
	for result in pool.imap(load_getter, load_args):
		load = result[0]
		acc = result[1]
		#Load will be None if the accession is in both query and target, but they still don't share even a single Kmer. Unlikely, but it happened once, so it WILL happen again.
		if load is not None:
			loads.append(load)
			ordered_accs.append(acc)
	
	pool.close()
	pool.join()
	
	loads = np.array(loads)
	ordered_accs = np.array(ordered_accs)
	
	order = loads.argsort()[::-1]
	
	loads = loads[order]
	ordered_accs = ordered_accs[order]
	
	load_balancer = {}
	accs_per_load = {}
	for i in range(0, threads):
		load_balancer[i] = 0
		accs_per_load[i] = []
	
	for i in range(0, loads.shape[0]):
		index = list(load_balancer.values()).index(min(list(load_balancer.values())))
		#print(index, load)
		load_balancer[index] += loads[i]
		accs_per_load[index].append(int(ordered_accs[i]))
		
	del loads
	del ordered_accs
		
	print("done!")
	if verbose:
		print("FastAAI has balanced the workload of calculating AAI from your data.")
		for index in accs_per_load:
			print("Thread", index, "will handle", len(accs_per_load[index]), "accessions.")
	print("FastAAI is beginning the calculation of AAI between your query and target genomes.")
	
	del load_balancer
	
	input_queue = multiprocessing.Queue()
	output_queue = multiprocessing.Queue()
	
	for thread in accs_per_load:
		input_queue.put(accs_per_load[thread])
		
	for i in range(0, threads):
		input_queue.put('STOP')
	
	for i in range(0, threads):
		multiprocessing.Process(target=accession_worker, args=(input_queue, output_queue, query, target, query_len, target_len, jacc_precision)).start()
	
	print("")
	
	results = np.zeros(shape = (query_len, target_len), dtype = jacc_precision)
	
	#Counter to keep the threads running until the whole process is done.
	donezo = threads
	while donezo > 0:
		row = output_queue.get()
		try:
			results[row[0]] += row[1]
		except:
			donezo -= 1
	
	print("AAI calculations complete. Formatting results for writing.")
	
	#global glob_prec
	#glob_prec = jacc_precision
	
	rdb_name = os.path.normpath(output+"/temp/aai_calc_db.db")
	rdb = calculation_database(rdb_name, precision)
	rdb.activate_connection()
	rdb.initialize_database()
	
	#Get the data ready for passing to children...
	
	results = np.split(results, query_len, axis = 0)
	
	insertable = []
	#iterate over results and turn them into tuples.
	for i in range(0, query_len):
		insertable.append((i, results[i].tobytes()))
		results[i] = None
	
	rdb.cursor.executemany("INSERT INTO jaccards VALUES (?, ?)", (insertable))
	rdb.connection.commit()
	
	rdb.close_connection()
	
	del insertable
	del results
	
	#Now we split the query genomes into chunk and have threads process each chunk in parallel with its respective shared prot counts.
	query_chunks = split_indicies(query_len, threads)
	query_args = [([rdb_name], query_chunks[i], output, query, target, precision) for i in range(0, threads)]
	
	print("Results formatted. Writing results starting at", curtime())
	
	pool = multiprocessing.Pool(threads)
	
	pool.map(finish_jaccards, query_args)
	
	pool.close()
	pool.join()
	
	os.remove(rdb_name)
	
	print("FastAAI complete! Results at:", os.path.normpath(output+"/results/"))
	
	return None
	
#Assess the number of comparisons that will have to be made to complete an accession so that balanced loads can be passed to threads
def load_getter(args):
	query, target, accession = args[0], args[1], args[2]
	query.activate_connection()
	target.activate_connection()

	original_index = generate_accessions_index()
	accession_inverter = {}
	for acc in original_index:
		sql_friendly_accession = acc.replace(".", "_")
		accession_inverter[original_index[acc]] = sql_friendly_accession

	sql_friendly_accession = accession_inverter[accession].replace('.', '_')
	sql = "SELECT kmer FROM "+ sql_friendly_accession
	query.cursor.row_factory = lambda cursor, row: row[0]
	#query_kmers = set(query.cursor.execute(sql).fetchall()).intersection()
	target.cursor.row_factory = lambda cursor, row: row[0]
	#target_kmers = target.cursor.execute(sql).fetchall()
	
	shared_kmers = list(set(query.cursor.execute(sql).fetchall()).intersection(target.cursor.execute(sql).fetchall()))
	query.cursor.row_factory = None
	target.cursor.row_factory = None
	
	bytes_sql = "SELECT sum(length(genomes)) FROM " + sql_friendly_accession + " WHERE kmer IN ({kmers})".format(kmers=','.join(['?']*len(shared_kmers)))
	
	if len(shared_kmers) > 0:
		tgt_res = target.cursor.execute(bytes_sql, shared_kmers).fetchone()[0]
		query_res = query.cursor.execute(bytes_sql, shared_kmers).fetchone()[0]
		#This if *should* always happen, if it gets checked.
		if tgt_res is not None and query_res is not None:
			load = int(tgt_res/(4096) * query_res/(4096))
		else:
			load = None
	else:
		load = None
	
	query.close_connection()
	target.close_connection()
		
	return [load, accession]
	
def accession_worker(in_queue, out_queue, query, target, qlen, tlen, prec):
	original_index = generate_accessions_index()
	accession_inverter = {}
	for acc in original_index:
		sql_friendly_accession = acc.replace(".", "_")
		accession_inverter[original_index[acc]] = sql_friendly_accession

	query.activate_connection()
	target.activate_connection()
	query.load_genome_index()
	target.load_genome_index()
	
	for my_accessions in iter(in_queue.get, 'STOP'):
	
		#print(my_accessions)
		
		target.load_accessions(permitted_accessions = my_accessions)
		query.load_accessions(permitted_accessions = my_accessions)
		
		query_data = {}
		target_data = {}
		
		for acc in my_accessions:

			sql_friendly_accession = accession_inverter[acc].replace('.', '_')
			
			query_data[acc] = dict(query.cursor.execute("SELECT * FROM "+sql_friendly_accession+"_genomes").fetchall())
			
			query.cursor.row_factory = lambda cursor, row: row[0]
			selected_kmers = list(query.cursor.execute("SELECT kmer FROM "+sql_friendly_accession).fetchall())
			query.cursor.row_factory = None
			
			target_sql = "SELECT * FROM " + sql_friendly_accession + " WHERE kmer in ({kmers})".format(kmers=','.join(['?']*len(selected_kmers)))
			target_data[acc] = dict(target.cursor.execute(target_sql, selected_kmers).fetchall())
			
		target_kmer_cts_by_acc = {}
		for acc in my_accessions:
			target_kmer_cts_by_acc[acc] = np.zeros(tlen, dtype = np.int16)
			
		for genome in target.gak:
			for acc in target.gak[genome]:
				target_kmer_cts_by_acc[acc][genome] = target.gak[genome][acc]
					
		#No longer needed.
		target.gak = None
		#We want each thread to report every single genome
		for genome in query.gak:
			#count += 1
			#print("Thread", my_thread, "genome", count, "of", total)
			these_jaccards = np.zeros(tlen, dtype = np.float64)
			for acc in query.gak[genome]:
				these_intersections = np.zeros(tlen, dtype = np.int16)
				query_kmers = query_data[acc][genome]
				query_kmer_ct = query_kmers.shape
				for kmer in query_kmers:
					if kmer in target_data[acc]:
						these_intersections[target_data[acc][kmer]] += 1
						
				these_jaccards += np.divide(these_intersections, np.subtract(np.add(query_kmer_ct, target_kmer_cts_by_acc[acc]), these_intersections))
				
			out_queue.put([genome, these_jaccards])
			
	target.close_connection()
	query.close_connection()
	out_queue.put("Based")

	return None
	
def finish_jaccards(args):
	partial_dbs, my_query_genomes, output, query, target, prec = args[0], args[1], args[2], args[3] ,args[4], args[5]
	#Load protein counts
	#for each genome, query each partial and sum matching genomes, then divide by shared counts.
	
	query.activate_connection()
	target.activate_connection()
	query.load_genome_index()
	target.load_genome_index()
	
	selected_query_genomes = range(my_query_genomes[0], my_query_genomes[1])

	offset = my_query_genomes[0]
	
	target_len = len(target.genome_index)
	query_len = my_query_genomes[1] - my_query_genomes[0]
	
	#get shared protein counts
	query.load_accessions(permitted_genomes = selected_query_genomes)
	
	max_acc = 122
	
	query_set = np.zeros(shape = (query_len, max_acc), dtype = np.int16)
	
	for g in query.gak:
		query_set[(g-offset), list(query.gak[g])] += 1
	
	target_set = np.zeros(shape = (max_acc, len(target.genome_index)), dtype = np.int16)
	
	target.load_accessions()
	
	target_protein_counts = np.zeros(target_len, dtype = np.int16)
	for t in target.gak:
		target_set[list(target.gak[t]), t] += 1
		target_protein_counts[t] = len(target.gak[t])
		
	#This will be used to divide the jaccs and such. If disk, then disk, tho...
	shared_prot_counts_by_genome = np.dot(query_set, target_set)
	
	del query_set
	del target_set
	
	target.gak = None
	
	query.close_connection()
	target.close_connection()
	
	activated_DBs = []
	idx = 0
	for db in partial_dbs:
		activated_DBs.append(calculation_database(db, prec))
		activated_DBs[idx].activate_connection()
		idx += 1

	
	for genome in selected_query_genomes:
		sql = "SELECT jaccards FROM jaccards WHERE genome="+str(genome)
		total_jaccs = np.zeros(target_len, dtype = np.float64)
		shared_acc_counts = shared_prot_counts_by_genome[genome - offset]
		for db in activated_DBs:
			result = db.cursor.execute(sql).fetchone()[0]
			total_jaccs += result
			
		total_jaccs = np.divide(total_jaccs, shared_acc_counts)
	
		aai_est = numpy_kaai_to_aai(total_jaccs)
		
		no_hit = np.where(shared_acc_counts == 0)
		#Actual hits is already stored in shared_acc_counts
		possible_hits = np.minimum(len(query.gak[genome]), target_protein_counts).astype(str)
		
		total_jaccs = np.round(total_jaccs, 4).astype(str)
		
		shared_acc_counts = shared_acc_counts.astype(str)
			
		total_jaccs[no_hit] = "N/A"
		aai_est[no_hit] = "N/A"
		shared_acc_counts[no_hit] = "N/A"
		possible_hits[no_hit] = "N/A"
		
		name = query.reverse_genome_index[genome]
		
		output_file = output +"/results/"+name+"_results.txt"
		fh = open(output_file, "w")

		for tgt in range(0, target_len):
			target_name = target.reverse_genome_index[tgt]
			if target_name == name:
				fh.write(name+"\t"+target_name+"\t"+"100.0"+"\t"+"0.0"+"\t"+shared_acc_counts[tgt]+"\t"+possible_hits[tgt]+"\t"+"100.0"+"\n")
			else:
				fh.write(name+"\t"+target_name+"\t"+total_jaccs[tgt]+"\t"+"N/A"+"\t"+shared_acc_counts[tgt]+"\t"+possible_hits[tgt]+"\t"+aai_est[tgt]+"\n")
		
		fh.close()
		
		#Write partial to file, here.
	
	for db in activated_DBs:
		db.close_connection()
	
	return None
	
	
#Here's the DB SQL querying functionality/limited version.
def do_query_vs_target_sql(query, target, threads, output, verbose, do_stdev):
	#Save the file paths.
	query_name, target_name = query, target
	
	query = fastaai_database(query_name)
	query.activate_connection()
	query.load_genome_index()
	query.just_accessions()
	
	converter = generate_accessions_index()
	acc_sql = "SELECT name FROM sqlite_master WHERE type='table'"
	tables = [item[0] for item in query.cursor.execute(acc_sql).fetchall()]
	cleaned_tables = []
	for table in tables:
		if table.endswith("_genomes"):
			acc_name = table.split("_genomes")[0]
			acc_name = acc_name.replace("_", ".")
			index = acc_name
			cleaned_tables.append((table, index))

	del tables
			
	#Go through tables and load data.
	query_acc_kmers = defaultdict(dict)
	
	sys.stdout.write("\n")
	sys.stdout.write("Loading query data at " + curtime() + " ...\n")
	sys.stdout.flush()
	
	for tab_idx in cleaned_tables:
		table = tab_idx[0]
		accession = tab_idx[1]
		for result in query.cursor.execute("SELECT * FROM " + table).fetchall():
			query_acc_kmers[result[0]][accession] = result[1]
		
	query.close_connection()
		
		
	sys.stdout.write("\n")
	sys.stdout.write("Loading target data at " + curtime() + " ...\n")
	sys.stdout.flush()
		
	target = fastaai_database(target_name)
	target.activate_connection()
	target.load_genome_index()
	target.load_accessions()
	target.close_connection()
		
	query_args = []
	for genome in query_acc_kmers:
		query_args.append((target, query.reverse_genome_index[genome], query_acc_kmers[genome], os.path.normpath(output+"/results")))
	
	detected_query_accs = query.accessions
	query_length = len(query.genome_index)
	
	#Cleanup
	del query
	del query_acc_kmers
	
	#global target_kmer_cts
	target_kmer_cts = {}
	
	target_len = len(target.gak)
	
	for accession in np.intersect1d(detected_query_accs, target.accessions):
		target_kmer_cts[accession] = np.zeros(target_len, dtype = np.int16)
		for g in target.gak:
			if accession in target.gak[g]:
				target_kmer_cts[accession][g] = target.gak[g][accession]
	
	#global target_protein_counts
	target_protein_counts = np.zeros(target_len, dtype = np.int16)
	for g in target.gak:
		target_protein_counts[g] = len(target.gak[g])
	
	target_length = len(target.gak)
	
	target.gak = None
	
	#Should just load the stuff then straightforward sql
	sys.stdout.write("\n")
	sys.stdout.write("FastAAI will search "+ str(query_length) + " query genomes against " + str(target_length) + " target genomes.\n")
	sys.stdout.write("\n")
	
	count = 0
	total = len(query_args)
	
	sys.stdout.write("Beginning AAI calculation at " + curtime())
	
	if verbose:
		print("")
	#progress bar - possible dangerous use of the return to line start sequence.
		try:
			percentage = 0
			sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% (Query genome ' + str(count) + " of " + str(total) + ' done at '+curtime()+')\n')
			sys.stdout.flush()
			last_pct = 0
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass

	pool = multiprocessing.Pool(threads, initializer = sql_query_thread_starter, initargs = (target_kmer_cts, target_protein_counts,))

	#Process as we go.
	if do_stdev:
		for file in pool.imap(do_sql_query, query_args):
			if verbose:
			#progress bar - possible dangerous use of the return to line start sequence.
				try:
					count += 1
					percentage = (count/total)*100
					if int(percentage/2) > last_pct or count == total:
						sys.stdout.write('\033[A')
						sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% (Query genome ' + str(count) + " of " + str(total) + ' done at '+curtime()+')\n')
						sys.stdout.flush()
					last_pct = int(percentage/2)
				except:
					#It's not really a big deal if the progress bar cannot be printed.
					pass
		
		pool.close()
		pool.join()
	else:
	
		for file in pool.imap(do_sql_query_no_SD, query_args):
			
			if verbose:
				#progress bar - possible dangerous use of the return to line start sequence.
				try:
					count += 1
					percentage = (count/total)*100
					if int(percentage/2) > last_pct or count == total:
						sys.stdout.write('\033[A')
						sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% (Query genome ' + str(count) + " of " + str(total) + ' done at '+curtime()+')\n')
						sys.stdout.flush()
					last_pct = int(percentage/2)
				except:
					#It's not really a big deal if the progress bar cannot be printed.
					pass
		
		pool.close()
		pool.join()
	
	print("AAI calculation complete! Results at:", os.path.normpath(output+"/results"))

	return None
	
#This can also take the genomes-first formatted prots in the DB and search them memory-efficiently, if not time efficiently.
def do_sql_query(args):
	kmer_index = create_kmer_index()
	accession_index = generate_accessions_index()
	#database, file.basename, file.best_hits_kmers, os.path.normpath(output+"/temp")
	database, name, acc_kmers, temp_out = args[0],args[1],args[2],args[3]
	
	database.activate_connection()
	
	res_ct = 0
	target_len = len(database.genome_index)
	
	results = np.zeros(shape = (len(acc_kmers), target_len), dtype = np.float64)
	row = 0
	
	shared_acc_counts = np.zeros(target_len, dtype = np.int16)
	
	for accession in acc_kmers:
		acc_index = accession_index[accession]
		sql_friendly_accession = accession.replace(".", "_")
		if acc_index in database.accessions:
			temp_name = name+"_"+sql_friendly_accession
			#The accession was found for this target genome, for each tgt genome.
			shared_acc_counts[np.nonzero(target_kmer_cts[acc_index])] += 1
			
			#Miguel issue redev.
			these_kmers = [(int(kmer),) for kmer in acc_kmers[accession]]
			temp_name = name+"_"+sql_friendly_accession
			temp_name = temp_name.replace(".", "_")
			
			temp_tab = "CREATE TEMP TABLE " + temp_name + " (kmer INTEGER)"
			database.cursor.execute(temp_tab)
			database.connection.commit()
			insert_table = "INSERT INTO " + temp_name + " VALUES (?)"
			database.cursor.executemany(insert_table, these_kmers)
			database.connection.commit()
			join_and_select_sql = "SELECT genomes FROM " + temp_name + " INNER JOIN " + sql_friendly_accession + " ON "+ temp_name+".kmer = " + sql_friendly_accession+".kmer;"
			
			these_intersections = np.zeros(target_len, dtype = np.int16)
			#sql_query = "SELECT genomes FROM " + sql_friendly_accession + " WHERE kmer in ({kmers})".format(kmers=','.join(['?']*len(these_kmers)))
			for result in database.cursor.execute(join_and_select_sql).fetchall():
				these_intersections[result] += 1
		
			results[row] = np.divide(these_intersections, np.subtract(np.add(acc_kmers[accession].shape[0], target_kmer_cts[acc_index]), these_intersections))
		
		row += 1
	
	database.close_connection()
	
	#These are the jacc averages
	jaccard_averages = np.divide(np.sum(results, axis = 0), shared_acc_counts)
	
	#Get the differences from the mean per hit
	results = results - jaccard_averages
	#Square them
	results = np.square(results)
	#Sum squares and divide by shared acc. count, the sqrt to get SD.
	jaccard_SDs = np.sqrt(np.divide(np.sum(results, axis = 0), shared_acc_counts))

	aai_est = numpy_kaai_to_aai(jaccard_averages)
	
	no_hit = np.where(shared_acc_counts == 0)
	#Actual hits is already stored in shared_acc_counts
	possible_hits = np.minimum(len(acc_kmers), target_protein_counts).astype(str)
	
	
	jaccard_averages = np.round(jaccard_averages, 4).astype(str)
	jaccard_SDs = np.round(jaccard_SDs, 4).astype(str)
	
	shared_acc_counts = shared_acc_counts.astype(str)
		
	jaccard_averages[no_hit] = "N/A"
	aai_est[no_hit] = "N/A"
	jaccard_SDs[no_hit] = "N/A"
	shared_acc_counts[no_hit] = "N/A"
	possible_hits[no_hit] = "N/A"
	
	output_file = temp_out +"/"+name+"_results.txt"
	fh = open(output_file, "w")

	for target in range(0, target_len):
		target_name = database.reverse_genome_index[target]
		if target_name == name:
			fh.write(name+"\t"+target_name+"\t"+"100.0"+"\t"+"0.0"+"\t"+shared_acc_counts[target]+"\t"+possible_hits[target]+"\t"+"100.0"+"\n")
		else:
			fh.write(name+"\t"+target_name+"\t"+jaccard_averages[target]+"\t"+jaccard_SDs[target]+"\t"+shared_acc_counts[target]+"\t"+possible_hits[target]+"\t"+aai_est[target]+"\n")
	
	fh.close()
	
	return output_file

#This can also take the genomes-first formatted prots in the DB and search them memory-efficiently, if not time efficiently.
def do_sql_query_no_SD(args):
	kmer_index = create_kmer_index()
	accession_index = generate_accessions_index()
	#database, file.basename, file.best_hits_kmers, os.path.normpath(output+"/temp")
	database, name, acc_kmers, temp_out = args[0],args[1],args[2],args[3]
	
	database.activate_connection()
	
	res_ct = 0
	target_len = len(database.genome_index)
	
	results = np.zeros(shape = target_len, dtype = np.float64)
	#row = 0
	
	shared_acc_counts = np.zeros(target_len, dtype = np.int16)
	
	for accession in acc_kmers:
		acc_index = accession_index[accession]
		sql_friendly_accession = accession.replace(".", "_")
		if acc_index in database.accessions:
			#The accession was found for this target genome, for each tgt genome.
			shared_acc_counts[np.nonzero(target_kmer_cts[acc_index])] += 1
			
			these_kmers = [(int(kmer),) for kmer in acc_kmers[accession]]
			#Miguel issue redev.
			temp_name = name+"_"+sql_friendly_accession
			temp_name = temp_name.replace(".", "_")
			
			temp_tab = "CREATE TEMP TABLE " + temp_name + " (kmer INTEGER)"
			database.cursor.execute(temp_tab)
			database.connection.commit()
			insert_table = "INSERT INTO " + temp_name + " VALUES (?)"
			database.cursor.executemany(insert_table, these_kmers)
			database.connection.commit()
			join_and_select_sql = "SELECT genomes FROM " + temp_name + " INNER JOIN " + sql_friendly_accession + " ON "+ temp_name+".kmer = " + sql_friendly_accession+".kmer;"
			
			these_intersections = np.zeros(target_len, dtype = np.int16)
			#sql_query = "SELECT genomes FROM " + sql_friendly_accession + " WHERE kmer in ({kmers})".format(kmers=','.join(['?']*len(these_kmers)))
			for result in database.cursor.execute(join_and_select_sql).fetchall():
				these_intersections[result] += 1
			
			results += np.divide(these_intersections, np.subtract(np.add(acc_kmers[accession].shape[0], target_kmer_cts[acc_index]), these_intersections))
	
	database.close_connection()
	
	#These are the jacc averages
	jaccard_averages = np.divide(results, shared_acc_counts)
	del results

	aai_est = numpy_kaai_to_aai(jaccard_averages)
	
	no_hit = np.where(shared_acc_counts == 0)
	
	possible_hits = np.minimum(len(acc_kmers), target_protein_counts).astype(str)
	
	jaccard_averages = np.round(jaccard_averages, 4).astype(str)
	
	shared_acc_counts = shared_acc_counts.astype(str)
		
	jaccard_averages[no_hit] = "N/A"
	aai_est[no_hit] = "N/A"
	shared_acc_counts[no_hit] = "N/A"
	possible_hits[no_hit] = "N/A"
	
	output_file = temp_out +"/"+name+"_results.txt"
	fh = open(output_file, "w")

	for target in range(0, target_len):
		target_name = database.reverse_genome_index[target]
		if target_name == name:
			fh.write(name+"\t"+target_name+"\t"+"100.0"+"\t"+"0.0"+"\t"+shared_acc_counts[target]+"\t"+possible_hits[target]+"\t"+"100.0"+"\n")
		else:
			fh.write(name+"\t"+target_name+"\t"+jaccard_averages[target]+"\t"+"N/A"+"\t"+shared_acc_counts[target]+"\t"+possible_hits[target]+"\t"+aai_est[target]+"\n")
	
	fh.close()
	
	return output_file

def numpy_kaai_to_aai(kaai_array):
	#aai_hat = (-0.3087057 + 1.810741 * (np.exp(-(-0.2607023 * np.log(kaai))**(1/3.435))))*100
	
	#Protect the original jaccard averages memory item
	aai_hat_array = kaai_array.copy()
	
	non_zero = np.where(aai_hat_array > 0)
	is_zero  = np.where(aai_hat_array <= 0)
	
	#I broke this down into its original components
	#Avoid zeroes in log - still actually works, but it produces warnings I don't want to see.
	aai_hat_array[non_zero] = np.log(aai_hat_array[non_zero])
	
	aai_hat_array = np.multiply(np.subtract(np.multiply(np.exp(np.negative(np.power(np.multiply(aai_hat_array, -0.2607023), (1/3.435)))), 1.810741), 0.3087057), 100)
	'''
	Same as the above, broken down into easier-to-follow steps.
	aai_hat_array = np.multiply(aai_hat_array, -0.2607023)
	aai_hat_array = np.power(aai_hat_array, (1/3.435))
	aai_hat_array = np.negative(aai_hat_array)
	aai_hat_array = np.exp(aai_hat_array)
	aai_hat_array = np.multiply(aai_hat_array, 1.810741)
	aai_hat_array = np.subtract(aai_hat_array, 0.3087057)
	aai_hat_array = np.multiply(aai_hat_array, 100)
	'''
	
	#<30 and >90 values
	smol = np.where(aai_hat_array < 30)
	big  = np.where(aai_hat_array > 90)
	
	aai_hat_array = np.round(aai_hat_array, 2)
	
	#Convert to final printables
	aai_hat_array = aai_hat_array.astype(str)
	aai_hat_array[smol] = "<30%"
	aai_hat_array[big]  = ">90%"
	#The math of the above ends up with zero values being big, so we fix those.
	aai_hat_array[is_zero] = "<30%"
	
	return aai_hat_array
	
def curtime():
	time_format = "%d/%m/%Y %H:%M:%S"
	timer = datetime.datetime.now()
	time = timer.strftime(time_format)
	return time

#Manages the query process.
def db_query_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module takes two FastAAI databases and searches all of the genomes in the QUERY against all of the genomes in the TARGET
	
	If you have many genomes (more than 1000), it will be faster to create the query database using FastAAI build_db, 
	then search it against an existing target using this module than it is to do the same thing with an SQL query.
	
	If you give the same database as query and target, a special all vs. all search of the genomes in the database will be done.
	''')
	parser.add_argument('-q', '--query',  dest = 'query',  default = None, help =  'Path to the query database. The genomes FROM the query will be searched against the genomes in the target database')
	parser.add_argument('-t', '--target', dest = 'target', default = None, help =  'Path to the target database.')
				
	parser.add_argument('-o', '--output', dest = 'output', default = "FastAAI", help = 'The directory where FastAAI will place the result of this query. By default, a directory named "FastAAI" will be created in the current working directory and results will be placed there.')

	parser.add_argument('--threads',      dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')
	parser.add_argument('--verbose',      dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')

	parser.add_argument('--do_stdev',   dest = "do_stdev", action='store_true',                 help = 'Off by default. Calculate std. deviations on Jaccard indicies. Increases memory usage and runtime slightly. Does NOT change estimated AAI values at all.')
	parser.add_argument('--unlimited_resources', dest = "large_mem",  action = 'store_true', help = 'Off by default. Use a faster algorithm that consumes more RAM. FastAAI cannot calculate std. deviations with this algorithm, so they will automatically be skipped.')
	parser.add_argument('--mem',        dest = "precision", default = "med",                 help = 'One of low/med/high. Medium by default. Save RAM in return for slightly rounded AAI estimates. Only affects FastAAI if you are also using the "--unlimited_resources" flag.')

	args, unknown = parser.parse_known_args()
	
	return parser, args
	
#Control the query process for any DB-first query.
def db_query(query, target, verbose, output, threads, do_stdev, precision, memory_efficient):
	print("")
	
	#Sanity checks.
	if not os.path.exists(target):
		print("Target database not found. Exiting FastAAI")
		sys.exit()
	
	if not os.path.exists(query):
		print("Query database not found. Exiting FastAAI")
		sys.exit()
		
	#status = "exists"
	query_ok = assess_db(query)
	target_ok = assess_db(target)
	
	if query_ok != "exists":
		print("Query database improperly formatted. Exiting FastAAI")
		sys.exit()
		
	if target_ok != "exists":
		print("Query database improperly formatted. Exiting FastAAI")
		sys.exit()
	
	#Check if the database is querying against itself.
	if target is None or query is None:
		print("I require both a query and a target database. FastAAI exiting.")
		sys.exit()
	
	if query == target:
		print("Performing an all vs. all query on", query)
		#all_vs_all = True
	else:
		print("Querying", query, "against", target)
		#all_vs_all = False
	
	#Ready the output directories as needed.	
	#The databases are already created, the only state they can be in in P+H
	good_to_go = prepare_directories(output, "protein and HMM", "query")
	if not good_to_go:
		print("Exiting FastAAI")
		sys.exit()
	
	if precision not in ["high", "med", "low"]:
		print("Selected memory usage setting not found. Defaulting to med. Select one with --mem high/med/low.")
		precision = 'med'
	
	#Default
	if (not memory_efficient) or do_stdev:
		do_query_vs_target_sql(query, target, threads, output, verbose, do_stdev)
	#Not default.
	else:
		do_query_vs_target_aai_only(query, target, threads, output, precision, verbose)
	
	print("")
	

#Perform a minimal-memory query of a target database from input files. Lighter weight function for low memory
def sql_query_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module takes one or many genomes, proteins, or proteins and HMMs as a QUERY and searches them against an existing FastAAI database TARGET using SQL
	If you only have a few genomes - or not enough RAM to hold the entire target database in memory - this is the probably the best option for you.
	
	If you provide FastAAI with genomes or only proteins (not proteins and HMMs), this FastAAI module will produce the required protein and HMM files as needed
	and place them in the output directory, just like it does while building a database. 
	
	Once these inputs are ready to be queried against the database (each has both a protein and HMM file), they will be processed independently, 1 per thread at a time.
	
	Note: Protein and HMM files generated during this query can be supplied to build a FastAAI database from proteins and HMMs using the build_db module, without redoing preprocessing.
	''')
			
	parser.add_argument('-g', '--genomes',  dest = 'genomes', default = None,  help =  'A directory containing genomes in FASTA format.')
	parser.add_argument('-p', '--proteins', dest = 'proteins', default = None, help = 'A directory containing protein amino acids in FASTA format.')
	parser.add_argument('-m', '--hmms',     dest = 'hmms', default = None,     help =     'A directory containing the results of an HMM search on a set of proteins.')
	
	parser.add_argument('--target',   dest = 'target', default = None, help =   'A path to the FastAAI database you wish to use as the target')
	
	parser.add_argument('-o', '--output',   dest = 'output', default = "FastAAI", help = 'The directory where FastAAI will place the result of this query and any protein or HMM files it has to generate. By default, a directory named "FastAAI" will be created in the current working directory and results will be placed there.')
	
	parser.add_argument('--genome_file',    dest = 'gf', default = None, help =  'Alternative way to supply genomes. A file containing paths to your genome files, 1 per line.')
	parser.add_argument('--protein_file',   dest = 'pf', default = None, help = 'Alternative way to supply proteins. A file containing paths to your protein files, 1 per line.')
	parser.add_argument('--hmm_file',       dest = 'hf', default = None, help =     'Alternative way to supply HMMs. A file containing paths to your HMM files, 1 per line.')
	
	parser.add_argument('--threads',  dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')
	parser.add_argument('--verbose',        dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
	
	parser.add_argument('--do_stdev',   dest = "do_stdev", action='store_true',   help = 'Off by default. Calculate std. deviations on Jaccard indicies. Increases memory usage and runtime slightly. Does NOT change estimated AAI values at all.')
	
	args, unknown = parser.parse_known_args()
	
	return parser, args

def sql_query_thread_starter(kmer_cts, protein_cts):
	global target_kmer_cts 
	global target_protein_counts
	target_kmer_cts = kmer_cts
	target_protein_counts = protein_cts

	
def sql_query(genomes, proteins, hmms, gf, pf, hf, db_name, output, threads, verbose, do_stdev):
	
	if not os.path.exists(db_name):
		print("")
		print("FastAAI can't find your database:", db_name)
		print("Are you sure that the path you've given to the database is correct and that the database exists?")
		print("FastAAI exiting.")
		print("")
		sys.exit()
		
	start = check_out_input_files(genomes, proteins, hmms, gf, pf, hf)
	
	#If something failed, we stop.
	if start is None:
		sys.exit()

	
	
	good_to_go = prepare_directories(output, start, "query")
	
	if not good_to_go:
		print("Exiting FastAAI")
		sys.exit()
	
	#global kmer_index
	#kmer_index = create_kmer_index()
	
	
	print("")
	print("Preparing inputs for querying...")
	
	prepared_files = advance_inputs(genomes = genomes, proteins = proteins, hmms = hmms, genomes_file = gf, proteins_file = pf, hmms_file = hf, output = output, threads = threads, verbose = verbose, db_name = db_name)
	
	if prepared_files is None:
		return None
	
	query_accessions_detected = set()
	for file in prepared_files:
		query_accessions_detected = query_accessions_detected.union(file.best_hits.values())
	
	#We don't want to get more than we have to.
	query_accessions_detected = list(query_accessions_detected)
	
	if prepared_files is None:
		print("Exiting FastAAI")
		sys.exit()
	
	if verbose:
		print("")	
	print("Gathering database information...")
			
	database = fastaai_database(db_name)
	database.activate_connection()
	database.load_genome_index()
	database.load_accessions()
	database.close_connection()
	
	#formatted_dataset = [(database, file.basename, file.best_hits_kmers, os.path.normpath(output+"/results")) for file in prepared_files] 
	
	#global accession_index
	accession_index = generate_accessions_index()
	
	#Translate to indicies.
	query_accessions_detected = [accession_index[a] for a in query_accessions_detected]
	
	#global target_kmer_cts
	target_kmer_cts = {}
	
	for accession in np.intersect1d(database.accessions, query_accessions_detected):
		target_kmer_cts[accession] = np.zeros(len(database.genome_index), dtype = np.int16)
		for g in database.gak:
			if accession in database.gak[g]:
				target_kmer_cts[accession][g] = database.gak[g][accession]
	
	#global target_protein_counts
	target_protein_counts = np.zeros(len(database.gak), dtype = np.int16)
	for g in database.gak:
		target_protein_counts[g] = len(database.gak[g])
	
	database.gak = None
	
	formatted_dataset = [(database, file.basename, file.best_hits_kmers, os.path.normpath(output+"/results")) for file in prepared_files] 
	
	if verbose:
		print("")
		print("-"*100)
		print("")
	
	count = 0
	total = len(formatted_dataset)
	
	print("Beginning AAI calculation")
	
	#globals to pass... target_kmer_cts target_protein_counts
	#Just remake these in the procs. kmer_index accession_index
	
	if verbose:
		print("")
	#progress bar - possible dangerous use of the return to line start sequence.
		try:
			percentage = 0
			sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% (Query genome ' + str(count) + " of " + str(total) + ' done at '+curtime()+' )\n')
			sys.stdout.flush()
			last_pct = 0
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass
	
	#If parallelized, do parallel

	pool = multiprocessing.Pool(threads, initializer = sql_query_thread_starter, initargs = (target_kmer_cts, target_protein_counts,))
		
	#Process as we go.
	if do_stdev:
		for file in pool.imap(do_sql_query, formatted_dataset):
			
			'''
			handle = open(file, "r")
			
			for line in handle:
				final_result.write(line)
				
			handle.close()
			os.remove(file)	
			'''
			if verbose:
			#progress bar - possible dangerous use of the return to line start sequence.
				try:
					count += 1
					percentage = (count/total)*100
					if int(percentage/2) > last_pct or count == total:
						sys.stdout.write('\033[A')
						sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% (Query genome ' + str(count) + " of " + str(total) + ' done at '+curtime()+' )\n')
						sys.stdout.flush()
					last_pct = int(percentage/2)
				except:
					#It's not really a big deal if the progress bar cannot be printed.
					pass
		
		pool.close()
		pool.join()
	else:

		for file in pool.imap(do_sql_query_no_SD, formatted_dataset):
			'''
			handle = open(file, "r")
			
			for line in handle:
				final_result.write(line)
				
			handle.close()
			os.remove(file)	
			'''
			if verbose:
			#progress bar - possible dangerous use of the return to line start sequence.
				try:
					count += 1
					percentage = (count/total)*100
					if int(percentage/2) > last_pct or count == total:
						sys.stdout.write('\033[A')
						sys.stdout.flush()
						sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% (Query genome ' + str(count) + " of " + str(total) + ' done at '+curtime()+' )\n')
						sys.stdout.flush
					last_pct = int(percentage/2)
				except:
					#It's not really a big deal if the progress bar cannot be printed.
					pass
		
		pool.close()
		pool.join()
	
	if verbose:
		print("")
		print("-"*100)
		print("")
		
	if os.path.exists(output+"/temp"):
		os.rmdir(output+"/temp")
		
	print("FastAAI query complete! Results at:", os.path.normpath(output + "/results"))
	return None
		

#Check to see if the file exists and is a valid fastAAI db
def assess_db(path):
	status = None
	if os.path.exists(path):
		db = fastaai_database(path)
		try:
			db.activate_connection()
			sql = "SELECT name FROM sqlite_master WHERE type='table'"
			
			db.cursor.row_factory = lambda cursor, row: row[0]
			tables = db.cursor.execute(sql).fetchall()
			db.cursor.row_factory = None
			
			db.close_connection()
						
			if len(tables) > 2 and "genome_index" in tables and "genome_acc_kmer_counts" in tables:		
				status = "exists"
			else:
				status = "wrong format"
			
		except:
			status = "wrong format"
		
	else:
		try:
			db = fastaai_database(path)
			db.activate_connection()
			db.initialize_parent_database()
			db.close_connection()
			status = "created"
		except:
			status = "unable to create"
			
	return status
	
#Add one FastAAI DB to another FastAAI DB
def merge_db_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module allows you to add the contents of one or more FastAAI databases to another.
	You must have at least two already-created FastAAI databases using the build_db module before this module can be used.
	
	Supply a comma-separated list of at least one donor database and a single recipient database. 
	If the recipient already exists, then genomes in all the donors will be added to the recipient.
	If the recipient does not already exist, a new database will be created, and the contents of all the donors will be added to it.
	
	Example:
	FastAAI.py merge_db --donors databases/db1.db,databases/db2.db -recipient databases/db3.db --threads 3
	This command will create a new database called "db3.db", merge the data in db1.db and db2.db, and then add the merged data into db3.db
	
	Only the recipient database will be modified; the donors will be left exactly as they were before running this module.
	''')
			
	parser.add_argument('-d', '--donors',      dest = 'donors',    default = None,   help =  'Comma-separated string of paths to one or more donor databases. The genomes FROM the donors will be added TO the recipient and the donors will be unaltered')

	parser.add_argument('-r', '--recipient',   dest = 'recipient', default = None, help =  'Path to the recipient database. Any genomes FROM the donor database not already in the recipient will be added to this database.')	
		
	parser.add_argument('--verbose',           dest = 'verbose',   action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
	
	parser.add_argument('--threads',      dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')

	args, unknown = parser.parse_known_args()
	
	return parser, args
	
def merge_db_thread_starter(rev_index, per_db_accs):
	global reverse_genome_indicies
	global accs_per_db	
	reverse_genome_indicies = rev_index
	accs_per_db = per_db_accs
	
	
	
def merge_db(recipient, donors, verbose, threads):
	#Prettier on the CLI

	if donors is None or recipient is None:
		print("Either donor or target not given. FastAAI is exiting.")
		return None
	
	print("")
	
	donors = donors.split(",")
	valid_donors = []
	for d in donors:
		if os.path.exists(d):
			if d == recipient:
				print("Donor database", d, "is the same as the recipient. This database will be skipped.")
			else:
				check = assess_db(d)
				if check == "exists":
					if d not in valid_donors:
						valid_donors.append(d)
					else:
						print("It appears that database", d, "was already added to the list of donors. Did you type it twice in the list of donors? Skipping it.")
				else:
					if check == "created":
						print("Donor database", d, "not found! Skipping.")
					else:
						print("Something was wrong with supplied database:", d+". A status check found:", check)
		else:
			print("Donor database", d, "not found! Are you sure the path is correct and this donor exists? This database will be skipped.")
	
	if len(valid_donors) == 0:
		print("None of the supplied donor databases were able to be accessed. FastAAI cannot continue if none of these databases are valid. Exiting.")
		sys.exit()
		
	recip_check = assess_db(recipient)
		
	if recip_check == "created" or recip_check == "exists":
		for donor in valid_donors:
			print("Donor database:", donor, "will be added to recipient database:", recipient)
			
		recipient = fastaai_database(recipient)
	else:
		print("I couldn't find or create the recipient database at", recipient+".", "Does the folder you're trying to place this database in exist, and do you have permission to write files to it? FastAAI exiting.")
		sys.exit()
		
	if recipient is None or len(valid_donors) == 0:
		print("I require both a valid donor and a recipient database. FastAAI exiting.")
		sys.exit()

	donor_dbs = []
	for d in valid_donors:
		donor_dbs.append(fastaai_database(d))
		
	all_accessions = set()
	#global joint_genome_index
	joint_genome_index = {}
	joint_genome_counts = {}
	max_index = 0
	#The idea here is to create a set of arrays whose values span the range of each donor's genomes and translate those into an overall list, in order.
	
	#global reverse_genome_indicies
	reverse_genome_indices = {}
	
	#global accs_per_db
	accs_per_db = {}
		
	#Load recipient data, if any.
	if recip_check == "exists":
		recipient.activate_connection()
		recipient.just_accessions()
		recipient.load_genome_index()
		recipient.close_connection()
		
		all_accessions = all_accessions.union(recipient.accessions)
		accs_per_db[recipient.path] = recipient.accessions
		recipient.accessions = None
		max_index = len(recipient.genome_index)
		
		joint_genome_index = dict(zip(recipient.genome_index.keys(), recipient.genome_index.values()))
		joint_genome_counts = dict(zip(recipient.protein_counts_by_genome.keys(), recipient.protein_counts_by_genome.values()))
		
		#reverse_genome_index = dict(zip(joint_genome_index.values(),joint_genome_index.keys()))
		#So... the keys are the genome indicies of the recip. These... shouldn't need any updates. Only the donors need to match.
		ct = 0
		path = recipient.path
		reverse_genome_indices[path] = []
		for idx in sorted(recipient.genome_index.values()):
			reverse_genome_indices[path].append(idx)
		reverse_genome_indices[path] = np.array(reverse_genome_indices[path], dtype = np.int32)
		recipient.genome_index = None
	
	#Donors should always exist, never be created.
	for d in donor_dbs:
		d.activate_connection()
		d.just_accessions()
		d.load_genome_index()
		d.close_connection()
		accs_per_db[d.path] = d.accessions
		all_accessions = all_accessions.union(d.accessions)
		d.accessions = None
		reverse_genome_indices[d.path] = []
		#Database construction indicates this should always be 0-COUNT
		for g in sorted(d.genome_index.keys()):
			if g not in joint_genome_index:
				reverse_genome_indices[d.path].append(max_index)
				joint_genome_index[g] = max_index
				#Map the counts on.
				joint_genome_counts[max_index] = d.protein_counts_by_genome[d.genome_index[g]]
				#reverse_genome_index[max_index] = g
				max_index += 1
			else:
				reverse_genome_indices[d.path].append(joint_genome_index[g])
		#Make it an array, now
		reverse_genome_indices[d.path] = np.array(reverse_genome_indices[d.path], dtype = np.int32)
		d.genome_index = None
	
	#global accession_index
	accession_index = generate_accessions_index()
	
	#global accession_inverter
	accession_inverter = {}
	for acc in accession_index:
		sql_friendly_accession = acc.replace(".", "_")
		accession_inverter[accession_index[acc]] = sql_friendly_accession
		
	all_accessions = list(all_accessions)
	
	#if not os.path.exists("FastAAI_temp"):
	#	os.mkdir("FastAAI_temp")
	
	print("")
	print("Formatting data to add to database. Started at", curtime())
	
	temp_dir = tempfile.mkdtemp()
	try:
		acc_args = [(acc, donor_dbs, recipient, temp_dir) for acc in all_accessions]
	
		if verbose:
			print("")
			count = 0
			total_counts = len(acc_args)
			try:
				percentage = (count/total_counts)*100
				sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(count) + " of " + str(total_counts) + ' done at ' + curtime() + " )\n")
				sys.stdout.flush()			
			except:
				#It's not really a big deal if the progress bar cannot be printed.
				pass
		
		last_pct = 0
		
		pool = multiprocessing.Pool(threads, initializer=merge_db_thread_starter, initargs = (reverse_genome_indices, accs_per_db,))
		
		quiverfull = []
		for result in pool.imap_unordered(pull_and_merge_accession, acc_args):
			acc = result[0]
			child = result[1]
			#sub_gak = result[2]
			
			quiverfull.append([acc, child])
			#gaks.extend(sub_gak)
			
			if verbose:
				count += 1
				try:
					percentage = (count/total_counts)*100
					log_time = curtime()
					sys.stdout.write('\033[A')
					sys.stdout.flush()
					sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(count) + " of " + str(total_counts) + ' done at ' + curtime() + " )\n")
					sys.stdout.flush()			
				except:
					#It's not really a big deal if the progress bar cannot be printed.
					pass
		
		pool.close()
		pool.join()
		
		print("")
		print("Adding data to final database. Started at", curtime())
		
		if verbose:
			print("")
			
			count = 0
			total_counts = len(acc_args)
			try:
				percentage = (count/total_counts)*100
				sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(count) + " of " + str(total_counts) + ' done at ' + curtime() + " )\n")
				sys.stdout.flush()
			except:
				#It's not really a big deal if the progress bar cannot be printed.
				pass
		
		last_pct = 0
		
		recipient.activate_connection()
		genome_list_update_sql = "INSERT OR REPLACE INTO genome_index VALUES (?, ?, ?)"
		genome_reindex = []
		for g in joint_genome_index:
			genome_reindex.append((g, joint_genome_index[g], joint_genome_counts[joint_genome_index[g]]))
			
		recipient.cursor.executemany(genome_list_update_sql, genome_reindex)
		recipient.connection.commit()
		
		del genome_reindex
		
		for result in quiverfull:
			acc = result[0]
			child = result[1]
			
			recipient.add_child_to_parent(acc, child, genomes_too = True, update_gak = True)
			
			if verbose:
				count += 1
				try:
					percentage = (count/total_counts)*100
					log_time = curtime()
					sys.stdout.write('\033[A')
					sys.stdout.flush()
					sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/2)).ljust(50)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(count) + " of " + str(total_counts) + ' done at ' + curtime() + " )\n")
					sys.stdout.flush()
				except:
					#It's not really a big deal if the progress bar cannot be printed.
					pass
	except:
		#Error
		shutil.rmtree(temp_dir)
	finally:
		#Success
		shutil.rmtree(temp_dir)
		
	print("\nDatabases merged!")
	
	return None
	
def pull_and_merge_accession(args):
	accession_index = generate_accessions_index()
	
	#global accession_inverter
	accession_inverter = {}
	for acc in accession_index:
		sql_friendly_accession = acc.replace(".", "_")
		accession_inverter[accession_index[acc]] = sql_friendly_accession
	
	#joint_genome_index, accession_index, accession_inverter, accs_per_db are global already.
	acc, donor_dbs, recipient, temp = args[0], args[1], args[2], args[3]
	
	acc_name = accession_inverter[acc]
	acc_name_gens = acc_name + "_genomes"
	
	query_sql = "SELECT * FROM " + acc_name
	
	temp_db = fastaai_database(os.path.normpath(temp+"/"+acc_name+".db"))
	temp_db.activate_connection()
	
	create_command = "CREATE TABLE IF NOT EXISTS " + acc_name + " (kmer INTEGER PRIMARY KEY, genomes array)"
	temp_db.cursor.execute(create_command)
	temp_db.connection.commit()
	
	create_command = "CREATE TABLE IF NOT EXISTS " + acc_name + "_genomes (genome INTEGER PRIMARY KEY, kmers array)"
	temp_db.cursor.execute(create_command)
	temp_db.connection.commit()
	
	query_lists = {}
	for db in donor_dbs:
		if acc in accs_per_db[db.path]:
			db.activate_connection()
			
			for result in db.cursor.execute(query_sql).fetchall():
				kmer = result[0]
				genomes = result[1]
				translated_genomes = reverse_genome_indicies[db.path][genomes]
				
				if kmer in query_lists:
					query_lists[kmer] = np.union1d(query_lists[kmer], translated_genomes)
				else:
					query_lists[kmer] = translated_genomes
			
			db.close_connection()
	
	#Recipient is not guaranteed to be in the accs per db - if it was created anew, it wouldn't be.
	if recipient.path in accs_per_db:
		if acc in accs_per_db[recipient.path]:
			recipient.activate_connection()	
			
			for result in recipient.cursor.execute(query_sql).fetchall():
				kmer = result[0]
				genomes = result[1]
				translated_genomes = reverse_genome_indicies[recipient.path][genomes]
				if kmer in query_lists:
					query_lists[kmer] = np.union1d(query_lists[kmer], translated_genomes)
				else:
					query_lists[kmer] = translated_genomes
			
			recipient.close_connection()
	
	#Byte-string these.
	for kmer in query_lists:
		query_lists[kmer] = query_lists[kmer].tobytes()
	
	temp_db.cursor.executemany("INSERT INTO " + acc_name + " VALUES (?,?)", zip(query_lists.keys(), query_lists.values()))
	temp_db.connection.commit()
	
	del query_lists
	
	#Reset. Do genomes
	query_genomes_sql = "SELECT * FROM " + acc_name_gens
	query_lists = {}
	for db in donor_dbs:
		if acc in accs_per_db[db.path]:
			db.activate_connection()
			
			for result in db.cursor.execute(query_genomes_sql).fetchall():
				genome = result[0]
				kmers = result[1]
				translated_genome = int(reverse_genome_indicies[db.path][genome])
				#Each genome gets added only once, no dupes.
				if translated_genome not in query_lists:
					query_lists[translated_genome] = kmers

			db.close_connection()
	
	if recipient.path in accs_per_db:
		if acc in accs_per_db[recipient.path]:
			recipient.activate_connection()	
			
			for result in recipient.cursor.execute(query_genomes_sql).fetchall():
				genome = result[0]
				kmers = result[1]
				translated_genome = int(reverse_genome_indicies[recipient.path][genome])
				#Each genome gets added only once, no dupes.
				if translated_genome not in query_lists:
					query_lists[translated_genome] = kmers
			
			recipient.close_connection()
	
	#Byte-string these.
	#gak = []
	for g in query_lists:
		#gak.append((g, acc, query_lists[g].shape[0]))
		query_lists[g] = query_lists[g].tobytes()
		
	
	temp_db.cursor.executemany("INSERT INTO " + acc_name_gens + " VALUES (?,?)", zip(query_lists.keys(), query_lists.values()))
	temp_db.connection.commit()
	
	temp_db.close_connection()
	
	return [acc_name, temp_db.path]

#Query 1 genome vs. 1 target using Carlos' method - just needs query, target, threads
def single_query_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module takes a single query genome, protein, or protein and HMM pair and a single target genome, protein, or protein and HMM pair as inputs and calculates AAI between the two.
	
	If you supply a genome as either query or target, a protein and HMM file will be made for the genome.
	If you supply a protein as either query or target, an HMM file will be made for it.
	If you supply both an HMM and protein, the search will start right away. You cannot provide only an HMM.
	
	No database will be built, and you cannot query multiple genomes with this module.
	
	If you wish to query multiple genomes against themselves in all vs. all AAI search, use aai_index instead.
	If you wish to query multiple genomes against multiple targets, use multi_query instead.
	''')
	parser.add_argument('-qg', '--query_genome',  dest = 'query_genome',  default = None, help =  'Query genome')
	parser.add_argument('-tg', '--target_genome', dest = 'target_genome', default = None, help =  'Target genome')
	
	parser.add_argument('-qp', '--query_protein',  dest = 'query_protein',  default = None, help =  'Query protein')
	parser.add_argument('-tp', '--target_protein', dest = 'target_protein', default = None, help =  'Target protein')
	
	parser.add_argument('-qh', '--query_hmm',  dest = 'query_hmm',  default = None, help =  'Query HMM')
	parser.add_argument('-th', '--target_hmm', dest = 'target_hmm', default = None, help =  'Target HMM')
	
	parser.add_argument('-o', '--output', dest = 'output', default = "FastAAI", help = 'The directory where FastAAI will place the result of this query. By default, a directory named "FastAAI" will be created in the current working directory and results will be placed there.')	
	
	parser.add_argument('--threads',  dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')
	parser.add_argument('--verbose',        dest = 'verbose', action='store_true', help =   'Print minor updates to console. Major updates are printed regardless.')

	#Alternative file input
	
	args, unknown = parser.parse_known_args()
	
	return parser, args

def do_single_query(input_file):
	input_file.preprocess()
	return input_file
	
def intersect_kmer_lists(pair):
	intersection = np.intersect1d(pair[0], pair[1]).shape[0]
	union = pair[0].shape[0] + pair[1].shape[0] - intersection
	return (intersection/union)
	
def kaai_to_aai(kaai):
	# Transform the kAAI into estimated AAI values
	aai_hat = (-0.3087057 + 1.810741 * (np.exp(-(-0.2607023 * np.log(kaai))**(1/3.435))))*100
	
	return aai_hat
	
#This one's unique. It doesn't do anything with the DB, which means it doesn't access any other functionality outside of the input_file class. It just advances a pair of inputs in parallel and does intersections.
def single_query(query_args, target_args, shared_args):
	
	output, threads, verbose = shared_args[0], shared_args[1], shared_args[2]

	genomes, proteins, hmms = query_args[0], query_args[1], query_args[2]
	
	if genomes is None and proteins is None and hmms is None:
		print("Please supply a query genome, protein, or protein and HMM pair.")
		sys.exit()
	
	query = None
	
	if genomes is not None:
		query = input_file(genomes, output, verbose)
		query.set_genome(genomes)
	if proteins is not None:
		if query is not None:
			print("If you supply a genome for either query or target, you must supply ONLY the genome, not a genome and either a protein or HMM.")
			sys.exit()
		else:
			query = input_file(proteins, output, verbose)
			query.set_protein(proteins)
	if hmms is not None:
		if query is None:
			print("If you supply an HMM for either query or target, you must also supply the protein from which the HMM was generated.")
			sys.exit()
		else:
			query.set_hmm(hmms)
		
	genomes, proteins, hmms = target_args[0], target_args[1], target_args[2]

	if genomes is None and proteins is None and hmms is None:
		print("Please supply a target genome, protein, or protein and HMM pair.")
		sys.exit()
	
	target = None
	
	if genomes is not None:
		target = input_file(genomes, output, verbose)
		target.set_genome(genomes)
	if proteins is not None:
		if target is not None:
			print("If you supply a genome for either target or target, you must supply ONLY the genome, not a genome and either a protein or HMM.")
			sys.exit()
		else:
			target = input_file(proteins, output, verbose)
			target.set_protein(proteins)
	if hmms is not None:
		if target is None:
			print("If you supply an HMM for either target or target, you must also supply the protein from which the HMM was generated.")
			sys.exit()
		else:
			target.set_hmm(hmms)
		
	if query.basename == target.basename:
		print("You've selected the same query and target genome. The AAI is 100%.")
		print("FastAAI exiting.")
		return None
		
	statuses = ["genome", "protein", "protein and hmm"]
	query_stat = statuses.index(query.status)
	target_stat = statuses.index(target.status)
	minimum_status = statuses[min(query_stat, target_stat)]
	
	start_printouts = ["[Genome] Protein Protein+HMM", " Genome [Protein] Protein+HMM", "Genome Protein [Protein+HMM]"]
	
	print("")
	print("Query start: ", start_printouts[query_stat])
	print("Target start:", start_printouts[target_stat])
	print("")
	
	good_to_go = prepare_directories(output, minimum_status, "build")
	
	if not good_to_go:
		print("Exiting FastAAI")
		sys.exit()	
	
	qname = query.basename
	tname = target.basename
	
	name = qname + "_vs_" + tname + ".aai.txt"
	print("Output will be located at", os.path.normpath(output) + "/results/"+name)
	
	#Give the data for kmer indexing to the parallel processes
	global kmer_index
	kmer_index = create_kmer_index()
	
	advance_me = [query, target]
	#All we need to do this.
	pool = multiprocessing.Pool(min(threads, 2))
	
	results = pool.map(do_single_query, advance_me)
	
	pool.close()
	pool.join()
	
	query = results[0]
	target = results[1]
	
	#One of the printouts
	max_poss_prots = max(len(query.best_hits_kmers), len(target.best_hits_kmers))
	
	accs_to_view = set(query.best_hits_kmers.keys()).intersection(set(target.best_hits_kmers.keys()))
	
	seq_pairs = [[query.best_hits_kmers[acc], target.best_hits_kmers[acc]] for acc in accs_to_view]
	
	pool = multiprocessing.Pool(min(threads, len(accs_to_view)))
	
	results = np.array(pool.map(intersect_kmer_lists, seq_pairs))
	
	pool.close()
	pool.join()
	
	jacc_mean = np.mean(results)
	jacc_std = np.std(results)
	actual_prots = len(results)
	aai_est = round(kaai_to_aai(jacc_mean), 2)
	
	if aai_est > 90:
		aai_est = "> 90%"
	else:
		if aai_est < 30:
			aai_est = "< 30%"
	
	output = open(name, "w")
	
	print(qname, tname, round(jacc_mean, 4), round(jacc_std, 4), actual_prots, aai_est, file = output)
	
	output.close()
	
	print("FastAAI single query done! Estimated AAI:", aai_est)
	
def aai_index_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module takes a set of genomes, proteins, or proteins and HMMs, creates a FastAAI database from them, and then executes an all vs. all AAI search of the genomes in the database
	''')
	
	parser.add_argument('-g', '--genomes',  dest = 'genomes', default = None, help =  'A directory containing genomes in FASTA format.')
	parser.add_argument('-p', '--proteins', dest = 'proteins', default = None, help = 'A directory containing protein amino acids in FASTA format.')
	parser.add_argument('-m', '--hmms',     dest = 'hmms', default = None, help =     'A directory containing the results of an HMM search on a set of proteins.')
	
	parser.add_argument('-d', '--database', dest = 'db_name', default = "FastAAI_database.sqlite.db", help =  'The name of the database you wish to create or add to. The database will be created if it doesn\'t already exist and placed in the output directory.')
	
	parser.add_argument('-o', '--output',   dest = 'output', default = "FastAAI", help = 'The directory to place the database and any protein or HMM files FastAAI creates. By default, a directory named "FastAAI" will be created in the current working directory and results will be placed there.')
		
	parser.add_argument('--genome_file',    dest = 'gf', default = None, help =  'Alternative way to supply genomes. A file containing paths to your genome files, 1 per line.')
	parser.add_argument('--protein_file',   dest = 'pf', default = None, help = 'Alternative way to supply proteins. A file containing paths to your protein files, 1 per line.')
	parser.add_argument('--hmm_file',       dest = 'hf', default = None, help =     'Alternative way to supply HMMs. A file containing paths to your HMM files, 1 per line.')
		
	parser.add_argument('--verbose',        dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
	parser.add_argument('--threads',  dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')

	
	parser.add_argument('--do_stdev',   dest = "do_stdev", action='store_true',                 help = 'Off by default. Calculate std. deviations on Jaccard indicies. Increases memory usage and runtime slightly. Does NOT change estimated AAI values at all.')
	parser.add_argument('--unlimited_resources', dest = "large_mem",  action = 'store_true', help = 'Off by default. Use a faster algorithm that consumes more RAM. FastAAI cannot calculate std. deviations with this algorithm, so they will automatically be skipped.')
	parser.add_argument('--mem',        dest = "precision", default = "med",                 help = 'One of low/med/high. Medium by default. Save RAM in return for slightly rounded AAI estimates. Only affects FastAAI if you are also using the "--unlimited_resources" flag.')
	
	args, unknown = parser.parse_known_args()
	
	return parser, args

#Build a DB and query a dataset vs. self
def aai_index(genomes, proteins, hmms, db_name, output, threads, gf, pf, hf, verbose, do_stdev, memory_use, unlimited_resources):
	#run build DB and then db_query with the fresh DB
	success = build_db(genomes, proteins, hmms, db_name, output, threads, gf, pf, hf, verbose)
	if success:
		accessible_name = os.path.normpath(output + "/database/" + db_name)
		db_query(accessible_name, accessible_name, verbose, output, threads, do_stdev, memory_use, unlimited_resources)
	else:
		print("Database could not be built. FastAAI exiting.")
		
	return None
	
#Build 2 DBs and query query DB vs target DB
def multi_query_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module takes a set of query genomes/proteins/proteins+HMMs and a set of target genomes/proteins/proteins+HMMs.
	Two FastAAI databases will be created, one for the query and one for the target, then the query database will have AAI calculated against the target database
	''')

	parser.add_argument('-qg', '--query_genomes',  dest = 'query_genomes', default = None, help =  'A directory containing query genomes in FASTA format.')
	parser.add_argument('-qp', '--query_proteins', dest = 'query_proteins', default = None, help = 'A directory containing query protein amino acids in FASTA format.')
	parser.add_argument('-qm', '--query_hmms',     dest = 'query_hmms', default = None, help =     'A directory containing the results of an HMM search on the set of query proteins.')
	
	parser.add_argument('-tg', '--target_genomes',  dest = 'target_genomes', default = None, help =  'A directory containing target genomes in FASTA format.')
	parser.add_argument('-tp', '--target_proteins', dest = 'target_proteins', default = None, help = 'A directory containing target protein amino acids in FASTA format.')
	parser.add_argument('-tm', '--target_hmms',     dest = 'target_hmms', default = None, help =     'A directory containing the results of an HMM search on the set of target proteins.')
	
	
	parser.add_argument('-qd', '--query_database', dest = 'query_db_name', default =   "FastAAI_query_database.sqlite.db", help =  'The name of the query database you wish to create or add to. The database will be created if it doesn\'t already exist and placed in the output directory.')
	parser.add_argument('-td', '--target_database', dest = 'target_db_name', default = "FastAAI_target_database.sqlite.db", help =  'The name of the target database you wish to create or add to. The database will be created if it doesn\'t already exist and placed in the output directory.')
	
	parser.add_argument('-o', '--output',   dest = 'output', default = "FastAAI", help = 'The directory to place the database and any protein or HMM files FastAAI creates. By default, a directory named "FastAAI" will be created in the current working directory and results will be placed there.')
						
	parser.add_argument('--query_genome_file',    dest = 'qgf', default = None, help =  'Alternative way to supply genomes. A file containing paths to your query genome files, 1 per line.')
	parser.add_argument('--query_protein_file',   dest = 'qpf', default = None, help = 'Alternative way to supply proteins. A file containing paths to your query protein files, 1 per line.')
	parser.add_argument('--query_hmm_file',       dest = 'qhf', default = None, help =     'Alternative way to supply HMMs. A file containing paths to your query HMM files, 1 per line.')
		
	parser.add_argument('--target_genome_file',    dest = 'tgf', default = None, help =  'Alternative way to supply genomes. A file containing paths to your target genome files, 1 per line.')
	parser.add_argument('--target_protein_file',   dest = 'tpf', default = None, help = 'Alternative way to supply proteins. A file containing paths to your target protein files, 1 per line.')
	parser.add_argument('--target_hmm_file',       dest = 'thf', default = None, help =     'Alternative way to supply HMMs. A file containing paths to your target HMM files, 1 per line.')
	
	parser.add_argument('--threads',  dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')
	parser.add_argument('--verbose',        dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')

	parser.add_argument('--do_stdev',   dest = "do_stdev", action='store_true',                 help = 'Off by default. Calculate std. deviations on Jaccard indicies. Increases memory usage and runtime slightly. Does NOT change estimated AAI values at all.')
	parser.add_argument('--unlimited_resources', dest = "large_mem",  action = 'store_true', help = 'Off by default. Use a faster algorithm that consumes more RAM. FastAAI cannot calculate std. deviations with this algorithm, so they will automatically be skipped.')
	parser.add_argument('--mem',        dest = "precision", default = "med",                 help = 'One of low/med/high. Medium by default. Save RAM in return for slightly rounded AAI estimates. Only affects FastAAI if you are also using the "--unlimited_resources" flag.')
	
	args, unknown = parser.parse_known_args()
	
	return parser, args

#Build 2 DBs and query query DB vs target DB
def multi_query(query_arg_list, target_arg_list, shared_args):
	pass
	output, threads, verbose, do_stdev, mem, efficient = shared_args[0], shared_args[1], shared_args[2], shared_args[3], shared_args[4], shared_args[5]
	
	genomes, proteins, hmms, gf, pf, hf, db_name = query_arg_list[0], query_arg_list[1], query_arg_list[2], query_arg_list[3], query_arg_list[4], query_arg_list[5], query_arg_list[6]
	accessible_name_query = os.path.normpath(output + "/database/" + db_name)
	build_db(genomes, proteins, hmms, db_name, output, threads, gf, pf, hf, verbose)
	
	genomes, proteins, hmms, gf, pf, hf, db_name = target_arg_list[0], target_arg_list[1], target_arg_list[2], target_arg_list[3], target_arg_list[4], target_arg_list[5], target_arg_list[6]
	accessible_name_target = os.path.normpath(output + "/database/" + db_name)
	build_db(genomes, proteins, hmms, db_name, output, threads, gf, pf, hf, verbose)
	
	db_query(accessible_name_query, accessible_name_target, verbose, output, threads, do_stdev, mem, efficient)

	
#Manages the query process.
def miga_merge_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module is intended for internal use by the Microbial Genomes Atlas (MiGA).
	
	It preprocesses a single genome and adds it to an existing database with some assumptions to make the process faster.
	''')

	parser.add_argument('--target',   dest = 'target', default = None,   help =   'A path to the FastAAI database you wish to use as the target')

	parser.add_argument('--genome_file',    dest = 'gf', default = None, help =  'A single FASTA format genome file')
	parser.add_argument('--protein_file',   dest = 'pf', default = None, help =  'A single FASTA format AA protein file')
	parser.add_argument('--hmm_file',       dest = 'hf', default = None, help =  'A single FastAAI HMM result file.')
	
	parser.add_argument('--verbose',        dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
		
	args, unknown = parser.parse_known_args()
	
	return parser, args
	
	
def miga_merge(gf, pf, hf, target, verbose):
	#Temp dir for the outputs
	output = tempfile.mkdtemp()
	manager = None
	if gf is not None:
		manager = input_file(gf, output, verbose)
		manager.set_genome(gf)
		
	if pf is not None and manager is not None:
		manager.set_protein(pf)
		
	if pf is not None and manager is None:
		manager = input_file(pf, output, verbose)
		manager.set_protein(pf)
		
	#We don't need to worry about the hmm file having a none manager - if it does, it's already broke.
	if hf is not None and manager is not None:
		manager.set_hmm(hf)
		
	if manager is None:
		print("Something went wrong! I didn't get any files!")
	else:
		global kmer_index
		kmer_index = create_kmer_index()
		manager.preprocess()
		
		acc_index = generate_accessions_index()

		tgt = fastaai_database(target)
		tgt.activate_connection()
		tgt.just_accessions()
		tgt.load_genome_index()
		
		if manager.basename not in tgt.genome_index:
			
			my_index = max(list(tgt.genome_index.values())) + 1
			
			genome_byte = np.array(my_index, dtype = np.int32)
			genome_byte = genome_byte.tobytes()
			
			#for adding to genome acc kmer count table
			gak = []
			for scg in manager.best_hits:
				scg_hit = manager.best_hits[scg]
				sql_friendly = scg_hit.replace('.', '_')
				acc_num = acc_index[scg_hit]
				
				gak.append((my_index, acc_num, manager.protein_kmer_count[scg_hit],))
				#Add to each row
				genome_row = {my_index : manager.best_hits_kmers[scg_hit]}
				
				tgt.add_genomes_first(scg_hit, genome_row)
				
				#Prepare for insertion
				my_kmer_list = []
				for kmer in manager.best_hits_kmers[scg_hit]:
					my_kmer_list.append((int(kmer), genome_byte, genome_byte,))
				
				if acc_num not in tgt.accessions:
					#If we need to, make new table
					tgt.cursor.execute("CREATE TABLE " + sql_friendly + " (kmer INTEGER PRIMARY KEY, genomes array)")
					tgt.connection.commit()
				
				modifier_sql = "INSERT INTO "+sql_friendly + " VALUES (?, ?) ON CONFLICT(kmer) DO UPDATE SET genomes=genomes || (?)"  
				tgt.cursor.executemany(modifier_sql, my_kmer_list)
				tgt.connection.commit()
			
			#And add to the GAK.
			gak_sql = "INSERT OR REPLACE INTO genome_acc_kmer_counts VALUES (?,?,?)"
			tgt.cursor.executemany(gak_sql, gak)
			tgt.connection.commit()
			
			#Update the genome table
			tgt.cursor.execute("INSERT OR REPLACE INTO genome_index VALUES (?, ?, ?)", (manager.basename, my_index, manager.protein_count))
			tgt.connection.commit()
		
		#We open either way, so we close either way.
		tgt.close_connection()
		
		return None
	
'''
Main
'''
def main():
	#The currently supported modules.
	modules = ["build_db", "merge_db", "simple_query", "db_query", "single_query", "aai_index", "multi_query", "miga_merge"]
	
	#Print modules if someone just types FastAAI
	if len(sys.argv) < 2:
		print("")
		print("    Welcome to FastAAI")
		print("")
		print("")
		print("    Please select one of the following modules:")
		print("")
		print("------------------------------------------- Quick Usage Options -------------------------------------------")
		print("")
		print("    single_query |" + " Quickly query ONE query genome against ONE target genome")
		print("    multi_query  |" + " Create a query DB and a target DB, then calculate query vs. target AAI")
		print("    aai_index    |" + " Create a database from multiple genomes and do an all vs. all AAI index of the genomes")
		print("")
		print("-------------------------------------- Database Construction Options --------------------------------------")
		print("")
		print("    build_db     |" + " Create or add to a FastAAI database from genomes, proteins, or proteins and HMMs")
		print("    merge_db     |" + " Add the contents of one FastAAI DB to another")
		print("")
		print("---------------------------------------------- Query Options ----------------------------------------------")
		print("")
		print("    simple_query |" + " Query a genome or protein (one or many) against an existing FastAAI database")
		print("    db_query     |" + " Query the genomes in one FastAAI database against the genomes in another FastAAI database")
		print("")
		print("-----------------------------------------------------------------------------------------------------------")
		print("")
		print("                    To select a module, enter 'FastAAI [module]' into the command line!")
		print("")
		sys.exit()
		
	#This is the module selection
	selection = sys.argv[1]
	
	if selection not in modules:
		print("")
		print("    I couldn't find the module you specified. Please select one of the following modules:")
		print("")
		print("------------------------------------------- Quick Usage Options -------------------------------------------")
		print("")
		print("    single_query |" + " Quickly query ONE query genome against ONE target genome")
		print("    multi_query  |" + " Create a query DB and a target DB, then calculate query vs. target AAI")
		print("    aai_index    |" + " Create a database from multiple genomes and do an all vs. all AAI index of the genomes")
		print("")
		print("-------------------------------------- Database Construction Options --------------------------------------")
		print("")
		print("    build_db     |" + " Create or add to a FastAAI database from genomes, proteins, or proteins and HMMs")
		print("    merge_db     |" + " Add the contents of one FastAAI DB to another")
		print("")
		print("---------------------------------------------- Query Options ----------------------------------------------")
		print("")
		print("    simple_query |" + " Query a genome or protein (one or many) against an existing FastAAI database")
		print("    db_query     |" + " Query the genomes in one FastAAI database against the genomes in another FastAAI database")
		print("")
		print("-----------------------------------------------------------------------------------------------------------")
		print("")
		print("    To select a module, enter 'FastAAI [module]' into the command line!")
		print("")
		sys.exit()
	
#################### Database build or add ########################
	
	if selection == "build_db":
		parser, opts = build_db_opts()
		
		#module name only
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
		
		#Directory based
		genomes, proteins, hmms = opts.genomes, opts.proteins, opts.hmms
		
		#Input list based
		gf, pf, hf = opts.gf, opts.pf, opts.hf
		
		output  = os.path.normpath(opts.output)
	
		threads = opts.threads
		verbose = opts.verbose
		
		#Database handle
		db_name = opts.db_name
		
		
		#genomes, proteins, hmms, db_name, output, threads, gf, pf, hf, verbose
		build_db(genomes, proteins, hmms, db_name, output, threads, gf, pf, hf, verbose)
		
#################### Add two DBs ########################	

	if selection == "merge_db":
		parser, opts = merge_db_opts()
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
			
		recipient = opts.recipient
		donors = opts.donors
		verbose = opts.verbose
		threads = opts.threads
		
		merge_db(recipient, donors, verbose, threads)
		
#################### Query files vs DB ########################	
		
	if selection == "simple_query":
		parser, opts = sql_query_opts()
		
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
			
		#directory based
		genomes, proteins, hmms = opts.genomes, opts.proteins, opts.hmms
		
		#Input list based
		gf, pf, hf = opts.gf, opts.pf, opts.hf
		
		db_name = opts.target
		
		output  = opts.output
		threads = opts.threads
		verbose = opts.verbose
		
		do_stdev = opts.do_stdev
				
		sql_query(genomes, proteins, hmms, gf, pf, hf, db_name, output, threads, verbose, do_stdev)
		
			
#################### Query DB vs DB ###########################		
	if selection == "db_query":
		parser, opts = db_query_opts()
		#module name only
		
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
			
		query = opts.query
		target = opts.target
		verbose = opts.verbose
		
		do_stdev = opts.do_stdev
		#massive = opts.massive
		
		mem = opts.precision
		efficient = opts.large_mem
		
		output = opts.output
		threads = opts.threads
		
		db_query(query, target, verbose, output, threads, do_stdev, mem, efficient)

#################### One-pass functions #######################
	if selection == "single_query":
		parser, opts = single_query_opts()
		#module name only
		
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
			
		shared_opts = []
		output = os.path.normpath(opts.output)
		threads = opts.threads
		verbose = opts.verbose
		
		shared_opts.append(output)
		
		shared_opts.append(threads)
		shared_opts.append(verbose)
		
		query_opts = []
		
		query_genome = opts.query_genome
		query_protein = opts.query_protein
		query_hmm = opts.query_hmm
		
		
		query_opts.append(query_genome)
		query_opts.append(query_protein)
		query_opts.append(query_hmm)
		
		target_opts = []
		
		target_genome = opts.target_genome
		target_protein = opts.target_protein
		target_hmm = opts.target_hmm
		
		#tg = opts.target_genome_file
		#tp = opts.target_protein_file
		#th = opts.target_hmm_file
		
		target_opts.append(target_genome)
		target_opts.append(target_protein)
		target_opts.append(target_hmm)
			
		single_query(query_opts, target_opts, shared_opts)
		
	if selection == "aai_index":
		parser, opts = aai_index_opts()
		#module name only
		
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
			
			
		genomes, proteins, hmms = opts.genomes, opts.proteins, opts.hmms
		#Text file versions of genomes/proteins/hmms
		gf, pf, hf = opts.gf, opts.pf, opts.hf
		
		db_name = opts.db_name
		
		output = opts.output
		threads = opts.threads
		verbose = opts.verbose
		
		do_stdev = opts.do_stdev
		#massive = opts.massive
		
		mem = opts.precision
		efficient = opts.large_mem
		
		aai_index(genomes, proteins, hmms, db_name, output, threads, gf, pf, hf, verbose, do_stdev, mem, efficient)
	
	if selection == "multi_query":
		parser, opts = multi_query_opts()
		#module name only
		
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
			
		shared_arg_list = []
		output = os.path.normpath(opts.output)
		threads = opts.threads
		verbose = opts.verbose
		
		do_stdev = opts.do_stdev
		#massive = opts.massive
		
		mem = opts.precision
		efficient = opts.large_mem
		
		shared_arg_list.append(output)
		shared_arg_list.append(threads)
		shared_arg_list.append(verbose)
		shared_arg_list.append(do_stdev)
		shared_arg_list.append(mem)
		shared_arg_list.append(efficient)

		query_arg_list = []
		genomes, proteins, hmms = opts.query_genomes, opts.query_proteins, opts.query_hmms
		#Text file versions of genomes/proteins/hmms
		gf, pf, hf = opts.qgf, opts.qpf, opts.qhf
		query_db_name = opts.query_db_name
		
		query_arg_list.append(genomes)
		query_arg_list.append(proteins)
		query_arg_list.append(hmms)
		query_arg_list.append(gf)
		query_arg_list.append(pf)
		query_arg_list.append(hf)
		query_arg_list.append(query_db_name)
		
		target_arg_list = []
		genomes, proteins, hmms = opts.target_genomes, opts.target_proteins, opts.target_hmms
		#Text file versions of genomes/proteins/hmms
		gf, pf, hf = opts.tgf, opts.tpf, opts.thf
		target_db_name = opts.target_db_name
		
		target_arg_list.append(genomes)
		target_arg_list.append(proteins)
		target_arg_list.append(hmms)
		target_arg_list.append(gf)
		target_arg_list.append(pf)
		target_arg_list.append(hf)
		target_arg_list.append(target_db_name)

		multi_query(query_arg_list, target_arg_list, shared_arg_list)
		
############## Super duper secret MiGA module #################
	if selection == "miga_merge":
		parser, opts = miga_merge_opts()
		
		#module name only
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
		
		#file based
		gf, pf, hf = opts.gf, opts.pf, opts.hf
		
		target = opts.target
		
		verbose = opts.verbose
		
		miga_merge(gf, pf, hf, target, verbose)
	
	return None
	

if __name__ == "__main__":
	main()

	