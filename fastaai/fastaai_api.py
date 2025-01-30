from fastaai import *
import os
import sys
import sqlite3
import numpy as np
import multiprocessing as mp

import datetime

def prepare_genome(genome_file, protein_file = None, hmm_results = None):
	hmm_file = find_hmm()
	hmm_preproc_initializer(hmm_file)
	
	gprep = input_file(genome_file, write_outputs = False)
	gprep.add_triplet(genome_file, protein_file, hmm_results)
	gprep.preprocess()
		
	return gprep

class fastaai_database:
	def __init__(self, db_path):
		self.path = db_path
	
		self.conn = None
		self.curs = None
		
		self.tables_dict = None
		self.genome_index = None
		self.reverse_gix = None
		self.scp_cts = None
		self.gak = None
		
		self.gak_array = None
		self.acc_idx = generate_accessions_index(forward = True)
		
	def open(self, in_mem = False):
		if os.path.exists(self.path):
			self.conn = sqlite3.connect(self.path)
			self.curs = self.conn.cursor()
			try:
				self.curs.execute("SELECT * FROM genome_index").fetchone()
			except:
				print("Error accessing", self.path)
				print("Is it a FastAAI database?")
				self.curs.close()
				self.conn.close()
				self.conn = None
				self.curs = None
				
			if self.curs is not None:
				if in_mem:
					print("Loading database into memory")
					memory_db = sqlite3.connect(":memory:")
					self.curs.close()
					self.conn.backup(memory_db)
					self.conn.close()
					self.conn = memory_db
					self.curs = self.conn.cursor()
					print("Database loaded!")
		else:
			print("No file found at", self.path)
				
	def close(self):
		if self.conn is not None:
			self.curs.close()
			self.conn.close()
			self.conn = None
			self.curs = None
	

	#Access database metadata, table kmer records for in-mem searches before queries are execute
	def prepare_db_for_search(self):
		self.open()
		
		self.tables_dict = None
		self.genome_index = None
		self.scp_cts = None
		self.gak = None
		
		self.genome_index = {}
		self.scp_cts = {}
		
		num_to_acc = generate_accessions_index(forward = False)
		
		table_sql = "SELECT name FROM sqlite_master"
		tables = self.curs.execute(table_sql).fetchall()
		for t in tables:
			actual_name = t[0]
			if actual_name == "genome_index":
				print("Loading genome index")
				res = self.curs.execute("SELECT * FROM genome_index").fetchall()
				for r in res:
					genome = r[0]
					gix = r[1]
					count = r[2]
					self.genome_index[genome] = gix
					self.scp_cts[genome] = count
					
				self.reverse_gix = dict((v,k) for k,v in self.genome_index.items())
				
				
			if actual_name == "genome_acc_kmer_counts":
				print("Loading metadata")
				res = self.curs.execute("SELECT * FROM genome_acc_kmer_counts").fetchall()
				self.gak = {}
				for r in res:
					genome_index = r[0]
					accession_index = r[1]
					acc_name = num_to_acc[accession_index]
					kmer_ct = r[2]
					if acc_name not in self.gak:
						self.gak[acc_name] = {}
					self.gak[acc_name][genome_index] = kmer_ct
				
			if not actual_name.endswith("_genomes") and not actual_name.endswith("_index") and actual_name != "genome_index" and actual_name != "genome_acc_kmer_counts":
				try:
					kmers = self.curs.execute('SELECT kmer FROM "{tab}"'.format(tab=actual_name)).fetchall()
					if len(kmers) > 0:
						if self.tables_dict is None:
							self.tables_dict = {}
						self.tables_dict[actual_name] = []
						for k in kmers:
							self.tables_dict[actual_name].append(k[0])
						self.tables_dict[actual_name] = np.array(self.tables_dict[actual_name], dtype = np.int32)
						self.tables_dict[actual_name] = np.sort(self.tables_dict[actual_name])
				except:
					print("Fail on table name", actual_name)
		
		self.close()
		
		#Ensures this happens after genome index has been populated
		self.gak_array = np.zeros(shape = (122, len(self.genome_index)), dtype = np.int32)
		for acc in self.gak:
			row_id = self.acc_idx[acc]
			for genome in self.gak[acc]:
				self.gak_array[row_id, genome] = self.gak[acc][genome]
	
	def search_prepped_genome_against_database(self, prepped_genome):
		#self.open()
		
		shared_acc_counts = np.zeros(len(self.genome_index), dtype = np.int32)
		jacc_sums = np.zeros(len(self.genome_index), dtype = np.float64)
		
		load_time = 0
		calc_time = 0
		total_kmers_accessed = 0
		total_genome_hits_recovered = 0
		
		for accession_name in prepped_genome.best_hits_kmers:
			start = datetime.datetime.now()
			acc_id = self.acc_idx[accession_name]
			shared_acc_counts[np.where(self.gak_array[acc_id] > 0)] += 1
			shared_kmers = np.intersect1d(prepped_genome.best_hits_kmers[accession_name], self.tables_dict[accession_name])
			query_k_size = shared_kmers.shape[0]
			
			if len(shared_kmers) > 998:
				#Ceiling
				num_grps = int(round((len(shared_kmers) / 998)+0.5))
				shared_k_grps = split_seq(shared_kmers, num_grps)
				
			else:
				shared_k_grps = [shared_kmers.tolist()]
			
			genome_collections = []
			for grp in shared_k_grps:
				fetch_binding = ', '.join(['?']*len(grp))
				next_set = self.curs.execute("SELECT genomes FROM {tab} WHERE kmer IN ({fetch_kmers})".format(tab = accession_name, fetch_kmers = fetch_binding), grp).fetchall()
				
				for r in next_set:
					genome_collections.append(r[0])
			
			end = datetime.datetime.now()
			
			load_time += (end-start).total_seconds()
			
			start = datetime.datetime.now()
			
			genome_collections = b''.join(genome_collections)
			genome_collections = np.frombuffer(genome_collections, dtype = np.int32)
			
			total_kmers_accessed += query_k_size
			total_genome_hits_recovered += genome_collections.shape[0]
			
			genome_collections = np.bincount(genome_collections, minlength = len(self.genome_index))
			union_sizes = self.gak_array[acc_id] + query_k_size - genome_collections
			
			jacc_sums += genome_collections/union_sizes
			
			end = datetime.datetime.now()
			
			calc_time += (end-start).total_seconds()
				
		#self.close()
		
		print(prepped_genome.name, "loaded in", load_time, "calc'd in", calc_time, "kmers queried:", total_kmers_accessed, "genomes_intersections_recovered:", total_genome_hits_recovered)
		
		final_jaccs = jacc_sums / shared_acc_counts
		aai = numpy_kaai_to_aai_just_nums(final_jaccs, as_float=True)
		
		return aai
	
def fastaai_batch_process(genomes_directory, fastaai_database_path, output_file = None, threads = 1, memory_database = False):
	genomes = os.listdir(genomes_directory)
	genomes = [os.path.normpath(genomes_directory+"/"+g) for g in genomes]
	
	genomes = genomes[0:100]
	
	database = fastaai_database(fastaai_database_path)
	database.prepare_db_for_search()
	
	if output_file is not None:
		handle = open(output_file, "w")
		
		print("query_genome", "database_genome", "estimated_AAI", sep = "\t", file = handle)
	else:
		print("query_genome", "database_genome", "estimated_AAI", sep = "\t")
		
	speed_check = []
	
	database.open(in_mem = memory_database)
	pool = multiprocessing.Pool(threads, maxtasksperchild = 1)
	for processed_genome in pool.imap_unordered(prepare_genome, genomes):
		#print("Genes predicted for", processed_genome.name)
		speed_check.append(processed_genome)
		
		aai_scores = database.search_prepped_genome_against_database(processed_genome)
		ct = 0
		for a in aai_scores:
			if output_file is not None:
				print(processed_genome.name, database.reverse_gix[ct], a, sep = "\t", file = handle)
			else:
				print(processed_genome.name, database.reverse_gix[ct], a, sep = "\t")
			ct += 1
		
		
	
	pool.close()
	pool.join()
	
	database.close()
	
	if output_file is not None:
		handle.close()
	
def quick_compare(processed_genome_1 = None, processed_genome_2 = None):
	jacc_avg = None
	if processed_genome_1 is None or processed_genome_2 is None:
		print("Provide two preprocessed genome files prepared with the 'prepare_genome' function.")
	else:
		if processed_genome_1.best_hits_kmers is None or processed_genome_2.best_hits_kmers is None:
			print("The genome files supplied do not appear to be processed.")
		else:
			shared_acc_ct = 0
			jacc_sum = 0
			
			for acc in processed_genome_1.best_hits_kmers:
				if acc in processed_genome_2.best_hits_kmers:
					shared_acc_ct += 1
					intersection_size = np.intersect1d(processed_genome_1.best_hits_kmers[acc], processed_genome_2.best_hits_kmers[acc]).shape[0]
					union_size = processed_genome_1.best_hits_kmers[acc].shape[0] + processed_genome_2.best_hits_kmers[acc].shape[0] - intersection_size
					jacc_sum += (intersection_size/union_size)
					
			jacc_avg = jacc_sum / shared_acc_ct
	
	return jacc_avg
