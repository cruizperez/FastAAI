import sys
import os

import gzip
import argparse 
import json
import gzip

import numpy as np

import multiprocessing
import sqlite3

def convert_array(bytestring):
	return np.frombuffer(bytestring, dtype = np.int32)

sqlite3.register_converter("array", convert_array)

class acc_indexer:
	def __init__(self):
		self.forward = None
		self.reverse = None
		self.generate_accessions_index()

	def list_to_index_dict(self, list):
		result = {}
		counter = 0
		for item in list:
			result[item] = counter
			counter += 1
		return result
		
	def rev_list_to_index_dict(self, list):
		result = {}
		counter = 0
		for item in list:
			result[counter] = item
			counter += 1
		return result
		
	def generate_accessions_index(self):
		acc_list = ['PF01780_19', 'PF03948_14', 'PF17144_4', 'PF00830_19', 'PF00347_23', 'PF16906_5', 'PF13393_6',
			'PF02565_15', 'PF01991_18', 'PF01984_20', 'PF00861_22', 'PF13656_6', 'PF00368_18', 'PF01142_18', 'PF00312_22', 'PF02367_17',
			'PF01951_16', 'PF00749_21', 'PF01655_18', 'PF00318_20', 'PF01813_17', 'PF01649_18', 'PF01025_19', 'PF00380_19', 'PF01282_19',
			'PF01864_17', 'PF01783_23', 'PF01808_18', 'PF01982_16', 'PF01715_17', 'PF00213_18', 'PF00119_20', 'PF00573_22', 'PF01981_16',
			'PF00281_19', 'PF00584_20', 'PF00825_18', 'PF00406_22', 'PF00177_21', 'PF01192_22', 'PF05833_11', 'PF02699_15', 'PF01016_19',
			'PF01765_19', 'PF00453_18', 'PF01193_24', 'PF05221_17', 'PF00231_19', 'PF00416_22', 'PF02033_18', 'PF01668_18', 'PF00886_19',
			'PF00252_18', 'PF00572_18', 'PF00366_20', 'PF04104_14', 'PF04919_12', 'PF01912_18', 'PF00276_20', 'PF00203_21', 'PF00889_19',
			'PF02996_17', 'PF00121_18', 'PF01990_17', 'PF00344_20', 'PF00297_22', 'PF01196_19', 'PF01194_17', 'PF01725_16', 'PF00750_19',
			'PF00338_22', 'PF00238_19', 'PF01200_18', 'PF00162_19', 'PF00181_23', 'PF01866_17', 'PF00709_21', 'PF02006_16', 'PF00164_25',
			'PF00237_19', 'PF01139_17', 'PF01351_18', 'PF04010_13', 'PF06093_13', 'PF00828_19', 'PF02410_15', 'PF01176_19', 'PF02130_17',
			'PF01948_18', 'PF01195_19', 'PF01746_21', 'PF01667_17', 'PF03874_16', 'PF01090_19', 'PF01198_19', 'PF01250_17', 'PF17136_4',
			'PF06026_14', 'PF03652_15', 'PF04019_12', 'PF01201_22', 'PF00832_20', 'PF01264_21', 'PF03840_14', 'PF00831_23', 'PF00189_20',
			'PF02601_15', 'PF01496_19', 'PF00411_19', 'PF00334_19', 'PF00687_21', 'PF01157_18', 'PF01245_20', 'PF01994_16', 'PF01632_19',
			'PF00827_17', 'PF01015_18', 'PF00829_21', 'PF00410_19', 'PF00833_18', 'PF00935_19', 'PF01992_16']
		
		self.forward = list_of_poss_accs = self.list_to_index_dict(acc_list)
		self.reverse = list_of_poss_accs = self.rev_list_to_index_dict(acc_list)
	
class miga_db_builder:
	def __init__(self, path):
		self.path = path
		self.conn = None
		self.curs = None
		
	def activate(self):
		self.conn = sqlite3.connect(self.path)
		self.curs = self.conn.cursor()
		
	def deactivate(self):
		self.curs.close()
		self.conn.close()
		self.conn = None
		self.curs = None	
		
	def initialize_metadata(self):
		self.curs.execute("CREATE TABLE IF NOT EXISTS genome_index (genome text, gen_id integer, protein_count integer)")
		self.curs.execute("CREATE TABLE IF NOT EXISTS genome_acc_kmer_counts (genome integer, accession integer, count integer)")
		self.curs.execute("CREATE INDEX IF NOT EXISTS kmer_acc ON genome_acc_kmer_counts (genome, accession);")
		
	def insert_genome_index(self, gi):
		self.curs.executemany("INSERT INTO genome_index VALUES (?, ?, ?)", gi)
		self.conn.commit()
		
	def insert_gak(self, gak):
		self.curs.executemany("INSERT INTO genome_acc_kmer_counts VALUES (?, ?, ?)", gak)
		self.conn.commit()
		
	def add_acc_genomes(self, acc, data):
		create_sql = "CREATE TABLE IF NOT EXISTS {acc}_genomes (genome INTEGER PRIMARY KEY, kmers array)"
		create_sql = create_sql.format(acc = acc)
		self.curs.execute(create_sql)
		insert_sql = "INSERT INTO {acc}_genomes VALUES (?, ?)"
		insert_sql = insert_sql.format(acc=acc)
		self.curs.executemany(insert_sql, data)
		
		self.conn.commit()
		
		
	def add_acc_kmers(self, acc, data):
		create_sql = "CREATE TABLE IF NOT EXISTS {acc} (kmer INTEGER PRIMARY KEY, genomes array)"
		create_sql = create_sql.format(acc = acc)
		self.curs.execute(create_sql)
		insert_sql = "INSERT INTO {acc} VALUES (?, ?)"
		insert_sql = insert_sql.format(acc=acc)
		self.curs.executemany(insert_sql, data)
		
		self.conn.commit()
		
		self.curs.execute("CREATE INDEX {acc}_index ON {acc} (kmer)".format(acc=acc))
		self.conn.commit()
	
#Class for loading crystal files and prepping them for consumption.
class ravenous_crystal_lizard:
	def __init__(self, crystal_list, database, overwrite = False):
		self.paths_file = crystal_list
		self.input_paths = None
		
		self.crystal_contents = None
		
		self.accession_index = acc_indexer()
		
		self.genome_index = None
		self.genome_prot_ct = None
		self.gak = None
		
		self.db_already_exists = os.path.exists(database)
		self.overwrite = overwrite
		self.og_db_path = database
		self.db = miga_db_builder(database)
		
	def consume_list(self):
		with open(self.paths_file) as fh:
			self.input_paths = fh.readlines()
			
		self.input_paths = [path.strip() for path in self.input_paths]
				
	def consume_crystal_data(self):
		self.crystal_contents = {}
		self.genome_index = []
		self.gak = []
		current_index = 0
		for crystal in self.input_paths:
			if crystal.endswith(".gz"):
				with gzip.open(crystal, "rb") as fh:
					next_crystal = fh.read()
					next_crystal = next_crystal.decode('utf-8')
					next_crystal = json.loads(next_crystal)
			else:
				with open(crystal, "r") as fh:
					next_crystal = json.load(fh)
				
			filename = next_crystal["filename"]
			next_crystal = next_crystal["protein_data"]
			protein_count = len(next_crystal)
			
			next_index = (filename, current_index, protein_count,)
			self.genome_index.append(next_index)
			
			for acc in next_crystal:
				acc_id = self.accession_index.forward[acc]
				
				if acc not in self.crystal_contents:
					self.crystal_contents[acc] = {}
				
				kmer_list = np.array(next_crystal[acc]["kmers"], dtype = np.int32)
				kmer_ct = kmer_list.shape[0]
				
				next_gak = (current_index, acc_id, kmer_ct, )
				self.gak.append(next_gak)
				
				self.crystal_contents[acc][current_index] = kmer_list
				
			current_index += 1
					
		self.db.activate()
		self.db.initialize_metadata()
		self.db.insert_genome_index(self.genome_index)
		self.db.insert_gak(self.gak)
		
		for acc in self.crystal_contents:
			#self.db.add_acc_genomes(acc, self.crystal_contents[acc])
			flipped_dataset = self.invert_to_kmer_first(self.crystal_contents[acc])
			self.db.add_acc_kmers(acc, flipped_dataset)
			flipped_dataset = None
			insertable_genomes = []
			for genome, kmer_array in self.crystal_contents[acc].items():
				next_row = (genome, kmer_array.tobytes(),)
				insertable_genomes.append(next_row)
				
			self.db.add_acc_genomes(acc, insertable_genomes)
			self.crystal_contents[acc] = None
					
		self.db.deactivate()
		
	#Take a set of genome : kmer_lists and flip them to an equivalent set of kmer : genome_lists
	def invert_to_kmer_first(self, dataset):
		genomes = []
		counts = []
		kmer_unlist = []
		for genome_index in dataset:
			genomes.append(genome_index)
			counts.append(dataset[genome_index].shape[0])
			kmer_unlist.append(dataset[genome_index])
		
		genomes = np.array(genomes, dtype = np.int32)
		counts = np.array(counts, dtype = np.int32)
		
		kmer_unlist = np.concatenate(kmer_unlist) #A 1-d array of all of the kmers for all of the genomes containing this SCP
		counted_gens = np.repeat(genomes, counts) #A 1-d array of the same length as kmer_unlist with the corresp. genome index for each kmer
		
		#This contains a list of kmers and genome indices repeated enough times to match their kmer collection 1 to 1 in the same order
		formatted_pairs = np.vstack([kmer_unlist, counted_gens]) 
		kmer_unlist = None
		counted_gens = None
				
		#Sort the list based on kmer, then genome
		sorted_indices = np.lexsort((formatted_pairs[1, :], formatted_pairs[0, :]))
		
		formatted_pairs = formatted_pairs[:, sorted_indices]

		#Collect an ordered list of unique kmers
		discovered_kmers = np.unique(formatted_pairs[0, :])
		
		#Collect a list of the genomes associated with each kmer
		formatted_pairs = np.split(formatted_pairs[1, :], np.unique(formatted_pairs[0, :], return_index=True)[1][1:])

		final_dataset = []
		for kmer, genomes in zip(discovered_kmers, formatted_pairs):
			genome_bytestring = genomes.tobytes()
			
			kmer = int(kmer)
			
			final_dataset.append((kmer, genome_bytestring,))
			
		return final_dataset
		
	def run(self):
		do_run = True
		if self.db_already_exists:
			if self.overwrite:
				os.remove(self.og_db_path)
			else:
				print("")
				print("Target database file already exists! I'm quitting.")
				print("Supply a different path or use --overwrite")
				do_run = False
				
		if do_run:
			self.consume_list()
			self.consume_crystal_data()
		
#Add options	
def options():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''Dedicated MiGA db builder module.
			
	Takes a file containing a list of paths to crystals and builds a database from those.
			
	Notes:
			
	Assumes the supplied database path does not exist. 
	Use --overwrite to delete an existing DB under the same path if you want to.''')
			
	parser.add_argument('--crystal_list',   dest = 'crystals', default = None, help = 'File containing a list of paths to FastAAI crystals')
	parser.add_argument('--database_path',   dest = 'db', default = None, help = 'Path to a NEW database to be built.')
	
	parser.add_argument('--overwrite',   dest = 'overwrite', action = 'store_true', help = 'Delete an existing database at --database_path and create a new one. Otw. quits to preserve existing db.')
		
	args, unknown_opts = parser.parse_known_args()
	
	return parser, args

def main():
	p, a = options()
	crystal_file = a.crystals
	db = a.db
	overwrite = a.overwrite
	
	if len(sys.argv) < 3:
		p.print_help()
	
	if crystal_file is None:
		print("I need a file containing a list of paths to crystals")
		sys.exit()
		
	if db is None:
		print("I need a path to an output database")
		sys.exit()
	
	mn = ravenous_crystal_lizard(crystal_list = crystal_file, 
								database = db, 
								overwrite = overwrite)
	mn.run()
		
if __name__ == "__main__":
	main()
