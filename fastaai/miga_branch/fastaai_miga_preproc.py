import sys
import os
import pyrodigal as pd
import pyhmmer

import gzip
from collections import namedtuple
import argparse 
import datetime 
import json

import numpy as np

class fasta_file:
	def __init__(self, file):
		self.file_path = os.path.abspath(file)
		
		self.contents = {}
		
		self.read_fasta()
		
	def read_fasta(self):
		cur_seq = ""
		cur_prot = ""
		
		contents = {}
		deflines = {}
		
		fasta = agnostic_reader(self.file_path)
		for line in fasta:
			if line.startswith(">"):
				if len(cur_seq) > 0:
					contents[cur_prot] = cur_seq
					deflines[cur_prot] = defline
					
				cur_seq = ""
				cur_prot = line.strip().split()[0][1:]
				defline = line.strip()[len(cur_prot)+1 :].strip()
				
			else:
				cur_seq += line.strip()
					
		fasta.close()
	
		#Final iter
		if len(cur_seq) > 0:
			contents[cur_prot] = cur_seq
			deflines[cur_prot] = defline
		
		self.contents = contents
		
		#return contents, deflines	
		return None	
		
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

class pyrodigal_manager:
	def __init__(self, sequences = None,
				trans_tables = [11, 4], #translation tables to use - only relevant if training input is None and is_meta is False
				aa = None, compress = False):
		
		
		self.meta = False
		self.gene_predictors = {}
		
		self.sequences = sequences
		self.seqlens = {}

		self.training_seq = [] 
		self.running_sum = 0 #training sequence current length
		self.training_needs_join = False #needed if more than 1 seq added to training sequence
		self.training_done = False
		
		self.trans_tables = trans_tables
		self.training_data = {}
		
		self.aa = aa
		
		self.predicted_genes = {}
		self.coding_densities = {}
		if self.meta:
			self.predicted_genes[-1] = {}
			self.coding_densities[-1] = {}
		else:
			for t in self.trans_tables:
				self.predicted_genes[t] = {}
				self.coding_densities[t] = {}
		
		self.log = []
		self.do_compress = compress
				
	def sequence_handler(self):
		self.training_seq = []
		self.training_needs_join = False
		if not self.meta:
			for seqid in self.sequences:
				current_seqlen = len(self.sequences[seqid])
				self.seqlens[seqid] = current_seqlen #get length
				self.sequences[seqid] = self.sequences[seqid].encode() #to binary
				
				self.training_seq.append(self.sequences[seqid]) #add to training set
				self.running_sum += current_seqlen #running total of 32 million or less
				
				if self.training_needs_join:
					self.running_sum += 12
					
				self.training_needs_join = True
				
				if self.running_sum > 32000000:
					self.train_manager()
			
			if not self.training_done:
				self.train_manager()
			
	def convert_seq(self, this_seqid):
		self.sequences[this_seqid] = ''.join(self.sequences[this_seqid])
		seqlen = len(self.sequences[this_seqid])
		self.sequences[this_seqid] = self.sequences[this_seqid].encode()
		return seqlen
	
	def train_manager(self):
		if not self.training_done: #Make sure I don't happen twice
			self.training_done = True #Make sure I don't happen twice
			if self.running_sum < 20000:
				self.log.append("Can't train on 20 thousand or fewer characters. Switching to meta mode.")
				self.gene_predictor = pd.OrfFinder(meta=True)
				self.is_meta = True	
			else:
				#Collect sequences into a prodigal-formatted string
				self.training_seq = b'TTAATTAATTAA'.join(self.sequences.values())
				
				#Truncate to 32 million bp if needed
				if self.running_sum > 32000000:
					self.log.append("Warning:  Sequence is long (max 32000000 for training).")
					self.log.append("Training on the first 32000000 bases.")
					self.training_seq = self.training_seq[0:32000000]
				
				#G is 71, C is 67; we're counting G + C and dividing by the total.
				gc = round(((self.training_seq.count(67) + self.training_seq.count(71))/ len(self.training_seq)) * 100, 2)
					
				self.log.append(str(len(self.training_seq)) + " bp training seq created, " + str(gc) + " pct GC")
				
				#Intialize orffinder
				self.gene_predictor = pd.OrfFinder(meta=False)

				#Create training data on each sequence
				for i in range(0, len(self.trans_tables)):
					next_table = self.trans_tables.pop(0)
					
					self.gene_predictors[next_table] = pd.OrfFinder(meta = False)
					self.gene_predictors[next_table].train(self.training_seq, translation_table = next_table)
					
					#Making another OrfFinder instance with this will allow quick swapping while comparing tables.
					self.training_data[next_table] = self.gene_predictors[next_table].training_info
							
				#Clean up afterwards
				self.training_seq = None
	
	def predict(self):
		#Eliminate sequence entries to prevent memory bloat.
		#Usually just grabs one sequence.
		remaining_sequence_ids = tuple(self.sequences.keys())
		for seqid in remaining_sequence_ids:
			sequence = self.sequences.pop(seqid)
			
			for tt in self.gene_predictors:
				#How do we get this working with the training data instances...
				next_genes = self.gene_predictors[tt].find_genes(sequence)
				
				#Keep internal copy if the goal is to reuse them in another program
				self.predicted_genes[tt][seqid] = next_genes #Easier to retain like this and call gene functions.
						
		self.compare_predicted_genes()
		
	def compare_predicted_genes(self):
		if len(self.predicted_genes) == 1:
			pass
		else:
			for tt in self.predicted_genes:
				total_seqlen = 0
				total_coding_bases = 0
				for seqid in self.predicted_genes[tt]:
					seqlen = self.seqlens[seqid]
					total_seqlen += seqlen
					for gene in self.predicted_genes[tt][seqid]:
						total_coding_bases += (gene.end - gene.begin + 1) #Sequence is 1 longer because it's inclusive
									
				self.coding_densities[tt] = total_coding_bases/total_seqlen
			
			tables_to_remove = list(self.coding_densities.keys())
			winning_table = None
			winning_density = 0
			for tt in self.coding_densities:
				if self.coding_densities[tt] > 1.1 * winning_density:
					winning_density = self.coding_densities[tt]
					winning_table = tt
				
			tables_to_remove.pop(tables_to_remove.index(winning_table)) #keep the winning table by removing all others
			
			
			self.log.append("Winning translation table was: " + str(winning_table) + " with coding density " +  str(round(winning_density, 4)))
			for t in tables_to_remove:
				self.log.append("Losing translation table: " + str(t) + " had coding density" + str(round(self.coding_densities[t], 4)))
		
			self.predicted_genes = self.predicted_genes[winning_table] #keep the winning set.
	
	def format_seq(self, seq, num_chars = 60):
		#ceiling funciton without the math module
		ceiling = int(round((len(seq)/num_chars)+0.5, 0))
		formatted = '\n'.join([seq[(i*num_chars):(i+1)*num_chars] for i in range(0, ceiling)])
		formatted = formatted.strip()
		
		return formatted
	
	def write_aa_file(self):
		if self.aa is not None:
			content = []
			
			seqnum = 1
			for seqid in self.predicted_genes:
				gene_num = 1
				for g in self.predicted_genes[seqid]:
					#print(g)
					protein_name = ">" + seqid + "_" + str(gene_num)
					#table = g.translation_table
					start = str(g.begin)
					end = str(g.end)
					strand = str(g.strand)
					annotation = g._gene_data(seqnum)
					translation = g.translate()
					writeable_trans = self.format_seq(translation)
					translation = None
					
					header = " # ".join([protein_name, start, end, strand, annotation])
					
					content.append(header)
					content.append(writeable_trans)
					
					gene_num += 1
				
				seqnum += 1
			
			content = "\n".join(content)
			content += "\n" #final newline
			
			if self.do_compress:
				if not self.aa.endswith(".gz"):
					self.aa += ".gz"
					
				content = content.encode()
					
				output_writer = gzip.open(self.aa, "wb")
			else:
				output_writer = open(self.aa, "w")
				
			output_writer.write(content)
			
			output_writer.close()
			
			content = None
					
	def convert_to_internal_rep(self): #go from pyrodigal objects to protein name:translation dict
		conversion = {}
		for seqid in self.predicted_genes:
			gene_num = 1
			for g in self.predicted_genes[seqid]:
				#print(g)
				protein_name = seqid + "_" + str(gene_num)
				translation = g.translate()
				conversion[protein_name] = translation
				gene_num += 1
			
		self.predicted_genes = conversion
		conversion = None
	
	def run(self):
		self.sequence_handler()
		self.predict()
		self.write_aa_file()
		self.convert_to_internal_rep()

class pyhmmer_manager:
	def __init__(self, do_compress):
		self.hmm_model = []
		
		self.proteins_to_search = []
		self.protein_descriptions = None
		
		self.hmm_result_proteins = []
		self.hmm_result_accessions = []
		self.hmm_result_scores = []
		
		self.printable_lines = []
		
		self.bacterial_SCPs = None
		self.archaeal_SCPs = None
		self.assign_hmm_sets()
		self.domain_counts = {"Bacteria" : 0, "Archaea": 0}
		self.voted_domain = {"Bacteria" : len(self.bacterial_SCPs), "Archaea" : len(self.archaeal_SCPs)}
		
		self.bacterial_fraction = None
		self.archaeal_fraction = None
		
		self.best_hits = None
		
		self.do_compress = do_compress
		
	#Load HMM
	def load_hmm_from_file(self, hmm_path):
		hmm_set = pyhmmer.plan7.HMMFile(hmm_path)
		for hmm in hmm_set:
			self.hmm_model.append(hmm)
			
	#Set archaeal and bacterial HMM sets.
	def assign_hmm_sets(self):
		self.bacterial_SCPs = {'PF00709_21': 'Adenylsucc_synt', 'PF00406_22': 'ADK', 'PF01808_18': 'AICARFT_IMPCHas', 'PF00231_19': 'ATP-synt',
		'PF00119_20': 'ATP-synt_A', 'PF01264_21': 'Chorismate_synt', 'PF00889_19': 'EF_TS', 'PF01176_19': 'eIF-1a',
		'PF02601_15': 'Exonuc_VII_L', 'PF01025_19': 'GrpE', 'PF01725_16': 'Ham1p_like', 'PF01715_17': 'IPPT',
		'PF00213_18': 'OSCP', 'PF01195_19': 'Pept_tRNA_hydro', 'PF00162_19': 'PGK', 'PF02033_18': 'RBFA', 'PF02565_15': 'RecO_C',
		'PF00825_18': 'Ribonuclease_P', 'PF00687_21': 'Ribosomal_L1', 'PF00572_18': 'Ribosomal_L13',
		'PF00238_19': 'Ribosomal_L14', 'PF00252_18': 'Ribosomal_L16', 'PF01196_19': 'Ribosomal_L17',
		'PF00861_22': 'Ribosomal_L18p', 'PF01245_20': 'Ribosomal_L19', 'PF00453_18': 'Ribosomal_L20',
		'PF00829_21': 'Ribosomal_L21p', 'PF00237_19': 'Ribosomal_L22', 'PF00276_20': 'Ribosomal_L23',
		'PF17136_4': 'ribosomal_L24', 'PF00189_20': 'Ribosomal_S3_C', 'PF00281_19': 'Ribosomal_L5', 'PF00181_23': 'Ribosomal_L2',
		'PF01016_19': 'Ribosomal_L27', 'PF00828_19': 'Ribosomal_L27A', 'PF00830_19': 'Ribosomal_L28',
		'PF00831_23': 'Ribosomal_L29', 'PF00297_22': 'Ribosomal_L3', 'PF01783_23': 'Ribosomal_L32p',
		'PF01632_19': 'Ribosomal_L35p', 'PF00573_22': 'Ribosomal_L4', 'PF00347_23': 'Ribosomal_L6',
		'PF03948_14': 'Ribosomal_L9_C', 'PF00338_22': 'Ribosomal_S10', 'PF00411_19': 'Ribosomal_S11',
		'PF00416_22': 'Ribosomal_S13', 'PF00312_22': 'Ribosomal_S15', 'PF00886_19': 'Ribosomal_S16',
		'PF00366_20': 'Ribosomal_S17', 'PF00203_21': 'Ribosomal_S19', 'PF00318_20': 'Ribosomal_S2',
		'PF01649_18': 'Ribosomal_S20p', 'PF01250_17': 'Ribosomal_S6', 'PF00177_21': 'Ribosomal_S7',
		'PF00410_19': 'Ribosomal_S8', 'PF00380_19': 'Ribosomal_S9', 'PF00164_25': 'Ribosom_S12_S23',
		'PF01193_24': 'RNA_pol_L', 'PF01192_22': 'RNA_pol_Rpb6', 'PF01765_19': 'RRF', 'PF02410_15': 'RsfS',
		'PF03652_15': 'RuvX', 'PF00584_20': 'SecE', 'PF03840_14': 'SecG', 'PF00344_20': 'SecY', 'PF01668_18': 'SmpB',
		'PF00750_19': 'tRNA-synt_1d', 'PF01746_21': 'tRNA_m1G_MT', 'PF02367_17': 'TsaE', 'PF02130_17': 'UPF0054',
		'PF02699_15': 'YajC'}
		
		self.archaeal_SCPs = {'PF00709_21': 'Adenylsucc_synt', 'PF05221_17': 'AdoHcyase', 'PF01951_16': 'Archease', 'PF01813_17': 'ATP-synt_D',
		'PF01990_17': 'ATP-synt_F', 'PF01864_17': 'CarS-like', 'PF01982_16': 'CTP-dep_RFKase', 'PF01866_17': 'Diphthamide_syn',
		'PF04104_14': 'DNA_primase_lrg', 'PF01984_20': 'dsDNA_bind', 'PF04010_13': 'DUF357', 'PF04019_12': 'DUF359',
		'PF04919_12': 'DUF655', 'PF01912_18': 'eIF-6', 'PF05833_11': 'FbpA', 'PF01725_16': 'Ham1p_like',
		'PF00368_18': 'HMG-CoA_red', 'PF00334_19': 'NDK', 'PF02006_16': 'PPS_PS', 'PF02996_17': 'Prefoldin',
		'PF01981_16': 'PTH2', 'PF01948_18': 'PyrI', 'PF00687_21': 'Ribosomal_L1', 'PF00572_18': 'Ribosomal_L13',
		'PF00238_19': 'Ribosomal_L14', 'PF00827_17': 'Ribosomal_L15e', 'PF00252_18': 'Ribosomal_L16',
		'PF01157_18': 'Ribosomal_L21e', 'PF00237_19': 'Ribosomal_L22', 'PF00276_20': 'Ribosomal_L23',
		'PF16906_5': 'Ribosomal_L26', 'PF00831_23': 'Ribosomal_L29', 'PF00297_22': 'Ribosomal_L3',
		'PF01198_19': 'Ribosomal_L31e', 'PF01655_18': 'Ribosomal_L32e', 'PF01780_19': 'Ribosomal_L37ae',
		'PF00832_20': 'Ribosomal_L39', 'PF00573_22': 'Ribosomal_L4', 'PF00935_19': 'Ribosomal_L44', 'PF17144_4': 'Ribosomal_L5e',
		'PF00347_23': 'Ribosomal_L6', 'PF00411_19': 'Ribosomal_S11', 'PF00416_22': 'Ribosomal_S13',
		'PF00312_22': 'Ribosomal_S15', 'PF00366_20': 'Ribosomal_S17', 'PF00833_18': 'Ribosomal_S17e',
		'PF00203_21': 'Ribosomal_S19', 'PF01090_19': 'Ribosomal_S19e', 'PF00318_20': 'Ribosomal_S2',
		'PF01282_19': 'Ribosomal_S24e', 'PF01667_17': 'Ribosomal_S27e', 'PF01200_18': 'Ribosomal_S28e',
		'PF01015_18': 'Ribosomal_S3Ae', 'PF00177_21': 'Ribosomal_S7', 'PF00410_19': 'Ribosomal_S8',
		'PF01201_22': 'Ribosomal_S8e', 'PF00380_19': 'Ribosomal_S9', 'PF00164_25': 'Ribosom_S12_S23',
		'PF06026_14': 'Rib_5-P_isom_A', 'PF01351_18': 'RNase_HII', 'PF13656_6': 'RNA_pol_L_2',
		'PF01194_17': 'RNA_pol_N', 'PF03874_16': 'RNA_pol_Rpb4', 'PF01192_22': 'RNA_pol_Rpb6',
		'PF01139_17': 'RtcB', 'PF00344_20': 'SecY', 'PF06093_13': 'Spt4', 'PF00121_18': 'TIM', 'PF01994_16': 'Trm56',
		'PF00749_21': 'tRNA-synt_1c', 'PF00750_19': 'tRNA-synt_1d', 'PF13393_6': 'tRNA-synt_His',
		'PF01142_18': 'TruD', 'PF01992_16': 'vATP-synt_AC39', 'PF01991_18': 'vATP-synt_E', 'PF01496_19': 'V_ATPase_I'}
		
	#Convert passed sequences.
	def convert_protein_seqs_in_mem(self, contents):
		#Clean up.
		self.proteins_to_search = []
		
		for protein in contents:
			#Skip a protein if it's longer than 100k AA.
			if len(contents[protein]) >= 100000:
				continue
			as_bytes = protein.encode()
			#Pyhmmer digitization of sequences for searching.
			easel_seq = pyhmmer.easel.TextSequence(name = as_bytes, sequence = contents[protein])
			easel_seq = easel_seq.digitize(pyhmmer.easel.Alphabet.amino())
			self.proteins_to_search.append(easel_seq)
			
		easel_seq = None		
			
	def execute_search(self):
		top_hits = list(pyhmmer.hmmsearch(self.hmm_model, self.proteins_to_search, cpus=1, bit_cutoffs="trusted"))

		self.printable_lines = []
		
		self.hmm_result_proteins = []
		self.hmm_result_accessions = []
		self.hmm_result_scores = []
		
		for model in top_hits:
			for hit in model:
				target_name = hit.name.decode()
				target_acc = hit.accession
				if target_acc is None:
					target_acc = "-"
				else:
					target_acc = target_acc.decode()
				
				query_name = hit.best_domain.alignment.hmm_name.decode()
				query_acc = hit.best_domain.alignment.hmm_accession.decode()
				
				full_seq_evalue = "%.2g" % hit.evalue
				full_seq_score = round(hit.score, 1)
				full_seq_bias = round(hit.bias, 1)
				
				best_dom_evalue = "%.2g" % hit.best_domain.alignment.domain.i_evalue
				best_dom_score = round(hit.best_domain.alignment.domain.score, 1)
				best_dom_bias = round(hit.best_domain.alignment.domain.bias, 1)

				#I don't know how to get most of these values.
				exp = 0
				reg = 0
				clu = 0
				ov  = 0
				env = 0
				dom = len(hit.domains)
				rep = 0
				inc = 0
				
				try:
					description = self.protein_descriptions[target_name]
				except:
					description = ""
				
				writeout = [target_name, target_acc, query_name, query_acc, full_seq_evalue, \
				full_seq_score, full_seq_bias, best_dom_evalue, best_dom_score, best_dom_bias, \
				exp, reg, clu, ov, env, dom, rep, inc, description]
				
				#Format and join.
				writeout = [str(i) for i in writeout]
				writeout = '\t'.join(writeout)
				
				self.printable_lines.append(writeout)
				
				self.hmm_result_proteins.append(target_name)
				self.hmm_result_accessions.append(query_acc)
				self.hmm_result_scores.append(best_dom_score)
		
	def filter_to_best_hits(self):
		hmm_file = np.transpose(np.array([self.hmm_result_proteins, self.hmm_result_accessions, self.hmm_result_scores]))
		
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
		
		sql_friendly_names = [i.replace(".", "_") for i in hmm_file[:,1]]
		
		self.best_hits = dict(zip(hmm_file[:,0], sql_friendly_names))
		
		hmm_file = None
	
	#Count per-dom occurs.
	def assign_domain(self):
		for prot in self.best_hits.values():
			if prot in self.bacterial_SCPs:
				self.domain_counts["Bacteria"] += 1
			if prot in self.archaeal_SCPs:
				self.domain_counts["Archaea"] += 1
			
		self.bacterial_fraction = self.domain_counts["Bacteria"] / self.voted_domain["Bacteria"]
		self.aechaeal_fraction = self.domain_counts["Archaea"]  / self.voted_domain["Archaea"]
			
		if self.bacterial_fraction >= self.aechaeal_fraction:
			self.voted_domain = "Bacteria"
		else:
			self.voted_domain = "Archaea"
			
		pop_keys = list(self.best_hits.keys())
		for key in pop_keys:
			if self.voted_domain == "Bacteria":
				if self.best_hits[key] not in self.bacterial_SCPs:
					self.best_hits.pop(key)
			if self.voted_domain == "Archaea":
				if self.best_hits[key] not in self.archaeal_SCPs:
					self.best_hits.pop(key)
				
	def to_hmm_file(self, output):
		if output is not None:
			#PyHMMER data is a bit hard to parse. For each result:
			content = '\n'.join(self.printable_lines) + '\n'
			
			if self.do_compress:
				if not output.endswith(".gz"):
					output += ".gz"
					
				content = content.encode()
							
				fh = gzip.open(output, "wb")
				fh.write(content)
				fh.close()
				content = None
				
			else:
				fh = open(output, "w")
				
				fh.write(content)
					
				fh.close()
				
		content = None
		
	#If we're doing this step at all, we've either loaded the seqs into mem by reading the prot file
	#or have them in mem thanks to pyrodigal.
	def run_for_fastaai(self, prots, hmm_output):
		#self.convert_protein_seqs_in_mem(prots)
		#self.execute_search()
		#self.filter_to_best_hits()
	
		try:
			self.convert_protein_seqs_in_mem(prots)
			self.execute_search()
			self.filter_to_best_hits()
			try:
				self.to_hmm_file(hmm_output)
			except:
				print(output, "cannot be created. HMM search failed. This file will be skipped.")

		except:
			print(hmm_output, "failed to run through HMMER!")
			self.best_hits = None

class mining_straight_down:
	def __init__(self, basename = None, protein_list = None, crystal_output = None):
		self.basename = basename
		self.proteins_to_format = protein_list
		self.output_file = crystal_output
		self.formatted_data = None
		
	#Translate tetramers to unique int32 indices.
	def unique_kmer_simple_key(self, seq):
		#num tetramers = len(seq) - 4 + 1, just make it -3.
		n_kmers = len(seq) - 3
		
		#Converts the characters in a sequence into their ascii int value
		as_ints = np.array([ord(i) for i in seq], dtype = np.int32)
		
		#create seq like 0,1,2,3; 1,2,3,4; 2,3,4,5... for each tetramer that needs a value
		kmers = np.arange(4*n_kmers)
		kmers = kmers % 4 + kmers // 4
		
		#Select the characters (as ints) corresponding to each tetramer all at once and reshape into rows of 4, 
		#each row corresp. to a successive tetramer
		kmers = as_ints[kmers].reshape((n_kmers, 4))
		
		#Given four 2-digit numbers, these multipliers work as offsets so that all digits are preserved in order when summed
		mult = np.array([1000000, 10000, 100, 1], dtype = np.int32)
		
		#the fixed values effectively offset the successive chars of the tetramer by 2 positions each time; 
		#practically, this is concatenation of numbers
		#Matrix mult does this for all values at once.
		return np.unique(np.dot(kmers, mult))

	def prepare_data(self):
		self.formatted_data = {"filename": self.basename, "protein_data":{}}
		for prot_acc_seq in self.proteins_to_format:
			prot = prot_acc_seq[0]
			acc = prot_acc_seq[1]
			kmerized_seq = self.unique_kmer_simple_key(prot_acc_seq[2])
			kmerized_seq = kmerized_seq.tolist()
			#print(kmerized_seq)
			
			self.formatted_data["protein_data"][acc] = {"protein_name":prot, "kmers":kmerized_seq}
			
	def to_json(self):
		with open(self.output_file, "w") as fh:
			json.dump(self.formatted_data, fh, indent = 4)
			
class input_file:
	def __init__(self, genome = None, protein = None, hmm = None, #data inputs
				output_protein = None, output_hmm = None, output_crystal = None, #data outputs
				output_log = None, verbose = False, compress_outputs = False):
		
		self.verbose = verbose
		self.do_compress = compress_outputs
		
		self.genome_input = genome
		self.protein_input = protein
		self.hmm_input = hmm
		
		self.protein_output = output_protein
		self.hmm_output = output_hmm
		self.crystal_output = output_crystal
		
		self.log_contents = []
		self.log_file = output_log
		
		self.initial_status = None
		self.current_status = "genome"
		
		self.basename = None
		
		self.genome = None
		self.proteins = None
		self.hmm_besthits = None
		
		current_datetime = datetime.datetime.now()
		self.timestamps = {"start":current_datetime,
							"protein_pred":current_datetime,
							"hmm_search":current_datetime,
							"crystal":current_datetime}
							
		self.runtimes = None
		
		self.hmm_file = None
	
	def curtime(self, step = None):
		if step is not None:
			self.timestamps[step] = datetime.datetime.now()
	
	def timediffs(self):
		self.runtimes = {}
		protein_pred_time = self.timestamps["protein_pred"] - self.timestamps["start"]
		protein_pred_time = round(protein_pred_time.total_seconds(), 2)
		
		hmm_search_time = self.timestamps["hmm_search"] - self.timestamps["protein_pred"]
		hmm_search_time = round(hmm_search_time.total_seconds(), 2)
		
		crystal_time = self.timestamps["crystal"] - self.timestamps["hmm_search"]
		crystal_time = round(crystal_time.total_seconds(), 2)
		
		self.runtimes["protein_pred"] = protein_pred_time
		self.runtimes["hmm_search"] = hmm_search_time
		self.runtimes["crystal"] = crystal_time
	
	def get_initial_status(self):
		if self.genome_input is not None:
			self.initial_status = "genome"
			
		if self.protein_input is not None:
			self.initial_status = "protein"
		
		if self.hmm_input is not None and self.protein_input is not None:
			self.initial_status = "hmm"
			
	def get_file_basename(self):
		if self.initial_status == "genome":
			self.basename = self.file_basename(self.genome_input)
		if self.initial_status == "protein":
			self.basename = self.file_basename(self.protein_input)
		if self.initial_status == "hmm":
			self.basename = self.file_basename(self.protein_input)
		
	#Not an input sanitizer - simply replaces characters that would throw SQLite for a loop.
	def sql_safe(self, string):
		#Sanitize for SQL
		#These are chars safe for sql
		sql_safe = set('_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
		current_chars = set(string)
		#self.sql_name = self.basename
		#Identify SQL-unsafe characters as those outside the permissible set and replace all with underscores.
		for char in current_chars - sql_safe:
			string = string.replace(char, "_")
			
		return string
		
	#Gonna have to go put this everywhere...
	#Consistent file basename behavior
	def file_basename(self, file):
		#Get the name after the final directory path
		name = os.path.basename(file)
		#Extract the portion of a filename prior to the first '.' separator.
		while name != os.path.splitext(name)[0]:
			name = os.path.splitext(name)[0]
		
		name = self.sql_safe(name)
		
		return name	
		
	def find_hmm(self):
		self.hmm_file = None
		try:
			#Look in the same dir as the script; old method/MiGA friendly
			script_path = os.path.dirname(__file__)
			if len(script_path) == 0:
				script_path = "."
			hmm_complete_model = os.path.abspath(os.path.normpath(script_path +"/"+ "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"))
			self.hmm_file  = str(hmm_complete_model)
		except:
			#Try to locate the data bundled as it would be with a pip/conda install.
			script_path = os.path.dirname(sys.modules['fastAAI_HMM_models'].__file__)
			if len(script_path) == 0:
				script_path = "."
			hmm_complete_model = os.path.abspath(os.path.normpath(script_path + '/00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm'))
			self.hmm_file = str(hmm_complete_model)
			#Check that the file exists or fail to the except.
			fh = open(self.hmm_file)
			fh.close()
	
	#Load existing files functions
	def read_genomes(self):
		if self.genome_input is not None:
			genome_seqs = fasta_file(self.genome_input)
			self.genome = genome_seqs.contents
			genome_seqs = None
			
	def read_proteins(self):
		if self.protein_input is not None:
			#Simple dict of seqid:sequence
			protein_seqs = fasta_file(self.protein_input)
			self.proteins = protein_seqs.contents
			protein_seqs = None
		
	def read_hmms(self):
		if self.hmm_input is not None:
			prots = []
			accs = []
			scores = []
			f = agnostic_reader(self.hmm_input)
			for line in f:
				if line.startswith("#"):
					continue
				else:
					segs = line.strip().split()
					
					if len(segs) < 9:
						continue
					
					prots.append(segs[0])
					accs.append(segs[3])
					scores.append(segs[8])
				
			f.close()
			
			if len(prots) < 1:
				self.best_hits = {}
			
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
			
			sql_friendly_names = [i.replace(".", "_") for i in hmm_file[:,1]]
			self.hmm_besthits = dict(zip(hmm_file[:,0], sql_friendly_names))
	
	#runner functions
	def predict_proteins(self):
		mn = pyrodigal_manager(sequences = self.genome,
								aa = self.protein_output,
								compress = self.do_compress)
		mn.run()
		self.proteins = mn.predicted_genes
		
		mn = None

	def hmm_search_and_BH(self):
		hmm_manager = pyhmmer_manager(self.do_compress)
		hmm_manager.load_hmm_from_file(self.hmm_file)
				
		hmm_manager.run_for_fastaai(prots = self.proteins, hmm_output = self.hmm_output)
		
		self.hmm_besthits = hmm_manager.best_hits		

	def filter_bh_prots(self):
		cleaned_prots = []
		for protein in self.proteins:
			if protein in self.hmm_besthits:
				accession = self.hmm_besthits[protein]
				
				next_item = (protein, accession, self.proteins[protein])
				
				cleaned_prots.append(next_item)
				
		self.proteins = cleaned_prots
		cleaned_prots = None
		
	def crystalize(self):
		mn = mining_straight_down(basename = self.basename, protein_list = self.proteins, crystal_output = self.crystal_output)
		mn.prepare_data()
		mn.to_json()
	
	def run(self):
		self.get_initial_status()
		self.get_file_basename()
		
		self.current_status = self.initial_status
		
		if self.current_status == "genome":
			self.read_genomes()
			self.predict_proteins()
			self.current_status = "protein"
			
		if self.initial_status == "protein":
			self.read_proteins()
		
		if self.verbose:		
			self.curtime("protein_pred")
				
		if self.current_status == "protein":
			self.find_hmm()
			self.hmm_search_and_BH()
			self.current_status = "hmm"
			
		if self.initial_status == "hmm":
			self.read_proteins()
			self.read_hmms()
			
		if self.verbose:
			self.curtime("hmm_search")
			
		if self.current_status == "hmm":
			self.filter_bh_prots()
			self.crystalize()
		
		if self.verbose:
			self.curtime("crystal")
			
		if self.verbose:
			self.timediffs()
			print(self.basename, "complete.")
			print("\tRuntimes: ", self.runtimes)

#Add options	
def options():
	'''
	genome = None, protein = None, hmm = None, #data inputs
	output_protein = None, output_hmm = None, output_crystal = None, #data outputs
	output_log = None, verbose = False, compress_outputs = False
	'''

	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''''')
			
	parser.add_argument('--genome',   dest = 'in_gen', default = None, help = 'Input genome in nt FASTA format.')
	parser.add_argument('--protein',   dest = 'in_prot', default = None, help = 'Input proteome for a genome in AA FASTA format')
	parser.add_argument('--hmm',   dest = 'in_hmm', default = None, help = 'Input FastAAI HMM search result for this proteome. Must be paired with --protein to work.')
	
	parser.add_argument('--output_protein',   dest = 'out_prot', default = None, help = 'An output containing predicted proteins for this genome in AA FASTA format. If omitted, no proteins file will be produced.')
	parser.add_argument('--output_hmm',   dest = 'out_hmm', default = None, help = 'An output containing the results of an HMM search of this proteome against FastAAIs SCPs. If omitted, no HMM file will be produced.')
	parser.add_argument('--output_crystal',   dest = 'out_crystal', default = None, required = True, help = 'Required. A JSON-format output representing the fully preprocessed input.')
	
	parser.add_argument('--compress',   dest = 'compress', action='store_true', help = 'GZIP protein and HMM outputs')
	parser.add_argument('--verbose',   dest = 'verbose', action='store_true', help = 'Print feedback to stdout')

	args, unknown_opts = parser.parse_known_args()
	
	return parser, args

def main():	
	parser, opts = options()
	
	if len(sys.argv) < 3:
		parser.print_help()
	
	ing = opts.in_gen
	inp = opts.in_prot
	inh = opts.in_hmm
	
	outp = opts.out_prot
	outh = opts.out_hmm
	outc = opts.out_crystal
	
	comp = opts.compress
	verb = opts.verbose
	
	mn = input_file(genome = ing, 
					protein = inp,
					hmm = inh,
					output_protein = outp, 
					output_hmm = outh, 
					output_crystal = outc, 
					compress_outputs = comp,
					verbose = verb)
	
	mn.run()
		
if __name__ == "__main__":
	main()
