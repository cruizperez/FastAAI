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
#This as well
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

import pyrodigal as pd
import pyhmmer

from collections import namedtuple

from math import ceil

import re


class progress_tracker:
	def __init__(self, total, step_size = 2, message = None, one_line = True):
		self.current_count = 0
		self.max_count = total
		#Book keeping.
		self.start_time = None
		self.end_time = None
		#Show progrexx every [step] percent
		self.step = step_size
		self.justify_size = ceil(100/self.step)
		self.last_percent = 0
		self.message = message
		
		self.pretty_print = one_line
		
		self.start()

	def curtime(self):
		time_format = "%d/%m/%Y %H:%M:%S"
		timer = datetime.datetime.now()
		time = timer.strftime(time_format)
		return time
		
	def start(self):
		print("")
		if self.message is not None:
			print(self.message)
		
		try:
			percentage = (self.current_count/self.max_count)*100
			sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/self.step)).ljust(self.justify_size)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(self.current_count) + " of " + str(self.max_count) + ' ) at ' + self.curtime() + "\n")
			sys.stdout.flush()
			
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass
	
	def update(self):
		self.current_count += 1
		percentage = (self.current_count/self.max_count)*100
		try:
			if percentage // self.step > self.last_percent:
				if self.pretty_print:
					sys.stdout.write('\033[A')
				sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/self.step)).ljust(self.justify_size)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(self.current_count) + " of " + str(self.max_count) + ' ) at ' + self.curtime() + "\n")
				sys.stdout.flush()
				self.last_percent = percentage // self.step
			#Bar is always full at the end.
			if count == self.max_count:
				if self.pretty_print:
					sys.stdout.write('\033[A')
				sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*self.justify_size).ljust(self.justify_size)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(self.current_count) + " of " + str(self.max_count) + ' ) at ' + self.curtime() + "\n")
				sys.stdout.flush()
				#Add space at end.
				print("")
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass
					
					
#Takes a bytestring from the SQL database and converts it to a numpy array.
def convert_array(bytestring):
	return np.frombuffer(bytestring, dtype = np.int32)

def convert_float_array_16(bytestring):
	return np.frombuffer(bytestring, dtype = np.float16)

def convert_float_array_32(bytestring):
	return np.frombuffer(bytestring, dtype = np.float32)
	
def convert_float_array_64(bytestring):
	return np.frombuffer(bytestring, dtype = np.float64)

def read_fasta(file):
	cur_seq = ""
	cur_prot = ""
	
	contents = {}
	deflines = {}
	
	fasta = agnostic_reader(file)
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
		
	return contents, deflines
	
class fasta_file:
	def __init__(self, file, type = "genome"):
		self.file_path = os.path.abspath(file)
		self.name = os.path.basename(file)
		self.no_ext = os.path.splitext(self.name)[0]
		self.type = type
		
		self.tuple_structure = namedtuple("fasta", ["seqid", "description", "sequence"])
		self.contents = {}
		
	def convert(self, contents, descriptions):
		for protein in contents:
			self.contents = self.tuple_structure(seqid = protein, description = descriptions[protein], sequence = contents[protein])
			
	
	def def_import_file(self):
		contents, descriptions = read_fasta(self.file_path)
		self.convert(contents, descriptions)
		
class pyhmmer_manager:
	def __init__(self, do_compress):
		self.hmm_model = []
		self.hmm_model_optimized = None
		
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
		
		self.do_compress = False
		
	def optimize_models(self):
		try:
			self.hmm_model_optimized = []
			
			for hmm in self.hmm_model:
				prof = pyhmmer.plan7.Profile(M = hmm.insert_emissions.shape[0], alphabet = pyhmmer.easel.Alphabet.amino())
				prof.configure(hmm = hmm, background = pyhmmer.plan7.Background(alphabet = pyhmmer.easel.Alphabet.amino()), L = hmm.insert_emissions.shape[0]-1)
				optim = prof.optimized()
				self.hmm_model_optimized.append(optim)
				
			#Clean up.
			self.hmm_model = None
		except:
			#Quiet fail condition - fall back on default model.
			self.hmm_model_optimized = None
		
	#Load HMM and try to optimize.
	def load_hmm_from_file(self, hmm_path):
		hmm_set = pyhmmer.plan7.HMMFile(hmm_path)
		for hmm in hmm_set:
			self.hmm_model.append(hmm)
			
		#This doesn't seem to be improving performance currently.
		self.optimize_models()
			
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
			
	def load_protein_seqs_from_file(self, prots_file):
		#Pyhmmer has a method for loading a fasta file, but we need to support gzipped inputs, so we do it manually.
		contents, deflines = read_fasta(prots_file)
		self.protein_descriptions = deflines
		self.convert_protein_seqs_in_mem(contents)
		
	def execute_search(self):
		if self.hmm_model_optimized is None:
			top_hits = list(pyhmmer.hmmsearch(self.hmm_model, self.proteins_to_search, cpus=1, bit_cutoffs="trusted"))
		else:
			top_hits = list(pyhmmer.hmmsearch(self.hmm_model_optimized, self.proteins_to_search, cpus=1, bit_cutoffs="trusted"))
		
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
		#PyHMMER data is a bit hard to parse. For each result:
		
		content = '\n'.join(self.printable_lines) + '\n'
		
		if self.do_compress:
			#Clean
			if os.path.exists(output):
				os.remove(output)
			#content = '\n'.join(self.printable_lines) + '\n'
			content = content.encode()
			
			#for formatted_line in self.printable_lines:	
			#	content += self.hmm_result_proteins[i].encode() + b'\t' + self.hmm_result_accessions[i].encode()+ b'\t' + str(round(self.hmm_result_scores[i], 1)).encode() + b'\n'
				
			fh = gzip.open(output+".gz", "wb")
			fh.write(content)
			fh.close()
			content = None
			
		else:
			#Clean
			if os.path.exists(output+".gz"):
				os.remove(output+".gz")
				
			fh = open(output, "w")
			#for formatted_line in self.printable_lines:
			
			fh.write(content)
				
			fh.close()
			
		content = None
		
	#If we're doing this step at all, we've either loaded the seqs into mem by reading the prot file
	#or have them in mem thanks to pyrodigal.
	def run_for_fastaai(self, prots, hmm_output):
		try:
			self.convert_protein_seqs_in_mem(prots)
			self.execute_search()
			self.filter_to_best_hits()
			self.to_hmm_file(hmm_output)
		except:
			print(output, "cannot be created. HMM search failed. This file will be skipped.")
			self.best_hits = None
		
def hmm_preproc_initializer(hmm_file, do_compress = False):
	global hmm_manager
	hmm_manager = pyhmmer_manager(do_compress)
	hmm_manager.load_hmm_from_file(hmm_file)
	
class pyrodigal_manager:
	def __init__(self, file = None, aa_out = None, nt_out = None, is_meta = False, full_headers = True, trans_table = 11,
				num_bp_fmt = True, verbose = True, do_compress = "0", compare_against = None):
		#Input NT sequences
		self.file = file
		
		#List of seqs read from input file.
		self.sequences = None
		#Concatenation of up to first 32 million bp in self.sequences - prodigal caps at this point.
		self.training_seq = None
		
		#Predicted genes go here
		self.predicted_genes = None
		#Record the translation table used.
		self.trans_table = trans_table
		
		#This is the pyrodigal manager - this does the gene predicting.
		self.manager = pd.OrfFinder(meta=is_meta)
		self.is_meta = is_meta
		
		#Full prodigal header information includes more than just a protein number.
		#If full_headers is true, protein deflines will match prodigal; else, just protein ID.
		self.full_headers = full_headers
		
		#Prodigal prints info to console. I enhanced the info and made printing default, but also allow them to be totally turned off.
		self.verbose = verbose
		
		#Prodigal formats outputs with 70 bases per line max
		self.num_bp_fmt = num_bp_fmt
		
		#File names for outputs
		self.aa_out = aa_out
		self.nt_out = nt_out
		
		#List of proteins in excess of 100K base pairs (HMMER's limit) and their lengths. This is also fastAAI specific.
		self.excluded_seqs = {}
		
		#Gzip outputs if asked.
		self.compress = do_compress
		
		self.labeled_proteins = None
		
		#Normally, we don't need to keep an input sequence after it's had proteins predicted for it - however
		#For FastAAI and MiGA's purposes, comparisons of two translation tables is necessary.
		#Rather than re-importing sequences and reconstructing the training sequences, 
		#keep them for faster repredict with less I/O
		self.compare_to = compare_against
		if self.compare_to is not None:
			self.keep_seqs = True
			self.keep_after_train = True
		else:
			self.keep_seqs = False
			self.keep_after_train = False
	
	#Imports a fasta as binary.
	def import_sequences(self):
		if self.sequences is None:
			self.sequences = {}
			
		#check for zipped and import as needed.
		with open(self.file, 'rb') as test_gz:
			#Gzip magic number
			is_gz = (test_gz.read(2) == b'\x1f\x8b')
		
		if is_gz:
			fh = gzip.open(self.file)
		else:
			fh = open(self.file, "rb")
		
		imp = fh.readlines()
		
		fh.close()
		
		cur_seq = None
		for s in imp:
			s = s.decode().strip()
			#> is 62 in ascii. This is asking if the first character is '>'
			if s.startswith(">"):
				#Skip first cycle, then do for each after
				if cur_seq is not None:
					self.sequences[cur_seq] = ''.join(self.sequences[cur_seq])
					self.sequences[cur_seq] = self.sequences[cur_seq].encode()
					#print(cur_seq, len(self.sequences[cur_seq]))
				cur_seq = s[1:]
				cur_seq = cur_seq.split()[0]
				cur_seq = cur_seq.encode('utf-8')
				self.sequences[cur_seq] = []
			else:
				#Remove the newline character.
				#bases = s[:-1]
				self.sequences[cur_seq].append(s)
		
		#Final set
		self.sequences[cur_seq] = ''.join(self.sequences[cur_seq])
		self.sequences[cur_seq] = self.sequences[cur_seq].encode()
		
		#Now we have the data, go to training.
		if not self.is_meta:
			self.train_manager()
		
	#Collect up to the first 32 million bases for use in training seq.
	def train_manager(self):
		running_sum = 0
		seqs_added = 0
		if self.training_seq is None:
			self.training_seq = []
			for seq in self.sequences:
				running_sum += len(self.sequences[seq])
				if seqs_added > 0:
					#Prodigal interleaving logic - add this breaker between sequences, starting at sequence 2
					self.training_seq.append(b'TTAATTAATTAA')
					running_sum += 12
					
				seqs_added += 1
					
				#Handle excessive size
				if running_sum >= 32000000:					
					print("Warning:  Sequence is long (max 32000000 for training).")
					print("Training on the first 32000000 bases.")
				
					to_remove = running_sum - 32000000
					
					#Remove excess characters
					cut_seq = self.sequences[seq][:-to_remove]
					#Add the partial seq
					self.training_seq.append(cut_seq)
					
					#Stop the loop and move to training
					break
				
				#add in a full sequence
				self.training_seq.append(self.sequences[seq])

			if seqs_added > 1:
				self.training_seq.append(b'TTAATTAATTAA')
				
			self.training_seq = b''.join(self.training_seq)
		
		if len(self.training_seq) < 20000:
			if self.verbose:
				print("Can't train on 20 thousand or fewer characters. Switching to meta mode.")
			self.manager = pd.OrfFinder(meta=True)
			self.is_meta = True
		else:
			if self.verbose:
				print("")
				#G is 71, C is 67; we're counting G + C and dividing by the total.
				gc = round(((self.training_seq.count(67) + self.training_seq.count(71))/ len(self.training_seq)) * 100, 2)
				print(len(self.training_seq), "bp seq created,", gc, "pct GC")
				
			#Train
			self.manager.train(self.training_seq, translation_table = self.trans_table)
		
		if not self.keep_after_train:
			#Clean up
			self.training_seq = None
		
	def predict_genes(self):
		if self.is_meta:
			if self.verbose:
				print("Finding genes in metagenomic mode")
		else:
			if self.verbose:
				print("Finding genes with translation table", self.trans_table)
				print("")
			
		self.predicted_genes = {}
		for seq in self.sequences:
			
			if self.verbose:
				print("Finding genes in sequence", seq.decode(), "("+str(len(self.sequences[seq]))+ " bp)... ", end = '')
				
			self.predicted_genes[seq] = self.manager.find_genes(self.sequences[seq])
				
			#If we're comparing multiple tables, then we want to keep these for re-prediction.
			if not self.keep_seqs:
				#Clean up
				self.sequences[seq] = None
			
			if self.verbose:
				print("done!")
			
	#Predict genes with an alternative table, compare results, and keep the winner.	
	def compare_alternative_table(self, table):
		if table == self.trans_table:
			print("You're trying to compare table", table, "with itself.")
		else:
			if self.verbose:
				print("Comparing translation table", self.trans_table, "against table", table)
			old_table = self.trans_table
			old_genes = self.predicted_genes
			old_size = 0
			for seq in self.predicted_genes:
				for gene in self.predicted_genes[seq]:
					old_size += (gene.end - gene.begin)
			
			self.trans_table = table
			self.train_manager()
			self.predict_genes()
				
			new_size = 0
			for seq in self.predicted_genes:
				for gene in self.predicted_genes[seq]:
					new_size += (gene.end - gene.begin)
			
			if (old_size / new_size) > 1.1:
				if self.verbose:
					print("Translation table", self.trans_table, "performed better than table", old_table, "and will be used instead.")
			else:
				if self.verbose:
					print("Translation table", self.trans_table, "did not perform significantly better than table", old_table, "and will not be used.")
				self.trans_table = old_table
				self.predicted_genes = old_genes
			
			#cleanup
			old_table = None
			old_genes = None
			old_size = None
			new_size = None
		
	def predict_and_compare(self):
		self.predict_genes()
	
		#Run alt comparisons in gene predict.
		if self.compare_to is not None:
			while len(self.compare_to) > 0:
				try:
					next_table = int(self.compare_to.pop(0))
					
					if len(self.compare_to) == 0:
						#Ready to clean up.
						self.keep_after_train = True
						self.keep_seqs = True
					
					self.compare_alternative_table(next_table)
				except:
					print("Alternative table comparison failed! Skipping.")
		
	#Break lines into size base pairs per line. Prodigal's default for bp is 70, aa is 60.
	def num_bp_line_format(self, string, size = 70):
		#ceiling funciton without the math module
		ceiling = int(round((len(string)/size)+0.5, 0))
		formatted = '\n'.join([string[(i*size):(i+1)*size] for i in range(0, ceiling)])
		return formatted
	
	#Writeouts
	def write_nt(self):
		if self.nt_out is not None:
			if self.verbose:
				print("Writing nucleotide sequences... ")
			if self.compress == '1' or self.compress == '2':
				out_writer = gzip.open(self.nt_out+".gz", "wb")
				
				content = b''
				
				for seq in self.predicted_genes:
					seqname = b">"+ seq + b"_"
					#Gene counter
					count = 1
					for gene in self.predicted_genes[seq]:
						#Full header lines
						if self.full_headers:
							content += b' # '.join([seqname + str(count).encode(), str(gene.begin).encode(), str(gene.end).encode(), str(gene.strand).encode(), gene._gene_data.encode()])
						else:
							#Reduced headers if we don't care.
							content += seqname + str(count).encode()
							
						content += b'\n'
							
						if self.num_bp_fmt:
							#60 bp cap per line
							content += self.num_bp_line_format(gene.sequence(), size = 70).encode()
						else:
							#One-line sequence.
							content += gene.sequence().encode()
							
						content += b'\n'
						count += 1
				
				out_writer.write(content)
				out_writer.close()
			
			if self.compress == '0' or self.compress == '2':
				out_writer = open(self.nt_out, "w")
			
				for seq in self.predicted_genes:
					#Only do this decode once.
					seqname = ">"+ seq.decode() +"_"
					#Gene counter
					count = 1
					
					for gene in self.predicted_genes[seq]:
						#Full header lines
						if self.full_headers:
							#Standard prodigal header
							print(seqname + str(count), gene.begin, gene.end, gene.strand, gene._gene_data, sep = " # ", file = out_writer)
						else:
							#Reduced headers if we don't care.
							print(seqname + str(count), file = out_writer)
							
						if self.num_bp_fmt:
							#60 bp cap per line
							print(self.num_bp_line_format(gene.sequence(), size = 70), file = out_writer)
						else:
							#One-line sequence.
							print(gene.sequence(), file = out_writer)
							
						count += 1
							
				out_writer.close()
		
	def write_aa(self):
		if self.aa_out is not None:
			if self.verbose:
				print("Writing amino acid sequences...")
				
			self.labeled_proteins = {}
			content = ''
			for seq in self.predicted_genes:
				count = 1
				seqname = ">"+ seq.decode() + "_"
				for gene in self.predicted_genes[seq]:
					prot_name = seqname + str(count)
					translation = gene.translate()
					self.labeled_proteins[prot_name[1:]] = translation
					defline = " # ".join([prot_name, str(gene.begin), str(gene.end), str(gene.strand), str(gene._gene_data)])
					content += defline
					content += "\n"
					count += 1
					content += self.num_bp_line_format(translation, size = 60)
					content += "\n"					
				
			if self.compress == '0' or self.compress == '2':
				out_writer = open(self.aa_out, "w")
				out_writer.write(content)
				out_writer.close()
				
			if self.compress == '1' or self.compress == '2':
				content = content.encode()
				out_writer = gzip.open(self.aa_out+".gz", "wb")
				out_writer.write(content)
				out_writer.close()
				
	def run_for_fastaai(self):
		self.verbose = False
		self.import_sequences()
		self.train_manager()
		self.predict_and_compare()
		self.write_aa()
	
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

'''
Class for handling all of the raw genome/protein/protein+HMM file inputs when building a database.

Takes a file or files and processes them from genome -> protein, protein -> hmm, prot+HMM -> kmerized protein best hits as numpy int arrays according to the kmer_index

'''

class input_file:
	def __init__(self, input_path, output = "", verbosity = False, do_compress = False):
		#starting path for the file; irrelevant for protein and hmm, but otherwise useful for keeping track.
		self.path = input_path
		#Output directory starts with this
		self.output = os.path.normpath(output + "/")
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
			
		#Sanitize for SQL
		#These are chars safe for sql
		sql_safe = set('_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
		current_chars = set(self.basename)
		#self.sql_name = self.basename
		#Identify SQL-unsafe characters as those outside the permissible set and replace all with underscores.
		for char in current_chars - sql_safe:
			self.basename = self.basename.replace(char, "_")
			
		#'genome' or 'protein' or 'protein and HMM' 
		self.status = None
		#These will keep track of paths for each stage of file for us.
		self.genome = None
		self.protein = None
		self.hmm = None
		
		self.ran_hmmer = False
		
		#If pyrodigal is run, then the protein sequences are already loaded into memory. 
		#We reuse them in kmer extraction instead of another I/O
		self.prepared_proteins = None
		
		self.intermediate = None
		
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
		
		#Check if the file failed to produce ANY SCP HMM hits.
		self.is_empty = False
		
		self.do_compress = do_compress
	
		self.crystal = None
		
		self.init_time = None
		#default to 0 time.
		self.prot_pred_time = None
		self.hmm_search_time = None
		self.besthits_time = None
	
	def curtime(self):
		time_format = "%d/%m/%Y %H:%M:%S"
		timer = datetime.datetime.now()
		time = timer.strftime(time_format)
		return time
		
	def partial_timings(self):
		protein_pred = self.prot_pred_time-self.init_time
		hmm_search = self.hmm_search_time-self.prot_pred_time
		besthits = self.besthits_time-self.hmm_search_time
		
		protein_pred = protein_pred.total_seconds()
		hmm_search = hmm_search.total_seconds()
		besthits = besthits.total_seconds()
		
		self.prot_pred_time = protein_pred
		self.hmm_search_time = hmm_search
		self.besthits_time = besthits
	
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
		
	def set_crystal(self, path):
		self.status = 'crystal'
		self.crystal = path

	#Runs prodigal, compares translation tables and stores faa files
	def genome_to_protein(self):
		if self.genome is None:
			print(self.name, "wasn't a declared as a genome! I can't make this into a protein!")
		else:		
			protein_output = os.path.normpath(self.output + "/predicted_proteins/" + self.basename + '.faa')
					
			if self.do_compress:
				compress_level = "1"
			else:
				compress_level = "0"
			
			mn = pyrodigal_manager(file = self.genome, aa_out = protein_output, compare_against = [4], do_compress = compress_level)
			mn.run_for_fastaai()
			
			self.trans_table = str(mn.trans_table)
			
			for prot in mn.excluded_seqs:
				self.err_log += "Protein " + prot + " was observed to have >100K amino acids ( " + str(mn.excluded_seqs[prot]) + " AA found ). It will not be included in predicted proteins for this genome;"
			
			self.prepared_proteins = mn.labeled_proteins
			
			del mn
			
			#If there are zipped files leftover and we didn't want them, clean them up.
			if self.do_compress:
				self.set_protein(str(protein_output)+".gz")
				#Clean up unzipped version on reruns
				if os.path.exists(str(protein_output)):
					os.remove(str(protein_output))
			else:
				self.set_protein(str(protein_output))
				#Clean up a zipped version on reruns
				if os.path.exists(str(protein_output)+".gz"):
					os.remove(str(protein_output)+".gz")

			self.prot_pred_time = datetime.datetime.now()
					
	#run hmmsearch on a protein	
	def protein_to_hmm(self):
		if self.protein is None:
			print(self.basename, "wasn't a declared as a protein! I can't make this into an HMM!")
		else:

			folder = os.path.normpath(self.output + "/hmms")
						
			hmm_output = os.path.normpath(folder +"/"+ self.basename + '.hmm')
			
			if self.prepared_proteins is None:
				self.prepared_proteins, deflines = read_fasta(self.protein)

			hmm_manager.run_for_fastaai(self.prepared_proteins, hmm_output)
			
			self.ran_hmmer = True
			
			if self.do_compress:
				self.set_hmm(str(hmm_output)+".gz")
				if os.path.exists(str(hmm_output)):
					os.remove(str(hmm_output))
			else:
				self.set_hmm(str(hmm_output))
				if os.path.exists(str(hmm_output)+".gz"):
					os.remove(str(hmm_output)+".gz")
					
			self.hmm_search_time = datetime.datetime.now()
	
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

	def load_hmm_and_filter_from_file(self):
		prots = []
		accs = []
		scores = []
		f = agnostic_reader(self.hmm)
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
		self.best_hits = dict(zip(hmm_file[:,0], sql_friendly_names))
			
		
	#This should consider the domain by majority vote...
	def prot_and_hmm_to_besthits(self):
		if self.ran_hmmer:
			#Manager has a filter built in.
			self.best_hits = hmm_manager.best_hits
		else:
			#Load the best hits file via old numpy method.
			self.load_hmm_and_filter_from_file()
		
		hit_count = 0
				
		#from pyrodigal predictions or HMM intermediate production, the sequences are already in mem and don't need read in.
		if self.prepared_proteins is None:
			#But otherwise, we need to read them in.
			self.prepared_proteins, deflines = read_fasta(self.protein)
		
		self.protein_kmer_count = {}
		self.best_hits_kmers = {}

		#Kmerize proteins and record metadata
		for protein in self.prepared_proteins:
			if protein in self.best_hits:
				accession = self.best_hits[protein]
				kmer_set = self.unique_kmer_simple_key(self.prepared_proteins[protein])
				self.protein_kmer_count[accession] = kmer_set.shape[0]
				self.protein_count += 1
				self.best_hits_kmers[accession] = kmer_set
				hit_count += 1
				
			#Free the space either way
			self.prepared_proteins[protein] = None
			
		#Final free.
		self.prepared_proteins = None
			
		#No HMM hits.
		if hit_count == 0:
			self.is_empty = True
			
		self.besthits_time = datetime.datetime.now()
		self.status = "best hits found"

	def preprocess(self):
		self.init_time = datetime.datetime.now()
		#default to 0 time.
		self.prot_pred_time = self.init_time
		self.hmm_search_time = self.init_time
		self.besthits_time = self.init_time
	
		#There's no advancement stage for protein and HMM
		if self.status == 'genome':
			start_time = self.curtime()
			#report = True
			if self.start_time is None:
				self.start_time = start_time
			
			if self.initial_state == "protein+HMM":
				self.initial_state = "genome"
			
			self.genome_to_protein()
			
		if self.status == 'protein':
			start_time = self.curtime()
			#report = True
			if self.start_time is None:
				self.start_time = start_time
				
			if self.initial_state == "protein+HMM":
				self.initial_state = "protein"
			
			self.protein_to_hmm()
			
		if self.status == 'protein and hmm':
			start_time = self.curtime()
			
			if self.start_time is None:
				self.start_time = start_time
				
			self.prot_and_hmm_to_besthits()
		
		#Add an end time if either genome -> protein -> HMM or protein -> HMM happened.
		if self.start_time is not None:
			end_time = self.curtime()
			self.end_time = end_time
		else:
			#Start was protein+HMM. There was no runtime, and intitial state is p+hmm
			#self.initial_state = "protein+HMM"
			self.start_time = "N/A"
			self.end_time = "N/A"
			
		#Protein not generated on this run.
		if self.trans_table is None:
			self.trans_table = "unknown"
		
		self.partial_timings()

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

def find_hmm():
	hmm_path = None
	try:
		#Try to locate the data bundled as it would be with a pip/conda install.
		script_path = os.path.dirname(sys.modules['fastAAI_HMM_models'].__file__)
		if len(script_path) == 0:
			script_path = "."
		hmm_complete_model = os.path.abspath(os.path.normpath(script_path + '/00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm'))
		hmm_path = str(hmm_complete_model)
		#Check that the file exists or fail to the except.
		fh = open(hmm_path)
		fh.close()
	except:
		#Look in the same dir as the script; old method/MiGA friendly
		script_path = os.path.dirname(__file__)
		if len(script_path) == 0:
			script_path = "."
		hmm_complete_model = os.path.abspath(os.path.normpath(script_path +"/"+ "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"))
		hmm_path = str(hmm_complete_model)
	
	return hmm_path

#Build DB from genomes
	
def unique_kmers(seq, ksize):
	n_kmers = len(seq) - ksize + 1
	kmers = []
	for i in range(n_kmers):
		kmers.append(kmer_index[seq[i:i + ksize]])
	#We care about the type because we're working with bytes later.
	return np.unique(kmers).astype(np.int32)

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
	
def split_seq_indices(seq, num_grps):
	newseq = []
	splitsize = 1.0/num_grps*len(seq)
	for i in range(num_grps):
		newseq.append((int(round(i*splitsize)), int(round((i+1)*splitsize)),))
	return newseq
	
	
def list_to_index_dict(list):
	result = {}
	counter = 0
	for item in list:
		result[item] = counter
		counter += 1
	return result
	
	
def rev_list_to_index_dict(list):
	result = {}
	counter = 0
	for item in list:
		result[counter] = item
		counter += 1
	return result
	
def generate_accessions_index(forward = True):
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
	if forward:
		list_of_poss_accs = list_to_index_dict(acc_list)
	else:
		list_of_poss_accs = rev_list_to_index_dict(acc_list)
		
	return list_of_poss_accs

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
						
	parser.add_argument('--threads',  dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')
	parser.add_argument('--verbose',        dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
	parser.add_argument('--compress', dest = "do_comp", action = 'store_true', help = 'Gzip compress generated proteins, HMMs. Off by default.')

		
	args, unknown = parser.parse_known_args()
	
	return parser, args
	
def run_build(input_file):
	input_file.preprocess()
	if len(input_file.best_hits_kmers) < 1:
		input_file.best_hits_kmers = None
		input_file.err_log += " This file did not successfully complete. No SCPs could be found."
		
	return input_file 
	
def acc_transformer_init(db, outpath):
	sqlite3.register_converter("array", convert_array)
	global indb
	indb = db
	global outbase
	outbase = outpath
	global ok
	ok = generate_accessions_index()
	
def acc_transformer(acc_name):
	source = sqlite3.connect(indb)
	scurs = source.cursor()
	
	data = scurs.execute("SELECT * FROM {acc}_genomes".format(acc=acc_name)).fetchall()
	
	scurs.close()
	source.close()
	
	reformat = {}
	
	for row in data:
		genome, kmers = row[0], np.frombuffer(row[1], dtype=np.int32)
		for k in kmers:
			if k not in reformat:
				reformat[k] = []
			reformat[k].append(genome)
	
	data = None
	
	to_add = []
	for k in reformat:
		as_bytes = np.array(reformat[k], dtype = np.int32)
		as_bytes = as_bytes.tobytes()
		reformat[k] = None
		to_add.append((int(k), as_bytes,))
	
	my_acc_db = os.path.normpath(outbase + "/temp/"+acc_name+".db")
	
	if os.path.exists(my_acc_db):
		os.remove(my_acc_db)
	
	my_db  = sqlite3.connect(my_acc_db)
	curs = my_db.cursor()
	curs.execute("CREATE TABLE {acc} (kmer INTEGER PRIMARY KEY, genomes array)".format(acc=acc_name))
	my_db.commit()
	
	curs.executemany("INSERT INTO {acc} VALUES (?, ?)".format(acc = acc_name), to_add)
	
	my_db.commit()
	
	to_add = None
	
	curs.execute("CREATE INDEX {acc}_index ON {acc} (kmer)".format(acc=acc_name))
	my_db.commit()
	
	curs.close()
	my_db.close()
	
	return [my_acc_db, acc_name]
	
def build_db(genomes, proteins, hmms, db_name, output, threads, verbose, do_compress):
	success = True
	
	imported_files = fastaai_file_importer(genomes = genomes, proteins = proteins, hmms = hmms, output = output, compress = do_compress)
	imported_files.determine_inputs()
	
	if imported_files.error:
		print("Exiting FastAAI due to input file error.")
		quit()

	good_to_go = prepare_directories(output, imported_files.status, "query")
	
	db_path = os.path.normpath(output + "/database")
	if not os.path.exists(db_path):
		os.mkdir(db_path)
	
	if not good_to_go:
		print("Exiting FastAAI")
		sys.exit()
	
	print("")

	hmm_path = find_hmm()
	
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
			parent = sqlite3.connect(final_database)
			curs = parent.cursor()
			
			existing_genome_IDs = {}
			sql_command = "SELECT genome, gen_id FROM genome_index"
			for result in curs.execute(sql_command).fetchall():
				genome = result[0]
				id = int(result[1])
				existing_genome_IDs[genome] = id
			
			curs.close()
			parent.close()
	except:
		print("You specified an existing file to be a database, but it does not appear to be a FastAAI database.")
		print("FastAAI will not be able to continue. Please give FastAAI a different database name and continue.")
		print("Exiting.")
		success = False
	
	if success:
		hmm_file = find_hmm()
		if existing_genome_IDs is not None:
			genome_idx = max(list(existing_genome_IDs.values()))+1
		else:
			existing_genome_IDs = {}
			genome_idx = 0
		
		td = os.path.normpath(output+"/temp")
		if not os.path.exists(td):
			os.mkdir(td)
		
		temp_db = os.path.normpath(td+"/FastAAI_temp_db.db")

		if os.path.exists(temp_db):
			os.remove(temp_db)
		
		sqlite3.register_converter("array", convert_array)
		worker = sqlite3.connect(temp_db)
		wcurs = worker.cursor()
		wcurs.execute("CREATE TABLE genome_index (genome text, gen_id integer, protein_count integer)")
		wcurs.execute("CREATE TABLE genome_acc_kmer_counts (genome integer, accession integer, count integer)")
		ok = generate_accessions_index()
		for t in ok:
			wcurs.execute("CREATE TABLE " + t + "_genomes (genome INTEGER PRIMARY KEY, kmers array)")
			
		worker.commit()

		new_gens = []
		new_gak = []
		accs_seen = {}
		if verbose:
			tracker = progress_tracker(total = len(imported_files.in_files), message = "Processing inputs")
		else:
			print("Processing inputs")
			
		#Only build_db makes a log.
		if not os.path.exists(os.path.normpath(output + "/" + "logs")):
			os.mkdir(os.path.normpath(output + "/" + "logs"))
		
		logger = open(os.path.normpath(output+"/logs/"+"FastAAI_preprocessing_log.txt"), "a")
		print("file", "start_date", "end_date", "starting_format", 
		"prot_prediction_time", "trans_table", "hmm_search_time", "besthits_time", 
		"errors", sep = "\t", file = logger)
		
		pool = multiprocessing.Pool(threads, initializer = hmm_preproc_initializer, initargs = (hmm_file, do_compress,))
		
		for result in pool.imap(run_build, imported_files.in_files):
			#log data, regardless of kind
			print(result.basename, result.start_time, result.end_time, result.initial_state, 
			result.prot_pred_time, result.trans_table, result.hmm_search_time, result.besthits_time,
			result.err_log, sep = "\t", file = logger)
			
			if result.best_hits_kmers is not None:
				genome_name = result.original_name
				
				if genome_name in existing_genome_IDs:
					print(genome_name, "Already present in final database and will be skipped.")
					print("")
				else:
					protein_count = result.protein_count
					for acc_name in result.best_hits_kmers:
						if acc_name not in accs_seen:
							accs_seen[acc_name] = 0
						acc_id = ok[acc_name]
						kmer_ct = result.protein_kmer_count[acc_name]
						kmers = result.best_hits_kmers[acc_name]
						kmers = kmers.tobytes()
						wcurs.execute("INSERT INTO {acc}_genomes VALUES (?, ?)".format(acc=acc_name), (genome_idx, kmers,))
						new_gak.append((genome_idx, acc_id, kmer_ct,))
						
					new_gens.append((genome_name, genome_idx, protein_count,))
					genome_idx += 1
					
					worker.commit()
					
			if verbose:	
				tracker.update()
	
		pool.close()
		
		logger.close()
		
		wcurs.executemany("INSERT INTO genome_index VALUES (?,?,?)", new_gens)
		wcurs.executemany("INSERT INTO genome_acc_kmer_counts VALUES (?,?,?)", new_gak)
		worker.commit()
		
		wcurs.close()
		worker.close()
		
		accs_seen = list(accs_seen.keys())
		
		parent = sqlite3.connect(final_database)
		curs = parent.cursor()
				
		curs.execute("attach '" + temp_db + "' as worker")
		#initialize if needed.
		curs.execute("CREATE TABLE IF NOT EXISTS genome_index (genome text, gen_id integer, protein_count integer)")
		curs.execute("CREATE TABLE IF NOT EXISTS genome_acc_kmer_counts (genome integer, accession integer, count integer)")
		
		curs.execute("INSERT INTO genome_index SELECT * FROM worker.genome_index")
		curs.execute("INSERT INTO genome_acc_kmer_counts SELECT * FROM worker.genome_acc_kmer_counts")
		curs.execute("CREATE INDEX IF NOT EXISTS kmer_acc ON genome_acc_kmer_counts (genome, accession);")
		parent.commit()
		
		if verbose:
			tracker = progress_tracker(total = len(accs_seen), message = "Collecting results")
		else:
			print("Collecting results")
		
		pool = multiprocessing.Pool(threads, initializer = acc_transformer_init, initargs = (temp_db, output,))
		
		for result in pool.imap_unordered(acc_transformer, accs_seen):
			database, accession = result[0], result[1]
			curs.execute("CREATE TABLE IF NOT EXISTS {acc} (kmer INTEGER PRIMARY KEY, genomes array)".format(acc=accession))
			curs.execute("CREATE TABLE IF NOT EXISTS {acc}_genomes (genome INTEGER PRIMARY KEY, kmers array)".format(acc=accession))
			curs.execute("CREATE INDEX IF NOT EXISTS {acc}_index ON {acc}(kmer)".format(acc=accession))
			
			#Get the genomes from worker db.
			curs.execute("INSERT INTO {acc}_genomes SELECT * FROM worker.{acc}_genomes".format(acc=accession))
			
			parent.commit()
			
			accdb = sqlite3.connect(database)
			acc_curs = accdb.cursor()
			
			to_update = acc_curs.execute("SELECT kmer, genomes, genomes FROM {acc}".format(acc=accession)).fetchall()
			
			acc_curs.close()
			accdb.close()
			
			update_concat_sql = "INSERT INTO {acc} VALUES (?,?) ON CONFLICT(kmer) DO UPDATE SET genomes=genomes || (?)".format(acc=accession)
			#ON CONFLICT(kmer) DO UPDATE SET genomes=genomes || acc.{acc}.genomes;".format(acc=accession)
			#print(update_concat_sql)
			curs.executemany(update_concat_sql, to_update)
			
			parent.commit()
			
			os.remove(database)
			if verbose:
				tracker.update()
		
		pool.close()
		
		curs.execute("detach worker")
		
		parent.commit()
		
		curs.close()
		parent.close()	
		
		os.remove(temp_db)
		
	shutil.rmtree(td)
		

	if success:
		print("Database build complete!")
		
	return success

def file_v_db_initializer(tgak, tgt_names, tgt_cts, hmm_file, do_compress, tgt_ct, sd, out, style, in_mem, build_q, tdb):
	#num_tgts, self.do_sd, self.output, self.style, self.as_mem_db, self.do_db_build
	global _tdb
	_tdb = tdb
	
	global _tgt_gak
	_tgt_gak = tgak
	
	global _tname
	_tname = tgt_names
	
	global _tct
	_tct = tgt_cts
	
	global hmm_manager
	hmm_manager = pyhmmer_manager(do_compress)
	hmm_manager.load_hmm_from_file(hmm_file)
	
	global num_tgts
	num_tgts = tgt_ct
	
	global _do_sd
	_do_sd = sd
	
	global out_style
	out_style = style
	
	global out_base
	out_base = out
	
	global db_is_in_mem
	db_is_in_mem = in_mem
	
	global make_query_db
	make_query_db = build_q
	
	return _tdb, _tgt_gak, _tname, _tct, hmm_manager, num_tgts, _do_sd, out_base, out_style, db_is_in_mem, make_query_db
	
def file_v_db_worker(query_args):
	#query info for this particular query
	in_file = query_args[0]
	
	in_file.preprocess()
	
	qname = in_file.basename
	
	do_sd = _do_sd
	
	#std dev. calcs are not meaningful with matrix style output.
	if out_style == "matrix":
		do_sd = False
	
	if do_sd:
		results = []
		shared_acc_counts = []
	else:
		results = np.zeros(shape = num_tgts, dtype = np.float_)
		shared_acc_counts = np.zeros(shape = num_tgts, dtype = np.int32)
	
	if db_is_in_mem:
		#The connection is already given as MDB if the db is in mem
		tconn = _tdb
	else:
		#db is on disk and the connection has to be established.
		tconn = sqlite3.connect(_tdb)
		
	tcurs = tconn.cursor()
		
	#This is a difference from the DB-first method.
	acc_idx = generate_accessions_index(forward = True)
	
	genome_lists = {}
	
	tcurs.row_factory = lambda cursor, row: row[0]
	
	
	if make_query_db:
		ret = [qname, None, []]
	else:
		ret = [qname, None, None]
	
	#We need to purge accsessions not in tgt.
	for acc in in_file.best_hits_kmers:
		one = in_file.best_hits_kmers[acc]
		acc_id = acc_idx[acc]
		
		if make_query_db:
			ret[2].append((qname, acc_id, one.tobytes(),))
		
		#Check working.
		if acc_id in _tgt_gak:

			kmer_ct = one.shape[0]

			if do_sd:
				hits = np.zeros(shape = num_tgts, dtype = np.int32)
				hits[np.nonzero(_tgt_gak[acc_id])] = 1
				shared_acc_counts.append(hits)
			else:
				shared_acc_counts[np.nonzero(_tgt_gak[acc_id])] += 1
			
			#SQL has a max binding size of 999, for some reason.
			if kmer_ct > 998:
				#Each kmer needs to be a tuple.
				these_kmers = [(int(kmer),) for kmer in one]
			
				temp_name = "_" + qname +"_" + acc
				temp_name = temp_name.replace(".", "_")
				
				tcurs.execute("CREATE TEMP TABLE " + temp_name + " (kmer INTEGER)")
				tconn.commit()
				insert_table = "INSERT INTO " + temp_name + " VALUES (?)"
				tcurs.executemany(insert_table, these_kmers)
				tconn.commit()
				join_and_select_sql = "SELECT genomes FROM " + temp_name + " INNER JOIN " + acc + " ON "+ temp_name+".kmer = " + acc+".kmer;"

				set = tcurs.execute(join_and_select_sql).fetchall()
			else:
				#kmers must be a list, not a tuple.
				these_kmers = [int(kmer) for kmer in one]
				select = "SELECT genomes FROM " + acc + " WHERE kmer IN ({kmers})".format(kmers=','.join(['?']*len(these_kmers)))
				
				set = tcurs.execute(select, these_kmers).fetchall()
				
			#join results into one bytestring.
			set = b''.join(set)
			
			these_intersections = np.bincount(np.frombuffer(set, dtype = np.int32), minlength = num_tgts)
			set = None
			#Add tgt kmer counts to query kmer counts, find union size based on intersection size, cald jacc
			jacc = np.divide(these_intersections, np.subtract(np.add(_tgt_gak[acc_id], kmer_ct), these_intersections))
			
			if do_sd:
				results.append(jacc)
			else:
				results += jacc
			
	tcurs.row_factory = None
	tcurs.close()
	
	if do_sd:
		results = np.vstack(results)
		has_accs = np.vstack(shared_acc_counts)
		
		shared_acc_counts = np.sum(has_accs, axis = 0)
		
		#final jacc_means
		jaccard_averages = np.divide(np.sum(results, axis = 0), shared_acc_counts)
		
		aai_ests = numpy_kaai_to_aai(jaccard_averages)
		
		#find diffs from means; this includes indicies corresponding to unshared SCPs that should not be included. 
		results = results - jaccard_averages
		
		#fix those corresponding indicies to not contribute to the final SD.
		results[np.nonzero(has_accs == 0)] = 0
		
		#Square them
		results = np.square(results)
		#Sum squares and divide by shared acc. count, the sqrt to get SD.
		jaccard_SDs = np.sqrt(np.divide(np.sum(results, axis = 0), shared_acc_counts))
		jaccard_SDs = np.round(jaccard_SDs, 4).astype(str)
				
	else:
		#other condition.
		jaccard_SDs = None
		jaccard_averages = np.divide(results, shared_acc_counts)
		#we don't want to pass char arrays to main, so skip this here and do it in main instead.
		if out_style != "matrix":
			aai_ests = numpy_kaai_to_aai(jaccard_averages)
		
	del results

	#Since the outputs go to separate files, it makes more sense to do them within the worker processes instead of in main.
	if out_style == "tsv":
		no_hit = np.where(shared_acc_counts == 0)
		
		possible_hits = np.minimum(len(in_file.best_hits_kmers), _tct).astype(str)
		jaccard_averages = np.round(jaccard_averages, 4).astype(str)
		shared_acc_counts = shared_acc_counts.astype(str)
		
		jaccard_averages[no_hit] = "N/A"
		aai_ests[no_hit] = "N/A"
		shared_acc_counts[no_hit] = "N/A"
		possible_hits[no_hit] = "N/A"
			
		output_name = os.path.normpath(out_base + "/"+qname+"_results.txt")
	
		out = open(output_name, "w")
		out.write("query\ttarget\tavg_jacc_sim\tjacc_SD\tnum_shared_SCPs\tposs_shared_SCPs\tAAI_estimate\n")
		if do_sd:
			jaccard_SDs[no_hit] = "N/A"
			for i in range(0, len(aai_ests)):
				out.write(qname+"\t"+_tname[i]+"\t"+jaccard_averages[i]+"\t"+jaccard_SDs[i]+"\t"+shared_acc_counts[i]+"\t"+possible_hits[i]+"\t"+aai_ests[i]+"\n")
		else:			
			for i in range(0, len(aai_ests)):
				out.write(qname+"\t"+_tname[i]+"\t"+jaccard_averages[i]+"\t"+"N/A"+"\t"+shared_acc_counts[i]+"\t"+possible_hits[i]+"\t"+aai_ests[i]+"\n")
		out.close()

		
	#We're just gonna pass this back to the main to print.
	if out_style == "matrix":
		ret[1] = jaccard_averages
	
	return ret
	
#Handles both query and target types for a db vs db query
class file_vs_db_query:
	def __init__(self, in_memory = False, input_file_objects = None,
	target = None, threads = 1, do_sd = False, output_base = "FastAAI", output_style = "tsv", 
	build_db_from_queries = True, qdb_name = "Query_FastAAI_database.db", hmm_path = None, 
	do_comp = True, verbose = True):
		#files to work with
		self.queries = input_file_objects
		self.do_db_build = build_db_from_queries
		self.dbname = qdb_name
		
		self.t = target
		self.valids = None
		
		#Originally this was made to be a memory database only block of code, but just if/else one change makes it work on disk and it doesn't need a redev, then.
		self.as_mem_db = in_memory
		
		self.t_conn = None
		self.t_curs = None
		
		self.threads = threads
		self.do_sd = do_sd

		self.output_base = output_base
		self.output = os.path.normpath(output_base + "/results")
		self.style = output_style
		
		if hmm_path is not None:
			self.hmm_path = hmm_path
		else:
			self.hmm_path = find_hmm()
			
		self.do_comp = do_comp
		
		self.verbose = verbose

	'''
	Workflow is:
		load target db as mem (optional)
		assess valid targets
		create query db output (optional)
		pass query args to workers
		preproc query args
		write results
		fill query_db_out (optional)
	'''
		
		
	def open(self):
		if self.as_mem_db:
			self.t_conn = sqlite3.connect(':memory:')
		else:
			self.t_conn = sqlite3.connect(self.t)
			
		self.t_curs = self.t_conn.cursor()	
		
		if self.as_mem_db:
			self.t_curs.execute("attach '" + self.t + "' as targets")
			
			self.t_curs.execute("CREATE TABLE genome_index AS SELECT * FROM targets.genome_index")
			self.t_curs.execute("CREATE TABLE genome_acc_kmer_counts AS SELECT * FROM targets.genome_acc_kmer_counts")
			self.t_curs.execute("CREATE INDEX t_gi ON genome_index (gen_id)")
			self.t_curs.execute("CREATE INDEX t_gak ON genome_acc_kmer_counts (accession)")
		
		if self.as_mem_db:
			table_sql = "SELECT name FROM targets.sqlite_master"
		else:
			table_sql = "SELECT name FROM sqlite_master"
		
		
		ok = generate_accessions_index()
		ok_names = set(list(ok.keys()))
		successful_tables = []
		
		for name in self.t_curs.execute(table_sql).fetchall():
			name = name[0]
			if name in ok_names:
				successful_tables.append(ok[name])
				if self.as_mem_db:
					self.t_curs.execute("CREATE TABLE " + name + " AS SELECT * FROM targets."+name)
					self.t_curs.execute("CREATE INDEX "+name+"_index ON " + name+" (kmer)" )
					
		if self.as_mem_db:
			self.t_conn.commit()
			self.t_curs.execute("detach targets")
			
		self.valids = tuple(successful_tables)
			
	def close(self):
		self.t_curs.close()
		self.t_curs = None
		
	def clean_up(self):
		self.t_conn.close()
		self.t_conn = None
		
	def sqlite_table_schema(self, conn, name):
		"""Return a string representing the table's CREATE"""
		cursor = conn.execute("SELECT sql FROM sqlite_master WHERE name=?;", [name])
		sql = cursor.fetchone()[0]
		cursor.close()
		return sql
		
	def execute(self):
		print("FastAAI is running.")
		tgt_id_res = self.t_curs.execute("SELECT * FROM genome_index ORDER BY gen_id").fetchall()
		
		tgt_ids = []
		tgt_naming = []
		tgt_counts = []
		for r in tgt_id_res:
			genome, id, prot_ct = r[0], r[1], r[2]
			tgt_ids.append(genome)
			tgt_naming.append(genome)
			tgt_counts.append(prot_ct)
			
		num_tgts = len(tgt_ids)
		tgt_counts = np.array(tgt_counts, dtype = np.int32)
			
		tgts_gak = {}
		gak_sql = "SELECT * FROM genome_acc_kmer_counts WHERE accession in ({accs})".format(accs=','.join(['?']*len(self.valids)))
		
		for result in self.t_curs.execute(gak_sql, self.valids).fetchall():
			genome, acc, ct = result[0], result[1], result[2]
			if acc not in tgts_gak:
				tgts_gak[acc] = np.zeros(num_tgts, dtype = np.int32)
			tgts_gak[acc][genome] += ct
			
		#If the DB is a memory DB, we need to maintain the connection, but neither needs to maintain the curor in main.
		self.close()
			
		query_groups = []
		
		for query_input in self.queries:
			query_groups.append((query_input,))
		
		#And if it's a physical database, we do want to close it.
		if not self.as_mem_db:
			self.t_conn.close()
		
		num_queries = len(query_groups)
		
		if self.do_db_build:
			sqlite3.register_converter("array", convert_array)
			qdb_path = os.path.normpath(self.output_base + "/database/"+self.dbname)
			if not os.path.exists(os.path.normpath(self.output_base + "/database")):
				try:
					os.mkdir(os.path.normpath(self.output_base + "/database"))
				except:
					print("Couldn't make database at", qdb_path)
					self.do_db_build = False
					
			if os.path.exists(qdb_path):
				print("Database for queries already exists. I can't make one at:", qdb_path)
				self.do_db_build = False
			else:
				query_db_conn = sqlite3.connect(qdb_path)
				q_curs = query_db_conn.cursor()
				q_curs.execute("CREATE TABLE storage (genome INTEGER, accession INTEGER, kmers array)")
				q_curs.execute("CREATE INDEX store_idx ON storage (genome, accession)")
				query_genome_index = []
				qgi_ct = 0
				qg_gak = []
		
		if self.verbose:
			tracker = progress_tracker(total = num_queries, message = "Calculating AAI...", one_line = True)
		
		if self.style == "matrix":
			output_name = os.path.normpath(self.output + "/FastAAI_matrix.txt")
			output = open(output_name, "w")
			#needs target names.
			print("query_genome", *tgt_ids, sep = "\t", file = output)
		
		#Need to pass these
		
		#both initializers will share this.
		shared_args = [tgts_gak, tgt_naming, tgt_counts, self.hmm_path, self.do_comp, num_tgts, self.do_sd, self.output, 
			self.style, self.as_mem_db, self.do_db_build]
		
		if self.as_mem_db:
			shared_args.append(self.t_conn)
			shared_args = tuple(shared_args)
			pool = multiprocessing.Pool(self.threads, initializer = file_v_db_initializer, initargs = shared_args)
		else:
			#db is on disk,
			shared_args.append(self.t)
			shared_args = tuple(shared_args)
			pool = multiprocessing.Pool(self.threads, initializer = file_v_db_initializer, initargs = shared_args)

		for result in pool.imap(file_v_db_worker, query_groups):
			if self.verbose:
				tracker.update()
			qname = result[0]
			if self.style == "matrix":
				printout = numpy_kaai_to_aai(result[1])
				print(qname, *printout, sep = "\t", file = output)
				
			if self.do_db_build:
				query_genome_index.append((qname, qgi_ct, len(result[2]),))
				for row in result[2]:
					num_kmers = int(len(row[2])/4)
					qg_gak.append((qgi_ct, row[1], num_kmers,))
				qgi_ct += 1
				q_curs.executemany("INSERT INTO storage VALUES (?, ?, ?)", result[2])
				query_db_conn.commit()
				
		pool.close()
		
		if self.style == "matrix":
			output.close()
			
		if self.do_db_build:
			q_curs.execute("CREATE TABLE genome_index (genome text, gen_id integer, protein_count integer)")
			q_curs.execute("CREATE TABLE genome_acc_kmer_counts (genome integer, accession integer, count integer)")
			q_curs.executemany("INSERT INTO genome_index VALUES (?,?,?)", query_genome_index)
			q_curs.executemany("INSERT INTO genome_acc_kmer_counts VALUES (?,?,?)", qg_gak)
			query_db_conn.commit()
			
			acc_id_to_name = generate_accessions_index(forward = False)
			qgi_dict = {}
			for tup in query_genome_index:
				qgi_dict[tup[0]] = tup[1]
			
			accs_in_db = q_curs.execute("SELECT DISTINCT(accession) FROM genome_acc_kmer_counts").fetchall()
			if self.verbose:
				tracker = progress_tracker(total = len(accs_in_db), message = "Crafting database from query outputs.", one_line = True)
			
			for acc in accs_in_db:
				acc = acc[0]
				acc_name = acc_id_to_name[acc]
				q_curs.execute("CREATE TABLE " + acc_name + " (kmer INTEGER PRIMARY KEY, genomes array)")
				q_curs.execute("CREATE TABLE " + acc_name + "_genomes (genome INTEGER PRIMARY KEY, kmers array)")
				data = q_curs.execute("SELECT genome, kmers FROM storage WHERE accession = ?", (acc,)).fetchall()
				
				ins = []
				#group by kmer
				kmers_by_gen = {}
				for row in data:
					gen = row[0]
					gen = qgi_dict[gen]
					kmers = np.frombuffer(row[1], dtype = np.int32)
					ins.append((gen, kmers,))
					for k in kmers:
						#typecast
						k = int(k)
						if k not in kmers_by_gen:
							kmers_by_gen[k] = []
						kmers_by_gen[k].append(gen)
						
				data = None
				
				q_curs.executemany("INSERT INTO "+ acc_name + "_genomes VALUES (?,?)", ins)
				
				ins = []
				for k in kmers_by_gen:
					dat = kmers_by_gen[k]
					dat = np.sort(np.array(dat, dtype = np.int32))
					ins.append((k, dat.tobytes()))
				
				q_curs.executemany("INSERT INTO "+ acc_name + " VALUES (?,?)", ins)
				
				ins = None
				
				query_db_conn.commit()
				
				q_curs.execute("CREATE INDEX IF NOT EXISTS " + acc_name + "_index ON " + acc_name + " (kmer)")
				
				if self.verbose:
					tracker.update()
				
				
			q_curs.execute("CREATE INDEX IF NOT EXISTS kmer_acc ON genome_acc_kmer_counts (genome, accession);")
			q_curs.execute("DROP INDEX store_idx")
			q_curs.execute("DROP TABLE storage")
			query_db_conn.commit()
			q_curs.execute("VACUUM")
			query_db_conn.commit()
			q_curs.close()
			query_db_conn.close()

	#Actually run the thing.
	def run(self):
		self.open()
		self.execute()
		#Clean up the db connections; free the mem.
		self.clean_up()
	
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
	
#Also includes a multiply by 100 and type conversion compared to original - this is some silliness for saving memory.
def numpy_kaai_to_aai_just_nums(kaai_array, as_float = False):
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
	
	aai_hat_array = np.round(aai_hat_array, 2)
	
	#<30 and >90 values
	smol = np.where(aai_hat_array < 30)
	big  = np.where(aai_hat_array > 90)
	
	#We can find these later.
	aai_hat_array[smol] = 15
	aai_hat_array[big]  = 95
	
	if as_float:
		aai_hat_array = np.round(aai_hat_array, 2)
	else:
		aai_hat_array = np.multiply(aai_hat_array, 100)
		aai_hat_array = np.round(aai_hat_array, 2)
		aai_hat_array = aai_hat_array.astype(np.int16)
	
	return aai_hat_array

	
def curtime():
	time_format = "%d/%m/%Y %H:%M:%S"
	timer = datetime.datetime.now()
	time = timer.strftime(time_format)
	return time

#Perform a minimal-memory query of a target database from input files. Lighter weight function for low memory
def sql_query_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	This FastAAI module takes one or many genomes, proteins, or proteins and HMMs as a QUERY and searches them against an existing FastAAI database TARGET using SQL
	If you only have a few genomes - or not enough RAM to hold the entire target database in memory - this is the probably the best option for you.
	
	To provide files, supply either a directory containing only one type of file (e.g. only genomes in FASTA format), a file containing paths to files of a type, 1 per line,
	or a comma-separated list of files of a single type (no spaces)
	
	If you provide FastAAI with genomes or only proteins (not proteins and HMMs), this FastAAI module will produce the required protein and HMM files as needed
	and place them in the output directory, just like it does while building a database. 
	
	Once these inputs are ready to be queried against the database (each has both a protein and HMM file), they will be processed independently, 1 per thread at a time.
	
	Note: Protein and HMM files generated during this query can be supplied to build a FastAAI database from proteins and HMMs using the build_db module, without redoing preprocessing.
	''')
			
	parser.add_argument('-g', '--genomes',  dest = 'genomes', default = None,  help = 'Genomes in FASTA format.')
	parser.add_argument('-p', '--proteins', dest = 'proteins', default = None, help = 'Protein amino acids in FASTA format.')
	parser.add_argument('-m', '--hmms',     dest = 'hmms', default = None,     help = 'HMM search files produced by FastAAI on a set of proteins.')
	
	parser.add_argument('--target',   dest = 'target', default = None, help =   'A path to the FastAAI database you wish to use as the target')
	
	parser.add_argument('-o', '--output',   dest = 'output', default = "FastAAI", help = 'The directory where FastAAI will place the result of this query and any protein or HMM files it has to generate. By default, a directory named "FastAAI" will be created in the current working directory and results will be placed there.')
	parser.add_argument('--output_style', dest = "style", default = 'tsv', help = "Either 'tsv' or 'matrix'. Matrix produces a simplified output of only AAI estimates.")
	parser.add_argument('--do_stdev',   dest = "do_stdev", action='store_true',   help = 'Off by default. Calculate std. deviations on Jaccard indicies. Increases memory usage and runtime slightly. Does NOT change estimated AAI values at all.')
	
	parser.add_argument('--threads',  dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')
	parser.add_argument('--verbose',        dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
	
	parser.add_argument('--in_memory', dest = "in_mem", action = 'store_true', help = 'Load the target database into memory before querying. Consumes more RAM, but is faster and reduces file I/O substantially.')
	
	parser.add_argument('--create_query_db', dest = "make_db", action = 'store_true', help = 'Create a query database from the genomes.')
	parser.add_argument('--query_db_name', dest = "qdb_name", default = "Query_FastAAI_db.db", help = 'Name the query database. This file must not already exist.')
	
	parser.add_argument('--compress', dest = "do_comp", action = 'store_true', help = 'Gzip compress generated proteins, HMMs. Off by default.')

	args, unknown = parser.parse_known_args()
	
	return parser, args

def sql_query_thread_starter(kmer_cts, protein_cts):
	global target_kmer_cts 
	global target_protein_counts
	target_kmer_cts = kmer_cts
	target_protein_counts = protein_cts

#took a function from fastaai 2.0
class fastaai_file_importer:
	def __init__(self, genomes = None, proteins = None, hmms = None, crystals = None, output = "FastAAI", compress = False):
		#genomes, prots, hmms can be supplied as either directory, a file with paths 1/line, or comma-sep paths. Type is determined automatically.
		self.genomes = genomes
		self.proteins = proteins
		self.hmms = hmms
		self.crystals = crystals
		
		self.genome_list = None
		self.protein_list = None
		self.hmm_list = None
		self.crystal_list = None
		
		#file base names.
		self.identifiers = None
		
		self.error = False
		
		self.in_files = None
		
		self.status = "genome"
		self.output = output
		
		self.do_comp = compress
	
	def retrieve_files(self, arg):	
		done = False
		files = []
		names = []
		#Case where a directory is supplied.
		if os.path.isdir(arg):
			for file in sorted(os.listdir(arg)):
				#Retrieve file name
				if file.endswith(".gz"):
					name = os.path.splitext(os.path.basename(file[:-3]))[0]
				else:
					name = os.path.splitext(os.path.basename(file))[0]
				
				names.append(name)				
				files.append(os.path.abspath(os.path.normpath(arg + '/' +file)))
				
			done = True
				
			
		#Case where a file containing paths is supplied.
		if os.path.isfile(arg):
			handle = agnostic_reader(arg)
			for line in handle:
				file = line.strip()
				if os.path.exists(file):
					if file.endswith(".gz"):
						name = os.path.splitext(os.path.basename(file[:-3]))[0]
					else:
						name = os.path.splitext(os.path.basename(file))[0]
					
					names.append(name)				
					files.append(os.path.abspath(os.path.normpath(file)))
					
			handle.close()
			done = True
			
			if len(names) == 0 and len(files) == 0:
				#Try interpreting the file as a singular path.
				done = False
							
		#Last check.
		if not done:
			for file in arg.split(","):
				if os.path.exists(file):
					if file.endswith(".gz"):
						name = os.path.splitext(os.path.basename(file[:-3]))[0]
					else:
						name = os.path.splitext(os.path.basename(file))[0]
					
					names.append(name)				
					files.append(os.path.abspath(os.path.normpath(file)))
				
		return files, names
	
	#Check if g/p/h
	def determine_inputs(self):
		if self.genomes is not None:
			self.genome_list, self.identifiers = self.retrieve_files(self.genomes)
			if self.proteins is not None or self.hmms is not None:
				print("You can supply genomes or proteins or proteins and HMMS, but not genomes and anything else.")
				self.error = True
			
		#Proteins, but no HMMs
		if self.proteins is not None and self.hmms is None:
			self.protein_list, self.identifiers = self.retrieve_files(self.proteins)
			
		if self.proteins is not None and self.hmms is not None:
			self.protein_list, prot_names = self.retrieve_files(self.proteins)
			self.hmm_list, hmm_names = self.retrieve_files(self.hmms)
			
			if len(self.protein_list) != len(self.hmm_list):
				print("Different number of proteins and HMMs supplied. You must supply the same number of each, and they must be matched pairs.")
				self.error = True
			else:
				all_same = True
				for p, h in zip(prot_names, hmm_names):
					if p != h:
						all_same = False
				
				if all_same:
					self.identifiers = prot_names
					prot_names = None
					hmm_names = None
				else:
					self.error = True
					
			if self.crystals is not None:
				self.crystal_list, self.identifiers = self.retrieve_files(self.crystals)
				
		if not self.error:
			self.prep_input_files()
			
	def prep_input_files(self):
		self.in_files = []
		if self.genome_list is not None:
			self.status = "genome"
			for g in self.genome_list:
				f = input_file(g, output = self.output, do_compress = self.do_comp)
				f.set_genome(g)
				self.in_files.append(f)
			
		if self.protein_list is not None:
			self.status = "protein"
			for p in self.protein_list:
				f = input_file(p, output = self.output, do_compress = self.do_comp)
				f.set_protein(p)
				self.in_files.append(f)
			
		if self.hmm_list is not None:
			self.status = "protein+HMM"
			for h, f in zip(self.hmm_list, self.in_files):
				f.set_hmm(h)
	
def sql_query(genomes, proteins, hmms, db_name, output, threads, verbose, do_stdev, style, in_mem, make_db, qdb_name, do_comp):
	
	if not os.path.exists(db_name):
		print("")
		print("FastAAI can't find your database:", db_name)
		print("Are you sure that the path you've given to the database is correct and that the database exists?")
		print("FastAAI exiting.")
		print("")
		sys.exit()
		
	#importer opts
	#genomes = None, proteins = None, hmms = None, crystals = None
	imported_files = fastaai_file_importer(genomes = genomes, proteins = proteins, hmms = hmms, output = output)
	imported_files.determine_inputs()

	if imported_files.error:
		print("Exiting FastAAI due to input file error.")
		quit()

	good_to_go = prepare_directories(output, imported_files.status, "query")
	
	if not good_to_go:
		print("Exiting FastAAI")
		sys.exit()
	
	print("")

	'''
	self, in_memory = False, input_file_objects = None,
	target = None, threads = 1, do_sd = False, output_base = "FastAAI", output_style = "tsv", 
	build_db_from_queries = True, qdb_name = "Query_FastAAI_database.db", hmm_path = "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm", 
	do_comp = True, verbose = True
	'''
	hmm_path = find_hmm()
	
	mdb = file_vs_db_query(in_memory = in_mem, input_file_objects = imported_files.in_files, target=db_name, 
	threads = threads, output_base = output, do_sd = do_stdev, output_style = style, do_comp = do_comp, 
	build_db_from_queries = make_db, qdb_name = qdb_name, verbose = verbose, hmm_path = hmm_path)
	
	mdb.run()
	
	#Here's where the querying db comes in

		
	print("FastAAI query complete! Results at:", os.path.normpath(output + "/results"))
	return None
	
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
	parser.add_argument('--output_style', dest = "style", default = 'tsv', help = "Either 'tsv' or 'matrix'. Matrix produces a simplified output of only AAI estimates.")
	parser.add_argument('--do_stdev',   dest = "do_stdev", action='store_true',                 help = 'Off by default. Calculate std. deviations on Jaccard indicies. Increases memory usage and runtime slightly. Does NOT change estimated AAI values at all.')

	parser.add_argument('--threads',      dest = 'threads', type=int, default = 1, help = 'The number of processors to use. Default 1.')
	parser.add_argument('--verbose',      dest = 'verbose', action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
	parser.add_argument('--in_memory', dest = "in_mem", action = 'store_true', help = 'Load both databases into memory before querying. Consumes more RAM, but is faster and reduces file I/O substantially. Consider reducing number of threads')
	parser.add_argument('--store_results', dest = "storage", action = 'store_true', help = 'Keep partial results in memory. Only works with --in_memory. Fewer writes, but more RAM. Default off.')
	
	args, unknown = parser.parse_known_args()
	
	return parser, args
	
	
#db-db query; in-mem		
def parse_db_init(query, target, outpath):
	global qdb
	qdb = query
	global tdb
	tdb = target
	global output_path
	output_path = outpath
	
	global query_gak
	global target_gak
	
	return qdb, tdb, output_path
	
def parse_accession(acc):
	tmp = sqlite3.connect(":memory:")
	curs = tmp.cursor()
	
	curs.execute("attach '" + qdb + "' as queries")
	curs.execute("attach '" + tdb + "' as targets")
	
	sql = '''
	SELECT queries.{acc}.genomes, targets.{acc}.genomes 
	FROM queries.{acc} INNER JOIN targets.{acc} 
	ON queries.{acc}.kmer=targets.{acc}.kmer
	'''.format(acc = acc)
	
	res = curs.execute(sql).fetchall()
	
	curs.execute("detach queries")
	curs.execute("detach targets")
	
	curs.close()
	tmp.close()
	
	tl = []
	ql = {}
	
	acc_id = generate_accessions_index()
	acc_id = acc_id[acc]
	
	indexer = 0
	for r in res:
		queries = np.frombuffer(r[0], dtype = np.int32)
		tgt = np.frombuffer(r[1], dtype = np.int32)
		tl.append(tgt)
		
		for q in queries:
			if q not in ql:
				ql[q] = {}
			if acc_id not in ql[q]:
				ql[q][acc_id] = []
				
			ql[q][acc_id].append(indexer)
		
		indexer += 1
		
	tl = np.array(tl, dtype = object)
			
	for q in ql:
		if acc_id in ql[q]:
			ql[q][acc_id] = np.array(ql[q][acc_id], dtype=np.int32)
	
	out_file = os.path.normpath(output_path+"/temp/"+acc+".pickle")
	
	with open(out_file, "wb") as out:
		pickle.dump([ql, tl], out)
		
	return([acc, out_file])
			
#all of this is exclusive to the in-mem approach for db db query
def one_init(ql, tl, num_tgt, qgak_queue, tgak, tpres, sd, sty, temp_path, store_results, progress_queue, qnames, tnames):
	global _ql
	_ql = ql
	global _tl
	_tl = tl
	global _nt
	_nt = num_tgt
	
	qgak_data = qgak_queue.get()
	
	global out_base
	out_base = temp_path
	
	global group_id
	group_id = os.path.normpath(temp_path + "/temp/partial_results_group_" + str(qgak_data[0])+ ".txt")
	
	global _qgak
	_qgak = qgak_data[1]
	
	global query_grouping
	query_grouping = qgak_data[2]
	
	qgak_data = None
	
	global _tgak
	_tgak = tgak
	
	global _tpres
	_tpres = tpres
	
	global _tct
	_tct = np.sum(_tpres, axis = 0)
	
	global do_sd
	do_sd = sd
	global style
	style = sty
	#Suppress div by zero warning - it's handled.
	np.seterr(divide='ignore')
	
	global store
	store = store_results
	if store:
		global holder
		holder = []
	else:
		global outwriter
		outwriter = open(group_id, "w")
		
	global prog_queue
	prog_queue = progress_queue
	
	global _qnames
	_qnames = qnames
	
	global _tnames
	_tnames = tnames
	
def one_work(placeholder):	
	for q in query_grouping:
		results = []
		#We also need to count the accs in the query genome, but which are not part of the inner join.
		for acc in _qgak[q][0]:
			if acc in _ql[q]:
				#the bincount is intersections.
				these_intersections = np.bincount(np.concatenate(_tl[acc][_ql[q][acc]]), minlength = _nt)
			else:
				#there are no intersections even though this accession is shared with at least one target
				#number of intersects is all zeros
				these_intersections = np.zeros(_nt, dtype = np.int32)
			
			#Append the counts or zeros, either way.
			results.append(these_intersections)
				
		results = np.vstack(results)
		
		target_kmer_counts = _tgak[_qgak[q][0], :]
		
		#unions = size(A) + size(B) - size(intersections(A, B)) 
		#unions = target_kmer_counts + query_kmers_by_acc - intersections
		unions = np.subtract(np.add(target_kmer_counts, _qgak[q][1][:, None]), results)
		
		#These are now jaccards, not #intersections
		results = np.divide(results, unions)
		
		shared_acc_counts = np.sum(_tpres[_qgak[q][0], :], axis = 0)
		
		no_hit = np.where(shared_acc_counts == 0)
		
		jaccard_averages = np.divide(np.sum(results, axis = 0), shared_acc_counts)

		#Skip SD if output is matrix
		if style == "tsv":
			aai_ests = numpy_kaai_to_aai(jaccard_averages)

			if do_sd:
				#find diffs from means; this includes indicies corresponding to unshared SCPs that should not be included. 
				results = results - jaccard_averages
				
				#fix those corresponding indicies to not contribute to the final SD.
				results[np.logical_not(_tpres[_qgak[q][0], :])] = 0
				#results[np.nonzero(has_accs == 0)] = 0
				
				#Square them; 0^2 = 0, so we don't have to think about the fixed indices any more.
				results = np.square(results)
				#Sum squares and divide by shared acc. count, the sqrt to get SD.
				jaccard_SDs = np.sqrt(np.divide(np.sum(results, axis = 0), shared_acc_counts))
				jaccard_SDs = np.round(jaccard_SDs, 4).astype(str)
				
			no_hit = np.where(shared_acc_counts == 0)
			
			#addtl.shape[0] is the query acc count
			possible_hits = np.minimum(_qgak[q][0].shape[0], _tct).astype(str)
			
			jaccard_averages = np.round(jaccard_averages, 4).astype(str)
			shared_acc_counts = shared_acc_counts.astype(str)
			
			jaccard_averages[no_hit] = "N/A"
			aai_ests[no_hit] = "N/A"
			shared_acc_counts[no_hit] = "N/A"
			possible_hits[no_hit] = "N/A"
				
			qname = _qnames[q]
				
			output_name = os.path.normpath(out_base + "/results/"+qname+"_results.txt")
		
			out = open(output_name, "w")
			out.write("query\ttarget\tavg_jacc_sim\tjacc_SD\tnum_shared_SCPs\tposs_shared_SCPs\tAAI_estimate\n")
			if do_sd:
				jaccard_SDs[no_hit] = "N/A"
				for i in range(0, len(aai_ests)):
					out.write(qname+"\t"+_tnames[i]+"\t"+jaccard_averages[i]+"\t"+jaccard_SDs[i]+"\t"+shared_acc_counts[i]+"\t"+possible_hits[i]+"\t"+aai_ests[i]+"\n")
			else:			
				for i in range(0, len(aai_ests)):
					out.write(qname+"\t"+_tnames[i]+"\t"+jaccard_averages[i]+"\t"+"N/A"+"\t"+shared_acc_counts[i]+"\t"+possible_hits[i]+"\t"+aai_ests[i]+"\n")
			out.close()
			
			
		else:
			if store:
				aai_ests = numpy_kaai_to_aai_just_nums(jaccard_averages, as_float = False)
				aai_ests[no_hit] = 0
				#add zeros at misses/NAs
				holder.append(aai_ests)
			else:
				aai_ests = numpy_kaai_to_aai_just_nums(jaccard_averages, as_float = True)
				aai_ests[no_hit] = 0
				print(*aai_ests, sep = "\t", file = outwriter)
		
		prog_queue.put(q)
	
	prog_queue.put("done")
	
	return None

def two_work(i):
	if store:
		hold_together = np.vstack(holder)
		np.savetxt(group_id, hold_together, delimiter = "\t", fmt='%4d')
	else:
		outwriter.close()
		
	return group_id
	
def on_disk_init(query_database_path, target_database_path, num_tgt, query_queue, target_gak, tpres, sd, sty, temp_path, progress_queue, qnames, tnames, valids):
	global database
	database = sqlite3.connect(":memory:")
	
	curs = database.cursor()
	curs.execute("attach '" + query_database_path + "' as queries")
	curs.execute("attach '" + target_database_path + "' as targets")
	curs.close()
	
	global _nt
	_nt = num_tgt
	
	qgak_data = query_queue.get()
	
	global out_base
	out_base = temp_path
	
	global group_id
	group_id = os.path.normpath(temp_path + "/temp/partial_results_group_" + str(qgak_data[0])+ ".txt")
	
	global _qgak
	_qgak = qgak_data[1]
	
	global query_grouping
	query_grouping = qgak_data[2]
	
	global _tgak
	_tgak = target_gak
	
	global _tpres
	_tpres = tpres
	
	global _tct
	_tct = np.sum(_tpres, axis = 0)
	
	global do_sd
	do_sd = sd
	global style
	style = sty
	#Suppress div by zero warning - it's handled.
	np.seterr(divide='ignore')
	
	if style == "matrix":
		global outwriter
		outwriter = open(group_id, "w")
		
	global prog_queue
	prog_queue = progress_queue
	
	global _qnames
	_qnames = qnames
	
	global _tnames
	_tnames = tnames
	
	global acc_indexer
	acc_indexer = generate_accessions_index(forward = False)
	
	global _valids
	_valids = valids
	
def on_disk_work_one(placeholder):
	curs = database.cursor()
	for q in query_grouping:
		results = []
		qname = _qnames[q]
		for acc in _qgak[q][0]:
			acc_name = acc_indexer[acc]
				
			if acc_name in _valids:
				
				one = curs.execute("SELECT kmers FROM queries."+acc_name+"_genomes WHERE genome=?", (str(q),)).fetchone()[0]
				one = np.frombuffer(one, dtype = np.int32)
				
				if one.shape[0] > 998:
					#Each kmer needs to be a tuple.
					these_kmers = [(int(kmer),) for kmer in one]
				
					temp_name = "_" + qname +"_" + acc_name
					temp_name = temp_name.replace(".", "_")
					
					curs.execute("CREATE TEMP TABLE " + temp_name + " (kmer INTEGER)")
					insert_table = "INSERT INTO " + temp_name + " VALUES (?)"
					curs.executemany(insert_table, these_kmers)
					
					join_and_select_sql = "SELECT genomes FROM " + temp_name + " INNER JOIN targets." + acc_name + " ON "+ temp_name+".kmer = targets." + acc_name + ".kmer;"
					
					matches = curs.execute(join_and_select_sql).fetchall()
				else:
					#kmers must be a list, not a tuple.
					these_kmers = [int(kmer) for kmer in one]
					select = "SELECT genomes FROM targets." + acc_name + " WHERE kmer IN ({kmers})".format(kmers=','.join(['?']*len(these_kmers)))
					matches = curs.execute(select, these_kmers).fetchall()
				
				set = []
				for row in matches:
					set.append(row[0])
				set = b''.join(set)
					
				matches = None
				these_intersections = np.bincount(np.frombuffer(set, dtype = np.int32), minlength = _nt)
				set = None
				results.append(these_intersections)
				
			else:
				results.append(np.zeros(_nt, dtype=np.int32))
			
		results = np.vstack(results)
		
		target_kmer_counts = _tgak[_qgak[q][0], :]
		
		#unions = size(A) + size(B) - size(intersections(A, B)) 
		#unions = target_kmer_counts + query_kmers_by_acc - intersections
		unions = np.subtract(np.add(target_kmer_counts, _qgak[q][1][:, None]), results)
		
		#These are now jaccards, not #intersections
		results = np.divide(results, unions)

		shared_acc_counts = np.sum(_tpres[_qgak[q][0], :], axis = 0)
		
		no_hit = np.where(shared_acc_counts == 0)
		
		jaccard_averages = np.divide(np.sum(results, axis = 0), shared_acc_counts)
		
		#Skip SD if output is matrix
		if style == "tsv":
			aai_ests = numpy_kaai_to_aai(jaccard_averages)
			
			if do_sd:
				#find diffs from means; this includes indicies corresponding to unshared SCPs that should not be included. 
				results = results - jaccard_averages
				
				#fix those corresponding indicies to not contribute to the final SD.
				results[np.logical_not(_tpres[_qgak[q][0], :])] = 0
				#results[np.nonzero(has_accs == 0)] = 0
				
				#Square them; 0^2 = 0, so we don't have to think about the fixed indices any more.
				results = np.square(results)
				#Sum squares and divide by shared acc. count, the sqrt to get SD.
				jaccard_SDs = np.sqrt(np.divide(np.sum(results, axis = 0), shared_acc_counts))
				jaccard_SDs = np.round(jaccard_SDs, 4).astype(str)
				
			no_hit = np.where(shared_acc_counts == 0)
			
			#_qgak[q][0] is the query acc count
			possible_hits = np.minimum(_qgak[q][0].shape[0], _tct).astype(str)
			
			jaccard_averages = np.round(jaccard_averages, 4).astype(str)
			shared_acc_counts = shared_acc_counts.astype(str)
			
			jaccard_averages[no_hit] = "N/A"
			aai_ests[no_hit] = "N/A"
			shared_acc_counts[no_hit] = "N/A"
			possible_hits[no_hit] = "N/A"
			
			output_name = os.path.normpath(out_base + "/results/"+qname+"_results.txt")
		
			out = open(output_name, "w")
			out.write("query\ttarget\tavg_jacc_sim\tjacc_SD\tnum_shared_SCPs\tposs_shared_SCPs\tAAI_estimate\n")
			if do_sd:
				jaccard_SDs[no_hit] = "N/A"
				for i in range(0, len(aai_ests)):
					out.write(qname+"\t"+_tnames[i]+"\t"+jaccard_averages[i]+"\t"+jaccard_SDs[i]+"\t"+shared_acc_counts[i]+"\t"+possible_hits[i]+"\t"+aai_ests[i]+"\n")
			else:			
				for i in range(0, len(aai_ests)):
					out.write(qname+"\t"+_tnames[i]+"\t"+jaccard_averages[i]+"\t"+"N/A"+"\t"+shared_acc_counts[i]+"\t"+possible_hits[i]+"\t"+aai_ests[i]+"\n")
			out.close()
			
		else:
			aai_ests = numpy_kaai_to_aai_just_nums(jaccard_averages, as_float = True)
			aai_ests[no_hit] = 0
			print(*aai_ests, sep = "\t", file = outwriter)
		
		prog_queue.put(q)
		
	curs.close()
	prog_queue.put("done")
	
def on_disk_work_two(i):
	outwriter.close()
	return group_id

def sorted_nicely(l):
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
	
class db_db_remake:
	def __init__(self, in_memory = False, store_mat_res = False, 
		query = None, target = None, threads = 1, do_sd = False, 
		output_base = "FastAAI", output_style = "tsv", verbose = True):
		
		#databases to eat
		self.q = query
		self.t = target
		
		#metadata
		self.ok = generate_accessions_index(forward = True)
		self.rev = generate_accessions_index(forward = False)
		self.valids = None
		
		#Originally this was made to be a memory database only block of code, but just if/else one change makes it work on disk and it doesn't need a redev, then.
		self.as_mem_db = in_memory
		self.store_mat = store_mat_res
		
		#in-mem stuff
		self.conn = None
		self.curs = None
		
		self.threads = threads
		self.do_sd = do_sd

		self.output_base = output_base
		self.output = os.path.normpath(output_base + "/results")
		self.style = output_style
			
		self.query_names = None
		self.target_names = None
		
		self.num_queries = None
		self.num_targets = None
		
		self.query_gak = None
		self.target_gak = None
		self.target_presence = None
		
		self.query_dict = None
		self.target_dict = None
		
		self.verbose = verbose
	
	#getting the db metadata happens the same way in every case
	def open(self):
		if self.verbose:
			print("Perusing database metadata")
			
		self.conn = sqlite3.connect(":memory:")			
		self.curs = self.conn.cursor()
		
		self.curs.execute("attach '" + self.q + "' as queries")
		self.curs.execute("attach '" + self.t + "' as targets")
			
		#Find the shared accessions for these databases
		shared_accs_sql = '''
		SELECT queries.sqlite_master.name 
		FROM queries.sqlite_master INNER JOIN targets.sqlite_master 
		ON queries.sqlite_master.name = targets.sqlite_master.name
		'''
		self.valids = {}
		for table in self.curs.execute(shared_accs_sql).fetchall():
			table = table[0]
			#Filter to 
			if table in self.ok:
				self.valids[table] = self.ok[table]
		
		self.query_names = []
		for r in self.curs.execute("SELECT genome FROM queries.genome_index ORDER BY gen_id").fetchall():
			self.query_names.append(r[0])
		
		self.target_names = []
		for r in self.curs.execute("SELECT genome FROM targets.genome_index ORDER BY gen_id").fetchall():
			self.target_names.append(r[0])
			
		self.num_queries = len(self.query_names)
		self.num_targets = len(self.target_names)
		
		gak_sql = '''
		SELECT * FROM {db}.genome_acc_kmer_counts
		WHERE accession in ({accs})
		ORDER BY genome
		'''
		
		acc_ids = list(self.valids.values())
		acc_ids.sort()
		acc_ids = tuple(acc_ids)
		
		#query genome-acc-kmers (gak) is ordered by genome first, then accession
		self.query_gak = {}
		#for result in self.curs.execute(gak_sql.format(db = "queries", accs=','.join(['?']*len(self.valids))), acc_ids).fetchall():
		for result in self.curs.execute("SELECT * FROM queries.genome_acc_kmer_counts ORDER BY genome").fetchall():
			genome, accession, kmer_ct = result[0], result[1], result[2]
			if genome not in self.query_gak:
				self.query_gak[genome] = [[],[]]
			self.query_gak[genome][0].append(accession) 
			self.query_gak[genome][1].append(kmer_ct) 
			
		#refigure into numpy arrays for quicker array access later.
		for genome in self.query_gak:
			self.query_gak[genome] = (np.array(self.query_gak[genome][0], dtype = np.int32), np.array(self.query_gak[genome][1], dtype = np.int32))
		
		#Split these into ordered groups - this makes joining results at the end easier.
		qgak_queue = multiprocessing.Queue()
		groupings = split_seq_indices(np.arange(self.num_queries), self.threads)
		group_id = 0
		for group in groupings:
			next_set = {}
			for i in range(group[0], group[1]):
				next_set[i] = self.query_gak[i]
				self.query_gak[i] = None
			#this ensures that the selection of qgak and the query index range match
			qgak_queue.put((group_id, next_set, np.arange(group[0], group[1]),))
			group_id += 1
		
		self.query_gak = qgak_queue
		qgak_queue = None
		
		#tgt gak is organized by accession first, then genome
		self.target_gak = np.zeros(shape = (122, self.num_targets), dtype = np.int32)
		for result in self.curs.execute(gak_sql.format(db = "targets", accs=','.join(['?']*len(self.valids))), acc_ids).fetchall():
			genome, accession, kmer_ct = result[0], result[1], result[2]
			self.target_gak[accession, genome] += kmer_ct
	
		self.target_presence = self.target_gak > 0
		self.target_presence = self.target_presence.astype(bool)
	
	#This needs to have a TSV write method
	def load_in_mem(self):
		tempdir_path = os.path.normpath(self.output_base+"/temp")
		if not os.path.exists(tempdir_path):
			os.mkdir(tempdir_path)
			
		ql = {}
		tl = {}
		for t in self.valids.values():
			tl[t] = None
		for i in range(0, self.num_queries):
			ql[i] = {}
		
		if self.verbose:
			tracker = progress_tracker(total = len(self.valids), message = "Loading data in memory.")
		else:
			print("\nLoading data in memory.")
			
		pool = multiprocessing.Pool(self.threads, initializer = parse_db_init, initargs = (self.q, self.t, self.output_base,))
		
		for result in pool.imap_unordered(parse_accession, self.valids.keys()):
			this_accession = result[0]
			
			this_acc_id = self.ok[this_accession]
			
			with open(result[1], "rb") as inp:
				this_acc_data = pickle.load(inp)
			os.remove(result[1])
			
			tl[this_acc_id] = this_acc_data[1]
			
			for q in this_acc_data[0]:
				#We know that this acc must be in every ql for this loaded data.
				ql[q][this_acc_id] = this_acc_data[0][q][this_acc_id]
			if self.verbose:
				tracker.update()
			
		pool.close()
		
		if self.verbose:
			tracker = progress_tracker(total = self.num_queries, message = "Calculating AAI")
		else:
			print("\nCalculating AAI.")
		
		query_groups = []
		for grouping in split_seq_indices(np.arange(self.num_queries), self.threads):
			query_groups.append(np.arange(grouping[0], grouping[1]))
		
		result_queue = multiprocessing.Queue()
		remaining_procs = self.threads
		still_going = True
		
		pool = multiprocessing.Pool(self.threads, initializer = one_init, initargs = (ql, tl, self.num_targets, self.query_gak, self.target_gak, self.target_presence, self.do_sd, 
		self.style, self.output_base, self.store_mat, result_queue, self.query_names, self.target_names,))
		
		some_results = pool.imap(one_work, query_groups)
		
		while still_going:
			item = result_queue.get() 
			if item == "done":
				remaining_procs -= 1
				if remaining_procs == 0:
					still_going = False
			else:
				if self.verbose:
					tracker.update()
				else:
					pass
		
		if self.style == "matrix":
			result_files = []
			
			for result in pool.map(two_work, range(0, self.threads)):
				result_files.append(result)
			
			pool.close()
		
			self.write_mat_from_files(result_files)
		else:
			pool.close()

	#This needs to be implemented from existing code.
	def db_on_disk(self):				
		if self.style == "matrix":
			tempdir_path = os.path.normpath(self.output_base+"/temp")
			if not os.path.exists(tempdir_path):
				os.mkdir(tempdir_path)
			#Safety check - this is pointless on the normal version.
			self.store_mat = False
			
		result_queue = multiprocessing.Queue()
		remaining_procs = self.threads
		still_going = True
		
		if self.verbose:
			tracker = progress_tracker(total = self.num_queries, message = "Calculating AAI")
		else:
			print("\nCalculating AAI")
		
		query_groups = []
		for grouping in split_seq_indices(np.arange(self.num_queries), self.threads):
			query_groups.append(np.arange(grouping[0], grouping[1]))
		
		pool = multiprocessing.Pool(self.threads, initializer = on_disk_init, initargs = (self.q, self.t, 
		self.num_targets, self.query_gak, self.target_gak, self.target_presence, self.do_sd, 
		self.style, self.output_base, result_queue, self.query_names, self.target_names, self.valids))
		
		some_results = pool.imap(on_disk_work_one, query_groups)
		
		while still_going:
			item = result_queue.get() 
			if item == "done":
				remaining_procs -= 1
				if remaining_procs == 0:
					still_going = False
			else:
				if self.verbose:
					tracker.update()
				else:
					pass
		
		if self.style == "matrix":
			result_files = []
			for result in pool.map(on_disk_work_two, range(0, self.threads)):
				result_files.append(result)
			
		pool.close()
			
		if self.style == "matrix":
			self.write_mat_from_files(result_files)
			
			
	def write_mat_from_files(self, result_files):
		tempdir_path = os.path.normpath(self.output_base+"/temp")
	
		result_files = sorted_nicely(result_files)
		
		#print("Combining:")
		#for f in result_files:
		#	print(f)
		
		if self.verbose:
			tracker = progress_tracker(total = self.threads, step_size = 2, message = "Finalizing results.")
		else:
			print("\nFinalizing results.")
		
		output_file = os.path.normpath(self.output+"/FastAAI_matrix.txt")
		final_outwriter = open(output_file, "w")
		print("query_genome\t"+'\t'.join(self.target_names), file = final_outwriter)
		
		row = 0
			
		for f in result_files:
			fh = open(f, "r")
			cur = fh.readlines()
			fh.close()
			
			for i in range(0, len(cur)):
				if self.store_mat:
					#Add the decimals - we don't need to do this is we've been writing line-wise.
					#values will ALWAYS be 4 digits in this method, so groups of 2 dec. works.
					cur[i] = re.sub("(\d{2})(\d{2})", "\\1.\\2", cur[i])
				#Add in the query name to the row
				cur[i] = self.query_names[row]+"\t"+cur[i]
				row += 1
			
			final_outwriter.write(''.join(cur))
			cur = None	
			
			if self.verbose:
				tracker.update()
		
		final_outwriter.close()
		
		shutil.rmtree(tempdir_path)
		
	def close(self):
		self.curs.close()
		self.curs = None
		
	def clean_up(self):
		self.conn.close()
		self.conn = None
		
	def run(self):
		self.open()
		
		#work
		if self.as_mem_db:
			self.load_in_mem()
		else:
			self.db_on_disk()
		
		self.close()
		self.clean_up()
	
	
#Control the query process for any DB-first query.
def db_query(query, target, verbose, output, threads, do_stdev, style, in_mem, store_results):
	print("")
	
	#Sanity checks.
	if target is None:
		print("You need to supply a databasae for --target")
		sys.exit()
		
	#Sanity checks.
	if query is None:
		print("You need to supply a databasae for --query")
		sys.exit()
	

	
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
	
	#todo
	mdb = db_db_remake(in_memory = in_mem, store_mat_res = store_results, query = query, target = target, threads = threads, do_sd = do_stdev, output_base = output, output_style = style, verbose = verbose)
	mdb.run()
	
	print("")
	

#Check to see if the file exists and is a valid fastAAI db
def assess_db(path):
	status = None
	if os.path.exists(path):
		conn = sqlite3.connect(path)
		curs = conn.cursor()
		try:
			sql = "SELECT name FROM sqlite_master WHERE type='table'"
			
			curs.row_factory = lambda cursor, row: row[0]
			tables = curs.execute(sql).fetchall()
			curs.row_factory = None
			
			curs.close()
			conn.close()
						
			if len(tables) > 2 and "genome_index" in tables and "genome_acc_kmer_counts" in tables:		
				status = "exists"
			else:
				status = "wrong format"
			
		except:
			status = "wrong format"
		
	else:
		try:
			conn = sqlite3.connect(path)
			conn.close()
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
	
def merge_db_init(indexer, table_record, donor_dbs, tempdir):
	global mgi
	mgi = indexer
	global accs_per_db
	accs_per_db = table_record
	global tdb_list
	tdb_list = donor_dbs
	global work_space
	work_space = tempdir
	
def acc_transformer_merge(acc_name_genomes):
	acc_name = acc_name_genomes.split("_genomes")[0]
	my_acc_db = os.path.normpath(work_space + "/"+acc_name+".db")
	if os.path.exists(my_acc_db):
		os.remove(my_acc_db)
		
	my_db  = sqlite3.connect(my_acc_db)
	curs = my_db.cursor()
	curs.execute("CREATE TABLE {acc} (kmer INTEGER PRIMARY KEY, genomes array)".format(acc=acc_name))
	curs.execute("CREATE TABLE {acc} (genome INTEGER PRIMARY KEY, kmers array)".format(acc=acc_name_genomes))
	my_db.commit()
	
	reformat = {}
	for d in tdb_list:
		simple_rows = []
		#do nothing if the acc is not in the donor.
		if acc_name_genomes in accs_per_db[d]:
			donor_conn = sqlite3.connect(d)
			dcurs = donor_conn.cursor()
			data = dcurs.execute("SELECT * FROM {acc}".format(acc=acc_name_genomes)).fetchall()
			dcurs.close()
			donor_conn.close()
			
			for row in data:
				genome, kmers = row[0], row[1]
				new_index = mgi[d][genome]
				#-1 is the value indicating an already-seen genome that should not be added.
				if new_index > -1:
					simple_rows.append((new_index, kmers,))
					kmers = np.frombuffer(kmers, dtype=np.int32)
					for k in kmers:
						if k not in reformat:
							reformat[k] = []
						reformat[k].append(new_index)
			
			if len(simple_rows) > 0:
				curs.executemany("INSERT INTO {acc} VALUES (?,?)".format(acc=acc_name_genomes), simple_rows)
				my_db.commit()
				
			simple_rows = None
			data = None
	
	to_add = []
	for k in reformat:
		as_bytes = np.array(reformat[k], dtype = np.int32)
		as_bytes = as_bytes.tobytes()
		reformat[k] = None
		to_add.append((int(k), as_bytes,))
	
	curs.executemany("INSERT INTO {acc} VALUES (?, ?)".format(acc = acc_name), to_add)
	
	my_db.commit()
	
	to_add = None
	
	curs.execute("CREATE INDEX {acc}_index ON {acc} (kmer)".format(acc=acc_name))
	my_db.commit()
	
	curs.close()
	my_db.close()
	
	return [my_acc_db, acc_name]

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
		print("Donor databases:")
		for donor in valid_donors:
			print("\t", donor)
		print("Will be added to recipient database:", recipient)
	else:
		print("I couldn't find or create the recipient database at", recipient+".", "Does the folder you're trying to place this database in exist, and do you have permission to write files to it? FastAAI exiting.")
		sys.exit()
		
	if recipient is None or len(valid_donors) == 0:
		print("I require both a valid donor and a recipient database. FastAAI exiting.")
		sys.exit()
	
	gen_counter = 0
	multi_gen_ids = {}
	all_gens = {}
	
	#Load recipient data, if any.
	if recip_check == "exists":
		conn = sqlite3.connect(recipient)
		curs = conn.cursor()
		data = curs.execute("SELECT genome, gen_id FROM genome_index").fetchall()
		tabs = curs.execute("SELECT name FROM sqlite_master").fetchall()
		curs.close()
		conn.close()
		
		multi_gen_ids[recipient] = {}
		for row in data:
			genome, index = row[0], row[1]
			all_gens[genome] = 0
			multi_gen_ids[recipient][genome] = index
		
		gen_counter = max(list(multi_gen_ids[recipient].values())) + 1
				
	genome_index_to_add = []
	gak_to_add = []
	tables = {}
	#Donors should always exist, never be created.
	for d in valid_donors:
		#load
		conn = sqlite3.connect(d)
		curs = conn.cursor()
		data = curs.execute("SELECT * FROM genome_index").fetchall()
		tabs = curs.execute("SELECT name FROM sqlite_master").fetchall()
		gak = curs.execute("SELECT * FROM genome_acc_kmer_counts").fetchall()
		curs.close()
		conn.close()
		multi_gen_ids[d] = {}
		for row in data:
			genome, index, prot_ct = row[0], row[1], row[2]
			if genome not in all_gens:
				all_gens[genome] = 0
				#We need to be able to convert number to number.
				multi_gen_ids[d][index] = gen_counter
				genome_index_to_add.append((genome, gen_counter, prot_ct,))
				gen_counter += 1
			else:
				#This is a remove condition for later.
				multi_gen_ids[d][index] = -1
		data = None
		
		for row in gak:
			genome_id, acc_id, kmer_ct = row[0], row[1], row[2]
			new_index = multi_gen_ids[d][genome_id]
			if new_index > -1:
				gak_to_add.append((new_index, acc_id, kmer_ct,))
		
		tables[d] = []
		for tab in tabs:
			tab = tab[0]
			if tab.endswith("_genomes"):
				tables[d].append(tab)
		tables[d] = set(tables[d])
	
	all_tabs = set()
	for t in tables:
		all_tabs = all_tabs.union(tables[t])
	
	all_tabs = list(all_tabs)
	
	temp_dir = tempfile.mkdtemp()
	try:
		if verbose:
			tracker = progress_tracker(len(all_tabs), message = "Formatting data to add to database")
		else:
			print("Formatting data to add to database")
		
		conn = sqlite3.connect(recipient)
		curs = conn.cursor()
		
		#indexer, table_record, donor_dbs, tempdir
		pool = multiprocessing.Pool(threads, initializer=merge_db_init, initargs = (multi_gen_ids, tables, valid_donors, temp_dir,))

		for result in pool.imap_unordered(acc_transformer_merge, all_tabs):
			db, accession = result[0], result[1]
			curs.execute("CREATE TABLE IF NOT EXISTS {acc} (kmer INTEGER PRIMARY KEY, genomes array)".format(acc=accession))
			curs.execute("CREATE TABLE IF NOT EXISTS {acc}_genomes (genome INTEGER PRIMARY KEY, kmers array)".format(acc=accession))
			curs.execute("CREATE INDEX IF NOT EXISTS {acc}_index ON {acc}(kmer)".format(acc=accession))
			conn.commit()
			
			curs.execute("attach '" + db + "' as acc")
			conn.commit()
			
			#Get the genomes from worker db.
			curs.execute("INSERT INTO {acc}_genomes SELECT * FROM acc.{acc}_genomes".format(acc=accession))
			to_update = curs.execute("SELECT kmer, genomes, genomes FROM acc.{acc}".format(acc=accession)).fetchall()
			update_concat_sql = "INSERT INTO {acc} VALUES (?,?) ON CONFLICT(kmer) DO UPDATE SET genomes=genomes || (?)".format(acc=accession)
			curs.executemany(update_concat_sql, to_update)
			conn.commit()
			
			curs.execute("detach acc")
			conn.commit()
			
			os.remove(db)
			
			if verbose:
				tracker.update()
		
		pool.close()
		pool.join()
		
		curs.execute("CREATE TABLE IF NOT EXISTS genome_index (genome text, gen_id integer, protein_count integer)")
		curs.execute("CREATE TABLE IF NOT EXISTS genome_acc_kmer_counts (genome integer, accession integer, count integer)")
		
		curs.executemany("INSERT INTO genome_index VALUES (?,?,?)", genome_index_to_add)
		curs.executemany("INSERT INTO genome_acc_kmer_counts VALUES (?,?,?)", gak_to_add)
		
		curs.execute("CREATE INDEX IF NOT EXISTS kmer_acc ON genome_acc_kmer_counts (genome, accession);")
		
		conn.commit()
			
	except:
		curs.close()
		conn.close()
		#Error
		shutil.rmtree(temp_dir)
	finally:
		curs.close()
		conn.close()
		#Success
		shutil.rmtree(temp_dir)
		
	valid_accs = generate_accessions_index()
	
	
	print("\nDatabases merged!")
	
	return None
	
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
	parser.add_argument('--compress', dest = "do_comp", action = 'store_true', help = 'Gzip compress generated proteins, HMMs. Off by default.')
	
	args, unknown = parser.parse_known_args()
	
	return parser, args

def kaai_to_aai(kaai):
	# Transform the kAAI into estimated AAI values
	aai_hat = (-0.3087057 + 1.810741 * (np.exp(-(-0.2607023 * np.log(kaai))**(1/3.435))))*100
	
	return aai_hat
	
#This one's unique. It doesn't do anything with the DB, which means it doesn't access any other functionality outside of the input_file class. It just advances a pair of inputs in parallel and does intersections.
def single_query(qf, tf, output, verbose, threads, do_compress):	
	
	if qf.identifiers[0] == tf.identifiers[0]:
		print("You've selected the same query and target genome. The AAI is 100%.")
		print("FastAAI exiting.")
		return None
	
	statuses = ["genome", "protein", "protein and hmm"]
	query_stat = statuses.index(qf.status)
	target_stat = statuses.index(tf.status)
	minimum_status = statuses[min(query_stat, target_stat)]
	
	start_printouts = ["[Genome] Protein Protein+HMM", " Genome [Protein] Protein+HMM", "Genome Protein [Protein+HMM]"]
	
	print("")
	print("Query start: ", start_printouts[query_stat])
	print("Target start:", start_printouts[target_stat])
	print("")
	
	
	qname = qf.identifiers[0]
	tname = tf.identifiers[0]
	
	name = os.path.normpath(output + "/results/" + qname + "_vs_" + tname + ".aai.txt")
	print("Output will be located at", name)
	
	advance_me = [qf.in_files[0], tf.in_files[0]]
	#All we need to do this.
	hmm_file = find_hmm()
	pool = multiprocessing.Pool(min(threads, 2), initializer = hmm_preproc_initializer, initargs = (hmm_file, do_compress,))
	
	results = pool.map(run_build, advance_me)
	
	pool.close()
	pool.join()
	
	query = results[0]
	target = results[1]
	
	print(query.partial_timings())
	print(target.partial_timings())
	
	#One of the printouts
	max_poss_prots = max(len(query.best_hits_kmers), len(target.best_hits_kmers))
	
	accs_to_view = set(query.best_hits_kmers.keys()).intersection(set(target.best_hits_kmers.keys()))
	
	results = []
	for acc in accs_to_view:
		intersect = np.intersect1d(query.best_hits_kmers[acc], target.best_hits_kmers[acc])
		intersect = intersect.shape[0]
		union = query.best_hits_kmers[acc].shape[0] + target.best_hits_kmers[acc].shape[0] - intersect
		jacc = intersect/union
		results.append(jacc)
		
	results = np.array(results, dtype = np.float_)
	
	jacc_mean = np.mean(results)
	jacc_std = np.std(results)
	actual_prots = len(results)
	poss_prots = max(len(query.best_hits_kmers), len(target.best_hits_kmers))
	aai_est = round(kaai_to_aai(jacc_mean), 2)
	
	if aai_est > 90:
		aai_est = ">90%"
	else:
		if aai_est < 30:
			aai_est = "<30%"
	
	output = open(name, "w")
	
	print("query\ttarget\tavg_jacc_sim\tjacc_SD\tnum_shared_SCPs\tposs_shared_SCPs\tAAI_estimate", file = output)
	print(qname, tname, round(jacc_mean, 4), round(jacc_std, 4), actual_prots, poss_prots, aai_est, sep = "\t", file = output)
	
	output.close()
	
	print("query\ttarget\tavg_jacc_sim\tjacc_SD\tnum_shared_SCPs\tposs_shared_SCPs\tAAI_estimate")
	print(qname, tname, round(jacc_mean, 4), round(jacc_std, 4), actual_prots, poss_prots, aai_est, sep = "\t")

	
	print("FastAAI single query done! Estimated AAI:", aai_est)

def miga_merge_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	Hello, Miguel.
	
	Give one genome in nt, aa, or aa+hmm format and a database to create or add to. 
	It'll add the genome as efficiently as possible.
	
	The normal merge command creates parallel processes and gathers data in 
	one-SCP databases to add to the main DB. Great for many genomes. A lot of extra 
	work for just one.
	
	This version skips the creation of subordinate DBs and just directly adds the genome. 
	Faster, fewer writes, no parallel overhead. Maybe we can even daemonize this one.
	''')
			
	parser.add_argument('--genome',      dest = 'gen',    default = None,   help =  'Path to one genome, FASTA format')
	parser.add_argument('--protein',      dest = 'prot',    default = None,   help =  'Path to one protein, AA FASTA format')
	parser.add_argument('--hmm',      dest = 'hmm',    default = None,   help =  'Path to one HMM file as predicted by FastAAI')

	parser.add_argument('--output', dest = 'output', default = "FastAAI", help = 'Place the partial output files into a directory with this base. Default "FastAAI"')
	parser.add_argument('--target',   dest = 'database', default = None, help =  'Path to the target database. The genome supplied will be added to this. The DB will be created if needed.')	
		
	parser.add_argument('--verbose',           dest = 'verbose',   action='store_true', help = 'Print minor updates to console. Major updates are printed regardless.')
	parser.add_argument('--compress',           dest = 'compress',   action='store_true', help = 'Compress generated file output')
	
	args, unknown = parser.parse_known_args()
	
	return parser, args

def miga_merge(infile, target_db, verbose, do_compress):
	status = assess_db(target_db)
	if status == "wrong format":
		print("The database", target_db, "exists, but appears to not be a FastAAI database.")
		print("FastAAI will not alter this file. Quitting.")
		return None
		
	if status == "unable to create":
		print("The database", target_db, "could not be created.")
		print("Are you sure that the path you gave is valid? Quitting.")
		return None

	if verbose:
		print("Processing genome")
	
	next_id = 0
	exist_gens = {}
	conn = sqlite3.connect(target_db)
	curs = conn.cursor()
	if status == 'exists':
		for row in curs.execute("SELECT * FROM genome_index ORDER BY gen_id").fetchall():
			genome, id, prot_ct = row[0], row[1], row[2]
			exist_gens[genome] = id
			next_id += 1
			
	if infile.basename in exist_gens:
		print("It looks like the file you're trying to add already exists in the database.")
		print("Adding it is too likely to corrupt the database. Quitting.")
		return None

	hmm_file = find_hmm()
	global hmm_manager

	hmm_manager = pyhmmer_manager(do_compress)
	hmm_manager.load_hmm_from_file(hmm_file)
	
	infile.preprocess()
	
	if len(infile.best_hits_kmers) > 0:
	
		ok = generate_accessions_index()
		gak_to_add = []
		
		gen_id = np.zeros(1, dtype = np.int32)
		gen_id[0] = next_id
		gen_id = gen_id.tobytes()
		
		for accession in infile.best_hits_kmers:
			acc_id = ok[accession]
			gak_to_add.append((next_id, acc_id, infile.best_hits_kmers[accession].shape[0],))
			
			curs.execute("CREATE TABLE IF NOT EXISTS {acc} (kmer INTEGER PRIMARY KEY, genomes array)".format(acc=accession))
			curs.execute("CREATE TABLE IF NOT EXISTS {acc}_genomes (genome INTEGER PRIMARY KEY, kmers array)".format(acc=accession))
			curs.execute("CREATE INDEX IF NOT EXISTS {acc}_index ON {acc}(kmer)".format(acc=accession))
			
			gen_first = (next_id, infile.best_hits_kmers[accession].tobytes(),)
			curs.execute("INSERT INTO {acc}_genomes VALUES (?,?)".format(acc=accession), gen_first)
			
			kmers_first = []
			for k in infile.best_hits_kmers[accession]:
				#we know there's only one genome in these cases.
				kmers_first.append((int(k), gen_id, gen_id, ))
			
			update_concat_sql = "INSERT INTO {acc} VALUES (?,?) ON CONFLICT(kmer) DO UPDATE SET genomes=genomes || (?)".format(acc=accession)
			
			curs.executemany(update_concat_sql, kmers_first)
			
		#Safety checks.
		curs.execute("CREATE TABLE IF NOT EXISTS genome_index (genome text, gen_id integer, protein_count integer)")
		curs.execute("CREATE TABLE IF NOT EXISTS genome_acc_kmer_counts (genome integer, accession integer, count integer)")
		
		gen_idx_to_add = (infile.basename, next_id, len(infile.best_hits_kmers))
		curs.execute("INSERT INTO genome_index VALUES (?, ?, ?)", gen_idx_to_add)
		#gak was made over the loops.
		curs.executemany("INSERT INTO genome_acc_kmer_counts VALUES (?,?,?)", gak_to_add)
		curs.execute("CREATE INDEX IF NOT EXISTS kmer_acc ON genome_acc_kmer_counts (genome, accession);")

		conn.commit()	
		
	else:
		print("No proteins to add for this genome:",infile.basename,"Database will be unaltered. Exiting.")
	
	curs.close()
	conn.close()
	
'''
Main
'''
def main():
	#The currently supported modules.
	modules = ["build_db", "merge_db", "simple_query", "db_query", "single_query", "miga_merge"]
	
	#Print modules if someone just types FastAAI
	if len(sys.argv) < 2:
		print("")
		print("    I couldn't find the module you specified. Please select one of the following modules:")
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
		print("------------------------------------------- Other Options -------------------------------------------")
		print("")
		print("    single_query |" + " Query ONE query genome against ONE target genome")
		#print("    multi_query  |" + " Create a query DB and a target DB, then calculate query vs. target AAI")
		#print("    aai_index    |" + " Create a database from multiple genomes and do an all vs. all AAI index of the genomes")
		print("")
		print("-----------------------------------------------------------------------------------------------------------")
		print("    To select a module, enter 'FastAAI [module]' into the command line!")
		print("")
		sys.exit()
		
	#This is the module selection
	selection = sys.argv[1]
	
	if selection not in modules:
		print("")
		print("    I couldn't find the module you specified. Please select one of the following modules:")
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
		print("------------------------------------------- Other Options -------------------------------------------")
		print("")
		print("    single_query |" + " Query ONE query genome against ONE target genome")
		#print("    multi_query  |" + " Create a query DB and a target DB, then calculate query vs. target AAI")
		#print("    aai_index    |" + " Create a database from multiple genomes and do an all vs. all AAI index of the genomes")
		print("")
		print("-----------------------------------------------------------------------------------------------------------")
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
		
		output  = os.path.normpath(opts.output)
	
		threads = opts.threads
		verbose = opts.verbose
		
		#Database handle
		db_name = opts.db_name
		
		do_comp = opts.do_comp
		
		build_db(genomes, proteins, hmms, db_name, output, threads, verbose, do_comp)
	
		
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
			
		genomes, proteins, hmms = opts.genomes, opts.proteins, opts.hmms
		
		db_name = opts.target
		
		output  = opts.output
		threads = opts.threads
		verbose = opts.verbose
		
		do_stdev = opts.do_stdev
		
		style, in_mem, make_db, qdb_name, do_comp = opts.style, opts.in_mem, opts.make_db, opts.qdb_name, opts.do_comp
				
		sql_query(genomes, proteins, hmms, db_name, output, threads, verbose, do_stdev, style, in_mem, make_db, qdb_name, do_comp)
		
			
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
		output = opts.output
		threads = opts.threads
		
		style, in_mem, store = opts.style, opts.in_mem, opts.storage
			
		
		db_query(query, target, verbose, output, threads, do_stdev, style, in_mem, store)

#################### One-pass functions #######################
	if selection == "single_query":
		parser, opts = single_query_opts()
		#module name only
		
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
			
		output = os.path.normpath(opts.output)
		try:
			threads = int(opts.threads)
		except:
			print("Couldn't interpret your threads. Defaulting to 1.")
			threads = 1
		verbose = opts.verbose
		do_compress = opts.do_comp
		
		query_genome = opts.query_genome
		query_protein = opts.query_protein
		query_hmm = opts.query_hmm
		
		query_file = fastaai_file_importer(genomes = query_genome, proteins = query_protein, hmms = query_hmm, output = output, compress = do_compress)
		query_file.determine_inputs()
		
		target_genome = opts.target_genome
		target_protein = opts.target_protein
		target_hmm = opts.target_hmm
		
		target_file = fastaai_file_importer(genomes = target_genome, proteins = target_protein, hmms = target_hmm, output = output, compress = do_compress)
		target_file.determine_inputs()
		
		is_ok = True
		if len(query_file.in_files) != 1:
			print("Query genome unacceptable. Check your inputs")
			is_ok = False
		
		if len(target_file.in_files) != 1:
			print("target genome unacceptable. Check your inputs")
			is_ok = False
		if is_ok:
			good_to_go = prepare_directories(output, query_file.status, "query")
			if good_to_go:
				good_to_go = prepare_directories(output, target_file.status, "query")
			if good_to_go:
				single_query(query_file, target_file, output, verbose, threads, do_compress)
		
############## MiGA module #################
	if selection == "miga_merge":
		parser, opts = miga_merge_opts()
		
		#module name only
		if len(sys.argv) < 3:
			print(parser.print_help())
			sys.exit()
		
		g,p,h = opts.gen, opts.prot, opts.hmm
		
		target = opts.database
		
		verbose = opts.verbose
		
		output_path = opts.output
		
		if target == None:
			target = os.path.normpath(output_path + "/database/FastAAI_database.sqlite.db")
		
		do_compress = opts.compress
		
		imported_files = fastaai_file_importer(genomes = g, proteins = p, hmms = h, 
		output = output_path, compress = do_compress)
		
		imported_files.determine_inputs()
		
		if len(imported_files.in_files) == 0:
			print("Something was wrong with your input file.")
		else:
			input_genome = imported_files.in_files[0]
			
			good_to_go = prepare_directories(output_path, imported_files.status, "build")
			
			miga_merge(input_genome, target, verbose, do_compress)
		
			#This is where a new db would normally be created, 
			#which is not what happens when the supplied target is some other sort of path.
			output_default = os.path.normpath(output_path + "/database")
			if len(os.listdir(output_default)) == 0:
				os.rmdir(output_default)
	
	return None
	
if __name__ == "__main__":
	main()

	