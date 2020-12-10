#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  0.8
# Date:		 March 02, 2020

# Description: Calculates the average amino acid identity using k-mers
from single copy genes. It is a faster version of the regular AAI (Blast
or Diamond) and the hAAI implemented in MiGA.
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
from functools import partial


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
    folder = file_path.parent
    protein_output = folder / (filename + '.faa')
    output_11 = folder / (filename + '.faa.11')
    temp_output = folder / (filename + '.temp')
    subprocess.call(["prodigal", "-i", str(file_path), "-a", str(output_11), 
                    "-p", "meta", "-q", "-o", str(temp_output)])
    output_4 = folder / (filename + '.faa.4')
    temp_output = folder / (filename + '.temp')
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
    
    if (length_4 / length_11) >= 1.1:
        shutil.copy(output_4, protein_output)
    else:
        shutil.copy(str(output_11), str(protein_output))
    
    # Remove intermediate files
    output_4.unlink()
    output_11.unlink()
    temp_output.unlink()

    # Remove stop '*' codons from protein sequences
    with open(protein_output, 'r') as final_protein, open(temp_output, 'w') as temporal_file:
        for line in final_protein:
            if line.startswith(">"):
                temporal_file.write("{}".format(line))
            else:
                line = line.replace('*', '')
                temporal_file.write("{}".format(line))
    shutil.copy(str(temp_output), str(protein_output))
    temp_output.unlink()

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
    folder = file_path.parent
    protein_output = folder / (filename + '.faa')
    temp_output = folder / (filename + '.temp')
    subprocess.call(["prodigal", "-i", str(file_path), "-a", str(protein_output), 
                    "-p", "meta", "-q", "-o", str(temp_output)])
    
    # Remove intermediate files
    temp_output.unlink()

    # Remove stop '*' codons from protein sequences
    with open(protein_output, 'r') as final_protein, open(temp_output, 'w') as temporal_file:
        for line in final_protein:
            if line.startswith(">"):
                temporal_file.write("{}".format(line))
            else:
                line = line.replace('*', '')
                temporal_file.write("{}".format(line))
    shutil.copy(str(temp_output), str(protein_output))
    temp_output.unlink()

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
    folder = file_path.parent
    name = file_path.name
    hmm_output = folder / (name + '.hmm')
    temp_output = folder / (name + '.temp')
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
def hmm_filter(scg_hmm_file, keep):
    """
    Filters HMM results for best hits per protein
    
    Arguments:
        SCG_HMM_file {file path} -- Path to HMM results file
        keep {bool} -- Keep HMM files
    
    Returns:
        outfile -- Path to filtered files
    """
    hmm_path = Path(scg_hmm_file)
    name = hmm_path.name
    folder = hmm_path.parent
    outfile = folder / (name + '.filt')
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
    with open(protein_file) as fasta_in:
        for line in fasta_in:
            if line.startswith(">"):
                if store_sequence == True:
                    kmers = build_kmers(protein_sequence, kmer_size)
                    kmers = set(kmers.split(","))
                    scg_kmers.update(kmers)
                    protein_sequence = ""
                else:
                    protein_sequence = ""
                    store_sequence = True
            else:
                protein_sequence += line.strip()
    genome_kmers = {final_filename : list(scg_kmers)}
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

# --- Convert kmers to indices ---
# ------------------------------------------------------
def convert_kmers_to_indices(kmer_dict):
    print("Converting kmers to indices")
    for genome in kmer_dict:
        for protein_marker in kmer_dict[genome]:
            kmer_index = []
            for kmer in kmer_dict[genome][protein_marker].split(','):
                kmer_index.append(global_kmer_index_dictionary[kmer])
            kmer_index = np.sort(np.unique(np.array(kmer_index, dtype=np.int32)))
            kmer_dict[genome][protein_marker] = kmer_index

    return kmer_dict
# ------------------------------------------------------

# --- Create global dictionary with unique kmers and indices for each one ---
# ------------------------------------------------------
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

# --- Transform kmer dictionaries to index dictionaries ---
# ------------------------------------------------------
def transform_kmer_dicts_to_arrays(kmer_dict, temporal_working_directory, single_dataset):
    kmer_dict = convert_kmers_to_indices(kmer_dict)
    #Get skip indices
    smartargs = []
    genome_ids = list(kmer_dict.keys())
    for i in range(0, len(genome_ids)):
        if single_dataset == True:
            smartargs.append((temporal_working_directory, genome_ids[i], i))
        else:
            smartargs.append((temporal_working_directory, genome_ids[i]))
        
    return kmer_dict, smartargs
# ------------------------------------------------------

# --- Parse kAAI when query == reference ---
# ------------------------------------------------------
def single_kaai_parser(arguments):
    """
    Calculates the Jaccard distances using single protein markers shared by two genomes
    
    Arguments:
        arguments {tuple} -- Tuple with the temporal folder, the query id and the index of said query_id
    
    Returns:
        [Path to output] -- Path to output file
    """
    temporal_folder = arguments[0]
    query_id = arguments[1]
    skip_first_n = arguments[2]
    
    temporal_folder = Path(str(temporal_folder.name))
    temporal_file = Path(query_id).name + '.aai.temp'
    temporal_output = temporal_folder / temporal_file
    
    query_scg_list = np.array(list(query_kmer_dictionary[query_id].keys()))

    with open(temporal_output, 'w') as out_file:
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
            
            # Allow for numpy in-builts; they're a little faster.
            jaccard_similarities = np.array(jaccard_similarities, dtype=np.float_)
            try:
                mean = np.mean(jaccard_similarities)
                var = np.std(jaccard_similarities)
                if mean >= 0.9:
                    aai_est = ">90%"
                else:
                    aai_est = kaai_to_aai(mean)
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           round(mean, 4), round(var, 4),
                           len(jaccard_similarities), shorter_genome, aai_est))
            except:
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           "NA", "NA", "NA", "NA", "NA"))
    return temporal_output
# ------------------------------------------------------

# --- Parse viral kAAI when query == reference ---
# ------------------------------------------------------
def single_virus_kaai_parser(query_id):
    """
    Calculates Jaccard distances on kmers from viral proteins
    
    Arguments:
        query_id {str} -- Id of the query genome
    
    Returns:
        [Path to output] -- Path to output file
    """
    file_path = Path(query_id)
	
	#Carlos, tempdir for safety
    tmp_folder = tempfile.TemporaryDirectory()
    running_folder = tmp_folder.name
	
	
    temp_output = running_folder / file_path.with_suffix('.aai.temp')
    # Start comparison with all genomes in the query dictionary
    with open(temp_output, 'w') as out_file:
        for target_genome, kmers_target in query_kmer_dictionary.items():
            jaccard_index = None
            kmers_query = set(query_kmer_dictionary[query_id])
            intersection = len(kmers_query.intersection(kmers_target))
            union = len(kmers_query.union(kmers_target))
            try:
                jaccard_index = intersection / union
                out_file.write("{}\t{}\t{}\n".format(query_id, target_genome, jaccard_index))
            except:
                out_file.write("{}\t{}\tNA\n".format(query_id, target_genome))
    return temp_output
# ------------------------------------------------------

# --- Parse kAAI when query != reference ---
# ------------------------------------------------------
def double_kaai_parser(arguments):
    """
    Calculates the Jaccard distances using single protein markers shared by two genomes
    
    Arguments:
        arguments {tuple} -- Tuple with the temporal folder, the query id and the index of said query_id
    
    Returns:
        [Path to output] -- Path to output file
    """
    temporal_folder = arguments[0]
    query_id = arguments[1]
    
    temporal_folder = Path(str(temporal_folder.name))
    temporal_file = Path(query_id).name + '.aai.temp'
    temporal_output = temporal_folder / temporal_file
    
    query_scg_list = np.array(list(query_kmer_dictionary[query_id].keys()))

    with open(temporal_output, 'w') as out_file:
        #for target_genome, scg_ids in query_kmer_dictionary.items():
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
            jaccard_similarities = np.array(jaccard_similarities, dtype=np.float_)
            try:
                mean = np.mean(jaccard_similarities)
                var = np.std(jaccard_similarities)
                if mean >= 0.9:
                    aai_est = ">90%"
                else:
                    aai_est = kaai_to_aai(mean)
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           round(mean, 4), round(var, 4),
                           len(jaccard_similarities), shorter_genome, aai_est))
            except:
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           "NA", "NA", "NA", "NA", "NA"))
    return temporal_output


    """
    Calculates Jaccard distances on kmers from proteins shared
    
    Arguments:
        query_id {str} -- Id of the query genome
    
    Returns:
        [Path to output] -- Path to output file
    """
    file_path = Path(query_id)
	
	#Carlos, tempdir for safety
    tmp_folder = tempfile.TemporaryDirectory()
    running_folder = tmp_folder.name
	
	
    temp_output = running_folder / file_path.with_suffix('.aai.temp')
    # Get number and list of SCG detected in query
    query_num_scg = len(query_kmer_dictionary[query_id])
    query_scg_list = query_kmer_dictionary[query_id].keys()
    # Start comparison with all genomes in the query dictionary
    with open(temp_output, 'w') as out_file:
        for target_genome, scg_ids in ref_kmer_dictionary.items():
            jaccard_similarities = []
            # Get number and list of SCG detected in reference
            target_num_scg = len(scg_ids)
            target_scg_list = scg_ids.keys()
            # Choose the smallest set of proteins
            if query_num_scg > target_num_scg:
                final_scg_list = target_scg_list
            else:
                final_scg_list = query_scg_list
            # Compare all the proteins in the final SCG list
            for accession in final_scg_list:
                if accession in query_scg_list and accession in target_scg_list:
                    # Get set and list for each SCG accession
                    kmers_query = set(query_kmer_dictionary[query_id][accession].split(','))
                    kmers_target = ref_kmer_dictionary[target_genome][accession].split(',')
                    # Calculate jaccard_similarity
                    intersection = len(kmers_query.intersection(kmers_target))
                    union = len(kmers_query.union(kmers_target))
                    jaccard_similarities.append(intersection / union)
                else:
                    continue
            try:
                n = len(jaccard_similarities)
                mean = sum(jaccard_similarities)/n
                var = sum([ (x - mean)**2 for x in jaccard_similarities ])/(n - 1)
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           round(mean, 4), round(var**0.5, 4),
                           len(jaccard_similarities), len(final_scg_list)))
            except:
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           "NA", "NA", "NA", "NA"))
    return temp_output
# ------------------------------------------------------

# --- Parse viral kAAI when query != reference ---
# ------------------------------------------------------
def double_viral_kaai_parser(query_id):
    """
    Calculates Jaccard distances on kmers from viral proteins
    
    Arguments:
        query_id {str} -- Id of the query genome
    
    Returns:
        [Path to output] -- Path to output file
    """
    file_path = Path(query_id)
	
	#Carlos, tempdir for safety
    tmp_folder = tempfile.TemporaryDirectory()
    running_folder = tmp_folder.name
	
	
    temp_output = running_folder / file_path.with_suffix('.aai.temp')
    # Start comparison with all genomes in the query dictionary
    with open(temp_output, 'w') as out_file:
        for target_genome, kmers_target in ref_kmer_dictionary.items():
            jaccard_index = None
            kmers_query = set(query_kmer_dictionary[query_id])
            intersection = len(kmers_query.intersection(kmers_target))
            union = len(kmers_query.union(kmers_target))
            try:
                jaccard_index = intersection / union
                out_file.write("{}\t{}\t{}\n".format(query_id, target_genome, jaccard_index))
            except:
                out_file.write("{}\t{}\tNA\n".format(query_id, target_genome))
    return temp_output
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

# --- Merge kmer dictionaries ---
# ------------------------------------------------------
def kaai_to_aai(kaai):
    # Transform the kAAI into estimated AAI values
    aai_hat = (-0.3087057 + 1.810741 * (np.exp(-(-0.2607023 * np.log(kaai))**(1/3.435))))*100
    return aai_hat
# ------------------------------------------------------
	

################################################################################
"""---2.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script calculates the average amino acid identity using k-mers\n'''
                        '''from single copy genes. It is a faster version of the regular AAI '''
                        '''(Blast or Diamond) and the hAAI implemented in MiGA.'''
            '''Usage: ''' + argv[0] + ''' -p [Protein Files] -t [Threads] -o [Output]\n'''
            '''Global mandatory parameters: -g [Genome Files] OR -p [Protein Files] OR -s [SCG HMM Results] -o [AAI Table Output]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')
    mandatory_options = parser.add_argument_group('Mandatory i/o options. You must select an option for the queries and one for the references.')
    mandatory_options.add_argument('--qg', dest='query_genomes', action='store', required=False,
                                    help='File with list of query genomes.')
    mandatory_options.add_argument('--qp', dest='query_proteins', action='store', required=False,
                                    help='File with list of query proteins.')
    mandatory_options.add_argument('--qh', dest='query_hmms', action='store', required=False,
                                    help=textwrap.dedent('''
                                    File with list of pre-computed query hmmsearch results.
                                    If you select this option you must also provide a file with 
                                    a list of protein files for the queries (with --qp).
                                    '''))
    mandatory_options.add_argument('--qd', dest='query_database', action='store', required=False,
                                    help='File with list of pre-indexed query databases.')
    mandatory_options.add_argument('--rg', dest='reference_genomes', action='store', required=False,
                                    help='File with list of reference genomes.')
    mandatory_options.add_argument('--rp', dest='reference_proteins', action='store', required=False,
                                    help='File with list of reference proteins.')
    mandatory_options.add_argument('--rh', dest='reference_hmms', action='store', required=False,
                                    help=textwrap.dedent('''
                                    File with list of pre-computed reference hmmsearch results.
                                    If you select this option you must also provide a file with 
                                    a list of protein files for the references (with --qp).
                                    '''))
    mandatory_options.add_argument('--rd', dest='reference_database', action='store', required=False,
                                    help='File with list of pre-indexed reference databases.')
    mandatory_options.add_argument('-o', '--output', dest='output', action='store', required=False, help='Output file. By default kaai_comparisons.txt')
    additional_input_options = parser.add_argument_group('Behavior modification options.')
    additional_input_options.add_argument('-e', '--ext', dest='extension', action='store', required=False, 
                                            help='Extension to remove from original filename, e.g. ".fasta"')
    additional_input_options.add_argument('-i', '--index', dest='index_db', action='store_true', required=False, 
                                            help='Only index and store databases, i.e., do not perform comparisons.')
    misc_options = parser.add_argument_group('Miscellaneous options')
    misc_options.add_argument('--virus', dest='virus', action='store_true', required=False,
                                help='Toggle virus-virus comparisons. Use only with viral genomes or proteins.')
    misc_options.add_argument('-t', '--threads', dest='threads', action='store', default=1, type=int, required=False,
                                help='Number of threads to use, by default 1')
    misc_options.add_argument('-k', '--keep', dest='keep', action='store_false', required=False,
                                help='Keep intermediate files, by default true')

    args = parser.parse_args()

    query_genomes = args.query_genomes
    reference_genomes = args.reference_genomes
    query_proteins = args.query_proteins
    reference_proteins = args.reference_proteins
    query_hmms = args.query_hmms
    reference_hmms = args.reference_hmms
    query_database = args.query_database
    reference_database = args.reference_database
    output = args.output
    if output == None:
        output == "kaai_comparisons.txt"
    extension = args.extension
    index_db = args.index_db
    threads = args.threads
    keep = args.keep
    virus = args.virus

    print("kAAI started on {}".format(datetime.datetime.now()))
    # Check user input
    # ------------------------------------------------------
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
    # ------------------------------------------------------

    # Create temporal working directory
    temporal_working_directory = TemporaryDirectory()
    # ------------------------------------------------------

    # Check if queries are the same as references (an all-vs-all comparison)
    # ------------------------------------------------------
    same_inputs = False
    if query_input == reference_input:
        same_inputs = True
    if same_inputs == True:
        print('You specified the same query and reference files.')
        print('I will perform an all vs all comparison :)')
    # ------------------------------------------------------

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
    # If starting from database and query == reference
    if same_inputs == True:
        if query_database != None:
            with open(query_database) as query_database_files:
                for db_location in query_database_files:
                    if Path(db_location.strip()).is_file():
                        with gzip.open(db_location.strip(), 'rb') as database_handle:
                            temp_dict = pickle.load(database_handle)
                            if isinstance(temp_dict,dict):
                                query_kmer_dict_list.append(temp_dict)
								#Carlos, this line serves no purpose but does take a bunch of time and mem.
                                #print(query_kmer_dict_list)
                            else:
                                exit("One of the database files appear to have the wrong format. Please provide a correctly formated databases.")
            query_kmer_dict = merge_dicts(query_kmer_dict_list)
    else:
    # If the inputs are not the same:
        # If query and ref are provided
        if query_database != None and reference_database != None:
            with open(query_database, 'r') as query_database_files:
                for db_location in query_database_files:
                    if Path(db_location.strip()).is_file():
                        with gzip.open(db_location.strip(), 'rb') as database_handle:
                            temp_dict = pickle.load(database_handle)
                            if isinstance(temp_dict,dict):
                                query_kmer_dict_list.append(temp_dict)
                            else:
                                exit("One of the query database files appear to have the wrong format. Please provide a correctly formated databases.")
            query_kmer_dict = merge_dicts(query_kmer_dict_list)
            with open(reference_database) as reference_database_files:
                for db_location in reference_database_files:
                    if Path(db_location.strip()).is_file():
                        with gzip.open(db_location.strip(), 'rb') as database_handle:
                            temp_dict = pickle.load(database_handle)
                            if isinstance(temp_dict,dict):
                                reference_kmer_dict_list.append(temp_dict)
                            else:
                                exit("One of the reference database files appear to have the wrong format. Please provide a correctly formated databases.")
            reference_kmer_dict = merge_dicts(reference_kmer_dict_list)
        # If only the query has a db
        elif query_database != None and reference_database == None:
            with open(query_database) as query_database_files:
                for db_location in query_database_files:
                    if Path(db_location.strip()).is_file():
                        with gzip.open(db_location.strip(), 'rb') as database_handle:
                            temp_dict = pickle.load(database_handle)
                            if isinstance(temp_dict,dict):
                                query_kmer_dict_list.append(temp_dict)
                            else:
                                exit("One of the query database files appear to have the wrong format. Please provide a correctly formated databases.")
            query_kmer_dict = merge_dicts(query_kmer_dict_list)
        # If only the reference has a db
        elif query_database == None and reference_database != None:
            with open(reference_database) as reference_database_files:
                for db_location in reference_database_files:
                    if Path(db_location.strip()).is_file():
                        with gzip.open(db_location.strip(), 'rb') as database_handle:
                            temp_dict = pickle.load(database_handle)
                            if isinstance(temp_dict,dict):
                                reference_kmer_dict_list.append(temp_dict)
                            else:
                                exit("One of the reference database files appear to have the wrong format. Please provide a correctly formated databases.")
            reference_kmer_dict = merge_dicts(reference_kmer_dict_list)
    # ------------------------------------------------------

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
            with open(query_input, 'r') as query_input_fh:
                for line in query_input_fh:
                    query_list.append(line.strip())
            for index, query in enumerate(query_list):
                query_name = str(Path(query).name)
                if extension != None:
                    query_name = query_name.replace(extension, "")
                if query_hmms != None:
                    query_protein_list = []
                    with open(query_proteins, 'r') as query_protein_fh:
                        for line in query_protein_fh:
                            query_protein_list.append(line.strip())
                    query_file_names[query_name] = [None, query_protein_list[index], query, query + '.filt']
                elif query_proteins != None:
                    query_file_names[query_name] = [None, query, query + '.hmm', query + '.hmm.filt']
                elif query_genomes != None:
                    query_file_names[query_name] = [query, query + '.faa', query + '.faa.hmm', query + '.faa.hmm.filt']
    # For viral genomes
    else:
        if query_database != None:
            pass
        else:
            with open(query_input, 'r') as query_input_fh:
                for line in query_input_fh:
                    query_list.append(line.strip())
            for index, query in enumerate(query_list):
                query_name = str(Path(query).name)
                if extension != None:
                    query_name = query_name.replace(extension, "")
                if query_proteins != None:
                    query_file_names[query_name] = [None, query]
                elif query_genomes != None:
                    query_file_names[query_name] = [query, query + '.faa']

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
                with open(reference_input, 'r') as reference_input_fh:
                    for line in reference_input_fh:
                        reference_list.append(line.strip())
                for index, reference in enumerate(reference_list):
                    reference_name = str(Path(reference).name)
                    if extension != None:
                        reference_name = reference_name.replace(extension, "")
                    if reference_hmms != None:
                        reference_protein_list = []
                        with open(reference_proteins, 'r') as reference_protein_fh:
                            for line in reference_protein_fh:
                                reference_protein_list.append(line.strip())
                        reference_file_names[reference_name] = [None, reference_protein_list[index], reference, reference + '.filt']
                    elif reference_proteins != None:
                        reference_file_names[reference_name] = [None, reference, reference + '.hmm', reference + '.hmm.filt']
                    elif query_genomes != None:
                        reference_file_names[reference_name] = [reference, reference + '.faa', reference + '.faa.hmm', reference + '.faa.hmm.filt']
        # For viral genomes
        else:
            if reference_database != None:
                pass
            else:
                with open(reference_input, 'r') as reference_input_fh:
                    for line in reference_input_fh:
                        reference_list.append(line.strip())
                for index, reference in enumerate(reference_list):
                    reference_name = str(Path(reference).name)
                    if extension != None:
                        reference_name = reference_name.replace(extension, "")
                    if reference_proteins != None:
                        reference_file_names[reference_name] = [None, reference]
                    elif query_genomes != None:
                        reference_file_names[reference_name] = [reference, reference + '.faa']
    # ------------------------------------------------------

    # Pre-index and store databases
    # ------------------------------------------------------
    # Pre-index queries
    if query_kmer_dict == None:
        print("Processing queries...")
        # If using bacterial genomes
        if virus == False:
            if query_hmms != None:
                query_hmm_results = query_list
            elif query_proteins != None:
                query_protein_files = query_list
                print("Searching against HMM models...")
                try:
                    pool = multiprocessing.Pool(threads)
                    query_hmm_results = pool.map(run_hmmsearch, query_protein_files)
                finally:
                    pool.close()
                    pool.join()
            elif query_genomes != None:
                print("Predicting proteins...")
                # Predict query proteins 
                try:
                    pool = multiprocessing.Pool(threads)
                    query_protein_files = pool.map(run_prodigal, query_list)
                finally:
                    pool.close()
                    pool.join()
                print("Done!")
                print("Searching against HMM models...")
                # Run hmmsearch against proteins predicted
                try:
                    pool = multiprocessing.Pool(threads)
                    query_hmm_results = pool.map(run_hmmsearch, query_protein_files)
                finally:
                    pool.close()
                    pool.join()
                print("Done!")
            print("Filtering query hmmsearch results...")
            # Filter query HMM search results
            try:
                pool = multiprocessing.Pool(threads)
                pool.map(partial(hmm_filter, keep=keep), query_hmm_results)
            finally:
                pool.close()
                pool.join()
            print("Extracting kmers from query proteins...")
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
            query_kmer_dict = merge_dicts(kmer_results)
            del kmer_results
        # If using viral genomes
        else:
            if query_genomes != None:
                print("Predicting proteins...")
                # Predict query proteins 
                try:
                    pool = multiprocessing.Pool(threads)
                    query_protein_files = pool.map(run_prodigal_virus, query_list)
                finally:
                    pool.close()
                    pool.join()
                print("Done!")
            elif query_proteins != None:
                query_protein_files = query_list
            print("Extracting kmers from query proteins...")
            query_information = []
            for name, values in query_file_names.items():
                query_information.append((name, values[1], 4))
            try:
                pool = multiprocessing.Pool(threads)
                kmer_results = pool.map(read_viral_kmers_from_file, query_information)
            finally:
                pool.close()
                pool.join()
            query_kmer_dict = merge_dicts(kmer_results)
            del kmer_results

    # Pre-index references (if different from queries)
    if same_inputs == False and reference_kmer_dict == None:
        print("Processing references...")
        # If using bacterial genomes
        if virus == False:
            if reference_hmms != None:
                reference_hmm_results = reference_list
            elif reference_proteins != None:
                reference_protein_files = reference_list
                print("Searching against HMM models...   ")
                try:
                    pool = multiprocessing.Pool(threads)
                    reference_hmm_results = pool.map(run_hmmsearch, reference_protein_files)
                finally:
                    pool.close()
                    pool.join()
            if reference_genomes != None:
                print("Predicting proteins...")
                # Predict reference proteins 
                try:
                    pool = multiprocessing.Pool(threads)
                    reference_protein_files = pool.map(run_prodigal, reference_list)
                finally:
                    pool.close()
                    pool.join()
                print("Done!")
                print("Searching against HMM models...")
                # Run hmmsearch against proteins predicted
                try:
                    pool = multiprocessing.Pool(threads)
                    reference_hmm_results = pool.map(run_hmmsearch, reference_protein_files)
                finally:
                    pool.close()
                    pool.join()
                print("Done!")
            print("Filtering reference hmmsearch results...")
            # Filter reference HMM search results
            try:
                pool = multiprocessing.Pool(threads)
                pool.map(partial(hmm_filter, keep=keep), reference_hmm_results)
            finally:
                pool.close()
                pool.join()
            print("Extracting kmers from reference proteins...") 
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
            reference_kmer_dict = merge_dicts(kmer_results)
            del kmer_results
        # If using viral genomes
        else:
            if query_genomes != None:
                print("Predicting proteins...")
                # Predict query proteins 
                try:
                    pool = multiprocessing.Pool(threads)
                    query_protein_files = pool.map(run_prodigal, query_list)
                finally:
                    pool.close()
                    pool.join()
                print("Done!")
            elif query_proteins != None:
                query_protein_files = query_list
            print("Extracting kmers from query proteins...")
            reference_information = []
            for name, values in reference_file_names.items():
                reference_information.append((name, values[1], 4))
            try:
                pool = multiprocessing.Pool(threads)
                kmer_results = pool.map(read_viral_kmers_from_file, reference_information)
            finally:
                pool.close()
                pool.join()
            query_kmer_dict = merge_dicts(kmer_results)
            del kmer_results
    # ------------------------------------------------------
    
    # Create or database(s) and compress it(them)
    # ------------------------------------------------------
    if same_inputs == True and query_database == None:
        print("Saving pre-indexed database...")
        query_database_name = query_input + '.db.gz'
        with gzip.open(query_database_name, 'wb') as database_handle:
            pickle.dump(query_kmer_dict, database_handle, protocol=4)
    if same_inputs == False and query_database == None and reference_database == None:
        print("Saving pre-indexed databases...")
        query_database_name = query_input + '.db.gz'
        reference_database_name = reference_input + '.db.gz'
        with gzip.open(query_database_name, 'wb') as database_handle:
            pickle.dump(query_kmer_dict, database_handle, protocol=4)
        with gzip.open(reference_database_name, 'wb') as database_handle:
            pickle.dump(reference_kmer_dict, database_handle, protocol=4)
    elif same_inputs == False and query_database == None:
        print("Saving pre-indexed query database...")
        query_database_name = query_input + '.db.gz'
        with gzip.open(query_database_name, 'wb') as database_handle:
            pickle.dump(query_kmer_dict, database_handle, protocol=4)
    elif same_inputs == False and reference_database == None:
        print("Saving pre-indexed reference database...")
        reference_database_name = reference_input + '.db.gz'
        with gzip.open(reference_database_name, 'wb') as database_handle:
            pickle.dump(reference_kmer_dict, database_handle, protocol=4)
    # ------------------------------------------------------
    # Calculate Jaccard distances
    # ------------------------------------------------------
    if index_db == True:
        print("Finished pre-indexing databases.")
        print("Next time you can run the program using only these files with --qd and(or) --rd.")
    else:
        print("Calculating shared kmer fraction...")
        if virus == False:
            if same_inputs == True:
                # Create global kmer index dictionary "global_kmer_index_dictionary"
                global_unique_kmers([query_kmer_dict])
                query_kmer_dict, query_smart_args_tempdir = transform_kmer_dicts_to_arrays(query_kmer_dict, temporal_working_directory, single_dataset=True)
                print("Beginning FastAAI pairwise calculations now.")
                try:
                    pool = multiprocessing.Pool(threads, initializer = single_dictionary_initializer, initargs = (query_kmer_dict,))
                    Fraction_Results = pool.map(single_kaai_parser, query_smart_args_tempdir)
                finally:
                    pool.close()
                    pool.join()
            else:
                global_unique_kmers([query_kmer_dict, reference_kmer_dict])
                query_kmer_dict, query_smart_args_tempdir = transform_kmer_dicts_to_arrays(query_kmer_dict, temporal_working_directory, single_dataset=False)
                reference_kmer_dict, _ref_smart_args_tempdir = transform_kmer_dicts_to_arrays(reference_kmer_dict, temporal_working_directory, single_dataset=False)
                print("Beginning FastAAI pairwise calculations now.")
                try:
                    pool = multiprocessing.Pool(threads, initializer = two_dictionary_initializer, initargs = (query_kmer_dict, reference_kmer_dict))
                    Fraction_Results = pool.map(double_kaai_parser, query_smart_args_tempdir)
                finally:
                    pool.close()
                    pool.join()
        else:
            if same_inputs == True:
                query_id_list = query_kmer_dict.keys()
                try:
                    pool = multiprocessing.Pool(threads, initializer = single_dictionary_initializer, initargs = (query_kmer_dict,))
                    Fraction_Results = pool.map(single_virus_kaai_parser, query_id_list)
                finally:
                    pool.close()
                    pool.join()
            else:
                query_id_list = query_kmer_dict.keys()
                try:
                    pool = multiprocessing.Pool(threads, initializer = two_dictionary_initializer, initargs = (query_kmer_dict, reference_kmer_dict))
                    Fraction_Results = pool.map(double_viral_kaai_parser, query_id_list)
                finally:
                    pool.close()
                    pool.join()
    # ------------------------------------------------------
    
    # Merge results into a single output
    # ------------------------------------------------------
        print("Merging results...")
        with open(output, 'w') as outfile:
            for file in Fraction_Results:
                with open(file) as Temp:
                    shutil.copyfileobj(Temp, outfile)
                file.unlink()
        print("kAAI finishied correctly on {}".format(datetime.datetime.now()))
    # ------------------------------------------------------
    # If comparing viral genomes


	
	
	
if __name__ == "__main__":
    main()
