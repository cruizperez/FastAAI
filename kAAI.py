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
import subprocess
import argparse
import multiprocessing
import datetime
import shutil
from numpy import random
from pathlib import Path
from sys import argv
from sys import exit
from functools import partial
import sys


################################################################################
"""---1.0 Define Functions---"""
# --- Run prodigal ---
def run_prodigal(input_file):
    """
    Runs prodigal, compares translation tables and stores faa files

    Arguments:
       input_file -- Path to genome FastA file
    
    Returns:
        output -- Path to amino acid fasta result
    """
    file_path = Path(input_file)
    prefix = Path(file_path.stem)
    folder = file_path.parent
    final_output = folder / prefix.with_suffix('.faa')
    output_11 = folder / prefix.with_suffix('.faa.11')
    temp_output = folder / prefix.with_suffix('.temp')
    subprocess.call(["prodigal", "-i", str(file_path), "-a", str(output_11), 
                    "-p", "meta", "-q", "-o", str(temp_output)])
    output_4 = folder / prefix.with_suffix('.faa.4')
    temp_output = folder / prefix.with_suffix('.temp')
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
        shutil.copy(output_4, final_output)
    else:
        shutil.copy(str(output_11), str(final_output))
    
    # Remove intermediate files
    output_4.unlink()
    output_11.unlink()
    temp_output.unlink()

    with open(final_output) as final_protein, open(temp_output) as temporal_file:
        for line in final_protein:
            if line.startswith(">"):
                temporal_file.write("{}".format(line))
            else:
                line.replace('*', '')
                temporal_file.write("{}".format(line))
    shutil.copy(str(temp_output), str(final_output))
    temp_output.unlink()

    return final_output

# --- Run hmmsearch ---
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
    name = Path(file_path.name)
    output = folder / name.with_suffix('.hmm')
    temp_output = folder / name.with_suffix('.temp')
    script_path = Path(__file__)
    script_dir = script_path.parent
    HMM_complete_model = script_dir / "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"
    subprocess.call(["hmmsearch", "--tblout", str(output), "-o", str(temp_output), "--cut_tc", "--cpu", "1",
                    str(HMM_complete_model), str(file_path)])
    temp_output.unlink()
    return output

# --- Filter HMM results for best matches ---
def HMM_filter(SCG_HMM_file, keep):
    """
    Filters HMM results for best hits per protein
    
    Arguments:
        SCG_HMM_file {file path} -- Path to HMM results file
        keep {bool} -- Keep HMM files
    
    Returns:
        outfile -- Path to filtered files
    """
    HMM_path = Path(SCG_HMM_file)
    name = HMM_path.name
    folder = HMM_path.parent
    outfile = folder / Path(name).with_suffix('.filt')
    HMM_hit_dict = {}
    with open(SCG_HMM_file, 'r') as hit_file:
        for line in hit_file:
            if line.startswith("#"):
                continue
            else:
                hit = line.strip().split()
                protein_name = hit[0]
                score = hit[8]
                if protein_name in HMM_hit_dict:
                    #! Attention
                    if score > HMM_hit_dict[protein_name][0] and score >= 50:
                    # if score > HMM_hit_dict[protein_name][0]:
                        HMM_hit_dict[protein_name] = [score, line]
                    elif score < HMM_hit_dict[protein_name][0]:
                        continue
                    else:
                        if random.randint(2) > 0 and score >= 50:
                        # if random.randint(2) > 0:
                            HMM_hit_dict[protein_name] = [score, line]
                else:
                    if score >= 50:
                        HMM_hit_dict[protein_name] = [score, line]
                    # HMM_hit_dict[protein_name] = [score, line]
    with open(outfile, 'w') as output:
        for hits in HMM_hit_dict.values():
            output.write("{}".format(hits[1]))
    if keep == False:
        HMM_path.unlink()
    return outfile



# --- Find Kmers from HMM results ---
def Kmer_Parser(SCG_HMM_file):
    """
    Extract kmers from protein files that have hits
    in the HMM searches.
    
    Arguments:
        SCG_HMM_file {file path} -- Path to filtered HMM results.
    
    Returns:
        [genome_kmers] -- Dictionary of kmers per gene. 
    """
    HMM_path = Path(SCG_HMM_file)
    name = HMM_path.name
    genome = HMM_path.stem
    folder = HMM_path.parent
    protein_file = folder / Path(name).with_suffix('.faa')
    positive_matches = {}
    positive_proteins = []
    with open(HMM_path, 'r') as HMM_Input:
        for line in HMM_Input:
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
    genome_kmers = {name : scg_kmers}
    #! Attention
    print("Dictionary size: {}".format(get_size(genome_kmers)))
    SCG_HMM_file.unlink()
    return genome_kmers

# --- Read Kmers from SCGs ---
def read_kmers_from_file(filename, positive_hits, ksize):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    scg_kmers = {}
    with open(filename) as Fasta_in:
        for title, sequence in SimpleFastaParser(Fasta_in):
            protein = title.split()[0]
            if protein in positive_hits:
                kmers = build_kmers(sequence, ksize)
                scg_kmers[protein] = kmers
    return scg_kmers

# --- Build Kmers ---
def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
    kmers_set = set(kmers)
    #!Attention
    # kmers_set = ','.join(set(kmers))
    return kmers_set

# --- Parse kAAI ---
def kAAI_Parser(query_id):
    """
    Calculates Jaccard distances on kmers from proteins shared
    
    Arguments:
        query_id {str} -- Id of the query genome
    
    Returns:
        [Path to output] -- Path to output file
    """
    file_path = Path(query_id)
    running_folder = Path.cwd()
    temp_output = running_folder / file_path.with_suffix('.aai.temp')
    query_num_scg = len(total_kmer_dictionary[query_id])
    query_scg_list = total_kmer_dictionary[query_id].keys()
    with open(temp_output, 'w') as out_file:
        for target_genome, scg_ids in total_kmer_dictionary.items():
            start = datetime.datetime.now().time()
            jaccard_similarities = []
            target_num_scg = len(scg_ids)
            target_scg_list = scg_ids.keys()
            if query_num_scg > target_num_scg:
                final_scg_list = target_scg_list
            else:
                final_scg_list = query_scg_list
            for accession in final_scg_list:
                if accession in query_scg_list and accession in target_scg_list:
                    #!Attention
                    kmers_query = total_kmer_dictionary[query_id][accession]
                    kmers_target = total_kmer_dictionary[target_genome][accession]
                    # kmers_query = set(total_kmer_dictionary[query_id][accession].split(','))
                    # kmers_target = total_kmer_dictionary[target_genome][accession].split(',')
                    intersection = len(kmers_query.intersection(kmers_target))
                    union = len(kmers_query.union(kmers_target))
                    jaccard_similarities.append(intersection / union)
                else:
                    continue
            try:
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           sum(jaccard_similarities)/len(jaccard_similarities),
                           len(jaccard_similarities), len(final_scg_list)))
            except:
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           "NA", "NA", "NA"))
            end = datetime.datetime.now().time()
            print("Comparison time: {}".format(datetime.datetime.combine(datetime.date.min, end) - datetime.datetime.combine(datetime.date.min, start)))
    return temp_output

# --- Initializer function ---
def child_initialize(_dictionary):
    """
    Make dictionary available for multiprocessing
    """
    global total_kmer_dictionary
    total_kmer_dictionary = _dictionary

def merge_dicts(Dictionaries):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in Dictionaries:
        result.update(dictionary)
    return result

def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size
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
    parser.add_argument('-g', '--genomes', dest='genome_list', action='store', nargs='+', required=False, help='List of input genomes.')
    parser.add_argument('-p', '--proteins', dest='protein_files', action='store', nargs='+', required=False, help='List of input protein files. They should have the .faa extension.')
    parser.add_argument('-s', '--scg_hmm', dest='HMM_files', action='store', nargs='+', required=False, help='List of hmm search results')
    parser.add_argument('-o', '--output', dest='outfile', action='store', required=True, help='Output File')
    parser.add_argument('-t', '--threads', dest='threads', action='store', default=1, type=int, required=False, help='Number of threads to use, by default 1')
    parser.add_argument('-k', '--keep', dest='keep', action='store_false', required=False, help='Keep intermediate files, by default true')
    args = parser.parse_args()

    genome_list = args.genome_list
    protein_files = args.protein_files
    HMM_files = args.HMM_files
    outfile = args.outfile
    threads = args.threads
    keep = args.keep

    print("kAAI started on {}".format(datetime.datetime.now())) # Remove after testing
    # Check input
    # ------------------------------------------------------
    if HMM_files != None and protein_files != None:
        exit('Please provide only one input. You provided Proteins and HMM results')
    elif HMM_files != None and genome_list != None:
        exit('Please provide only one input. You provided HMM results and Genomes')
    elif protein_files != None and genome_list != None:
        exit('Please provide only one input. You provided Proteins and Genomes')
    elif protein_files == None and genome_list == None and HMM_files == None:
        exit('No input provided, please provide genomes "-g", protein "-p", or scg hmm searches "-s"')
    # ------------------------------------------------------

    # Predict proteins and perform HMM searches
    # ------------------------------------------------------
    print("kAAI started on {}".format(datetime.datetime.now()))
    if genome_list != None:
        print("Starting from Genomes.")
        print("Predicting proteins...   ", end="")
        protein_files = []
        try:
            pool = multiprocessing.Pool(threads)
            protein_files = pool.map(run_prodigal, genome_list)
        finally:
            pool.close()
            pool.join()
        print("Done")
        print("Searching HMM models...   ", end="")
        try:
            pool = multiprocessing.Pool(threads)
            HMM_Search_Files = pool.map(run_hmmsearch, protein_files)
        finally:
            pool.close()
            pool.join()
        print("Done")
    elif protein_files != None:
        print("Starting from Proteins.")
        print("Searching HMM models...   ", end="")
        try:
            pool = multiprocessing.Pool(threads)
            HMM_Search_Files = pool.map(run_hmmsearch, protein_files)
        finally:
            pool.close()
            pool.join()
        print("Done")
    elif HMM_files != None:
        print("Starting from HMM searches.")
        HMM_Search_Files = HMM_files
    # ------------------------------------------------------
    
    # Filter HMM results, retaining best hit per protein
    # ------------------------------------------------------
    print("Filtering HMM results...")
    print(datetime.datetime.now())
    try:
        pool = multiprocessing.Pool(threads)
        filtered_files = pool.map(partial(HMM_filter, keep=keep), HMM_Search_Files)
    finally:
        pool.close()
        pool.join()
    # ------------------------------------------------------

    # Find kmers per SCG per genome
    # ------------------------------------------------------
    print("Parsing HMM results...")
    print(datetime.datetime.now())
    try:
        pool = multiprocessing.Pool(threads)
        Kmer_Results = pool.map(Kmer_Parser, filtered_files)
    finally:
        pool.close()
        pool.join()
    Final_Kmer_Dict = merge_dicts(Kmer_Results)
    del Kmer_Results
    # ------------------------------------------------------
    
    # Calculate Jaccard distances
    # ------------------------------------------------------
    print("Calculating shared Kmer fraction...")
    print(datetime.datetime.now())
    ID_List = Final_Kmer_Dict.keys()
    try:
        pool = multiprocessing.Pool(threads, initializer = child_initialize, initargs = (Final_Kmer_Dict,))
        Fraction_Results = pool.map(kAAI_Parser, ID_List)
    finally:
        pool.close()
        pool.join()
    # ------------------------------------------------------
    
    # Merge results into a single output
    # ------------------------------------------------------
    print("Merging results...")
    print(datetime.datetime.now())
    with open(outfile, 'w') as output:
        for file in Fraction_Results:
            with open(file) as Temp:
                shutil.copyfileobj(Temp, output)
            file.unlink()
    print("kAAI finishied correctly!!")
    # ------------------------------------------------------

if __name__ == "__main__":
    main()
