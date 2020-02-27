#!/usr/bin/env python

"""
########################################################################
# Author:	   Carlos Ruiz
# Intitution:   Georgia Institute of Technology
# Version:	  1.0
# Date:		 Nov 27 2019

# Description: Calculates the average amino acid identity using k-mers
from single copy genes. It is a faster version of the regular AAI (Blast
or Diamond) and the hAAI implemented in MiGA.
########################################################################
"""

################################################################################
"""---0.0 Import Modules---"""
import subprocess
import numpy as np
import argparse
import multiprocessing
import datetime
import shutil
from pathlib import Path
from sys import argv
from sys import exit
from functools import partial


################################################################################
"""---1.0 Define Functions---"""
# --- Run prodigal ---
def run_prodigal(input_file):
    """
    Runs prodigal and stores faa files

    Arguments:
       input_file -- Path to genome FastA file
    
    Returns:
        output -- Path to amino acid fasta result
    """
    file_path = Path(input_file)
    prefix = Path(file_path.stem)
    folder = file_path.parent
    output = folder / prefix.with_suffix('.faa')
    temp_output = folder / prefix.with_suffix('.temp')
    subprocess.call(["prodigal", "-i", str(file_path), "-a", str(output), "-p", "meta", "-q", "-o", str(temp_output)])
    temp_output.unlink()
    return output

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
    subprocess.call(["hmmsearch", "--tblout", str(output), "-o", str(temp_output), "--cut_ga", "--cpu", "1",
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
    genome = HMM_path.stem
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
                    if score > HMM_hit_dict[protein_name][0]:
                        HMM_hit_dict[protein_name] = [score, line]
                    elif score < HMM_hit_dict[protein_name][0]:
                        continue
                    else:
                        if np.random.randint(2) > 0:
                            HMM_hit_dict[protein_name] = [score, line]
                else:
                    HMM_hit_dict[protein_name] = [score, line]
    with open(outfile, 'w') as output:
        for hits in HMM_hit_dict.values():
            output.write("{}".format(hits[1]))
    if keep == False:
        HMM_Path.unlink()
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
    kmer_dic = {}
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
    genome_kmers = {genome : scg_kmers}
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
    return kmers_set

# --- Read Kmers from files ---
def read_total_kmers_from_file(filename, positive_hits, ksize):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    all_kmers = []
    with open(filename) as Fasta_in:
        for _, sequence in SimpleFastaParser(Fasta_in):
            kmers = build_kmers(sequence, ksize)
            all_kmers += kmers
    return all_kmers

# --- Parse kAAI ---
def kAAI_Parser(query_id):

    file_path = Path(query_id)
    running_folder = Path.cwd()
    temp_output = running_folder / file_path.with_suffix('.aai.temp')
    query_num_scg = len(total_kmer_dictionary[query_id])
    query_scg_list = total_kmer_dictionary[query_id].keys()
    with open(temp_output, 'w') as out_file:
        for target_genome, scg_ids in total_kmer_dictionary.items():
            jaccard_similarities = []
            target_num_scg = len(scg_ids)
            target_scg_list = scg_ids.keys()
            if query_num_scg > target_num_scg:
                final_scg_list = target_scg_list
            else:
                final_scg_list = query_scg_list
            for accession in final_scg_list:
                if accession in query_scg_list and accession in target_scg_list:
                    kmers_query = total_kmer_dictionary[query_id][accession]
                    kmers_target = total_kmer_dictionary[target_genome][accession]
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
    return temp_output

# --- Initialize function ---
def child_initialize(_dictionary):
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
    parser.add_argument('-g', '--genomes', dest='Genome_List', action='store', nargs='+', required=False, help='List of input genomes. Implies step 1')
    parser.add_argument('-p', '--proteins', dest='Protein_Files', action='store', nargs='+', required=False, help='List of input protein files')
    parser.add_argument('-s', '--scg_hmm', dest='HMM_Files', action='store', nargs='+', required=False, help='List of hmm search results')
    parser.add_argument('-o', '--output', dest='Output', action='store', required=True, help='Output File')
    parser.add_argument('-t', '--threads', dest='Threads', action='store', default=1, type=int, required=False, help='Number of threads to use, by default 1')
    parser.add_argument('-k', '--keep', dest='keep', action='store_false', required=False, help='Keep intermediate files, by default true')
    args = parser.parse_args()

    Genome_List = args.Genome_List
    Protein_Files = args.Protein_Files
    HMM_Files = args.HMM_Files
    Output = args.Output
    Threads = args.Threads
    keep = args.keep

    # Predict proteins and perform HMM searches
    print("kAAI started on {}".format(datetime.datetime.now())) # Remove after testing
    if Genome_List != None:
        print("Starting from Genomes...")
        print("Predicting proteins...")
        Protein_Files = []
        try:
            pool = multiprocessing.Pool(Threads)
            Protein_Files = pool.map(run_prodigal, Genome_List)
        finally:
            pool.close()
            pool.join()
        print("Searching HMM models...")
        try:
            pool = multiprocessing.Pool(Threads)
            HMM_Search_Files = pool.map(run_hmmsearch, Protein_Files)
        finally:
            pool.close()
            pool.join()
    elif Protein_Files != None:
        print("Starting from Proteins...")
        print("Searching HMM models...")
        try:
            pool = multiprocessing.Pool(Threads)
            HMM_Search_Files = pool.map(run_hmmsearch, Protein_Files)
        finally:
            pool.close()
            pool.join()
    elif HMM_Files != None:
        print("Starting from HMM searches...")
        HMM_Search_Files = HMM_Files
    elif HMM_Files != None and Protein_Files != None:
        exit('Please provide only one input. You provided Proteins and HMM results')
    elif HMM_Files != None and Genome_List != None:
        exit('Please provide only one input. You provided HMM results and Genomes')
    elif Protein_Files != None and Genome_List != None:
        exit('Please provide only one input. You provided Proteins and Genomes')
    else:
        exit('No input provided, please provide genomes "-g", protein "-p", or scg hmm searches "-s"')
    # ---------------------------------------------------------------

    # Filter HMM results, retaining best hit per protein
    print("Filtering HMM results...")
    print(datetime.datetime.now())
    try:
        pool = multiprocessing.Pool(Threads)
        filtered_files = pool.map(partial(HMM_filter, keep=keep), HMM_Search_Files)
    finally:
        pool.close()
        pool.join()
    
    # Parse HMM results, calculate distances and compile results
    print("Parsing HMM results...")
    print(datetime.datetime.now())
    try:
        pool = multiprocessing.Pool(Threads)
        Kmer_Results = pool.map(Kmer_Parser, filtered_files)
    finally:
        pool.close()
        pool.join()
    Final_Kmer_Dict = merge_dicts(Kmer_Results)

    # Calculate shared Kmer fraction
    print("Calculating shared Kmer fraction...")
    print(datetime.datetime.now())
    ID_List = Final_Kmer_Dict.keys()
    try:
        pool = multiprocessing.Pool(Threads, initializer = child_initialize, initargs = (Final_Kmer_Dict,))
        Fraction_Results = pool.map(kAAI_Parser, ID_List)
    finally:
        pool.close()
        pool.join()

    # Merge results into a single output
    print("Merging results...")
    print(datetime.datetime.now())
    with open(Output, 'w') as OutFile:
        for file in Fraction_Results:
            with open(file) as Temp:
                shutil.copyfileobj(Temp, OutFile)
            file.unlink()


if __name__ == "__main__":
    main()
