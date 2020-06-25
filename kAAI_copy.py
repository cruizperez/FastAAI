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
import sys, textwrap



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
    output_11 = folder / (filename + '.faa11')
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
                        if random.randint(2) > 0:
                            hmm_hit_dict[protein_name] = [score, line]
                else:
                    hmm_hit_dict[protein_name] = [score, line]
    with open(outfile, 'w') as output:
        for hits in hmm_hit_dict.values():
            output.write("{}".format(hits[1]))
    if keep == False:
        hmm_path.unlink()
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
    # hmm_path = Path(scg_hmm_file)
    # name = hmm_path.name
    # folder = hmm_path.parent
    # protein_file = folder / Path(name).with_suffix('.faa')
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
    Path(scg_hmm_file).unlink()
    return genome_kmers
# ------------------------------------------------------

# --- Extract kmers from protein sequences ---
# ------------------------------------------------------
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

# --- Parse kAAI when query == reference ---
# ------------------------------------------------------
def single_kaai_parser(query_id):
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
    # Get number and list of SCG detected in query
    query_num_scg = len(query_kmer_dictionary[query_id])
    query_scg_list = query_kmer_dictionary[query_id].keys()
    # Start comparison with all genomes in the query dictionary
    with open(temp_output, 'w') as out_file:
        for target_genome, scg_ids in query_kmer_dictionary.items():
            start = datetime.datetime.now().time()
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
                    kmers_target = query_kmer_dictionary[target_genome][accession].split(',')
                    # Calculate jaccard_similarity
                    intersection = len(kmers_query.intersection(kmers_target))
                    union = len(kmers_query.union(kmers_target))
                    jaccard_similarities.append(intersection / union)
                else:
                    continue
            try:
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           round(sum(jaccard_similarities)/len(jaccard_similarities), 4),
                           len(jaccard_similarities), len(final_scg_list)))
            except:
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           "NA", "NA", "NA"))
            end = datetime.datetime.now().time()
            print("Comparison time: {}".format(datetime.datetime.combine(datetime.date.min, end) - datetime.datetime.combine(datetime.date.min, start)))
    return temp_output
# ------------------------------------------------------

# --- Parse kAAI when query == reference ---
# ------------------------------------------------------
def double_kaai_parser(query_id):
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
    # Get number and list of SCG detected in query
    query_num_scg = len(query_kmer_dictionary[query_id])
    query_scg_list = query_kmer_dictionary[query_id].keys()
    # Start comparison with all genomes in the query dictionary
    with open(temp_output, 'w') as out_file:
        for target_genome, scg_ids in ref_kmer_dictionary.items():
            start = datetime.datetime.now().time()
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
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           round(sum(jaccard_similarities)/len(jaccard_similarities), 4),
                           len(jaccard_similarities), len(final_scg_list)))
            except:
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
                           "NA", "NA", "NA"))
            end = datetime.datetime.now().time()
            print("Comparison time: {}".format(datetime.datetime.combine(datetime.date.min, end) - datetime.datetime.combine(datetime.date.min, start)))
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
def two_dictionary_initializer(_query_dictionary, _ref_dictionary):
    """
    Make dictionary available for multiprocessing
    """
    global query_kmer_dictionary
    global ref_kmer_dictionary
    query_kmer_dictionary = _query_dictionary
    ref_kmer_dictionary = _ref_dictionary
# ------------------------------------------------------

# --- Merge kmer dictionaries ---
# ------------------------------------------------------
def merge_dicts(dictionaries):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dictionaries:
        result.update(dictionary)
    return result
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
    mandatory_options = parser.add_argument_group('Mandatory i/o options.')
    mandatory_options.add_argument('-q', '--queries', dest='queries', action='store', required=True,
                                    help='File with list of queries.')
    mandatory_options.add_argument('-r', '--references', dest='references', action='store', required=True,
                                    help='File with list of references.')
    mandatory_options.add_argument('-g', '--genomes', dest='genomes', action='store_true', required=False,
                                    help='Start from genome files.')
    mandatory_options.add_argument('-p', '--proteins', dest='proteins', action='store_true', required=False,
                                    help='Start from protein files.')
    mandatory_options.add_argument('-s', '--hmms', dest='hmms', action='store_true', required=False, 
                                    help=textwrap.dedent('''
                                    Start from HMM search result files.
                                    If you select this option you must also provide two files with 
                                    a list of protein files for the query and the references.
                                    (in the same order).
                                    '''))
    mandatory_options.add_argument('-o', '--output', dest='outfile', action='store', required=True, help='Output File')
    additional_input_options = parser.add_argument_group('Additional input options.')
    additional_input_options.add_argument('-e', '--ext', dest='extension', action='store', required=False, 
                                            help='Extension to remove from original filename, e.g. ".fasta"')
    additional_input_options.add_argument('--prot_query', dest='proteins_query', action='store', required=False,
                                            help=textwrap.dedent('''
                                            File with the list of proteins for the queries (same order as queries).
                                            Only required if starting from hmmsearch results (-s).
                                            '''))
    additional_input_options.add_argument('--prot_ref', dest='proteins_reference', action='store', required=False,
                                            help=textwrap.dedent('''
                                            File with the list of proteins for the references (same order as references).
                                            Only required if starting from hmmsearch results (-s).
                                            '''))
    misc_options = parser.add_argument_group('Miscellaneous options')
    misc_options.add_argument('-t', '--threads', dest='threads', action='store', default=1, type=int, required=False, help='Number of threads to use, by default 1')
    misc_options.add_argument('-k', '--keep', dest='keep', action='store_false', required=False, help='Keep intermediate files, by default true')

    args = parser.parse_args()

    queries = args.queries
    references = args.references
    genomes = args.genomes
    proteins = args.proteins
    hmms = args.hmms
    outfile = args.outfile
    extension = args.extension
    proteins_query = args.proteins_query
    proteins_reference = args.proteins_reference
    threads = args.threads
    keep = args.keep

    print("kAAI started on {}".format(datetime.datetime.now())) # Remove after testing
    # Check user input
    # ------------------------------------------------------
    if queries == None or references == None:
        exit('Please prove a query and reference list files.')
    if (genomes + proteins + hmms) > 1 or (genomes + proteins + hmms) == 0:
        exit('Please specify if the program should start from genomes (-g), proteins (-p) OR hmm search results (-s)')
    if hmms == True:
        if proteins_query == None or proteins_reference == None:
            print("You chose to start from hmmsearch results (-s).")
            print("However, I also need the location of the proteins used for the search.")
            exit("Please provide it with --prot_query for the queries and with --prot_ref for the references.")

    # ------------------------------------------------------

    # Get files from the query and reference lists
    # ------------------------------------------------------
    same_genomes = False
    if queries == references:
        same_genomes = True
        print('You specified the same query and reference files.')
        print('I will perform an all vs all comparison :)')
    query_list = []
    reference_list = []
    with open(queries, 'r') as query_file:
        for line in query_file:
            query_list.append(line.strip())
    if same_genomes == False:
        with open(references, 'r') as reference_file:
            for line in reference_file:
                reference_list.append(line.strip())
    # ------------------------------------------------------

    # Create a dictionary with resulting filenames and a list with dictionary keys
    # ------------------------------------------------------
    query_file_names = {}
    reference_file_names = {}
    query_key_names = []
    reference_key_names = []
    for query in query_list:
        query_name = str(Path(query).name)
        # Add final name to query list
        if extension != None:
            query_name = query_name.replace(extension, "")
            query_key_names.append(query_name)
        else:
            query_key_names.append(query_name)
        if genomes == True:
            query_file_names[query_name] = [query, query + '.faa', query + '.faa.hmm', query + '.faa.hmm.filt']
        elif proteins == True:
            query_file_names[query_name] = [None, query, query + '.hmm', query + '.hmm.filt']
        elif hmms == True:
            query_file_names[query_name] = [None, None, query, query + '.filt']
    if same_genomes == False:
        for reference in reference_list:
            reference_name = str(Path(reference).name)
            # Add final name to query list
            if extension != None:
                reference_name = reference_name.replace(extension, "")
                reference_key_names.append(reference_name)
            else:
                reference_key_names.append(reference_name)
            if genomes == True:
                reference_file_names[reference_name] = [reference, reference + '.faa', reference + '.faa.hmm', reference + '.faa.hmm.filt']
            elif proteins == True:
                reference_file_names[reference_name] = [None, reference, reference + '.hmm', reference + '.hmm.filt']
            elif hmms == True:
                reference_file_names[reference_name] = [None, None, reference, reference + '.filt']
    print(query_file_names) #! Remove after testing
    # ------------------------------------------------------
    #TODO: Add support for databases here (future)

    # Predict proteins and perform HMM searches
    # ------------------------------------------------------
    if genomes == True:
        print("Starting from genome files.")
        print("Predicting proteins...   ", end="")
        # Predict query proteins 
        try:
            pool = multiprocessing.Pool(threads)
            query_protein_files = pool.map(run_prodigal, query_list)
        finally:
            pool.close()
            pool.join()
        # Predict reference proteins
        if same_genomes == False:
            try:
                pool = multiprocessing.Pool(threads)
                reference_protein_files = pool.map(run_prodigal, reference_list)
            finally:
                pool.close()
                pool.join()
        print(query_file_names)#! Remove after testing
        print("Done")
        print("Searching against HMM models...   ", end="")
        # Search queries agains HMM SCG models
        try:
            pool = multiprocessing.Pool(threads)
            query_hmm_results = pool.map(run_hmmsearch, query_protein_files)
        finally:
            pool.close()
            pool.join()
        # Search references agains HMM SCG models
        if same_genomes == False:
            try:
                pool = multiprocessing.Pool(threads)
                reference_hmm_results = pool.map(run_hmmsearch, reference_protein_files)
            finally:
                pool.close()
                pool.join()
        print(query_file_names)#! Remove after testing
        print("Done")
    elif proteins == True:
        query_protein_files = query_list
        reference_protein_files = reference_list
        print("Starting from protein files.")
        print("Searching against HMM models...   ", end="")
        # Search queries against HMM models
        try:
            pool = multiprocessing.Pool(threads)
            query_hmm_results = pool.map(run_hmmsearch, query_protein_files)
        finally:
            pool.close()
            pool.join()
        # Search references against HMM models
        if same_genomes == False:
            try:
                pool = multiprocessing.Pool(threads)
                reference_hmm_results = pool.map(run_hmmsearch, reference_protein_files)
            finally:
                pool.close()
                pool.join()
        print(query_file_names)#! Remove after testing
        print("Done")
    elif hmms == True:
        print("Starting from HMM searches.")
        query_hmm_results = query_list
        reference_hmm_results = reference_list
    # ------------------------------------------------------
    
    # Filter HMM results, retaining best hit per protein
    # ------------------------------------------------------
    print("Filtering HMM results...", end="")
    print(datetime.datetime.now()) #! Remove after testing
    # Filter query HMM search results
    try:
        pool = multiprocessing.Pool(threads)
        query_filtered_hmms = pool.map(partial(hmm_filter, keep=keep), query_hmm_results)
    finally:
        pool.close()
        pool.join()
    # Filter reference HMM search results
    if same_genomes == False:
        try:
            pool = multiprocessing.Pool(threads)
            query_filtered_hmms = pool.map(partial(hmm_filter, keep=keep), reference_hmm_results)
        finally:
            pool.close()
            pool.join()
    print('Done')

    # ------------------------------------------------------

    # Find kmers per SCG per genome
    # ------------------------------------------------------
    print("Parsing HMM results...")
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
    # Finding kmers for all references
    if same_genomes == False:
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
    # ------------------------------------------------------
    
    # Calculate Jaccard distances
    # ------------------------------------------------------
    print("Calculating shared Kmer fraction...")
    print(datetime.datetime.now())
    if same_genomes == True:
        query_id_list = query_kmer_dict.keys()
        try:
            pool = multiprocessing.Pool(threads, initializer = single_dictionary_initializer, initargs = (query_kmer_dict,))
            Fraction_Results = pool.map(single_kaai_parser, query_id_list)
        finally:
            pool.close()
            pool.join()
    else:
        query_id_list = query_kmer_dict.keys()
        try:
            pool = multiprocessing.Pool(threads, initializer = two_dictionary_initializer, initargs = (query_kmer_dict, reference_kmer_dict))
            Fraction_Results = pool.map(double_kaai_parser, query_id_list)
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
