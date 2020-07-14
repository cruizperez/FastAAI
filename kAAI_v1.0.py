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
                                    help='File with pre-indexed query database.')
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
                                    help='File with pre-indexed reference database.')
    mandatory_options.add_argument('-o', '--output', dest='output', action='store', required=True, help='Output file')
    additional_input_options = parser.add_argument_group('Behavior modification options.')
    additional_input_options.add_argument('-e', '--ext', dest='extension', action='store', required=False, 
                                            help='Extension to remove from original filename, e.g. ".fasta"')
    additional_input_options.add_argument('-i', '--index', dest='index_db', action='store_true', required=False, 
                                            help='Only index and store databases, i.e., do not perform comparisons.')
    misc_options = parser.add_argument_group('Miscellaneous options')
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
    extension = args.extension
    index_db = args.index_db
    threads = args.threads
    keep = args.keep

    print("kAAI started on {}".format(datetime.datetime.now()))
    # Check user input
    # ------------------------------------------------------
    # Check if no query was provided
    if query_genomes == None and query_proteins == None and query_hmms == None and query_database == None:
        exit('Please prove a file with a list of queries, e.g., --qg, --qp, --qh, or --qd)')
    # Check query inputs
    query_input = None
    if query_hmms != None:
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
        print("Starting from the pre-indexed query database.")
    # Check if no reference was provided
    if reference_genomes == None and reference_proteins == None and reference_hmms == None and reference_database == None:
        exit('Please prove a file with a list of references, e.g., --rg, --rp, --rh, or --rd)')
    # Check reference inputs
    reference_input = None
    if reference_hmms != None:
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
        print("Starting from the pre-indexed reference database.")
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

    # If using pre-indexed databases, check if they are valid files.
    # ------------------------------------------------------
    # If any of the starting points is from database, then store the
    # kmer structures in the corresponding dictionaries.
    # Otherwise get read the file list and get the filenames
    query_kmer_dict = None
    reference_kmer_dict = None
    # If starting from database and query == reference
    if query_database != None and same_inputs == True:
        if Path(query_database).is_file():
            with gzip.open(query_database, 'rb') as database_handle:
                query_kmer_dict = pickle.load(database_handle)
            if isinstance(query_kmer_dict,dict):
                pass
            else:
                exit("The database appears to have the wrong format. Please provide a correctly formated database.")
    # If starting from database and query != reference
    elif query_database != None and reference_database == None:
        # First check database
        if Path(query_database).is_file():
            with gzip.open(query_database, 'rb') as database_handle:
                query_kmer_dict = pickle.load(database_handle)
            if isinstance(query_kmer_dict,dict):
                pass
            else:
                exit("The query database appears to have the wrong format. Please provide a correctly formated database.")
        else:
            exit("I cannot locate the query database you proveded: {}". format(query_database))
    elif reference_database != None:
        if Path(reference_database).is_file():
            with gzip.open(reference_database, 'rb') as database_handle:
                reference_kmer_dict = pickle.load(database_handle)
            if isinstance(reference_kmer_dict,dict):
                pass
            else:
                exit("The reference database appears to have the wrong format. Please provide a correctly formated database.")
    # ------------------------------------------------------

    # Get files from the query and reference lists and then
    # create a dictionary with resulting filenames and a list with dictionary keys
    # The structure of the dictionary is:
    # original_query, proteins, hmms, filtered_hmms
    # ------------------------------------------------------
    # First parse the query:
    query_list = []
    query_file_names = {}
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
    # Then parse the references:
    reference_list = []
    reference_file_names = {}
    if same_inputs == True:
        pass
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
    # ------------------------------------------------------

    # Pre-index and store databases
    # ------------------------------------------------------
    # Pre-index queries
    if query_kmer_dict == None:
        print("Processing queries...")
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
    # Pre-index references (if different from queries)
    if same_inputs == False and reference_kmer_dict == None:
        print("Processing references...")
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
    print("Done!")
    # ------------------------------------------------------

    # Calculate Jaccard distances
    # ------------------------------------------------------
    if index_db == True:
        print("Finished pre-indexing databases.")
        print("Next time you can run the program using only these files with --qd and(or) --rd.")
    else:
        print("Calculating shared Kmer fraction...")
        if same_inputs == True:
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
        with open(output, 'w') as outfile:
            for file in Fraction_Results:
                with open(file) as Temp:
                    shutil.copyfileobj(Temp, outfile)
                file.unlink()
        print("kAAI finishied correctly on {}".format(datetime.datetime.now()))
    # ------------------------------------------------------

if __name__ == "__main__":
    main()
