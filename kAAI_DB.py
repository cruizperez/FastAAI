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
import sqlite3
import pickle
from contextlib import closing


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

    return final_output

# --- Run hmmsearch ---
def run_hmmsearch(input_file):
    """
    Runs hmmsearch on the set of SCGs and stores the hmm results
    
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
        HMM_path.unlink()
    return outfile



# --- Find Kmers from HMM results ---
def Kmer_parser(SCG_HMM_file_filtered, database):
    """
    Extract kmers from protein files that have hits
    in the HMM searches.
    
    Arguments:
        SCG_HMM_file_filtered {file path} -- Path to filtered HMM results.
    
    Returns:
        [genome_kmers] -- Dictionary of kmers per gene. 
    """
    HMM_path = Path(SCG_HMM_file_filtered)
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
    SCG_HMM_file_filtered.unlink()
    pick_dic = pickle.dumps(genome_kmers, pickle.HIGHEST_PROTOCOL)
    return (genome, pick_dic)

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

# --- Look in DB from kmer info and calculate kAAI ---
def kAAI_calculator(query_id, list_and_db):
    genomes = list_and_db[0]
    database = list_and_db[1]
    file_path = Path(query_id)
    running_folder = Path.cwd()
    temp_output = running_folder / file_path.with_suffix('.aai.temp')
    conn = sqlite3.connect(database)
    conn.text_factory = str
    cur = conn.cursor()
    out_file = open(temp_output, 'w')

    sql_comm = 'SELECT kmer_info FROM Genome_Kmers WHERE genome = ?;'
    cur.execute(sql_comm, (query_id,))
    data = cur.fetchone()[0]
    query_id_kmer = pickle.loads(data)
    query_id_kmer = pickle.loads(query_id_kmer)
    query_num_scg = len(query_id_kmer[query_id])
    query_scg_list = query_id_kmer[query_id].keys()
    for target_id in genomes:
        comp = query_id + "___" + target_id
        sql_comm = 'SELECT COUNT(*) FROM Genome_Comparisons WHERE comparison = ?;'
        cur.execute(sql_comm, (comp,))
        result = cur.fetchone()[0]
        if result == 0:
            jaccard_similarities = []
            sql_comm = 'SELECT kmer_info FROM Genome_Kmers WHERE genome = ?;'
            cur.execute(sql_comm, (target_id,))
            data = cur.fetchone()[0]
            target_id_kmer = pickle.loads(data)
            target_id_kmer = pickle.loads(target_id_kmer)
            target_num_scg = len(target_id_kmer[target_id])
            target_scg_list = target_id_kmer[target_id].keys()
            if query_num_scg > target_num_scg:
                final_scg_list = target_scg_list
            else:
                final_scg_list = query_scg_list
            for accession in final_scg_list:
                if accession in query_scg_list and accession in target_scg_list:
                    kmers_query = query_id_kmer[query_id][accession]
                    kmers_target = target_id_kmer[target_id][accession]
                    intersection = len(kmers_query.intersection(kmers_target))
                    union = len(kmers_query.union(kmers_target))
                    jaccard_similarities.append(intersection / union)
                else:
                    continue
            try:
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_id,
                           sum(jaccard_similarities)/len(jaccard_similarities),
                           len(jaccard_similarities), len(final_scg_list)))
            except:
                out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_id,
                           "NA", "NA", "NA"))
        else:
            continue


    #     for target_genome, scg_ids in total_kmer_dictionary.items():
    #         jaccard_similarities = []
    #         target_num_scg = len(scg_ids)
    #         target_scg_list = scg_ids.keys()
    #         if query_num_scg > target_num_scg:
    #             final_scg_list = target_scg_list
    #         else:
    #             final_scg_list = query_scg_list
    #         for accession in final_scg_list:
    #             if accession in query_scg_list and accession in target_scg_list:
    #                 kmers_query = total_kmer_dictionary[query_id][accession]
    #                 kmers_target = total_kmer_dictionary[target_genome][accession]
    #                 intersection = len(kmers_query.intersection(kmers_target))
    #                 union = len(kmers_query.union(kmers_target))
    #                 jaccard_similarities.append(intersection / union)
    #             else:
    #                 continue
    #         try:
    #             out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
    #                        sum(jaccard_similarities)/len(jaccard_similarities),
    #                        len(jaccard_similarities), len(final_scg_list)))
    #         except:
    #             out_file.write("{}\t{}\t{}\t{}\t{}\n".format(query_id, target_genome,
    #                        "NA", "NA", "NA"))
    return temp_output

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


################################################################################
"""---2.0 Main Function---"""

def main():
    # Setup parser for arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script calculates the average amino acid identity using k-mers\n'''
                        '''from single copy genes. It is a faster version of the regular AAI\n'''
                        '''(Blast or Diamond) and the hAAI implemented in MiGA.\n\n'''
            '''Usage: ''' + argv[0] + ''' -p [Protein Files] -t [Threads] -o [Output]\n'''
            '''Global mandatory parameters: -g [Genome Files] OR -p [Protein Files] OR -s [SCG HMM Results] -o [AAI Table Output]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')

    input_options = parser.add_argument_group('General input/output options')
    input_options.add_argument('-g', '--genomes', dest='genome_list', action='store', nargs='+', required=False,
                        help='List of input genomes. Implies step 1')
    input_options.add_argument('--gql', dest='Genome_List', action='store', nargs='+', required=False,
                        help='List of input genomes. Implies step 1')
    input_options.add_argument('-p', '--proteins', dest='protein_files', action='store', nargs='+', required=False,
                        help='List of input protein files')
    input_options.add_argument('-s', '--scg_hmm', dest='HMM_files', action='store', nargs='+', required=False,
                        help='List of hmm search results')
    input_options.add_argument('-o', '--output', dest='output', action='store', required=True, help='Output File')
    
    database_options = parser.add_argument_group('SQL database options')
    database_options.add_argument('--database', dest='database', action='store', default=None, required=False,
                        help='''Database with pre-computed distances. Note that the entry names\n'''
                             '''in the database are determined by the file name, e.g., "GCA_003019315.1.faa"\n'''
                             '''will be stored as "GCA_003019315.1", regardless of the extension. Choice between:\n'''
                             '''refseq: set of complete refseq bacteria/archaea genomes\n'''
                             '''[db file]: custom database file location''')
    database_options.add_argument('--db_output', dest='database_output', action='store', required=False,
                        help='''File to save updated database file. By default don't save\n''')
    database_options.add_argument('--recalculate', dest='recalculate', action='store_true', required=False,
                        help='''Recalculate distances already stored in the database. By default False\n''')

    
    misc_options = parser.add_argument_group('Miscellaneous options')
    misc_options.add_argument('-t', '--threads', dest='threads', action='store', default=1, type=int, required=False,
                        help='Number of threads to use, by default 1')
    misc_options.add_argument('-k', '--keep', dest='keep', action='store_false', required=False,
                        help='Keep intermediate files, by default true')

    args = parser.parse_args()
    
    genome_list = args.genome_list
    protein_files = args.protein_files
    HMM_files = args.HMM_files
    output = args.output
    database = args.database
    database_output = args.database_output
    recalculate = args.recalculate
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

    # Use requested database or create new one
    # ------------------------------------------------------
    # Table names
    tab1_name = "Genome_Kmers"
    tab2_name = "Genome_Comparisons"
    tab1_index = "genome_index"
    tab2_index = "comparison_index"
    col_names_tab1 = ("genome", "kmer_info")
    col_names_tab2 = ("comparison", "jaccard_dist", "prots_used", "prots_short")
    conn = None
    cur = None

    if database == None:
        print("Building temporal database. {}".format(datetime.datetime.now()))
        conn = sqlite3.connect("Temporal_DB.sql")
        conn.text_factory = str  # allows utf-8 data to be stored
        cur = conn.cursor()
        # Create table for genome information
        sql = 'DROP TABLE IF EXISTS ' + tab1_name
        cur.execute(sql)
        sql = 'CREATE TABLE Genome_Kmers ("genome", "kmer_info")'
        cur.execute(sql)
        # Create index for faster search
        index = tab1_index
        sql = 'CREATE UNIQUE INDEX genome_index ON Genome_Kmers (genome);'
        cur.execute(sql)
        # Create table for comparisons
        sql = 'DROP TABLE IF EXISTS ' + tab2_name
        cur.execute(sql)
        sql = 'CREATE TABLE Genome_Comparisons ("comparison", "jaccard_dist", "prots_used", "prots_short")'
        cur.execute(sql)
        # Create index for faster search
        index = "%s" % (tab2_index)
        sql = 'CREATE UNIQUE INDEX comparison_index ON Genome_Comparisons (comparison);'
        cur.execute(sql)
    elif database == "refseq":
        script_path = Path(__file__)
        script_dir = script_path.parent
        db_location = str(script_dir / "00.Libraries/02.SB/refseq.sql")
        conn = sqlite3.connect(db_location)
        conn.text_factory = str
        cur = conn.cursor()
    else:
        conn = sqlite3.connect(database)
        conn.text_factory = str
        cur = conn.cursor()
    # ------------------------------------------------------

    # Predict proteins and perform HMM searches
    # ------------------------------------------------------
    original_input_files = []
    HMM_search_files = []
    if database == None or recalculate == True:
        if genome_list != None:
            print("Starting from Genomes...")
            print("Predicting proteins...")
            protein_files = []
            original_input_files = genome_list
            try:
                pool = multiprocessing.Pool(threads)
                protein_files = pool.map(run_prodigal, genome_list)
            finally:
                pool.close()
                pool.join()
            print("Searching HMM models...")
            try:
                pool = multiprocessing.Pool(threads)
                HMM_search_files = pool.map(run_hmmsearch, protein_files)
            finally:
                pool.close()
                pool.join()
        elif protein_files != None:
            print("Starting from Proteins...")
            print("Searching HMM models...")
            original_input_files = protein_files
            try:
                pool = multiprocessing.Pool(threads)
                HMM_search_files = pool.map(run_hmmsearch, protein_files)
            finally:
                pool.close()
                pool.join()
        elif HMM_files != None:
            original_input_files = HMM_files
            print("Starting from HMM searches...")
            HMM_search_files = HMM_files
    else:
        if genome_list != None:
            original_input_files = genome_list
            print("Starting from Genomes...")
            print("Predicting proteins...")
            missing_genomes = []
            for genome in genome_list:
                genome_name = str(Path(genome).with_suffix(''))
                sql_comm = 'SELECT COUNT(*) FROM Genome_Kmers WHERE genome = ?;'
                cur.execute(sql_comm, (genome_name,))
                result = cur.fetchone()[0]
                if result == 0:
                    missing_genomes.append(genome)
                else:
                    continue
            protein_files = []
            try:
                pool = multiprocessing.Pool(threads)
                protein_files = pool.map(run_prodigal, missing_genomes)
            finally:
                pool.close()
                pool.join()
            print("Searching HMM models...")
            try:
                pool = multiprocessing.Pool(threads)
                HMM_search_files = pool.map(run_hmmsearch, protein_files)
            finally:
                pool.close()
                pool.join()
        elif protein_files != None:
            print("Starting from Proteins...")
            print("Searching HMM models...")
            original_input_files = protein_files
            missing_proteins = []
            for protein in protein_files:
                protein_name = str(Path(protein).with_suffix(''))
                sql_comm = 'SELECT COUNT(*) FROM Genome_Kmers WHERE genome = ?;'
                cur.execute(sql_comm, (protein_name,))
                result = cur.fetchone()[0]
                if result == 0:
                    print(protein_name, "TEST")
                    missing_proteins.append(protein)
                else:
                    continue
            try:
                pool = multiprocessing.Pool(threads)
                HMM_search_files = pool.map(run_hmmsearch, protein_files)
            finally:
                pool.close()
                pool.join()
        elif HMM_files != None:
            print("Starting from HMM searches...")
            original_input_files = HMM_files
            for hmm_search in HMM_files:
                name = Path(hmm_search).name
                hmm_name = str(Path(name).with_suffix(''))
                sql_comm = 'SELECT COUNT(*) FROM Genome_Kmers WHERE genome = ?;'
                cur.execute(sql_comm, (hmm_name,))
                result = cur.fetchone()[0]
                if result == 0:
                    HMM_search_files.append(hmm_search)
                else:
                    continue

    # ------------------------------------------------------

    # Filter HMM results, retaining best hit per protein
    # ------------------------------------------------------
    print("Filtering HMM results...")
    print(datetime.datetime.now())
    filtered_files = []
    try:
        pool = multiprocessing.Pool(threads)
        filtered_files = pool.map(partial(HMM_filter, keep=keep), HMM_search_files)
    finally:
        pool.close()
        pool.join()
    # ------------------------------------------------------
    
    # Find kmers per SCG per genome and populate genome information table in DB
    # ------------------------------------------------------
    print("Parsing HMM results...")
    print(datetime.datetime.now())
    if database == None:
        database = "Temporal_DB.sql"
    try:
        pool = multiprocessing.Pool(threads)
        kmer_results = pool.map(partial(Kmer_parser, database=database), filtered_files)
    finally:
        pool.close()
        pool.join()
    print("Populating database")
    for genome in kmer_results:
        pick_dic = pickle.dumps(genome[1], pickle.HIGHEST_PROTOCOL)
        insertsql = "INSERT INTO " + tab1_name + " VALUES(?,?);"
        cur.execute(insertsql, (genome[0], sqlite3.Binary(pick_dic)))
    # ------------------------------------------------------
    
    conn.commit()
    cur.close()
    conn.close()
    
    # Calculate shared Kmer fraction
    # ------------------------------------------------------
    print("Calculating shared Kmer fraction...")
    print(datetime.datetime.now())
    genomes_to_compare = []
    for genome in original_input_files:
        name = Path(genome).name
        file_name = str(Path(name).with_suffix(''))
        genomes_to_compare.append(file_name)
    list_and_db = (genomes_to_compare, database)
    try:
        pool = multiprocessing.Pool(threads)
        comparison_results = pool.map(partial(kAAI_calculator, list_and_db=list_and_db), genomes_to_compare)
    finally:
        pool.close()
        pool.join()

    # for i in genomes_to_compare:
    #     kAAI_calculator(i, (genomes_to_compare, database))
    
    # ------------------------------------------------------

    # # Calculate shared Kmer fraction
    # print("Calculating shared Kmer fraction...")
    # print(datetime.datetime.now())
    # ID_List = Final_Kmer_Dict.keys()
    # try:
    #     pool = multiprocessing.Pool(Threads, initializer = child_initialize, initargs = (Final_Kmer_Dict,))
    #     Fraction_Results = pool.map(kAAI_Parser, ID_List)
    # finally:
    #     pool.close()
    #     pool.join()

    # Merge results into a single output
    print("Merging results...")
    print(datetime.datetime.now())
    with open(output, 'w') as outFile:
        for file in comparison_results:
            with open(file) as temp:
                shutil.copyfileobj(temp, outFile)
            file.unlink()

if __name__ == "__main__":
    main()
