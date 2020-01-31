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
from pathlib import Path

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
    output_bacteria = folder / name.with_suffix('.bacteria.hmm')
    output_archaea = folder / name.with_suffix('.archaea.hmm')
    temp_output_bacteria = folder / name.with_suffix('.bacteria.temp')
    temp_output_archaea = folder / name.with_suffix('.archaea.temp')
    script_path = Path(__file__)
    script_dir = script_path.parent
    HMM_bacteria_model = script_dir / "00.Libraries/01.SCG_HMMs/Bacteria_SCG.hmm"
    HMM_archaea_model = script_dir / "00.Libraries/01.SCG_HMMs/Archaea_SCG.hmm"
    HMM_complete_model = script_dir / "00.Libraries/01.SCG_HMMs/Complete_CSG_DB.hmm"
    subprocess.call(["hmmsearch", "--tblout", str(output_bacteria), "-o", str(temp_output_bacteria), "--cut_ga", "--cpu", "1",
                    str(HMM_bacteria_model), str(file_path)])
    subprocess.call(["hmmsearch", "--tblout", str(output_archaea), "-o", str(temp_output_archaea), "--cut_ga", "--cpu", "1",
                    str(HMM_bacteria_model), str(file_path)])
    temp_output_bacteria.unlink()
    temp_output_archaea.unlink()

    counter_bacteria = 0
    counter_archaea = 0
    with open(output_bacteria) as hmm_bac:
        for line in hmm_bac:
            if line.startswith("#"):
                continue
            else:
                counter_bacteria += 1
    with open(output_archaea) as hmm_arc:
        for line in hmm_arc:
            if line.startswith("#"):
                continue
            else:
                counter_archaea += 1
    if counter_bacteria > counter_archaea:
        output_archaea.unlink()
        output = output_bacteria
    elif counter_bacteria < counter_archaea:
        output_bacteria.unlink()
        output = output_archaea
    else:
        print("Not sure which model to use. Default to complete dataset")
        output_archaea.unlink()
        output_bacteria.unlink()
        output_complete = folder / name.with_suffix('.complete.hmm')
        temp_output_complete = folder / name.with_suffix('.complete.temp')
        subprocess.call(["hmmsearch", "--tblout", str(output_complete), "-o", str(temp_output_complete), "--cut_ga", "--cpu", "1",
                    str(HMM_complete_model), str(file_path)])
        output = output_complete
    return output

# --- Find Kmers from HMM results ---
def Kmer_Parser(SCG_HMM_file, Keep):
    from pathlib import Path

    Kmer_Dic = {}
    HMM_Path = Path(SCG_HMM_file)
    Name = HMM_Path.name
    Folder = HMM_Path.parent
    Protein_File = Folder / HMM_Path.with_suffix('.faa')
    Positive_Matches = []
    with open(HMM_Path, 'r') as HMM_Input:
        for line in HMM_Input:
            if line.startswith("#"):
                continue
            else:
                Positive_Matches.append(line.strip().split()[0])
    if Keep == False:
        HMM_Path.unlink()
    kmers = read_kmers_from_file(Protein_File, Positive_Matches, 4)
    Kmer_Dic[Name] = set(kmers)

    return Kmer_Dic

# --- Build Kmers ---
def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

# --- Read Kmers from SCGs ---
def read_kmers_from_file(filename, positive_hits, ksize):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    all_kmers = []
    with open(filename) as Fasta_in:
        for title, sequence in SimpleFastaParser(Fasta_in):
            if title.split()[0] in positive_hits:
                kmers = build_kmers(sequence, ksize)
                all_kmers += kmers
    return all_kmers

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
def kAAI_Parser(ID):
    from pathlib import Path

    FilePath = Path(ID)
    Folder = Path.cwd()
    Output = Folder / FilePath.with_suffix('.aai.temp')
    with open(Output, 'w') as OutFile:
        for key2, value2 in Kmer_Dictionary.items():
            intersection = len(Kmer_Dictionary[ID].intersection(value2))
            shorter = min(len(list(Kmer_Dictionary[ID])), len(list(value2)))
            fraction = round(intersection/shorter, 3)
            OutFile.write("{}\t{}\t{}\t{}\t{}\n".format(ID, key2, intersection, shorter, fraction))
    return Output

# --- Initialize function ---
def child_initialize(_dictionary):
     global Kmer_Dictionary
     Kmer_Dictionary = _dictionary

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
    import argparse
    from sys import argv
    from sys import exit
    from pathlib import Path
    import subprocess
    import multiprocessing
    from functools import partial
    import datetime
    import shutil

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
    parser.add_argument('-k', '--keep', dest='Keep', action='store_false', required=False, help='Keep intermediate files, by default true')
    args = parser.parse_args()

    Genome_List = args.Genome_List
    Protein_Files = args.Protein_Files
    HMM_Files = args.HMM_Files
    Output = args.Output
    Threads = args.Threads
    Keep = args.Keep

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


    # Parse HMM results, calculate distances and compile results
    print("Parsing HMM results...")
    print(datetime.datetime.now()) # Remove after testing
    try:
        pool = multiprocessing.Pool(Threads)
        Kmer_Results = pool.map(partial(Kmer_Parser, Keep=Keep), HMM_Search_Files)
    finally:
        pool.close()
        pool.join()

    Final_Kmer_Dict = merge_dicts(Kmer_Results)

    # Calculate shared Kmer fraction
    print("Calculating shared Kmer fraction...")
    print(datetime.datetime.now()) # Remove after testing
    ID_List = Final_Kmer_Dict.keys()
    try:
        pool = multiprocessing.Pool(Threads, initializer = child_initialize, initargs = (Final_Kmer_Dict,))
        Fraction_Results = pool.map(kAAI_Parser, ID_List)
    finally:
        pool.close()
        pool.join()

     # Merge results into a single output
    print(datetime.datetime.now()) # Remove after testing
    with open(Output, 'w') as OutFile:
        for file in Fraction_Results:
            with open(file) as Temp:
                shutil.copyfileobj(Temp, OutFile)
            file.unlink()


if __name__ == "__main__":
    main()