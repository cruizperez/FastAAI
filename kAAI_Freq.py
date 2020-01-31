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
"""---1.0 Define Functions---"""
# --- Run prodigal ---
def run_prodigal(InputFile):
    import subprocess
    from pathlib import Path

    FilePath = Path(InputFile)
    Prefix = Path(FilePath.stem)
    Folder = FilePath.parent
    Output = Folder / Prefix.with_suffix('.faa')
    Temp_Output = Folder / Prefix.with_suffix('.temp')
    subprocess.call(["prodigal", "-i", str(FilePath), "-a", str(Output), "-p", "meta", "-q", "-o", str(Temp_Output)])
    Temp_Output.unlink()
    return Output

# --- Run hmmsearch ---
def run_hmmsearch(InputFile):
    import subprocess
    from pathlib import Path

    FilePath = Path(InputFile)
    Folder = FilePath.parent
    Name = Path(FilePath.name)
    Output = Folder / Name.with_suffix('.hmm')
    Temp_Output = Folder / Name.with_suffix('.temp')
    Script_path = Path(__file__)
    Script_dir = Script_path.parent
    HMM_Model = Script_dir / "00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm"
    subprocess.call(["hmmsearch", "--tblout", str(Output), "-o", str(Temp_Output), "--cut_ga", "--cpu", "1", str(HMM_Model), str(FilePath)])
    Temp_Output.unlink()
    return Output

# --- Find Kmers from HMM results ---
def Kmer_Parser(SCG_HMM_file, Keep):
    from pathlib import Path

    Kmer_Dic = {}
    HMM_Path = Path(SCG_HMM_file)
    Name = HMM_Path.name
    Folder = HMM_Path.parent
    Protein_File = Folder / Path(Name).with_suffix('.faa')
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
    Name = Name.replace('.hmm', '')
    Kmer_Dic[Name] = kmers

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
                sequence = sequence.replace('*','') # Remove asterisk from prodigal
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

#! New Function to Filter Abundant Kmers
def Kmer_Filter(Kmer_Dictionary, Output):
    """
    Counts frequency of kmers, fits an inverse
    gamma distribution and calculates one stdev
    to filter kmers higher than this value
    
    Arguments:
        Kmer_Dictionary {[Dictionary]} -- Dictionary of Kmers
        Output {[String]} -- Output for histogram plot
    """
    from collections import Counter
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import invgauss
    import pandas as pd
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Century Gothic'
    plt.rcParams.update({'font.size': 25})

    Total_Kmers = []
    for kmers in Kmer_Dictionary.values():
        Total_Kmers += kmers
    Total_Frequency_Count = Counter(Total_Kmers)
    Total_Frequency = pd.DataFrame.from_dict(Total_Frequency_Count, orient='index')
    Freqs = Total_Frequency.iloc[:,0]
    Freqs = np.asarray(Freqs)
    #* Fit and plot data
    mu, loc, scale = invgauss.fit(Freqs)
    X_values = np.linspace(Freqs.min(), Freqs.max(), 500)
    Y_Predicted = invgauss.pdf(X_values, mu, loc, scale)
    Std = invgauss.var(mu, loc, scale)**0.5
    Mean = invgauss.mean(mu, loc, scale)
    Figure, Axis = plt.subplots(1,1,figsize=(10,10),dpi=300,constrained_layout=True)
    Axis.hist(Freqs, bins=50, density=True, color="#4E8983", label="Kmer Frequencies")
    Axis.plot(X_values, Y_Predicted, linewidth=2, color="#217B55", label="InvGamma Distribuion")
    Axis.axvline(x=Mean, linewidth=2, color='#960D41', label="Kmer Frequency Cutoff")
    Axis.legend(loc="best")
    Axis.set_xlabel("Number of occurrences in dataset")
    Axis.set_ylabel("Frequency")
    Figure.savefig(Output + "_plot.png")
    Low_Abundance_Kmers = []
    print(Mean, Std)
    for key, value in Total_Frequency_Count.items():
        #if value >= (int(Mean) + int(Std)*2):
        # TESTING PURPOSES
        if ((value >= (int(Mean) + int(Std)*2)) and (value <= (int(Mean) + int(Std)*3))):
            Low_Abundance_Kmers.append(key)
    Low_Abundance_Kmers.sort()
    
    return Low_Abundance_Kmers


#! Calculate frequencies individually and return a dataframe
def Frequency_Calculator(Kmer_Dictionary, Low_Abundance_Kmers):
    """Calculates the frequencies of Kmers in each genome
       and returns a table
    
    Arguments:
        Kmer_Dictionary {Dictionary} -- Dictionary of kmers
    """
    import pandas as pd
    from collections import Counter
    import numpy as np

    Genome_ID = None
    for key, value in Kmer_Dictionary.items():
        Counts = dict(Counter(value))
        Genome_ID = key

    Freq_Table = pd.DataFrame.from_records([Counts])
    Freq_Table.rename(index={0:Genome_ID}, inplace=True)
    Freq_Table.filter(items=Low_Abundance_Kmers)
    Freq_Table = Freq_Table.reindex(Low_Abundance_Kmers, axis=1).fillna(0)
    Frequencies = Freq_Table.to_numpy()[0]
    Frequencies = Frequencies.astype(int)

    return Genome_ID, Frequencies

#! Calculate bray curtis distances?
def Distance_Calculator(Kmer_Freqs, Genomes):
    """
    Calculates bray curtis distances between
    each pair of genomes
    
    Arguments:
        Kmer_Freqs {Array} -- Dataframe with kmer frequencies
        Genomes {List} -- List of genome names
    """
    from skbio.diversity import beta_diversity

    BC_Distances = beta_diversity("braycurtis", Kmer_Freqs, Genomes)
#    BC_Distances = beta_diversity("euclidean", Kmer_Freqs, Genomes)
    BC_Dataframe = BC_Distances.to_data_frame()
    BC_Dataframe = 1-BC_Dataframe
    
    return BC_Dataframe

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
    import pandas as pd
    import numpy as np

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
    print("kAAI started on {}".format(datetime.datetime.now()))
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
    print(datetime.datetime.now())
    try:
        pool = multiprocessing.Pool(Threads)
        Kmer_Results = pool.map(partial(Kmer_Parser, Keep=Keep), HMM_Search_Files)
    finally:
        pool.close()
        pool.join()
    Final_Kmer_Dict = merge_dicts(Kmer_Results)

    # Estimate total frequency of kmers, plot, and return threshold for filter
    print("Filtering kmers by frequency...")
    print(datetime.datetime.now())
    Low_Abundance_Kmers = Kmer_Filter(Final_Kmer_Dict, Output)


    print("Calculating kmer frequencies per genome...")
    print(datetime.datetime.now())
    try:
        pool = multiprocessing.Pool(Threads)
        Frequency_Results = pool.map(partial(Frequency_Calculator, Low_Abundance_Kmers=Low_Abundance_Kmers),
         Kmer_Results)
    finally:
        pool.close()
        pool.join()
    
    print("Merging tables and filtering...")
    print(datetime.datetime.now())
    Genomes = []
    Arrays = []
    for i in Frequency_Results:
        Genomes.append(i[0])
        Arrays.append(i[1])
    Array_Kmers = np.vstack(Arrays)
    # Final_Dataframe = pd.DataFrame(Array_Kmers, index=Genomes)

    print("Calculating distances...")
    print(datetime.datetime.now())
    print(Array_Kmers)
    BC_Distances = Distance_Calculator(Array_Kmers, Genomes)
    BC_Distances.to_csv(Output + "_dist.mat", sep="\t", header=True, index=True)

    with open(Output + "_dist.tab", 'w') as Outfile:
        for i in BC_Distances.columns:
            for j in BC_Distances.columns:
                Outfile.write("{}\t{}\t{}\n".format(i, j, BC_Distances.loc[i,j]))

    print("kAAI finished on {}".format(datetime.datetime.now())) # Remove after testing

if __name__ == "__main__":
    main()
