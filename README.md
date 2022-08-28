# FastAAI
Fast estimation of Average Amino Acid Identities (AAI) for bacterial and archaeal genomes.

## Content Table
  * [Features](#features)
  * [Citation](#citation)
  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [FAQs](#faqs)
  * [License](#license)

## Features
Coming soon

## Citation
Coming soon

## Requirements:
- Python >=3.6 (3.9+ recommended)
- Additional Python Modules:
   - numpy
   - pyrodigal - https://github.com/althonos/pyrodigal/
   - pyhmmer - https://github.com/althonos/pyhmmer

## Installation

FastAAI and all its dependencies can be installed through pip with the following command:

```bash
pip install FastAAI
```

## Usage

FastAAI executes its behaviors through commands. A list of commands and their behaviors can be seen through simply calling FastAAI on the command line:

```bash
fastaai
```

The various commands each contain their own usage instructions, which can be accessed through calling fastaai [command], e.g.

```bash
fastaai build_db
```

The five FastAAI commands are

* build_db - Input a set of genomes and predict proteins, identify single-copy proteins, and construct (or add to) a FastAAI database.
* merge_db - Merge two or more FastAAI databases. Can create a new database or modify an existing one.
* simple_query - Input a set of genomes as a query and a prebuilt FastAAI database as a target; calculate AAI for each query against each target
* db_query - Query the genomes in one FastAAI database against the genomes in another (or itself). Calculate AAI for each genome pair between the two.
* single_query - Input exactly two genomes; preprocess as needed and calculate AAI between the pair of genomes.

## General Topics

#### Input files and their formats
FastAAI takes genomes, proteins, and tabular HMM search files (see below) as its basic inputs. Genomes are expected to be supplied in nucleotide FASTA format, with each genome (even if they are collections of multiple contigs) to be in a single, separate file. Each protein file is expected to contain the predicted proteome of a single genome in amino acid FASTA format. Each HMM file is expected to be the tabular output resulting from a search of a single genome's proteome against FastAAI's reference set of HMM models.

Inputs of each type (genome, protein, HMM) can be given in one of three ways:
* As a directory: a path to a directory containing only files of a particular type
* As a file: a text file containing paths to only files of a given type, one path per line
* As a string: a comma-separated string of paths to files; note that Python will give up if there are too many files given this way.

FastAAI will automatically detect the way you supply inputs.

#### HMM Files
In this repositoriy, you can find FastAAI's single copy protein HMM models under the heading of 00.Libraries/01.SCG_HMMs/Complete_SCG_DB.hmm. An HMM file produced by FastAAI will be the result of a search against this collection of models. While FastAAI uses [PyHMMER](https://github.com/althonos/pyhmmer) to implement its HMM search, any FastAAI HMM search can also be replicated with the following command:

```bash
hmmsearch --tblout [FastAAI_HMM_file] -o [file_to_dsiscard] --cut_tc Complete_SCG_DB.hmm [input_proteome_file]
```

#### Databases
FastAAI databases are SQLite 3 databases containing collections of genomes. Within every FastAAI database there will be two metadata tables describing the genomes the database contains and an additional two tables for each SCP observed within the set of genomes the database contains, up to 244 possible tables. The two metadata tables will always be named "genome_index" and "genome_acc_kmer_counts" with the paired, per-SCP tables named as "SCP_Accession_ID" and "SCP_Accession_ID_genomes", with all characters that have syntactic meaning in SQL replaced with underscores. An example of an SCP table pair would be for the SCP PF00119.20, which becomes PF00119_20 after character changes and will produce the tables PF00119_20 and PF00119_20_genomes. Schemas follow:

* genome_index (genome TEXT, gen_id INTEGER, protein_count INTEGER)
* genome_acc_kmer_counts (genome INTEGER, accession INTEGER, count INTEGER)
* [accession_ID] (kmer INTEGER, genomes BLOB)
* [accession_ID_genomes] (genome INTEGER, kmers BLOB)

Additional notes:

* Accession IDs are numbered according to an internal scheme used by FastAAI. Numbering for the IDs is available on this github in the metadata folder.
* Genomes are always numbered from 0 to (#genomes - 1) within a FastAAI database. No indices may be absent; FastAAI will refuse to add a genome if the genome contains no data FastAAI can use.
* Both the genomes blob and kmers blob are bytestrings of 32 bit integers.

### Database creation

The build_db module of FastAAI is designed to take a set of inputs, preprocess them and build a database from those inputs. Input formatting is discussed above.

#### Preprocessing:
Preprocessing consists of detecting the input format and (for genome inputs) predicting the proteome for each genome using [Pyrodigal](https://github.com/althonos/pyrodigal), (for genome and protein inputs) searching proteomes against FastAAI's set of SCP HMMs using PyHMMER, and (for all input types) extracting the unique tetramer sets of each SCP protein identified as a bidirectional best-match by HMMER. The final result of preprocessing for each genome is a list SCPs, each with the unique amino acid tetramer set for the corresponding protein.

#### Identification of best-matching SCPs
An HMM search consists of searching protein sequences against a prebuilt model. Proteins matching models below a cutoff are not reported by FastAAI. Where FastAAI is concerned, three remaining features are important: the SCP accession, protein ID, and HMM score. The score of a protein indicates the quality of match to the associated SCP, where higher scores indicate better matches. 

Among the proteins that pass the initial filter, FastAAI identifies the highest scoring SCP assignment for each protein and the highest scoring assignment for each SCP. If a protein appears multiple times (that is, it matched to more than one SCP), then only the protein's highest scoring match is considered and the others are discarded. Likewise, if an SCP appears multiple times in the results (that is, multiple proteins matched to it), then only the highest scoring protein is retained for that SCP. Bidirectional best-matches are the remaining protein-SCP pairs. This all means that for each best-match between a protein and SCP, the protein's highest scoring match and the SCP's highest scoring match must be their counterpart in the pair.

A consequence of this approach is that each SCP can appear only once for each genome and each protein in a proteome can only be the representative for one SCP. A genome can have as little as one SCP or as many as 122 - more typically, a genome will have 50-90.

#### Database construction:



## FAQs
Coming soon


## License

See LICENSE
