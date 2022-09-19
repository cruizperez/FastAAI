# FastAAI
Fast estimation of Average Amino Acid Identities (AAI) for bacterial and archaeal genomes.

## Content Table
  * [Features](#features)
  * [Citation](#citation)
  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Example](#example)
  * [Additional Information](#info)
  * [License](#license)

## Features
* 

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

## Example

Let's say we have a collection of genomes in a folder labeled "example_genomes" (which you can find in this respository). Each genome in the folder is in its own nucleotide FASTA-format file. The files can be gzipped or uncompressed - FastAAI doesn't care. The ones in the folder here are gzipped.

The first step is building a database. Here's an example command to do so:

```bash
fastaai build_db --genomes example_genomes/ --threads 4 --verbose --output example_build --database my_example_db.db --compress
```

This will create a folder called "example_build" which contains subfolders named "predicted_proteins," "hmms", "database", and "logs." The logs folder will contain a file named "FastAAI_preprocessing_log.txt," recording information about the protein prediction and HMM search results for each queryy genome. Finally, the database folder will contain "my_example_db.db," which is your completed FastAAI database.

Because we used the --compress flag, files in the predicted_proteins and hmms folders will be gzipped upon output, and because we used the --verbose flag, we'll get a progress report as FastAAI works that will look like so:

Completion |###############                               |   30.00% ( 3 of 10 ) at 19/09/2022 13:56:41

The report only updates every 2% completion, so it may be some time between updates if you're running hundreds or thousands of genomes. A build_db command will have two progress bars, one for preprocessing and one for database formatting, but they'll all look like so.

Next, we can calculate AAI:

```bash
fastaai db_query --query example_build/database/my_example_db.db --target example_build/database/my_example_db.db --threads 4 --verbose --output example_build
```

By supplying the same database as query and target, we'll be calculating an all vs. all AAI estimation for the genomes in the database. This will our all vs. all estimate for the genomes we had in our "example_genomes" folder.

We didn't supply --output_style matrix, so we'll be getting tabular output files. We also didn't tell FastAAI to calculate standard deviations with --do_stdev, so the fourth column will be all N/A. These files will be in example_build/results/, since we gave the same directory base as the output location.

When it's done (which should take less than a second), you'll find files that look like this:

query	target	avg_jacc_sim	jacc_SD	num_shared_SCPs	poss_shared_SCPs	AAI_estimate

_Pseudomonas__cissicola_GCA_002019225_1.fna.gz	Xanthomonas_albilineans_GCA_000962915_1.fna.gz	0.5199	N/A	79	79	68.75

_Pseudomonas__cissicola_GCA_002019225_1.fna.gz	Xanthomonas_albilineans_GCA_000962925_1.fna.gz	0.5176	N/A	79	79	68.63

_Pseudomonas__cissicola_GCA_002019225_1.fna.gz	Xanthomonas_albilineans_GCA_000962935_1.fna.gz	0.5193	N/A	79	79	68.72

_Pseudomonas__cissicola_GCA_002019225_1.fna.gz	Xanthomonas_albilineans_GCA_000962945_1.fna.gz	0.5189	N/A	79	79	68.7

...

That's it!

## Additional Information

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
* Kmers are represented in the database using 32-bit integers rather than as text. The integer representation of a particular tetramer can be found by finding the ASCII decimal values of each character and concatenating them in ordr, e.g. "KCMK" has ASCII values K = 75, C = 67, M = 77, K = 75, and thus is represented by the integer 75677775 in a database.

### Building a FastAAI Database

The build_db module of FastAAI is designed to take a set of inputs, preprocess them and build a database from those inputs. Input formatting is discussed above. This section discusses the components of a build.

#### Preprocessing:
Preprocessing consists of detecting the input format and (for genome inputs) predicting the proteome for each genome using [Pyrodigal](https://github.com/althonos/pyrodigal), (for genome and protein inputs) searching proteomes against FastAAI's set of SCP HMMs using PyHMMER (https://github.com/althonos/pyhmmer), and (for all input types) extracting the unique tetramer sets of each SCP protein identified as a bidirectional best-match by HMMER. The final result of preprocessing for each genome is a list SCPs, each with the unique amino acid tetramer set for the corresponding protein.

#### Identification of best-matching SCPs
An HMM search consists of searching protein sequences against a prebuilt model. Proteins matching models below a cutoff are not reported by FastAAI. Where FastAAI is concerned, three remaining features are important: the SCP accession, protein ID, and HMM score. The score of a protein indicates the quality of match to the associated SCP, where higher scores indicate better matches. 

Among the proteins that pass the initial filter, FastAAI identifies the highest scoring SCP assignment for each protein and the highest scoring assignment for each SCP. If a protein appears multiple times (that is, it matched to more than one SCP), then only the protein's highest scoring match is considered and the others are discarded. Likewise, if an SCP appears multiple times in the results (that is, multiple proteins matched to it), then only the highest scoring protein is retained for that SCP. Bidirectional best-matches are the remaining protein-SCP pairs. This all means that for each best-match between a protein and SCP, the protein's highest scoring match and the SCP's highest scoring match must be their counterpart in the pair.

A consequence of this approach is that each SCP can appear only once for each genome and each protein in a proteome can only be the representative for one SCP. A genome can have as little as one SCP or as many as 122 - more typically, a genome will have 50-90.

#### Database Construction:

After the input genomes are preprocessed, they are ready to be added to a database. The databases' genome_index table will be created or updated as needed to provide a numerical index of each genome, and metadata assosciated with the presence of each SCP and the count of kmers associated with it will be added to the genome_acc_kmer_counts table. Genome names are only represented as text in the genome_index table; in all other places, they are represented with integers according to the genome index.

FastAAI will add a record for each SCP in each genome to the corresponding SCP_genomes table, providing a record of the set of kmers associated with each genome in the database that is directly accessible using genome ID as a key. This genome-first representation is used by FastAAI when the database is used as a query.

FastAAI will then reorganize the data for each SCP into a kmer-first structure, listing each genome that contained a particular tetramer (e.g. KCMK) on the protein assosciated with a particular SCP (e.g. PF00119_20). This results in each tetramer being the key to a list of genome indices (e.g. table PF00119_20 , tetramer = KCMK, genomes = (0, 2, 5, 13, ...), where genomes 0, 2, 5, 13, etc., all have a representative protein that matched the HMM for SCP PF00119_20, and all of these representatives proteins contained the tetramer KCMK). This allows access to the set of genomes where a tetramer intersection would occur using a tetramer as a key. This tetramer-first representation is used by FastAAI when the database is used as a target.

Finally, tetramer tables are indexed if an index does not already exist. This is done to speed up the retrieval of tetramers during AAI calculation.

#### AAI Calculation

All of FastAAI's queries proceed essentially according to the same logic: the set of SCPs shared in common between a pair of genomes are selected, and the Jaccard index of unique tetramers is calculated for each shared SCP in the pairing. The average Jaccard index is calculated from the individual SCP pairings, unweighted, and the average Jaccard index is then transformed into an estimated AAI through an equation (see the FastAAI paper.)

To calculate the Jaccard index for a particular SCP, FastAAI selects one genome at a time as a query and sequentially searches each SCP in that genome against a target database. The organization of data discussed above allows FastAAI to select the tetramers associated with each of a query genomes' SCPs using the genome-first representation and then request all of the target genomes associated with each tetramer in the query using the tetramer-first representation. A tabulation of the number of appearances of each target genome produces the size of the tetramer intersection between the query and every target genome. 

The total number of tetramers associated with each SCP of each target genome are stored in the genome_acc_kmer_counts table, so the calculation of union size for each query and target pair is simply the sum of the number of the current query's tetramers and the number of tetramers in each target, minus the size of the intersection for the query and each target. Calculation of Jaccard index is trivial from here.

* size(union(Q, T)) = size(Q) + size(Q) - size(intersection(Q, T))
* Jaccard(Q, T) = size(intersection(Q, T))/size(union(Q, T))

#### Outputs

FastAAI allows for two primary output formats: tabular and matrix. The format of outputs is set using the "--output_style" argument. "--output_style tsv" produces a tab-separated output file for each query genome containing results for that query genome against all target genomes. These files have column headers which report:

* Query genome name
* Target genome name
* Average Jaccard index
* Jaccard index std. deviation ("N/A" unless --do_stdev is used)
* Count of shared SCPs
* Number of possibly shared SCPs (max number of SCPs in either member of the query-target genome pair)
* Estimated AAI.

All values other than query name and target name are "N/A" if a query-target paring share no SCPs.

The matrix format ("--output_style matrix") produces a tab-separated matrix containing query names in the first column, target names in the first row, and the final pairwise AAI estimate for each query-target pairing in the appropriate cells of the matrix. The matrix is complete, so there is a duplicate of all AAI estimates off of the main diagonal of the matrix. Further, there are some differences in the reporting of AAI when compared to the tabular format:

* The tabular format lists the AAI estimate for query-target pairs with no shared SCPs as "N/A." The matrix format reports these estimates as 0.0.
* The tabular format constrains AAI estimates between 30 <= AAI <= 90. AAI estimates <30% AAI are reported as "<30% AAI," rather than with a number, as are AAI estimates >90% AAI with the ">90% AAI" label. The matrix format reports these categorical estimates with 15.0 and 95.0 AAI, respectively.

These changes to avoid text labelling were done to make working with the (often quite large) tabular results in subsequent analyses using R, Python, or another language easier, as all of the results are already in numerical format. 

None of the other data aside from the estimated AAI is available in the matrix-formatted output.

## License

See LICENSE
