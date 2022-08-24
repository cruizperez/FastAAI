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

### Conda Installation
It appears we need a bunch of pre-requisites to run FastAAI. No worries, their installation using Conda is quite easy. If you don't have Conda, you can install it as follows:
1. Download Anaconda from https://www.anaconda.com/products/individual.
2. Run `bash Anaconda-latest-Linux-x86_64.sh` and follow the installation instructions.
3. Once installed you can run `conda -V`. You should get the version of conda that you installed.

### Pip Installation

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

### Database creation

The build_db module of FastAAI is designed to take a (1) set of inputs, (2) preprocess them and (3) build a database from those inputs.

#### (1) Inputs and their formats



## FAQs
Coming soon


## License

See LICENSE
