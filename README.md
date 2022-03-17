# FastAAI
Fast estimation of Average Amino Acid Identities (AAI) for bacterial and viral genomes.
Includes a module for the classification of viral genomes.

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
- Programs:
   - [HMMER](http://hmmer.org/) >= 3.1
- Python >=3.6,<3.9
- Base Python Modules:
   - argparse
   - datetime
   - pathlib
   - shutil
   - subprocess
   - gzip
   - multiprocessing
   - textwrap
   - pickle
   - tempfile
   - sys
   - functools
- Additional Python Modules:
   - numpy

## Installation

### Conda Installation
It appears we need a bunch of pre-requisites to run FastAAI. No worries, their installation using Conda is quite easy. If you don't have Conda, you can install it as follows:
1. Download Anaconda from https://www.anaconda.com/products/individual.
2. Run `bash Anaconda-latest-Linux-x86_64.sh` and follow the installation instructions.
3. Once installed you can run `conda -V`. You should get the version of conda that you installed.


Now, let's add the conda channels required to install the pre-requisites:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
```

Then, create an environment for FastAAI:

```bash
conda create -n fastaai hmmer prodigal numpy python=3.7
```

And activate it:

```bash
conda activate fastaai
```

### Pip Installation

Final installation of FastAAI should be done with pip. Once you have the  FastAAI with the following command:

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

### Database creation


## FAQs
Coming soon


## License

See LICENSE
