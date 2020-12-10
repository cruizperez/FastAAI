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
- Python Modules:
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
   - numpy
   - sys
   - functools

## Installation
### Conda Installation
FastAAIIt appears we need a bunch of pre-requisites to run FastAAI No worries, their installation using Conda is quite easy. If you don't have Conda, you can install it as follows:
1. Download Anaconda from https://www.anaconda.com/products/individual.
2. Run `bash Anaconda-latest-Linux-x86_64.sh` and follow the installation instructions.
3. Once installed you can run `conda -V`. You should get the version of conda that you installed.

Now, let's add the conda channels required to install the pre-requisites:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels cruizperez
```

Then, create an environment for MicrobeAnnotator:

```bash
conda create -n fastaai hmmer prodigal numpy python=3.7 fastaai
```

And activate it:

```bash
conda activate microbeannotator
```

Both main scripts (microbeannotator and microbeannotator_db_builder) should be in your path ready for use!
This should take care of most of the requirements except for Aspera Connect and KofamScan, which are a little more involved. Let's install those.

### Pip Installation
#Once you have installed the pre-requisites to run MicrobeAnnotator, or if you already had them and you are not using Conda, you can install MicrobeAnnotator using pip:


## Usage
### Database creation


## FAQs



## License
See LICENSE