# Recycler2

Recycler2 assembles plasmids from metagenomic assembly graphs.

## Installation

Recycler2 is written in Python3. Recycler2 uses NumPy, NetworkX, pySAM, and nose. The required versions of these required dependencies will all be installed by the setup.py script.

Recycler2 uses [BWA](https://github.com/lh3/bwa), [NCBI BLAST+ tools](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), and [samtools](https://github.com/samtools/samtools). The executables of these programs should be available on the system on which Recycler2 is run.

The [PlasClass classifier](https://github.com/Shamir-Lab/PlasClass) should also be installed in order to use the full functionality of Recycler2.

We recommend using a virtual environment. For example, in Linux, before running setup.py:
```
python -m venv recycler2-env
source recycler2-env/bin/activate
```
To install, download and run setup.py:
```
    git clone https://github.com/Shamir-Lab/Recycler2.git
    cd Recycler2
    python setup.py install
```
It is possible to install as a user without root permissions:
```
python setup.py install --user
```
#### Configuring paths to required executables
The BWA, samtools, and BLAST+ executables must be available to Recycler2. They can either be added to your `PATH` environment variable, or you can specify the paths to each of them in the file `bin/config.json`.
