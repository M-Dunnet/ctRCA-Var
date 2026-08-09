## ctRCA - Part 1: Software Installation

### Overview
The ctDNA Rolling Circle Amplification (ctRCA) pipeline requires the following publicly available software:

- [SAMTools](https://github.com/samtools/samtools)
- [Dorado](https://dorado-docs.readthedocs.io/en/latest/#__tabbed_1_2)
- [BEDTools](https://bedtools.readthedocs.io/en/latest/)
- [BB-8](https://github.com/M-Dunnet/BB-8)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [minimap2](https://github.com/lh3/minimap2)
- [ctRCA-Var](https://github.com/M-Dunnet/ctRCA-Var)

### Directory Setup
Install all software into a single directory:
```bash=
mkdir -p ~/software/ctRCA
cd ~/software/ctRCA
```
For the rest of this guide, this directory will be referred to as:
```bash=
export CTRCA_SOFT=~/software/ctRCA
```

### Python Virtual Environment
Several components require Python dependencies. Create and activate a virtual environment:

```bash=
cd $CTRCA_SOFT
python3 -m venv venv_ctRCA
source venv_ctRCA/bin/activate
```
Keep this environment active during installation of Python-based tools.

---

### Installing SAMtools (v1.20)
SAMTools is pre-installed on the OTAGO OnDemand Server, and often pre-installed on University-supplied computers. Check if it is already installed with:
```bash=
samtools --version
```

If not installed:
```bash=
cd $CTRCA_SOFT
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar -xjf samtools-1.20.tar.bz2
cd samtools-1.20

./configure --prefix=$CTRCA_SOFT/samtools-1.20
make
make install
```
Add to PATH:
```bash=
export PATH=$CTRCA_SOFT/samtools-1.20/bin:$PATH
```

---

### Installing Dorado (v1.4.0)

We are using Dorado v1.4.0 (current as of 20-Apr-2026).

Dorado can be installed by downloading a pre-built binary. Ensure the correct version for your operating system is selected. Download links are available on the [Dorado documentation page](https://software-docs.nanoporetech.com/dorado/latest/). 

To install from the command line:
```bash=
cd $CTRCA_SOFT
curl -L https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz -o dorado.tar.gz
tar -xzf dorado.tar.gz
mv dorado-1.4.0-linux-x64 dorado-1.4.0
```

Add to PATH:
```bash=
export PATH=$CTRCA_SOFT/dorado-1.4.0/bin:$PATH
```
Test installation:
```bash=
dorado --help
dorado --version
```

---

### Installing BEDTools (v2.29.1)

We are using BEDTools v.2.28.0. (Current as of 20-Apr-2026). All usage information can be found in the [GitHub repository](https://github.com/arq5x/bedtools2) or [documentation](https://bedtools.readthedocs.io/en/latest/).

```bash=
cd $CTRCA_SOFT
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -xzf bedtools-2.29.1.tar.gz
cd bedtools2
make
```
Add to PATH:
```bash=
export PATH=$CTRCA_SOFT/bedtools2/bin:$PATH
```

Test installation:
```bash=
bedtools --help
```

---

### Installing BB-8 (v1.0)
We are using BB-8 v1.0 (current as of 20-Apr-2026). BB-8 is a Consensus caller for R2C2 concatemeric reads - adapted directly from C3POa by R-Volden and C-Vollmers by the Guilford lab for ctDNA analyses. 

Ensure the virtual environment is active:
```bash=
source $CTRCA_SOFT/venv_ctRCA/bin/activate
```

Then clone the GitHub repository:
```bash=
cd $CTRCA_SOFT
git clone https://github.com/M-Dunnet/BB-8.git
cd BB-8

chmod +x setup.sh
./setup.sh
```
This installs required Python modules and dependencies (including `conk` and `blat`).

---

### Installing Cutadapt
We are using Cutadapt v.5.0. (Current as of 20-Apr-2026).

Install cutadapt into the virtual environment:
```bash=
source $CTRCA_SOFT/venv_ctRCA/bin/activate

pip install --upgrade pip
pip install cutadapt
```
Test installation
```bash=
cutadapt --version
```

---

### Installing minimap2

We are using Minimap2-2.30. (Current as of 20-Apr-2026). Documentation and build instructions can be found [here](https://lh3.github.io/minimap2/minimap2.html) and [here](https://github.com/lh3/minimap2/releases). Samtools is required for minimap2.

To Download:
```bash=
cd $CTRCA_SOFT
git clone https://github.com/lh3/minimap2.git
cd minimap2
make
```
Add to PATH:
```bash=
export PATH=$CTRCA_SOFT/minimap2:$PATH
```
**Indexing a Genome**
Indexing the reference genome is required before mapping reads. This step creates a minimap2 index file (`.mmi`), which allows for fast and efficient alignment.
```bash=
minimap2 -d <genome_index.mmi> <genome.fasta>
```

#### Custom Mapping Script
I have Written a script to map, filter, and index multiple files: `MapFiles.sh`. This is found in the preprocessing sub-directory of ctRCA-Var. 

---

### Installing ctRCA dependancies

Ensure the virtual environment is active:
```bash=
source $CTRCA_SOFT/venv_ctRCA/bin/activate
```

Clone the GitHub repository:
```bash=
cd $CTRCA_SOFT
git clone https://github.com/M-Dunnet/ctRCA-Var.git
```

---

### Final Notes

To make the software easily available in all future terminal sessions, add the PATH exports to your shell startup file. On many HPC systems, this is `~/.profile`.

Open the `.profile` file in your home directory:

```bash=
nano ~/.profile
```

This will either create a new file or open an existing one. Add the following lines:
```bash=
export CTRCA_SOFT=~/software/ctRCA
export PATH=$CTRCA_SOFT/samtools-1.20/bin:$CTRCA_SOFT/dorado-1.4.0/bin:$CTRCA_SOFT/bedtools2/bin:$CTRCA_SOFT/minimap2:$PATH
```
Save the file with `CTRL + O`, press `Enter`, then exit with `CTRL + X`. These settings will be applied automatically each time you open a new terminal.

To apply the changes immediately in your current session, run:
```bash
source ~/.profile
```
Test this has worked correctly:
```bash=
echo $PATH
```
You should see a long string of paths. It should begin with something similar to the following (displayed as a single line in the terminal):
```bash=
/home/{yourprofile}/software/ctRCA/samtools-1.20/bin:
/home/{yourprofile}/software/ctRCA/dorado-1.4.0/bin:
/home/{yourprofile}/software/ctRCA/bedtools2/bin:
/home/{yourprofile}/software/ctRCA/minimap2:
```
