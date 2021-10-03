viralFlye
=========

viralFlye is a pipeline to recover high-quality viral genomes from long-read metagenomic sequencing.

Installation
------------

You need the `conda` package manager for the installation process.
To install viralFlye locally, run: 

```
git clone https://github.com/Dmitry-Antipov/viralFlye
cd viralFlye
install.sh
```

This will create a conda environment viralFlye which contain all dependencies.

You also need the Pham HMM database for viral genome identification. If you don't have it yet, download using:

```
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz
```

Running viralFlye
-----------------

To run, activate the conda environemnt first, then invoke the pipeline:

```
conda activate viralFlye
./viralFlye.py
```

viraFlye takes an existing metaFlye assembly directory as input. Use [metaFlye v2.9+](https://github.com/fenderglass/Flye) with `--meta` option.
You will also need the original reads that were used to assemble.

Then, pipeline could be invoked as follows:

```
./viralFlye.py --dir flye_assembly_dir --hmm path_to_Pfam-A.hmm.gz --reads path_to_reads --outdir output_dir
```

Additional options are shown when running viralFlye script without parameters.

Output
------

viralFlye output can be found in the outdir directory (by default, equal to the input directory)
It consists of 3 fasta files, `linears_viralFlye.fasta`, `components_viralFlye.fasta` and `circulars_viralFlye.fasta`,
and a txt file that lists all erroneously circularized components.


Prediction of hosts within the sample is performed by a separate script `crispr_host_match.py`. 
It takes metaFlye result as an input, extracts circular and linear isolated contigs,predicts viruses and CRISPR spacers and matches them using BLAST. 
Result (BLAST output format 6) can be found in `blast.out` file in the output folder.

Dependencies
-----------

viralFlye package depends on the following software

* [viralVerify](https://github.com/ablab/viralVerify)
* [viralComplete](https://github.com/ablab/viralComplete)
* prodigal 
* samtools 
* seqtk 
* minced 
* minimap2 
* biopython 
* pysam 
* samtools
* freebayes
* bcftools
* numpy
* scipy
* BLAST


License
-------

viralFlye is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Code contributors:

* Dmitry Antipov (Center for AlgorithmicBiotechnology, Saint PetersburgState University)
* Mikhail Rayko  (Center for AlgorithmicBiotechnology, Saint PetersburgState University)
* Mikhail Kolmogorov (University of California, Santa Cruz)
