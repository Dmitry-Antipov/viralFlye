This is a support repo for "Assembling viruses and identifying their hosts from long-read metagenomics data"

To install you should run ./install.sh , which installs all required packages via conda. Then enviroment viralFlye should be activated.
---
Input
To run viralFlye.py, you'll need a directory with metaflye(any version) output, set of HMMs used for viral verification (can be downloaded from Pfam-A, ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/) and original long reads.
Additional options are shown when running viralFlye script without parameters.

---
Output
viralFlye output can be found in the outdir directory (by default, equal to the input directory)
It consists of 3 fasta files, linears_viralFlye.fasta, components_viralFlye.fasta and circulars_viralFlye.fasta, and a txt file that lists all erroneously circularized components.


Prediction of hosts within the sample is performed by a separate script `crispr_host_match.py`. It takes metaFlye result as an input, extracts circular and linear isolated contigs,predicts viruses and CRISPR spacers and matches them using BLAST. Result (BLAST output format 6) can be found in `blast.out` file in the output folder.
