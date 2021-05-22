This is a support repo for "Assembling viruses and identifying their hosts from long-read metagenomics data"

Requirements
----------
* python 3.*
* hmmsearch
* prodigal
* samtools
* seqtk
* minced
* minimap2
* [Flye](https://github.com/fenderglass/Flye) genome assembler 
* [viralVerify](https://github.com/ablab/viralVerify)
* [viralComplete](https://github.com/ablab/viralComplete)


All mentioned tools should be in your $PATH.


Output
---------
Information about circular and linear complete viruses can be found in `vc_circulars` and `vc_linears` subfolders in metaflye output folder.
Viral components can be found in `vv_components` subfolder. All viruses (including incomplete) can be found in `vv_circulars` and `vv_linears` subfolders. Coordinates of falsely circularized viruses detected by CircularDisconnector are in `<metaflye output folder>/vv_circulars/Prediction_results_fasta/linear_check/linears.txt`. Details can be found in viralVerify and viralComplete user manuals.

Prediction of hosts within the sample is performed by a separate script `crispr_host_match.py`. It takes metaFlye result as an input, extracts circular and linear isolated contigs,predicts viruses and CRISPR spacers and matches them using BLAST. Result (BLAST output format 6) can be found in `blast.out` file in the output folder.
