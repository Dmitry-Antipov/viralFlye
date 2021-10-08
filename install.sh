conda create -n viralFlye -c bioconda -c conda-forge -c mikeraiko "python>=3.6" prodigal viralverify vcflib seqtk minced minimap2 biopython pysam tabix samtools freebayes bcftools numpy scipy blast bwa viralcomplete

#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh
#git clone  https://github.com/Dmitry-Antipov/viralFlye
#export PATH=$PATH:$(pwd)/viralFlye
#export PATH=$PATH:$(pwd)
#export PATH=$PATH:$(pwd)/viralComplete
#chmod +x viralFlye.py
