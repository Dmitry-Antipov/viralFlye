conda create -n viralFlye -c bioconda -c conda-forge  "python>=3.6" prodigal viralverify samtools seqtk minced minimap2 biopython pysam tabix samtools freebayes bcftools numpy scipy blast

#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh
#git clone  https://github.com/Dmitry-Antipov/viralFlye
#export PATH=$PATH:$(pwd)/viralFlye
#export PATH=$PATH:$(pwd)
#export PATH=$PATH:$(pwd)/viralComplete
#chmod +x viralFlye.py

git clone  https://github.com/ablab/viralComplete
chmod +x viralComplete/viralcomplete.py
