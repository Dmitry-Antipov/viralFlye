conda create -n viralflye -c bioconda -c conda-forge  "python>=3.6" prodigal viralverify samtools seqtk minced minimap2 biopython pysam
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate viralflye
git clone  https://github.com/Dmitry-Antipov/viralFlye
git clone  https://github.com/ablab/viralComplete
export PATH=$PATH:$(pwd)/viralFlye
export PATH=$PATH:$(pwd)/viralComplete
chmod +x viralFlye/viralFlye.py
chmod +x viralComplete/viralcomplete.py
viralFlye/viralFlye.py