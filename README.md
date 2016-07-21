# MACSE_VR smaller verified db
# BIOCODETMEPLATE Complete database

python arms.py preprocess --name index2 -f data/Index2_R1.fastq \
       -r data/Index2_R2.fastq -b data/Adapter2.txt \
        --db data/BIOCODETEMPLATE --outDir outDir

This is partition?

python arms.py -t 10 clean -n Index2 -s outDir/samples/ -p macse -o outDir/clenaedSamples/ \
        --db data/BIOCODE_MACSE_VR.fasta

python arms.py removeChimeras -n Index2 -i outDir/clenaedSamples/MACSEOUT_MERGED.fasta -f outDir/clenaedSamples/Current_updated.names

python arms.py dropShort -n Index2 -i outDir/clenaedSamples/MACSEOUT_MERGED.fasta -f outDir/clenaedSamples/Current_updated.names  -l 90 -o SOMEOUT.fasta -s SOMECLEAN.names



