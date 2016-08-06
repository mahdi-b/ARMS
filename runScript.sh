#=================================================================================================
# Step 1: assembling and renaming sequences
# time	2.0s
# req	PEAR
#=================================================================================================


# assemble the forward and reverse reads. -o BALI is the prefix that will be assigned to the newly created files
#mkdir 1_assembled && cd 1_assembled
#pear -f ../testData/20K_R1.fq -r ../testData/20K_R2.fq -o BALI -j 10

'python chewbacca.py assemble -p pear -n BALI -f /home/greg/ARMS/testARMS/testData/20K_forward.fq  \
    -r /home/greg/ARMS/testARMS/testData/20K_reverse.fq  -o 1_assembled -t 1 -m 550
'
#=================================================================================================
# Step 2: demux the file into independent samples
# time	0.247s
#=================================================================================================
#cd ../
#mkdir 2_split && cd 2_split

# splits the assembled file based on the samples defined in the barcodes file.
# see testData/barcodes.txt for an example
# Note: this creates an unmatches file... those are the sequences that do not match any barcode.
# In what follows, I ignore the unmatched sequences
#cat ../1_assembled/BALI.assembled_renamed.fastq | fastx_barcode_splitter.pl \
#    --bcfile ../testData/barcodes.txt --prefix out_ \
#    	     --suffix .fastq --bol --mismatches 1

python chewbacca.py demux_samples -i 1_assembled -b ~/ARMS/testARMS/testData/barcodes.txt -o 2_split
python chewbacca.py rename -i 2_split -o 3_renamed -f fastq

#=================================================================================================
# Step 3: Trim the (a)dapter and the (b)arcode sequences from thedata.
# time:	2.0s + 1.9s
# req	flexbar
#=================================================================================================
#cd ../
#mkdir 3_trim_ab && cd 3_trim_ab
# First remove the barcodes from the left site,
#ls ../2_split/*.fastq | grep -v unmatched | parallel  flexbar \
#     -r {} -t {/.}_temp_out \
#     -ae LEFT -a ../testData/adapters.fasta

# then remove the adapter from the right side
##ls ../2_split/*.fastq | grep -v unmatched | parallel flexbar \
#     -r {/.}_temp_out.fastq -t {/.}_debarcoded \
#     -ae RIGHT -a ../testData/adapters_RC.fa

python chewbacca.py trim_adapters  -p flexbar -i 3_renamed -o 4_trim_ab -b ~/ARMS/testARMS/testData/adapters.fasta -a ~/ARMS/testARMS/testData/adapters_RC.fa

#=================================================================================================
# Step 4: clean the low quality of the reads using a sliding window
# Params for this step need to be explored more
# time	.836s
# req	default-jre
#=================================================================================================
#cd ../
#mkdir 4_cleaned && cd 4_cleaned

#ls ../3_trim_ab/*debarcoded.fastq | parallel "java -jar ~/ARMS/programs/Trimmomatic-0.33/trimmomatic-0.33.jar \
#   SE -phred33 {} {/.}_cleaned.fastq \
#      SLIDINGWINDOW:5:25 MINLEN:200"
python chewbacca.py clean_seqs -i 4_trim_ab -o 5_cleaned

#=================================================================================================
# Step 5: dereplicate the data set
# time	0.262s + 0.266s + 0.141s + 0.417s = 1.086s
# req	emboss NOTE: Usearch is not open licensed
#=================================================================================================
#cd ../
#mkdir 5_dereplicate && cd 5_dereplicate

# dereplication requires that we first convert the dataset to fasta
#ls ../4_cleaned/out_*fastq  | parallel  seqret {} {/.}.fasta
#ls *cleaned.fasta  | parallel  "~/ARMS/bin/usearch7.0.1090 -derep_fulllength {} -output {/.}_derep.fa -uc {/.}_uc.out"

# parse the number of identical reads included in each sequence and write them to the {sample_file}_parsed.out
#ls *uc.out | parallel "python ~/ARMS/bin/getSeedSequences.py {} {.}_parsed.out"

# rename the sequences to include the the number of identical reads. Ex. in the fasta file {sample_file}_derep_renamed.fa,
# read 123_10 indicate that for sequence which id is 123, there 10 sequences that identical to it and
# which were discarded.
#ls *cleaned.fasta | parallel "python ~/ARMS/bin/renameSequences.py {/.}_derep.fa {/.}_uc_parsed.out {/.}_derep_renamed.fa"
python chewbacca.py dereplicate_fasta -i 5_cleaned/ -o 6_derep

#
##=================================================================================================
## Step 6: Align the reads against the reference to infer orientation
## this will be most likely change in the future since we don't need to align to
## infer orientation (can be inferred from barcode)
## time	7.82s
##=================================================================================================
## run mother with the flip param
## make it from the same directory as align.seqs can cause errors if run with relative path.
#ls *debarcoded_cleaned_derep_renamed.fa | parallel  \
#   mothur '"#align.seqs(candidate={}, template=~/ARMS/data/BIOCODETEMPLATE, flip=t)"'
#
## create a new directory and move results to it after aligning the reads
#cd ..
#mkdir 6_aligned && cd 6_aligned
#mv ../5_dereplicate/*align* .
#mv ../5_dereplicate/*.flip.accnos .
#
##=================================================================================================
## Step 7: Split the file to prepare to run MACSE
## NOTE: We split all the files into smaller chunks so that we can run them on
## multiple core in parallle.
## IMPORTANT, once you split files, you do not want to run ls on that directory, since it
## could contain a large number of files. This would substancially slow you down
## time	0.365s
##=================================================================================================
#cd ../
#mkdir 7_split_aligned && cd 7_split_aligned
#
## IMPORTNAT, make sure you adjust the number of sequences of files to obtain a relatively small number
## of chunks. Ex: no more than 300 files.
## here, for the sake of simplicity, we are splitting in files of 1000 sequences each.
## note, the last parameter is the format
#ls ../6_aligned/*renamed.align | parallel "python ~/ARMS/bin/splitKperFasta.py {} {/.} 1000 fasta"
#
#
##=================================================================================================
## Step 8: run MACSE on all the chunks we created
## -t 20 indicates that we will be running on 20 cores. Adjust this to fit the number
## of cores available on the machine.
##=================================================================================================
#cd ..
#mkdir 8_macse_out && cd 8_macse_out
#
## Chewbacca runs enrichAlignment, then exportAlignment. In the last stage, it cleans the alignments and concatenates them
## into a single file MACSEOUT_MERGED.fasta, which can be clustered
#python ~/ARMS/SRC/ARMS/chewbacca.py -t 20 clean -n BALI -s ../7_split_aligned \
#       -p macse -o . --db ~/ARMS/data/BIOCODE_MACSE_VR.fasta
#
#
#
#
#
#
##=================================================================================================
## Step 9: Cluster the sequences
##=================================================================================================
#cd ../
#mkdir 9_cluster && cd 9_cluster
#
## prepare the file for clustering. Dereplicates across smaples and renames the resulting sequences with hashes
#~/bin/vsearch/bin/vsearch --derep_fulllength \
#  ../8_macse_out/MACSEOUT_MERGED.fasta --sizeout \
#   --fasta_width 0 --output amplicons_linearized_dereplicated.fasta -uc uc.out
#
# python ~/bin/getSeedSequences.py uc.out uc_parsed.out
#
# python ~/bin/formatReadWithCounts.py ../8_macse_out/MACSEOUT_MERGED.fasta uc_parsed.out dereplicated_renamed.fasta
#
#
#
## IMPORTANT: make sure that the abundances in dereplicated_renamed.fasta are sorted in decreasing order
#
#
## We need to explore this more. Run by default for now....
## Nore that any program can be used here as long as two files
## 1- seeds: contains the cluster centroids. This file contains the updates counts for each cluster.
## ex. a seq 97_2 from the cluster, if selected as a seed, would not be for example 97_100. This indicates
## that 98 sequences are now assigned to the cluster for which 97 is a seed.
## 2-clustering.out: contains the clustering results. (see file for sample format)
#~/bin/swarm/src/swarm dereplicated_renamed.fasta \
#  		      				       --output-file clustering.out -u uclust_file -w seeds
#
## We convert th seeds files to uppercase (strangely, it is output in lowercase by swarm)
## this might not be necessary with other clustering programs
#
#tr  '[:lower:]' '[:upper:]' < seeds > seeds.fasta
#rm seeds
#
#
##=================================================================================================
## Step 9_P: Find Chimeras
## we run this step after the clustering in order to speed up search since we align the sequences
## against themselves.
##=================================================================================================
#cd ..
#mkdir 9_p_uchime && cd 9_p_uchime
#cp ../9_cluster/seeds.fasta .
#grep ">"  seeds.fasta | sed 's/>//' | awk '{print $1"\t"$1}' > seeds.names
#mothur "#chimera.uchime(fasta=seeds.fasta, name=seeds.names)"
#python /tmp/mothur_test/removeChimeras.py seeds.fasta seeds.uchime.accnos
#
##=================================================================================================
## Step 10: Identify hits against biocode
##=================================================================================================
#cd ../
#mkdir 10_search_biocode && cd 10_search_biocode
#
#~/bin/vsearch/bin/vsearch-1.1.1-linux-x86_64 --usearch_global \
#  					     ../9_p_uchime/seeds.pick.fasta  --db ../data/BiocodePASSED_SAP.txt --id 0.9 \
#					     				     --userfields query+target+id+alnlen+qcov --userout out \
#									     		  --alnout alnout.txt
#
## Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in parsed_BIOCODE.out
## parameters can be changed and this command can be rerun as many times as necessary
#python ~/bin/parseVSearchout.py ../data/BiocodePASSED_SAP_tax_info.txt  out 97 85
#       > parsed_BIOCODE.out
#
##=================================================================================================
## Step 11: Run the alignment against the sequences extracted from the NCBI's NT.
## extracting the sequences and building the necessary database needs to be done once.
## I make the DB available on the server
##=================================================================================================
#cd ../
#mkdir 11_search_nt && cd 11_search_nt
#
#~/bin/vsearch/bin/vsearch-1.1.1-linux-x86_64 --usearch_global \
#  ../9_p_uchime/seeds.pick.fasta  --db /home/mahdi/refs/COI_DOWNLOADED/COI.fasta \
#      --id 0.9 \
#          --userfields query+target+id+alnlen+qcov --userout out \
#	      --alnout alnout.txt --userfields query+target+id+alnlen+qcov
#
## This needs the taxonomy to be available. See below...
## Parse the alignment results and put those that pass the criterion (97 similarity, 85 coverage) in parsed_nt.out
## parameters can be changed and this command can be rerun as many times as necessary
#python ~/bin/getTaxFromId.py out /home/mahdi/refs/taxonomy/ncbi.db 97 85 \
#       > parsed_nt.out
#
##=================================================================================================
## Step 12: Build the abundance matrix
##=================================================================================================
#
#cd ../
#mkdir 12_matrix && cd 12_matrix
#
#python ~/bin/builMatrix.py ../5_dereplicate/ ../9_cluster/uc_parsed.out ../9_cluster/clustering.out matrix.out
#
#
## FOR BLAST
## DO THIS ONCE TO EXTRACT SEQUENCES AND Format them
## Uncomment the code below to run only once. Not that this
## will take a very long time!!!!
#
## blastdbcmd -db ~/refs/nt/nt -entry all -outfmt "%g;;%t" > COI.ids
## awk -F ";;" '{print $1}' COI.ids | blastdbcmd -db ~/refs/nt/nt -dbtype nucl  -entry_batch - -out COI.fasta
## perl -ne 'if ($_ =~ />/){@data= split(/\|/, $_); print ">".$data[1]."\n"}else{print $_}' COI.fasta  > COI.fasta_
## mv COI.fasta_ COI.fasta
## sqlite3 ncbi.db
## create table gi_taxid(gi integer, taxid integer);
## .tables
## .mode list
## .separator \t
## .import gi_taxid_nucl.dmp gi_taxid
## CREATE UNIQUE INDEX gi_idx_on_gi_taxid ON gi_taxid(gi);
