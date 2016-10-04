#=================================================================================================
# Step 0: clean spades
#=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 preclean -f ~/ARMS/data/20K_R1.fq \
-r ~/ARMS/data/20K_R2.fq  -o 0_preclean -j 1

#=================================================================================================
# Step 1: assembling and renaming sequences
#=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 assemble -n BALI -f 0_preclean  -r 0_preclean  -o 1_assembled

#=================================================================================================
# Step 2: demux the file into independent samples and rename the sequences with sample names and serial ID#s
#=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 demux_samples -i 1_assembled -b ~/ARMS/data/barcodes.txt -o 2a_split
nice -10 python ~/ARMS/src/ARMS/chewbacca.py rename -i 2a_split -o 2b_renamed -f fastq

#=================================================================================================
# Step 3: Trim the (a)dapter and the (b)arcode sequences from thedata.
#======================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 trim_adapters  -i 2b_renamed -o 3_trim_ab -a ~/ARMS/data/adapters.fasta \
-arc ~/ARMS/data/adapters_RC.fa

#=================================================================================================
# Step 4: Clean the low quality of the reads using a sliding window
#=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 clean_seqs -i 3_trim_ab -o 4a_cleaned -m 200 -w 5 -q 25
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 convert_fastq_to_fasta -i 4a_cleaned -o 4b_fasta

#=================================================================================================
# Step 5: Dereplicate the data set
#=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 ungap_fasta -i 4b_fasta -o 4c_ungap -f fasta -g .-
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 merge_files -i 4c_ungap -o 4d_merged -f fasta -n BALI
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 dereplicate_fasta -i 4d_merged -o 5_derep

##=================================================================================================
## Step 6: Split the file to prepare to run MACSEx
##=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 partition -i 5_derep -o 6_partitioned -c 200 -f fasta

##=================================================================================================
## Step 7: run MACSE
##=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 macseAlign -i 6_partitioned -o 7_macseAligned -d ~/ARMS/data/BIOCODETEMPLATE 

##=====================================================================================================
## Step 8: Filter the data to be clustered
##+========================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 macseClean -i 7_macseAligned -o 8_macseCleaned -d ~/ARMS/data/BIOCODETEMPLATE -s 6_partitioned

##=================================================================================================
## Step 9: Cluster the sequences
##=================================================================================================
# SWARM
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 cluster_seqs -i 8_macseCleaned -o 10_clustered -g 5_derep_groups_files

# CROP
# nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 cluster_seqs -i 8_macseCleaned -o 10_clustered -g 5_derep_groups_files -p crop -z 50 -b 173

# Vsearch
# nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 cluster_seqs -i 8_macseCleaned -o 10_clustered -g 5_derep_groups_files -p vsearch


##=================================================================================================
## Step 10: Build the abundance matrix
##=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py build_matrix -b ~/ARMS/data/barcodes.txt -g 10_clustered_groups_files -s 2b_renamed_samples -o 11_buildmatrix


##=================================================================================================
## Step 11: Identify hits against biocode/Bold/ncbi
##=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 query_fasta  -j 4 -i 10_clustered -o 12_biocode -r ~/ARMS/data/BiocodePASSED_SAP.txt -x ~/ARMS/data/BiocodePASSED_SAP_tax_info_formatted.txt -s 97 -c 85
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 query_fasta  -j 4 -i 10_clustered  -o 13_bold -r ~/ARMS/data/bold100k.fna -x ~/ARMS/data/seq_lin.mapping -s 97 -c 85
nice -10 python ~/ARMS/src/ARMS/chewbacca.py -t 2 query_db -i 10_clustered -o 14_ncbi -r ~/ARMS/refs/COI.fasta -d ~/ARMS/refs/ncbi.db



##=================================================================================================
## Step 12: Annotate the abundance matrix with OTU identities
##=================================================================================================
nice -10 python ~/ARMS/src/ARMS/chewbacca.py annotate_matrix -i 11_buildmatrix -a 12_biocode -o 15_annotate_biocode
nice -10 python ~/ARMS/src/ARMS/chewbacca.py annotate_matrix -i 15_annotate_biocode -a 13_bold -o 16_annotate_bold
nice -10 python ~/ARMS/src/ARMS/chewbacca.py annotate_matrix -i 16_annotate_bold -a 14_ncbi -o 17_annotate_ncbi



