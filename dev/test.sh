#!/usr/bin/env bash
# get prot seqs
# grep '>' BIOCODE_071216_MACSE_RENAMED_CORRECTED_FILTERD |cut -f 1 | sed -s 's/>//' | shuf| head -100000 >names.txt
# seqtk subseq bold_coi_11_05_2015_corrected_filtered.faa names.txt > bold_sample.faa

# get corresponding nucleotides
# cat names.txt|sed 's/_/ /'|cut -f 1 -d ' '>nuc_names.txt
# seqtk subseq bold_coi_11_05_2015_corrected_filtered.fna nuc_names.txt > bold_sample.fna

# python ~/ARMS/src/ARMS/dev/makeBadSeqs.py ~/ARMS/data/biocode_sample.fa 1000  bads.fasta
python ~/ARMS/src/ARMS/chewbacca.py minhash -i bads.fasta -o . -d ~/ARMS/data/bold_sample.fna -m 2g -s 20
python ~/ARMS/src/ARMS/dev/alignReadsProt.py ~/ARMS/data/bold_sample.fa bads.fasta bads_20.out > rslt.out
