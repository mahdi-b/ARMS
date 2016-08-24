#!/usr/bin/env bash
python ~/ARMS/src/ARMS/dev/makeBadSeqs.py ~/ARMS/data/bold_sample.fasta 1000  bads.fasta
python ~/ARMS/src/ARMS/chewbacca.py minhash -i bads.fasta -o . -d ~/ARMS/data/bold_sample.fasta -m 2g -s 20
python ~/ARMS/src/ARMS/dev/alignReadsProt.py ~/ARMS/data/BIOCODE_071216_MACSE_RENAMED_CORRECTED_FILTERD.fa bads.fasta bads_20.out > rslt.out
