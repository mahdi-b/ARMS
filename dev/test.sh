#!/usr/bin/env bash
python ~/ARMS/src/ARMS/dev/makeBadSeqs.py ~/ARMS/data/BIOCODE_MACSE_VR.fasta 1000  bads.fasta
python ~/ARMS/src/ARMS/chewbacca.py minhash -i bads.fasta -o . -d ~/ARMS/data/BIOCODE_MACSE_VR.fasta -m 2g -s 20
python ~/ARMS/src/ARMS/dev/alignReadsProt.py ~/ARMS/data/BIOCODE_MACSE_VR.fasta bads.fasta bads_20.out > rslt.out
