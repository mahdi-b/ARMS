#!/usr/bin/env bash
# get prot seqs
# grep '>' ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.faa  |cut -f 1 | sed -s 's/>//' | shuf| head -600000 > nuc_names.txt
# pull them from the file
#seqtk subseq ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.faa nuc_names.txt > bold_sample.faa

# 100k's are DBS
# 10k's go into make bads

# Have BOLD.faa
# from BOLD.faa -> sample bold_110k.faa , bold110k.faa.names
    # get prot seqs
grep '>' ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.faa  |cut -f 1 | sed -s 's/>//' | shuf| head -110000 > ~/ARMS/data/bold110k.faa.names
    # pull them from the file
seqtk subseq ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.faa ~/ARMS/data/bold110k.faa.names > ~/ARMS/data/bold110k.faa
cat ~/ARMS/data/bold110k.faa.names |cut -f 1 -d '_' > ~/ARMS/data/bold110k.fna.names
paste ~/ARMS/data/bold110k.fna.names ~/ARMS/data/bold110k.faa.names > ~/ARMS/data/bold110k_name_pairs.txt
# from bold110k.faa -> sample 10k.faa, bold10k.faa.names
grep '>' ~/ARMS/data/bold110k.faa  |cut -f 1 | sed -s 's/>//' | shuf| head -10000 > ~/ARMS/data/bold10k.faa.names
seqtk subseq ~/ARMS/data/bold110k.faa ~/ARMS/data/bold10k.faa.names > ~/ARMS/data/bold10k.faa

# diffFasta bold10k from bold110k -> bold100k, bold 100k.faa.names
python ~/ARMS/src/ARMS/utils/diffFasta.py ~/ARMS/data/bold110k.faa ~/ARMS/data/bold10k.faa ~/ARMS/data/bold100k.faa
grep '>' ~/ARMS/data/bold100k.faa  |cut -f 1 | sed -s 's/>//' | shuf > ~/ARMS/data/bold100k.faa.names

# Convert faa names to fna names
cat ~/ARMS/data/bold100k.faa.names|cut -f 1 -d '_' > ~/ARMS/data/bold100k.fna.names
cat ~/ARMS/data/bold10k.faa.names|cut -f 1 -d '_' > ~/ARMS/data/bold10k.fna.names

# Have BOLD.fna, bold 100k.faa.names, bold10k.faa.names
# from BOLD.fna, bold100k.faa.names -> bold100k.fna, bold100k.fna.names
seqtk subseq ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.fna ~/ARMS/data/bold10k.fna.names > ~/ARMS/data/bold10k.fna
# from BOLD.fna, bold10k.faa.names -> bold10k.fna, bold10k.fna.names
seqtk subseq ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.fna ~/ARMS/data/bold100k.fna.names > ~/ARMS/data/bold100k.fna

paste ~/ARMS/data/bold100k.fna.names ~/ARMS/data/bold100k.faa.names > ~/ARMS/data/bold100k_name_pairs.txt

nice -10 python ~/ARMS/src/ARMS/dev/makeBadSeqs.py ~/ARMS/data/bold10k.fna 1000 bads.fasta
nice -10 python ~/ARMS/src/ARMS/chewbacca.py minhash -i bads.fasta -o . -d ~/ARMS/data/bold100k.fna -m 32g -s 40,20,5
nice -10 python ~/ARMS/src/ARMS/dev/alignReadsProt.py ~/ARMS/data/bold100k.fna bads.fasta bads_40.out  ~/ARMS/data/bold110k_name_pairs.txt ~/ARMS/data/bold100k.faa 1 > h_rslt.40.out

#grep Error.*_i.* -B 2 h_rslt.out > h_errors.txt

#nice -10 python ~/ARMS/src/ARMS/dev/indelPositionBadnessAnalysis.py ~/ARMS/testARMS/dev h_bads_40.out.tab
