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
nice -10 grep '>' ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.faa  |cut -f 1 | sed -s 's/>//' | shuf| head -10000 > ~/ARMS/data/bold10k.faa.names
    # pull them from the file
nice -10 seqtk subseq ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.faa ~/ARMS/data/bold10k.faa.names > ~/ARMS/data/bold10k.faa
nice -10 cat ~/ARMS/data/bold10k.faa.names |cut -f 1 -d '_' > ~/ARMS/data/bold10k.fna.names
nice -10 seqtk subseq ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.fna ~/ARMS/data/bold10k.fna.names > ~/ARMS/data/bold10k.fna
nice -10 paste ~/ARMS/data/bold10k.fna.names ~/ARMS/data/bold10k.faa.names > ~/ARMS/data/bold10k_name_pairs.txt

# diffFasta bold10k from bold110k -> bold100k, bold 100k.faa.names
nice -10 python ~/ARMS/src/ARMS/utils/diffFasta.py ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.faa ~/ARMS/data/bold10k.faa ~/ARMS/data/bold_ref.faa
nice -10 python ~/ARMS/src/ARMS/utils/diffFasta.py ~/ARMS/data/bold_coi_11_05_2015_corrected_filtered.fna ~/ARMS/data/bold10k.fna ~/ARMS/data/bold_ref.fna

nice -10 python ~/ARMS/src/ARMS/dev/makeBadSeqs.py ~/ARMS/data/bold10k.fna 10000 bads.fasta
nice -10 python ~/ARMS/src/ARMS/chewbacca.py minhash -i bads.fasta -o ~/ARMS/testARMS/dev -d ~/ARMS/data/bold_ref.fna -m 32g -s 40,20,5
nice -10 python ~/ARMS/src/ARMS/dev/alignReadsProt.py ~/ARMS/data/bold_ref.fna bads.fasta bads_40.out  ~/ARMS/data/bold10k_name_pairs.txt ~/ARMS/data/bold_ref.faa 1 > h_rslt.40.out

#grep Error.*_i.* -B 2 h_rslt.out > h_errors.txt

#nice -10 python ~/ARMS/src/ARMS/dev/indelPositionBadnessAnalysis.py ~/ARMS/testARMS/dev h_bads_40.out.tab
