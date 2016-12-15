File Types
==========
Chewbacca uses several common file types. Those are described below.

Fasta Files
-----------
**common extensions: .fa, .fasta, .FASTA**
Read more here: https://en.wikipedia.org/wiki/FASTA_format


FastQ Files
-----------
**extensions: .fq, .fastq, .FASTQ**
Read more here: https://en.wikipedia.org/wiki/FASTQ_format

Groups Files
------------
**extensions: .groups**

Groups files are used by Chewbacca to keep track of groups, or clusters, of relatively similar sequences.
Group files are generated or updated after each dereplication or clustering step.
A Groups file consists of one or more lines in the following format:

  *GROUPID <tab> SequenceID (<space> SequenceID)* ...
 
As an example:

| Hawaii_site1_0_ID111  Hawaii_site1_0_ID111 Hawaii_site1_0_ID112 Hawaii_site1_0_ID113
| Indonedia_site1_0_ID115       Indonedia_site1_0_ID115 Indonedia_site1_0_ID117
| Philippines_site1_0_ID119     philippines_site1_0_ID119

*Notes:*

1. The GROUPID for a group/cluster is a representative sequence from that cluster.
        This means that a sequenceId  will likely appear twice on a line (once as a GROUPID, and once in the sequence SequenceIds list).


2. See the "naming conventions" section for more info on chewbacca sequence naming standards.
   SequenceId are created using the a combination of sameple name, file offset,and the sequential number
   ex. Hawaii_site1_0_ID119.
   - Hawaii_site1: This sequence is from the Rodent_gut sample.
   - 0 file offset. When more than one sequence file is used, the files are annotated using different offesets.
     This makes it easy to track which sequences came from which file, which could potentially represent different
     sequencing runs, or other things of interest. 

Samples Files
-------------
**extensions: .samples**

Samples files are used by Chewbacca to map sequence names to the the name of their respective sample names.
This file is generally written once, early on in the anylitical process, at the time of sequence renaming.
The primary purposes for writing this file are for annotation and construction of an OTU table at the end of the analysis.

A Samples file consists of one or more lines in the following format:

  *SequenceID <tab> SampleID*

As an example:

|  Hawaii_site1_0_ID111 GUT_SAMPLE_21
|  Hawaii_site1_0_ID112 GUT_SAMPLE_21
|  Rodent_gutID113 GUT_SAMPLE_22
|  Rodent_noseID115     NOSE_SAMPLE_1
|  Rodent_stomachID115  STOMACH_SAMPLE_2


Note that more than one sample file is generated when sequences form asample are present in more than one files. In such case,
Each file is assigned the a different offet. Ex.: Rodent_gut_0.samples, Rodent_gut_1.samples, etc...

Barcodes Files
--------------
**extensions: .barcodes, .txt**

User provided file containing the barcode information. The format is

** sampleID <tab> DNA_barcode_Sequence**

|  BALI_site_1          agacgc
|  BALI_site_2          agtgta
|  Hawaii_site_1        actagc
|  etc..

If the fasta files have already been tagged using barcode (i.e., sample name was appended to sequence id).
The sequence information can be any string of valid DNA. Ex. AAAA.

Adapters files
--------------
**extensions: .adapters, .txt**

fasta file containing sequencing barcodes.
This is used in the cleaning up stage to remove any occurrence of the barcode sequences in the fasta reads.
