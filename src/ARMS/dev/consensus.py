from Bio.Seq import Seq

"""FROM http://biopython.org/DIST/docs/api/Bio.Align.AlignInfo-pysrc.html"""
def dumb_consensus(msa, threshold=.7, ambiguous="-", require_multiple=True):
    """Output a fast consensus sequence of the alignment.

 This doesn't do anything fancy at all. It will just go through the
 sequence residue by residue and count up the number of each type
 of residue (ie. A or G or T or C for DNA) in all sequences in the
 alignment. If the percentage of the most common residue type is
 greater then the passed threshold, then we will add that residue type,
 otherwise an ambiguous character will be added.

 This could be made a lot fancier (ie. to take a substitution matrix
 into account), but it just meant for a quick and dirty consensus.

 Arguments:
     - threshold - The threshold value that is required to add a particular
       atom.
     - ambiguous - The ambiguous character to be added when the threshold is
       not reached.
     - consensus_alpha - The alphabet to return for the consensus sequence.
       If this is None, then we will try to guess the alphabet.
     - require_multiple - If set as 1, this will require that more than
sequence be part of an alignment to put it in the consensus (ie.
       not just 1 sequence and gaps).
 """
    # Iddo Friedberg, 1-JUL-2004: changed ambiguous default to "X"
    consensus = ''

    # find the length of the consensus we are creating
    con_len = len(msa[0])

    # go through each seq item
    for n in range(con_len):
        # keep track of the counts of the different atoms we get
        atom_dict = {}
        num_atoms = 0

        for record in msa:
            # make sure we haven't run past the end of any sequences
            # if they are of different lengths
            if n < len(record.seq):
                if record.seq[n] != '-' and record.seq[n] != '.':
                    if record.seq[n] not in atom_dict:
                        atom_dict[record.seq[n]] = 1
                    else:
                        atom_dict[record.seq[n]] += 1

                    num_atoms = num_atoms + 1

        max_atoms = []
        max_size = 0

        for atom in atom_dict:
            if atom_dict[atom] > max_size:
                max_atoms = [atom]
                max_size = atom_dict[atom]
            elif atom_dict[atom] == max_size:
                max_atoms.append(atom)

        if require_multiple and num_atoms == 1:
            consensus += ambiguous
        elif (len(max_atoms) == 1) and ((float(max_size) /
                                             float(num_atoms)) >= threshold):
            consensus += max_atoms[0]
        else:
            consensus += ambiguous
    return consensus
