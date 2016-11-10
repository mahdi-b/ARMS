from Bio import SeqIO
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from classes.BufferedWriter import BufferedSeqWriter
from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import getInputFiles, debugPrintInputInfo, init_pool, run_parallel, printVerbose, strip_ixes, \
    makeDirOrdie, bulk_move_to_dir, cleanup_pool, makeAuxDir
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from classes.PythonRunner import PythonRunner
from Cluster_Helpers import handle_groups_file_update
from parse.parseUCtoGroups import parseUCtoGroups
from rename.renameWithoutCount import removeCountsFromGroupsFile


class Cluster_Program_Swarm(ChewbaccaProgram):
    name = "swarm"

    def execute_program(self):
        args = self.args
        self.cluster_swarm(args.input_f, args.outdir, args.groupsfile, args.processes, args.extraargstring)

    # TODO IMPORTANT: make sure that the abundances in dereplicated_renamed.fasta are sorted in decreasing order
    # We need to explore this more. Run by default for now....
    # Nore that any program can be used here as long as two files
    # 1- seeds: contains the cluster centroids. This file contains the updates counts for each cluster.
    # ex. a seq 97_2 from the cluster, if selected as a seed, would not be for example 97_100. This indicates
    # that 98 sequences are now assigned to the cluster for which 97 is a seed.
    # 2-clustering.out: contains the clustering results. (see file for sample format)

    # ~/bin/swarm/src/swarm dereplicated_renamed.fasta \

    def cluster_swarm(self, input_f, outdir, groupsfile, processes, extraargstring):
        """Clusters sequences using SWARM.
        :param input_f: A file or folder containing fasta files to cluster.
        :param outdir: The output directory results will be written to.
        :param groupsfile: A groups file or folder containing groups files that describe the input. Note: if no groups
                            file is supplied, then entries in the fasta file are assumed to be singleton sequences.
        :param processes: The maximum number of processes to use.
        :param extraargstring: Advanced program parameter string.
        """
        # Grab the fasta file(s) to cluster
        inputs = getInputFiles(input_f)
        debugPrintInputInfo(inputs, "clustered")
        pool = init_pool(min(len(inputs), processes))

        # RUN CLUSTERING
        run_parallel([ProgramRunner(ProgramRunnerCommands.CLUSTER_SWARM,
                                    [input_, "%s/%s_clustered" % (outdir, strip_ixes(input_)),
                                     "%s/%s_clustered_uc" % (outdir, strip_ixes(input_)),
                                     "%s/%s_clustered_seeds" % (outdir, strip_ixes(input_))],
                                    {"exists": [input_]}, extraargstring)
                      for input_ in inputs], pool)

        # PARSE UC FILE TO GROUPS FILE
        printVerbose("Parsing the clustered uc files to groups files")
        clustered_uc_files = getInputFiles(outdir, "*_clustered_uc")
        debugPrintInputInfo(clustered_uc_files, "parsed to groups")
        run_parallel([PythonRunner(parseUCtoGroups, [input_, "%s/%s.groups" % (outdir, strip_ixes(input_))],
                                   {"exists": [input_]})
                      for input_ in clustered_uc_files], pool)
        printVerbose("Done parsing groups files.")

        # REMOVE COUNTS FROM CLUSTERING GROUPS FILE
        printVerbose("Cleaning the .groups file from clustering")
        # Grab the current groups file and the new clustered groups file (which needs to be cleaned)
        clustered_groups_files = getInputFiles(outdir, "*_clustered.groups")
        debugPrintInputInfo(clustered_groups_files, "cleaned")
        run_parallel([PythonRunner(removeCountsFromGroupsFile,
                                   [input_, "%s/%s_uncount.groups" % (outdir, strip_ixes(input_))],
                                   {"exists": [input_]})
                      for input_ in clustered_groups_files], pool)
        printVerbose("Done cleaning groups files.")

        printVerbose("Capitalizing sequences")
        # Convert the seeds files to uppercase (swarm writes in lowercase)
        inputs = getInputFiles(outdir, "*_seeds")
        run_parallel([PythonRunner(capitalize_seqs, [input_, "%s.fasta" % input_], {"exists": [input_]})
                      for input_ in inputs], pool)
        printVerbose("Done capitalizing sequences")

        # Collect the groups file from clustering with counts removed
        cleaned_clustered_groups_files = getInputFiles(outdir, "*_uncount.groups", ignore_empty_files=False)

        # Resolve the user specified names file if necessary
        final_groups_files = handle_groups_file_update(outdir, groupsfile, cleaned_clustered_groups_files)

        # Move the final groups file(s) to the groups dir
        groups_dir = makeDirOrdie("%s_groups_files" % outdir)
        bulk_move_to_dir(final_groups_files, groups_dir)


        # Move aux files to the aux dir
        aux_files = getInputFiles(outdir, "*", "*_seeds.fasta", ignore_empty_files=False)
        aux_dir = makeAuxDir(outdir)
        bulk_move_to_dir(aux_files, aux_dir)

        # Cleanup the pool
        cleanup_pool(pool)


def capitalize_seqs(input_fasta, output_fasta, filetype='fasta'):
    """Capitalizes the ATGC sequence in a fasta file and writes it to a new file.

    :param input_fasta: Filepath to the input fasta file to capitalize.
    :param output_fasta: Filepath to the output fasta file.
    :param filetype: The file format to read and write.  Either 'fasta' or 'fastq'
    :return: Filepath to the output fasta file.
    """
    capitalized_output_file = BufferedSeqWriter(output_fasta, filetype)

    for sequence in SeqIO.parse(open(input_fasta, 'rU'), "fasta"):
        sequence.seq = Seq(str(sequence.seq).upper(), SingleLetterAlphabet())
        capitalized_output_file.write(sequence)
    capitalized_output_file.flush()
    return output_fasta
