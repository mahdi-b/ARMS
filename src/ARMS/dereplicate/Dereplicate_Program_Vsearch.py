from classes.ChewbaccaProgram import ChewbaccaProgram
from classes.Helpers import getInputFiles, debugPrintInputInfo, init_pool, run_parallel, printVerbose, strip_ixes, \
    makeAuxDir, bulk_move_to_dir, cleanup_pool, makeDirOrdie
from classes.ProgramRunner import ProgramRunner, ProgramRunnerCommands
from classes.PythonRunner import PythonRunner
from parse.parseUCtoGroups import parseUCtoGroups
from rename.renameWithReplicantCounts import renameWithReplicantCounts
from rename.renameWithoutCount import removeCountsFromFastFile
from util.updateGroups import update_groups


class Dereplicate_Program_Vsearch(ChewbaccaProgram):
    """Dereplicates a fasta file by grouping identical reads together under one representative sequence.
        The number of replicant sequences each representative represents is given by a replicant count at the end of
        the sequence name in output fasta file.  If a .groups file is provided, then replicant counts will take into
        account previous dereplication counts (e.g. a replicant sequence that represents 3 sequences will add 3 to its
        representative sequence's replicant count).  Replicant counts are denoted with a suffix of '_X' on the sequence
        name, where X is the dereplication count.

        e.g. the sequence below is named 'sequence_24' and has a replicant count of 5.
        >sequence_24_5
        AAACG
    """
    name = "vsearch"

    def execute_program(self):
        args = self.args
        self.dereplicate_vsearch(args.input_f, args.outdir, args.groupsfile, args.processes, args.stripcounts,
                                 args.extraargstring)

    def dereplicate_vsearch(self, input_f, outdir, groupsfile, processes, stripcounts, extraargstring):
        """Dereplicates with vsearch.

        :param input_f: Filepath to the file or folder of files to dereplicate.
        :param outdir: Filepath to the output directory.
        :param groupsfile: A groups file to use as a reference for replicant counting.  If no groups file is
                            provided, input sequences are conidered singletons (regardless of their name-annotated
                            dereplication count).
        :param processes: The number of processes to use to dereplicate the fileset.
        :param stripcounts: If True, strips the trailing dereplication counts from a file before dereplication.
        :param extraargstring: Advanced program parameter string.
        """
        inputs = getInputFiles(input_f)
        pool = init_pool(min(len(inputs), processes))
        # REMOVES COUNTS FROM SEQUENCE NAMES IN ORDER TO CLUSTER PROPERLY
        # strip counts if we need to.
        if stripcounts:
            printVerbose("Removing counts from sequence names...")
            debugPrintInputInfo(inputs, "renamed")
            run_parallel([PythonRunner(removeCountsFromFastFile,
                                       [input_, "%s/%s_uncount.fasta" % (outdir, strip_ixes(input_)), 'fasta'],
                                       {"exists": input_})
                          for input_ in inputs], pool)
            printVerbose("Done removing counts.")

            # Grab the cleaned files as input for the next step
            inputs = getInputFiles(outdir, "*_uncount.fasta")

        # DEREPLICATE
        debugPrintInputInfo(inputs, "dereplicated")
        printVerbose("Dereplicating...")
        run_parallel([ProgramRunner(ProgramRunnerCommands.DEREP_VSEARCH,
                                    [processes, input_,
                                     "%s/%s_derep.fasta" % (outdir, strip_ixes(input_)),
                                     "%s/%s_uc.out" % (outdir, strip_ixes(input_))],
                                    {"exists": [input_], "positive": [processes]},
                                    extraargstring)
                      for input_ in inputs], pool)
        printVerbose("Done dereplicating")

        # LOG DEREPLICATED SEQUENCES INTO A .GROUPS FILE
        # generates a .groups file named _uc_parsed.out
        # python parseUCtoGroups.py uc.out uc_parsed.out
        input_ucs = getInputFiles(outdir, "*_uc.out")
        printVerbose("Generating a groups file from dereplication.")
        debugPrintInputInfo(inputs, "parsed (into a .groups file)")
        run_parallel([PythonRunner(parseUCtoGroups, [input_, "%s/%s_derep.groups" % (outdir, strip_ixes(input_))],
                                   {"exists": [input_]})
                      for input_ in input_ucs], pool)

        most_recent_groups_files = getInputFiles(outdir, "*_derep.groups", ignore_empty_files=False)

        # UPDATE THE MOST CURRENT GROUPS FILES WITH DEREPLICATION COUNTS
        if groupsfile is not None:
            # Grab the oldgroups file and the dereplicated groups file
            old_groups_files = getInputFiles(groupsfile)
            derep_groups_files = getInputFiles(outdir, "*_derep.groups")

            printVerbose("Updating .groups files with dereplicated data")
            printVerbose("%d Reference (old) groups files to be read:" % len(old_groups_files))
            printVerbose(str(old_groups_files))
            printVerbose("%d Dereplicated (new) groups files to be read:" % len(derep_groups_files))
            printVerbose(str(derep_groups_files))

            update_groups(old_groups_files, derep_groups_files, outdir, "dereplicated")
            most_recent_groups_files = getInputFiles(outdir, "dereplicated*", ignore_empty_files=False)
            printVerbose("Done updating .groups files.")

        if len(inputs) != len(most_recent_groups_files):
            print ("Error: Number of input fastas (%d) is not equal to the number ofgroups files (%d)." %
                   (len(inputs), len(most_recent_groups_files)))
            exit()
        fasta_groups_pairs = zip(inputs, most_recent_groups_files)
        # ADD COUNT TO SEQUENCE NAMES AND SORT BY COUNT
        # python renameWithReplicantCounts.py
        #               8_macse_out/MACSEOUT_MERGED.fasta uc_parsed.out dereplicated_renamed.fasta
        printVerbose("Adding dereplication data to unique fasta")
        run_parallel([PythonRunner(renameWithReplicantCounts,
                                   [fasta, groups, "%s/%s_counts.fasta" % (outdir, strip_ixes(fasta)), 'fasta'],
                                   {"exists": [fasta, groups]})
                      for fasta, groups in fasta_groups_pairs], pool)
        printVerbose("Done adding data")

        aux_dir = makeAuxDir(outdir)
        groups_dir = makeDirOrdie("%s_groups_files" % outdir)
        bulk_move_to_dir(most_recent_groups_files, groups_dir)
        aux_files = getInputFiles(outdir, '*', "*_counts.fasta", ignore_empty_files=False)
        bulk_move_to_dir(aux_files, aux_dir)
        cleanup_pool(pool)
