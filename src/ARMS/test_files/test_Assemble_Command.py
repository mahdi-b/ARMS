from assemble.Assemble_Command import Assemble_Command
from classes.Helpers import getInputFiles, makeDirOrdie, makeAuxDir
from test_files.testHelpers import cleanup_files, assert_auxdir, assert_outdir,  fasta_to_list
from nose.tools import assert_equals, assert_true

def test_assemble():
    """
    Test the Assemble_Command to ensure that: \n
    1. command creates just one output file in outdir. \n
    2. output file is correct (1 assembled read). \n
    3. command creates aux dir.

    parser_assemble.add_argument('-f', '--input_f', required=True, help="Forward Fastq Reads file or folder.")
    parser_assemble.add_argument('-r', '--input_r', required=True, help="Reverse Fastq Reads file or folder.")
    parser_assemble.add_argument('-n', '--name', required=True, help="Assembled File Prefix.")
    parser_assemble.add_argument('-o', '--outdir', required=True, help="Directory where outputs will be saved.")
    """

    params = lambda: 0
    test_reads = "assemble_data/test_%s.fq"
    params.input_f = test_reads % "R1"
    params.input_r = test_reads % "R2"
    params.outdir = "rslt"
    params.name = "test"
    params.processes = 1
    params.pearthreads = 1
    params.extraargstring = ""

    for program in Assemble_Command.supported_programs:
        cleanup_files(params.outdir)
        makeDirOrdie(params.outdir)
        params.program = program.name
        Assemble_Command(params).execute_command()
        assert_outdir(params.outdir)
        output_files = getInputFiles(params.outdir, "*assembled*")
        assert_equals(len(output_files), 1)
        assert_auxdir(params.outdir)
        seqs = fasta_to_list(output_files[0], 'fastq')
        assert_equals(len(seqs), 1)
        for seq in seqs:
            assert_true("good" in seq.id)
    return True