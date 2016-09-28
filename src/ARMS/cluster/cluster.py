import sys
from cluster_crop import cluster_crop
from cluster_swarm import cluster_swarm
from cluster_vsearch import cluster_vsearch
from classes.Helpers import *


# TODO doc
def cluster_main(input_f, outdir, program, groupsfile, threads, aux_params):
    """Switchboard for clustering options.

    :param input_f:
    :param outdir:
    :param program:
    :param groupsfile:
    :param threads:
    :param aux_params:
    :return:
    """
    # Make the output directory, failing if it already exists

    makeDirOrdie(outdir)
    # Grab the fasta file(s) to cluster_main
    inputs = getInputFiles(input_f)
    debugPrintInputInfo(inputs, "clustered")
    pool = init_pool(min(len(inputs), threads))
    printVerbose("\nClustering with: %s\n" % program)
    rslt = (None,None)
    if program == "crop":
        required = (inputs, outdir, groupsfile, pool, )
        keys = ["blocksize", "clustpct", "maxmcmc", "maxsm", "rare"]
        missing_param =False
        for key in keys:
            val = aux_params.get(key, None)
            if val is None:
                print "Error: '%s' parameter is requried for crop clustering" % key
                missing_param = True
        if not missing_param:
            params = list(required)
            for key in keys:
                params.append(aux_params.get(key, None))
            print params
            rslt = cluster_crop(*params)

    elif program == "vsearch":
        vsearch_id_pct = aux_params.get('idpct', None)
        if vsearch_id_pct is None:
            print "Error: 'idpct' parameter is required for vsearch clustering"
        rslt = cluster_vsearch(inputs, outdir, groupsfile, pool, vsearch_id_pct)

    else: #"swarm"
        rslt = cluster_swarm(inputs, outdir, groupsfile, pool)

    # Grab the resulting file lists from clustering
    aux_file_list, groups_file_list = rslt

    # Move the final groups file(s) to the groups dir
    groups_dir = makeDirOrdie("%s_groups_files" % outdir)
    bulk_move_to_dir(groups_file_list, groups_dir)

    # Move aux files to the aux dir
    aux_dir = makeAuxDir(outdir)
    bulk_move_to_dir(aux_file_list, aux_dir)

    # Cleanup the pool
    cleanup_pool(pool)




if __name__ == "__main__":
    if len(sys.argv) < 5:
        print "Usage: file_to_clean    output_file    list_of_chars_to_remove   filetype "
    else:
        cluster_main(*sys.argv[1:])
