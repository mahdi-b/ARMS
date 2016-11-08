from classes.Helpers import getInputFiles, printVerbose, debugPrintInputInfo
from util.updateGroups import update_groups


def handle_groups_file_update(outdir, groupsfile, clustering_groups_files_uncount):
    """Checks if the user specified groups file exists, and updates the groupings with clustering data.
        Returns a list of the most up to date groups files.
    :param outdir: Filepath to the directory where outputfiles will be written.
    :param groupsfile: Optional filepath to the .groups file, or a folder of .groups files to use as a reference.
    :param clustering_groups_files_uncount: The output groups file from clustering, with trailing replication counts
        removed from sequence names.  Names in this file should match those used in the user-specified groups file
        groupsfile.
    :return: A list of filenames pointing to the most up to date groups files.
    """
    most_recent_groups_files = clustering_groups_files_uncount
    if groupsfile:
        # Try to grab groups files
        user_specified_groups_files = getInputFiles(groupsfile, critical=False)
        # If we have files at the given location
        if len(user_specified_groups_files) != 0:
            most_recent_groups_files = user_specified_groups_files
            printVerbose("Updating .groups files with clustering data")
            debugPrintInputInfo(most_recent_groups_files, "used as groups references")
            update_groups(most_recent_groups_files, clustering_groups_files_uncount, outdir, "postcluster")
            printVerbose("Done updating .groups files.")
            most_recent_groups_files = getInputFiles(outdir, "postcluster*.groups")
    else:
        printVerbose("No name files provided, assuming singletons...\n")
    return most_recent_groups_files
