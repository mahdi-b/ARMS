
.. _`dev_guide`:

Developer's Guide
=================

Chewbacca Design Philosophies:
------------------------------
1. **Make everything parallel!**
	The goal of Chewbacca is to be a high-throughput pipeline, capable of scaling from 1 to 1000 processors.  As such, everything should be as parallel as possible. \
	Run multiple processes! Run mulitple threads! Split files into chunks if you have to!  

2. **Output goes in folders!** (Keep the working directory clean)
	Every Chewbacca command should make an output folder.  If that folder already exists, die.  Don't overwrite results.  If a Command writes files that will be neeed downstream,
	make a new folder for them (separately).  
	* \*_aux folders are for files that we don't need downstream, but someone might want to look at one day.
	* .groups files go in a \*_groups directory
	* .samples go in a \*_samples directory
	* etc.

3. **Commands do one thing.**
	To make Chewbacca as useful as possible, keep the work units small.  Do one thing (split, merge, cluster, dereplicate, etc.).  One big command may be useful to you, but its probably not applicable to everyone else.
	Before writing a command, ask yourself: "Will this be useful to someone else?"

4. **Programs do all the work.**
	Commands represent high-level operations (cluster, dereplicate, rename).  \
	But under the hood, this could be an involved, messy, complicated process.  Programs should take care of their own buisiness.  \
	If you need to move a file, rename a file, parse a file, or whatever, do it within the program class (side scripts and helpers are encouraged).

Chewbacca Design Pattern:
-------------------------
**ChewbaccaCommands (classes/ChewbaccaCommand.py)**

ChewbaccaCommands represent high-level commands like "cluster" (Cluster_Command), "partition" (Partition_Command), "dereplicate" (Dereplicate_Command), \
and know which programs are supported (ChewbaccaCommand.supported_programs), but otherwise oblivious as to HOW those programs operate.

**ChewbaccaPrograms (classes/ChewbaccaPrograms.py)**

ChewbaccaPrograms contain implementation-specific functionality required to complete a ChewbaccaCommand using a particular program. \
For example, Assemble_Program_Pear implements the functionality requried to complete the Assemble_Command using the program Pear.

**ProgramRunner (classes/ProgramRunner.py)**

ProgramRunner objects handle the execution of a command line program (like Pear).  \
The class contains a dictionary of unformatted command line strings, into which arguments are inserted.

**argparse**

Chewbacca uses argparse to take in all command line options.

Chewbacca Naming Conventions:
-----------------------------
* Command classes should end in "_Command"
* Program classes should end in "<Command_Name>_Program_<X>", where:
	 <Command_Name> is the name of the Command your program implements (conceptually, not programatically), and
	<X> is the sub program name you're using (mothur, Qiime, PEAR, etc), or 'Chewbacca' if you're using a custom script.


Adding New Programs:
--------------------
1. Add your command line invocation to ProgramRunner
	1.1. Add a ProgramRunnerProgram enum.
		This allows users to define paths to the executable.
	1.2. Add the executable path to ProgramRunner.program_paths.
		This provides a default path to the executable.
	1.3. Add a ProgramRunnerCommand enum.
		This exposes your command to the world.
	1.4. Add your command string to the ProgramRunner.commandTemplates dictionary
		This makes your command actually do something.

2. Add a parser for your command.
	Add a subparser to the argparse object, and fill in the variables and help messages you'll expose to users.
3. Create a Command class.
	Copy one of the other Command classes and rename it with your new command!
4. Create a Program class
	Create a program class and fill in implementation details (this is where your custom code goes).
	The classes/Helpers module provides some useful functions for selecting, moving, and creating files and folders.
	NOTE: If you create a multiprocessing.pool object, make sure you call classes/Helpers.cleanup_pool() on it.
5. Add your new Program Class to your Command class.
	Add your Program class to the list of supported programs.
	Choose a default Program for your Command (probably your new Program Class).
6. Test it!
	Run chewbacca.py with your command name and see how it goes!  Diagnose errors as appropriate.



