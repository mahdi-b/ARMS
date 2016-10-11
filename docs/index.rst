.. Chewbacca documentation master file, created by
   sphinx-quickstart on Tue Sep 13 17:53:16 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index:

Chewbacca - A Toolkit for COI Analysis
=======================================
Chewbacca is a bioinformatics toolkit for COI analysis, meshing commonly used programs together to create a framework for automated analysis.
Chewbacca currently supports the cleaning, assembly, demultiplexing, clustering, Aligning, and Identification of COI data.
Chewbacca also allows users to build OTU tables and visualize their data.

.. _`quick_start`:

Quick Start
-----------
1. Grab the Docker image version x.x.x.
2. In your docker shell:

::

	$ docker load -it chewbacca_vx.x.x
	# cd ~/ARMS/testARMS
	# python ~/ARMS/src/ARMS/chewbacca.py --help

.. _faq:

FAQ
===
**What is Chewbacca?**

Chewbacca is a command line bioinformatics toolkit for COI analysis, meshing commonly used programs together to create a framework for automated analysis.
Chewbacca currently supports the cleaning, assembly, demultiplexing, clustering, Aligning, and Identification of COI data.
Chewbacca also allows users to build OTU tables and visualize their data.

**What does it do?**

Lots of stuff.  Its like mothur and Qiime, but for COI analysis.  Have a look at our available :ref:`commands`.  

**Who should use Chewbacca?**

Researchers who use COI data (at any stage of processing) for abundance/distribution questions.

**I have some fasta files.  I need to clean them.  Can you help?**

Yes!  Chewbacca comes with a 'default' set of steps that will take in raw reads (or reads at varying levels of assembly/cleaning) and give you back an OTU table and some nice graphs.

**What makes Chewbacca Different?**

1. Chewbacca is a toolkit designed with run_parallel processing in mind.  Chewbacca's operations are as run_parallel as possible.
2. Chewbacca is modular.  Different research requires different analytical pipelines. Chewbacca's subprograms are each designed to tackle one small problem.  Odds are good that you'll find some parts of the toolkit useful.
3. Chewbacca remembers.  Did you mess up in one of your steps?  Chewbacca saves the output of each step to a directory, meaning you don't have to start from scratch if you change part of your pipeline.

**Don't like what you see?**

Adding new processes/programs to Chewbacca is easy!  Take a look at our :ref:`dev_guide`.


.. _reference:

Reference
=========
* :ref:`file_types`
* :ref:`commands`
* :ref:`API`
* :ref:`dev_guide`

