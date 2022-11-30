
This directory contains code associated with the manuscript:
Escherichia coli Transcriptome Assembly from a Compendium of RNA-seq Data Sets

The code here implements an analysis pipeline that uses a number of third-party software tools. These tools must be installed and be available for execution from the current path for the pipeline to operate properly. The tools are:

esearch from NCBI's E-utilities
prefetch and fasterq-dump from SRA-Toolkit
HISAT2
SAMtools
StringTie
GffCompare

The code is provided under the GNU General Public License (see file COPYING for details).

*******************
*** Quick Start ***
*******************

0. Make sure the abovementioned third-party software tools are executable from the current path.

1. The provided file SraAccList.txt contains a list of accession numbers downloaded from SRA and used in this study. You can use this list (to reproduce results in this study) or download your own set of accession numbers from SRA.

2. Obtain meta data associated with each accesion number:
   python getMetaData.py SraAccList.txt

3. Execute the pipeline. This will process a batch of RNA-seq datasets. The pipeline downloads FASTQ files, aligns them to the genome, converts them to BAM format, assembles them, and merges the assemblies. Results are written to files in the directory rna-seq-data.
   python pipeline.py metaData.txt

4. Identification of operons based on the transcript assemblies and comparison to experimentally confirmed operons reported by RegulonDB:
   python operons.py rna-seq-data

