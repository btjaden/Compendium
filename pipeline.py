#!/usr/bin/python

# Copyright 2022 Brian Tjaden
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

import sys, os
import subprocess, shutil
import time


#####################
#####   USAGE   #####
#####################

if len(sys.argv) < 2:
	sys.stderr.write("\nUSAGE: pipeline.py <metaData.txt> <OPT batch_size>" + "\n\n")
	sys.stderr.write("pipeline takes a list of SRA Accession Numbers and executes a pipeline with the following steps. (1) Download FASTQ files from SRA using SRATOOLKIT. (2) Run HISAT2 on FASTQ files to determine alignments. SAM file is output. (3) Run SAMTOOLS to generate sorted BAM file. (4) Run STRINGTIE to determine assembly. GTF file is output. (5) The GTF files are then fed in to GFFCOMPARE using another Python script, run_gffcompare.py. The pipeline is set up to process a batch of *batch_size* files each time the pipeline is executed, where *batch_size* is an optional command line argument. Experiments that have been processed previously will not be processed again, so repeated execution of the pipeline should eventually process all experiments. Information about each file/experiment/alignment in stored in the file *rna-seq-dir/results.pipeline* and several additional useful files are output with in files with the prefix *rna-seq-dir/E_coli*.\n\n")
	sys.exit(1)



#########################
#####   FUNCTIONS   #####
#########################

# RETURNS A DICTIONARY OF PREVIOUSLY PROCESSED EXPERIMENTS
def prepareOutputFile():
        if (not os.path.exists(OUT_FILENAME)):
                out_file = open(OUT_FILENAME, 'w')
                out_file.close()
                return {}

        experiments = {}
        with open(OUT_FILENAME, 'r') as in_file:
                line = in_file.readline()
                while (line != ''):
                        if (line.startswith(HEADER)):
                                parse_line = line.strip().split()
                                experiments[parse_line[1]] = True
                        line = in_file.readline()
        sys.stderr.write('Number of experiments processed previously: ' + str(len(experiments)) + '\n')
        return experiments



def downloadFile(in_filename, experiments):
        name = ''
        isPaired = False
        totalNum = 0
        numFailed = 0
        numDownloaded = 0
        with open(in_filename, 'r') as f:
                line = f.readline()  # Ignore header line
                line = f.readline()
                while (line != '') and (numDownloaded == 0):
                        parse_line = line.split('\t')
                        name = parse_line[0]
                        isPaired = False
                        if (parse_line[15].lower() == 'paired'): isPaired = True
                        if ((parse_line[12] == 'RNA-Seq') or (parse_line[12] == 'ncRNA-Seq')) and (parse_line[14] == 'TRANSCRIPTOMIC') and (parse_line[18] == 'ILLUMINA') and (parse_line[28] == 'Escherichia coli str. K-12 substr. MG1655'):
                                totalNum += 1
                                if (name in experiments):
                                        line = f.readline()
                                        continue
                                prefetch_run = subprocess.run('./prefetch -q ' + name, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                                if ('err:' in prefetch_run.stderr.decode()):
                                        numFailed += 1
                                        line = f.readline()
                                        continue
                                fasterq_run = subprocess.run('./fasterq-dump -q ' + name, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                                if ('err:' in fasterq_run.stderr.decode()):
                                        numFailed += 1
                                        line = f.readline()
                                        continue
                                numDownloaded += 1
                                experiments[name] = True  # Add to dictionary
                        line = f.readline()
        return name, isPaired, numDownloaded



def getFileSize(name, isPaired):
        if (not isPaired) and (os.path.exists(name + '.fastq')):
                return os.path.getsize(name + '.fastq')
        if (isPaired) and (os.path.exists(name + '_1.fastq')) and (os.path.exists(name + '_2.fastq')):
                return os.path.getsize(name + '_1.fastq') + os.path.getsize(name + '_2.fastq')
        return -1


def run_HISAT(name, isPaired):
        fastq_files = '-U ' + name + '.fastq'
        if (isPaired): fastq_files = '-1 ' + name + '_1.fastq ' + '-2 ' + name + '_2.fastq'
        hisat_run = subprocess.run('./hisat2 -p 8 -x genome/E_coli ' + fastq_files + ' -S ' + name + '.sam', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        return hisat_run.stderr.decode()



def run_SAMTOOLS(name):
        temp_filename = 'mqixunecodkeju.trash'
        samtools_run = subprocess.run('./samtools sort --threads 8 -o ' + name + '.bam ' + name + '.sam', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        samtools_run = subprocess.run('./samtools stats --threads 8 ' + name + '.bam > ' + temp_filename, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        bam_stats = ''
        with open(temp_filename, 'r') as in_file:
                line = in_file.readline()
                while (line != ''):
                        if (line.startswith('SN')): bam_stats += line.strip() + '\n'
                        line = in_file.readline()
        os.remove(temp_filename)
        return bam_stats



def run_STRINGTIE(name):
        stringtie_run = subprocess.run('./stringtie -o ' + WORKING_DIR + name + '.gtf -p 8 -G genome/GCF_000005845.2_ASM584v2_genomic.gff -m 50 -s 1 -c 1 -g 5 -a 9999999 ' + name + '.bam', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)



##############################
##########   MAIN   ##########
##############################

BATCH_SIZE = 1  # Number of RNA-seq experiments to process
if (len(sys.argv) >= 3): BATCH_SIZE = int(sys.argv[2])
WORKING_DIR = 'rna-seq-data/'
OUT_FILENAME = WORKING_DIR + 'results.pipeline'
HEADER = '+++++'

# SET UP FILE TO OUTPUT RESULTS. AFTER FILE IS CREATED, WE ONLY
# APPEND TO IT. DETERMINE EXPERIMENTS PREVIOUSLY PROCESSED.
experiments = prepareOutputFile()

# PROCESS A BATCH OF RNA-SEQ EXPERIMENTS
for i in range(BATCH_SIZE):
        sys.stderr.write('\t' + str(i+1) + ' out of ' + str(BATCH_SIZE) + '\n')

        # DOWNLOAD FASTQ FILE FROM SRA
        t0 = time.time()
        name, isPaired, numDownloaded = downloadFile(sys.argv[1], experiments)
        t1 = time.time()
        if (numDownloaded > 0):
                sys.stderr.write('\tDownloaded:\t' + name + '\n')
                out_file = open(OUT_FILENAME, 'a')  # Add to file
                out_file.write(HEADER + ' ' + name + '\n')
                out_file.write('TIME (DOWNLOAD): ' + str(int(t1-t0)) + '\n')
                out_file.write('FILE SIZE: ' + str(getFileSize(name, isPaired)) + '\n')
                out_file.write('IS_PAIRED: ' + str(isPaired) + '\n')
                out_file.close()
        else:  continue  # Error downloading fastq file from SRA

        # EXECUTE HISAT2 ON FASTQ FILE
        t0 = time.time()
        hisat_results = run_HISAT(name, isPaired)
        if (hisat_results.startswith('Error') or hisat_results.startswith('Warning') or hisat_results.startswith('(ERR)')): hisat_results = '\n'
        t1 = time.time()
        sys.stderr.write('\tAligned:\t' + name + '\n')
        out_file = open(OUT_FILENAME, 'a')  # Add to file
        out_file.write('TIME (HISAT): ' + str(int(t1-t0)) + '\n')
        out_file.write('HISAT START\n' + hisat_results + '\nHISAT END\n')
        out_file.close()

        # EXECUTE SAMTOOLS ON SAM FILE TO CREATE SORTED BAM FILE
        t0 = time.time()
        samtools_results = run_SAMTOOLS(name)
        t1 = time.time()
        sys.stderr.write('\tBAM file:\t' + name + '\n')
        out_file = open(OUT_FILENAME, 'a')  # Add to file
        out_file.write('TIME (SAMTOOLS): ' + str(int(t1-t0)) + '\n')
        out_file.write(samtools_results + '\n')
        out_file.close()

        # EXECUTE STRINGTIE ON BAM FILE TO CREATE GTF FILE
        t0 = time.time()
        run_STRINGTIE(name)
        t1 = time.time()
        sys.stderr.write('\tAssembly:\t' + name + '\n')
        out_file = open(OUT_FILENAME, 'a')  # Add to file
        out_file.write('TIME (STRINGTIE): ' + str(int(t1-t0)) + '\n\n')
        out_file.close()

        # REMOVE TEMP FILES
        if (isPaired):
                if (os.path.exists(name + '_1.fastq')): os.remove(name + '_1.fastq')
                if (os.path.exists(name + '_2.fastq')): os.remove(name + '_2.fastq')
        elif (os.path.exists(name + '.fastq')): os.remove(name + '.fastq')
        if (os.path.exists(name + '.sam')): os.remove(name + '.sam')
        if (os.path.exists(name + '.bam')): os.remove(name + '.bam')
        sys.stderr.write('\n')  # Output blank line between each experiment



# CLEAN UP TEMP FILES
filelist = os.listdir()
for f in filelist:
        if (f.startswith('fasterq.tmp.tempest.')): os.rmdir(f)
        if (f.endswith('.fastq')): os.remove(f)

# Execute GFFCOMPARE
subprocess.run('python run_gffcompare.py ' + WORKING_DIR, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

