#!/usr/bin/python

# Copyright 2022 Brian Tjaden
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

import sys, os
import subprocess
import random
from statistics import mean


#####################
#####   USAGE   #####
#####################

if len(sys.argv) < 2:
	sys.stderr.write("\nUSAGE: run_gffcompare.py <rna-seq-data/>" + "\n\n")
	sys.stderr.write("run_gffcompare combines assemblies and estimates accuracy. Output is to several files with prefix *E_coli* in the input directory.\n\n")
	sys.exit(1)



#########################
#####   FUNCTIONS   #####
#########################

def write_gtf_files_to_file(working_dir):
        gtf_list = []
        file_list = os.listdir(working_dir)
        for f in file_list:
                if (f.endswith('.gtf')) and (not f.startswith(BASE)) and (not f.startswith('Merged.annotated.')): gtf_list.append(working_dir + f)
        random.shuffle(gtf_list)  # Randomly shuffle list of GTF files

        # Write list of GTF files out to file to be used as GFFCOMPARE input
        with open(working_dir + BASE + '.filelist', 'w') as out_file:
                for f in gtf_list[:RANDOM_FILES]: out_file.write(f + '\n')



# Execute gffcompare on each file individually in order to get stats
def execute_gffcompare_once_per_file(working_dir):
        base_sensitivity, base_precision, exon_sensitivity, exon_precision, transcript_sensitivity, transcript_precision, locus_sensitivity, locus_precision, missed_exons, novel_exons, missed_loci, novel_loci = [], [], [], [], [], [], [], [], [], [], [], []
        with open(working_dir + BASE + '.filelist', 'r') as in_file:
                gtf_list = in_file.readlines()
        with open(working_dir + BASE + '.stats.all', 'w') as out_file:
                for f in gtf_list:
                        subprocess.run('./gffcompare -T -o ' + working_dir + BASE + ' -r genome/GCF_000005845.2_ASM584v2_genomic.gff ' + f, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                        with open(working_dir + BASE + '.stats', 'r') as in_file:
                                lines = in_file.readlines()
                        if (not lines[9].startswith('#-----------------|')): continue
                        for line in lines: out_file.write(line.strip() + '\n')
                        base_sensitivity.append(float(lines[10].split()[2]))
                        base_precision.append(float(lines[10].split()[4]))
                        exon_sensitivity.append(float(lines[11].split()[2]))
                        exon_precision.append(float(lines[11].split()[4]))
                        transcript_sensitivity.append(float(lines[14].split()[2]))
                        transcript_precision.append(float(lines[14].split()[4]))
                        locus_sensitivity.append(float(lines[15].split()[2]))
                        locus_precision.append(float(lines[15].split()[4]))
                        missed_exons.append(int(lines[21].split()[2].split('/')[0]))
                        novel_exons.append(int(lines[22].split()[2].split('/')[0]))
                        index = 25
                        if (lines[24].strip().startswith('Missed loci:')): index = 24
                        missed_loci.append(int(lines[index].split()[2].split('/')[0]))
                        novel_loci.append(int(lines[index+1].split()[2].split('/')[0]))
        os.remove(working_dir + BASE + '.annotated.gtf')

        # Output summary stats
        with open(working_dir + BASE + '.stats.summary', 'w') as out_file:
                out_file.write('Average base level sensitivity: ' + '{:.1f}'.format(mean(base_sensitivity)) + '\n')
                out_file.write('Average base level precision: ' + '{:.1f}'.format(mean(base_precision)) + '\n')
                out_file.write('Average exon level sensitivity: ' + '{:.1f}'.format(mean(exon_sensitivity)) + '\n')
                out_file.write('Average exon level precision: ' + '{:.1f}'.format(mean(exon_precision)) + '\n')
                out_file.write('Average transcript level sensitivity: ' + '{:.1f}'.format(mean(transcript_sensitivity)) + '\n')
                out_file.write('Average transcript level precision: ' + '{:.1f}'.format(mean(transcript_precision)) + '\n')
                out_file.write('Average locus level sensitivity: ' + '{:.1f}'.format(mean(locus_sensitivity)) + '\n')
                out_file.write('Average locus level precision: ' + '{:.1f}'.format(mean(locus_precision)) + '\n')
                out_file.write('Average number missed exons: ' + str(int(mean(missed_exons))) + '\n')
                out_file.write('Total number missed exons: ' + str(int(sum(missed_exons))) + '\n')
                out_file.write('Average number novel exons: ' + str(int(mean(novel_exons))) + '\n')
                out_file.write('Total number novel exons: ' + str(int(sum(novel_exons))) + '\n')
                out_file.write('Average number missed loci: ' + str(int(mean(missed_loci))) + '\n')
                out_file.write('Total number missed loci: ' + str(int(sum(missed_loci))) + '\n')
                out_file.write('Average number novel loci: ' + str(int(mean(novel_loci))) + '\n')
                out_file.write('Total number novel loci: ' + str(int(sum(novel_loci))) + '\n')



# Execute on all GTF files to get merged results
def execute_gffcompare_all(working_dir):
        subprocess.run('./gffcompare -T -o ' + working_dir + BASE + ' -r genome/GCF_000005845.2_ASM584v2_genomic.gff -i ' + working_dir + BASE + '.filelist', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        subprocess.run('./gffcompare -T -o ' + working_dir + BASE2 + ' -r genome/GCF_000005845.2_ASM584v2_genomic.gff ' + working_dir + BASE + '.combined.gtf', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)



##############################
##########   MAIN   ##########
##############################

RANDOM_FILES = 10000  # Randomly select this many GTF files to be used by GFFCOMPARE
BASE = 'E_coli'
BASE2 = 'Merged'

working_dir = sys.argv[1]
if (working_dir[-1] != '/'): working_dir += '/'
write_gtf_files_to_file(working_dir)
execute_gffcompare_once_per_file(working_dir)
execute_gffcompare_all(working_dir)


