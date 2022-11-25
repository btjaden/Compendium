#!/usr/bin/python

# Copyright 2022 Brian Tjaden
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

import sys, os
import subprocess


#####################
#####   USAGE   #####
#####################

if len(sys.argv) < 2:
	sys.stderr.write("\nUSAGE: getMetaData.py <SraAccList.txt>" + "\n\n")
	sys.stderr.write("getMetaData takes a list of SRA Accession Numbers from the file data/SraAccList.txt given as a command line argument and downloads the metadata for each, outputting results to the files *metaData.txt* and metaDataSummary.txt*.\n\n")
	sys.exit(1)



#########################
#####   FUNCTIONS   #####
#########################

def getMetadata(filename):
        out_file = open('metaData.txt', 'w')
        header = 'Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash'
        out_file.write('\t'.join(header.split(',')) + '\n')
        headers = ['strategy', 'selection', 'source', 'layout', 'platform', 'model', 'scientificname']
        dictionaries = [{}, {}, {}, {}, {}, {}, {}]
        indices = [12, 13, 14, 15, 18, 19, 28]
        temp_filename = 'temp6789876.txt'
        count = 1
        with open(filename, 'r') as f:
                accession = f.readline().strip()
                while (accession != ''):
                        subprocess.run('esearch -db sra -query ' + accession + ' | efetch -format runinfo > ' + temp_filename, shell=True)
                        # Process metadata
                        with open(temp_filename, 'r') as f_temp:
                                foundData = False
                                line = f_temp.readline().strip()
                                line = f_temp.readline().strip()
                                while (line != ''):
                                        parse_line = line.split(',')
                                        if (parse_line[0] == accession) and (len(parse_line) > indices[-1]):
                                                foundData = True
                                                out_file.write('\t'.join(parse_line) + '\n')
                                                for i in range(len(indices)):
                                                        token = parse_line[indices[i]]
                                                        if (len(token) > 0):
                                                                dictionaries[i][token] = dictionaries[i].get(token, 0) + 1
                                        line = f_temp.readline().strip()
                                if (not foundData): sys.stderr.write('Could not obtain meta data for accession ' + accession + '\n')
                        accession = f.readline().strip()
                        if (count % 10 == 0): sys.stdout.write(str(count) + '\n')
                        count += 1
                        #if (count >=100): break
        out_file.close()
        os.remove(temp_filename)

        # Output metadata summary
        out_file = open('metaDataSummary.txt', 'w')
        for i in range(len(dictionaries)):
                out_file.write('Dictionary: ' + headers[i] + '\n')
                for key in dictionaries[i]:
                        out_file.write(key + '\t' + str(dictionaries[i][key]) + '\n')
                out_file.write('\n')
        out_file.close()



##############################
##########   MAIN   ##########
##############################

getMetadata(sys.argv[1])

