#!/usr/bin/python

# Copyright 2022 Brian Tjaden
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

import sys, os
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn import metrics
from sklearn import linear_model
from sklearn import ensemble


#####################
#####   USAGE   #####
#####################

if len(sys.argv) < 2:
	sys.stderr.write("\nUSAGE: operons.py <rna-seq-data/>" + "\n\n")
	sys.stderr.write("operons predicts operons based on all experiments. Summary info is output is to stdout and two files are created in the specified input directory: operonPairs.tsv and operonsComplete.tsv.\n\n")
	sys.exit(1)



#########################
#####   FUNCTIONS   #####
#########################

def readInGenes(gene_filename):
        genes = {}
        features = {'CDS':True, 'ncRNA':True, 'rRNA':True, 'tRNA':True, 'pseudogene':True, 'sequence_feature':True}
        with open(gene_filename, 'r') as in_file:
                line = in_file.readline()
                while (line != ''):
                        if (line.startswith('#')):  # Ignore comments
                                line = in_file.readline()
                                continue
                        parse_line = line.strip().split('\t')
                        feature = parse_line[2]
                        start = parse_line[3]
                        stop = parse_line[4]
                        strand = parse_line[6]
                        fields = parse_line[8]
                        parse_fields = fields.split(';')
                        ID, name, product = '?', '?', '?'
                        if (feature in features):
                                for p in parse_fields:
                                        ID_prefix = 'Parent='
                                        if (feature == 'pseudogene') or (feature == 'sequence_feature'): ID_prefix = 'ID='
                                        if (p.startswith(ID_prefix)): ID = p[len(ID_prefix):]
                                        if (p.startswith('gene=')): name = p[5:]
                                        if (p.startswith('product=')): product = p[8:]
                                if (len(ID) > 1): genes[ID] = [name, feature, start, stop, strand, product]
                        line = in_file.readline()
        genes2 = {}  # Dictionary keyed by gene name rather than Blattner number
        for g in genes: genes2[genes[g][0]] = g
        return genes, genes2



# Start coordinate is element at index 3 in list
def startCoord(myList):
        return myList[3]



# Start coordinate is element at index 2 in list
def startCoord2(myList):
        return myList[2]



# Create list of genes sorted by their start coordinate
def sortGenes(genes):
        list_of_genes = []
        for g in genes:
                list_of_genes.append([g, genes[g][0], genes[g][1], int(genes[g][2]), int(genes[g][3]), genes[g][4], genes[g][5]])
        list_of_genes.sort(key=startCoord)
        return list_of_genes



# Get list of consecutive genes on the same strand, i.e., candidate operons
# And consecutive genes on the opposite strand, i.e., negative operons
def getDirectonPairs(genes):
        list_of_genes = sortGenes(genes)

        # Extract pairs of consecutive genes on same strand and opposite strand
        genePairs, genePairs_negative = [], []
        for i in range(0, len(list_of_genes)-1):
                gene1, gene2 = list_of_genes[i], list_of_genes[i+1]
                if (gene1[1] == 'hokC') or (gene2[1] == 'hokC'): continue
                if (gene1[2] == 'CDS') and (gene2[2] == 'CDS'):
                        if (gene1[5] == gene2[5]): genePairs.append((gene1[1], gene2[1]))
                        else: genePairs_negative.append((gene1[1], gene2[1]))
        return genePairs, genePairs_negative



# Return a list of genes between the specified start and stop coordinates
def getGenesInOperon(list_of_genes, start, stop, strand):
        current_strand = '-'
        if (strand.lower() == 'forward'): current_strand = '+'

        # Binary Search
        mid = -1
        low, high = 0, len(list_of_genes) - 1
        while (low < high):
                mid = int((low + high) / 2)
                if (list_of_genes[mid][3] == start): break
                elif (list_of_genes[mid][3] < start): low = mid + 1
                else: high = mid - 1
        while (list_of_genes[mid][3] >= start) and (mid > 0): mid -= 1
        if (list_of_genes[mid][3] < start): mid += 1

        genesInOperon = []
        while (start <= list_of_genes[mid][3]) and (list_of_genes[mid][4] <= stop):
                if (list_of_genes[mid][5] == current_strand):
                        genesInOperon.append(list_of_genes[mid][1])
                mid += 1
        return genesInOperon



def readInRegulonDB(in_filename, genes):

        # Sort genes by start coordinate
        list_of_genes = sortGenes(genes)

        regulonDB = []
        with open(in_filename, 'r') as in_file:
                line = in_file.readline()
                while (line.startswith('#')): line = in_file.readline()
                while (line != ''):
                        parse_line = line.strip().split('\t')
                        numGenesInOperon = int(parse_line[4])
                        if (numGenesInOperon > 1):
                                start = int(parse_line[1])
                                stop = int(parse_line[2])
                                strand = parse_line[3]
                                genesInOperon1 = getGenesInOperon(list_of_genes, start, stop, strand)
                                genesInOperon2 = parse_line[5].split(',')  # Unused
                                evidence = '?'
                                if (len(parse_line) > 6): evidence = parse_line[6]
                                confidence = '?'
                                if (len(parse_line) > 7): confidence = parse_line[7]
                                if (len(genesInOperon1) > 1):
                                        regulonDB.append([start, stop, strand, genesInOperon1, evidence, confidence])
                        line = in_file.readline()
        return regulonDB



def getExpression(parse_line):
        numExpressed = 0
        expressionVector = []
        for i in range(4, len(parse_line)):
                if (len(parse_line[i]) > 1):
                        parse_expression = parse_line[i].split('|')
                        FPKM = float(parse_expression[3])
                        TPM = float(parse_expression[4])
                        numExpressed += 1
                        expressionVector.append(TPM)
                else: expressionVector.append(0.0)
        return numExpressed, expressionVector



def readInTranscripts(in_filename):
        annotated_transcripts = {}
        with open(in_filename, 'r') as in_file:
                line = in_file.readline()
                while (line != ''):
                        parse_line = line.strip().split('\t')
                        class_code = parse_line[3]
                        numExpressed, expressionVector = getExpression(parse_line)
                        if (parse_line[2] != '-') and (class_code != 'x'):  # Annotated gene
                                ID = parse_line[2].split('|')[1]
                                name = parse_line[2].split('|')[0]
                                if (ID.startswith('rna-')):
                                        ID = name
                                        name = 'rna'
                                if (ID not in annotated_transcripts) or (numExpressed > annotated_transcripts[ID][2]):
                                        annotated_transcripts[ID] = [name, class_code, numExpressed, expressionVector]
                        line = in_file.readline()
        return annotated_transcripts



# Retuen a list of distances between each pair of genes
def getFeatures(genes, genes2, transcripts, operonPairs, label):
        features = []
        for pair in operonPairs:
                ID1, ID2 = genes2[pair[0]], genes2[pair[1]]
                if (ID1 not in transcripts) or (ID2 not in transcripts): continue
                stop1 = int(genes[ID1][3])
                start2 = int(genes[ID2][2])
                stop2 = int(genes[ID2][3])
                distance = start2 - stop1
                if (stop2 <= stop1): distance = 0  # One gene is contained in the other
                expressionVector1 = transcripts[ID1][3]
                expressionVector2 = transcripts[ID2][3]
                corr = np.corrcoef(expressionVector1, expressionVector2)[0,1]
                features.append((distance, corr, label))
        return features



def outputNovelOperonPredictions(genes, genes2, regulonDB, transcripts):
        # Create dictionary mapping gene names to the operon (ID) in which they reside
        regulonDB_operons = {}  # Key is gene name, value is operon ID
        regulonDB_operons2 = []  # Index is operon ID, values are lists of genes in operon
        operonID = 0
        for op in regulonDB:
                geneList = op[3]
                regulonDB_operons2.append([])
                for g in geneList:
                        regulonDB_operons[g] = operonID
                        regulonDB_operons2[-1].append(g)
                operonID += 1

        directonPairs, directonPairs_negative = getDirectonPairs(genes)

        examplesPositive, examplesNegative = [], []
        for pair in directonPairs:
                ID1, ID2 = genes2[pair[0]], genes2[pair[1]]
                if (ID1 not in transcripts) or (ID2 not in transcripts): continue
                if (pair[0] in regulonDB_operons) and (pair[1] in regulonDB_operons) and (regulonDB_operons[pair[0]] == regulonDB_operons[pair[1]]):
                        examplesPositive.append(pair)
                else: examplesNegative.append(pair)
        featuresPositive = getFeatures(genes, genes2, transcripts, examplesPositive, 1)
        featuresNegative = getFeatures(genes, genes2, transcripts, examplesNegative, 0)

        # Predict operons
        DATA = np.array(featuresPositive + featuresNegative)
        X = DATA[:,:-1]
        y = DATA[:,-1]
        X = preprocessing.StandardScaler().fit_transform(X)  # Feature scaling
        model = linear_model.LogisticRegression(solver='lbfgs', C=0.0011)
        model.fit(X, y)
        y_pred = model.predict(X)
        y_pred_probs = model.predict_proba(X)
        sys.stdout.write('\nAccuracy:\t' + str(model.score(X, y)) + '\n')
        sys.stdout.write('F1 score:\t' + str(metrics.f1_score(y, y_pred)) + '\n')
        sys.stdout.write('Sensitivity:\t' + str(metrics.recall_score(y, y_pred)) + '\n')
        sys.stdout.write('Precision:\t' + str(metrics.precision_score(y, y_pred)) + '\n')
        sys.stdout.write('Operon pairs in RegulonDB:\t' + str(len(featuresPositive)) + '\n')
        sys.stdout.write('Novel predicted operon pairs:\t' + str(metrics.confusion_matrix(y, y_pred)[0][1]) + '\n')
        sys.stdout.write('\n')

        # Output operon *pairs*
        out_file = open(results_dir + 'operonPairs.tsv', 'w')
        out_file.write('Co-transcription Evidence' + '\t' + 'Probability of Co-transcription' + '\t' + 'RegulonDB' + '\t' + 'Gene1 name' + '\t' + 'Gene1 synonym' + '\t' + 'Gene1 start' + '\t' + 'Gene1 stop' + '\t' + 'Gene1 strand' + '\t' + 'Gene1 product' + '\t' + 'Gene2 name' + '\t' + 'Gene2 synonym' + '\t' + 'Gene2 start' + '\t' + 'Gene2 stop' + '\t' + 'Gene2 strand' + '\t' + 'Gene2 product' + '\n')
        genesInPairs = {}
        for i in range(len(examplesPositive)):
                predicted = 'NO'
                if (int(y_pred[i]) == 1): predicted = 'YES'
                ID1, ID2 = genes2[examplesPositive[i][0]], genes2[examplesPositive[i][1]]
                if (int(y_pred[i]) == 1): genesInPairs[ID1], genesInPairs[ID2] = True, True
                out_file.write(predicted + '\t' + str(y_pred_probs[i][1]) + '\t' + 'YES' + '\t' + ID1[5:] + '\t' + genes[ID1][0] + '\t' + genes[ID1][2] + '\t' + genes[ID1][3] + '\t' + genes[ID1][4] + '\t' + genes[ID1][5] + '\t' + ID2[5:] + '\t' + genes[ID2][0] + '\t' + genes[ID2][2] + '\t' + genes[ID2][3] + '\t' + genes[ID2][4] + '\t' + genes[ID2][5] + '\n')
        for i in range(len(examplesNegative)):
                predicted = 'NO'
                if (int(y_pred[len(examplesPositive) + i]) == 1): predicted = 'YES'
                ID1, ID2 = genes2[examplesNegative[i][0]], genes2[examplesNegative[i][1]]
                if (int(y_pred[len(examplesPositive) + i]) == 1): genesInPairs[ID1], genesInPairs[ID2] = True, True
                out_file.write(predicted + '\t' + str(y_pred_probs[len(examplesPositive) + i][1]) + '\t' + 'NO' + '\t' + ID1[5:] + '\t' + genes[ID1][0] + '\t' + genes[ID1][2] + '\t' + genes[ID1][3] + '\t' + genes[ID1][4] + '\t' + genes[ID1][5] + '\t' + ID2[5:] + '\t' + genes[ID2][0] + '\t' + genes[ID2][2] + '\t' + genes[ID2][3] + '\t' + genes[ID2][4] + '\t' + genes[ID2][5] + '\n')
        out_file.close()
        sys.stdout.write('Genes in operon pairs:\t' + str(len(genesInPairs)) + '\n')

        # Output *complete* operons
        operonPairs = []
        examplesAll = examplesPositive + examplesNegative
        for i in range(len(examplesAll)):
                if (int(y_pred[i]) == 0): continue
                ID1, ID2 = genes2[examplesAll[i][0]], genes2[examplesAll[i][1]]
                name1, start1, stop1, strand1, product1 = genes[ID1][0], genes[ID1][2], genes[ID1][3], genes[ID1][4], genes[ID1][5]
                name2, start2, stop2, strand2, product2 = genes[ID2][0], genes[ID2][2], genes[ID2][3], genes[ID2][4], genes[ID2][5]
                operonPairs.append((ID1, name1, int(start1), int(stop1), strand1, product1, ID2, name2, int(start2), int(stop2), strand2, product2))
                operonPairs.sort(key=startCoord2)  # Sort operon pairs by first gene coord
        genesToOperons = {}
        operonsToGenes = []  # Each list item is dictionary of genes in the operon
        operonIndex = 0
        for pair in operonPairs:
                name1, name2 = pair[1], pair[7]
                if (name1 not in genesToOperons) and (name2 not in genesToOperons):
                        operonsToGenes.append({name1:True, name2:True})
                        genesToOperons[name1], genesToOperons[name2] = operonIndex, operonIndex
                        operonIndex += 1
                elif (name1 not in genesToOperons) and (name2 in genesToOperons):
                        currentOperonIndex = genesToOperons[name2]
                        genesToOperons[name1] = currentOperonIndex
                        operonsToGenes[currentOperonIndex][name1] = True
                elif (name1 in genesToOperons) and (name2 not in genesToOperons):
                        currentOperonIndex = genesToOperons[name1]
                        genesToOperons[name2] = currentOperonIndex
                        operonsToGenes[currentOperonIndex][name2] = True
                elif (name1 in genesToOperons) and (name2 in genesToOperons):
                        sys.stdout.write('Error - not expecting both genes to already be in operons.\n')
        sys.stdout.write('Number of multi-gene operons:\t' + str(len(operonsToGenes)) + '\n')
        sys.stdout.write('Genes in multi-gene operons:\t' + str(len(genesToOperons)) + '\n')

        out_file = open(results_dir + 'operonsComplete.tsv', 'w')
        out_file.write('Start' + '\t' + 'Stop' + '\t' + 'Strand' + '\t' + 'Gene names' + '\t' + 'Gene synonyms' + '\t' + 'Gene products' + '\n')
        for genes_in_operon in operonsToGenes:
                list_of_genes = []
                for g in genes_in_operon: list_of_genes.append([genes2[g][5:], g, genes[genes2[g]][4], int(genes[genes2[g]][2]), int(genes[genes2[g]][3]), genes[genes2[g]][5]])
                list_of_genes.sort(key=startCoord)
                start = list_of_genes[0][3]
                stop = list_of_genes[-1][4]
                strand = list_of_genes[0][2]
                nameList, synonymList, productList = [], [], []
                for gene in list_of_genes:
                        nameList.append(gene[0])
                        synonymList.append(gene[1])
                        productList.append(gene[5])
                out_file.write(str(start) + '\t' + str(stop) + '\t' + strand + '\t' + ','.join(nameList) + '\t' + ','.join(synonymList) + '\t' + '... '.join(productList) + '\n')
        out_file.close()



##############################
##########   MAIN   ##########
##############################

working_dir = sys.argv[1]
if (working_dir[-1] != '/'): working_dir += '/'
results_dir = working_dir

genes, genes2 = readInGenes('genome/GCF_000005845.2_ASM584v2_genomic.gff')
transcripts = readInTranscripts(working_dir + 'E_coli.tracking')

# RegulonDB experimentally confirmed operons
regulonDB = readInRegulonDB('genome/RegulonDB.experiments', genes)
outputNovelOperonPredictions(genes, genes2, regulonDB, transcripts)

