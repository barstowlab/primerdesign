#!/sw/bin/python3.4



import pdb
import re

import sys
import datetime
import os

import subprocess

from sudokuutils6 import get_input2, ImportGenBankSequence, ImportFastaSequence, \
FindATandTAPositions2, ImportInsertionLocationFile2, FindClosestRealDisruptionIndex, \
ReverseComplement, WriteFastaSequence, ensure_dir, WritePrimer3InputForDetectionPrimerSelection, \
GenerateTrimmedGenomeSequence, ParsePrimer3DesignOutputFile

from vectorOutput2 import generateOutputMatrixWithHeaders, writeOutputMatrix


# ------------------------------------------------------------------------------------------------ #
# Design of primers for verification of Sudoku multiple and single occupancy predictions
# Last updated by Buz Barstow 2015-8-07
# ------------------------------------------------------------------------------------------------ #


####################################################################################################
####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Import the input parameter file

inputParameters = ['verificationPrimerBaseDir', 'shewyGenomeGenBankFileName', \
'primerLocationFile', 'transposonSequenceFileName']

# argv = sys.argv
# inputParameterValues = get_input2(argv)

inputParameterValues = get_input2(['','../../Input Files/PrimerDesign.inp'], inputParameters)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Parse the input data
verificationPrimerBaseDir = inputParameterValues['verificationPrimerBaseDir']
shewyGenomeGenBankFileName = inputParameterValues['shewyGenomeGenBankFileName']
primerLocationFile = inputParameterValues['primerLocationFile']
transposonSequenceFileName = inputParameterValues['transposonSequenceFileName']

chromosomeStart = 1
chromosomeEnd = 4969811
megaplasmidStart = 4969812
megaplasmidEnd = 5131424

chromosomeStartIndex = chromosomeStart - 1
chromosomeEndIndex = chromosomeEnd - 1
megaplasmidStartIndex = megaplasmidStart - 1
megaplasmidEndIndex = megaplasmidEnd - 1


primerThermodynamicParameterPath = '/Users/buz/Dropbox (BarstowLab)/bin/primer3-2.3.6/primer3_config/'

ensure_dir(verificationPrimerBaseDir)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
sequence = ImportGenBankSequence(shewyGenomeGenBankFileName)
transposonSequence = ImportFastaSequence(transposonSequenceFileName)
rcTransposonSequence = ReverseComplement(transposonSequence)
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Import the insertion location to check file
locDict = ImportInsertionLocationFile2(primerLocationFile)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Make a sequence with an inserted transposon and design primers to amplify this sequence

disruptionKeys = locDict.keys()
primerDesignOutputFilesDict = {}

for key in disruptionKeys:
	
	disruptionLocations = locDict[key][0].split(',')
	disruptionNames = locDict[key][1].split(',')
	
	i = 0
	while i < len(disruptionLocations):
	
		disruptionLoc = int(disruptionLocations[i])
		disruptionLocIndex = disruptionLoc - 1


		closestDisruptionIndex = FindClosestRealDisruptionIndex(disruptionLocIndex, sequence, \
		chromosomeStartIndex, chromosomeEndIndex, megaplasmidStartIndex, megaplasmidEndIndex)

		# Trim the genome sequence to the region around the transposon

		trimmedSequence, sequenceTarget, productSizeRange = \
		GenerateTrimmedGenomeSequence(closestDisruptionIndex, sequence, \
		chromosomeStartIndex, chromosomeEndIndex, megaplasmidStartIndex, megaplasmidEndIndex, \
		transposonSequence, rcTransposonSequence)


		name = str(disruptionLoc)

		primerDesignFileName = verificationPrimerBaseDir + '/' + str(disruptionLoc)+'.txt'
		primerOutputFile = verificationPrimerBaseDir + '/' + str(disruptionLoc)+'.primer'

		WritePrimer3InputForDetectionPrimerSelection(primerDesignFileName, trimmedSequence, name, \
		sequenceTarget, productSizeRange, primerThermodynamicParameterPath)

		callStr = ['primer3_core', '-output', primerOutputFile, primerDesignFileName]
	
		primerDesignOutputFilesDict[str(disruptionLoc)] = [primerOutputFile, disruptionNames[i]]

		subprocess.check_output(callStr)
	
		# Remove the temporary files placed in the working directory
		os.remove(name + '.for')
		os.remove(name + '.rev')
		
		i += 1


# ------------------------------------------------------------------------------------------------ #
# Parse the output files and make a list of primers pairs, melting temperatures and GC contents

rightPrimers = []
leftPrimers = []
rightPrimerTms = []
leftPrimerTms = []
rightPrimerGCContents = []
leftPrimerGCContents = []
disruptionLocations = []
disruptionNames = []

primerDesignOutputFileKeys = primerDesignOutputFilesDict.keys()

for key in primerDesignOutputFileKeys:

	file = primerDesignOutputFilesDict[key][0]
	disruptionName = primerDesignOutputFilesDict[key][1]
	
	rightPrimer, leftPrimer, rightPrimerTm, leftPrimerTm, rightPrimerGC, leftPrimerGC = \
	ParsePrimer3DesignOutputFile(file) 

	rightPrimers.append(rightPrimer)
	leftPrimers.append(leftPrimer)
	rightPrimerTms.append(rightPrimerTm)
	leftPrimerTms.append(leftPrimerTm)
	rightPrimerGCContents.append(rightPrimerGC)
	leftPrimerGCContents.append(leftPrimerGC)
	disruptionLocations.append(key)
	disruptionNames.append(disruptionName)


vectorList = [disruptionLocations, leftPrimers, rightPrimers, leftPrimerTms, rightPrimerTms, \
leftPrimerGCContents, rightPrimerGCContents, disruptionNames]

headers = ["Disruption Location", "Left Primer", "Right Primer", "Left Primer Tm", "Right Primer Tm", \
"Left Primer GC", "Right Primer GC", "Disruption Name"]

primerSummaryOutputFile = verificationPrimerBaseDir + '/' + 'primerSummary.csv'

oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')

writeOutputMatrix(primerSummaryOutputFile, oMatrix)

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Make an ordering list for IDT
fHandle = open(verificationPrimerBaseDir + '/' + 'Primer Summary for IDT.csv', 'w')

fHandle.write('Name,Sequence,Scale,Purification\n')

i = 0

while i < len(leftPrimers):
	fHandle.write(disruptionNames[i] + '_L' + ',' + leftPrimers[i] + ',' + '25nm' + ',' + 'STD' \
	+ '\n')
	
	fHandle.write(disruptionNames[i] + '_R' + ',' + rightPrimers[i] + ',' + '25nm' + ',' + 'STD' \
	+ '\n')
	
	i += 1

fHandle.close()
# ------------------------------------------------------------------------------------------------ #
