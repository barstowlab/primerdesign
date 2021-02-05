# ------------------------------------------------------------------------------------------------ #
# sudokuUtils.py
# Created by Buz Barstow 2015-6-14
# Last updated by Buz Barstow 2015-6-23
# Utility code for analysis of Sudoku sequencing data
# ------------------------------------------------------------------------------------------------ #

gSeparatorString = '# ----------------------------------------------------------------------------------------------- #\n'



# ------------------------------------------------------------------------------------------------ #
# Input parser
def get_input(argv, inputParameters):
	import re
	import pdb
	
	if len(argv) != 2:
		print("IndexSummarize takes 1 argument, the input file name")
		return
	
	file = open(argv[1],'r')
	inputData = file.readlines()
	
	# Go through all lines and search for input parameters	
	
	inputParameterValues = {}
	
	for item in inputParameters:
		result = inputParameterSearch(inputData, item)
		if result != None:
			inputParameterValues[item] = result
		#else:
		#	print "Paramter " + item + \
		#	" not found. All parameters need to be defined in the input file"
		#	return None
	
	return inputParameterValues
# ------------------------------------------------------------------------------------------------ #
			
# ------------------------------------------------------------------------------------------------ #
# Function to search input data for a given key
def inputParameterSearch(inputData, inputItem):
	import re
	import pdb
	
	#from string import rstrip
	regex = re.compile('\s*' + inputItem + '\s*' + '=', re.IGNORECASE)
	
	inputFound = False
	
# 	pdb.set_trace()
	
	for line in inputData:
		if regex.match(line) != None:
			if inputFound == False:
				inputLine = line.split('=')
				inputParameter = inputLine[1].strip()
				inputFound = True
			elif inputFound == True:
				print("An input parameter can only be defined once")
				inputParameter = None
			
	if inputFound == False:
		inputParameter = None
	
# 	pdb.set_trace()
	
	return inputParameter
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Function to search input data for a given key
def inputParameterSearch2(inputData, inputItem):
	import re
	import pdb
	import numpy
	
	#from string import rstrip
	regex = re.compile('\s*' + inputItem + '\s*' + '='+ '\s*', re.IGNORECASE)
	
	inputFound = False
	continuation = False
	
	for line in inputData:
		if line[0] != '#':
			if continuation == False:
				if regex.match(line) != None:
					if inputFound == False:
						inputLine = line.split('=')
					
						inputFound = True
					
						if line.strip()[-1] == '\\':
	# 						pdb.set_trace()
							continuation = True
							inputParameter = inputLine[1].strip()[0:-1]
						else:
							continuation = False
							inputParameter = inputLine[1].strip()
				
					elif inputFound == True:
						print("An input parameter can only be defined once")
						inputParameter = None
			
			elif continuation == True:
				if line.strip()[-1] == '\\':
					continuation = True
					inputParameter += line.strip()[0:-1]
				else:
					continuation = False
					inputParameter += line.strip()
			
	if inputFound == False:
		inputParameter = None
	

	return inputParameter
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
# Input parser
def get_input2(argv, inputParameters):
	import re
	import pdb
	
	if len(argv) != 2:
		print("IndexSummarize takes 1 argument, the input file name")
		return
	
	file = open(argv[1],'r')
	inputData = file.readlines()
	
	# Go through all lines and search for input parameters	
	
	inputParameterValues = {}
	
	for item in inputParameters:
		result = inputParameterSearch2(inputData, item)
		if result != None:
			inputParameterValues[item] = result
	
	return inputParameterValues
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
class Alignment:

	def __init__(self, readID, readAlignmentCoord, alignmentQuality, alignmentFound, \
	multipleAlignmentsFound, strangeFlagsSum):
		self.readID = readID
		self.readAlignmentCoord = readAlignmentCoord
		self.alignmentQuality = alignmentQuality
		self.alignmentFound = alignmentFound
		self.multipleAlignmentsFound = multipleAlignmentsFound
		self.strangeFlagsSum = strangeFlagsSum


# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ParseBowtie2Output(alignmentLines, reverseCorrection=-1):
# The key update to version 2 of the AlignSequencesToGenome code is in here. 
# We take a look at the second field of the alignment line

# Note, I've added the option to add 1 bp to the readAlignmentCoord if the read aligns to the 
# reverse strand. I can easily take this out in the future if I redesign the primers to avoid 
# binding to the ITR at the opposite end of the transposon. 
	
	import re
	import pdb
	
	xsRe = re.compile('XS:i:')
	
	lowScoreCount = 0
	noAlignCount = 0
	perfectAlignmentCount = 0
	multipleAlignmentCount = 0
	xsCount = 0
	xsIndices = []
	lowScoreIndices = []
	noAlignIndices = []
	readsWithStrangeFlagsSums = 0
	
	genomeAlignmentDict = {}
	
	i = 0
	while i < len(alignmentLines):
		line = alignmentLines[i]
		if len(line) > 0:
			if line[0] != '@':
				if xsRe.search(line) != None:
					xsCount += 1
					xsIndices.append(i)
					multipleAlignmentsFound = True
				else:
					multipleAlignmentsFound = False
					

				lineData = line.split('\t')
				readID = lineData[0]
				flagsSum = int(lineData[1])
				readAlignmentCoord = int(lineData[3])
				alignmentQuality = int(lineData[4])
				alignmentSequenceLength = len(lineData[9])
				
# 				pdb.set_trace()
				
				
				if flagsSum == 0 or flagsSum == 16:
					# Flags sum should be zero as we are analyzing non paired reads, if 
					# read aligns to forward strand.
					if flagsSum == 0:
						readAlignmentCoord = readAlignmentCoord
						strangeFlagsSum = False
					# Flags sum should be 16 as we are analyzing non paired reads, if 
					# read aligns to reverse strand.
					elif flagsSum == 16:
						readAlignmentCoord = readAlignmentCoord + alignmentSequenceLength \
						+ reverseCorrection
						# Note, I've added the option to add 1 bp here to normalize the result with 
						# the Sanger read. 
						strangeFlagsSum = False
				else:
					# This really shouldn't be needed, but in case it is, we have a way to 
					# report a strange flags sum. 
					readAlignmentCoord = readAlignmentCoord
					readsWithStrangeFlagsSums += 1
					strangeFlagsSum = True
				
				if alignmentQuality <= 1:
					lowScoreCount += 1
					lowScoreIndices.append(i)
							
				if readAlignmentCoord == 0:
					noAlignIndices.append(i)
					noAlignCount += 1
					alignmentFound = False
				else:
					alignmentFound = True
					if multipleAlignmentsFound == False:
						perfectAlignmentCount += 1
					elif multipleAlignmentsFound == True:
						multipleAlignmentCount += 1

			
				genomeAlignmentDict[readID] = Alignment(readID, readAlignmentCoord, \
				alignmentQuality, alignmentFound, multipleAlignmentsFound, strangeFlagsSum)				
		
		i += 1
	return [genomeAlignmentDict, perfectAlignmentCount, multipleAlignmentCount, lowScoreCount, \
	noAlignCount, readsWithStrangeFlagsSums]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def OutputGenomeAlignmentSummary(outputFileName, genomeAlignmentDict):
# Updated to include strange flags sum
	keys = genomeAlignmentDict.keys()

	outputHandle = open(outputFileName, 'w')
	
	for key in keys:
	
		readAlignmentCoord = str(genomeAlignmentDict[key].readAlignmentCoord)
		alignmentQuality = str(genomeAlignmentDict[key].alignmentQuality)
		alignmentFound = str(genomeAlignmentDict[key].alignmentFound)
		multipleAlignmentsFound = str(genomeAlignmentDict[key].multipleAlignmentsFound)
		strangeFlagsSum = str(genomeAlignmentDict[key].strangeFlagsSum)

		outputString = str(key) + ',' + readAlignmentCoord + ',' + alignmentQuality + ',' \
		+ alignmentFound + ',' + multipleAlignmentsFound + ',' + strangeFlagsSum + '\n'

		outputHandle.write(outputString)

	outputHandle.close()
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GenerateLogStringForIndividualMfaFile(startTime, endTime, currentFileNumber, totalFiles, \
perfectAlignmentCount, multipleAlignmentCount, lowScoreCount, noAlignCount, strangeFlagsSumCount):
# Updated to include strange flags sum
	duration = (endTime - startTime).total_seconds()
	
	summary = "Reads with Perfect Genome Alignments: " + str(perfectAlignmentCount) + '\n'
	summary += "Reads with Multiple Genome Alignnments: " + str(multipleAlignmentCount) + '\n'
	summary += "Reads with Low Alignment Scores: " + str(lowScoreCount) + '\n'
	summary += "Reads with No Alignment Scores: " + str(noAlignCount) + '\n'
	summary += "Reads with Strange Flags Sum: " + str(strangeFlagsSumCount) + '\n'

	lastOperationDurationStr = "Last operation duration: " + str(duration) + " seconds" + '\n'
		
	timeRemainingEstimate = duration*(totalFiles-currentFileNumber)/60
	timeRemainingStr = "Time remaining estimation: " + str(timeRemainingEstimate) + " minutes"  + '\n'
	
	logStr = summary + lastOperationDurationStr + timeRemainingStr

	return logStr
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateLogFileSummary(totalPerfectAlignmentCount, totalMultipleAlignmentCount, \
totalLowScoreCount, totalNoAlignCount, totalStrangeFlagsSumCount, initialStartTime, completeEndTime):
# Updated to include strange flags sum
	jobDuration = (completeEndTime - initialStartTime).total_seconds()

	summary = '# ----------------------------------------------------------------------------------------------- #\n'
	summary = "Total Reads with Perfect Genome Alignments: " + str(totalPerfectAlignmentCount) + '\n'
	summary += "Total Reads with Multiple Genome Alignnments: " + str(totalMultipleAlignmentCount) + '\n'
	summary += "Total Reads with Low Alignment Scores: " + str(totalLowScoreCount) + '\n'
	summary += "Total Reads with No Alignment Scores: " + str(totalNoAlignCount) + '\n'
	summary += "Total Reads with Strange Flags Sum: " + str(totalStrangeFlagsSumCount) + '\n'
	
	summary += "Complete Job Time: " + str(jobDuration/60) + " minutes\n"
	summary += '# ----------------------------------------------------------------------------------------------- #\n'
	
	return summary
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def UpdateLogFileData(outputLog, outputStr):
	
	logFileHandle = open(outputLog, 'a')
	logFileHandle.write(outputStr)
	logFileHandle.close()
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CountEntriesInSummaryFiles(alignmentFiles, calculateMaxReadIDLength=False):
	import gc
	
	alignmentFileEntries = 0
	
	maxReadIDLength = 10
	
	for file in alignmentFiles:
# 		print("Processing file: " + file)
	
		with open(file) as infile:
			for line in infile:
				alignmentFileEntries += 1
				if calculateMaxReadIDLength == True:
					lineArray = line.strip().split(',')	
					readIDLength = len(lineArray[0])
					if readIDLength > maxReadIDLength:
						maxReadIDLength = readIDLength
		
		gc.collect()
	
	if calculateMaxReadIDLength:
		return [alignmentFileEntries, maxReadIDLength]
	else:
		return alignmentFileEntries
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CountEntriesInAllSummaryFiles(genomeAlignmentFiles, himarRecognitionFiles, indexAlignmentFiles):
	
	print("Counting entries in genome alignment files and calculating maximum length of read ids")
	[genomeAlignmentEntries, maxReadIDLength] = CountEntriesInSummaryFiles(genomeAlignmentFiles, \
	calculateMaxReadIDLength=True)
	print("Counting entries in himar alignment files")
	himarRecognitionEntries = CountEntriesInSummaryFiles(himarRecognitionFiles)
	print("Counting entries in index alignment files")
	indexAlignmentEntries = CountEntriesInSummaryFiles(indexAlignmentFiles)

	return [genomeAlignmentEntries, himarRecognitionEntries, indexAlignmentEntries, maxReadIDLength]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ParseGenomeAlignmentFiles(readIDFieldCode, genomeAlignmentFiles, genomeAlignmentEntries):
# Updated to import flag that indicates if bowtie alignment flags sum was wrong
	import numpy
	import ast
	
	genomeArray = numpy.zeros(genomeAlignmentEntries, \
	dtype={'names':['readID', 'readAlignmentCoord', 'alignmentQuality', 'alignmentFound', \
	'multipleAlignmentsFound', 'himarRecognized', 'index', 'strangeFlagsSum'], \
	'formats':[readIDFieldCode, 'i32', 'i16', 'i8', 'i8', 'i8', 'i8', 'i8']})
	
	i = 0
	for file in genomeAlignmentFiles:
		print("Processing file: " + file)
	
		with open(file) as infile:
			for line in infile:
				lineArray = line.strip().split(',')
		
				genomeArray[i]['readID'] = lineArray[0]
				genomeArray[i]['readAlignmentCoord'] = int(lineArray[1])
				genomeArray[i]['alignmentQuality'] = int(lineArray[2])
			
				if ast.literal_eval(lineArray[3]):
					genomeArray[i]['alignmentFound'] = 1
				else:
					genomeArray[i]['alignmentFound'] = 0
			
				if ast.literal_eval(lineArray[4]):
					genomeArray[i]['multipleAlignmentsFound'] = 1
				else:
					genomeArray[i]['multipleAlignmentsFound'] = 0
					
				if ast.literal_eval(lineArray[5]):
					genomeArray[i]['strangeFlagsSum'] = 1
				else:
					genomeArray[i]['strangeFlagsSum'] = 0
			
				i += 1
	
	return genomeArray
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ParseIndexAlignmentFiles(readIDFieldCode, indexAlignmentFiles, indexAlignmentEntries):
	
	import numpy
	import ast
	
	indexArray = numpy.zeros(indexAlignmentEntries, dtype={'names':['readID', 'index'], \
	'formats':[readIDFieldCode,'i8']})

	
	i = 0
	for file in indexAlignmentFiles:
		print("Processing file: " + file)
	
		with open(file) as infile:
			for line in infile:
				lineArray = line.strip().split(',')	
				indexArray[i]['readID'] = lineArray[0]
			
				# If the index is not recognized, set the index to zero 
				if ast.literal_eval(lineArray[2]) == False:
					index = 0
				else:
					indexArray[i]['index'] = lineArray[1]	
			
				i += 1
		
	return indexArray
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ParseHimarAlignmentFiles(readIDFieldCode, himarRecognitionFiles, himarRecognitionEntries):
	import numpy
	import ast
	
	himarArray = numpy.zeros(himarRecognitionEntries, dtype={'names':['readID', 'himarRecognized'], \
	'formats':[readIDFieldCode,'i8']})

	
	i = 0
	for file in himarRecognitionFiles:
		print("Processing file: " + file)

		with open(file) as infile:
			for line in infile:
				lineArray = line.strip().split(',')		
				readID = lineArray[0]
				himarRecognized = ast.literal_eval(lineArray[1])

				himarArray[i]['readID'] = readID
			
				if himarRecognized:
					himarArray[i]['himarRecognized'] = 1
				else:
					himarArray[i]['himarRecognized'] = 0
			
				i += 1
	return himarArray
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def UpdateGenomeArrayWithHimarAndIndex(genomeArray, indexArray, himarArray):
	
	from pdb import set_trace
	import numpy
	
	sortedGenomeArray = numpy.sort(genomeArray, order='readID')
	sortedHimarArray = numpy.sort(himarArray, order='readID')
	sortedIndexArray = numpy.sort(indexArray, order='readID')
	
# 	set_trace()
	
	i = 0
	while i < len(sortedGenomeArray):
		if sortedGenomeArray[i]['readID'] == sortedHimarArray[i]['readID'] \
		and sortedGenomeArray[i]['readID'] == sortedIndexArray[i]['readID']:
		
			sortedGenomeArray[i]['himarRecognized'] = sortedHimarArray[i]['himarRecognized']
			sortedGenomeArray[i]['index'] = sortedIndexArray[i]['index']
		else:
			print("Read ids don't match")
		i += 1

	return sortedGenomeArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def WriteCompiledGenomeArray(outputFile, sortedGenomeArray):
# Updated to print out strangeFlagsSum flag

	fileHandle = open(outputFile, 'w')
	
	i = 0
	while i < len(sortedGenomeArray):
		
		outputStr = str(sortedGenomeArray[i]['readID'].tostring().decode('utf-8')) + ','
		outputStr += str(sortedGenomeArray[i]['readAlignmentCoord']) + ','
		outputStr += str(sortedGenomeArray[i]['alignmentQuality']) + ','
		outputStr += str(sortedGenomeArray[i]['alignmentFound']) + ','
		outputStr += str(sortedGenomeArray[i]['multipleAlignmentsFound']) + ','
		outputStr += str(sortedGenomeArray[i]['himarRecognized']) + ','
		outputStr += str(sortedGenomeArray[i]['index']) + ','
		outputStr += str(sortedGenomeArray[i]['strangeFlagsSum']) + '\n'
		
		fileHandle.write(outputStr)
		
		i += 1
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateSudokuReadTaxonomy(genomeArray, outputLog, operationTimingInterval=100000):
# Updated to include strange flags flag in read taxonomy	
	
	import datetime
	
	logStr = gSeparatorString
	logStr += 'Generating Sudoku Read Taxonomy'
	
	print(logStr)
	
	UpdateLogFileData(outputLog, logStr)
	
	taxDict = {'totalReads':0, 'totalReadsWithHimar':0,	'totalReadsWithIndex':0,\
	'totalReadsWithNoStrangeFlags':0, 'totalReadsWithNoHimar':0, 'totalReadsWithNoIndex':0,\
	'totalReadsWithStrangeFlags':0, 'totalReadsWithPerfectGenome':0, \
	'totalReadsWithMultipleAlignment':0, 'totalReadsWithNoAlignment':0,	\
	'totalReadsWithLowScoreAlignment':0, \
	'totalReadsWithIndexAndHimar':0, 'totalReadsWithIndexAndHimarAndPerfectGenome':0, \
	'totalReadsWithIndexAndHimarAndMultipleAlignment':0, \
	'totalReadsWithIndexAndHimarButNoAlignment':0,	\
	'totalReadsWithIndexAndHimarAndLowScoreAlignment':0,\
	'totalReadsWithIndexAndHimarAndPerfectGenomeAndNoStrangeFlags':0,	
	'totalReadsWithIndexAndNoHimar':0, 'totalReadsWithIndexNoHimarAndPerfectGenome':0,\
	'totalReadsWithIndexNoHimarAndMultipleAlignment':0,	\
	'totalReadsWithIndexNoHimarButNoAlignment':0, \
	'totalReadsWithIndexNoHimarAndLowScoreAlignment':0, 'totalReadsWithNoIndexAndHimar':0, \
	'totalReadsWithNoIndexAndHimarAndPerfectGenome':0,\
	'totalReadsWithNoIndexAndHimarAndMultipleAlignment':0,\
	'totalReadsWithNoIndexAndHimarButNoAlignment':0,\
	'totalReadsWithNoIndexAndHimarAndLowScoreAlignment':0, 'totalReadsWithNoIndexAndNoHimar':0,\
	'totalReadsWithNoIndexNoHimarAndPerfectGenome':0,\
	'totalReadsWithNoIndexNoHimarAndMultipleAlignment':0,\
	'totalReadsWithNoIndexNoHimarButNoAlignment':0,	\
	'totalReadsWithNoIndexNoHimarAndLowScoreAlignment':0}
	

	i = 0
	
	startTime = datetime.datetime.now()

	while i < len(genomeArray):
	
		sudokuRead = genomeArray[i]
		
		taxDict['totalReads'] += 1
	
		if sudokuRead['himarRecognized'] == 1:
			taxDict['totalReadsWithHimar'] += 1
		else:
			taxDict['totalReadsWithNoHimar'] += 1
		
		if sudokuRead['index'] != 0:
			taxDict['totalReadsWithIndex'] += 1
		else:	
			taxDict['totalReadsWithNoIndex'] += 1
			
		if sudokuRead['strangeFlagsSum'] == 0:
			taxDict['totalReadsWithNoStrangeFlags'] += 1
		else:
			taxDict['totalReadsWithStrangeFlags'] += 1
		
		if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
			taxDict['totalReadsWithPerfectGenome'] += 1
		elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] and sudokuRead['alignmentQuality'] > 1:
			taxDict['totalReadsWithMultipleAlignment'] += 1
		elif sudokuRead['alignmentFound'] == 0:
			taxDict['totalReadsWithNoAlignment'] += 1
		elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
			taxDict['totalReadsWithLowScoreAlignment'] += 1
	
	
		if sudokuRead['himarRecognized'] == 1 and sudokuRead['index'] != 0:
			taxDict['totalReadsWithIndexAndHimar'] += 1
			if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
				taxDict['totalReadsWithIndexAndHimarAndPerfectGenome'] += 1
				if sudokuRead['strangeFlagsSum'] == 0:
					taxDict['totalReadsWithIndexAndHimarAndPerfectGenomeAndNoStrangeFlags'] += 1
			elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 1 and sudokuRead['alignmentQuality'] > 1:
				taxDict['totalReadsWithIndexAndHimarAndMultipleAlignment'] += 1
			elif sudokuRead['alignmentFound'] == 0:
				taxDict['totalReadsWithIndexAndHimarButNoAlignment'] += 1
			elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
				taxDict['totalReadsWithIndexAndHimarAndLowScoreAlignment'] += 1
	
	
		elif sudokuRead['himarRecognized'] == 0 and sudokuRead['index'] != 0:
			taxDict['totalReadsWithIndexAndNoHimar'] += 1
			if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
				taxDict['totalReadsWithIndexNoHimarAndPerfectGenome'] += 1
			elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 1 and sudokuRead['alignmentQuality'] > 1:
				taxDict['totalReadsWithIndexNoHimarAndMultipleAlignment'] += 1
			elif sudokuRead['alignmentFound'] == 0:
				taxDict['totalReadsWithIndexNoHimarButNoAlignment'] += 1
			elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
				taxDict['totalReadsWithIndexNoHimarAndLowScoreAlignment'] += 1
	
	
		elif sudokuRead['himarRecognized'] and sudokuRead['index'] == 0:
			taxDict['totalReadsWithNoIndexAndHimar'] += 1
			if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
				taxDict['totalReadsWithNoIndexAndHimarAndPerfectGenome'] += 1
			elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 1 and sudokuRead['alignmentQuality'] > 1:
				taxDict['totalReadsWithNoIndexAndHimarAndMultipleAlignment'] += 1
			elif sudokuRead['alignmentFound'] == 0:
				taxDict['totalReadsWithNoIndexAndHimarButNoAlignment'] += 1
			elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
				taxDict['totalReadsWithNoIndexAndHimarAndLowScoreAlignment'] += 1


		elif sudokuRead['himarRecognized'] == 0 and sudokuRead['index'] == 0:
			taxDict['totalReadsWithNoIndexAndNoHimar'] += 1
			if sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 0 and sudokuRead['alignmentQuality'] > 1:
				taxDict['totalReadsWithNoIndexNoHimarAndPerfectGenome'] += 1
			elif sudokuRead['alignmentFound'] == 1 and sudokuRead['multipleAlignmentsFound'] == 1 and sudokuRead['alignmentQuality'] > 1:
				taxDict['totalReadsWithNoIndexNoHimarAndMultipleAlignment'] += 1
			elif sudokuRead['alignmentFound'] == 0:
				taxDict['totalReadsWithNoIndexNoHimarButNoAlignment'] += 1
			elif sudokuRead['alignmentQuality'] <= 1 and sudokuRead['alignmentFound'] == 1:
				taxDict['totalReadsWithNoIndexNoHimarAndLowScoreAlignment'] += 1
	
		if i % operationTimingInterval == 0:
			endTime = datetime.datetime.now()
			duration = (endTime - startTime).total_seconds()
			timeRemainingEstimate = (duration/operationTimingInterval)*(len(genomeArray)-i)/60
			
			logStr = "Time for last " + str(operationTimingInterval) + " operations: " \
			+ str(duration) + " seconds\n"
			logStr += "Time remaining estimation: " + str(timeRemainingEstimate) + " minutes\n"
			
			print(logStr)
			UpdateLogFileData(outputLog, logStr)
			
			startTime = datetime.datetime.now()
	
		i += 1
		
	return taxDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def OutputReadTaxonomy(taxonomyDict, outputLog):
	
	outputLogHandle  = open(outputLog, 'w')
	
	taxonomyDictKeys = taxonomyDict.keys()

	outputStr = ''
	
	for key in taxonomyDictKeys:
		outputStr += key + ': ' + str(taxonomyDict[key]) + '\n'
	
	print(outputStr)
	
	UpdateLogFileData(outputLog, outputStr)
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportGenomeCompilationFile(genomeCompilationFile, includesb=False):
	import gc
	import numpy
	import pdb
	
	readIDRe = re.compile('\d+:\d+:\d+:\d+')
	
	# Determine length of genomeCompilationFile
	
	fHandle = open(genomeCompilationFile)
	compilationEntries = 0
	maxReadIDLength = 15
	
	with open(genomeCompilationFile) as infile:
		for line in infile:
			compilationEntries += 1
			lineArray = line.strip().split(',')	
			readIDLength = len(lineArray[0])
			if readIDLength > maxReadIDLength:
				maxReadIDLength = readIDLength
	
	gc.collect()
# 	pdb.set_trace()
	
	readIDFieldCode = 'a' + str(maxReadIDLength+2)
	
	genomeArray = numpy.zeros(compilationEntries, \
	dtype={'names':['readID', 'readAlignmentCoord', 'alignmentQuality', 'alignmentFound', \
	'multipleAlignmentsFound', 'himarRecognized', 'index', 'strangeFlagsSum'], \
	'formats':[readIDFieldCode, 'i32', 'i16', 'i8', 'i8', 'i8', 'i8', 'i8']})
	
	i = 0
	with open(genomeCompilationFile) as infile:
		for line in infile:
			lineArray = line.strip().split(',')
			if includesb:
				genomeArray[i]['readID'] = readIDRe.search(lineArray[0]).group()
			else:
				genomeArray[i]['readID'] = lineArray[0]
			genomeArray[i]['readAlignmentCoord'] = int(lineArray[1])
			genomeArray[i]['alignmentQuality'] = int(lineArray[2])
			genomeArray[i]['alignmentFound'] = int(lineArray[3])
			genomeArray[i]['multipleAlignmentsFound'] = int(lineArray[4])
			genomeArray[i]['himarRecognized'] = int(lineArray[5])
			genomeArray[i]['index'] = int(lineArray[6])
			genomeArray[i]['strangeFlagsSum'] = int(lineArray[7])
			

			i += 1
	
	return genomeArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def GenerateBarcodeLookupTable(barcodeFile):
	barcodeFileHandle = open(barcodeFile, 'r')
	barcodeFileData = barcodeFileHandle.readlines()
	barcodeFileHandle.close()

	barcodeLookupTable = {}

	for line in barcodeFileData:
		if line[0] != '#':
			lineData = line.strip().split(',')
			revCompl = lineData[2]
			barcodeNumber = lineData[3]
			barcodeLookupTable[revCompl] = int(barcodeNumber)
		
	barcodes = barcodeLookupTable.keys()
	return [barcodeLookupTable, barcodes]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ParseGenomeArrayToPoolFiles(genomeArray, barcodeLookupTable, poolFileBaseDir, poolFilePrefix):	
	import re
	import gc
	
	# Generate pool file handles
	poolFileHandleDict = {}
	poolFiles = []
	
	indices = barcodeLookupTable.values()
	
	for index in indices:
		poolFileName = poolFileBaseDir + poolFilePrefix + str(index) + '.csv'
		poolFileHandleDict[index] = open(poolFileName, 'w')
		poolFiles.append(poolFileName)
	
	# Run through genome array and parse good reads into the pool files
	i = 0
	while i < len(genomeArray):
		read = genomeArray[i]
		readID = genomeArray[i]['readID'].tostring().decode('utf-8')
		readAlignmentCoord = genomeArray[i]['readAlignmentCoord']
		alignmentQuality = genomeArray[i]['alignmentQuality']
		alignmentFound = genomeArray[i]['alignmentFound']
		multipleAlignmentsFound= genomeArray[i]['multipleAlignmentsFound']
		himarRecognized = genomeArray[i]['himarRecognized']
		index = genomeArray[i]['index']
		strangeFlagsSum = genomeArray[i]['strangeFlagsSum']
		
		if alignmentQuality > 1 and multipleAlignmentsFound == 0 and himarRecognized == 1 \
		and index != 0 and strangeFlagsSum == 0:
			outputStr = str(readID) + ',' + str(readAlignmentCoord) + '\n'
			poolFileHandleDict[index].write(outputStr)
		i += 1
	
	# Close up all of the pool file handles
	for index in indices:
		poolFileHandleDict[index].close()
	
	return poolFiles
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateFileNamesFromFileListAndBaseDir(baseDir, files):
	
	fullFileNames = []
	
	for file in files:
		fullFileNames.append(baseDir + '/' + file)
		
	return fullFileNames
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateValidCoordsList(genomeArray):
	
	import numpy
	
	validCoords = 0
	
	i = 0
	maxReadIDLength = 15
	
	while i < len(genomeArray):
		
		coord = genomeArray[i]
		
		if coord['alignmentQuality'] > 1 and coord['alignmentFound'] == 1 \
		and coord['multipleAlignmentsFound'] == 0 and coord['himarRecognized'] == 1 \
		and coord['index'] > 0:
			validCoords += 1
			readIDLength = len(coord[0])
			if readIDLength > maxReadIDLength:
				maxReadIDLength = readIDLength
			
		i += 1
	
	readIDFieldCode = 'a' + str(maxReadIDLength+2)
	
	validGenomeArray = numpy.zeros(validCoords, \
	dtype={'names':['readID', 'readAlignmentCoord', 'alignmentQuality', 'index'], \
	'formats':[readIDFieldCode, 'i32', 'i16', 'i8', 'i8']})
	
	i = 0
	j = 0
	
	while i < len(genomeArray) and j < validCoords:
		
		coord = genomeArray[i]
		
		if coord['alignmentQuality'] > 1 and coord['alignmentFound'] == 1 \
		and coord['multipleAlignmentsFound'] == 0 and coord['himarRecognized'] == 1 \
		and coord['index'] > 0 and coord['strangeFlagsSum'] == 0:
			validGenomeArray[j]['readID'] = coord['readID']
			validGenomeArray[j]['readAlignmentCoord'] = coord['readAlignmentCoord']
			validGenomeArray[j]['alignmentQuality'] = coord['alignmentQuality']
			validGenomeArray[j]['index'] = coord['index']
			j += 1
		
		i += 1
	
	return validGenomeArray
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportPoolFiles(barcodeFile, poolFiles, onlyExportUniqueCoords=True):
	
	import os
	import gc
	from scipy import unique
	
	# Count lines in each pool file
	
	lineCountDict = {}
	maxReadIDLength = 15
	
	for poolFile in poolFiles:
		index = GetIndexCodeFromPoolFileName(poolFile)
		lineCount = 0
		with open(poolFile) as infile:
			for line in infile:
				lineCount += 1
				lineArray = line.strip().split(',')
				readIDLength = len(lineArray[0])
				if readIDLength > maxReadIDLength:
					maxReadIDLength = readIDLength
		
		lineCountDict[index] = [lineCount, poolFile]
		gc.collect()
	
	# Import pools
	
	pooledCoordArrayDict = {}
	lineCountKeys = lineCountDict.keys()
	
	readIDFieldCode = 'a' + str(maxReadIDLength+2)
	
	for key in lineCountKeys:
		lineCountArray = lineCountDict[key]
		lineCount = lineCountArray[0]
		fileName = lineCountArray[1]
		
		if onlyExportUniqueCoords == True:		
			coordArray = numpy.zeros(lineCount, dtype='i32')
			i = 0
			with open(fileName) as infile:
				for line in infile:
					coordArray[i] = int(line.strip().split(',')[1])
					i += 1
			
			pooledCoordArrayDict[key] = unique(coordArray)
			
		elif onlyExportUniqueCoords == False:
			coordArray = numpy.zeros(lineCount, dtype={'names':['readID', 'readAlignmentCoord'], \
			'formats':[readIDFieldCode, 'i32']})
			with open(fileName) as infile:
				i = 0
				for line in infile:
					lineArray = line.strip().split(',')
					coordArray[i]['readAlignmentCoord'] = int(lineArray[1])
					coordArray[i]['readID'] = str(lineArray[0])
					i += 1
					
				pooledCoordArrayDict[key] = coordArray

	return pooledCoordArrayDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GetIndexCodeFromPoolFileName(poolFileName, poolFileNamePrefix='pool_'):
	
	import re
	
	import pdb
	
	poolNameRe = re.compile(poolFileNamePrefix + '(\d+)')

# 	pdb.set_trace()
	[fileName, fileExtension] = ProcessFileNameExtension(poolFileName)
	
	poolNameMatch = poolNameRe.match(fileName)
	
	if poolNameRe.match(fileName) != None:
		indexCode = poolNameMatch.group(1)
	else:
		print("Error. Pool file name does not match regular expression: " + poolFileName)
	
	return indexCode
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GeneratePoolNameToPoolCodeLookupTable(barcodeFile):
# Note that this
	barcodeFileHandle = open(barcodeFile, 'r')
	barcodeFileData = barcodeFileHandle.readlines()
	barcodeFileHandle.close()

	indexLookupTable = {}

	for line in barcodeFileData:
		if line[0] != '#':
			lineData = line.strip().split(',')
			poolName = lineData[0]
			forwardSeq = lineData[1]
			revCompl = lineData[2]
			barcodeNumber = lineData[3]
			indexLookupTable[str(barcodeNumber)] = poolName
		
	
	return indexLookupTable
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateDTypeDictForPoolPresenceDict(indexLookupTable, numberPossibleAddresses=10):
	
	import pdb
	
	poolNames = indexLookupTable.values()
	poolNames = sorted(poolNames)
	
	dtypeDict = {}
	dtypeDict['names'] = []
	dtypeDict['formats'] = []
	
	# Part of string for genomic coordinate
	dtypeDict['names'].append('readAlignmentCoord')
	dtypeDict['formats'].append('i32')
	
	# Add in columns for pools
	i = 0
	while i < len(poolNames):
		dtypeDict['names'].append(poolNames[i])
		dtypeDict['formats'].append('i32')
		i += 1
	
	dtypeDict['names'].append('addresses')
	dtypeDict['formats'].append('(' + str(numberPossibleAddresses) + ',4)int32')
	
	return dtypeDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GeneratePoolPresenceTable(uniqueCoords, sortedValidGenomeArray, indexLookupTable, dtypeDict):
	import numpy
	import pdb
	
	poolPresenceTable = numpy.zeros(len(uniqueCoords), dtype=dtypeDict)
	
	i = 0
	while i < len(uniqueCoords):
		poolPresenceTable[i] = uniqueCoords[i]
		i += 1
	
	
	i = 0
	j = 0
	
	while j < len(uniqueCoords):
		while(i < len(sortedValidGenomeArray) and \
		sortedValidGenomeArray[i]['readAlignmentCoord'] == poolPresenceTable[j]['readAlignmentCoord']):
		
			index = str(sortedValidGenomeArray[i]['index'])
			column = indexLookupTable[index]
			poolPresenceTable[j][column] += 1
		
			i += 1
		
		j += 1
		
	return poolPresenceTable
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def WritePoolPresenceTable(filename, poolPresenceTable):

	fhandle = open(filename, 'w')
	
	poolColumns = [\
	'readAlignmentCoord',\
	'A','B','C','D','E','F','G','H', \
	'1','2','3','4','5','6','7','8','9','10','11','12',\
	'PR01','PR02','PR03','PR04','PR05','PR06','PR07','PR08','PR09','PR10',\
	'PR11','PR12','PR13','PR14','PR15','PR16','PR17','PR18','PR19','PR20',\
	'PC01','PC02','PC03','PC04','PC05','PC06','PC07','PC08','PC09','PC10','PC11','PC12','PC13',\
	'PC14','PC15','PC16','PC17','PC18','PC19','PC20','PC21',\
	'SoG','pta','Blank']
	
	addressColumns = [\
	'address1_r','address1_c','address1_pr','address1_pc',\
	'address2_r','address2_c','address2_pr','address2_pc', \
	'address3_r','address3_c','address3_pr','address3_pc', \
	'address4_r','address4_c','address4_pr','address4_pc', \
	'address5_r','address5_c','address5_pr','address5_pc', \
	'address6_r','address6_c','address6_pr','address6_pc',\
	'address7_r','address7_c','address7_pr','address7_pc',\
	'address8_r','address8_c','address8_pr','address8_pc',\
	'address9_r','address9_c','address9_pr','address9_pc',\
	'address10_r','address10_c','address10_pr','address10_pc']
	
		
	# Write out the header line
	
	headerLine = ''
	
	i = 0
	while i < len(poolColumns):
		headerLine += poolColumns[i] + ','
		i += 1
	
	i = 0
	while i < len(addressColumns):
		headerLine += addressColumns[i]
		if i < len(addressColumns) -1:
			headerLine += ','
		i += 1
	headerLine += '\n'
	
	
	
	totalStr = ''
	
	i = 0
	while i < poolPresenceTable.shape[0]:
		outputStr = ''
		j = 0
		while j < len(poolColumns):
			outputStr += str(poolPresenceTable[i][poolColumns[j]]) + ','
			j += 1
		
		j = 0
		addresses = poolPresenceTable[i]['addresses']
		while j < len(addresses):
			k = 0
			address = addresses[j]
			while k < len(address):
				outputStr += str(address[k])
				if j == (len(addresses) - 1) and k == (len(address) - 1):
					outputStr += ''
				else:
					outputStr += ','
				k += 1
			j += 1
# 		print(outputStr)
		outputStr += '\n'
		totalStr += outputStr
		i += 1
		
	writeStr = headerLine + totalStr
	
	fhandle.write(writeStr)
	fhandle.close()

	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def WritePoolPresenceTable2(filename, poolPresenceTable):

	fhandle = open(filename, 'w')
	
	poolColumns = [\
	'readAlignmentCoord',\
	'A','B','C','D','E','F','G','H', \
	'1','2','3','4','5','6','7','8','9','10','11','12',\
	'PR01','PR02','PR03','PR04','PR05','PR06','PR07','PR08','PR09','PR10',\
	'PR11','PR12','PR13','PR14','PR15','PR16','PR17','PR18','PR19','PR20',\
	'PC01','PC02','PC03','PC04','PC05','PC06','PC07','PC08','PC09','PC10','PC11','PC12','PC13',\
	'PC14','PC15','PC16','PC17','PC18','PC19','PC20','PC21',\
	'SoG','pta','Blank']
		
	# Write out the header line
	
	headerLine = ''
	
	i = 0
	while i < len(poolColumns):
		headerLine += poolColumns[i] + ','
		i += 1
	
	headerLine += '\n'
	
	totalStr = ''
	
	i = 0
	while i < poolPresenceTable.shape[0]:
		outputStr = ''
		j = 0
		while j < len(poolColumns):
			outputStr += str(poolPresenceTable[i][poolColumns[j]]) + ','
			j += 1
		
		outputStr += '\n'
		totalStr += outputStr
		i += 1
		
	writeStr = headerLine + totalStr
	
	fhandle.write(writeStr)
	fhandle.close()

	return
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def ProcessFileNameExtension(fileName):
	import re
	import os
	import pdb
	
# 	pdb.set_trace()
	
	baseName = os.path.basename(fileName)
	[fileName, fileExtension] = os.path.splitext(baseName)
		
	return [fileName, fileExtension]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateUniqueCoordsList(genomeArray):
	from scipy import unique
	
	coords = genomeArray['readAlignmentCoord']
	
	uniqueCoords = unique(coords)
	
	return uniqueCoords
# ------------------------------------------------------------------------------------------------ #






####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Functions for handling pool presence table analysis
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def ImportPoolPresenceTable(poolPresenceTableFileName, poolColumns):
	# Import the pool presence table
	import numpy
	import gc
	from pdb import set_trace
	
	
	fHandle = open(poolPresenceTableFileName, 'r')
	data = fHandle.readlines()
	
	headers = data[0].strip().split(',')
	
	dtypeDict = GenerateDTypeDictForPoolPresenceTable(poolColumns)
		
	poolColumnToHeaderIndexDict = {}
	for col in poolColumns:
		poolColumnToHeaderIndexDict[col] = headers.index(col)
		
	# Initialize a numpy array for the pool presence table
	poolPresenceTable = numpy.zeros(len(data)-1, dtype=dtypeDict)
	
	# Populate the poolPresenceTable
	i = 1
	while i < len(data):
		lineData = data[i].strip().split(',')
		for col in poolColumns:
			poolPresenceTable[i-1][col] = lineData[poolColumnToHeaderIndexDict[col]]
		i += 1
		
		
	return poolPresenceTable
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateDTypeDictForPoolPresenceTable(poolColumns):

	dtypeDict = {}
	dtypeDict['names'] = []
	dtypeDict['formats'] = []
	
	for col in poolColumns:
		if col == 'readAlignmentCoord':
			dtypeDict['names'].append(col)
			dtypeDict['formats'].append('i32')
		else:
			dtypeDict['names'].append(col)
			dtypeDict['formats'].append('i8')

	return dtypeDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def GenerateBlankPhysicalAddressDict(poolPresenceTable):

	uniqueCoords = poolPresenceTable['readAlignmentCoord']

	physicalAddressDict = {}

	i = 0

	while i < len(uniqueCoords):
		physicalAddressDict[int(uniqueCoords[i])] = []
		i += 1

	return physicalAddressDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def FindAddressCoords(poolPresenceTableLine, axisPools, threshold=5):
	
	axisAddresses = []
	
	i = 0
	while i < len(axisPools):
		if poolPresenceTableLine[axisPools[i]] >= threshold:
			axisAddresses.append([axisPools[i], poolPresenceTableLine[axisPools[i]]])
		i += 1
	

	return axisAddresses
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindAddressCoords2(poolPresenceTableLine, axisPools, threshold=5):
	# Very much the same as FindAddressCoords, but returns a list of just axis addresses and a 
	# dict of the read counts for each of these axis addresses, rather than combining this dict
	# into the axisAddresses array
	
	axisAddresses = []
	axisAddressScoreDict = {}
	
	i = 0
	while i < len(axisPools):
		if poolPresenceTableLine[axisPools[i]] >= threshold:
			axisAddresses.append(axisPools[i])
			axisAddressScoreDict[axisPools[i]] = poolPresenceTableLine[axisPools[i]]
		i += 1
	

	return [axisAddresses, axisAddressScoreDict]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculatePossibleAddresses(addresses_r, addresses_c, addresses_pr, addresses_pc):
# Calculate the unambiguous library addresses that that a read alignment coord maps to
# The score is the total number of reads that are associated with the library address assignment
	possibleAddresses = []
	
	for address_r in addresses_r:
		for address_c in addresses_c:
			for address_pr in addresses_pr:
				for address_pc in addresses_pc:
					
					row = address_r[0]
					row_score = address_r[1]
					col = address_c[0]
					col_score = address_c[1]
					pr = address_pr[0]
					pr_score = address_pr[1]
					pc = address_pc[0]
					pc_score = address_pc[1]
					
					possibleAddress = row + '_' + col + '_' + pr + '_' + pc
					possibleAddressScore = row_score + col_score + pr_score + pc_score
					
					possibleAddresses.append([possibleAddress, possibleAddressScore])

	return possibleAddresses
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateLibraryAddressesForPoolPresenceTableLine2(poolPresenceTableLine, rowPools, colPools, \
prPools, pcPools, logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict, \
threshold):

# Used in generation of pool presence table taxonomy
# Used to calculate possible addresses for pool presence table line
	
	coord = int(poolPresenceTableLine['readAlignmentCoord'])
	
	addresses_r = FindAddressCoords(poolPresenceTableLine, rowPools, threshold=threshold)
	addresses_c = FindAddressCoords(poolPresenceTableLine, colPools, threshold=threshold)
	addresses_pr = FindAddressCoords(poolPresenceTableLine, prPools, threshold=threshold)
	addresses_pc = FindAddressCoords(poolPresenceTableLine, pcPools, threshold=threshold)
	
	possibleAddressesAndScores = \
	CalculatePossibleAddresses2(addresses_r, addresses_c, addresses_pr, addresses_pc, \
	logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict)
	
	return possibleAddressesAndScores
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Functions for calculating a Voigt function to a read ratio histogram
def voigtFunction(x, p):
	
	from scipy.special import erfc
	from numpy import exp
	from numpy import sqrt
	from numpy import pi, float64

	a = p[0]
	c = p[1]
	delta = p[2]
	sigma = p[3]
	
	firstArg = ((-1.j)*(-c + x) + delta)/(sqrt(2)*sigma)
	secondArg = ((1.j)*(-c + x) + delta)/(sqrt(2)*sigma)

	voigtEquation = a\
	*(exp(firstArg**2)*erfc(firstArg) \
	+ exp(secondArg**2)*erfc(secondArg) ) \
	/ (2*sqrt(2*pi)*sigma)
	
	voigtEquation = float64(voigtEquation)
	
	return voigtEquation

	
def voigtResiduals(p, y, x): 
	err = y - voigtFunction(x,p) 
	return err

def voigtFit(x,y, p0):
	import scipy
	from scipy.optimize import leastsq
	
	plsq = leastsq(voigtResiduals, p0, args=(y, x), maxfev=2000)
	
	return [plsq[0], voigtFunction(x, plsq[0])]

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ConvertBinEdgesToCenters(binEdges):
	
	from numpy import zeros
	
	binCenters = zeros(len(binEdges)-1)
	
	i = 0
	
	while i < len(binEdges) - 1:
		
		binCenter = binEdges[i] + (binEdges[i + 1] - binEdges[i])
		
		binCenters[i] = binCenter
		i += 1
		
	
	return binCenters
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateVoigtScore(logNRatio, plsq, normalizationFactor):
	
	voigtScore = voigtFunction(logNRatio, plsq)/normalizationFactor

	return voigtScore
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def CalculatePossibleAddresses2(addresses_r, addresses_c, addresses_pr, addresses_pc, plsqDict, \
normalizationFactorDict):
# Calculate the unambiguous library addresses that that a read alignment coord maps to
# Same as CalculatePossibleAddresses
# However, the reported score is slightly different. It is a collection of the ratios of pool axes
# entries. 

# Bookmark: this is where I'll add the Voigt function
	
	import numpy
	import pdb

	possibleAddressesAndScores = []
	
	for address_r in addresses_r:
		for address_c in addresses_c:
			for address_pr in addresses_pr:
				for address_pc in addresses_pc:
					
					row = address_r[0]
					col = address_c[0]
					pr = address_pr[0]
					pc = address_pc[0]
					
					possibleAddress = row + '_' + col + '_' + pr + '_' + pc
					
					nRowReads = address_r[1]
					nColReads = address_c[1]
					nPRReads = address_pr[1]
					nPCReads = address_pc[1]
					
					totalReads = nRowReads + nColReads + nPRReads + nPCReads
					
					nr2nc = nRowReads/nColReads
					nr2npr = nRowReads/nPRReads
					nr2npc = nRowReads/nPCReads
					nc2npr = nColReads/nPRReads
					nc2npc = nColReads/nPCReads
					npr2npc = nPRReads/nPCReads
					
					logNr2nc = numpy.log(nr2nc)
					logNr2npr = numpy.log(nr2npr)
					logNr2npc = numpy.log(nr2npc)
					logNc2npr = numpy.log(nc2npr)
					logNc2npc = numpy.log(nc2npc)
					logNpr2npc = numpy.log(npr2npc)
										
					voigtScoreNr2nc = CalculateVoigtScore(logNr2nc, plsqDict['nr2nc'], \
					normalizationFactorDict['nr2nc'])
					
					voigtScoreNr2npr = CalculateVoigtScore(logNr2npr, plsqDict['nr2npr'], \
					normalizationFactorDict['nr2npr'])
					
					voigtScoreNr2npc = CalculateVoigtScore(logNr2npc, plsqDict['nr2npc'], \
					normalizationFactorDict['nr2npc'])
					
					voigtScoreNc2npr = CalculateVoigtScore(logNc2npr, plsqDict['nc2npr'], \
					normalizationFactorDict['nc2npr'])
					
					voigtScoreNc2npc = CalculateVoigtScore(logNc2npc, plsqDict['nc2npc'], \
					normalizationFactorDict['nc2npc'])
					
					voigtScoreNpr2npc = CalculateVoigtScore(logNpr2npc, plsqDict['npr2npc'], \
					normalizationFactorDict['npr2npc'])

					scoreDict = {'nr2nc':voigtScoreNr2nc, 'nr2npr':voigtScoreNr2npr, \
					'nr2npc':voigtScoreNr2npc, 'nc2npr':voigtScoreNc2npr, \
					'nc2npc':voigtScoreNc2npc, 'npr2npc':voigtScoreNpr2npc}
					
					score = voigtScoreNr2nc * voigtScoreNr2npr * voigtScoreNr2npc \
					* voigtScoreNc2npr * voigtScoreNc2npc * voigtScoreNpr2npc

					possibleAddressesAndScores.append([possibleAddress, scoreDict, score, \
					totalReads])

	return possibleAddressesAndScores
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FitReadNumberRatiosAndPlot(logNr2ncArray, logBins, xAxisLabel, plsq0, integrate=False):

	import numpy
	import matplotlib.pyplot as plt
	import pdb
	from scipy.integrate import simps

	valueslogNr2nc, baselogNr2nc = numpy.histogram(logNr2ncArray, bins=logBins)
	binCentersLogNr2nc = ConvertBinEdgesToCenters(baselogNr2nc)
	
	

	[plsqlogNr2nc, voigtlogNr2nc] = voigtFit(binCentersLogNr2nc, valueslogNr2nc, plsq0)
	
	plt.figure()
	
	
	plt.plot(binCentersLogNr2nc, voigtlogNr2nc, linewidth=1)
	plt.scatter(binCentersLogNr2nc, valueslogNr2nc)
	
	plt.xlabel(xAxisLabel)
	plt.ylabel("Number of Pool Presence Table Lines")
	
	plt.grid()
	
	plt.xlim(min(logBins),max(logBins))
	
# 	pdb.set_trace()
	
	plt.ylim(-100,2500)

	plt.show()
	
	if integrate == True:
		integral = simps(voigtlogNr2nc, binCentersLogNr2nc)
		return [valueslogNr2nc, baselogNr2nc, integral, plsqlogNr2nc, \
		voigtlogNr2nc, binCentersLogNr2nc]
	else:
		return [valueslogNr2nc, baselogNr2nc] 
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def CalculatePhysicalAddresses(poolPresenceTable, rowPools, colPools, prPools, pcPools, \
threshold=5):
	
	physicalAddressDict = GenerateBlankPhysicalAddressDict(poolPresenceTable)
	
	i = 0
	while i < len(poolPresenceTable):
		entry = poolPresenceTable[i]
		coord = int(poolPresenceTable[i]['readAlignmentCoord'])
		
		addresses_r = FindAddressCoords(poolPresenceTable[i], rowPools, threshold=threshold)
		addresses_c = FindAddressCoords(poolPresenceTable[i], colPools, threshold=threshold)
		addresses_pr = FindAddressCoords(poolPresenceTable[i], prPools, threshold=threshold)
		addresses_pc = FindAddressCoords(poolPresenceTable[i], pcPools, threshold=threshold)
		
		possibleAddresses = CalculatePossibleAddresses(addresses_r, addresses_c, addresses_pr, \
		addresses_pc)
		
		physicalAddressDict[coord] = possibleAddresses
		
		i += 1
	
	return physicalAddressDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculatePhysicalAddressCoordsForLine(poolPresenceTableLine, rowPools, colPools, prPools, \
pcPools, threshold=5):

	addresses_r = FindAddressCoords2(poolPresenceTableLine, rowPools, threshold=threshold)
	addresses_c = FindAddressCoords2(poolPresenceTableLine, colPools, threshold=threshold)
	addresses_pr = FindAddressCoords2(poolPresenceTableLine, prPools, threshold=threshold)
	addresses_pc = FindAddressCoords2(poolPresenceTableLine, pcPools, threshold=threshold)
	
	# Remember, each line in the addresses array is a 2 element list, the first containing the pool
	# name, and the second containing the number of reads associated with it. 
	
	return [addresses_r, addresses_c, addresses_pr, addresses_pc]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def GroupPoolPresenceTable(poolPresenceTable, rowPools, colPools, prPools, pcPools):
	
	run = []
	result = [run]
	expect = None
	CurrentLineMatchesPrevious = True

	i = 0
	while i < len(poolPresenceTable):
		if CurrentLineMatchesPrevious:
			run.append(poolPresenceTable[i])
		else:
			run = [poolPresenceTable[i]]
			result.append(run)
		
		CurrentLineMatchesPrevious = False
		
		if i < len(poolPresenceTable) - 1:
			CurrentLineMatchesPrevious = \
			DoesNextLineGroupWithCurrentLine(poolPresenceTable[i+1], poolPresenceTable[i], \
			rowPools, colPools, prPools, pcPools)
			# print(CurrentLineMatchesPrevious)

		i += 1
		
    # Assign likely real genomic coordinates to the groups of read alignment coordinates
    # Do something with AT and TA locations in genome 
    
	return result
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def DoesNextLineGroupWithCurrentLine(nextLine, currentLine, \
rowPools, colPools, prPools, pcPools, maxGap=1):
	
	# First, test if the two read alignment coordinates are adjacent to one another
	currentCoord = currentLine['readAlignmentCoord']
	expect = currentCoord + maxGap
	nextCoord = nextLine['readAlignmentCoord']
	
	nextLineMatchesWithCurrent = False
	
	if nextCoord <= expect:
		# print("coords are consecutive to within gap")
		poolPresenceProfileSimilar = \
		ArePoolPresenceProfilesSimilar(nextLine, currentLine, rowPools, colPools, prPools, pcPools,\
		readUseThreshold=1)
		if poolPresenceProfileSimilar:
			nextLineMatchesWithCurrent = True
			
	return nextLineMatchesWithCurrent
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def GroupPoolPresenceTable2(poolPresenceTable, rowPools, colPools, prPools, pcPools):
	
	import pdb
	
	run = []
	result = [run]
	expect = None
	CurrentLineMatchesPrevious = True

	i = 0
	while i < len(poolPresenceTable):
		if CurrentLineMatchesPrevious:
			run.append(poolPresenceTable[i])
		else:
			run = [poolPresenceTable[i]]
			result.append(run)
		
		CurrentLineMatchesPrevious = False
		
		if i < len(poolPresenceTable) - 1:
			CurrentLineMatchesPrevious = \
			DoesNextLineGroupWithCurrentLine2(poolPresenceTable[i+1], poolPresenceTable[i], \
			rowPools, colPools, prPools, pcPools)
			# print(CurrentLineMatchesPrevious)

		i += 1
		
    # Assign likely real genomic coordinates to the groups of read alignment coordinates
    # Do something with AT and TA locations in genome 
    
	return result
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def DoesNextLineGroupWithCurrentLine2(nextLine, currentLine, \
rowPools, colPools, prPools, pcPools, maxGap=1):
	
	# First, test if the two read alignment coordinates are adjacent to one another
	currentCoord = currentLine['readAlignmentCoord']
	expect = currentCoord + maxGap
	nextCoord = nextLine['readAlignmentCoord']
	
	nextLineMatchesWithCurrent = False
	
	if nextCoord <= expect:
		# print("coords are consecutive to within gap")
		poolPresenceProfileSimilar = \
		ArePoolPresenceProfilesSimilar2(nextLine, currentLine, rowPools, colPools, prPools, pcPools,\
		readUseThreshold=1)
		if poolPresenceProfileSimilar:
			nextLineMatchesWithCurrent = True
			
	return nextLineMatchesWithCurrent
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ArePoolPresenceProfilesSimilar2(nextLine, currentLine, rowPools, colPools, prPools, pcPools,\
readUseThreshold=1, diagnostics=False):
	
	from pdb import set_trace
	from scipy import intersect1d, setdiff1d, setxor1d, union1d

	
	currentLineCoords = CalculatePhysicalAddressCoordsForLine(currentLine, rowPools, colPools, \
	prPools, pcPools, threshold=readUseThreshold)
	
	nextLineCoords = CalculatePhysicalAddressCoordsForLine(nextLine, rowPools, colPools, prPools, \
	pcPools, threshold=readUseThreshold)
	
	# For each axis, test if there is first any intersection of coordinates between the two lines, 
	# and calculate the number of reads associated with this intersection.
	# Then, calculate the number of reads associated with the exclusive or of the coordinates sets
	# for each axis, and the number of reads associated here.
	
	similarityGreaterThanDifference = [False, False, False, False]
	
	i = 0
	while i < 4:
		axisCoordIntersections = intersect1d(nextLineCoords[i][0], currentLineCoords[i][0])
		axisCoordXOR = setxor1d(nextLineCoords[i][0], currentLineCoords[i][0])
		axisCoordUnion = union1d(nextLineCoords[i][0], currentLineCoords[i][0])
		axisCoordDiffNC = setdiff1d(nextLineCoords[i][0], currentLineCoords[i][0])
		axisCoordDiffCN = setdiff1d(currentLineCoords[i][0], nextLineCoords[i][0])
		
		
		axisCoordIntersectionsReads = \
		CountReadsAssociatedWithLineBooleanOperation(axisCoordIntersections, nextLineCoords[i], \
		currentLineCoords[i])
		
		axisCoordXORReads = \
		CountReadsAssociatedWithLineBooleanOperation(axisCoordXOR, nextLineCoords[i], \
		currentLineCoords[i])
		
		axisCoordUnionReads = \
		CountReadsAssociatedWithLineBooleanOperation(axisCoordUnion, nextLineCoords[i], \
		currentLineCoords[i])
		
		axisCoordDiffReadsNC = \
		CountReadsAssociatedWithLineBooleanOperation(axisCoordDiffNC, nextLineCoords[i], \
		currentLineCoords[i])
		
		axisCoordDiffReadsCN = \
		CountReadsAssociatedWithLineBooleanOperation(axisCoordDiffCN, nextLineCoords[i], \
		currentLineCoords[i])
		
		axisCoordDiffReads = axisCoordDiffReadsCN + axisCoordDiffReadsNC
		
		
		if diagnostics == True:
			print("currentLine: " + str(currentLine[0]))
			print("nextLine: " + str(nextLine[0]))
			print("i: " + str(i))
			print("axisCoordIntersectionsReads: " + str(axisCoordIntersectionsReads))
			print("axisCoordXORReads: " + str(axisCoordXORReads))
			print("axisCoordUnionReads: " + str(axisCoordUnionReads))
			print("axisCoordDiffReadsNC: " + str(axisCoordDiffReadsNC))
			print("axisCoordDiffReadsCN: " + str(axisCoordDiffReadsCN))
			print("axisCoordDiffReads: " + str(axisCoordDiffReads))


			print("\n")
				
		# Note, making a small change from version 1 of this algorithm, rather than
		# requiring axisCoordUnionReads >= axisCoordXORReads, I'm requiring 
		# axisCoordUnionReads > axisCoordXORReads
		# and I'm adding in the additional requirement that 
		# (axisCoordUnionReads > axisCoordDiffReads)

		if ((axisCoordIntersectionsReads > axisCoordXORReads) \
		or (axisCoordUnionReads > axisCoordXORReads)) and (axisCoordUnionReads > axisCoordDiffReads):
			similarityGreaterThanDifference[i] = True
	
		i += 1
	
# 	set_trace()
	i = 0
	axesWhereIntersectionGreaterThanXOR = 0
	while i < 4:
		if similarityGreaterThanDifference[i]:
			axesWhereIntersectionGreaterThanXOR += 1
		i += 1
	
	# Changing from >= 2 to > 2
	if axesWhereIntersectionGreaterThanXOR > 2:
		poolPresenceProfilesAreSimilar = True
	else:
		poolPresenceProfilesAreSimilar = False


	return poolPresenceProfilesAreSimilar
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateAverageReadAlignmentCoordinate(group, rowPools, colPools, prPools, pcPools, \
averagingType='median'):
	
	import numpy
	from pdb import set_trace
	import scipy.stats
	
	# Calculate the median and mean read alignment coordinate
	[readAlignmentCoords, totalReads] = \
	GenerateHistogramOfReadsVersusReadAlignmentCoord(group, rowPools, colPools, prPools, pcPools)
	
	readAlignmentCoordList = []
	i = 0
	while i < len(totalReads):
		j = 0
		while j < totalReads[i]:
			readAlignmentCoordList.append(readAlignmentCoords[i])
			j += 1
		i += 1
	
# 	set_trace()
	
# 	print(str(readAlignmentCoordList))
	
	if len(readAlignmentCoordList) == 0:
		averageReadAlignmentCoord = 0
		includeInSummedTable = False
	elif averagingType == 'median':
		averageReadAlignmentCoord = int(numpy.median(readAlignmentCoordList))
		includeInSummedTable = True
	elif averagingType == 'mode':
		averageReadAlignmentCoord = int(scipy.stats.mode(readAlignmentCoordList))
		includeInSummedTable = True
	elif averagingType == 'mean':
		averageReadAlignmentCoord = int(numpy.mean(readAlignmentCoordList))
		includeInSummedTable = True
	else:
		averageReadAlignmentCoord = int(numpy.median(readAlignmentCoordList))
		includeInSummedTable = True
	
	return [averageReadAlignmentCoord, includeInSummedTable]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def SumPoolPresenceTableGroup(group, dTypeDictForPoolPresenceTable, \
rowPools, colPools, prPools, pcPools, controlPools):
	
	import numpy
	
	from pdb import set_trace
	
	[consensusReadAlignmentCoord, includeInSummedTable] = \
	CalculateAverageReadAlignmentCoordinate(group, rowPools, colPools, prPools, pcPools, \
	averagingType='median')
	
	summedPoolPresenceTableEntry = numpy.zeros(1, dTypeDictForPoolPresenceTable)
	
# 	set_trace()
	
	if includeInSummedTable == True:	
		j = 0
	
		summedPoolPresenceTableEntry['readAlignmentCoord'] = consensusReadAlignmentCoord
	
		while j < len(group):
		
			for rowPool in rowPools:
				summedPoolPresenceTableEntry[rowPool] += int(group[j][rowPool])
		
			for colPool in colPools:
				summedPoolPresenceTableEntry[colPool] += int(group[j][colPool])
		
			for prPool in prPools:
				summedPoolPresenceTableEntry[prPool] += int(group[j][prPool])
		
			for pcPool in pcPools:
				summedPoolPresenceTableEntry[pcPool] += int(group[j][pcPool])
		
			for controlPool in controlPools:
				summedPoolPresenceTableEntry[controlPool] += int(group[j][controlPool])
		
		
			j += 1
		
# 	set_trace()
	
	return [summedPoolPresenceTableEntry, includeInSummedTable]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def SumPoolPresenceTable(groupedPoolPresenceTable, dTypeDictForPoolPresenceTable, \
rowPools, colPools, prPools, pcPools, controlPools):
	
	import numpy
	
	i = 0
	
	summedPoolPresenceTable = []
	
	while i < len(groupedPoolPresenceTable):
		
		group = groupedPoolPresenceTable[i]
		
		[summedPoolPresenceTableEntry, includeInSummedTable] = \
		SumPoolPresenceTableGroup(group, dTypeDictForPoolPresenceTable, \
		rowPools, colPools, prPools, pcPools, controlPools)
		
		if includeInSummedTable == True:
			summedPoolPresenceTable.append(summedPoolPresenceTableEntry)
		
		i += 1
	
	return summedPoolPresenceTable
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def SumPoolPresenceTableGroup2(group, dTypeDictForPoolPresenceTable, \
rowPools, colPools, prPools, pcPools, controlPools):
	
	import numpy
	
	from pdb import set_trace
	
	[consensusReadAlignmentCoord, includeInSummedTable] = \
	CalculateAverageReadAlignmentCoordinate(group, rowPools, colPools, prPools, pcPools, \
	averagingType='median')
	
	summedPoolPresenceTableEntry = numpy.zeros(1, dTypeDictForPoolPresenceTable)
	
# 	set_trace()
	
	if includeInSummedTable == True:	
		j = 0
	
		summedPoolPresenceTableEntry['readAlignmentCoord'] = consensusReadAlignmentCoord
	
		while j < len(group):
		
			for rowPool in rowPools:
				summedPoolPresenceTableEntry[rowPool] += int(group[j][rowPool])
		
			for colPool in colPools:
				summedPoolPresenceTableEntry[colPool] += int(group[j][colPool])
		
			for prPool in prPools:
				summedPoolPresenceTableEntry[prPool] += int(group[j][prPool])
		
			for pcPool in pcPools:
				summedPoolPresenceTableEntry[pcPool] += int(group[j][pcPool])
		
			for controlPool in controlPools:
				summedPoolPresenceTableEntry[controlPool] += int(group[j][controlPool])
		
		
			j += 1
	
	groupMemberCount = len(group)
	
	summedPoolPresenceTableEntry['groupMemberCount'] = groupMemberCount
	
# 	set_trace()
	
	return [summedPoolPresenceTableEntry, includeInSummedTable]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def SumPoolPresenceTable2(groupedPoolPresenceTable, dTypeDictForPoolPresenceTable, \
rowPools, colPools, prPools, pcPools, controlPools):
	
	import numpy
	
	i = 0
	
	summedPoolPresenceTable = []
	
	while i < len(groupedPoolPresenceTable):
		
		group = groupedPoolPresenceTable[i]
		
		[summedPoolPresenceTableEntry, includeInSummedTable] = \
		SumPoolPresenceTableGroup2(group, dTypeDictForPoolPresenceTable, \
		rowPools, colPools, prPools, pcPools, controlPools)
		
		if includeInSummedTable == True:
			summedPoolPresenceTable.append(summedPoolPresenceTableEntry)
		
		i += 1
	
	return summedPoolPresenceTable
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GenerateHistogramOfReadsVersusReadAlignmentCoord(groupedPoolPresenceTableGroup, \
rowPools, colPools, prPools, pcPools):

	i = 0
	
	readAlignmentCoords = []
	totalReads = []
	
	while i < len(groupedPoolPresenceTableGroup):
		
		readAlignmentCoord = groupedPoolPresenceTableGroup[i]['readAlignmentCoord']
		
		readAlignmentCoords.append(readAlignmentCoord)
		
		readCount = \
		CountReadsAssociatedWithCoordinateThatAreInLocationPools(groupedPoolPresenceTableGroup[i],\
		rowPools, colPools, prPools, pcPools)
		
		totalReads.append(readCount)
	
		i += 1
	return [readAlignmentCoords, totalReads]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def PlotHistogramOfReadsVersusReadAlignmentCoord(groupedPoolPresenceTableGroup, \
rowPools, colPools, prPools, pcPools):
	
	import matplotlib.pyplot as plt
	from matplotlib.pyplot import figure, plot
	
	[readAlignmentCoords, totalReads] = \
	GenerateHistogramOfReadsVersusReadAlignmentCoord(groupedPoolPresenceTableGroup, \
	rowPools, colPools, prPools, pcPools)
	
	figure()
	
	plot(readAlignmentCoords, totalReads, marker='o')
	
	return
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CountReadsAssociatedWithCoordinateThatAreInLocationPools(groupedPoolPresenceTableGroupLine,\
rowPools, colPools, prPools, pcPools):
	totalReads = 0
	
	i = 0
	while i < len(rowPools):
		totalReads += groupedPoolPresenceTableGroupLine[rowPools[i]]
		i += 1
	
	i = 0
	while i < len(colPools):
		totalReads += groupedPoolPresenceTableGroupLine[colPools[i]]
		i += 1
	
	i = 0
	while i < len(prPools):
		totalReads += groupedPoolPresenceTableGroupLine[prPools[i]]
		i += 1
		
	i = 0
	while i < len(pcPools):
		totalReads += groupedPoolPresenceTableGroupLine[pcPools[i]]
		i += 1

	return totalReads
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CountReadsAssociatedWithLineBooleanOperation(axisCoordIntersections, nextLineCoords, \
currentLineCoords):
	
	readsCount = 0
	
	for coord in axisCoordIntersections:
		if coord in nextLineCoords[0]:
			readsCount += nextLineCoords[1][coord]
		if coord in currentLineCoords[0]:
			readsCount += currentLineCoords[1][coord]
	
	return readsCount
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def ArePoolPresenceProfilesSimilar(nextLine, currentLine, rowPools, colPools, prPools, pcPools,\
readUseThreshold=1):
	
	from pdb import set_trace
	from scipy import intersect1d, setdiff1d, setxor1d, union1d

	
	currentLineCoords = CalculatePhysicalAddressCoordsForLine(currentLine, rowPools, colPools, \
	prPools, pcPools, threshold=readUseThreshold)
	
	nextLineCoords = CalculatePhysicalAddressCoordsForLine(nextLine, rowPools, colPools, prPools, \
	pcPools, threshold=readUseThreshold)
	
	# For each axis, test if there is first any intersection of coordinates between the two lines, 
	# and calculate the number of reads associated with this intersection.
	# Then, calculate the number of reads associated with the exclusive or of the coordinates sets
	# for each axis, and the number of reads associated here.
	
	similarityGreaterThanDifference = [False, False, False, False]
	
	i = 0
	while i < 4:
		axisCoordIntersections = intersect1d(nextLineCoords[i][0], currentLineCoords[i][0])
		axisCoordXOR = setxor1d(nextLineCoords[i][0], currentLineCoords[i][0])
		axisCoordUnion = union1d(nextLineCoords[i][0], currentLineCoords[i][0])
		
		axisCoordIntersectionsReads = \
		CountReadsAssociatedWithLineBooleanOperation(axisCoordIntersections, nextLineCoords[i], \
		currentLineCoords[i])
		
		axisCoordXORReads = \
		CountReadsAssociatedWithLineBooleanOperation(axisCoordXOR, nextLineCoords[i], \
		currentLineCoords[i])
		
		axisCoordUnionReads = \
		CountReadsAssociatedWithLineBooleanOperation(axisCoordUnion, nextLineCoords[i], \
		currentLineCoords[i])
		
		
		# print(str(nextLineCoords[i]) + ', ' + str(currentLineCoords[i]))
# 		print('axisCoordIntersectionsReads: ' + str(axisCoordIntersectionsReads))
# 		print('axisCoordXORReads: ' + str(axisCoordXORReads))
# 		
		if (axisCoordIntersectionsReads > axisCoordXORReads) \
		or (axisCoordUnionReads >= axisCoordXORReads):
			similarityGreaterThanDifference[i] = True
	
		i += 1
	
	i = 0
	axesWhereIntersectionGreaterThanXOR = 0
	while i < 4:
		if similarityGreaterThanDifference[i]:
			axesWhereIntersectionGreaterThanXOR += 1
		i += 1
	
	if axesWhereIntersectionGreaterThanXOR >= 2:
		poolPresenceProfilesAreSimilar = True
	else:
		poolPresenceProfilesAreSimilar = False


	return poolPresenceProfilesAreSimilar
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def WriteSummedPoolPresenceTable(filename, poolColumns, poolPresenceTable):
	
	from pdb import set_trace
	
	fhandle = open(filename, 'w')
	
	# Write out the header line
	
	headerLine = ''
	
	i = 0
	while i < len(poolColumns):
		headerLine += poolColumns[i] + ','
		i += 1
	headerLine += '\n'

	totalStr = ''
	
	i = 0
	while i < len(poolPresenceTable):
		outputStr = ''
		j = 0
# 		set_trace()
		while j < len(poolColumns):
			outputStr += str(poolPresenceTable[i][poolColumns[j]][0]) + ','
			j += 1
		
		outputStr += '\n'
		totalStr += outputStr
		i += 1
		
	writeStr = headerLine + totalStr
	
	fhandle.write(writeStr)
	fhandle.close()

	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def FindATandTAPositions(genomeFile):
	
	import re
	from pdb import set_trace
	
	
	fileHandle = open(genomeFile, 'r')
	fileData = fileHandle.readlines()
	
	commentRe = re.compile('>(\w+)')
	sequence = ''
	
	i = 0
	while i < len(fileData):
		line = fileData[i].strip()
		commentMatch = commentRe.match(line)
		
		if commentMatch != None:
			sequenceName = commentMatch.group(1)
		else:
			sequence += line
		
		i += 1

	
	ATandTAPositions = []
	
	atRegex = re.compile('at|ta', re.IGNORECASE)
	
# 	set_trace()
	
	i = 0
	while i < len(sequence) - 1:
		atMatch = atRegex.match(sequence[i:i+2])
		
		if atMatch != None:
			ATandTAPositions.append(i+1)
		
		i += 1
	
	
	
	return [ATandTAPositions, sequence]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportGenBankSequence(genBankFile):
	import re
	from pdb import set_trace
	
	originRe = re.compile(r'ORIGIN')
	endRe = re.compile(r'//')
	
	fileHandle = open(genBankFile, 'r')
	fileData = fileHandle.readlines()
	
	i = 0
	
	inSequence = False
	sequence = ''

	
	while i < len(fileData):
		
		line = fileData[i].strip()
		
		originMatch = originRe.match(line)
		
		if originMatch != None:
			inSequence = True
		
		endMatch = endRe.search(line)
		
		if endMatch != None:
			inSequence = False
		
		if (originMatch == None) and (inSequence == True):
			sequenceMatch = sequenceRe.search(line)
			
			if sequenceMatch != None:
				sequenceText = sequenceMatch.group(2)
				sequence += sequenceText
			else:
				print("Line not recognized")
	
		i += 1
	
	return sequence
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def WriteFastaSequence(sequence, name, fileName, width=100):
	
	fileHandle = open(fileName, 'w')
	
	fileHandle.write('>'+name+'\n')
	
	i = 0
	while i < len(sequence):
		
		if i > 0:
			if i%width == 0:
				fileHandle.write('\n')
		
		fileHandle.write(sequence[i])
		
		i += 1
	
	fileHandle.close()
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def WritePrimer3InputForDetectionPrimerSelection(fileName, sequence, sequenceName, sequenceTarget, \
productSizeRange, primerThermodynamicParameterPath, primerPairsToReturn=1):

	fileHandle = open(fileName, 'w')
	
	outputStr = 'SEQUENCE_ID=' + sequenceName + '\n'
	outputStr += 'SEQUENCE_TEMPLATE=' 
	outputStr += sequence
	outputStr += '\n'
	
	outputStr += 'SEQUENCE_TARGET=' + str(sequenceTarget[0]) + ',' + str(sequenceTarget[1]) + '\n'
	outputStr += 'PRIMER_TASK=pick_detection_primers\n'
	outputStr += 'PRIMER_PICK_LEFT_PRIMER=1\n'
	outputStr += 'PRIMER_PICK_INTERNAL_OLIGO=0\n'
	outputStr += 'PRIMER_PICK_RIGHT_PRIMER=1\n'
	outputStr += 'PRIMER_OPT_SIZE=21\n'
	outputStr += 'PRIMER_MIN_SIZE=15\n'
	outputStr += 'PRIMER_MAX_SIZE=30\n'
	outputStr += 'PRIMER_MAX_NS_ACCEPTED=1\n'
	
	outputStr += 'PRIMER_PRODUCT_SIZE_RANGE=' \
	+ str(productSizeRange[0]) + '-' + str(productSizeRange[1]) + '\n'
	
	outputStr += 'P3_FILE_FLAG=1\n'
	outputStr += 'SEQUENCE_INTERNAL_EXCLUDED_REGION=' \
	+ str(sequenceTarget[0]) + ',' + str(sequenceTarget[1]) + '\n'
	outputStr += 'PRIMER_EXPLAIN_FLAG=1\n'
	outputStr += 'PRIMER_LIBERAL_BASE=1\n'
	
	outputStr += 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH=' + primerThermodynamicParameterPath + '\n'
	
	outputStr += 'PRIMER_NUM_RETURN=' + str(primerPairsToReturn) + '\n'
	
	outputStr += '=\n'
	
	fileHandle.write(outputStr)

	return
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ImportFastaSequence(fastaFile):
	import re
	from pdb import set_trace
	
	
	fileHandle = open(fastaFile, 'r')
	fileData = fileHandle.readlines()
	
	commentRe = re.compile('>(\w+)')
	sequence = ''
	
	i = 0
	while i < len(fileData):
		line = fileData[i].strip()
		commentMatch = commentRe.match(line)
		
		if commentMatch != None:
			sequenceName = commentMatch.group(1)
		else:
			sequence += line
		
		i += 1
	
	return sequence
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ImportGenBankSequence(genBankFile):
	import re
	from pdb import set_trace
	
	originRe = re.compile(r'ORIGIN')
	endRe = re.compile(r'//')
	sequenceRe = re.compile(r'(\d+)\s+((\w+\s*)+)')
	
	fileHandle = open(genBankFile, 'r')
	fileData = fileHandle.readlines()
	
	i = 0
	
	inSequence = False
	sequence = ''

	sequenceLineCount = 0
	
	while i < len(fileData):
		
		line = fileData[i].strip()
		
		originMatch = originRe.match(line)
		
		if originMatch != None:
			inSequence = True
		
		endMatch = endRe.search(line)
		
		if endMatch != None:
			inSequence = False
		
		if (originMatch == None) and (inSequence == True):
			sequenceMatch = sequenceRe.search(line)
			
			if sequenceMatch != None:
				sequenceLineCount += 1				
				sequenceText = sequenceMatch.group(2)				
				sequence += sequenceText
			else:
				print("Line not recognized")
	
		i += 1
		

	sequenceNoSpaces = sequence.replace(" ", "")
	
	return sequenceNoSpaces
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def FindATandTAPositions2(genomeFile, format='genbank'):
# Does the same thing as FindATandTAPositions but can work with a GenBank or a Fasta file, \
# so you only  need one file format
	
	import re
	from pdb import set_trace
	
	if format == 'genbank':
		sequence = ImportGenBankSequence(genomeFile)
	elif format == 'fasta':
		sequence = ImportFastaSequence(genomeFile)
	
	ATandTAPositions = []
	
	atRegex = re.compile('(at|ta)', re.IGNORECASE)
	
# 	set_trace()
	
	i = 0
	while i < len(sequence) - 1:
		atMatch = atRegex.match(sequence[i:i+2])
		
		if atMatch != None:
			ATandTAPositions.append(i+1)
		
		i += 1
	
	
	
	return [ATandTAPositions, sequence]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CountEntriesInPoolSet(summedPoolPresenceTableLine, pools, countThreshold):
	
	import pdb
	
	readsInPoolGreaterThanThresholdCount = 0
	
	for pool in pools:
		print(pool + ': ' + str(summedPoolPresenceTableLine[pool][0]))
		if summedPoolPresenceTableLine[pool][0] > countThreshold:
			readsInPoolGreaterThanThresholdCount += 1

	return readsInPoolGreaterThanThresholdCount
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def FindUnambiguousEntriesInPoolPresenceTable(summedPoolPresenceTable, rowPools, colPools, \
prPools, pcPools, controlPools, countThresholdForPools=5, countThresholdForControls=1):
	
	import pdb
	
	unambiguousPoolPresenceTableEntries = []
	i = 0
	
	while i < len(summedPoolPresenceTable):
		
		print('i: ' + str(i))
		
		rowPoolEntriesAboveThresholdCount = CountEntriesInPoolSet(summedPoolPresenceTable[i], \
		rowPools, countThresholdForPools)
		
		print('rowPoolEntriesAboveThresholdCount: ' + str(rowPoolEntriesAboveThresholdCount))
	
		columnPoolEntriesAboveThresholdCount = CountEntriesInPoolSet(summedPoolPresenceTable[i], \
		colPools, countThresholdForPools)
		
		print('columnPoolEntriesAboveThresholdCount: ' + str(columnPoolEntriesAboveThresholdCount))
	
	
		plateRowPoolEntriesAboveThresholdCount = CountEntriesInPoolSet(summedPoolPresenceTable[i], \
		prPools, countThresholdForPools)
	
		print('plateRowPoolEntriesAboveThresholdCount: ' + str(plateRowPoolEntriesAboveThresholdCount))
	
		
		plateColPoolEntriesAboveThresholdCount = CountEntriesInPoolSet(summedPoolPresenceTable[i], \
		pcPools, countThresholdForPools)
		
		print('plateColPoolEntriesAboveThresholdCount: ' + str(plateColPoolEntriesAboveThresholdCount))
	
		controlsEntriesAboveThresholdCount = CountEntriesInPoolSet(summedPoolPresenceTable[i], \
		controlPools, countThresholdForControls)
	
		print('controlsEntriesAboveThresholdCount: ' + str(controlsEntriesAboveThresholdCount))
	
		
		if (rowPoolEntriesAboveThresholdCount == 1) and (columnPoolEntriesAboveThresholdCount == 1) \
		and (plateRowPoolEntriesAboveThresholdCount == 1) \
		and (plateColPoolEntriesAboveThresholdCount == 1) and \
		(controlsEntriesAboveThresholdCount == 0):
			unambiguousPoolPresenceTableEntries.append(summedPoolPresenceTable[i])
			print('Adding to unambiguousPoolPresenceTableEntries')
			
		
		i += 1
		
	return unambiguousPoolPresenceTableEntries
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CountFeaturesWithMoreThanOneDisruption(unambiguousPoolPresenceTableEntries, featureArray):
	i = 0

	while i < len(unambiguousPoolPresenceTableEntries):
	
		readAlignmentCoord = unambiguousPoolPresenceTableEntries[i]['readAlignmentCoord']
	
		j = 0
		while j < len(featureArray):
			featureStart = featureArray[j].startCoord
			featureEnd = featureArray[j].endCoord
		
			if featureStart < readAlignmentCoord < featureEnd:
				featureArray[j].poolPresenceEntries.append(unambiguousPoolPresenceTableEntries[i])
			j += 1
		
		i += 1


	i = 0
	featuresWithMoreThanOneDisruption = 0
	while i < len(featureArray):
		if len(featureArray[i].poolPresenceEntries) > 0:
			featuresWithMoreThanOneDisruption += 1
		i += 1
	return featuresWithMoreThanOneDisruption
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculatePhysicalAddressForPoolPresenceTableLine(poolPresenceTableLine, rowPools, colPools, prPools, pcPools, \
threshold=10):
	
	physicalAddressDict = GenerateBlankPhysicalAddressDict(poolPresenceTableLine)
	
	i = 0
	while i < len(poolPresenceTableLine):
		entry = poolPresenceTableLine[i]
		coord = int(poolPresenceTableLine[i]['readAlignmentCoord'])
		
		addresses_r = FindAddressCoords(poolPresenceTableLine, rowPools, threshold=threshold)
		addresses_c = FindAddressCoords(poolPresenceTableLine, colPools, threshold=threshold)
		addresses_pr = FindAddressCoords(poolPresenceTableLine, prPools, threshold=threshold)
		addresses_pc = FindAddressCoords(poolPresenceTableLine, pcPools, threshold=threshold)
		
		possibleAddresses = CalculatePossibleAddresses(addresses_r, addresses_c, addresses_pr, \
		addresses_pc)
		
		physicalAddressDict[coord] = possibleAddresses
		
		i += 1
	
	return physicalAddressDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ReturnClosestDisruptionToTranslationStart(feature):
	if len(feature.poolPresenceEntries) > 0:
# 		print('More than 0 pool presence entries')
		i = 0
		closestI = 0
		
		end = feature.endTranslation
		start = feature.startTranslation
		
		distanceOfClosestIToTranslationStart = abs(end - start)
# 		print(distanceOfClosestIToTranslationStart)
		
# 		print(str(feature.endTranslation))
		
		while i < len(feature.poolPresenceEntries):
			readAlignmentCoord = feature.poolPresenceEntries[i]['readAlignmentCoord'][0]
			startTranslation = feature.startTranslation 
# 			print(readAlignmentCoord)
# 			print(startTranslation)
			distanceToTranslationStart = abs(readAlignmentCoord - startTranslation)
# 			print(str(distanceToTranslationStart))
			
			if distanceToTranslationStart < distanceOfClosestIToTranslationStart:
				closestI = i
				distanceOfClosestIToTranslationStart = distanceToTranslationStart
			
			# print(closestI)
# 			print(distanceOfClosestIToTranslationStart)
# 			
			i += 1
		
		closestDisruption = feature.poolPresenceEntries[closestI]
		
		return closestDisruption
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def UpdateFeatureArrayWithPoolPresenceTableEntries(unambiguousPoolPresenceTableEntries, \
featureArray):

	i = 0

	while i < len(unambiguousPoolPresenceTableEntries):
	
		readAlignmentCoord = unambiguousPoolPresenceTableEntries[i]['readAlignmentCoord']
	
		j = 0
		while j < len(featureArray):
			featureStart = featureArray[j].startCoord
			featureEnd = featureArray[j].endCoord
		
			if featureStart < readAlignmentCoord < featureEnd:
				featureArray[j].poolPresenceEntries.append(unambiguousPoolPresenceTableEntries[i])
			j += 1
		
		i += 1


# 	i = 0
# 	featuresWithMoreThanOneDisruption = 0
# 	while i < len(featureArray):
# 		if len(featureArray[i].poolPresenceEntries) > 0:
# 			featuresWithMoreThanOneDisruption += 1
# 		i += 1
		
	return featureArray
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CalculateAddressesOfClosestDisruptionsToFeatureStarts(featuresWithDisruptions, \
rowPools, colPools, prPools, pcPools):
	
	from numpy import arange
	import numpy
	import pdb
	
	# This is a big kludge the generate the layout of the plate grid.
	# In future versions, this will be read in from a file rather than being hard coded
	# in like this. 
	
	i = 0
	plateRowLookupDict = {}
	while i < len(prPools):
		plateRowLookupDict[prPools[i]] = list(arange(21*i+1,21*i+22,1,dtype=int))
		i += 1
	

	plateRowLookupDict['PR20'][18] = 36
	plateRowLookupDict['PR20'][19] = 191
	
	i = 0

	featuresWithDisruptionsNamesAndAddresses = []

	while i < len(featuresWithDisruptions):
		name = featuresWithDisruptions[i].featureName
		start = featuresWithDisruptions[i].startTranslation
		end = featuresWithDisruptions[i].endTranslation
	
		closestDisruption = ReturnClosestDisruptionToTranslationStart(featuresWithDisruptions[i])
	
		readAlign = closestDisruption['readAlignmentCoord']
	
		address = CalculatePhysicalAddressForPoolPresenceTableLine(closestDisruption, rowPools, \
		colPools, prPools, pcPools, threshold=10)
	
		addressString = list(address.values())[0][0][0].split('_')
		
		humanReadableWellAddress = addressString[0] + addressString[1]
	
		plateRow = addressString[2]
		platCol = addressString[3]
	
		plate = \
		plateRowLookupDict[plateRow][pcPools.index(platCol)]

		featuresWithDisruptionsNamesAndAddresses.append([name, start, end, readAlign, address, \
		humanReadableWellAddress,  plate])
	
		i += 1
	
	return featuresWithDisruptionsNamesAndAddresses
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def WriteFeaturesWithDisruptionsAndAddresses(featuresWithDisruptionsNamesAndAddresses, \
geneDisruptionLocationsFileName):
	
	i = 0
	geneDisruptionListHandle = open(geneDisruptionLocationsFileName , 'w')

	geneDisruptionListHandle.write('Gene name, Translation Start, Translation End, Read Alignment Coord, Well, Plate\n')

	while i < len(featuresWithDisruptionsNamesAndAddresses):
		outputStr = ''
		outputStr += str(featuresWithDisruptionsNamesAndAddresses[i][0]) + ','
		outputStr += str(featuresWithDisruptionsNamesAndAddresses[i][1]) + ','
		outputStr += str(featuresWithDisruptionsNamesAndAddresses[i][2]) + ','
		outputStr += str(featuresWithDisruptionsNamesAndAddresses[i][3][0]) + ','
		outputStr += str(featuresWithDisruptionsNamesAndAddresses[i][5]) + ','
		outputStr += str(featuresWithDisruptionsNamesAndAddresses[i][6]) + '\n'
		geneDisruptionListHandle.write(outputStr)
		i += 1

	geneDisruptionListHandle.close()
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculatePoolCoordsForLine(poolPresenceTableLine, rowPools, colPools, prPools, \
pcPools, controlPools, threshold=5):
# Very much like CalculatePhysicalAddressCoordsForLine but also reports contents of control pools
# as well.  
# Used in generation of pool presence table taxonomy

	addresses_r = FindAddressCoords2(poolPresenceTableLine, rowPools, threshold=threshold)
	addresses_c = FindAddressCoords2(poolPresenceTableLine, colPools, threshold=threshold)
	addresses_pr = FindAddressCoords2(poolPresenceTableLine, prPools, threshold=threshold)
	addresses_pc = FindAddressCoords2(poolPresenceTableLine, pcPools, threshold=threshold)
	addresses_control = FindAddressCoords2(poolPresenceTableLine, controlPools, threshold=threshold)
	
	# Remember, each line in the addresses array is a 2 element list, the first containing the pool
	# name, and the second containing the number of reads associated with it. 
	
	return [addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def FindAddressCoords3(poolPresenceTableLine, axisPools, threshold=5):
	# Very much the same as FindAddressCoords, but returns a list of just axis addresses and a 
	# dict of the read counts for each of these axis addresses, rather than combining this dict
	# into the axisAddresses array
	
	
	axisAddressScoreDict = {}
	
	i = 0
	while i < len(axisPools):
		if poolPresenceTableLine[axisPools[i]] >= threshold:
			axisAddressScoreDict[axisPools[i]] = poolPresenceTableLine[axisPools[i]]
		i += 1
	

	return axisAddressScoreDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindLineNumberForCoord(poolPresenceTable, coord):
	i = 0
	found = False
	while (i < len(poolPresenceTable)) and (found == False):
		if poolPresenceTable[i]['readAlignmentCoord'] == coord:
			found = True
		else:
			i += 1
			
	index = i
	indexFound = found
	
	return [index, indexFound]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculatePoolCoordsForLine2(poolPresenceTableLine, spatialPools, controlPools, \
spatialThreshold=1, controlThreshold=0):
# Very much like CalculatePhysicalAddressCoordsForLine but also reports contents of control pools
# as well.  
# Used in generation of pool presence table taxonomy

	addresses_spatial = FindAddressCoords3(poolPresenceTableLine, spatialPools, \
	threshold=spatialThreshold)
	
	addresses_control = FindAddressCoords3(poolPresenceTableLine, controlPools, \
	threshold=controlThreshold)
	
	
	addressDict = {}
	
	for key in addresses_spatial.keys():
		addressDict[key] = addresses_spatial[key]
	
	for key in controlPools:
		addressDict[key] = addresses_control[key]
	

	
	return addressDict
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def CalculateLibraryAddressesForPoolPresenceTableLine(poolPresenceTableLine, rowPools, colPools, \
prPools, pcPools, threshold=5):
# Used in generation of pool presence table taxonomy
# Used to calculate possible addresses for pool presence table line
	
	coord = int(poolPresenceTableLine['readAlignmentCoord'])
	
	addresses_r = FindAddressCoords(poolPresenceTableLine, rowPools, threshold=threshold)
	addresses_c = FindAddressCoords(poolPresenceTableLine, colPools, threshold=threshold)
	addresses_pr = FindAddressCoords(poolPresenceTableLine, prPools, threshold=threshold)
	addresses_pc = FindAddressCoords(poolPresenceTableLine, pcPools, threshold=threshold)
	
	possibleAddresses = CalculatePossibleAddresses(addresses_r, addresses_c, addresses_pr, \
	addresses_pc)
	
	return possibleAddresses
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, addresses_pc, \
addresses_control):
# Used in generation of pool presence table taxonomy
# Used to calculate how many lines can be used to calculate library addresses
	
	nRowPoolAxisEntries = len(addresses_r[0])
	nColPoolAxisEntries = len(addresses_c[0])
	nPRPoolAxisEntries = len(addresses_pr[0])
	nPCPoolAxisEntries = len(addresses_pc[0])
	
	nAddressPoolAxesWithEntries = 0
	if nRowPoolAxisEntries > 0:
		nAddressPoolAxesWithEntries += 1
	if nColPoolAxisEntries > 0:
		nAddressPoolAxesWithEntries += 1
	if nPRPoolAxisEntries > 0:
		nAddressPoolAxesWithEntries += 1
	if nPCPoolAxisEntries > 0:
		nAddressPoolAxesWithEntries += 1
		
	
	nControlPoolsWithEntries = len(addresses_control[0])

	return [nAddressPoolAxesWithEntries, nControlPoolsWithEntries]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, \
addresses_pr, addresses_pc, addresses_control):
	
	
	nRowPoolAxisEntries = len(addresses_r[0])
	nColPoolAxisEntries = len(addresses_c[0])
	nPRPoolAxisEntries = len(addresses_pr[0])
	nPCPoolAxisEntries = len(addresses_pc[0])
	
	nAddressPoolAxesWithMoreThanOneEntry = 0
	if nRowPoolAxisEntries > 1:
		nAddressPoolAxesWithMoreThanOneEntry += 1
	if nColPoolAxisEntries > 1:
		nAddressPoolAxesWithMoreThanOneEntry += 1
	if nPRPoolAxisEntries > 1:
		nAddressPoolAxesWithMoreThanOneEntry += 1
	if nPCPoolAxisEntries > 1:
		nAddressPoolAxesWithMoreThanOneEntry += 1
	

	return nAddressPoolAxesWithMoreThanOneEntry
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateMaxNumberOfEntriesInSinglePoolAxis(addresses_r, addresses_c, \
addresses_pr, addresses_pc):
# This function finds the pool axis with the most coordinate entries and reports this number
# This number is important for guessing how many library addresses an ambiguous line might map to
# It is the minimum number of addresses that the line could map to. 
	
	nRowPoolAxisEntries = len(addresses_r[0])
	nColPoolAxisEntries = len(addresses_c[0])
	nPRPoolAxisEntries = len(addresses_pr[0])
	nPCPoolAxisEntries = len(addresses_pc[0])
	
	poolAxisEntries = [nRowPoolAxisEntries, nColPoolAxisEntries, nPRPoolAxisEntries, \
	nPCPoolAxisEntries]
	
	maxEntriesInSingleAddressPoolAxis = max(poolAxisEntries)
	

	return maxEntriesInSingleAddressPoolAxis
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def PrintSummedPoolPresenceLineTaxonomyDict(summedPoolPresenceLineTaxonomyDict, keysInPrintOrder):
	
	outputStr = ''
	
	for key in keysInPrintOrder:
		outputStr += key + ': ' + str(summedPoolPresenceLineTaxonomyDict[key]) + '\n'
	
	print(outputStr)
	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateSummedPoolPresenceTableTaxonomyDict(summedPoolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools):

	keysInPrintOrder = [\
	'totalGroups',\
	'totalLines',\
	'linesThatDoNotGroup',\
	'linesThatGroup',\
	'linesThatGroupThatMapToLibraryAddresses',\
	'linesThatGroupThatMapToSingleLibraryAddresses',\
	'linesThatGroupThatMapToMultipleLibraryAddresses',\
	'linesThatGroupThatMapToUnambiguousLibraryAddresses',\
	'linesThatGroupThatMapToAmbiguousLibraryAddresses', \
	'linesThatGroupThatDoNotMapToLibraryAddresses',\
	'linesThatGroupThatHaveNoReadsAboveThresholdInAnyPool',\
	'linesThatGroupThatMapToControlIndexesOnly',\
	'linesThatGroupThatHaveCoordinatesInNoPoolAxis',\
	'linesThatGroupThatHaveCoordinatesInOnlyOnePoolAxis',\
	'linesThatGroupThatHaveCoordinatesInOnlyTwoPoolAxes',\
	'linesThatGroupThatHaveCoordinatesInOnlyThreePoolAxes',\
	'linesThatDoNotGroupThatMapToLibraryAddresses',\
	'linesThatDoNotGroupThatMapToSingleLibraryAddresses',\
	'linesThatDoNotGroupThatMapToMultipleLibraryAddresses',\
	'linesThatDoNotGroupThatMapToUnambiguousLibraryAddresses',\
	'linesThatDoNotGroupThatMapToAmbiguousLibraryAddresses', \
	'linesThatDoNotGroupThatDoNotMapToLibraryAddresses',\
	'linesThatDoNotGroupThatHaveNoReadsAboveThresholdInAnyPool',\
	'linesThatDoNotGroupThatMapToControlIndexesOnly',\
	'linesThatDoNotGroupThatHaveCoordinatesInNoPoolAxis',\
	'linesThatDoNotGroupThatHaveCoordinatesInOnlyOnePoolAxis',\
	'linesThatDoNotGroupThatHaveCoordinatesInOnlyTwoPoolAxes',\
	'linesThatDoNotGroupThatHaveCoordinatesInOnlyThreePoolAxes']



	# Define the summed pool presence line taxonomy dict 
	summedPoolPresenceLineTaxonomyDict = {}

	for key in keysInPrintOrder:
		summedPoolPresenceLineTaxonomyDict[key] = 0


	summedPoolPresenceTablePoolsCoordsList = {}
	possibleAddressesDict = {}
	numberPoolAxesWithEntriesForLine = {}
	numberPoolAxesWithMoreThanOneEntry = {}

	i = 0
	while i < len(summedPoolPresenceTable):
	
		coord = int(summedPoolPresenceTable[i]['readAlignmentCoord'])
		
		summedPoolPresenceLineTaxonomyDict['totalLines'] += 1
		groupMemberCount = int(summedPoolPresenceTable[i]['groupMemberCount'])
	
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		summedPoolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=5)
	
		possibleAddressesDict[coord] = possibleAddresses
		lenPossibleAddresses = len(possibleAddresses)
	
		[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
		CalculatePoolCoordsForLine(\
		summedPoolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, threshold=5)
	
		poolCoordsForLine = [addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] 
	
		summedPoolPresenceTablePoolsCoordsList[coord] = poolCoordsForLine
	
		[nAddressPoolAxesWithEntries, nControlPoolsWithEntries] = \
		CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, addresses_pc, \
		addresses_control)
	
		numberPoolAxesWithEntriesForLine[coord] = [nAddressPoolAxesWithEntries, nControlPoolsWithEntries]
	
	
		# Calculate if there will be an ambiguous library address calculation by calculating the number
		# of pool axes with more than one entry. One axis with multiple (even if possible coords are
		# filled) is fine, but more axis with more than one entry leads to cross terms that are 
		# ambiguous. 
		nAddressPoolAxesWithMoreThanOneEntry = \
		CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
		addresses_pc, addresses_control)
	
		numberPoolAxesWithMoreThanOneEntry[coord] = nAddressPoolAxesWithMoreThanOneEntry
	
		outputStr = str(coord) + ',' + str(len(possibleAddresses)) + ',' + str(nAddressPoolAxesWithEntries) + ',' + str(nControlPoolsWithEntries) + '\n'
		#print(outputStr)
	
		
		if groupMemberCount == 1:
			summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroup'] += 1
		
			if (nAddressPoolAxesWithEntries == 0) and (nControlPoolsWithEntries == 0):
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatHaveNoReadsAboveThresholdInAnyPool'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
			if (nAddressPoolAxesWithEntries == 0) and (nControlPoolsWithEntries > 0):
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatMapToControlIndexesOnly'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
			if (nAddressPoolAxesWithEntries == 1):
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatHaveCoordinatesInOnlyOnePoolAxis'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
			if (nAddressPoolAxesWithEntries == 2):
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatHaveCoordinatesInOnlyTwoPoolAxes'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
			if (nAddressPoolAxesWithEntries == 3):
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatHaveCoordinatesInOnlyThreePoolAxes'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord))
				
			if (nAddressPoolAxesWithEntries == 0):
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatHaveCoordinatesInNoPoolAxis'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
		
			if nAddressPoolAxesWithMoreThanOneEntry > 1:
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatMapToAmbiguousLibraryAddresses'] += 1
			
			if nAddressPoolAxesWithMoreThanOneEntry <= 1 and lenPossibleAddresses >= 1:
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatMapToUnambiguousLibraryAddresses'] += 1
		
		
			if lenPossibleAddresses == 0:
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatDoNotMapToLibraryAddresses'] += 1
			
			elif lenPossibleAddresses == 1:
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatMapToSingleLibraryAddresses'] += 1
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatMapToLibraryAddresses'] += 1
			elif lenPossibleAddresses > 1:
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatMapToMultipleLibraryAddresses'] += 1
				summedPoolPresenceLineTaxonomyDict['linesThatDoNotGroupThatMapToLibraryAddresses'] += 1
	
	
	
		elif groupMemberCount > 1:
			summedPoolPresenceLineTaxonomyDict['totalGroups'] += 1
			summedPoolPresenceLineTaxonomyDict['linesThatGroup'] += groupMemberCount

			if (nAddressPoolAxesWithEntries == 0) and (nControlPoolsWithEntries == 0):
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatHaveNoReadsAboveThresholdInAnyPool'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
			if (nAddressPoolAxesWithEntries == 0) and (nControlPoolsWithEntries > 0):
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatMapToControlIndexesOnly'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
			if (nAddressPoolAxesWithEntries == 1):
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatHaveCoordinatesInOnlyOnePoolAxis'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
			if (nAddressPoolAxesWithEntries == 2):
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatHaveCoordinatesInOnlyTwoPoolAxes'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
			if (nAddressPoolAxesWithEntries == 3):
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatHaveCoordinatesInOnlyThreePoolAxes'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord)) 
		
			if (nAddressPoolAxesWithEntries == 0):
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatHaveCoordinatesInNoPoolAxis'] += 1
				if lenPossibleAddresses != 0:
					print(str(coord))
		
		
			if nAddressPoolAxesWithMoreThanOneEntry > 1:
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatMapToAmbiguousLibraryAddresses'] += 1


			if nAddressPoolAxesWithMoreThanOneEntry <= 1 and lenPossibleAddresses >= 1:
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatMapToUnambiguousLibraryAddresses'] += 1
		
			
		
			if lenPossibleAddresses == 0:
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatDoNotMapToLibraryAddresses'] += 1
			elif lenPossibleAddresses == 1:
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatMapToSingleLibraryAddresses'] += 1
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatMapToLibraryAddresses'] += 1
			elif lenPossibleAddresses > 1:
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatMapToMultipleLibraryAddresses'] += 1
				summedPoolPresenceLineTaxonomyDict['linesThatGroupThatMapToLibraryAddresses'] += 1
		i += 1
	
	PrintSummedPoolPresenceLineTaxonomyDict(summedPoolPresenceLineTaxonomyDict, keysInPrintOrder)
	
	return summedPoolPresenceLineTaxonomyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateTotalReadsInLine(poolPresenceTableLine, rowPools, colPools, prPools, pcPools, controlPools):

	i = 0
	totalReads = 0
	while i < len(rowPools):
		totalReads += poolPresenceTableLine[rowPools[i]]
		i +=1
	
	i = 0
	while i < len(colPools):
		totalReads += poolPresenceTableLine[colPools[i]]
		i +=1
	
	i = 0
	while i < len(prPools):
		totalReads += poolPresenceTableLine[prPools[i]]
		i +=1

	
	i = 0
	while i < len(pcPools):
		totalReads += poolPresenceTableLine[pcPools[i]]
		i +=1	

	i = 0
	while i < len(controlPools):
		totalReads += poolPresenceTableLine[controlPools[i]]
		i +=1
	
	return totalReads
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindEntriesInPoolPresenceLine(poolPresenceTableLine, poolsList):
	
	entries = []
	
	j = 0
	while j < len(poolsList):
		if 	poolPresenceTableLine[poolsList[j]] > 0:
			entries.append(poolPresenceTableLine[poolsList[j]])
		j += 1

	return entries
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculatePoolPresenceTableTotalReadCountHistogram(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools):
	
	poolPresenceTableHistogram = []
	poolPresenceTableHistogramNoControls = []
	
	i = 0
	while i < len(poolPresenceTable):
	
		lineEntries = []

		poolPresenceTableLine = poolPresenceTable[i]
		
		rowEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, rowPools)
		colEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, colPools)
		prEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, prPools)
		pcEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, pcPools)
		controlEntries = FindEntriesInPoolPresenceLine(poolPresenceTableLine, controlPools)
		
		lineEntries = rowEntries + colEntries + prEntries + pcEntries + controlEntries
		lineEntriesNoControls = rowEntries + colEntries + prEntries + pcEntries
		
		poolPresenceTableHistogram += lineEntries
		poolPresenceTableHistogramNoControls += lineEntriesNoControls
		
		i += 1

		
	
	return [poolPresenceTableHistogram, poolPresenceTableHistogramNoControls]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def PlotCumulativeDensityOfTotalReadCounts(base, cumulative, baseNC, cumulativeNC):
	
	import matplotlib.pyplot as plt
	
	plt.figure()
	plt.semilogx(base[:-1],cumulative/max(cumulative))
	plt.title("Cumulative Distribution of Total Read Count for All Pool Presence Table Entries")
	plt.xlabel("Total Read Count for Line")
	plt.ylabel("Fraction of All Non-Zero Pool Presence Table Entries") 
	plt.show()


	plt.figure()
	plt.semilogx(baseNC[:-1],cumulativeNC/max(cumulativeNC))
	plt.title("Cumulative Distribution of Total Read Count for Non-Control Pool Presence Table Entries")
	plt.xlabel("Total Read Count for Line")
	plt.ylabel("Fraction of All Non-Zero Pool Presence Table Entries") 
	plt.show()
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def PlotHistogramsOfTotalReadCountsPerLine(totalReadCountPerLineData, \
totalReadCountPerLineDataNoControls):
	
	import matplotlib.pyplot as plt
	
	plt.figure()
	plt.hist(totalReadCountPerLineData)
	plt.xlabel("Total Read Count Per Line")
	plt.ylabel("Number of Pool Presence Table Lines")
	plt.title("Histogram of Total Read Count Distribution in Pool Presence Table")
	


	plt.figure()
	plt.hist(totalReadCountPerLineDataNoControls)
	plt.xlabel("Total Read Count Per Line Excluding Control Counts")
	plt.ylabel("Number of Pool Presence Table Lines")
	plt.title("Histogram of Total Read Count Distribution without Controls in Pool Presence Table")
	
	
	plt.show()

# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CalculatePlotAndSavePoolPresenceTableTotalReadCountHistogram(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, cumulativePlotFileName, cumulativePlotFileNameNC):
# Generate, plot and save the pool presence table total read count per line histogram

	from vectorOutput2 import generateOutputMatrixWithHeaders, writeOutputMatrix
	import numpy
	import pdb

	[totalReadCountPerLineData, totalReadCountPerLineDataNoControls] = \
	CalculatePoolPresenceTableTotalReadCountHistogram(\
	poolPresenceTable, rowPools, colPools, prPools, pcPools, controlPools)
	
	# Plot histograms
	PlotHistogramsOfTotalReadCountsPerLine(totalReadCountPerLineData, \
	totalReadCountPerLineDataNoControls)
	
	
	# Cumulative density plots 

	bins = numpy.arange(0, max(totalReadCountPerLineData), 1)
	binsNC = numpy.arange(0, max(totalReadCountPerLineDataNoControls), 1)

	values, base = numpy.histogram(totalReadCountPerLineData, bins=bins)
	valuesNC, baseNC = numpy.histogram(totalReadCountPerLineDataNoControls, bins=binsNC)

	cumulative = numpy.cumsum(values)
	cumulativeNC = numpy.cumsum(valuesNC)
	
	PlotCumulativeDensityOfTotalReadCounts(base, cumulative, baseNC, cumulativeNC)
	
	# Save cumulative density plot data
	
	vectorList = [base[:-1], cumulative, cumulative/max(cumulative)]
	vectorListNC = [baseNC[:-1], cumulativeNC, cumulativeNC/max(cumulativeNC)]
	
	headers = ["Total Reads", "Cum. Lines", "Frac Cum. Lines"]
	headersNC = ["Total Reads NC", "Cum. Lines NC", "Frac Cum. Lines NC"]
	
# 	pdb.set_trace()

	oMatrix = generateOutputMatrixWithHeaders(vectorList, headers)
	oMatrixNC = generateOutputMatrixWithHeaders(vectorListNC, headersNC)
	
	writeOutputMatrix(cumulativePlotFileName, oMatrix)
	writeOutputMatrix(cumulativePlotFileNameNC, oMatrixNC)
	
	return
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
addresses_pr, addresses_pc):
	
	import pdb
	
	nRowCoords = len(addresses_r[0])
	nColCoords = len(addresses_c[0])
	nPRCoords = len(addresses_pr[0])
	nPCCoords = len(addresses_pc[0])
	
# 	pdb.set_trace()
	
	return [nRowCoords, nColCoords, nPRCoords, nPCCoords]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateCoordNumbersForAmbiguousAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, threshold):

	import pdb
	
	# Initialize the number arrays
	
	nRowCoordsArray = []
	nColCoordsArray = []
	nPRCoordsArray = []
	nPCCoordsArray = []
	
	i = 0
	
	while i < len(poolPresenceTable):
	
		# Figure out if the line is ambiguous or not
		
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		lenPossibleAddresses = len(possibleAddresses)
	
		[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
		CalculatePoolCoordsForLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, threshold=threshold)
	
		nAddressPoolAxesWithMoreThanOneEntry = \
		CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
		addresses_pc, addresses_control)
		
		# Here's the test for ambiguous mapping to multiple library addresses
		
		if nAddressPoolAxesWithMoreThanOneEntry > 1 and lenPossibleAddresses >= 1:
			
			[nRowCoords, nColCoords, nPRCoords, nPCCoords] = \
			CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
			addresses_pr, addresses_pc)
			
			nRowCoordsArray.append(nRowCoords)
			nColCoordsArray.append(nColCoords)
			nPRCoordsArray.append(nPRCoords)
			nPCCoordsArray.append(nPCCoords)
		
		i += 1

	return [nRowCoordsArray, nColCoordsArray, nPRCoordsArray, nPCCoordsArray]

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculatePossibleAddressNumbersAndCoordsForAmbiguousAddressLines(poolPresenceTable, rowPools, \
colPools, prPools, pcPools, controlPools, threshold):

	import pdb
	
	# Initialize the number arrays
	
	nPossibleAddressesArray = []
	nTotalCoordsArray = []
	maxNumberCoordsArray = []
	
	i = 0
	
	while i < len(poolPresenceTable):
	
		# Figure out if the line is ambiguous or not
		
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		lenPossibleAddresses = len(possibleAddresses)
	
		[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
		CalculatePoolCoordsForLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, threshold=threshold)
	
		nAddressPoolAxesWithMoreThanOneEntry = \
		CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
		addresses_pc, addresses_control)
		
		# Here's the test for ambiguous mapping to multiple library addresses
		
		if nAddressPoolAxesWithMoreThanOneEntry > 1 and lenPossibleAddresses >= 1:
			
			nPossibleAddressesArray.append(lenPossibleAddresses)
			
			[nRowCoords, nColCoords, nPRCoords, nPCCoords] = \
			CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
			addresses_pr, addresses_pc)
			
			nTotalCoordsArray.append(nRowCoords + nColCoords + nPRCoords + nPCCoords)
			
			maxNumberCoordsArray.append(max([nRowCoords, nColCoords, nPRCoords, nPCCoords]))
		
		i += 1

	return [nPossibleAddressesArray, nTotalCoordsArray, maxNumberCoordsArray]

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateAndPlotPossibleAddressNumbersAndCoordsForAmbiguousAddressLines(poolPresenceTable, rowPools, \
colPools, prPools, pcPools, controlPools, threshold):

	import matplotlib.pyplot as pyplot
	from numpy import arange, histogram, cumsum

	[nPossibleAddressesArray, nTotalCoordsArray, maxNumberCoordsArray] = \
	CalculatePossibleAddressNumbersAndCoordsForAmbiguousAddressLines(poolPresenceTable, rowPools, \
	colPools, prPools, pcPools, controlPools, threshold)
	
	nPossibleAddressesBins = arange(0,max(nPossibleAddressesArray),1)
	nTotalCoordsBins = arange(0,max(nTotalCoordsArray),1)
	maxNumberCoordsBins = arange(0,max(maxNumberCoordsArray),1)


	pyplot.figure()
	pyplot.hist(nPossibleAddressesArray, bins=nPossibleAddressesBins)
	pyplot.xlabel("Number Possible Addresses")
	pyplot.ylabel("Number of Ambiguous Lines")

	pyplot.figure()
	pyplot.hist(nTotalCoordsArray, bins=nTotalCoordsBins)
	pyplot.xlabel("Total Number of Coords")
	pyplot.ylabel("Number of Ambiguous Lines")
	
	
	pyplot.figure()
	pyplot.hist(maxNumberCoordsArray, bins=maxNumberCoordsBins)
	pyplot.xlabel("Maximum Number of Coords in a Pool Axis")
	pyplot.ylabel("Number of Ambiguous Lines")

	pyplot.show()
	
	
	binsPossAdd = arange(0, max(nPossibleAddressesArray), 1)
	valuesPossAdd, basePossAdd = histogram(nPossibleAddressesArray, bins=binsPossAdd)
	cumulativePossAdd = cumsum(valuesPossAdd)
	
	binsTotalCoords = arange(0, max(nTotalCoordsArray), 1)
	valuesTotalCoords, baseTotalCoords = histogram(nTotalCoordsArray, bins=binsTotalCoords)
	cumulativeTotalCoords = cumsum(valuesTotalCoords)
	
	binsMaxCoords = arange(0, max(maxNumberCoordsArray), 1)
	valuesMaxCoords, baseMaxCoords = histogram(maxNumberCoordsArray, bins=binsMaxCoords)
	cumulativeMaxCoords = cumsum(valuesMaxCoords)

	pyplot.figure()
	pyplot.plot(basePossAdd[:-1], cumulativePossAdd/max(cumulativePossAdd))
	pyplot.title("Cumulative Distribution of Possible Addresses for Ambiguous Pool Presence Table Entries")
	pyplot.xlabel("Possible Addresses")
	pyplot.ylabel("Fraction of Ambiguous Pool Presence Table Entries") 
	pyplot.show()
	
	pyplot.figure()
	pyplot.plot(baseTotalCoords[:-1], cumulativeTotalCoords/max(cumulativeTotalCoords))
	pyplot.title("Cumulative Distribution of Total Coords for Ambiguous Pool Presence Table Entries")
	pyplot.xlabel("Total Coordinates")
	pyplot.ylabel("Fraction of Ambiguous Pool Presence Table Entries") 
	pyplot.show()
	
	pyplot.figure()
	pyplot.plot(baseMaxCoords[:-1], cumulativeMaxCoords/max(cumulativeMaxCoords))
	pyplot.title("Cumulative Distribution of Maximum Pool Coordinates for Ambiguous Pool Presence Table Entries")
	pyplot.xlabel("Maximum Number of Coordinates in a Pool")
	pyplot.ylabel("Fraction of Ambiguous Pool Presence Table Entries") 
	pyplot.show()
	

	return
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def ReportAmbiguousPoolPresenceTableLines(poolPresenceTable, rowPools, \
colPools, prPools, pcPools, controlPools, threshold, nTotalCoordsQuery):
	
	import pdb
	import numpy
	
	# Initialize the number arrays
	
	reportedPoolPresenceTable = poolPresenceTable[0]
	
	i = 0
	
	while i < len(poolPresenceTable):
	
		# Figure out if the line is ambiguous or not
		
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		lenPossibleAddresses = len(possibleAddresses)
	
		[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
		CalculatePoolCoordsForLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, threshold=threshold)
	
		nAddressPoolAxesWithMoreThanOneEntry = \
		CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
		addresses_pc, addresses_control)
		
		# Here's the test for ambiguous mapping to multiple library addresses
		
		if nAddressPoolAxesWithMoreThanOneEntry > 1 and lenPossibleAddresses >= 1:
			
			[nRowCoords, nColCoords, nPRCoords, nPCCoords] = \
			CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
			addresses_pr, addresses_pc)
			
			
			
			nTotalCoords = nRowCoords + nColCoords + nPRCoords + nPCCoords
			
			if nTotalCoords == nTotalCoordsQuery:
# 				pdb.set_trace()
# 				print(poolPresenceTable[i])
				reportedPoolPresenceTable = numpy.append(reportedPoolPresenceTable, \
				poolPresenceTable[i])
			
		i += 1
	
	reportedPoolPresenceTable = numpy.delete(reportedPoolPresenceTable, 0)

	return reportedPoolPresenceTable
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CalculateReadNumberRatiosForSingleAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, threshold):
	
	import pdb
	from numpy import log
	
	# Initialize the ratio arrays
	
	nr2ncArray = []
	nr2nprArray = []
	nr2npcArray = []
	nc2nprArray = []
	nc2npcArray = []
	npr2npcArray = []
	
	logNr2ncArray = []
	logNr2nprArray = []
	logNr2npcArray = []
	logNc2nprArray = []
	logNc2npcArray = []
	logNpr2npcArray = []
	
	
	i = 0
	
	while i < len(poolPresenceTable):
	
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		# Calculate if the line is single occupancy
		
		lenPossibleAddresses = len(possibleAddresses)
	
		if lenPossibleAddresses == 1:
	
			[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
			CalculatePoolCoordsForLine(\
			poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, \
			threshold=threshold)
			
			
			nr = addresses_r[1][addresses_r[0][0]]
			nc = addresses_c[1][addresses_c[0][0]]
			npr = addresses_pr[1][addresses_pr[0][0]]
			npc = addresses_pc[1][addresses_pc[0][0]]
					
			nr2nc = nr/nc
			nr2npr = nr/npr
			nr2npc = nr/npc
			nc2npr = nc/npr
			nc2npc = nc/npc
			npr2npc = npr/npc

			nr2ncArray.append(nr2nc)
			nr2nprArray.append(nr2npr) 
			nr2npcArray.append(nr2npc)
			nc2nprArray.append(nc2npr)
			nc2npcArray.append(nc2npc)
			npr2npcArray.append(npr2npc)
			
			logNr2ncArray.append(log(nr2nc))
			logNr2nprArray.append(log(nr2npr)) 
			logNr2npcArray.append(log(nr2npc))
			logNc2nprArray.append(log(nc2npr))
			logNc2npcArray.append(log(nc2npc))
			logNpr2npcArray.append(log(npr2npc))

			nr2ncArray.append(nr2nc)
			nr2nprArray.append(nr2npr) 
			nr2npcArray.append(nr2npc)
			nc2nprArray.append(nc2npr)
			nc2npcArray.append(nc2npc)
			npr2npcArray.append(npr2npc)
		
		i += 1

	return [nr2ncArray, nr2nprArray, nr2npcArray, nc2nprArray, nc2npcArray, npr2npcArray, \
	logNr2ncArray, logNr2nprArray, logNr2npcArray, logNc2nprArray, logNc2npcArray, logNpr2npcArray]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculatePoolPresenceTableTaxonomyDict(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, threshold):

	keysInPrintOrder = [\
	'totalLines',\
	'linesThatMapToLibraryAddresses',\
	'linesThatMapToSingleLibraryAddresses',\
	'linesThatMapToMultipleLibraryAddresses',\
	'linesThatMapToUnambiguousLibraryAddresses',\
	'linesThatMapToAmbiguousLibraryAddresses', \
	'linesThatDoNotMapToLibraryAddresses',\
	'linesThatHaveNoReadsAboveThresholdInAnyPool',\
	'linesThatMapToControlIndexesOnly',\
	'linesThatHaveCoordinatesInNoPoolAxis',\
	'linesThatHaveCoordinatesInOnlyOnePoolAxis',\
	'linesThatHaveCoordinatesInOnlyTwoPoolAxes',\
	'linesThatHaveCoordinatesInOnlyThreePoolAxes',\
	'totalReads']

	# Define the pool presence line taxonomy dict 
	poolPresenceTableTaxonomyDict = {}

	for key in keysInPrintOrder:
		poolPresenceTableTaxonomyDict[key] = 0

	poolPresenceTablePoolsCoordsList = {}
	possibleAddressesDict = {}
	numberPoolAxesWithEntriesForLine = {}
	numberPoolAxesWithMoreThanOneEntry = {}

	i = 0
	while i < len(poolPresenceTable):
		
		readsInLine = CalculateTotalReadsInLine(poolPresenceTable[i], \
		rowPools, colPools, prPools, pcPools, controlPools)
		
		poolPresenceTableTaxonomyDict['totalReads'] += readsInLine
		
		coord = int(poolPresenceTable[i]['readAlignmentCoord'])
		
		poolPresenceTableTaxonomyDict['totalLines'] += 1
		
		possibleAddresses = CalculateLibraryAddressesForPoolPresenceTableLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		lenPossibleAddresses = len(possibleAddresses)
	
		[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
		CalculatePoolCoordsForLine(\
		poolPresenceTable[i], rowPools, colPools, prPools, pcPools, controlPools, threshold=threshold)
	
		poolCoordsForLine = [addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] 
		
		[nAddressPoolAxesWithEntries, nControlPoolsWithEntries] = \
		CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, addresses_pc, \
		addresses_control)
	
		
		# Calculate if there will be an ambiguous library address calculation by calculating the number
		# of pool axes with more than one entry. One axis with multiple (even if possible coords are
		# filled) is fine, but more axis with more than one entry leads to cross terms that are 
		# ambiguous. 
		nAddressPoolAxesWithMoreThanOneEntry = \
		CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
		addresses_pc, addresses_control)
	

		if (nAddressPoolAxesWithEntries == 0) and (nControlPoolsWithEntries == 0):
			poolPresenceTableTaxonomyDict['linesThatHaveNoReadsAboveThresholdInAnyPool'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
		if (nAddressPoolAxesWithEntries == 0) and (nControlPoolsWithEntries > 0):
			poolPresenceTableTaxonomyDict['linesThatMapToControlIndexesOnly'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
		if (nAddressPoolAxesWithEntries == 1):
			poolPresenceTableTaxonomyDict['linesThatHaveCoordinatesInOnlyOnePoolAxis'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
		if (nAddressPoolAxesWithEntries == 2):
			poolPresenceTableTaxonomyDict['linesThatHaveCoordinatesInOnlyTwoPoolAxes'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
		if (nAddressPoolAxesWithEntries == 3):
			poolPresenceTableTaxonomyDict['linesThatHaveCoordinatesInOnlyThreePoolAxes'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord)) 
	
		if (nAddressPoolAxesWithEntries == 0):
			poolPresenceTableTaxonomyDict['linesThatHaveCoordinatesInNoPoolAxis'] += 1
			if lenPossibleAddresses != 0:
				print(str(coord))
	
	
		if nAddressPoolAxesWithMoreThanOneEntry > 1 and lenPossibleAddresses >= 1:
			poolPresenceTableTaxonomyDict['linesThatMapToAmbiguousLibraryAddresses'] += 1


		if nAddressPoolAxesWithMoreThanOneEntry <= 1 and lenPossibleAddresses >= 1:
			poolPresenceTableTaxonomyDict['linesThatMapToUnambiguousLibraryAddresses'] += 1
	
		
	
		if lenPossibleAddresses == 0:
			poolPresenceTableTaxonomyDict['linesThatDoNotMapToLibraryAddresses'] += 1
		elif lenPossibleAddresses == 1:
			poolPresenceTableTaxonomyDict['linesThatMapToSingleLibraryAddresses'] += 1
			poolPresenceTableTaxonomyDict['linesThatMapToLibraryAddresses'] += 1
		elif lenPossibleAddresses > 1:
			poolPresenceTableTaxonomyDict['linesThatMapToMultipleLibraryAddresses'] += 1
			poolPresenceTableTaxonomyDict['linesThatMapToLibraryAddresses'] += 1
		
		i += 1
	
	PrintPoolPresenceTableTaxonomyDict(threshold, poolPresenceTableTaxonomyDict, keysInPrintOrder)
	
	return poolPresenceTableTaxonomyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def PrintPoolPresenceTableTaxonomyDict(threshold, poolPresenceTableTaxonomyDict, keysInPrintOrder):
	
	outputStr = ''
	
	outputStr += 'threshold: ' + str(threshold) + '\n'
	
	for key in keysInPrintOrder:
		outputStr += key + ': ' + str(poolPresenceTableTaxonomyDict[key]) + '\n'
	
	print(outputStr)
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
class SudokuPlate:

	def __init__(self, plateName, plateRow, plateCol, rowPools, colPools):
		
		self.plateName = plateName
		self.rowPools = rowPools
		self.colPools = colPools
		
		self.wellGrid = self.generateWellGrid(plateName, plateRow, plateCol, rowPools, colPools)
	
	def generateWellGrid(self, plateName, plateRow, plateCol, rowPools, colPools):
		
		wellGrid = {}
		
		i = 0
		while i < len(rowPools):
			wellGrid[rowPools[i]] = {}
			j = 0
			while j < len(colPools):
				wellGrid[rowPools[i]][colPools[j]] = SudokuWell(plateName, plateRow, plateCol, \
				rowPools[i], colPools[j])
				j += 1
			i += 1
		
		return wellGrid
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class SudokuWell:
	def __init__(self, plateName, plateRow, plateCol, row, col, OD=1.0):
		self.plateName = plateName
		self.row = row
		self.col = col
		self.libraryAddress = str(plateRow) + '_' + str(plateCol) + '_' +  str(row) + '_' + str(col)
		self.OD = OD
		self.readAlignmentCoords = []

 
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
class SudokuGenomicCoord:
	def __init__(self, coord, locatability, locatabilityScore, readCount):
		self.coord = coord
		self.locatability = locatability
		self.locatabilityScore = locatabilityScore
		self.readCount = readCount

# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ImportSudokuGridLayout(sudokuGridLayoutFile, rowPools, colPools):
	
	fileHandle = open(sudokuGridLayoutFile, 'r')
	
	data = fileHandle.readlines()
	
	pcPools = data[0].strip().split(',')
	pcPools.pop(0)

	
	prPools = []
	plateNameGrid = []
	
	i = 1
	while i < len(data):
		
		dataLine = data[i].strip().split(',')
		
		prPools.append(dataLine[0])
		
		plateNameGrid.append([])
		
		j = 1
		while j < len(dataLine):
			 plateNameGrid[i-1].append(dataLine[j])
			 j += 1
			 
		i += 1
		
	# The sudokuGridLookupDict is very similar to the plateRowLookupDict that was hard-coded into
	# earlier versions of the sudoku code. It is initially indexed by the plate row, followed by 
	# plate columns
	
	sudokuGridLookupDict = {}
	
	i = 0
	while i < len(prPools):
		sudokuGridLookupDict[prPools[i]] = {}
		j = 0
		while j < len(pcPools):
		
			plateName = plateNameGrid[i][j]
			
			sudokuGridLookupDict[prPools[i]][pcPools[j]] = SudokuPlate(plateName, \
			prPools[i], pcPools[j], rowPools, colPools)

			j += 1
		i += 1

	
	return [sudokuGridLookupDict, prPools, pcPools]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def PopulateSudokuGrid(sudokuGridLookupDict, summedPoolPresenceTable, rowPools, colPools, \
prPools, pcPools, controlPools, threshold):

	import pdb
	
	i = 0
	while i < len(summedPoolPresenceTable):
		
		summedPoolPresenceTableLine = summedPoolPresenceTable[i]
		
		coord = int(summedPoolPresenceTableLine['readAlignmentCoord'][0])
		
		# Calculate if the library address will be unambiguous, guessable, ambiguous, or un- 
		# locatable  by calculating the 
		# number of pool axes with more than one entry. One axis with multiple (even if all 
		# possible coords are filled) is fine, but more axis with more than one entry leads 
		# to cross terms that are ambiguous. 
		
		addresses = CalculateLibraryAddressesForPoolPresenceTableLine(summedPoolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, threshold=threshold)
		
		if len(addresses) > 0:
			

			[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
			CalculatePoolCoordsForLine(\
			summedPoolPresenceTableLine, rowPools, colPools, prPools, pcPools, controlPools, \
			threshold=threshold)
			
			[nAddressPoolAxesWithEntries, nControlPoolsWithEntries] = \
			CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, addresses_pc, \
			addresses_control)
	
			nAddressPoolAxesWithMoreThanOneEntry = \
			CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, addresses_pr, \
			addresses_pc, addresses_control)


			if nAddressPoolAxesWithEntries < 3:
				locatability = 'unlocatable'
			elif nAddressPoolAxesWithEntries == 3:
				if nAddressPoolAxesWithMoreThanOneEntry > 1:
					locatability = 'ambiguous'
				elif nAddressPoolAxesWithMoreThanOneEntry <= 1:
					locatability = 'guessable'
			elif nAddressPoolAxesWithEntries == 4:
				if nAddressPoolAxesWithMoreThanOneEntry > 1:
					locatability = 'ambiguous'
				elif nAddressPoolAxesWithMoreThanOneEntry <= 1:
					locatability = 'unambiguous'
			else:
				locatability = 'unlocatable'
		
			if locatability == 'unambiguous':
				j = 0
				while j < len(addresses):
# 					pdb.set_trace()
					addressCoords = addresses[j][0].split('_')
					
					row = addressCoords[0]
					col = addressCoords[1]
					pr = addressCoords[2]
					pc = addressCoords[3]
					
					outputStr =  str(coord) + ': [' + str(pr) + ',' + str(pc) + ',' + str(row) + ',' + str(col) + ']'
					print(outputStr)
										
					sudokuGridLookupDict[pr][pc].wellGrid[row][col].readAlignmentCoords.append(coord)
					
					j += 1
					
		
		i += 1
	
	return 
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def AddAddressToSudokuGrid(coord, locatability, possibleAddressesAndScoresEntry, \
sudokuGridLookupDict):
	
	addressCoords = possibleAddressesAndScoresEntry[0].split('_')
	locatabilityScore = possibleAddressesAndScoresEntry[2]
	readCount = possibleAddressesAndScoresEntry[3]
					 					
	row = addressCoords[0]
	col = addressCoords[1]
	pr = addressCoords[2]
	pc = addressCoords[3]
		
	sudokuCoord = SudokuGenomicCoord(coord, locatability, locatabilityScore, readCount)
	
	sudokuGridLookupDict[pr][pc].wellGrid[row][col].readAlignmentCoords.append(sudokuCoord)
	
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def PopulateSudokuGrid2(sudokuGridLookupDict, poolPresenceTable, rowPools, colPools, \
prPools, pcPools, controlPools, fitDict, areaDict, threshold, scoreThreshold):

	import pdb
	import operator
	
	i = 0
	while i < len(poolPresenceTable):
		
		poolPresenceTableLine = poolPresenceTable[i]
		
		coord = poolPresenceTableLine['readAlignmentCoord']
				
		
		possibleAddressesAndScores = \
		CalculateLibraryAddressesForPoolPresenceTableLine2(poolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, fitDict, \
		areaDict, threshold)
		
		if len(possibleAddressesAndScores) > 0:
		
			[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
			CalculatePoolCoordsForLine(\
			poolPresenceTableLine, rowPools, colPools, prPools, pcPools, controlPools, \
			threshold=threshold)
			
			[nAddressPoolAxesWithEntries, nControlPoolsWithEntries] = \
			CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, \
			addresses_pc, \
			addresses_control)
	
			nAddressPoolAxesWithMoreThanOneEntry = \
			CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, \
			addresses_pr, addresses_pc, addresses_control)
			
			maxEntriesInSingleAddressPoolAxis = \
			CalculateMaxNumberOfEntriesInSinglePoolAxis(addresses_r, addresses_c, \
			addresses_pr, addresses_pc)

			# Decide on the locatability of the genomic coordinate
			if nAddressPoolAxesWithEntries < 3:
				locatability = 'unlocatable'
			elif nAddressPoolAxesWithEntries == 3:
				if nAddressPoolAxesWithMoreThanOneEntry > 1:
					locatability = 'unlocatable'
				elif nAddressPoolAxesWithMoreThanOneEntry <= 1:
					locatability = 'guessable'
			elif nAddressPoolAxesWithEntries == 4:
				if nAddressPoolAxesWithMoreThanOneEntry > 1:
					locatability = 'ambiguous'
				elif nAddressPoolAxesWithMoreThanOneEntry <= 1:
					locatability = 'unambiguous'
			else:
				locatability = 'unlocatable'
		
		
			# If the genomic coordinate address is ambiguous or unambiguous, add it to the 
			# Sudoku grid
			if locatability == 'unambiguous':
				j = 0
				while j < len(possibleAddressesAndScores):
					AddAddressToSudokuGrid(coord, locatability, possibleAddressesAndScores[j], \
					sudokuGridLookupDict)
					j += 1
				
			elif locatability == 'ambiguous':
				
				sortedPossibleAddressesAndScores = sorted(possibleAddressesAndScores, \
				key=operator.itemgetter(2), reverse=True)
				
# 				pdb.set_trace()
				
				j = 0
				continueAddingAddressesToSudokuGrid = True
				
				while continueAddingAddressesToSudokuGrid == True:
					
					AddAddressToSudokuGrid(coord, locatability, \
					sortedPossibleAddressesAndScores[j], sudokuGridLookupDict)
					
					j += 1
					
					if (j < maxEntriesInSingleAddressPoolAxis) \
					and j < len(sortedPossibleAddressesAndScores):
						continueAddingAddressesToSudokuGrid = True
					else:
						continueAddingAddressesToSudokuGrid = False
					
					# elif j >= maxEntriesInSingleAddressPoolAxis:
# 						if j < len(sortedPossibleAddressesAndScores):
# 							nextScore = sortedPossibleAddressesAndScores[j][2]
# 							if nextScore < scoreThreshold:
# 								continueAddingAddressesToSudokuGrid = False
# 							else:
# 								continueAddingAddressesToSudokuGrid = True
# 						else:
# 							continueAddingAddressesToSudokuGrid = False
					

		i += 1
	
	return 
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateLibraryAddressLocatabilityForPoolPresenceTableLine(poolPresenceTableLine, \
rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, threshold, \
maxEntriesInSingleAddressPoolAxis, maxTotalCoords):
	
	import pdb
	
	[addresses_r, addresses_c, addresses_pr, addresses_pc, addresses_control] = \
	CalculatePoolCoordsForLine(\
	poolPresenceTableLine, rowPools, colPools, prPools, pcPools, controlPools, \
	threshold=threshold)
	
	[nRowCoords, nColCoords, nPRCoords, nPCCoords] = \
	CalculateNumberOfEntriesInPoolAxes(addresses_r, addresses_c, \
	addresses_pr, addresses_pc)
	
	nTotalCoords = nRowCoords + nColCoords + nPRCoords + nPCCoords
	
	[nAddressPoolAxesWithEntries, nControlPoolsWithEntries] = \
	CalculateNumberOfPoolAxesThatHaveEntries(addresses_r, addresses_c, addresses_pr, \
	addresses_pc, addresses_control)

	nAddressPoolAxesWithMoreThanOneEntry = \
	CalculateNumberOfPoolAxesThatHaveMoreThanOneEntry(addresses_r, addresses_c, \
	addresses_pr, addresses_pc, addresses_control)
	
	maxSinglePoolCoordNumber = \
	CalculateMaxNumberOfEntriesInSinglePoolAxis(addresses_r, addresses_c, \
	addresses_pr, addresses_pc)

	# Decide on the locatability of the genomic coordinate
	if nAddressPoolAxesWithEntries < 3:
		locatability = 'unlocatable'
		possibleAddressesAndScores = None
		
	
	elif maxSinglePoolCoordNumber > maxEntriesInSingleAddressPoolAxis or \
	nTotalCoords > maxTotalCoords:
		locatability = 'unlocatable'
		possibleAddressesAndScores = None
	
	elif nAddressPoolAxesWithEntries == 3:
		if nAddressPoolAxesWithMoreThanOneEntry > 1:
			locatability = 'unlocatable'
			possibleAddressesAndScores = None
		elif nAddressPoolAxesWithMoreThanOneEntry <= 1:
			locatability = 'guessable'
			possibleAddressesAndScores = None
	
	elif nAddressPoolAxesWithEntries == 4:
		
		possibleAddressesAndScores = \
		CalculateLibraryAddressesForPoolPresenceTableLine2(poolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, fitDict, \
		areaDict, threshold)
		
		if len(possibleAddressesAndScores) == 0:
			pdb.set_trace()
	
		if nAddressPoolAxesWithMoreThanOneEntry > 1:
			locatability = 'ambiguous'
		elif nAddressPoolAxesWithMoreThanOneEntry <= 1:
			locatability = 'unambiguous'
	
	else:
		locatability = 'unlocatable'
		possibleAddressesAndScores = None

	return [possibleAddressesAndScores, locatability, maxSinglePoolCoordNumber]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def AssignUnambiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
locatability):
	j = 0
	while j < len(possibleAddressesAndScores):
		AddAddressToSudokuGrid(coord, locatability, possibleAddressesAndScores[j], \
		sudokuGridLookupDict)
		j += 1

	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def AssignAmbiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
locatability, maxEntriesInSingleAddressPoolAxis):

	import operator
	import pdb
	
	sortedPossibleAddressesAndScores = sorted(possibleAddressesAndScores, \
	key=operator.itemgetter(2), reverse=True)
	
	# if  coord==5038711:
# 		pdb.set_trace()
	
	j = 0
	continueAddingAddressesToSudokuGrid = True
	
	while continueAddingAddressesToSudokuGrid == True:
		
		try:
			addressAndScore = sortedPossibleAddressesAndScores[j]
		except IndexError:
			pdb.set_trace()
		
		
		AddAddressToSudokuGrid(coord, locatability, addressAndScore, sudokuGridLookupDict)
		
		j += 1
		
		if (j < maxEntriesInSingleAddressPoolAxis) \
		and j < len(sortedPossibleAddressesAndScores):
			continueAddingAddressesToSudokuGrid = True
		else:
			continueAddingAddressesToSudokuGrid = False

	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AssignGuessableAddresses(sudokuGridLookupDict, coord, poolPresenceTableLine, \
rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, readCountThreshold, \
maxSinglePoolCoordNumber, maxTotalCoords):
	
	temporaryThreshold = 1
	validAddressesFound = False
	
	while validAddressesFound == False and temporaryThreshold <= readCountThreshold:
		
		[possibleAddressesAndScores, locatability, maxSinglePoolCoordNumber] = \
		CalculateLibraryAddressLocatabilityForPoolPresenceTableLine(poolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, temporaryThreshold, \
		maxSinglePoolCoordNumber, maxTotalCoords)
		
		if locatability == 'unambiguous':
			validAddressesFound = True
		else:
			temporaryThreshold += 1
	
	if validAddressesFound == True:
		AssignUnambiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
		'unambiguous')
		
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def PopulateSudokuGrid3(sudokuGridLookupDict, poolPresenceTable, rowPools, colPools, \
prPools, pcPools, controlPools, fitDict, areaDict, readCountThreshold, scoreThreshold, \
maxEntriesInSingleAddressPoolAxis, maxTotalCoords):

	import pdb
	
	i = 0
	while i < len(poolPresenceTable):
		
		poolPresenceTableLine = poolPresenceTable[i]
		
# 		pdb.set_trace()
		
		coord = poolPresenceTableLine['readAlignmentCoord']
		
		[possibleAddressesAndScores, locatability, maxSinglePoolCoordNumber] = \
		CalculateLibraryAddressLocatabilityForPoolPresenceTableLine(poolPresenceTableLine, \
		rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, readCountThreshold, \
		maxEntriesInSingleAddressPoolAxis, maxTotalCoords)
		
		# if  (coord==5038710) or (coord==5038711) or (coord==5038712):
# 			pdb.set_trace()
	
		if locatability == 'unambiguous':
			AssignUnambiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
			locatability)
			
		elif locatability == 'ambiguous':
			AssignAmbiguousAddresses(sudokuGridLookupDict, possibleAddressesAndScores, coord, \
			locatability, maxSinglePoolCoordNumber)
		
		elif locatability == 'guessable':
			AssignGuessableAddresses(sudokuGridLookupDict, coord, poolPresenceTableLine, \
			rowPools, colPools, prPools, pcPools, controlPools, fitDict, areaDict, \
			readCountThreshold, maxEntriesInSingleAddressPoolAxis, maxTotalCoords)
							
		i += 1
	
	return 
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GroupSudokuGenomicCoords(coordArray, maxGap=1):
	
	import pdb
	import numpy
	import operator
	
	i = 0
	expect = None
	run = []
	result = [run]
	currentLineMatchesPrevious = True
	
	while i < len(coordArray):
		if currentLineMatchesPrevious:
			run.append(coordArray[i])
		else:
			run = [coordArray[i]]
			result.append(run)
			
		currentLineMatchesPrevious = False
		
		if i < len(coordArray) - 1:
			currentCoord = coordArray[i].coord
			nextCoord = coordArray[i+1].coord
			expect = currentCoord + maxGap
			if nextCoord <= expect:
				currentLineMatchesPrevious = True
		i += 1
		
		
	groupedCoords = []
	
	i = 0
	while i < len(result):
		
		if len(result[i]) == 1:
# 			pdb.set_trace()
			groupedCoords.append(result[i][0])
		
		elif len(result[i]) > 1:
			coords = sorted(result[i], key=operator.attrgetter('readCount'), reverse=True)
			
			j = 0
			readCountList = []
			while j < len(coords):
				k = 0
				while k < coords[j].readCount:
					readCountList.append(coords[j].coord)
					k += 1
				j += 1
			
			representativeCoord = numpy.median(readCountList)
			
			j = 0 
			totalReadCount = 0
			locatabilityArray = []
			while j < len(coords):
				totalReadCount += coords[j].readCount
				locatabilityArray.append(coords[j].locatability)
				j +=1
						
			if 'unambiguous' in locatabilityArray:
				locatability = 'unambiguous'
			elif 'ambiguous' in locatabilityArray:
				locatability = 'ambiguous'
			elif 'guessable' in locatabilityArray:
				locatability = 'guessable'
			else:
				locatability = 'merged'
			
			
			locatabilityScore = 'merged'
			
			groupedCoords.append(SudokuGenomicCoord(representativeCoord, locatability, \
			locatabilityScore, totalReadCount))
		
		i += 1
	
	return groupedCoords
	
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GroupReadAlignmentCoordsInSudokuGrid(sudokuGridLookupDict, prPools, pcPools, rowPools, \
colPools, maxGap=4):
	
	import operator
	
	for prPool in prPools:
		
		for pcPool in pcPools:
		
			sudokuPlate = None
			try:
				sudokuPlate = sudokuGridLookupDict[prPool][pcPool]
			except IndexError:
				print('No plate at: ' + rowPool + '_' + colPool)
				pass
			
		
			if sudokuPlate != None:
				plateName = sudokuPlate.plateName
		
				for colPool in colPools:
					for rowPool in rowPools:
						sudokuWell = sudokuPlate.wellGrid[rowPool][colPool]
						readAlignmentCoords = sudokuWell.readAlignmentCoords
						
						readAlignmentCoords = sorted(readAlignmentCoords, \
						key=operator.attrgetter('coord'))
						
						groupedReadAlignmentCoords = GroupSudokuGenomicCoords(readAlignmentCoords, \
						maxGap=maxGap)
						sudokuWell.readAlignmentCoords = groupedReadAlignmentCoords
						
						
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GroupReadAlignmentIntegerCoordsArray(coordArray, maxGap=1):
	
	i = 0
	expect = None
	run = []
	result = [run]
	currentLineMatchesPrevious = True
	
	while i < len(coordArray):
		if currentLineMatchesPrevious:
			run.append(coordArray[i])
		else:
			run = [coordArray[i]]
			result.append(run)
			
		currentLineMatchesPrevious = False
		
		if i < len(coordArray) - 1:
			currentCoord = coordArray[i]
			nextCoord = coordArray[i+1]
			expect = currentCoord + maxGap
			if nextCoord <= expect:
				currentLineMatchesPrevious = True
		i += 1
		
	
	i = 0
	maxGroupSize = 1
	while i < len(result):
		groupSize = len(result[i])
		if groupSize > maxGroupSize:
			maxGroupSize = groupSize
			
		i += 1
	
	if maxGroupSize > 1:
		groupingNeeded = True
	else:
		groupingNeeded = False
	
	return [result, groupingNeeded]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def TestIfGroupingIsNeededInSudokuGrid(sudokuGridLookupDict, prPools, pcPools, rowPools, colPools,\
maxGap=4):

	import pdb

	groupingNeededInSinglyOccupiedWells = 0

	for prPool in prPools:
		
		for pcPool in pcPools:
		
			sudokuPlate = None
			try:
				sudokuPlate = sudokuGridLookupDict[prPool][pcPool]
			except IndexError:
				print('No plate at: ' + rowPool + '_' + colPool)
				pass
			
		
			if sudokuPlate != None:
				plateName = sudokuPlate.plateName
		
				for colPool in colPools:
					for rowPool in rowPools:
						sudokuWell = sudokuPlate.wellGrid[rowPool][colPool]
						readAlignmentCoords = sudokuWell.readAlignmentCoords
					
						coordArray = []
						
# 						pdb.set_trace()
					
						for coord in readAlignmentCoords:
							try:
								coordArray.append(coord.coord)
							except:
								pdb.set_trace()
					
						try:
							coordArray = sorted(coordArray)
						except:
							pdb.set_trace()
					
						[groupedCoords, groupingNeeded] =  \
						GroupReadAlignmentIntegerCoordsArray(coordArray, maxGap=maxGap)
					
						lenGroupedCoords = len(groupedCoords)
					
						if (groupingNeeded == True) and (lenGroupedCoords == 1):
					
							outputStr = plateName + ' ' + rowPool + colPool + ': ' \
							+ str(groupingNeeded) + ': '+ str(groupedCoords)
					
# 							print(outputStr)
						
							groupingNeededInSinglyOccupiedWells += 1
						
						
	print('groupingNeededInSinglyOccupiedWells: '+ str(groupingNeededInSinglyOccupiedWells))

# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateSudokuGridOccupancyTaxonomy(sudokuGridLookupDict, rowPools, colPools, \
prPools, pcPools):

	from scipy import unique

	taxonomyDict = {'emptyWells':0, 'singlyOccupiedWells':0, 'multiplyOccupiedWells':0}

	i = 0
	while i < len(prPools):
		j = 0
		while j < len(pcPools):
		
			k = 0
			while k < len(rowPools):
				l = 0
				while l < len(colPools):
				
					pr = prPools[i]
					pc = pcPools[j]
					row = rowPools[k]
					col = colPools[l]
					
					well = sudokuGridLookupDict[pr][pc].wellGrid[row][col]
					
					readAlignmentCoords = well.readAlignmentCoords
					
					# tempCoords = []
# 					for coord in readAlignmentCoords:
# 						tempCoords.append(int(coord))
# 					
# 					uniqueCoords = unique(tempCoords)
					
					nGenomicCoords = len(readAlignmentCoords)
					
					if nGenomicCoords == 0:
						taxonomyDict['emptyWells'] += 1
					elif nGenomicCoords == 1:
						taxonomyDict['singlyOccupiedWells'] += 1
					elif nGenomicCoords > 1:
						taxonomyDict['multiplyOccupiedWells'] += 1
					
					l += 1
				k += 1
			j += 1
		i += 1


	return taxonomyDict
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def PrintSudokuGridOccupancyTaxonomy(taxonomyDict):
	
	keysInPrintOrder = [\
	'emptyWells',\
	'singlyOccupiedWells',\
	'multiplyOccupiedWells']
	
	outputStr = ''
	
	for key in keysInPrintOrder:
		outputStr += key + ': ' + str(taxonomyDict[key]) + '\n'
	
	print(outputStr)
	return


# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class Feature:

	def __init__(self, coordinates, featureName):
		self.coordinates = coordinates
		
		[self.startCoord, self.endCoord, self.startTranslation, self.endTranslation] \
		= self.calculateStartAndEndCoords(coordinates)
		
		self.featureName = featureName
		self.poolPresenceEntries = []
		self.sudokuGridEntries = []
	
	def calculateStartAndEndCoords(self, coordinates):
		import re
# 		print(coordinates)
		coordsRe = re.compile(r'(\d+)\.\.(\d+)')
		complementRe = re.compile(r'complement\(')
		
		coordsMatch = coordsRe.search(coordinates)
		complementMatch = complementRe.search(coordinates)
		
		if coordsMatch == None:
			print('Error in matching coordinates!')
		else:
			if complementMatch == None:
				startCoord = int(coordsMatch.group(1))
				endCoord = int(coordsMatch.group(2))
				startTranslation = int(coordsMatch.group(1))
				endTranslation = int(coordsMatch.group(2))
			elif complementMatch != None:
				startCoord = int(coordsMatch.group(1))
				endCoord = int(coordsMatch.group(2))
				startTranslation = int(coordsMatch.group(2))
				endTranslation = int(coordsMatch.group(1))
		return [startCoord, endCoord, startTranslation, endTranslation]
	
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class Feature2:
# An updated version of the feature class.
	def __init__(self, featureType, coordinates):
		self.coordinates = coordinates
		self.featureType = featureType
		self.poolPresenceEntries = []
		self.sudokuGridEntries = []
		self.tagDict = {}
		self.cherryPicks = []
		self.featureName = ''
		
		
		self.sudokuWellWithClosestDisruptionToTransStart = None
		self.sudokuWellWithClosestDisruptionToTransStartPicked = False
		
		self.singlyOccupiedSudokuWellWithClosestDisruptionToTransStart = None
		self.singlyOccupiedSudokuWellWithClosestDisruptionToTransStartPicked = False
		
		self.anySingleOccupancySudokuWellPicked = False
		
		
		[self.startCoord, self.endCoord, self.startTranslation, self.endTranslation] \
		= self.calculateStartAndEndCoords(coordinates)
		

	
	def updateTagDict(self, tag, tagData):
		
		tagDictKeys = list(self.tagDict.keys())
		
		if tag not in tagDictKeys:
			self.tagDict[tag] = []
			self.tagDict[tag].append(tagData)
		else:
			self.tagDict[tag].append(tagData)
		
		return
		
	def updateName(self):
		featureTagDictKeys = list(self.tagDict.keys())
	
		if 'gene' in featureTagDictKeys:
			geneName = self.tagDict['gene'][0]
		elif 'note' in featureTagDictKeys:
			geneName = self.tagDict['note'][0]
		elif 'locus_tag' in featureTagDictKeys:
			geneName = self.tagDict['locus_tag'][0]
		else:
			geneName = ''
	
		self.featureName = geneName
		
		return
	
	def calculateStartAndEndCoords(self, coordinates):
		import re
		coordsRe = re.compile(r'(\d+)\.\.(\d+)')
		complementRe = re.compile(r'complement\(')
		joinRe = re.compile(r'join\(')
		joinedCoordsRe = re.compile(r'(\d+)\.\.(\d+)\,(\d+)\.\.(\d+)')
		
		coordsMatch = coordsRe.search(coordinates)
		complementMatch = complementRe.search(coordinates)
		joinMatch = joinRe.search(coordinates)
		joinedCoordsMatch = joinedCoordsRe.search(coordinates)
		
		if (coordsMatch == None) and (joinedCoordsMatch == None):
			print('Error in matching coordinates!')
		else:
			if joinMatch != None and complementMatch == None:
				startCoord = int(joinedCoordsMatch.group(1))
				endCoord = int(joinedCoordsMatch.group(4))
				startTranslation = int(joinedCoordsMatch.group(1))
				endTranslation = int(joinedCoordsMatch.group(4))
			elif joinMatch != None and complementMatch != None:
				startCoord = int(joinedCoordsMatch.group(1))
				endCoord = int(joinedCoordsMatch.group(4))
				startTranslation = int(joinedCoordsMatch.group(4))
				endTranslation = int(joinedCoordsMatch.group(1))
			elif complementMatch == None and joinMatch == None:
				startCoord = int(coordsMatch.group(1))
				endCoord = int(coordsMatch.group(2))
				startTranslation = int(coordsMatch.group(1))
				endTranslation = int(coordsMatch.group(2))
			elif complementMatch != None and joinMatch == None:
				startCoord = int(coordsMatch.group(1))
				endCoord = int(coordsMatch.group(2))
				startTranslation = int(coordsMatch.group(2))
				endTranslation = int(coordsMatch.group(1))
		return [startCoord, endCoord, startTranslation, endTranslation]
	
	def assignDisruptionsClosestToTranslationStart(self):
		
		singlyOccupiedDisrutpionNeeded = False
		
		[closestWell, closestWellPicked] = \
		FindClosestGrowingDisruptionToTranslationStart(self, False)
		
		self.sudokuWellWithClosestDisruptionToTransStart = closestWell
		self.sudokuWellWithClosestDisruptionToTransStartPicked = closestWellPicked
		
		
		[closestWell, closestWellPicked] = \
		FindClosestGrowingDisruptionToTranslationStart(self, True)
		
		self.singlyOccupiedSudokuWellWithClosestDisruptionToTransStart = closestWell
		self.singlyOccupiedSudokuWellWithClosestDisruptionToTransStartPicked = closestWellPicked
								
		return
	
	def findIfAnySingleOccupancyWellHasBeenPicked(self):
		
		anySingleOccupancyPicked = False
		
		i = 0
		while i < len(self.cherryPicks):
			
			if self.cherryPicks[i].multipleOccupancy == False:
				anySingleOccupancyPicked = True
			
			i += 1
		
		self.anySingleOccupancySudokuWellPicked = anySingleOccupancyPicked
		
		return
			
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindClosestGrowingDisruptionToTranslationStart(sudokuFeature, singlyOccupiedDisrutpionNeeded):
	
	minDistanceToTranslationStart = abs(sudokuFeature.startTranslation \
	- sudokuFeature.endTranslation)
	
	closestWell = None
	closestWellPicked = False
	
	for sudokuWellObj in sudokuFeature.sudokuGridEntries:
		
		if singlyOccupiedDisrutpionNeeded == True:
			if len(sudokuWellObj.readAlignmentCoords) == 1:
				evaluateEntryOK = True
			else:
				evaluateEntryOK = True
		else:
			evaluateEntryOK = True
		
		if evaluateEntryOK == True:
			for coordObj in sudokuWellObj.readAlignmentCoords:
				coord = coordObj.coord
				
				distanceToTranslationStart = abs(coord - sudokuFeature.startTranslation)
				
				growthAssumedForWellObj = True
				wellObjPicked = False
				
				for cherryPick in sudokuFeature.cherryPicks:
					if IsCherryPickEquivalentToSudokuWell(cherryPick, sudokuWellObj):
						growthAssumedForWellObj = cherryPick.growthReported
						wellObjPicked = True


				if (distanceToTranslationStart < minDistanceToTranslationStart) \
				and growthAssumedForWellObj == True:
					
					minDistanceToTranslationStart = distanceToTranslationStart
				
					closestWell = sudokuWellObj
					closestWellPicked = wellObjPicked
		
	return [closestWell, closestWellPicked]
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def IsCherryPickEquivalentToSudokuWell(cherryPickObj, sudokuWellObj):
	
	equivalent = False
	
	queryPlate = str(sudokuWellObj.plateName)
	queryRow = str(sudokuWellObj.row)
	queryCol = str(sudokuWellObj.col)
	
	testRow = str(cherryPickObj.sourceRow)
	testCol = str(cherryPickObj.sourceCol)
	testPlate = str(cherryPickObj.sourcePlate)
	
	if (testRow == queryRow) and (testCol == queryCol) and (testPlate == queryPlate):
		equivalent = True
	
	return equivalent
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def FindFeaturesInGenome(genBankFile):
	
	import re
	import pdb
	
	featureNameRe = re.compile(r'\s+/gene=[\'|"](\w+)[\'|"]')
	featureNoteRe = re.compile(r'\s+/note=[\'|"](\w+)[\'|"]')
	
	# Works!
	featureStartRe = re.compile(r'\s+(\w+)\s+([complement\(]*\d+\.\.\d+[\)]*)')

	
	fileHandle = open(genBankFile, 'r')
	
	data = fileHandle.readlines()
	fileHandle.close()
	
# 	pdb.set_trace()
	
	featureArray = []
	
	i = 0
	currentCoordinates = None
	currentName = None
	currentNote = None
	inInterestingRecord = False
	
	
	while i < len(data):
		line = data[i]
		
		featureStartMatch = featureStartRe.search(line)

		if featureStartMatch != None:
			if currentCoordinates != None and inInterestingRecord == True:
				if currentName == None:
					if currentNote != None:
						currentName = currentNote
					else:
						currentName = ''
					
				featureToStore = Feature(currentCoordinates, currentName)
				# print(featureToStore.featureName + ': ' + str(featureToStore.startCoord) \
# 				+ '..' + str(featureToStore.endCoord))
				featureArray.append(Feature(currentCoordinates, currentName))
			
			featureType = featureStartMatch.group(1)
			if featureType == 'gene':
# 				print(featureType)
				currentCoordinates = featureStartMatch.group(2)
# 				print(currentCoordinates)
				currentName = None
				inInterestingRecord = True
				currentNote = None
			else:
				currentCoordinates = None
				currentName = None
				inInterestingRecord = False

		if inInterestingRecord == True:
			featureNameMatch = featureNameRe.search(line)
		
			if featureNameMatch != None:
				currentName = featureNameMatch.group(1)

			featureNoteMatch = featureNoteRe.search(line)
			if featureNoteMatch != None:
				currentNote = featureNoteMatch.group(1)
# 				print(currentNote)
				
		i += 1
		
	
	return featureArray
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def FindFeaturesInGenome2(genBankFile):
# A heavily updated version of FindFeaturesInGenome. Can deal with all sorts of tags, rather than
# just genes and can deal with joined CDSs. 

	import re
	import pdb
	
	featureNameRe = re.compile(r'\s+/gene=[\'|"](\w+)[\'|"]')
	
# 	featureStartRe = re.compile(r'\s+(\w+)\s+([complement\(]*\d+\.\.\d+[\)]*)')
	featureStartRe = re.compile(r'\s+(\w+)\s+([complement\(|join\(]*\d+\.\.\d+[\,\d+\.\.\d+]*[\)]*)')
	featureTagRe = \
	re.compile(r'\s+/([a-zA-Z0-9_\-]+)=([\'|"]*)([\[a\]-zA-Z0-9_\-\.\,\: \t\&\*\(\)/;\+\/`\']+)([\'|"]*)')
	featureTagRe2 = \
	re.compile(r'\s+/([a-zA-Z0-9_\-]+)')
	
	sequenceSectionStartRe = re.compile(r'\s*ORIGIN')
	featureSectionStartRe = re.compile(r'\s*FEATURES')

	
	fileHandle = open(genBankFile, 'r')
	
	data = fileHandle.readlines()
	fileHandle.close()
	
	
	featureArray = []
	
	i = 0
	currentCoordinates = None
	currentFeature = None
	
	inFeatureEntry = False
	
	inFeatureSection = False
	inSequenceSection = False
	
	continuationNeeded = False
	
	
	while i < len(data):
			
		line = data[i]
		
		# Start looking for the start of a gene record
		
		featureSectionStartMatch = featureSectionStartRe.match(line)
		sequenceSectionStartMatch = sequenceSectionStartRe.match(line)
		
		if featureSectionStartMatch != None:
			inFeatureSection = True
			inSequenceSection = False
			inFeatureEntry = False
			
		if sequenceSectionStartMatch != None:
			inFeatureSection = False
			inSequenceSection = True
			inFeatureEntry = False
			
			if currentFeature != None:
				featureArray.append(currentFeature)
			
		
		if inFeatureSection and not inSequenceSection:
			
			featureStartMatch = featureStartRe.search(line)
					
			if (featureStartMatch != None):
				if currentFeature != None:
					featureArray.append(currentFeature)
				
				# Now, start getting the new feature
				inInterestingRecord = True
				currentFeatureType = featureStartMatch.group(1)
				currentCoordinates = featureStartMatch.group(2)
			
				# Initialize a new current feature
				currentFeature = Feature2(currentFeatureType, currentCoordinates)
				inFeatureEntry = True
				
			else:
				if inFeatureEntry == True:
					
					featureTagMatch = featureTagRe.search(line)
					featureTagMatch2 = featureTagRe2.match(line)
					
					if (featureTagMatch != None):
						
						if continuationNeeded == True:
							print("continuation error step 1!")
							return
						
						
						currentTag = featureTagMatch.group(1)
						currentTagData = featureTagMatch.group(3).strip()
						
						currentStartQuote = featureTagMatch.group(2).strip()
						currentEndQuote = featureTagMatch.group(4).strip()
						
						# If the tag data starts with a quote but doesn't end with a quote, 
						# then continue onto the next line
						if ((currentStartQuote == '"') or (currentStartQuote == '\'')) \
						and (currentEndQuote == ''):
							continuationNeeded = True
					
					if (featureTagMatch2 != None) and (featureTagMatch == None) \
					and (continuationNeeded == False):
						# print(featureTagMatch2.group(1))
# 						print(currentCoordinates)
						currentTag = featureTagMatch2.group(1)
						currentTagData = True
					
							
					if continuationNeeded == True and featureTagMatch == None:
						currentTagData += line.strip()
						currentTagEnding = currentTagData[-1]
						
						if ((currentTagEnding == '"') or (currentTagEnding == '\'')):
							currentTagData = currentTagData[0:-1]
							continuationNeeded = False
						
					if continuationNeeded == False:
# 						print(currentTag + str(currentTagData))
						currentFeature.updateTagDict(currentTag, currentTagData)
				
		i += 1
	
	for feature in featureArray:
		feature.updateName()
		
	
	return featureArray
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def FindFeaturesWithDisruptions(featureArray): 
	i = 0

	featuresWithDisruptions = []

	while i < len(featureArray):
		if len(featureArray[i].poolPresenceEntries) > 0:
			featureName = featureArray[i].featureName
		
			startCoord = featureArray[i].startCoord
			endCoord = featureArray[i].endCoord
		
			startTranslation = featureArray[i].startTranslation
			endTranslation = featureArray[i].endTranslation
		
			featuresWithDisruptions.append(featureArray[i])
		
		i += 1

	return featuresWithDisruptions
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def FindWellsWithDisruptionBetweenStartAndEndCoords(startCoord, endCoord, sudokuGridLookupDict, \
rowPools, colPools, prPools, pcPools):
	
	i = 0
	
	wellsToReturn = []
	
	while i < len(prPools):
		j = 0
		while j < len(pcPools):
		
			k = 0
			while k < len(rowPools):
				l = 0
				while l < len(colPools):
				
					pr = prPools[i]
					pc = pcPools[j]
					row = rowPools[k]
					col = colPools[l]
					
					well = sudokuGridLookupDict[pr][pc].wellGrid[row][col]
					
					readAlignmentCoordsArray = well.readAlignmentCoords
					
					tempWellsToReturn = []
					
					# For now, I'll take even ambiguous or guessed coordinates, but in the future
					# might only take un-ambiguous ones. 
					for sudokuCoord in readAlignmentCoordsArray:
						
						coord = sudokuCoord.coord
						
						if startCoord < coord < endCoord:
							tempWellsToReturn.append({'plateName':well.plateName, \
							'multipleOccupancy':False, \
							'readAlignmentCoord':coord, 'row':row, 'col':col, \
							'plateRow':pr, 'plateCol':pc})
					
					if len(readAlignmentCoordsArray) > 1:
						for tempWell in tempWellsToReturn:
							tempWell['multipleOccupancy'] = True
					
					for tempWell in tempWellsToReturn:
						wellsToReturn.append(tempWell)
					
					l += 1
				k += 1
			j += 1
		i += 1
		
	
	return wellsToReturn
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def UpdateFeatureArrayWithSudokuGridEntries(featureArray, sudokuGridLookupDict, \
rowPools, colPools, prPools, pcPools):

	import pdb

	i = 0 
	while i < len(featureArray):
		feature = featureArray[i]		
		startCoord = feature.startCoord
		endCoord = feature.endCoord
		
		wellsWithDisruptionsMatchingFeature = \
		FindWellsWithDisruptionBetweenStartAndEndCoords(startCoord, endCoord, \
		sudokuGridLookupDict, rowPools, colPools, prPools, pcPools)
		
		for well in wellsWithDisruptionsMatchingFeature:
			feature.sudokuGridEntries.append(well)
		
		i += 1
			
		
	return featureArray
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindWellsWithDisruptionBetweenStartAndEndCoords2(startCoord, endCoord, sudokuGridLookupDict, \
rowPools, colPools, prPools, pcPools):
	
	i = 0
	
	wellsToReturn = []
	
	while i < len(prPools):
		j = 0
		while j < len(pcPools):
		
			k = 0
			while k < len(rowPools):
				l = 0
				while l < len(colPools):
				
					pr = prPools[i]
					pc = pcPools[j]
					row = rowPools[k]
					col = colPools[l]
					
					well = sudokuGridLookupDict[pr][pc].wellGrid[row][col]
					readAlignmentCoordsArray = well.readAlignmentCoords
					
					wellContainsDesiredMutant = False
					
					for sudokuCoord in readAlignmentCoordsArray:
						coord = sudokuCoord.coord
						if startCoord < coord < endCoord:
							wellContainsDesiredMutant = True
					
					if wellContainsDesiredMutant:
						wellsToReturn.append(well)
					
					l += 1
				k += 1
			j += 1
		i += 1
		
	
	return wellsToReturn
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def UpdateFeatureArrayWithSudokuGridEntries2(featureArray, sudokuGridLookupDict, \
rowPools, colPools, prPools, pcPools):

	import pdb
	from copy import deepcopy

	i = 0 
	while i < len(featureArray):
		feature = featureArray[i]		
		startCoord = feature.startCoord
		endCoord = feature.endCoord
		
		wellsWithDisruptionsMatchingFeature = \
		FindWellsWithDisruptionBetweenStartAndEndCoords2(startCoord, endCoord, \
		sudokuGridLookupDict, rowPools, colPools, prPools, pcPools)
		
		for well in wellsWithDisruptionsMatchingFeature:
			feature.sudokuGridEntries.append(deepcopy(well))
			
					
		i += 1
			
		
	return featureArray
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CalculateFeatureArrayTaxonomy2(featureArray):
	
	keysInPrintOrder = [\
	'totalFeatures',\
	'featuresWithAtLeastOneDisruption',\
	'featuresWithNoDisruptions',\
	'featuresWithNoDisruptionsInSinglyOccupiedWells', \
	'featuresWithDisruptionsInSinglyOccupiedWells']
	
	taxonomyDict = {}
	
	for key in keysInPrintOrder:
		taxonomyDict[key] = 0
	
	for feature in featureArray:
		taxonomyDict['totalFeatures'] += 1
		
		if len(feature.sudokuGridEntries) == 0:
			taxonomyDict['featuresWithNoDisruptions'] += 1
		
		if len(feature.sudokuGridEntries) > 0:
			taxonomyDict['featuresWithAtLeastOneDisruption'] += 1
			
			entriesInSinglyOccupiedWells = 0
			
			for entry in feature.sudokuGridEntries:
				if len(entry.readAlignmentCoords) > 1:
					entriesInSinglyOccupiedWells += 1
			
			if entriesInSinglyOccupiedWells == 0:
				taxonomyDict['featuresWithNoDisruptionsInSinglyOccupiedWells'] += 1
			elif entriesInSinglyOccupiedWells > 0:
				taxonomyDict['featuresWithDisruptionsInSinglyOccupiedWells'] += 1
	
	
	PrintFeatureArrayTaxonomyDict(taxonomyDict, keysInPrintOrder)
	
	return taxonomyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateFeatureArrayTaxonomy(featureArray):
	
	keysInPrintOrder = [\
	'totalFeatures',\
	'featuresWithAtLeastOneDisruption',\
	'featuresWithNoDisruptions',\
	'featuresWithNoDisruptionsInSinglyOccupiedWells', \
	'featuresWithDisruptionsInSinglyOccupiedWells']
	
	taxonomyDict = {}
	
	for key in keysInPrintOrder:
		taxonomyDict[key] = 0
	
	for feature in featureArray:
		taxonomyDict['totalFeatures'] += 1
		
		if len(feature.sudokuGridEntries) == 0:
			taxonomyDict['featuresWithNoDisruptions'] += 1
		
		if len(feature.sudokuGridEntries) > 0:
			taxonomyDict['featuresWithAtLeastOneDisruption'] += 1
			
			entriesInSinglyOccupiedWells = 0
			
			for entry in feature.sudokuGridEntries:
				if entry['multipleOccupancy'] == False:
					entriesInSinglyOccupiedWells += 1
			
			if entriesInSinglyOccupiedWells == 0:
				taxonomyDict['featuresWithNoDisruptionsInSinglyOccupiedWells'] += 1
			elif entriesInSinglyOccupiedWells > 0:
				taxonomyDict['featuresWithDisruptionsInSinglyOccupiedWells'] += 1
	
	
	PrintFeatureArrayTaxonomyDict(taxonomyDict, keysInPrintOrder)
	
	return taxonomyDict
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def PrintFeatureArrayTaxonomyDict(taxonomyDict, keysInPrintOrder):
	
	outputStr = ''
	
	for key in keysInPrintOrder:
		outputStr += key + ': ' + str(taxonomyDict[key]) + '\n'
	
	print(outputStr)
	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ReturnClosestDisruptionToTranslationStartInSinglyOccupiedWell(feature):
	
	import pdb
	
# 	pdb.set_trace()
	
	closestDisruption = None
	
	entriesInSinglyOccupiedWells = 0
	
	for entry in feature.sudokuGridEntries:
		if entry['multipleOccupancy'] == False:
			entriesInSinglyOccupiedWells += 1
	
	if len(feature.sudokuGridEntries) > 0 and entriesInSinglyOccupiedWells > 0:
		i = 0
		closestI = 0
		
		end = feature.endTranslation
		start = feature.startTranslation
		distanceOfClosestIToTranslationStart = abs(end - start)
		
		while i < len(feature.sudokuGridEntries):
			readAlignmentCoord = feature.sudokuGridEntries[i]['readAlignmentCoord']
			startTranslation = feature.startTranslation 
			distanceToTranslationStart = abs(readAlignmentCoord - startTranslation)
			multipleOccupancy = feature.sudokuGridEntries[i]['multipleOccupancy']
			
			if (distanceToTranslationStart < distanceOfClosestIToTranslationStart) \
			and (multipleOccupancy == False):
				closestI = i
				distanceOfClosestIToTranslationStart = distanceToTranslationStart
			
			i += 1
		
		closestDisruption = feature.sudokuGridEntries[closestI]
		
	return closestDisruption
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ConvertWellIDToRowAndColumn(well):

	import re

	wellIDRe = re.compile(r'([A-H])(\d+)', re.IGNORECASE)
	
	wellMatch = wellIDRe.match(well)
	
	row = wellMatch.group(1)
	col = wellMatch.group(2)
	
	return [row, col]
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def IncrementWellInColumnsOrder(previousWell, previousPlate):
	
	import re
	
	
	gRows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	gColumns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
		
	wellIDRe = re.compile(r'([A-H])(\d+)', re.IGNORECASE)
	
	wellMatch = wellIDRe.search(previousWell)
	
	row = wellMatch.group(1)
	col = wellMatch.group(2)
	
	rowIndex = gRows.index(row)
	colIndex = gColumns.index(col)
	
	if rowIndex == len(gRows) - 1:
		nextRowIndex = 0
		nextColIndex = colIndex+1
	else:
		nextRowIndex = rowIndex + 1
		nextRow = gRows[nextRowIndex]
		nextColIndex = colIndex
	
	
	if nextColIndex == len(gColumns):
		nextPlate = previousPlate + 1
		nextColIndex = 0
	else:
		nextPlate = previousPlate
	
	nextWell = gRows[nextRowIndex] + gColumns[nextColIndex].zfill(2)
		
	return [nextWell, nextPlate]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def IncrementWellInRowsOrder(previousWell, previousPlate):
	
	import re
	
	
	gRows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	gColumns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
	
	wellIDRe = re.compile(r'([A-H])(\d+)', re.IGNORECASE)
	
	wellMatch = wellIDRe.search(previousWell)
	
	row = wellMatch.group(1)
	col = wellMatch.group(2)
	
	rowIndex = gRows.index(row)
	colIndex = gColumns.index(col)
	
	if colIndex == len(gColumns) - 1:
		nextColIndex = 0
		nextRowIndex = rowIndex+1
	else:
		nextRowIndex = rowIndex
		nextColIndex = colIndex + 1
		nextCol = gColumns[nextColIndex]
		
	
	if nextRowIndex == len(gRows):
		nextPlate = previousPlate + 1
		nextRowIndex = 0
	else:
		nextPlate = previousPlate
	
	nextWell = gRows[nextRowIndex] + gColumns[nextColIndex].zfill(2)
		
	return [nextWell, nextPlate]
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def GetNextDestinationWell(previousWell, previousPlate, controlWellsDict, fillOrder='columns'):
	
	import re
	
	nextDestinationWellOK = False
	
	while nextDestinationWellOK == False:
	
		if fillOrder == 'columns':
			[nextWell, nextPlate] = IncrementWellInColumnsOrder(previousWell, previousPlate)
		elif fillOrder == 'rows':
			[nextWell, nextPlate] = IncrementWellInRowsOrder(previousWell, previousPlate)
		
		controlWells = controlWellsDict[nextPlate]
# 		print(controlWells)
		
		if nextWell in controlWells:
			nextDestinationWellOK = False
			previousWell = nextWell
			previousPlate = nextPlate
		else:
			nextDestinationWellOK = True
	
	
	return [nextWell, nextPlate]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateControlWellsDict(numberPlates):
	controlWellsDict = {}
	
	if numberPlates > 96:
		print("Can not deal with more than 96 plates!")
		return 
	
	i = 1
	controlWellsDict[1] = 'A01'
	i = 2
	previousWell = 'A01'
	previousPlate = 1
	while i <= numberPlates:
		[nextWell, nextPlate] = IncrementWellInRowsOrder(previousWell, previousPlate)
		controlWellsDict[i] = nextWell
		previousWell = nextWell
		i += 1
	
	return controlWellsDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ensure_dir(f):
    import os
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
# ------------------------------------------------------------------------------------------------ #







####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 0: Initialize folders
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ----------------------------------------------------------------------------------------------- #
def InitializeOutputFoldersAndFiles(outputLog, himarRecognitionBaseDir, indexSummaryBaseDir, \
trimmedSequencesBaseDir, genomeAlignmentBaseDir, genomeArrayFileName, poolFileBaseDir, \
poolPresenceTableFileName):
	
	import os.path
	
	import pdb
	
	ensure_dir(himarRecognitionBaseDir)
	ensure_dir(indexSummaryBaseDir)
	ensure_dir(genomeAlignmentBaseDir)
	ensure_dir(poolFileBaseDir)
	ensure_dir(trimmedSequencesBaseDir)
	
	outputLogDir = os.path.dirname(outputLog)
	genomeArrayDir = os.path.dirname(genomeArrayFileName)
	poolPresenceTableDir = os.path.dirname(poolPresenceTableFileName)
	
# 	pdb.set_trace()
	
	ensure_dir(outputLogDir + '/')
	ensure_dir(genomeArrayDir + '/')
	ensure_dir(poolPresenceTableDir + '/')
	
	fileHandle = open(outputLog, 'w')
	fileHandle.close()

# ----------------------------------------------------------------------------------------------- #






####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 1: Functions for Checking Himar Content and Trimming Sequence Reads
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ----------------------------------------------------------------------------------------------- #
def Mismatch(Search,n):
	import itertools
	SearchL = list(Search)
	List     = [] # hold output
	# print list of indices to replace with '$'
	idxs = itertools.combinations(range(len(SearchL)),n)        
	# for each combination `idx` in idxs, replace str[idx] with '$':
	for idx in idxs:
		str = SearchL[:] # make a copy
		for i in idx:
			str[i]='[ACTGN]'
		List.append( ''.join(str) ) # convert back to string
	return List
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GenerateHimarMismatchRegexes(himarSequence, mismatchNumber=2):
	import re
	
	himarMismatchRegexList = []
	fullSequenceHimarMismatchRegexList = []
	
	searchList = Mismatch(himarSequence,mismatchNumber)
	searchListRegexes = []
	for searchTerm in searchList:
		himarMismatchRegexList.append(re.compile(r'' + searchTerm))
		fullSequenceHimarMismatchRegexList.append(re.compile(r'([NACTG]{0,8})' + r'(' + searchTerm + r')' + r'([ACTGN]+)'))
	
	return [himarMismatchRegexList, fullSequenceHimarMismatchRegexList]
	
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def CheckAndTrimHimarReadSequence(readSeq, idString, himarRe, himarMismatchRegexList, \
fullSequenceHimarMismatchRegexList):
	
	from pdb import set_trace
	
	himarReMatch = himarRe.match(readSeq)
	
		
	if himarReMatch != None:
		himarAssignment = True
		trimmedSequence = himarReMatch.group(3)
		assignmentType = 'perfect'
	else:
		himarMismatchMatches = 0
		
		trimmedSequence = ''
		
		i = 0
		while i < len(himarMismatchRegexList):
			regex = himarMismatchRegexList[i]
			if regex.search(readSeq) != None:
				himarMismatchMatches += 1
				fullSequenceRe = fullSequenceHimarMismatchRegexList[i].search(readSeq)
				if fullSequenceRe == None:
# 					set_trace()
					trimmedSequence = ''
				else:
					trimmedSequence = fullSequenceRe.group(3)
			i += 1
		
		if trimmedSequence == '':
			himarAssignment = False
			assignmentType = 'unmatchable'
			trimmedSequence = readSeq
		elif himarMismatchMatches >= 1:
			himarAssignment = True
			assignmentType = 'matchable'	
		else:
			himarAssignment = False
			assignmentType = 'unmatchable'
			trimmedSequence = readSeq
		
		
	return [himarAssignment, assignmentType, trimmedSequence]
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def ParseHimarReadFile(file, himarMismatchRegexList, fullSequenceHimarMismatchRegexList, himarRe):
	# Import the fastq files and summarize
	
	from pdb import set_trace
	import re
	
	# We'll assume the following format for the read id

	# @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	# EAS139 is the unique instrument name
	# 136 is the run id
	# FC706VJ is the flowcell id
	# 2	is the flowcell lane
	# 2104 is the tile number within the flowcell lane
	# 15343	is the 'x'-coordinate of the cluster within the tile
	# 197393 is the 'y'-coordinate of the cluster within the tile
	# 1	is the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
	# Y	is Y if the read is filtered, N otherwise
	# 18i is 0 when none of the control bits are on, otherwise it is an even number
	# ATCACG index sequence

	# '@JLK5VL1:704:H3M72BCXX:2:1101:1209:2029 2:N:0:'
	identifierLineRe = re.compile('@(\w+):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)')

	# In these dicts, the key is the read number, and the associated value is true or false for 
	# himar recognition

	himarRecognitionDict = {}
	perfectHimarReads = 0
	imperfectButMatchableHimarReads = 0
	unMatchableHimarReads = 0
	himarTrimmedSequenceDict = {}

	
	fileHandle = open(file, 'r')
	fileData = fileHandle.readlines()
	
	i = 0
	while i < len(fileData):
		line = fileData[i]
		
		identifierMatch = identifierLineRe.match(line)
		
		if identifierMatch != None:
				
			lane = identifierMatch.group(4)
			tileNumber = identifierMatch.group(5)
			xCoord = identifierMatch.group(6)
			yCoord = identifierMatch.group(7)
			idString = lane + ':' + tileNumber + ':' + xCoord + ':' + yCoord			
			readSeq = fileData[i+1].strip()
			
			
			[himarAssignment, assignmentType, trimmedSequence] = \
			CheckAndTrimHimarReadSequence(readSeq, idString, himarRe, himarMismatchRegexList, \
			fullSequenceHimarMismatchRegexList)
			
			himarTrimmedSequenceDict[idString] = trimmedSequence
			himarRecognitionDict[idString] = himarAssignment
			
			if assignmentType == 'perfect':
				perfectHimarReads += 1
			elif assignmentType == 'matchable':
				imperfectButMatchableHimarReads += 1
			elif assignmentType == 'unmatchable':
				unMatchableHimarReads += 1
			
		i += 1

	return [himarRecognitionDict, himarTrimmedSequenceDict, perfectHimarReads, \
	imperfectButMatchableHimarReads, unMatchableHimarReads]
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def WriteHimarSequenceSummaries(himarRecognitionOutputFileName, himarRecognitionDict):
	
	summaryFileHandle = open(himarRecognitionOutputFileName, 'w')
	
	himarRecognitionDictKeys = himarRecognitionDict.keys()
	
	for key in himarRecognitionDictKeys:		
		summaryFileHandle.write(str(key) + ',' + str(himarRecognitionDict[key]) + '\n')

	summaryFileHandle.close()
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def WriteTrimmedSequences(trimmedFileName, himarTrimmedSequenceDict):
	
	mfastaHandle = open(trimmedFileName, 'w')
	
	himarSequenceDictKeys = himarTrimmedSequenceDict.keys()
	
	for key in himarSequenceDictKeys:
		mfastaHandle.write('>' + str(key) + '\n')
		mfastaHandle.write(str(himarTrimmedSequenceDict[key]) + '\n')

	mfastaHandle.close()
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GenerateLogStringForIndividualFastQFileForHimarRecognition(startTime, endTime, \
currentFileNumber, totalFiles, perfectIndexReads, mismatchedButMatchableIndexReads, \
unMatchableIndexReads ):

	duration = (endTime - startTime).total_seconds()
	
	indexReadSummary = "Perfect Himar Reads: " + str(perfectIndexReads) + '\n'
	indexReadSummary += "Mismatched But Readable Himar Reads: " + \
	str(mismatchedButMatchableIndexReads) + '\n'
	indexReadSummary += "Unmatchable Himar Reads: " + str(unMatchableIndexReads) + '\n'

	lastOperationDurationStr = "Last operation duration: " + str(duration) + " seconds" + '\n'
	
	timeRemainingEstimate = duration*(totalFiles-currentFileNumber)/60
	timeRemainingStr = "Time remaining estimation: " + str(timeRemainingEstimate) + " minutes"  \
	+ '\n'
	
	logStr = indexReadSummary + lastOperationDurationStr + timeRemainingStr

	return logStr
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
def GenerateLogFileSummaryForHimarRecognition(totalPerfectHimarReads, \
totalImperfectButMatchableHimarReads, totalUnMatchableHimarReads, initialStartTime, \
completeEndTime):

	jobDuration = (completeEndTime - initialStartTime).total_seconds()

	summary = gSeparatorString
	summary += "Total Perfect Himar Reads: " + str(totalPerfectHimarReads) + '\n'
	summary += "Total Mismatched But Readable Himar Reads: " + \
	str(totalImperfectButMatchableHimarReads) + '\n'
	summary += "Total Unmatchable Himar Reads: " + str(totalUnMatchableHimarReads) + '\n'
	
	summary += "Complete Job Time: " + str(jobDuration/60) + " minutes\n"
	summary += gSeparatorString
	
	return summary
# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
def CheckHimarContentAndTrimSequences(fastqFiles, \
himarSequence, himarRecognitionFilePrefix, himarRecognitionBaseDir, \
trimmedSequencesFilePrefix, trimmedSequencesBaseDir, outputLog):

# Import the fastq files, check the himar and summarize the himar content and trim the sequence
# to only the genomic material

# In these dicts, the key is the read number, and the associated value is true or false for himar
# recognition

# Note, as of the current revision this is a bit of kludgy hold over from earlier versions of the 
# code. This function writes out the summary files and returns the file names, only for the 
# next function in the analysis chain to read them back in. 

	import datetime
	import re
	
	# Generate himar regular expressions and mismatch regular expressions
	himarRe = re.compile(r'([NACTG]{0,8})' + r'(' + himarSequence + r')' + r'([ACTGN]+)')

	[himarMismatchRegexList, fullSequenceHimarMismatchRegexList] = \
	GenerateHimarMismatchRegexes(himarSequence, mismatchNumber=2)


	logFileHandle = open(outputLog, 'a')
	initialStartTime = datetime.datetime.now()

	j = 1
	totalPerfectHimarReads = 0
	totalImperfectButMatchableHimarReads = 0
	totalUnMatchableHimarReads = 0
	
	himarRecognitionFiles = []
	trimmedMfaFiles = []

	for file in fastqFiles:

		startTime = datetime.datetime.now()
	
		# Parse the individual himar read fastq file
		processStr = 'Processing ' + str(j) + ' of ' + str(len(fastqFiles)) + ': ' + file + '\n'
		print(processStr)
		
		[himarRecognitionDict, himarTrimmedSequenceDict, perfectHimarReads, \
		imperfectButMatchableHimarReads, unMatchableHimarReads] \
		= ParseHimarReadFile(file, himarMismatchRegexList, fullSequenceHimarMismatchRegexList, \
		himarRe)
	
	
		totalPerfectHimarReads += perfectHimarReads
		totalImperfectButMatchableHimarReads += imperfectButMatchableHimarReads
		totalUnMatchableHimarReads += unMatchableHimarReads

		# Write out the himar sequence summaries
		himarRecognitionOutputFileName = himarRecognitionBaseDir + '/' \
		+ himarRecognitionFilePrefix + '_' + str(j) + '.csv'
		
		WriteHimarSequenceSummaries(himarRecognitionOutputFileName, himarRecognitionDict)
		
		himarRecognitionFiles.append(himarRecognitionOutputFileName)
	
		# Write out the trimmed sequences
		trimmedFileName = trimmedSequencesBaseDir + '/' \
		+ trimmedSequencesFilePrefix + '_' + str(j) + '.mfa'
		
		WriteTrimmedSequences(trimmedFileName, himarTrimmedSequenceDict)
		
		trimmedMfaFiles.append(trimmedFileName)

		endTime = datetime.datetime.now()
		
		summaryLogStr = GenerateLogStringForIndividualFastQFileForHimarRecognition(startTime, \
		endTime, j, len(fastqFiles), perfectHimarReads, imperfectButMatchableHimarReads, \
		unMatchableHimarReads)
	
		print(summaryLogStr)
	
		# Write data to log file
		UpdateLogFileData(outputLog, processStr+summaryLogStr)

	
		j += 1
	
	completeEndTime = datetime.datetime.now()

	summaryStr = GenerateLogFileSummaryForHimarRecognition(totalPerfectHimarReads, \
	totalImperfectButMatchableHimarReads, \
	totalUnMatchableHimarReads, initialStartTime, completeEndTime)

	print(summaryStr)

	UpdateLogFileData(outputLog, summaryStr)
	
	return [himarRecognitionFiles, trimmedMfaFiles]
# ------------------------------------------------------------------------------------------------ #





####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 2: Functions for Index Summarization
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ----------------------------------------------------------------------------------------------- #
def AssignBarcodeReadToPool(indexSeq, mismatchRegexListDict, barcodeLookupTable, barcodeList):
	
	from scipy import unique
	from pdb import set_trace
	
	
	if indexSeq in barcodeList:
		barcodeAssignment = [barcodeLookupTable[indexSeq], True]
		assignmentType = 'perfect'
	else:	
		mismatchedSequencePossibleMatches = []
	
		for barcode in barcodeList:
			searchListRegexes = mismatchRegexListDict[barcode]
			for regex in searchListRegexes:
				if regex.search(indexSeq) != None:
					mismatchedSequencePossibleMatches.append(barcodeLookupTable[barcode])
				
		uniqueMismatchedSequencePossibleMatches = unique(mismatchedSequencePossibleMatches)
			
		if len(uniqueMismatchedSequencePossibleMatches) == 1:
			barcodeAssignment = [mismatchedSequencePossibleMatches[0], True]
			assignmentType = 'matchable'
		else:
			barcodeAssignment = [indexSeq, False]
			assignmentType = 'unmatchable'
	

	return [barcodeAssignment, assignmentType]
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def GenerateIndexMismatchRegexes(barcodeList, mismatchNumber=2):
	
	import re
	
	mismatchRegexListDict = {}
	
	for barcode in barcodeList:
	
		searchList = Mismatch(barcode,mismatchNumber)
		searchListRegexes = []
	
		for searchTerm in searchList:
			searchListRegexes.append(re.compile(r'' + searchTerm))

		mismatchRegexListDict[barcode] = searchListRegexes
	
	return mismatchRegexListDict
	
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def GenerateLogFileSummaryForIndexSummary(totalPerfectIndexReads, \
totalMismatchedButMatchableIndexReads, \
totalUnMatchableIndexReads, initialStartTime, completeEndTime):

	jobDuration = (completeEndTime - initialStartTime).total_seconds()

	summary = gSeparatorString
	summary += "Total Perfect Index Reads: " + str(totalPerfectIndexReads) + '\n'
	summary += "Total Mismatched But Readable Index Reads: " + \
	str(totalMismatchedButMatchableIndexReads) + '\n'
	summary += "Total Unmatchable Index Reads: " + str(totalUnMatchableIndexReads) + '\n'
	
	summary += "Complete Job Time: " + str(jobDuration/60) + " minutes\n"
	summary += gSeparatorString
	
	return summary
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def GenerateLogStringForIndividualFastQFileForIndexSummary(startTime, endTime, currentFileNumber, \
totalFiles, perfectIndexReads, mismatchedButMatchableIndexReads, unMatchableIndexReads ):

	duration = (endTime - startTime).total_seconds()
	
	indexReadSummary = "Perfect Index Reads: " + str(perfectIndexReads) + '\n'
	indexReadSummary += "Mismatched But Readable Index Reads: " + \
	str(mismatchedButMatchableIndexReads) + '\n'
	indexReadSummary += "Unmatchable Index Reads: " + str(unMatchableIndexReads) + '\n'

	lastOperationDurationStr = "Last operation duration: " + str(duration) + " seconds" + '\n'
	
	timeRemainingEstimate = duration*(totalFiles-currentFileNumber)/60
	timeRemainingStr = "Time remaining estimation: " + str(timeRemainingEstimate) + " minutes"  + '\n'
	
	logStr = indexReadSummary + lastOperationDurationStr + timeRemainingStr

	return logStr
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def ParseIndexReadFile(file, mismatchRegexListDict, barcodeLookupTable, barcodeList):
	
	import re
	
	# We'll assume the following format for the read id

	# @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	# EAS139 is the unique instrument name
	# 136 is the run id
	# FC706VJ is the flowcell id
	# 2	is the flowcell lane
	# 2104 is the tile number within the flowcell lane
	# 15343	is the 'x'-coordinate of the cluster within the tile
	# 197393 is the 'y'-coordinate of the cluster within the tile
	# 1	is the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
	# Y	is Y if the read is filtered, N otherwise
	# 18i is 0 when none of the control bits are on, otherwise it is an even number
	# ATCACG index sequence

	# '@JLK5VL1:704:H3M72BCXX:2:1101:1209:2029 2:N:0:'
	identifierLineRe = re.compile('@(\w+):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)')
	
	barcodeRecognitionDict = {}
	
	fileHandle = open(file, 'r')
	fileData = fileHandle.readlines()
	
	perfectIndexReads = 0
	mismatchedButMatchableIndexReads = 0
	unMatchableIndexReads = 0
	
	i = 0
	while i < len(fileData):
		line = fileData[i]
		
		identifierMatch = identifierLineRe.match(line)
		
		if identifierMatch != None:
				
			lane = identifierMatch.group(4)
			tileNumber = identifierMatch.group(5)
			xCoord = identifierMatch.group(6)
			yCoord = identifierMatch.group(7)
			idString = lane + ':' + tileNumber + ':' + xCoord + ':' + yCoord
			
			indexSeq = fileData[i+1].strip()
			
			[barcodeAssignment, assignmentType] \
			= AssignBarcodeReadToPool(indexSeq, mismatchRegexListDict, barcodeLookupTable, barcodeList)
			
			barcodeRecognitionDict[idString] = barcodeAssignment
			
			if assignmentType == 'perfect':
				perfectIndexReads += 1
			elif assignmentType == 'matchable':
				mismatchedButMatchableIndexReads += 1
			elif assignmentType == 'unmatchable':
				unMatchableIndexReads += 1
				
		i += 1

	return [barcodeRecognitionDict, perfectIndexReads, mismatchedButMatchableIndexReads, \
	unMatchableIndexReads]
# ----------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------- #
def OutputBarcodeRecognitionTable(barcodeRecognitionDict, barcodeRecognitionOutputFileName):

	barcodeRecognitionKeys = list(barcodeRecognitionDict.keys())

	fileHandle = open(barcodeRecognitionOutputFileName, 'w')

	for key in barcodeRecognitionKeys:
		
		outputStr = str(key) + ',' + str(barcodeRecognitionDict[key][0]) \
		+ ',' + str(barcodeRecognitionDict[key][1]) + '\n'
	
		fileHandle.write(outputStr)
	
	fileHandle.close()
# ----------------------------------------------------------------------------------------------- #




# ------------------------------------------------------------------------------------------------ #
def SummarizeIndexReads(indexFastqFiles, barcodeFile, indexSummaryFilePrefix, indexSummaryBaseDir, \
outputLog):
# Import the fastq files and summarize

# In these dicts, the key is the read number, and the associated value is the code for the pool
# barcode

	import datetime
	
	
	# Generate the barcode lookup table and  the mismatch regexes

	[barcodeLookupTable, barcodes] = GenerateBarcodeLookupTable(barcodeFile)

	barcodeList = list(barcodes)
	mismatchRegexListDict = GenerateIndexMismatchRegexes(barcodeList, mismatchNumber=2)

	
	

	j = 1
	totalPerfectIndexReads = 0
	totalMismatchedButMatchableIndexReads = 0
	totalUnMatchableIndexReads = 0
	
	indexAlignmentFiles = []

	logFileHandle = open(outputLog, 'a')

	initialStartTime = datetime.datetime.now()

	for file in indexFastqFiles:
	
		startTime = datetime.datetime.now()
	
		# Parse the individual index read fastq file
		processStr = 'Processing ' + str(j) + ' of ' + str(len(indexFastqFiles)) + ': ' \
		+ file + '\n'
		
		print(processStr)
	
		[barcodeRecognitionDict, perfectIndexReads, mismatchedButMatchableIndexReads, \
		unMatchableIndexReads] = \
		ParseIndexReadFile(file, mismatchRegexListDict, barcodeLookupTable, barcodeList)
	
		totalPerfectIndexReads += perfectIndexReads
		totalMismatchedButMatchableIndexReads += mismatchedButMatchableIndexReads
		totalUnMatchableIndexReads += unMatchableIndexReads
		
		# Output the barcode recognition table
		barcodeRecognitionOutputFileName = indexSummaryBaseDir + '/' \
		+ indexSummaryFilePrefix + '_' + str(j) + '.csv'
		
		
		OutputBarcodeRecognitionTable(barcodeRecognitionDict, barcodeRecognitionOutputFileName)
		
		indexAlignmentFiles.append(barcodeRecognitionOutputFileName)
	
		endTime = datetime.datetime.now()
	
		summaryLogStr = GenerateLogStringForIndividualFastQFileForIndexSummary(startTime, \
		endTime, j, len(indexFastqFiles), perfectIndexReads, mismatchedButMatchableIndexReads, \
		unMatchableIndexReads)
	
		print(summaryLogStr)
	
		# Write data to log file
		UpdateLogFileData(outputLog, processStr+summaryLogStr)
	
		j += 1


	completeEndTime = datetime.datetime.now()

	summaryStr = GenerateLogFileSummaryForIndexSummary(totalPerfectIndexReads, \
	totalMismatchedButMatchableIndexReads, totalUnMatchableIndexReads, initialStartTime, \
	completeEndTime)

	print(summaryStr)

	UpdateLogFileData(outputLog, summaryStr)


	return indexAlignmentFiles
# ------------------------------------------------------------------------------------------------ #




####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 3: Functions for Genome Alignment
# ------------------------------------------------------------------------------------------------ #
####################################################################################################


# ------------------------------------------------------------------------------------------------ #
def AlignReadsToGenome(mfaFiles, outputLog, \
genomeAlignmentBaseDir, genomeAlignmentFilePrefix, referenceIndexPrefix, \
bowtieMode='--end-to-end'):
# Code to align illumina reads to reference genome and return alignment point

	import gc
	import subprocess
	import datetime
	from pdb import set_trace
	import re
	
	totalLowScoreCount = 0
	totalNoAlignCount = 0
	totalMultipleAlignmentCount = 0
	totalPerfectAlignmentCount = 0
	totalStrangeFlagsSumCount = 0
	
	genomeAlignmentFiles = []

	xsRe = re.compile('XS:i:')

	logFileHandle = open(outputLog, 'a')
	initialStartTime = datetime.datetime.now()

	j = 1

	for file in mfaFiles:
	
		startTime = datetime.datetime.now()
	
		processStr = 'Processing ' + str(j) + ' of ' + str(len(mfaFiles)) + ': ' + file + '\n'
		print(processStr)

# 		fullFilePath = mfaBaseDir + '/' + file
		
# 		set_trace()
		
		alignmentData = \
		subprocess.check_output(['bowtie2', '-x', referenceIndexPrefix, '-f', file, \
		bowtieMode, '--threads', '4'])

		alignmentStr = alignmentData.decode("utf-8")
		alignmentLines = alignmentStr.split('\n')

		[genomeAlignmentDict, perfectAlignmentCount, multipleAlignmentCount, lowScoreCount, \
		noAlignCount, strangeFlagsSumCount] = ParseBowtie2Output(alignmentLines)
		
		# Update the read summary counts
		totalLowScoreCount += lowScoreCount
		totalNoAlignCount += noAlignCount
		totalMultipleAlignmentCount += multipleAlignmentCount
		totalPerfectAlignmentCount += perfectAlignmentCount
		totalStrangeFlagsSumCount += strangeFlagsSumCount

		# Write out the genome alignment files
		outputFileName = \
		genomeAlignmentBaseDir + '/' + genomeAlignmentFilePrefix + '_' + str(j) + '.csv'
		OutputGenomeAlignmentSummary(outputFileName, genomeAlignmentDict)
		
		genomeAlignmentFiles.append(outputFileName)
	
		# End timing
		endTime = datetime.datetime.now()

		summaryLogStr = GenerateLogStringForIndividualMfaFile(startTime, endTime, j, \
		len(mfaFiles), perfectAlignmentCount, multipleAlignmentCount, lowScoreCount, noAlignCount, \
		strangeFlagsSumCount)
	
		print(summaryLogStr)
	
		# Write data to log file
		UpdateLogFileData(outputLog, processStr+summaryLogStr)

		j += 1

	completeEndTime = datetime.datetime.now()

	summaryStr = GenerateLogFileSummary(totalPerfectAlignmentCount, totalMultipleAlignmentCount, \
	totalLowScoreCount, totalNoAlignCount, totalStrangeFlagsSumCount, initialStartTime, completeEndTime)

	print(summaryStr)

	UpdateLogFileData(outputLog, summaryStr)
	
	gc.collect()
		
	return genomeAlignmentFiles
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def StatisticsForCoordinateGroupsWithMoreThanOneMember(\
groupedPoolPresenceTable, rowPools, colPools, prPools, pcPools, controlPools, \
dTypeDictForPoolPresenceTable, \
threshold=5, loAddEntThresh=0,\
hiAddEntThresh_r=None, hiAddEntThresh_c=None, hiAddEntThresh_pr=None, hiAddEntThresh_pc=None):

# For the groups with more than one member, count how many can be assigned to library addresses
# Note to be include, address entries must be greater than loAddEntThresh, not equal to or greater
# than.

	groupsWithMoreThanOneEntry = 0
	groupsWithMoreThanOneEntryNotToIncludeInSummedTable = 0
	groupsWithMoreThanOneEntryToIncludeInSummedTable = 0
	groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress = 0

	if hiAddEntThresh_r == None:
		hiAddEntThresh_r = len(rowPools)
	if hiAddEntThresh_c == None:
		hiAddEntThresh_c = len(colPools)
	if hiAddEntThresh_pr == None:
		hiAddEntThresh_pr = len(prPools)
	if hiAddEntThresh_pc == None:
		hiAddEntThresh_pc = len(pcPools)

	i = 0

	while i < len(groupedPoolPresenceTable):
	
		if len(groupedPoolPresenceTable[i]) > 1:
			
			groupsWithMoreThanOneEntry += 1
			
			[summedPoolPresenceTableGroup, includeInSummedTable] = \
			SumPoolPresenceTableGroup(groupedPoolPresenceTable[i], \
			dTypeDictForPoolPresenceTable, rowPools, colPools, prPools, pcPools, controlPools)
		
			if includeInSummedTable == True:
			
				groupsWithMoreThanOneEntryToIncludeInSummedTable += 1
			
				addresses_r = FindAddressCoords(summedPoolPresenceTableGroup, rowPools, threshold=threshold)
				addresses_c = FindAddressCoords(summedPoolPresenceTableGroup, colPools, threshold=threshold)
				addresses_pr = FindAddressCoords(summedPoolPresenceTableGroup, prPools, threshold=threshold)
				addresses_pc = FindAddressCoords(summedPoolPresenceTableGroup, pcPools, threshold=threshold)
		
				nEntries_r = len(addresses_r)
				nEntries_c = len(addresses_c)
				nEntries_pr = len(addresses_pr)
				nEntries_pc = len(addresses_pc)
		
				if (loAddEntThresh < nEntries_r < hiAddEntThresh_r) \
				and (loAddEntThresh < nEntries_c < hiAddEntThresh_c) \
				and (loAddEntThresh < nEntries_pr < hiAddEntThresh_pr) \
				and (loAddEntThresh < nEntries_pc < hiAddEntThresh_pc):
					groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress += 1
			else:
				groupsWithMoreThanOneEntryNotToIncludeInSummedTable += 1
			
		
	
		i += 1
	
	return [groupsWithMoreThanOneEntry, groupsWithMoreThanOneEntryToIncludeInSummedTable, \
	groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress, \
	groupsWithMoreThanOneEntryNotToIncludeInSummedTable]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def StatisticsForCoordinateGroupsWithOnlyOneMember(\
groupedPoolPresenceTable, rowPools, colPools, prPools, pcPools, \
dTypeDictForPoolPresenceTable, \
threshold=5, loAddEntThresh=0,\
hiAddEntThresh_r=None, hiAddEntThresh_c=None, hiAddEntThresh_pr=None, hiAddEntThresh_pc=None):


	groupsWithOnlyOneEntry = 0
	groupsWithOnlyOneEntryWithValidLibraryAddress = 0
	
	i = 0

	if hiAddEntThresh_r == None:
		hiAddEntThresh_r = len(rowPools)
	if hiAddEntThresh_c == None:
		hiAddEntThresh_c = len(colPools)
	if hiAddEntThresh_pr == None:
		hiAddEntThresh_pr = len(prPools)
	if hiAddEntThresh_pc == None:
		hiAddEntThresh_pc = len(pcPools)


	while i < len(groupedPoolPresenceTable):
	
		if len(groupedPoolPresenceTable[i]) <= 1:
			groupsWithOnlyOneEntry += 1
		
			addresses_r = FindAddressCoords(groupedPoolPresenceTable[i][0], rowPools, threshold=threshold)
			addresses_c = FindAddressCoords(groupedPoolPresenceTable[i][0], colPools, threshold=threshold)
			addresses_pr = FindAddressCoords(groupedPoolPresenceTable[i][0], prPools, threshold=threshold)
			addresses_pc = FindAddressCoords(groupedPoolPresenceTable[i][0], pcPools, threshold=threshold)
		
			nEntries_r = len(addresses_r)
			nEntries_c = len(addresses_c)
			nEntries_pr = len(addresses_pr)
			nEntries_pc = len(addresses_pc)
		
			if (loAddEntThresh < nEntries_r < hiAddEntThresh_r) \
			and (loAddEntThresh < nEntries_c < hiAddEntThresh_c) \
			and (loAddEntThresh < nEntries_pr < hiAddEntThresh_pr) \
			and (loAddEntThresh < nEntries_pc < hiAddEntThresh_pc):
				groupsWithOnlyOneEntryWithValidLibraryAddress += 1
			
		
	
		i += 1
		
	return [groupsWithOnlyOneEntry, groupsWithOnlyOneEntryWithValidLibraryAddress]
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def AnalyzeGroupedPoolPresenceTable(groupedPoolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, dTypeDictForPoolPresenceTable, \
threshold=5, loAddEntThresh=0):

	[groupsWithOnlyOneEntry, groupsWithOnlyOneEntryWithValidLibraryAddress] = \
	StatisticsForCoordinateGroupsWithOnlyOneMember(groupedPoolPresenceTable, \
	rowPools, colPools, prPools, pcPools, dTypeDictForPoolPresenceTable, \
	threshold=threshold, loAddEntThresh=loAddEntThresh)
	
	
	[groupsWithMoreThanOneEntry, groupsWithMoreThanOneEntryToIncludeInSummedTable, \
	groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress, \
	groupsWithMoreThanOneEntryNotToIncludeInSummedTable] = \
	StatisticsForCoordinateGroupsWithMoreThanOneMember(groupedPoolPresenceTable, \
	rowPools, colPools, prPools, pcPools, controlPools, dTypeDictForPoolPresenceTable, \
	threshold=threshold, loAddEntThresh=loAddEntThresh)
	
	
	outputStr = "groupsWithOnlyOneEntry: " + str(groupsWithOnlyOneEntry) + '\n'
	
	outputStr += "groupsWithOnlyOneEntryWithValidLibraryAddress: " \
	+ str(groupsWithOnlyOneEntryWithValidLibraryAddress) + '\n'
	
	outputStr += "groupsWithMoreThanOneEntry" + str(groupsWithMoreThanOneEntry) + '\n'
	
	outputStr += "groupsWithMoreThanOneEntryNotToIncludeInSummedTable: " \
	+ str(groupsWithMoreThanOneEntryNotToIncludeInSummedTable) + '\n'
	
	outputStr += "groupsWithMoreThanOneEntryToIncludeInSummedTable: " \
	+ str(groupsWithMoreThanOneEntryToIncludeInSummedTable) + '\n'
	
	outputStr += "groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress: " \
	+ str(groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress)
	
	groupedPoolPresenceAnalysisDict = {}
	groupedPoolPresenceAnalysisDict['groupsWithOnlyOneEntry'] = groupsWithOnlyOneEntry
	groupedPoolPresenceAnalysisDict['groupsWithOnlyOneEntryWithValidLibraryAddress'] = groupsWithOnlyOneEntryWithValidLibraryAddress
	groupedPoolPresenceAnalysisDict['groupsWithMoreThanOneEntry: '] = groupsWithMoreThanOneEntry
	groupedPoolPresenceAnalysisDict['groupsWithMoreThanOneEntryNotToIncludeInSummedTable'] = groupsWithMoreThanOneEntryNotToIncludeInSummedTable
	groupedPoolPresenceAnalysisDict['groupsWithMoreThanOneEntryToIncludeInSummedTable'] = groupsWithMoreThanOneEntryToIncludeInSummedTable
	groupedPoolPresenceAnalysisDict['groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress'] = groupsWithMoreThanOneEntryToIncludeInSummedTableWithValidLibraryAddress
	
	
	
	return [outputStr, groupedPoolPresenceAnalysisDict]

# ------------------------------------------------------------------------------------------------ #






####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 4: Functions for Alignment Compilation
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def CompileAlignmentsWithIndexAndHimarSummaries(genomeAlignmentFiles, indexAlignmentFiles, \
himarRecognitionFiles, outputLog, genomeArrayFileName):

	import gc
	
	outputStr = gSeparatorString
	outputStr += 'Compiling Genome Alignments \n'
	UpdateLogFileData(outputLog, outputStr)
	
	# Count lines in the summary files
	[genomeAlignmentEntries, himarRecognitionEntries, indexAlignmentEntries, maxReadIDLength] \
	= CountEntriesInAllSummaryFiles(genomeAlignmentFiles, indexAlignmentFiles, himarRecognitionFiles)

	genomeHimarIndexLinesSummaryEntry = "Genome alignment entries: " + str(genomeAlignmentEntries) \
	+ '\n'
	genomeHimarIndexLinesSummaryEntry += "Index alignment entries: " + str(indexAlignmentEntries) \
	+ '\n'
	genomeHimarIndexLinesSummaryEntry += "Himar alignment entries: " + str(himarRecognitionEntries) \
	+ '\n'
	genomeHimarIndexLinesSummaryEntry += "Maximum read id length: " + str(maxReadIDLength) + '\n'
	print(genomeHimarIndexLinesSummaryEntry)
	
	
	UpdateLogFileData(outputLog, genomeHimarIndexLinesSummaryEntry)
	
	
	# Allocate arrays for import of data
	readIDFieldCode = 'a' + str(maxReadIDLength+2)

	genomeArray = ParseGenomeAlignmentFiles(readIDFieldCode, genomeAlignmentFiles, \
	genomeAlignmentEntries)
	indexArray = ParseIndexAlignmentFiles(readIDFieldCode, indexAlignmentFiles, \
	indexAlignmentEntries)
	himarArray =  ParseHimarAlignmentFiles(readIDFieldCode, himarRecognitionFiles, \
	himarRecognitionEntries)

	genomeArray = UpdateGenomeArrayWithHimarAndIndex(genomeArray, indexArray, himarArray)
	
	# Write compiled genome array
	WriteCompiledGenomeArray(genomeArrayFileName, genomeArray)
	
	
	# Generate the read taxonomy
	taxonomyDict = GenerateSudokuReadTaxonomy(genomeArray, outputLog)
	OutputReadTaxonomy(taxonomyDict, outputLog)
	
	UpdateLogFileData(outputLog, gSeparatorString)
	print(gSeparatorString)
	
	gc.collect()
	
	return genomeArray
# ------------------------------------------------------------------------------------------------ #


####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 5: Functions for Filling Pool Files
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def FillPoolFiles(barcodeFile, poolFileBaseDir, poolFilePrefix, outputLog, genomeArray=None, \
genomeArrayFileName=None, importGenomeFile=False):
	
	import gc
	
	outputStr = gSeparatorString
	outputStr += 'Filling Pool Files \n'
	outputStr += "Importing Barcodes\n"
	UpdateLogFileData(outputLog, outputStr)
	print(outputStr)
	
	[barcodeLookupTable, barcodes] = GenerateBarcodeLookupTable(barcodeFile)
	barcodeList = list(barcodes)
	
	# Either import the compiled genome array from a file, or deal with the one that was passed
	# as an argument
	if importGenomeFile == True:	
		outputStr = "Importing Genome Compilation File\n"
		print(outputStr)
		UpdateLogFileData(outputLog, outputStr)
	
		genomeArray = ImportGenomeCompilationFile(genomeCompilationFile)


	# Parse out the genome array into pool files
	outputStr = "Parsing Genome Array to Pools\n"
	print(outputStr)
	UpdateLogFileData(outputLog, outputStr)
	
	poolFiles = \
	ParseGenomeArrayToPoolFiles(genomeArray, barcodeLookupTable, poolFileBaseDir, poolFilePrefix)

	print(gSeparatorString)
	UpdateLogFileData(outputLog, gSeparatorString)
	
	gc.collect()
	
	if importGenomeFile == True:
		return poolFiles
	else:
		return [poolFiles, genomeArray]
# ------------------------------------------------------------------------------------------------ #


####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 6: Functions for initially populating the pool presence table
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def BuildInitialPoolPresenceTable(genomeArray, outputLog, barcodeFile, poolFiles, \
poolPresenceTableFileName):
	
	import numpy
	import pdb
	
	outputStr = gSeparatorString
	outputStr += 'Building Initial Pool Presence Table\n'
	UpdateLogFileData(outputLog, outputStr)
	print(outputStr)
	
	
	# Compile the valid reads
	outputStr = "Making Valid Genome Array\n"
	UpdateLogFileData(outputLog, outputStr)
	print(outputStr)
	
	validGenomeArray = GenerateValidCoordsList(genomeArray)
	validGenomeArray = numpy.sort(validGenomeArray, order='readAlignmentCoord')
	
	# Import the pool files
# 	outputStr = "Importing Pool Files\n"
# 	UpdateLogFileData(outputLog, outputStr)
# 	print(outputStr)	
# 	poolFilesDict = ImportPoolFiles(barcodeFile, poolFiles)

	# Generate the unique coordinates list
	outputStr = "Generating Unique Coordinates List\n"
	UpdateLogFileData(outputLog, outputStr)
	print(outputStr)
	
	uniqueCoords = GenerateUniqueCoordsList(validGenomeArray)

	# Make the first round of the pool presence table
	indexLookupTable = GeneratePoolNameToPoolCodeLookupTable(barcodeFile)

	dtypeDict = GenerateDTypeDictForPoolPresenceDict(indexLookupTable)
		

	print("Generating Pool Presence Table")
	poolPresenceTable = GeneratePoolPresenceTable(uniqueCoords, validGenomeArray, \
	indexLookupTable, dtypeDict)

	WritePoolPresenceTable(poolPresenceTableFileName, poolPresenceTable)

	return
# ------------------------------------------------------------------------------------------------ #


####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Step 7: Functions for analyzing the pool presence table
# ------------------------------------------------------------------------------------------------ #
####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def PlotEffectOfReadThresholdOnPoolPresenceTableTaxonomy(poolPresenceTable, \
rowPools, colPools, prPools, pcPools, controlPools, \
startThreshold, maxThreshold, poolPresenceTaxonomyTableFileName):

	import matplotlib.pyplot as plt
	from vectorOutput2 import generateOutputMatrixWithHeaders, writeOutputMatrix
	from pdb import set_trace

	# Count the number of pool presence table lines with only control wells present

	thresholds = []
	linesThatMapToLibraryAddresses = []
	linesThatMapToSingleLibraryAddresses = []
	linesThatMapToMultipleLibraryAddresses = []
	linesThatMapToUnambiguousLibraryAddresses = []
	linesThatMapToAmbiguousLibraryAddresses = []
	linesThatDoNotMapToLibraryAddresses = []
	linesThatHaveNoReadsAboveThresholdInAnyPool = []
	linesThatMapToControlIndexesOnly = []
	linesThatHaveCoordinatesInNoPoolAxis = []
	linesThatHaveCoordinatesInOnlyOnePoolAxis = []
	linesThatHaveCoordinatesInOnlyTwoPoolAxes = []
	linesThatHaveCoordinatesInOnlyThreePoolAxes = []

	threshold = startThreshold
	
	while threshold <= maxThreshold:
		
		poolPresenceTaxonomyDict = CalculatePoolPresenceTableTaxonomyDict(poolPresenceTable, \
		rowPools, colPools, prPools, pcPools, controlPools, threshold)
	
		thresholds.append(threshold)
		linesThatMapToLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToLibraryAddresses'])
		
		linesThatMapToSingleLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToSingleLibraryAddresses'])
		
		linesThatMapToMultipleLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToMultipleLibraryAddresses'])
		
		linesThatMapToUnambiguousLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToUnambiguousLibraryAddresses'])
		
		linesThatMapToAmbiguousLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatMapToAmbiguousLibraryAddresses'])
		
		linesThatDoNotMapToLibraryAddresses.append(\
		poolPresenceTaxonomyDict['linesThatDoNotMapToLibraryAddresses'])
		
		linesThatHaveNoReadsAboveThresholdInAnyPool.append(\
		poolPresenceTaxonomyDict['linesThatHaveNoReadsAboveThresholdInAnyPool'])
		
		linesThatMapToControlIndexesOnly.append(\
		poolPresenceTaxonomyDict['linesThatMapToControlIndexesOnly'])
		
		linesThatHaveCoordinatesInNoPoolAxis.append(\
		poolPresenceTaxonomyDict['linesThatHaveCoordinatesInNoPoolAxis'])
		
		linesThatHaveCoordinatesInOnlyOnePoolAxis.append(\
		poolPresenceTaxonomyDict['linesThatHaveCoordinatesInOnlyOnePoolAxis'])
		
		linesThatHaveCoordinatesInOnlyTwoPoolAxes.append(\
		poolPresenceTaxonomyDict['linesThatHaveCoordinatesInOnlyTwoPoolAxes'])
		
		linesThatHaveCoordinatesInOnlyThreePoolAxes.append(\
		poolPresenceTaxonomyDict['linesThatHaveCoordinatesInOnlyThreePoolAxes'])
	
	
		threshold += 1


	headers = ['threshold', 
	'linesThatMapToLibraryAddresses', \
	'linesThatMapToSingleLibraryAddresses', \
	'linesThatMapToMultipleLibraryAddresses', \
	'linesThatMapToUnambiguousLibraryAddresses', \
	'linesThatMapToAmbiguousLibraryAddresses', \
	'linesThatDoNotMapToLibraryAddresses', \
	'linesThatHaveNoReadsAboveThresholdInAnyPool', \
	'linesThatMapToControlIndexesOnly', \
	'linesThatHaveCoordinatesInNoPoolAxis', \
	'linesThatHaveCoordinatesInOnlyOnePoolAxis', \
	'linesThatHaveCoordinatesInOnlyTwoPoolAxes', \
	'linesThatHaveCoordinatesInOnlyThreePoolAxes']

	vectorList = [\
	thresholds, \
	linesThatMapToLibraryAddresses, \
	linesThatMapToSingleLibraryAddresses, \
	linesThatMapToMultipleLibraryAddresses, \
	linesThatMapToUnambiguousLibraryAddresses, \
	linesThatMapToAmbiguousLibraryAddresses, \
	linesThatDoNotMapToLibraryAddresses, \
	linesThatHaveNoReadsAboveThresholdInAnyPool, \
	linesThatMapToControlIndexesOnly, \
	linesThatHaveCoordinatesInNoPoolAxis, \
	linesThatHaveCoordinatesInOnlyOnePoolAxis, \
	linesThatHaveCoordinatesInOnlyTwoPoolAxes, \
	linesThatHaveCoordinatesInOnlyThreePoolAxes]


	oMatrix = generateOutputMatrixWithHeaders(vectorList, headers)
	writeOutputMatrix(poolPresenceTaxonomyTableFileName, oMatrix)
	
	plt.figure()
	plt.plot(thresholds, linesThatMapToLibraryAddresses, \
	label='linesThatMapToLibraryAddresses')
	
	plt.plot(thresholds, linesThatMapToSingleLibraryAddresses, \
	label='linesThatMapToSingleLibraryAddresses')
	
	plt.plot(thresholds, linesThatMapToMultipleLibraryAddresses, \
	label='linesThatMapToMultipleLibraryAddresses')
	
	plt.plot(thresholds, linesThatMapToUnambiguousLibraryAddresses, \
	label='linesThatMapToUnambiguousLibraryAddresses')
	
	plt.plot(thresholds, linesThatMapToAmbiguousLibraryAddresses, \
	label='linesThatMapToAmbiguousLibraryAddresses')
	
	plt.plot(thresholds, linesThatDoNotMapToLibraryAddresses, \
	label='linesThatDoNotMapToLibraryAddresses')
	
	plt.plot(thresholds, linesThatHaveNoReadsAboveThresholdInAnyPool, \
	label='linesThatHaveNoReadsAboveThresholdInAnyPool')
	
	plt.plot(thresholds, linesThatMapToControlIndexesOnly, \
	label='linesThatMapToControlIndexesOnly')
	
	plt.plot(thresholds, linesThatHaveCoordinatesInNoPoolAxis, \
	label='linesThatHaveCoordinatesInNoPoolAxis')
	
	plt.plot(thresholds, linesThatHaveCoordinatesInOnlyOnePoolAxis, \
	label='linesThatHaveCoordinatesInOnlyOnePoolAxis')
	
	plt.plot(thresholds, linesThatHaveCoordinatesInOnlyTwoPoolAxes, \
	label='linesThatHaveCoordinatesInOnlyTwoPoolAxes')
	
	plt.plot(thresholds, linesThatHaveCoordinatesInOnlyThreePoolAxes, \
	label='linesThatHaveCoordinatesInOnlyThreePoolAxes')
	
	plt.grid()
	plt.legend()
	plt.xlabel('Read Count Threshold')
	plt.show()
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateRatioOfReadNumbersForSingleAddressLines(poolPresenceTable, rowPools, colPools, \
prPools, pcPools, controlPools, threshold, logReadNumberRatioHistogramsFileName, \
fitFileName):
	
	from numpy import array, arange
	from vectorOutput2 import generateOutputMatrixWithHeaders, writeOutputMatrix
	
	import pdb
	
	[nr2ncArray, nr2nprArray, nr2npcArray, nc2nprArray, nc2npcArray, npr2npcArray, \
	logNr2ncArray, logNr2nprArray, logNr2npcArray, logNc2nprArray, logNc2npcArray, \
	logNpr2npcArray] = \
	CalculateReadNumberRatiosForSingleAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
	pcPools, controlPools, threshold)

	a0 = 4000
	c0 = -0.5
	delta0 = 0.35
	sigma0 = 0.5

	plsq0 = array([a0, c0, delta0, sigma0])

	logBins = arange(-10,10.2,0.2)
	
	
	[valueslogNr2nc, baselogNr2nc, integralNr2nc, plsqNr2nc, voigtlogNr2nc, binCentersLogNr2nc] = \
	FitReadNumberRatiosAndPlot(logNr2ncArray, logBins, "log(n_row/n_col)", plsq0, integrate=True)

	[valueslogNr2npr, baselogNr2npr, integralNr2npr, plsqNr2npr, voigtlogNr2npr, \
	binCentersLogNr2npr] = \
	FitReadNumberRatiosAndPlot(logNr2nprArray, logBins, "log(n_row/n_pr)", plsq0, integrate=True)

	[valueslogNr2npc, baselogNr2npc, integralNr2npc, plsqNr2npc, voigtlogNr2npc, \
	binCentersLogNr2npc] = \
	FitReadNumberRatiosAndPlot(logNr2npcArray, logBins, "log(n_row/n_pc)", plsq0, integrate=True)
 
	[valueslogNc2npr, baselogNc2npr, integralNc2npr, plsqNc2npr, voigtlogNc2npr, \
	binCentersLogNc2npr] = \
	FitReadNumberRatiosAndPlot(logNc2nprArray, logBins, "log(n_col/n_pr)", plsq0, integrate=True)
 
	[valueslogNc2npc, baselogNc2npc, integralNc2npc, plsqNc2npc, voigtlogNc2npc, \
	binCentersLogNc2npc] = \
	FitReadNumberRatiosAndPlot(logNc2npcArray, logBins, "log(n_col/n_pc)", plsq0, integrate=True)

	[valueslogNpr2npc, baselogNpr2npc, integralNpr2npc, plsqNpr2npc, voigtlogNpr2npc, \
	binCentersLogNpr2npc] = \
	FitReadNumberRatiosAndPlot(logNpr2npcArray, logBins, "log(n_pr/n_pc)", plsq0, integrate=True)

	
	headers = [\
	'baselogNr2nc', 'valueslogNr2nc', \
	'baselogNr2npr', 'valueslogNr2npr', \
	'baselogNr2npc', 'valueslogNr2npc', \
	'baselogNc2npr', 'valueslogNc2npr', \
	'baselogNc2npc', 'valueslogNc2npc', \
	'baselogNpr2npc', 'valueslogNpr2npc' ]
	
	fitHeaders = [\
	'binCenterNr2nc', 'fitValueslogNr2nc', \
	'binCenterNr2npr', 'fitValueslogNr2npr', \
	'binCenterNr2npc', 'fitValueslogNr2npc', \
	'binCenterNc2npr', 'fitValueslogNc2npr', \
	'binCenterNc2npc', 'fitValueslogNc2npc', \
	'binCenterNpr2npc', 'fitValueslogNpr2npc' ]


	vectorList = [\
	binCentersLogNr2nc, valueslogNr2nc, \
	binCentersLogNr2npr, valueslogNr2npr, \
	binCentersLogNr2npc, valueslogNr2npc, \
	binCentersLogNc2npr, valueslogNc2npr, \
	binCentersLogNc2npc, valueslogNc2npc, \
	binCentersLogNpr2npc, valueslogNpr2npc]
	
	fitVectorList = [\
	binCentersLogNr2nc, voigtlogNr2nc, \
	binCentersLogNr2npr, voigtlogNr2npr, \
	binCentersLogNr2npc, voigtlogNr2npc, \
	binCentersLogNc2npr, voigtlogNc2npr, \
	binCentersLogNc2npc, voigtlogNc2npc, \
	binCentersLogNpr2npc, voigtlogNpr2npc ]
	
# 	pdb.set_trace()

	oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
	writeOutputMatrix(logReadNumberRatioHistogramsFileName, oMatrix)
	
	oMatrixFit = generateOutputMatrixWithHeaders(fitVectorList, fitHeaders, delimeter=',')
	writeOutputMatrix(fitFileName, oMatrixFit)
	
	
	logReadNumberRatioHistogramFitDict = {\
	'nr2nc':plsqNr2nc, 'nr2npr':plsqNr2npr, \
	'nr2npc':plsqNr2npc, 'nc2npr':plsqNc2npr, \
	'nc2npc':plsqNc2npc, 'npr2npc':plsqNpr2npc }
	
	logReadNumberRatioHistogramIntegralDict = {\
	'nr2nc':integralNr2nc, 'nr2npr':integralNr2npr, \
	'nr2npc':integralNr2npc, 'nc2npr':integralNc2npr, \
	'nc2npc':integralNc2npc, 'npr2npc':integralNpr2npc }
	
	return [logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateScoreHistogramForSingleAddressLines(poolPresenceTable, rowPools, colPools, prPools, \
pcPools, controlPools, logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict,\
threshold, cumulativeScoreDistributionFileName):
	
	from numpy import histogram, log, arange, cumsum
	from matplotlib.pyplot import hist
	import matplotlib.pyplot as plt
	
	i = 0
	scoreArray = []
	
	while i < len(poolPresenceTable):
		
		poolPresenceTableLine = poolPresenceTable[i]
		
		possibleAddressesAndScores = \
		CalculateLibraryAddressesForPoolPresenceTableLine2(poolPresenceTableLine, rowPools, \
		colPools, prPools, pcPools, logReadNumberRatioHistogramFitDict, \
		logReadNumberRatioHistogramIntegralDict, threshold)
		
		# Calculate if the line is single occupancy
		lenPossibleAddresses = len(possibleAddressesAndScores)
	
		if lenPossibleAddresses == 1:	
			scoreArray.append(possibleAddressesAndScores[0][2])

		i += 1
		
	logScoreArray = log(scoreArray)
	
	logBins = arange(min(logScoreArray),max(scoreArray),0.1)
	bins = arange(0,max(scoreArray),0.001)
	
	values, base = histogram(logScoreArray, bins=logBins)
	cumulative = cumsum(values)
	
	plt.figure()
	hist(scoreArray, bins)
	plt.xlabel("Score")
	plt.ylabel("Number of Single Address Lines")
	
	plt.figure()
	plt.plot(base[:-1],cumulative/max(cumulative))
	plt.title("Cumulative Score Distribution of Single Address Pool Presence Table Entries")
	plt.xlabel("Natural Log of Score")
	plt.ylabel("Fraction of Single Address Pool Presence Table Entries") 
	plt.show()
	

	return [scoreArray, cumulative, base[:-1]]
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ImportCherryPickedLibrary(fileName):
# Import a file listing picked wells from library
	
	import pdb
	import re
	
	wellIDRe = re.compile(r'([A-H])(\d+)', re.IGNORECASE)
	noGrowthReportedRe = re.compile(r'No Growth', re.IGNORECASE)
	predictedGenomicCoordRe = re.compile('\d+')
	
	fileHandle = open(fileName, 'r')
	data = fileHandle.readlines()
	header = data[0].strip().split(',')
	
	cherryPickedList = []
	
	i = 1
	while i < len(data):
		dataLine = data[i].strip().split(',')
		
		featureName = dataLine[0]
		
		predictedGenomicCoordMatch = predictedGenomicCoordRe.match(dataLine[1])
		
		if predictedGenomicCoordMatch != None:
			predictedGenomicCoord = int(predictedGenomicCoordMatch.group())
		else:
			predictedGenomicCoord = None
		
		
		try:
			sourcePlate = int(dataLine[2])
		except:
			pdb.set_trace()
		
		sourceWell = dataLine[3]
		condensedPlate = dataLine[4]		
		condensedWell = dataLine[5]
		comment = dataLine[6]
		
		noGrowthReportedMatch = noGrowthReportedRe.search(comment)
		if noGrowthReportedMatch != None:
			growthReported = False
		else:
			growthReported = True
		
		
		sourceWellMatch = wellIDRe.search(sourceWell)
		
		if sourceWellMatch != None:
			sourceRow = sourceWellMatch.group(1)
			sourceCol = sourceWellMatch.group(2)
		else:
			pdb.set_trace()

		
		condensedWellMatch = wellIDRe.search(condensedWell)
		condensedRow = condensedWellMatch.group(1)
		condensedCol = condensedWellMatch.group(2)
		
		
		sudokuCherryPick = SudokuCherryPick(featureName, predictedGenomicCoord, sourcePlate, \
		sourceRow, sourceCol, condensedPlate, condensedRow, condensedCol, comment, growthReported)
		
		
		cherryPickedList.append(sudokuCherryPick)
		
		i += 1
	

	return cherryPickedList
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindClosestGenomicCoord(readAlignmentCoords, reportedGenomicCoord):
	
	import pdb
	
	if len(readAlignmentCoords) == 1:
		latestCalculatedGenomicCoord = readAlignmentCoords[0].coord
	elif len(readAlignmentCoords) == 0:
		latestCalculatedGenomicCoord = 0
	elif len(readAlignmentCoords) > 1:
		j = 0
		diffList = []
		coordList = []
		
		while j < len(readAlignmentCoords):
			
			if reportedGenomicCoord == None:
				reportedGenomicCoord = 0
			
			diffList.append(abs(readAlignmentCoords[j].coord - reportedGenomicCoord))
			coordList.append(readAlignmentCoords[j].coord)
			j += 1
		sortedDiffs = sorted(diffList)
		if len(sortedDiffs) == 0:
			pdb.set_trace()
		latestCalculatedGenomicCoord = coordList[diffList.index(sortedDiffs[0])]
		

	return latestCalculatedGenomicCoord

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CalculateDifferenceBetweenReportedAndLatestCalculatedCoords(multipleOccupancy, \
readAlignmentCoords, reportedGenomicCoord):
	
	import pdb
	
	if reportedGenomicCoord == None:
		reportedGenomicCoord = 0
	
	if multipleOccupancy == False:
# 		pdb.set_trace()
		diff = readAlignmentCoords[0].coord - reportedGenomicCoord
		
	elif multipleOccupancy == True:
		closestGenomicCoord = FindClosestGenomicCoord(readAlignmentCoords, reportedGenomicCoord)
		diff = closestGenomicCoord - reportedGenomicCoord
	elif multipleOccupancy == None:
		diff = reportedGenomicCoord
		
	return diff
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def CheckSudokuCherryPickAgainstSudokuWell(plateRow, plateCol, cherryPick, sudokuGridLookupDict):

	cherryPick.sourcePlateRow = plateRow
	cherryPick.sourcePlateCol = plateCol

	sourceRow = cherryPick.sourceRow
	sourceCol = cherryPick.sourceCol
	
	sourceWell = sourceRow + sourceCol

	sudokuWell = \
	sudokuGridLookupDict[plateRow][plateCol].wellGrid[sourceRow][sourceCol]

	cherryPickSourcePlate = sudokuWell.plateName
	
	if len(sudokuWell.readAlignmentCoords) == 1:
		print("Single occupancy for " + str(cherryPickSourcePlate) + " well " + sourceWell)
		cherryPick.multipleOccupancy = False
		returnCode = 'singleOccupancy'
	
	elif len(sudokuWell.readAlignmentCoords) == 0:
		print("No calculated contents for " + str(cherryPickSourcePlate) + " well " + sourceWell)
		cherryPick.latestCalculatedGenomicCoord = 'Empty'
		returnCode = 'noOccupancy'
	
	elif len(sudokuWell.readAlignmentCoords) > 1:
		print("Multiple calculated contents for " + str(cherryPickSourcePlate) + " well " \
		+ sourceWell)
		cherryPick.multipleOccupancy = True
		returnCode = 'multipleOccupancy'
		
	
	cherryPick.latestCalculatedGenomicCoord = FindClosestGenomicCoord(\
	sudokuWell.readAlignmentCoords, cherryPick.reportedGenomicCoord)
	
	
	cherryPick.diffLatestReportedGenomicCoord = \
	CalculateDifferenceBetweenReportedAndLatestCalculatedCoords(cherryPick.multipleOccupancy, \
	sudokuWell.readAlignmentCoords, cherryPick.reportedGenomicCoord)

	return returnCode
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def AssignLatestCoordsOccupancySourcePlateRowsAndColumnsToCherryPickedList(cherryPickedList, \
sudokuGridLookupDict, prPools, pcPools):

	import re

	multipleOccupancies = 0
	singleOccupancies = 0
	noOccupancies = 0

	
	plateNameRe = re.compile(r'\d+')
	
	i = 0
	while i < len(cherryPickedList):
		
		cherryPick = cherryPickedList[i]
		cherryPickSourcePlate = cherryPickedList[i].sourcePlate
		
		matchingPlateFound = False
		
		for plateRow in prPools:
			for plateCol in pcPools:
				plate =  sudokuGridLookupDict[plateRow][plateCol].plateName
				
				plateNameMatch = plateNameRe.match(plate)
				
				if plateNameMatch != None:
					
					if int(plate) == cherryPickSourcePlate:
						
						matchingPlateFound = True
						
						returnCode = CheckSudokuCherryPickAgainstSudokuWell(plateRow, plateCol, \
						cherryPick, sudokuGridLookupDict)
						
						if returnCode == 'multipleOccupancy':
							multipleOccupancies += 1
						elif returnCode == 'singleOccupancy':
							singleOccupancies += 1
						elif returnCode == 'noOccupancy':
							noOccupancies += 1
						
		if matchingPlateFound == False:
			print("No matching plate found in sudoku look up grid for plate " \
			+ str(cherryPickSourcePlate))
		i += 1
	
	print("No occupancies: " + str(noOccupancies))
	print("Single occupancies: " + str(singleOccupancies))
	print("Multiple occupancies: " + str(multipleOccupancies))

	
	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def WriteSudokuCherryPickList(filename, cherryPickedList):
	
	fileHandle = open(filename, 'w')
	
	headerStr = ''
	headerStr += 'Reported Feature Name,'
	headerStr += 'ReportedGenomicCoord,'
	headerStr += 'Source Plate,'
	headerStr += 'Source Row,'
	headerStr += 'Source Col,'
	headerStr += 'Condensed Plate,'
	headerStr += 'Condensed Row,'
	headerStr += 'Condensed Col,'
	headerStr += 'Comment,'
	headerStr += 'Multiple Occupancy,'
	headerStr += 'Growth Reported,'
	headerStr += 'Source Plate Row,'
	headerStr += 'Source Plate Col,'
	headerStr += 'Latest Calculated Genomic Coord,'
	headerStr += 'Latest - Reported Coord,'
	headerStr += 'Latest Feature Name\n'
	
	fileHandle.write(headerStr)
	
	i = 0
	while i < len(cherryPickedList):
		
		pick = cherryPickedList[i]
		
		outputStr = ''
		outputStr += pick.reportedFeatureName + ','
		outputStr += str(pick.reportedGenomicCoord)  + ','
		outputStr += str(pick.sourcePlate)  + ','
		outputStr += pick.sourceRow  + ','
		outputStr += str(pick.sourceCol)  + ','
		outputStr += pick.condensedPlate  + ','
		outputStr += pick.condensedRow  + ','
		outputStr += pick.condensedCol  + ','
		outputStr += '"' + pick.comment + '"'  + ','
		outputStr += str(pick.multipleOccupancy)  + ','
		outputStr += str(pick.growthReported)  + ','
		outputStr += pick.sourcePlateRow + ','
		outputStr += pick.sourcePlateCol + ','
		outputStr += str(pick.latestCalculatedGenomicCoord)  + ','
		outputStr += str(pick.diffLatestReportedGenomicCoord) + ','
		outputStr += str(pick.latestFeatureName) + '\n'
		
		fileHandle.write(outputStr)
		
		i += 1

	fileHandle.close()

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindCherryPickWellsWithDisruptionBetweenStartAndEndCoords(startCoord, endCoord, featureName, \
cherryPickList):
	
	i = 0
	
	wellsToReturn = []
	
	while i < len(cherryPickList):
							
		cherryPick = cherryPickList[i]		
		coord = cherryPick.latestCalculatedGenomicCoord
		
		if startCoord < coord < endCoord:
			wellsToReturn.append(cherryPick)
			
			cherryPick.latestFeatureName = featureName
				
		i += 1
		
	return wellsToReturn
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindFeatureName(feature):
	
	featureTagDictKeys = list(feature.tagDict.keys())
	
	if 'gene' in featureTagDictKeys:
		geneName = feature.tagDict['gene'][0]
	elif 'note' in featureTagDictKeys:
		geneName = feature.tagDict['note'][0]
	elif 'locus_tag' in featureTagDictKeys:
		geneName = feature.tagDict['locus_tag'][0]
	else:
		geneName = ''
	
	return geneName
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def UpdateFeatureArrayWithCondensedLibraryEntries(featureArray, cherryPickList, \
rowPools, colPools, prPools, pcPools):

	import pdb

	i = 0 
	while i < len(featureArray):
	
		feature = featureArray[i]		
		startCoord = feature.startCoord
		endCoord = feature.endCoord
		
		featureName = FindFeatureName(feature)
		
		wellsWithDisruptionsMatchingFeature = \
		FindCherryPickWellsWithDisruptionBetweenStartAndEndCoords(startCoord, endCoord, \
		featureName, cherryPickList)
		
		for well in wellsWithDisruptionsMatchingFeature:
			feature.cherryPicks.append(well)

		
		featureArray[i].assignDisruptionsClosestToTranslationStart()
		featureArray[i].findIfAnySingleOccupancyWellHasBeenPicked()

		
		i += 1

	return featureArray
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def TestIfReportedFeatureNameEqualsLatestCalculatedFeatureName(cherryPickList):
	
	matches = 0
	noMatches = 0
	
	i = 0
	while i < len(cherryPickList):
		cherryPick = cherryPickList[i]
		
		if cherryPick.reportedFeatureName == cherryPick.latestFeatureName:
			matches += 1
		else:
			noMatches += 1
		
		i += 1
	
	print("Cherry Picks With Feature Name That Matches Old Feature Name: " + str(matches))
	print("Cherry Picks With Feature Name That Does Not Match Old Feature Name: " + str(noMatches))
	
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class SudokuCherryPick:
	
	def __init__(self, reportedFeatureName, reportedGenomicCoord, sourcePlate, sourceRow, \
	sourceCol, condensedPlate, condensedRow, condensedCol, comment, growthReported):
		
		self.reportedFeatureName = reportedFeatureName
		self.reportedGenomicCoord = reportedGenomicCoord
		self.sourcePlate = sourcePlate
		self.sourceRow = sourceRow
		self.sourceCol = sourceCol
		self.condensedPlate = condensedPlate
		self.condensedRow = condensedRow
		self.condensedCol = condensedCol
		self.comment = comment
		self.multipleOccupancy = None
		self.growthReported = growthReported
		self.sourcePlateRow = None
		self.sourcePlateCol = None
		self.latestCalculatedGenomicCoord = None
		self.diffLatestReportedGenomicCoord = 0
		self.latestFeatureName = ''
# ------------------------------------------------------------------------------------------------ #



####################################################################################################
# ------------------------------------------------------------------------------------------------ #
# Sanger verification code

# ------------------------------------------------------------------------------------------------ #
class BlastAlignment:

	def __init__(self, score, expectation):
		
		import re
		
		self.score = score
		self.expectation = expectation
		self.strandOrientation = ''
		self.alignmentLines = []
		self.readAlignmentCoord = 0
		self.strandOrientationRe = re.compile(r"(Plus|Minus)/(Minus|Plus)")
		self.alignmentSubjectLineRe = re.compile(r'Sbjct\s+(\d+)\s+\w+\s+(\d+)')
	
	def CalculateReadAlignmentCoord(self):
		
		import pdb
		
# 		print(self.strandOrientation)
		strandOrientationMatch = self.strandOrientationRe.search(self.strandOrientation)
		
		strandOrientationQuery = strandOrientationMatch.group(1)
		strandOrientationSubject = strandOrientationMatch.group(2)
		
		if strandOrientationQuery == 'Plus':
			if strandOrientationSubject == 'Minus':
# 				pdb.set_trace()
				readAlignmentCoord = self.alignmentSubjectLineRe.search(self.alignmentLines[0]).group(1)
			elif strandOrientationSubject == 'Plus':
				readAlignmentCoord = self.alignmentSubjectLineRe.search(self.alignmentLines[0]).group(1)
		elif strandOrientationQuery == 'Minus':
			print("Error!")
			readAlignmentCoord = 0
			
		self.readAlignmentCoord = readAlignmentCoord
		
		return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ParseBlastOutput3(alignmentLines):
	
	import re
	import pdb
	
	strandOrientationRe = re.compile(r"Strand=((Plus|Minus)/(Minus|Plus))")
	alignmentSubjectLineRe = re.compile(r'Sbjct\s+(\d+)\s+\w+\s+(\d+)')
	noHitsRe = re.compile(r'\*\*\*\*\* No hits found \*\*\*\*\*')
	newAlignmentRe = re.compile(r'Score\s+=\s+(\d+)\s+bits\s+\(\d+\),\s+Expect\s+=\s+(\d+e*\.*[+|-]*\d+)')
	
	alignments = []
	noHitsFound = False
	currentAlignment = None
	
	for line in alignmentLines:
		if noHitsRe.search(line) != None:
			noHitsFound = True
			readAlignmentCoords = None
	
	
	if noHitsFound == False:
		for line in alignmentLines:
		
			newAlignmentMatch = newAlignmentRe.search(line)
			strandOrientationMatch = strandOrientationRe.search(line)
			alignmentSubjectLineMatch = alignmentSubjectLineRe.search(line)
		
			if newAlignmentMatch != None:
			
				if currentAlignment != None:
# 					print(currentAlignment.score)
					currentAlignment.CalculateReadAlignmentCoord()
					alignments.append(currentAlignment)
			
				score = newAlignmentMatch.group(1)
				expectation = newAlignmentMatch.group(2)
				currentAlignment = BlastAlignment(score, expectation)
		
			elif (newAlignmentMatch == None) and (currentAlignment != None):
				if strandOrientationMatch != None:
					currentAlignment.strandOrientation = strandOrientationMatch.group(1)
				if alignmentSubjectLineMatch != None:
					currentAlignment.alignmentLines.append(line)
			
	if noHitsFound == False:
# 		pdb.set_trace()
# 		print(currentAlignment.score)	
		currentAlignment.CalculateReadAlignmentCoord()
		alignments.append(currentAlignment)	
	

	return alignments
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AlignSequenceToReferenceByBlast(referenceFile, queryFilePath):

	import subprocess
	
	alignmentData = \
	subprocess.check_output(['blastn', '-subject', referenceFile, '-query', queryFilePath])

	alignmentStr = alignmentData.decode("utf-8")
	alignmentLines = alignmentStr.split('\n')
	
	# Find read alignment coordinates
	readAlignmentCoords = ParseBlastOutput3(alignmentLines)

	return readAlignmentCoords

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportVerificationDataFile(sangerVerificationFile):
	
	fileHandle = open(sangerVerificationFile, 'r')
	
	data = fileHandle.readlines()
	
	verificationDict = {}
	
	i = 1
	while i < len(data):
		
		dataLine = data[i].strip().split(',')
		
		sequenceName = dataLine[0]
		fileDir = dataLine[1]
		fileName = dataLine[2]
		well = dataLine[3]
		plate = dataLine[4]
		predictionRound = dataLine[5]
		
		
		verificationDict[sequenceName] = \
		{'fileDir':fileDir, 'fileName':fileName, 'well':well, 'plate':plate, \
		'predictionRound':predictionRound}
		
		i += 1
	
	return verificationDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GenerateSudokuReadAlignmentOutputStrings(sudokuReadAlignmentCoords):
	
	import pdb
	
	sudokuReadAlignmentCoordsArray = []
	
	j = 0
	while j < len(sudokuReadAlignmentCoords):
	
# 		pdb.set_trace()
		
		coord = sudokuReadAlignmentCoords[j].coord
		sudokuReadAlignmentCoordsArray.append(coord)
		j += 1
	
	sudokuReadAlignmentCoordsStr = '"'
	
	j = 0
	
	while j < len(sudokuReadAlignmentCoordsArray):
		
		sudokuReadAlignmentCoordsStr += str(sudokuReadAlignmentCoordsArray[j])
		
		if j < len(sudokuReadAlignmentCoordsArray) - 1:
			sudokuReadAlignmentCoordsStr += ","
	
		j += 1
	
	
	sudokuReadAlignmentCoordsStr += '"'

	return sudokuReadAlignmentCoordsStr

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateSangerReadAlignmentOutputStrings(sangerReadAlignmentCoords):
	
	sangerReadAlignmentCoordsArray = []
	sangerScoresArray = []
	sangerExpectationsArray = []
	sangerOrientationsArray = []
	
	j = 0
	while j < len(sangerReadAlignmentCoords):
		
		coord = sangerReadAlignmentCoords[j].readAlignmentCoord
		expectation = sangerReadAlignmentCoords[j].expectation
		score = sangerReadAlignmentCoords[j].score
		orientation = sangerReadAlignmentCoords[j].strandOrientation
		
		sangerReadAlignmentCoordsArray.append(coord)
		sangerScoresArray.append(expectation)
		sangerExpectationsArray.append(score)
		sangerOrientationsArray.append(orientation)
		
		j += 1
	
	sangerReadAlignmentCoordsStr = '"'
	sangerScoresStr = '"'
	sangerExpectationsStr = '"'
	sangerOrientationsStr = '"'
	
	j = 0
	
	while j < len(sangerReadAlignmentCoordsArray):
		
		sangerReadAlignmentCoordsStr += str(sangerReadAlignmentCoordsArray[j])
		sangerScoresStr += str(sangerScoresArray[j])
		sangerExpectationsStr += str(sangerExpectationsArray[j])
		sangerOrientationsStr += str(sangerOrientationsArray[j])
		
		if j < len(sangerReadAlignmentCoordsArray) - 1:
			sangerReadAlignmentCoordsStr += ","
			sangerScoresStr += ","
			sangerExpectationsStr += ","
			sangerOrientationsStr += ","
	
		j += 1
	
	
	sangerReadAlignmentCoordsStr += '"'
	sangerScoresStr += '"'
	sangerExpectationsStr += '"'
	sangerOrientationsStr += '"'

	return [sangerReadAlignmentCoordsStr, sangerScoresStr, sangerExpectationsStr, \
	sangerOrientationsStr]

# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def OutputVerificationDict(verificationDict, outputLog):
	
	import pdb
	
	fileHandle = open(outputLog, 'w')

	verificationKeys = sorted(list(verificationDict.keys()))
	
	headerStr = "Nominal Genetic Locus,File Directory,File Name,Well,Plate,Prediction Round," \
	+ "Sanger Read Alignment Coordinates,Sanger Read Alignment Score," \
	+ "Sanger Read Alignment Expectation,Sanger Read Alignment Orientation," \
	+ "Sudoku Read Alignment Coordinates\n"
	
	fileHandle.write(headerStr)
	
	outputStr = ""
	
	i = 0
	while i < len(verificationKeys):
		key = verificationKeys[i]
		
		verificationDataDict = verificationDict[key]
		
		fileDir = verificationDataDict['fileDir']
		fileName = verificationDataDict['fileName']
		well = verificationDataDict['well']
		plate = verificationDataDict['plate']
		predictionRound = verificationDataDict['predictionRound']
		sangerReadAlignmentCoords = verificationDataDict['sangerReadAlignmentCoords']
		sudokuReadAlignmentCoords = verificationDataDict['sudokuReadAlignmentCoords']
		
		[sangerReadAlignmentCoordsStr, sangerScoresStr, sangerExpectationsStr, \
		sangerOrientationsStr] = GenerateSangerReadAlignmentOutputStrings(sangerReadAlignmentCoords)
		
		sudokuReadAlignmentOutputStr = \
		GenerateSudokuReadAlignmentOutputStrings(sudokuReadAlignmentCoords)
		
		
		outputStr += key + ',' + fileDir + ',' + fileName + ',' + well + ',' + plate \
		+ ',' + predictionRound + ',' + sangerReadAlignmentCoordsStr + ',' + sangerScoresStr \
		+ ',' + sangerExpectationsStr + ',' + sangerOrientationsStr + ',' \
		+ sudokuReadAlignmentOutputStr + '\n'
		
		
		i += 1
	
	fileHandle.write(outputStr)

	fileHandle.close()

# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def ConvertPlateNumberToPlateRowAndPlateColCoords(plateNumber, sudokuGridLookupDict, prPools, \
pcPools):
	import re	

	plateNameRe = re.compile(r'\d+')
	
	matchingPlateFound = False
	matchingPlateCol = None
	matchingPlateRow = None
	
	i = 0
	while (i < len(prPools)) and (matchingPlateFound == False):
		plateRow = prPools[i]
		j = 0
		while (j < len(pcPools)) and (matchingPlateFound == False):
			plateCol = pcPools[j]
			plate =  sudokuGridLookupDict[plateRow][plateCol].plateName
			plateNameMatch = plateNameRe.match(plate)

			if plateNameMatch != None:
				if int(plate) == int(plateNumber):
					matchingPlateFound = True
					matchingPlateCol = plateCol
					matchingPlateRow = plateRow
			j += 1
		i += 1
	
	return [matchingPlateRow, matchingPlateCol]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
class GenomeSequenceIndexFailure(Exception):
	pass
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindClosestRealDisruptionIndex(disruptionLocIndex, sequence, chromosomeStartIndex, \
chromosomeEndIndex, megaplasmidStartIndex, megaplasmidEndIndex):
	
	import re
	
	atRegex = re.compile('at', re.IGNORECASE)
	taRegex = re.compile('ta', re.IGNORECASE)
	
	
	possibleRealLocationIndices = []
	
	# This will cover almost all cases, except for a very small number.
	# I'm going to write a reporter function to shout out if it doesn't work.
	if (chromosomeStartIndex + 2 <= disruptionLocIndex <= chromosomeEndIndex - 2) \
	or (megaplasmidStartIndex + 2 <= disruptionLocIndex <= megaplasmidEndIndex - 2):
		
		searchStartIndex = disruptionLocIndex - 2
		searchEndIndex = disruptionLocIndex + 2
		
		i = searchStartIndex
		
		while i <= searchEndIndex:
			atMatch = atRegex.match(sequence[i:i+2])
			taMatch = taRegex.match(sequence[i:i+2])
	
			if atMatch != None:
				possibleRealLocationIndices.append([i, 'AT'])
			elif taMatch != None:
				possibleRealLocationIndices.append([i, 'TA'])

			i += 1
		
		# From the list of possible disruptions, choose the one closest to the measured
		# location
		if len(possibleRealLocationIndices) == 0:
			closestDisruptionIndex = [disruptionLocIndex, 'AT']
		elif len(possibleRealLocationIndices) > 0:
			differences = []
		
			for loc in possibleRealLocationIndices:
				differences.append(abs(loc[0] - disruptionLocIndex))
		
			minIndex = differences.index(min(differences))
			closestDisruptionIndex = possibleRealLocationIndices[minIndex]
		
		
	else:
		ex = GenomeSequenceIndexFailure('Genome Sequence Index Failure')
		print("Disruption location index " + str(disruptionLocIndex) + " is too close to " \
		+ "chromosome and megaplasmid edges.")
		raise ex 
	

	
	
	return closestDisruptionIndex
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def GenerateTrimmedGenomeSequence(closestDisruptionIndex, sequence, \
chromosomeStartIndex, chromosomeEndIndex, megaplasmidStartIndex, megaplasmidEndIndex, \
transposonSequence, rcTransposonSequence, startTrim=100, endTrim=100):

	trimStart = closestDisruptionIndex[0] - startTrim
	trimEnd = closestDisruptionIndex[0] + endTrim + 1
	
	if ((chromosomeStartIndex <= trimStart <= chromosomeEndIndex) \
	and (chromosomeStartIndex <= trimEnd <= chromosomeEndIndex)) \
	or ((megaplasmidStartIndex <= trimStart <= megaplasmidEndIndex) \
	and (megaplasmidStartIndex <= trimEnd <= megaplasmidEndIndex)):
	
		
		
		trimmedSequence = sequence[trimStart:closestDisruptionIndex[0]]
		
		if closestDisruptionIndex[1] == 'AT':
			trimmedSequence += transposonSequence
		elif closestDisruptionIndex[1] == 'TA':
			trimmedSequence += rcTransposonSequence
		
		trimmedSequence += sequence[closestDisruptionIndex[0]:trimEnd]
		
		sequenceTarget = [startTrim, len(transposonSequence)]
		productSizeRange = [len(transposonSequence), \
		int(len(transposonSequence) + startTrim/2 + endTrim/2)]
		
	
	else:
		ex = GenomeSequenceIndexFailure('Genome Sequence Index Failure')
		raise ex
	
	
	
	return trimmedSequence, sequenceTarget, productSizeRange
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportInsertionLocationFile(primerLocationFile):
	
	import re
	
	# Import the insertion location to check file
	locFileHandle = open(primerLocationFile, 'r')
	locData = locFileHandle.readlines()
	locDict = {}
	i = 1

	locationsRe = re.compile(r'(\d+),([A-H]\d+),\"*\d+')
	locRe = re.compile(r'\d+')

	while i < len(locData):
		line  = locData[i]
	
		lineData = line.strip().split(',')
		plate = lineData[0]
		well = lineData[1]
	
		locations = []
		j = 2
		while j < len(lineData):
			location = locRe.search(lineData[j]).group(0)
			#print(location)
			locations.append(int(location))
			j += 1
	
		wellKey = plate + '_' + well
		locDict[wellKey] = locations
		i += 1

	locFileHandle.close()
	
	return locDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ParseCSVLine(line):
	
	import pdb
	
	line = line.strip()
	
	i = 0
	columns = []
	currentColumn = ''
	inQuotedBlock = False
	
	while i < len(line):
		letter = line[i]
		
		if letter == '"':
			if inQuotedBlock == True:
				inQuotedBlock = False
			elif inQuotedBlock == False:
				inQuotedBlock = True
		
		
		elif letter == ',' and inQuotedBlock == False:
			columns.append(currentColumn)
			currentColumn = ''
		
		else:
			currentColumn += letter
		
		i += 1
	
	columns.append(currentColumn)
	
# 	pdb.set_trace()
		
	
	return columns
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ImportInsertionLocationFile2(primerLocationFile):
	
	import re
	import pdb
	
	# Import the insertion location to check file
	locFileHandle = open(primerLocationFile, 'r')
	locData = locFileHandle.readlines()
	locDict = {}
	i = 1

	
	locRe = re.compile(r'\d+')

	while i < len(locData):
		line  = locData[i]
	
		lineData = ParseCSVLine(line.strip())
		
		plate = lineData[0]
		well = lineData[1]
		locations = lineData[2]
		try:
			primerCodes = lineData[3]
		except:
			pdb.set_trace()
			
			
		wellKey = plate + '_' + well
		locDict[wellKey] = [locations, primerCodes]
		i += 1

	locFileHandle.close()
	
	return locDict
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def ReverseComplement(sequence):
	
	import pdb
	
	revComplement = ''
	
	i = len(sequence) - 1
	while i >= 0:
		
		sequenceLetter = sequence[i]
		
		if sequenceLetter == 'a' or sequenceLetter == 'A':
			revComplementLetter = 'T'
		elif sequenceLetter == 't' or sequenceLetter == 'T':
			revComplementLetter = 'A'
		elif sequenceLetter == 'c' or sequenceLetter == 'C':
			revComplementLetter = 'G'
		elif sequenceLetter == 'g' or sequenceLetter == 'G':
			revComplementLetter = 'C'
		
		revComplement += revComplementLetter
		
		
		i -= 1

	return revComplement
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ParsePrimer3DesignOutputFile(file) :
	
	import re
	
	primerLeft0Re = re.compile('PRIMER_LEFT_0_SEQUENCE=([a-zA-Z]+)', re.IGNORECASE)
	primerRight0Re = re.compile('PRIMER_RIGHT_0_SEQUENCE=([a-zA-Z]+)', re.IGNORECASE)
	
	primertLeft0TmRe = re.compile('PRIMER_LEFT_0_TM=(\d+\.\d+)')
	primertRight0TmRe = re.compile('PRIMER_RIGHT_0_TM=(\d+\.\d+)')
	
	primertLeft0GCRe = re.compile('PRIMER_LEFT_0_GC_PERCENT=(\d+\.\d+)')
	primertRight0GCRe = re.compile('PRIMER_RIGHT_0_GC_PERCENT=(\d+\.\d+)')
	
	
	fileHandle = open(file, 'r')
	data = fileHandle.readlines()
	
	i = 0
	while i < len(data):
	
		line = data[i].strip()
	
		primerLeft0Match = primerLeft0Re.search(line)
		primerRight0Match = primerRight0Re.search(line)

		primertLeft0TmMatch = primertLeft0TmRe.search(line)  
		primertRight0TmMatch = primertRight0TmRe.search(line)

		primertLeft0GCMatch = primertLeft0GCRe.search(line)
		primertRight0GCMatch = primertRight0GCRe.search(line)
		
		if primerLeft0Match != None:
			leftPrimer = primerLeft0Match.group(1)
		elif primerRight0Match != None:
			rightPrimer = primerRight0Match.group(1)
		elif primertLeft0TmMatch != None:
			leftPrimerTm = float(primertLeft0TmMatch.group(1))
		elif primertRight0TmMatch != None:
			rightPrimerTm = float(primertRight0TmMatch.group(1))
		elif primertLeft0GCMatch != None:
			leftPrimerGC = float(primertLeft0GCMatch.group(1))
		elif primertRight0GCMatch != None:
			rightPrimerGC = float(primertRight0GCMatch.group(1))
		
		
		i += 1

	
	
	return rightPrimer, leftPrimer, rightPrimerTm, leftPrimerTm, rightPrimerGC, leftPrimerGC
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
####################################################################################################



