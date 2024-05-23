#!/usr/bin/env python

import sys, getopt, errno
import time
import glob
import math
start = time.time()

def ReadFastaFile(_inputFile):
	header = _inputFile.readline()
	if header[0] != '>':
		print("Error: " + _inputFile.name + " does not appear to be in fasta format.")
		sys.exit(2)

	contigs = list()
	line = _inputFile.readline()
	read = ''
	while line != '':
		if line[0] != '>':
			read += line.strip()
		else:
			contigs.append(read.strip())
			read = ''
		line = _inputFile.readline()
	contigs.append(read.strip())
	return contigs

def ReadFastqFile(_inputFile, _makePlots):
	header = _inputFile.readline()
	if header[0] != '@':
		print("Error: " + _inputFile.name + " does not appear to be in fastq format.")
		sys.exit(2)	
	contigs = list()

	#initalize _qualityDistribution array
	baseQualityDistribution = [0] * 100
	readQualityDistribution = [0] * 100

	totalBaseQuality = 0
	totalReadQuality = 0
	totalReadAccuracy = 0
	while header != '':
		read 	= _inputFile.readline()
		qHeader = _inputFile.readline()
		quality = _inputFile.readline()

		assert qHeader.strip()[0] == '+'
		contigs.append(read.strip())

		quality = quality.strip()
		readAccuracy = 0
		for char in quality:
			score = ord(char) - 33
			accuracy = pow(10, -score/10.0)
			readAccuracy += accuracy
			if _makePlots:
				baseQualityDistribution[score] += 1

		avgReadAccuracy = -10 * math.log(float(readAccuracy) / len(quality), 10)
		totalReadAccuracy += avgReadAccuracy
		readQualityDistribution[int(avgReadAccuracy)] += 1
		header 	= _inputFile.readline()

	return contigs, baseQualityDistribution, readQualityDistribution, totalReadAccuracy

inputfile = ''
isFastq = False
makePlots = False
countBases = False
ACount = 0
CCount = 0
GCount = 0
TCount = 0
NCount = 0
GCComposition = 0
ATComposition = 0
NComposition = 0
genomeSize = 0
doCoverages = False
numge500 = 0
histPrefix = ""

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:qg:bp:",["ifile="])
except getopt.GetoptError:
	print("Option not recognised.")
	print("python get_stats.py -i <inputfile>")
	print("python get_stats.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("get_contig_stats.py -i <inputfile> -g <genomesize> - q -b -p")
		print("-i <inputfile>\t\t Name of fasta or fastq file to analyse")
		print("-g <genomesize>\t\t Size of genome to compare to for coverage stats. Use M for megabases, G for gigabases.")
		print("-q             \t\t Specify that input is in fastq format.")
		print("-b             \t\t Get base counts and GC composition (can be slow!).")
		print("-p <prefix>    \t\t Write length and quality histograms.")
		sys.exit()
	elif opt in ("-i", "--ifile"):
		inputfile = arg
	elif opt in("-q"):
		isFastq = True
	elif opt in ("-g"):
		doCoverages = True
		i = 0
		numString = ""
		unit = ""
		for char in arg:
			if not char.isdigit():
				numString = arg[:i]
				unit = arg[i:]
			i += 1
		unit = unit.upper()
		if unit == "":
			genomeSize = int(arg)
		elif unit == 'M':
			genomeSize = int(numString) * 1000
		elif unit == 'G':
			genomeSize = int(numString) * 1000000
		else:
			print("Did not understand unit " + unit + ". Use M or G.")

	elif opt in("-p"):
		makePlots = True
		histPrefix = arg
	elif opt in("-b"):
		countBases = True


if inputfile == '':
	print("You must specify -i")
	print("get_contig_stats.py -i <inputfile>")
	sys.exit(2)

s = inputfile.split(".")
prefix = ''
for i in range(0, len(s)-1):
	prefix += s[i]
if s[-1] == "fastq":
	isFastq = True

allContigs = list()
allTotalReadAccuracy = 0
allBaseQualityDistribution = [0] * 100
allReadQualityDistribution = [0] * 100

filenames = glob.glob(inputfile)
if len(filenames) > 1:
	prefix = "all"

for filename in filenames:
	try:
		with open(filename, 'r') as infile:
			if isFastq:
				contigs, baseQualityDistribution, readQualityDistribution, totalReadAccuracy = ReadFastqFile(infile, makePlots)
				allTotalReadAccuracy += totalReadAccuracy
				if makePlots:
					for i in range(0, len(baseQualityDistribution)):
						allBaseQualityDistribution[i] += baseQualityDistribution[i]
						allReadQualityDistribution[i] += readQualityDistribution[i]
			else:
				contigs = ReadFastaFile(infile)
			allContigs.extend(contigs)
	except (OSError, IOError) as e: 
		if getattr(e, 'errno', 0) == errno.ENOENT:
			print("Could not find file " + inputfile)
		print(e)
		sys.exit(2)

allContigs.sort(key=len, reverse=True)

numContigs = len(allContigs)
longestContig = 0
shortestContig = 0
if len(allContigs) > 0:
	longestContig = len(allContigs[0])
	shortestContig = len(allContigs[numContigs - 1])
contigSum = 0
lengths = list()

for contig in allContigs:
	contigSum += len(contig)
	if len(contig) >= 500:
		numge500 += 1
	if makePlots:
		lengths.append(len(contig))
	if countBases:
		for char in contig:
			if char == 'A':
				ACount += 1
			elif char == 'C':
				CCount += 1
			elif char == 'G':
				GCount += 1
			elif char == 'T':
				TCount += 1
			else:
				NCount += 1

averageContigLength = float(contigSum) / numContigs
if isFastq:
	averageReadAccuracy = float(allTotalReadAccuracy) / numContigs
if countBases:
	GCComposition = float(GCount + CCount) / contigSum
	ATComposition = float(ACount + TCount) / contigSum
	NComposition = float(NCount) / contigSum

#calculate N50 and N90
if doCoverages:
	numCoverages = contigSum / genomeSize
	genomeCoverages = [None] * (int)(numCoverages)
	genomeCoverageLengths = [None] * (int)(numCoverages)
N50 		= 0
N90 		= 0
N50Count 	= 0
N90Count 	= 0
N50Length = contigSum * 0.5
N90Length = contigSum * 0.9
cumulativeLength = 0
count = 0
i = 1
median = 0
even = False
if len(allContigs) % 2 == 0:
	even = True

for contig in allContigs:
	cumulativeLength += len(contig)
	count += 1
	if N90 == 0 and cumulativeLength >= N90Length:
		N90 = len(contig)
		N90Count = count
	if N50 == 0 and cumulativeLength >= N50Length:
		N50 = len(contig)
		N50Count = count
	if doCoverages:
		if cumulativeLength >= i * genomeSize:
			genomeCoverages[i-1] = count
			genomeCoverageLengths[i-1] = len(contig)
			i+=1
	if not even and count == int((len(allContigs) + 1)/2):
		median = len(contig)
	if even and count == int(len(allContigs)/2):
		median += len(contig)
	if even and count == int((len(allContigs)/2) + 1):
		median += len(contig)

if even:
	median = float(median)/2


print("NumContigs:\t" 			+ str(numContigs))
print("TotalSum:\t" 			+ str(contigSum))
print("MeanLength:\t"			+ str(averageContigLength))
print("MedianLength:\t"			+ str(median))
print("Shortest:\t"				+ str(shortestContig))
print("Longest:\t"				+ str(longestContig))
print("N50Length:\t"			+ str(N50))
print("N50Count:\t"				+ str(N50Count))
print("N90Length:\t"			+ str(N90))
print("N90Count:\t"				+ str(N90Count))
if isFastq:
	print("Mean read accuracy:\t" + str(averageReadAccuracy))
print("NumContigs >= 500: \t" 	+ str(numge500))
if countBases:
	print("No. A:\t\t"			+ str(ACount))
	print("No. T:\t\t"			+ str(TCount))
	print("No. C:\t\t"			+ str(CCount))
	print("No. G:\t\t"			+ str(GCount))
	print("No. N:\t\t"			+ str(NCount))
	print("GC-Composition:\t"	+ str(GCComposition))
	print("AT-Composition:\t"	+ str(ATComposition))
	print("N-Composition:\t"	+ str(NComposition))

if doCoverages:
	print("Coverage counts for genome of size " + str(genomeSize) + "b:")
	print("\tCoverage\tContig Count\tContig Length")
	for i in range(0, len(genomeCoverages)):
		print("\t" + str(i+1) + "x\t\t" + str(genomeCoverages[i]) + "\t\t" + str(genomeCoverageLengths[i]))

if makePlots:
	with open(histPrefix + ".lengths.txt", 'w') as outfile:
		for length in lengths:
			outfile.write(str(length) + "\n")

	if isFastq:
		with open(histPrefix + ".basequality.txt", 'w') as outfile:
			for i in range(0, len(allBaseQualityDistribution)):
				outfile.write(str(i) + "\t" + str(allBaseQualityDistribution[i]) + "\n")
		with open(histPrefix + ".readquality.txt", 'w') as outfile:
			for i in range(0, len(allReadQualityDistribution)):
				outfile.write(str(i) + "\t" + str(allReadQualityDistribution[i]) + "\n")

end = time.time()
print("get_contig_starts.py end time: " + str(end - start))

