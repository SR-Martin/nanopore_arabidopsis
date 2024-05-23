#!/usr/bin/python

import sys, getopt, math

def getReadQScore(qualityString):
	totalScore = 0
	qualityString = qualityString.strip()
	for i in xrange(0, len(qualityString)):
		qscore = ord(qualityString[i]) - 33
		totalScore += math.pow(10,-qscore/10.)
	meanScore = totalScore / len(qualityString)
	return -10 * math.log10(meanScore)

def writeGoodScores(filename_in, filename_out, minScore):
	infile = open(filename_in, 'r')
	outfile = open(filename_out, 'w+')
	numReads = 0
	numReadsKept = 0
	firstLine 	= infile.readline()
	secondLine 	= infile.readline()
	thirdLine 	= infile.readline()
	fourthLine 	= infile.readline()

	while fourthLine != '':
		numReads += 1
		assert(thirdLine.strip() == "+")
		if getReadQScore(fourthLine) >= minScore:
			numReadsKept += 1
			outfile.write(firstLine)
			outfile.write(secondLine)
			outfile.write(thirdLine)
			outfile.write(fourthLine)
		firstLine 	= infile.readline()
		secondLine 	= infile.readline()
		thirdLine 	= infile.readline()
		fourthLine 	= infile.readline()

	infile.close()
	outfile.close()
	print "Removed " + str(numReads - numReadsKept) + "/" + str(numReads) + " reads."

inputfile = ''
outputfile = ''
quality = 9
try:
	opts, args = getopt.getopt(sys.argv[1:],"hq:i:o:",["ifile=","ofile="])
except getopt.GetoptError:
	print 'Option not recognised.'
	print 'QScoreFilter.py -i <inputfile> -o <outputfile> -q <minquality>'
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print 'QScoreFilter.py -i <inputfile> -o <outputfile> -q <minquality>'
		sys.exit()
	elif opt in ("-i", "--ifile"):
		inputfile = arg
	elif opt in ("-o", "--ofile"):
		outputfile = arg
	elif opt in ("-q"):
		quality = arg

if inputfile == '' or outputfile =='':
	print 'Please provide an input file and an output file.'
	print 'QScoreFilter.py -i <inputfile> -o <outputfile> -q <minquality>'
	sys.exit(2)

print "Filtering out reads with Q score less than " + str(quality)
writeGoodScores(inputfile, outputfile, int(quality))