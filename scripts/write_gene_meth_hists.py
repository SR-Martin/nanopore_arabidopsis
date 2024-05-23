#!/usr/bin/env python
import sys, getopt, errno

class GenePos:
	def __init__(self, start, end):
		self.start = start
		self.end = end


class Gene:
	def __init__(self, name, strand, chrom, start, end):
		self.name = name
		self.strand = strand
		self.chrom = chrom
		self.coords = []
		self.coords.append(GenePos(start, end))


	def __lt__(self, other):
		if self.chrom < other.chrom:
			return True
		elif self.chrom > other.chrom:
			return False
		elif self.strand == "+" and other.strand == "-":
			return True
		elif self.strand == "-" and other.strand == "+":
			return False
		else:
			return self.coords[0].start < other.coords[0].start

	def __str__(self):
		return self.name + "\t" + str(self.chrom) + "\t" + self.strand + "\t" + str(self.start) + "\t" + str(self.end)

	def appendCoords(self, start, end):
		self.coords.append(GenePos(start, end))

	def getBin(self, pos, n):
		genePos = 0
		length = 0
		for coord in self.coords:
			interval = coord.end - coord.start
			if pos >= coord.start and pos <= coord.end:
				genePos += pos - coord.start
			elif pos > coord.end:
				genePos += interval
			length += interval
		return int(n * float(genePos)/length)

numBins = 100
binSize = 10
methycallFilename = "../new_analysis/LSK109/CpG_5mC_counts_all.methylCall"

geneNames = []
genes = dict()

try:
	opts, args = getopt.getopt(sys.argv[1:],"hg:f:n:m:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python write_gene_meth_hists.py -g <gene_list.txt> -f <gene_sites.fasta> -n <bin_size> -m <methylation file>")
	print("python write_gene_meth_hists.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("python write_gene_meth_hists.py -g <gene_list.txt> -f <gene_sites.fasta> - n <bin_size>")
		print("-g <gene_list.txt>\t\t List of genes")
		print("-f <gene_sites.fasta>\t\t Fasta file of gene sites")
		print("-n <bin_size>\t\t Size of overlapping bins (%) for histogram")
		print("-m <methylation file>\t\t File containing methylation calls")
		sys.exit()
	elif opt in ("-g"):
		geneListFilename = arg
	elif opt in ("-f"):
		fastaFilename = arg
	elif opt in ("-n"):
		binSize = int(arg)
	elif opt in ("-m"):
		methycallFilename = arg


methylatedFreqBins = [.0] * numBins
totalSites = [0] * numBins

if not geneListFilename or not fastaFilename:
	print("You must specify -g and -f")
	sys.exit(2)


with open(geneListFilename, 'r') as infile:
	for line in infile:
		geneNames.append(line.strip())

with open(fastaFilename, 'r') as infile:
	lastGeneName = ""
	for line in infile:
		if line[0] == '>':
			fields = line.split("|")
			if "." in fields[0]:
				geneName = fields[0][1:].split(".")[0]
			else:
				geneName = fields[0][1:].strip()

			if geneName in geneNames:
				coordString = fields[-1].strip()
				chromString = coordString.split(":")[0]
				chrom = int(chromString[3])
				coords = (coordString.split()[0]).split(":")[1]
				start = int(coords.split("-")[0])
				end = int (coords.split("-")[1])
				strand = "+"
				if coordString.split()[1] == "REVERSE":
					strand = "-"

				if chrom not in genes.keys():
					genes[chrom] = list()

				if lastGeneName != geneName:
					lastGeneName = geneName
					genes[chrom].append(Gene(geneName, strand, chromString[3], start, end))
				else:
					genes[chrom][-1].appendCoords(start, end)


for key in genes.keys():
	genes[key].sort()
	#for gene in genes[key]:
	#	print(gene)

with open(methycallFilename, 'r') as infile:
	currentChrom = min(genes.keys())
	currentStrand = "+"
	geneCount = 0
	geneList = genes[currentChrom]
	infile.readline()

	line = infile.readline()
	fields = line.split()
	pos = int(fields[1])
	strand = fields[2]

	for key in genes.keys():
		geneList = genes[key]
		for gene in geneList:
			geneStart = gene.coords[0].start
			geneEnd = gene.coords[-1].end
			while pos < geneStart or strand != gene.strand or key != chrom:
				line = infile.readline()
				fields = line.split()
				pos = int(fields[1])
				strand = fields[2]
				chrom = int(fields[0])

			while pos < geneEnd and strand == gene.strand and key == chrom:
				for coord in gene.coords:
					if pos >= coord.start and pos <= coord.end:
						methylatedFreq = float(fields[7])
						lastBin = gene.getBin(pos, numBins)
						firstBin = max(lastBin - binSize, 0)
						for i in range(firstBin, lastBin + 1):
							methylatedFreqBins[i] += methylatedFreq
							totalSites[i] += 1

				line = infile.readline()
				fields = line.split()
				pos = int(fields[1])
				strand = fields[2]
				chrom = int(fields[0])


for b in range(numBins):
	if totalSites[b] > 0:
		print(str(b) + "\t" + str(methylatedFreqBins[b]/totalSites[b]))
	else:
		print("0\t0")















