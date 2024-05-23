import sys, getopt, errno

def printHelp():
	print("bed_to_methylcall.py -i <bed TSV> -o <methylation call TSV")

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:")
except getopt.GetoptError:
	print("Option not recognised.")
	printHelp()
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		printHelp()
		sys.exit()
	elif opt in ("-i"):
		inputfile = arg
	elif opt in ("-o"):
		outputfile = arg

if inputfile == '' or outputfile == '':
	print("You must specify -i and -o")
	printHelp()
	sys.exit(2)

class MethylCall:
	def __init__(self, pos, strand, chrom, methylatedBases, totalBases):
		self.strand = strand
		self.chr = chrom
		self.pos = pos + 1
		self.methylatedBases = methylatedBases
		self.totalBases = totalBases

chromosome = dict()
for strand in ["+", "-"]:
	chromosome[strand] = list()

try:
	with open(inputfile, 'r') as infile:
		for line in infile:
			fields = line.split("\t")
			chrName = fields[0]
			if chrName[:3] == "Chr":
				chrNumber = chrName[-1]
				pos = int(fields[1]) - 1				
				strand = fields[5]
				methylatedBases = int(fields[12])
				totalBases = int(fields[11]) + int(fields[12]) + int(fields[15])
				if totalBases > 0:
					chromosome[strand].append(MethylCall(pos,strand,chrNumber,methylatedBases,totalBases))

except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + inputfile)
	sys.exit(2)

with open(outputfile, 'w') as outfile:
	outfile.write("chrom\tpos\tstrand\tmc_class\tmethylated_bases\ttotal_bases\tmethylation_call\tmethyl_freq\n")

	for strand in ["+", "-"]:
		methylCalls = chromosome[strand] 
		for mc in methylCalls:
			freq = float(mc.methylatedBases) / mc.totalBases
			call = 0
			if freq > 0.5 or (freq == 0.5 and mc.totalBases > 2):
				call = 1
			outfile.write(str(mc.chr) + "\t" + str(mc.pos) + "\t" + strand + "\tCCG\t" + str(mc.methylatedBases) + "\t" + str(mc.totalBases) + "\t" + str(call) + "\t" + str(freq) + "\n")





