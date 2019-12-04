import sys, getopt, os

def camiToKrona(argv):
	
	rankArr = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain", "dummy"]
	camiOutput = ''
	outFile = ''
	
	try:
		opts, args = getopt.getopt(argv, "i:o:", [])
	except getopt.GetoptError:
		print('-i <input> -o <output>')
	
	for opt, arg in opts:
		if opt in ("-i",):
			camiOutput = open(arg)
		elif opt in ("-o",):
			outFile = open(arg, 'w')
	
	#look for smallest taxonomic rank
	rank = ""
	rankIdx = 0
	for line in camiOutput:
		if "@" in line or "#" in line:
			continue
		if rankArr[rankIdx] in line:
			rank = rankArr[rankIdx]
		else:
			rankIdx += 1
	
	sumOfValues = 0.0
	camiOutput.seek(0)
	for line in camiOutput:
		if "@" in line or "#" in line:
			continue
		line = line.rstrip("\n")
		if line == "":
			continue
		line = line.split("\t")
		if rank == line[1]:
			sumOfValues += float(line[4])
			outFile.write(line[4] + "\t" + line[3].replace("|","\t") + "\n")
		else:
			outFile.write("0.0" + "\t" + line[3].replace("|","\t") + "\n")
	outFile.write(str(100 - sumOfValues))
	outFile.close()

camiToKrona(sys.argv[1:])