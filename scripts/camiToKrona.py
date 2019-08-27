import sys, getopt, os

def camiToKrona(argv):
	
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
	
	for line in camiOutput:
		if "@" in line or "#" in line:
			continue
		line = line.rstrip("\n")
		if line == "":
			continue
		line = line.split("\t")
		outFile.write(line[4] + "\t" + line[3].replace("|","\t") + "\n")
	outFile.close()

camiToKrona(sys.argv[1:])