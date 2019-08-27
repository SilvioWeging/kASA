import sys, getopt, os, json

def jasonToCAMIBin(argv):
	try:
		opts, args = getopt.getopt(argv, "i:o:", [])
	except getopt.GetoptError:
		print('-i <input> -o <output>')
	
	for opt, arg in opts:
		if opt in ("-i",):
			kASAOutput = json.load(open(arg))
		elif opt in ("-o",):
			outFile = open(arg, 'w')
	
	outFile.write("#CAMI Format for Binning created from kASA json output\n@Version:0.9.0\n@SEQUENCEID\tTAXID")
	
	for read in kASAOutput:
		taxa = read["Matched taxa"]
		outFile.write("\n" + read["Specifier from input file"] + "\t" + taxa[0]["tax ID"])

jasonToCAMIBin(sys.argv[1:])