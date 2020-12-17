import urllib.request
import sys, getopt
import os.path


try:
	opts, args = getopt.getopt(sys.argv[1:], "i:o:", [])
except getopt.GetoptError:
	print('-i <content file> -o <output path>')

for opt, arg in opts:
	if opt in ("-i",):
		contentFile = open(arg)
	elif opt in ("-o",):
		outPath = arg

for line in contentFile:
	line = line.rstrip("\r\n")
	if line == "":
		continue
	accnrs = (line.split("\t"))[3]
	for accnr in accnrs.split(";"):
		if accnr != "":
			if os.path.isfile(outPath + accnr + ".fasta"):
				print("File already exists")
				continue
			print("Downloading file: ", accnr + ".fasta")
			outfileForGenome = open(outPath + accnr + ".fasta", 'wb')
			callToNCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + accnr + "&rettype=fasta&retmode=text"
			outfileForGenome.write((urllib.request.urlopen(callToNCBI)).read())#genome
			outfileForGenome.close()
