import sys, getopt, os

def sumFreqsOnTaxLvl(argv):

	kASAOutput = ''
	nodes = ''
	outFile = ''
	nameFile = ''
	rank = ''
	
	#python /home/weging/scripts/fqkASASP.py -k /work-local/weging/out/bacteriaTestSetART_8.rtg -m /work/weging/bacteriaTestSet.txt -f /gpfs1/work/weging/bacteriaTestSetART.fq
	try:
		opts, args = getopt.getopt(argv, "i:n:m:r:o:", [])
	except getopt.GetoptError:
		print('you failed')
	
	for opt, arg in opts:
		if opt in ("-i",):
			frequencyFile = open(arg)
		elif opt in ("-n",):
			nodes = open(arg)
		elif opt in ("-m",):
			nameFile = open(arg)
		elif opt in ("-r",):
			rank = arg
		elif opt in ("-o",):
			outFile = open(arg, 'w')
	
	
	nodesDict = {}
	for line in nodes:
		line = line.split('|')
		#2	|	131567	|	superkingdom	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
		nodesDict[(line[0]).rstrip('\t')] = ( ((line[1]).rstrip('\t')).lstrip('\t') , ((line[2]).rstrip('\t')).lstrip('\t') )
	nodes.close()
	
	nameDict = {}
	for line in nameFile:
		line = line.split('|')
		if line[3] == "\t" + "scientific name\t":
			nameDict[line[0].rstrip('\t')] = (line[1].rstrip('\t')).lstrip('\t')
		#6	|	Azorhizobium	|		|	scientific name	|
		#6	|	Azorhizobium Dreyfus et al. 1988 emend. Lang et al. 2013	|		|	authority	|
		#6	|	Azotirhizobium	|		|	misspelling	|
	nameFile.close()
	
	#frequencyFile.next()
	
	resultDict = {}
	for line in frequencyFile:
		line = line.rstrip("\r\n")
		if "" == line:
			continue
		line = line.split('\t')
		specID = line[0]
		quantity = float(line[3])
		
		#print(specID, quantity)
		#break
		
		if specID in nodesDict:
			nextID = nodesDict[specID][0]
			nextRank = nodesDict[specID][1]
			#testParam = 0
			while nextRank != rank and nextID != "1":
				specID = nextID
				nextID = nodesDict[specID][0]
				nextRank = nodesDict[specID][1]
			
				#print specID, nextRank
				#if testParam == 5:
				#	return
				#testParam += 1
			
			if specID in resultDict:
				resultDict[specID] += quantity
			else:
				resultDict[specID] = quantity
		#else:
			#print specID, quantity
	#print(resultDict)
	for entry in resultDict:
		outStr = nameDict[entry] + "\t" + nodesDict[entry][1]
		outStr += "\t" + str(resultDict[entry]) +"\n"
		outFile.write(outStr)
	outFile.close()

#end of sumFreqsOnTaxLvl

sumFreqsOnTaxLvl(sys.argv[1:])