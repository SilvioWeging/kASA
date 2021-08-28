import sys, getopt, os, json

def jsonReadQuantity(argv):
	threshold = 0.0
	
	try:
		opts, args = getopt.getopt(argv, "i:o:t:", [])
	except getopt.GetoptError:
		print('-i <input> -o <output>')
	
	for opt, arg in opts:
		if opt in ("-i",):
			kASAOutput = json.load(open(arg))
		elif opt in ("-o",):
			outFile = open(arg, 'w')
		elif opt in ("-t",):
			threshold = float(arg)
	
	resultDict = {}
	readCount = 0
	
	
	for read in kASAOutput:
		readCount += 1
		taxa = read["Top hits"]
		#print(),read["Read number"],read["Specifier from input file"])
		#return
		if len(taxa) == 0:
			continue
		#taxIDs = (line[2]).split(";")
		#namesAndScores = (line[3]).split(";")
		#scores = (line[4]).split(";")
		
		until = 0
		startingScore = taxa[0]["Relative Score"]
		if startingScore < threshold:
		 continue
		
		for elem in taxa:
			if elem["Relative Score"] >= startingScore:
				until += 1
			else:
				break
		
		for i in range(until):
			taxID = taxa[i]["tax ID"]
			if taxID in resultDict:
				res = resultDict[taxID]
				resultDict[taxID] = (res[0], res[1] + 1.0/until)
			else:
				resultDict[taxID] = ( taxa[i]["Name"], 1.0/until )
	
	resultList = []
	
	for entry in resultDict:
		resultList.append((entry , (resultDict[entry])[0] , (resultDict[entry])[1] , (resultDict[entry])[1] / readCount))
	
	resultList.sort(key=lambda x: int(x[2]), reverse=True)
	
	for entry in resultList:
		outFile.write(entry[0] + "\t" + entry[1] + "\t" + str(entry[2]) + "\t" + str(entry[3]) + "\n")
	

jsonReadQuantity(sys.argv[1:])