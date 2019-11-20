import sys, getopt, os

def humanReadableReadQuantity(argv):
	threshold = 0.0
	
	try:
		opts, args = getopt.getopt(argv, "i:o:t:", [])
	except getopt.GetoptError:
		print('-i <input> -o <output>')
	
	for opt, arg in opts:
		if opt in ("-i",):
			kASAOutput = open(arg)
		elif opt in ("-o",):
			outFile = open(arg, 'w')
		elif opt in ("-t",):
			threshold = float(arg)
	
	resultDict = {}
	readCount = 0
	
	next(kASAOutput)
	for read in kASAOutput:
		read = (read.rstrip("\r\n")).split("\t")
		readCount += 1

		if read[2] == "-":
			continue
		taxIDs = (read[2]).split(";")
		names = (read[3]).split(";")
		scores = (read[4]).split(";")
		
		until = 0
		startingScore = float((scores[0]).split(",")[0])
		if startingScore < threshold:
		 continue
		
		for i in range(len(taxIDs)):
			if float((scores[i]).split(",")[0]) >= startingScore:
				until += 1
			else:
				break
		
		for i in range(until):
			taxID = taxIDs[i]
			if taxID in resultDict:
				res = resultDict[taxID]
				resultDict[taxID] = (res[0], res[1] + 1.0/until)
			else:
				resultDict[taxID] = ( names[i], 1.0/until )
	
	resultList = []
	
	for entry in resultDict:
		resultList.append((entry , (resultDict[entry])[0] , (resultDict[entry])[1] , (resultDict[entry])[1] / readCount))
	
	resultList.sort(key=lambda x: int(x[2]), reverse=True)
	
	for entry in resultList:
		outFile.write(entry[0] + "\t" + entry[1] + "\t" + str(entry[2]) + "\t" + str(entry[3]) + "\n")
	

humanReadableReadQuantity(sys.argv[1:])