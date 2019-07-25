import sys, getopt, re, os

def generateContentFile(argv):
	
	inPath = ''
	taxonomyFile = ''
	nameFile = ''
	accToTaxFile = ''
	contentFile = ''
	taxonomicUnit = ''
	
	try:
		opts, args = getopt.getopt(argv, "i:c:a:t:n:o:u:", [])
	except getopt.GetoptError:
		print('you failed')
	
	for opt, arg in opts:
		if opt in ("-i",):
			inPath = arg
		elif opt in ("-c",):
			contentFile = open(arg, 'w')
		elif opt in ("-a",):
			accToTaxFile = open(arg)
		elif opt in ("-t",):
			nodesFile = open(arg)
		elif opt in ("-n",):
			nameFile = open(arg)
		elif opt in ("-u",):
			taxonomicUnit = arg
	
	print("go through fastas and gather all accs")
	listOfAccs = []
	listOfFiles = []
	if inPath[-1] == '/':
		listOfFiles = os.listdir(inPath)
	else:
		listOfFiles.append("")
	
	dummys = []
	for filename in listOfFiles:
		fastaFile = open(inPath+filename)
		for line in fastaFile:
			if '>' in line:
				try:
					acc = ""
					if " " in line:
						numbers = ((((line.rstrip("\n")).split(' '))[0]).lstrip('>')).split('|')
						for e in numbers:
							if "." in e:
								acc = e
								break
					else:
						acc = (line.rstrip("\n")).lstrip(">")
					if acc != "":
						listOfAccs.append(acc)
					else:
						print("No accession number found in: ", line, "Providing dummy ID later..")
						dummys.append((line.rstrip("\n")).lstrip('>'))
				except:
					print(line)
					print(inPath+filename)
					fastaFile.close()
					break
		fastaFile.close()
	
	
	listOfAccs = list(set(listOfAccs))
	print(listOfAccs,dummys)
	
	print(str(len(listOfAccs))+" accession numbers found")
	
	
	print("sorting and translating to taxIDs")
	listOfAccs.sort()
	listOfAccsAndTax = []
	iIndexOflistOfAccs = 0
	for line in accToTaxFile:
		if iIndexOflistOfAccs >= len(listOfAccs):
			break
		line = line.split('\t')
		
		
		while listOfAccs[iIndexOflistOfAccs] < line[1]:
			iIndexOflistOfAccs += 1
			if iIndexOflistOfAccs >= len(listOfAccs):
				#print(listOfAccs[iIndexOflistOfAccs-1], line[1])
				break
		
		if iIndexOflistOfAccs >= len(listOfAccs):
			break
		
		if listOfAccs[iIndexOflistOfAccs] == line[1]:
			listOfAccsAndTax.append( (listOfAccs[iIndexOflistOfAccs], (line[2]).rstrip('\n\r')) )
			iIndexOflistOfAccs += 1
			if iIndexOflistOfAccs >= len(listOfAccs):
				break
	accToTaxFile.close()
	
	
	
	if len(listOfAccsAndTax) != len(listOfAccs):
		print("The following accession numbers have no tax ID:")
		listIdx = 0
		for acc in listOfAccs:
			if acc != listOfAccsAndTax[listIdx][0]:
				print(acc)
				dummys.append(acc)
			else:
				listIdx += 1
	
	
	print("Providing dummy tax IDs...")
	pool = 2147483646
	dummyTaxIDs = []
	if len(dummys) > pool:
		print("Too many dummys, I can't handle this!")
		return
	for entries in dummys:
		dummyTaxIDs.append( (entries,pool) )
		pool -= 1
	
	print("create a dictionary of the names")
	nameDict = {}
	for line in nameFile:
		line = line.split('|')
		if line[3] == "\t" + "scientific name\t":
			nameDict[line[0].rstrip('\t')] = (line[1].rstrip('\t')).lstrip('\t')
		#6	|	Azorhizobium	|		|	scientific name	|
		#6	|	Azorhizobium Dreyfus et al. 1988 emend. Lang et al. 2013	|		|	authority	|
		#6	|	Azotirhizobium	|		|	misspelling	|
	nameFile.close()
	
	print("create dictionary for lookup to map taxids to their genus")
	nodesDict = {}
	for line in nodesFile:
		line = line.split('|')
		#2	|	131567	|	superkingdom	|		|	0	|	0	|	11	|	0	|	0	|	0	|	0	|	0	|		|
		nodesDict[(line[0]).rstrip('\t')] = ( ((line[1]).rstrip('\t')).lstrip('\t') , ((line[2]).rstrip('\t')).lstrip('\t') )
	nodesFile.close()
	
	
	print("link at taxonomicUnit level")
	tempList = {}
	debugIter = 0
	for (acc, tax) in listOfAccsAndTax:
		debugIter += 1
		upperTax = tax
		entry = ""
		if upperTax in nodesDict:
			entry = nodesDict[upperTax]
		else:
			entry = ( "1", "" )
		#print entry
		if taxonomicUnit != '':
			while entry[1] != taxonomicUnit and entry[0] != "1":
				upperTax = entry[0]
				entry = nodesDict[upperTax]
				#print entry
			if entry[0] == "1": #no genus to be found
				upperTax = tax
		else:
			upperTax = tax
		if upperTax in tempList:
			tempList[upperTax] = (tempList[upperTax][0] + ";" + tax , tempList[upperTax][1] + ";" + acc)
		else:
			tempList[upperTax] = (tax , acc)
	
	
	print("create contentFile")
	unnamedCounter = 0
	tempListSize = len(tempList)
	tempListIndex = 0
	for entry in tempList:
		nonredundantTaxList = list(set((tempList[entry][0]).split(';')))
		nonredundantEntry = nonredundantTaxList[0]
		for i in range(1,len(nonredundantTaxList)):
			nonredundantEntry += ";" + nonredundantTaxList[i]
		if entry in nameDict:
			contentFile.write(nameDict[entry] + "\t" + entry + "\t" + nonredundantEntry + "\t" + tempList[entry][1] + "\n")
		else:
			contentFile.write("unnamed_"+str(unnamedCounter)+ "\t" + entry + "\t" + nonredundantEntry + "\t" + tempList[entry][1] + "\n")
			unnamedCounter += 1
		tempListIndex += 1
	
	dummyCounter = 0
	for entry in dummyTaxIDs:
		contentFile.write("Dummy_"+ str(dummyCounter) + "\t" + str(entry[1]) + "\t" + str(entry[1]) + "\t" + entry[0] + "\n")
		dummyCounter += 1
	contentFile.close()
	
	
#end of generateContentFile

generateContentFile(sys.argv[1:])