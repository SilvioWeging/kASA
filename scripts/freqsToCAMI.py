import sys, getopt, os

def csvToCAMI(argv):
	kASAOutput = ''
	nodes = ''
	outFile = ''
	nameFile = ''
	threshold = 0.0
	
	try:
		opts, args = getopt.getopt(argv, "i:n:m:o:t:", [])
	except getopt.GetoptError:
		print('Wrong input!')
	
	for opt, arg in opts:
		if opt in ("-i",):
			kASAOutput = open(arg)
		elif opt in ("-n",):
			nodes = open(arg)
		elif opt in ("-m",):
			nameFile = open(arg)
		elif opt in ("-o",):
			outFile = open(arg, 'w')
		elif opt in ("-t",):
			threshold = float(arg)

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

	taxPaths = {}
	for line in kASAOutput:
		line = line.rstrip("\r\n")
		if line == "":
			break
		line = line.split('\t')
		specID = line[0]
		quantity = float(line[3])*100.0
		
		if quantity > threshold:
			taxIDPath = specID
			if specID in nameDict:
				taxNamePath = nameDict[specID]
			else:
				taxNamePath = "unnamed"
			
			# go upwards the taxonomic tree
			if specID in nodesDict:
				taxRankPath = nodesDict[specID][1]
				nextID = nodesDict[specID][0]
				currRank = nodesDict[specID][1]
				
				while currRank != "superkingdom" and specID != "1":
					specID = nextID
					currRank = nodesDict[specID][1]
					if currRank != "no rank":
						taxIDPath = nextID + "|" + taxIDPath
						if nextID in nameDict:
							taxNamePath = nameDict[nextID] + "|" + taxNamePath
						else:
							taxNamePath = "unnamed" + "|" + taxNamePath
					else:
						taxIDPath = "|" + taxIDPath
						taxNamePath = "|" + taxNamePath
					taxRankPath = currRank + "|" + taxRankPath
					nextID = nodesDict[specID][0]

				taxIDs = taxIDPath.split("|")
				taxNames = taxNamePath.split("|")
				taxRanks = taxRankPath.split("|")

				# sum frequencies
				for i in range(len(taxIDs)-1, -1, -1):
					tID = taxIDs[i]
					if tID != "":
						if tID in taxPaths:
							taxPath = (taxPaths[tID]).split("\t")
							newQuantitity = taxPath[4]
							taxPaths[tID] = taxPath[0] + "\t" + taxPath[1] + "\t" + taxPath[2] + "\t" + taxPath[3] + "\t" + str(float(newQuantitity) + quantity)
						else:
							shortenedIDPath = taxIDs[0]
							shortenedNamePath = taxNames[0]
							for j in range(1,i+1,1):
								shortenedIDPath += "|" + taxIDs[j]
								shortenedNamePath += "|" + taxNames[j]
							taxPaths[tID] = tID + "\t" + taxRanks[i] + "\t" + shortenedIDPath + "\t" + shortenedNamePath + "\t" + str(quantity)
			#else:
			#	print specID
	outFile.write("#CAMI Submission for Taxonomic Profiling\n@SampleID:\n@Version:0.9.2\n@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n@TaxonomyID:?\n@__program__:kASA\n@@TAXID	RANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
	rankArr = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
	while len(rankArr) > 0:
		for path in taxPaths:
			splitPath = (taxPaths[path]).split("\t")
			amount = splitPath[4]
			rank = splitPath[1]
			if rank == rankArr[0]:
				outFile.write(taxPaths[path] + "\n")
		rankArr = rankArr[1:]

	outFile.close()

csvToCAMI(sys.argv[1:])