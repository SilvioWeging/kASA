/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*
*  Copyright (C) 2019 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
**************************************************************************/
#include "kASA.hpp"
#include "Read.hpp"
#include "Shrink.hpp"
#include "Compare.hpp"
#include "Update.hpp"

//#include "Unittests.h"

struct SCompareForSortIDsMAP {
	bool operator() (const uint64_t& a, const uint64_t& b) const {
		return a < b;
	}
	static uint64_t min_value() { return numeric_limits<uint64_t>::min(); };
	static uint64_t max_value() { return numeric_limits<uint64_t>::max(); };
};


int main(int argc, char* argv[]) {
	try {
		vector<string> vParameters(argv, argv + argc);

		string cMode = "", sDBPathOut = "", sTempPath = "", sInput = "", contentFileIn = "", readToTaxaFile = "", tableFile = "", indexFile = "", delnodesFile = "", codonTable = "", sTaxonomyPath = "", sAccToTaxFiles = "", sTaxLevel = "", sStxxlMode = "", sCodonID = "1";
		bool bSpaced = false, bVerbose = false, bTranslated = false, bHumanReadable = false, bRAM = false, bUnique = false;
		kASA::Shrink::ShrinkingStrategy eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::TrieHalf;
		int32_t iNumOfThreads = 1, iHigherK = 12, iLowerK = 7, iNumOfCall = 0, iNumOfBeasts = 3;
		uint64_t iMemorySizeAvail = 0;
		float fPercentageOfThrowAway = 0.f;
		uint8_t iTrieDepth = 6;

		auto timeRightNow = chrono::system_clock::to_time_t(chrono::system_clock::now());
		cout << "OUT: " << "kASA version " << kASA_VERSION << " ran on " <<
#if  _WIN64
			"Windows "
#elif __linux__
			"Linux "
#elif __APPLE__
			"macOS "
#else
			"Something other than Win/Linux/macOS "
#endif
			<< "at " << ctime(&timeRightNow) << "OUT: ";
		for (int32_t i = 0; i < argc; ++i) {
			 cout << vParameters[i] << " ";
		}
		cout << endl;

		if (argc == 1) {
			throw runtime_error("No Parameters given!");
		}
		cMode = vParameters[1];

		// still available parameters: jw

		for (int32_t i = 2; i < argc; ++i) {
			string sParameter = vParameters[i];
			if (sParameter == "-o" || sParameter == "--outgoing") {
				sDBPathOut = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-t" || sParameter == "--temp") {
				sTempPath = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-u" || sParameter == "--level") {
				sTaxLevel = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-e" || sParameter == "--unique") {
				bUnique = true;
			}
			else if (sParameter == "-f" || sParameter == "--acc2tax") {
				sAccToTaxFiles = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-y" || sParameter == "--taxonomy") {
				sTaxonomyPath = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-v" || sParameter == "--verbose") {
				bVerbose = true;
			}
			else if (sParameter == "-z" || sParameter == "--translated") {
				bTranslated = true;
			}
			else if (sParameter == "-d" || sParameter == "--database") {
				indexFile = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (cMode != "build") {
					if (!ifstream(indexFile)) {
						throw runtime_error("Index file not found");
					}
					if (!ifstream(indexFile + "_info.txt")) {
						throw runtime_error("Info file not found");
					}
					if (!ifstream(indexFile + "_f.txt") && cMode == "identify") {
						throw runtime_error("Frequency file not found");
					}
					if (!ifstream(indexFile + "_trie") && cMode == "identify") {
						throw runtime_error("Trie file not found");
					}
					if (!ifstream(indexFile + "_trie.txt") && cMode == "identify") {
						throw runtime_error("Trie info file not found");
					}
				}
			}
			else if (sParameter == "-a" || sParameter == "--alphabet") {
				codonTable = Utilities::removeSpaceAndEndline(vParameters[++i]);
				sCodonID = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-b" || sParameter == "--beasts") {
				iNumOfBeasts = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				if (iNumOfBeasts == 0) {
					iNumOfBeasts = 1;
				}
			}
			else if (sParameter == "-r" || sParameter == "--ram") {
				bRAM = true;
			}
			else if (sParameter == "-g" || sParameter == "--percentage") {
				fPercentageOfThrowAway = stof(Utilities::removeSpaceAndEndline(vParameters[++i]));
			}
			else if (sParameter == "-x" || sParameter == "--callidx") {
				iNumOfCall = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
			}
			else if (sParameter == "-n" || sParameter == "--threads") {
				iNumOfThreads = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
			}
			else if (sParameter == "-k") {
				iHigherK = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				iLowerK = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				iHigherK = (iHigherK > 12) ? 12 : iHigherK;
				iLowerK = (iLowerK < 1) ? 1 : iLowerK;
				if (iLowerK > iHigherK) {
					swap(iLowerK, iHigherK);
				}
			}
			else if (sParameter == "--kH") {
				iHigherK = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				iHigherK = (iHigherK > 12) ? 12 : iHigherK;
			}
			else if (sParameter == "--kL") {
				iLowerK = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				iLowerK = (iLowerK < 1) ? 1 : iLowerK;
			}
			else if (sParameter == "-i" || sParameter == "--input") {
				sInput = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (!ifstream(sInput) && sInput.back() != '/') {
					throw runtime_error("Input file not found");
				}
			}
			else if (sParameter == "-q" || sParameter == "--rtt") {
				readToTaxaFile = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-p" || sParameter == "--profile") {
				tableFile = Utilities::removeSpaceAndEndline(vParameters[++i]);
			}
			else if (sParameter == "-m" || sParameter == "--memory") {
				const auto& userGivenMemory = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (userGivenMemory == "inf") {
					iMemorySizeAvail = numeric_limits<uint64_t>::max() / (1024ull * 1024ull);
				}
				else {
					iMemorySizeAvail = 1024ull * stoi(userGivenMemory);
				}
			}
			else if (sParameter == "-s" || sParameter == "--strategy") {
				int32_t iChoice = stoi(Utilities::removeSpaceAndEndline(vParameters[++i]));
				switch (iChoice) {
				case 1:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::EveryNth;
					break;
				case 2:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::TrieHalf;
					break;
				case 3:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::Entropy;
					break;
				case 4:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::Overrepresented;
					break;
				default:
					eShrinkingStrategy = kASA::Shrink::ShrinkingStrategy::TrieHalf;
					break;
				}
			}
			else if (sParameter == "-c" || sParameter == "--content") {
				contentFileIn = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (!ifstream(contentFileIn) && cMode != "build" && cMode != "generateCF") {
					throw runtime_error("Content file not found");
				}
			}
			else if (sParameter == "-l" || sParameter == "--deleted") {
				delnodesFile = Utilities::removeSpaceAndEndline(vParameters[++i]);
				if (!ifstream(delnodesFile)) {
					throw runtime_error("Deleted nodes file not found");
				}
			}
			else if (sParameter == "-h" || sParameter == "--human") {
				bHumanReadable = true;
			}
			else if (sParameter == "--stxxl") {
				sStxxlMode = vParameters[++i];
			}
			else {
				throw runtime_error("Some unknown parameter has been inserted, please check your command line.");
			}
		}
#ifdef ENVIRONMENT32
		iMemorySizeAvail = (iMemorySizeAvail == 0 || iMemorySizeAvail >= numeric_limits<uint32_t>::max() || iMemorySizeAvail >= 2048) ? 1024 : iMemorySizeAvail; // TODO: wenn kleiner als 1024, auf 1024 setzen, sonst bug in stxxl::sort
#else
		iMemorySizeAvail = (iMemorySizeAvail == 0) ? 5120 : iMemorySizeAvail;
#endif
		iMemorySizeAvail *= 1024ull * 1024ull;


		if (cMode == "build") {
			kASA::Read kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, bTranslated, sStxxlMode);
			if (codonTable != "") {
				kASAObj.setCodonTable(codonTable, sCodonID);
			}
			auto start = std::chrono::high_resolution_clock::now();

			// No content file yet created
			if (contentFileIn == "") {
				if (sTaxLevel != "lowest" && (sAccToTaxFiles == "" || sTaxonomyPath == "")) {
					throw runtime_error("No acc2Tax file or taxonomy path given...");
				}
				else {
					contentFileIn = indexFile + "_content.txt";
					if (bVerbose) {
						cout << "OUT: Creating content file: " << contentFileIn << endl;
					}
					kASAObj.generateContentFile(sTaxonomyPath, sAccToTaxFiles, sInput, contentFileIn, sTaxLevel);
				}
			}
			if (iMemorySizeAvail*0.9 < 1024ull * 1024ull * 1024ull) {
				throw runtime_error("Not enough memory given!");
			}
			kASAObj.BuildAll(contentFileIn, sInput, indexFile, static_cast<uint64_t>(iMemorySizeAvail*0.9 - 1024ull * 1024ull * 1024ull), fPercentageOfThrowAway);

			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
		}
		else if (cMode == "generateCF") {
			if (sDBPathOut == "") {
				throw runtime_error("Where should I put the content file?");
			}
			if (sTaxLevel != "lowest" && (sAccToTaxFiles == "" || sTaxonomyPath == "")) {
				throw runtime_error("No acc2Tax file or taxonomy path given...");
			}
			kASA::kASA kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose);
			kASAObj.generateContentFile(sTaxonomyPath, sAccToTaxFiles, sInput, sDBPathOut, sTaxLevel);
		}
		else if (cMode == "update") {
			kASA::Update kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, bTranslated, sStxxlMode);
			if (codonTable != "") {
				kASAObj.setCodonTable(codonTable, sCodonID);
			}

			// Content file was created together with the index
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
			}

			// get taxIds that are new and map the accession numbers to those
			if (sAccToTaxFiles != "" && sTaxonomyPath != "") {
				kASAObj.addToContentFile(sTaxonomyPath, sAccToTaxFiles, sInput, sTaxLevel, contentFileIn);
			}

			auto start = std::chrono::high_resolution_clock::now();
			kASAObj.UpdateFromFasta(contentFileIn, indexFile, sInput, sDBPathOut, (indexFile == sDBPathOut) || (sDBPathOut == ""), static_cast<uint64_t>(iMemorySizeAvail*0.9 - 1024ull * 1024ull * 1024ull), fPercentageOfThrowAway);
			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
		}
		else if (cMode == "delete") {
			kASA::Update kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, bTranslated, sStxxlMode);
			auto start = std::chrono::high_resolution_clock::now();
			kASAObj.DeleteFromLib(indexFile, sDBPathOut, delnodesFile, (indexFile == sDBPathOut), static_cast<uint64_t>(iMemorySizeAvail*0.9 - 1024ull * 1024ull * 1024ull));
			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
		}
		else if (cMode == "shrink") {
			if (indexFile == sDBPathOut) {
				throw runtime_error("Paths and names of input and output are the same!");
			}

			// Default value
			bool bCopyContentFile = false;
			if (sDBPathOut == "") {
				sDBPathOut = indexFile + "_s";
			}

			// Content file was created together with the index
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
				bCopyContentFile = true;
			}

			kASA::Shrink kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);
			kASAObj.ShrinkLib(indexFile, sDBPathOut, eShrinkingStrategy, contentFileIn, fPercentageOfThrowAway);

			if (bCopyContentFile) {
				// Copy file so that both indices can be used by default parameters for -c
				ifstream CFile(contentFileIn, std::ios::binary);
				ofstream SCFile(sDBPathOut + "_content.txt", std::ios::binary);
				SCFile << CFile.rdbuf();
			}
		}
		else if (cMode == "identify") {
			kASA::Compare kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, iNumOfBeasts, bVerbose, bTranslated, sStxxlMode);
			if (codonTable != "") {
				kASAObj.setCodonTable(codonTable, sCodonID);
			}

			// Content file was created together with the index
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
			}

			kASAObj.bHumanReadable = bHumanReadable;
			auto start = std::chrono::high_resolution_clock::now();
			kASAObj.CompareWithLib_partialSort(contentFileIn, indexFile, sInput, readToTaxaFile, tableFile, iTrieDepth, iMemorySizeAvail, bSpaced, bRAM, bUnique);
			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
#if _WIN32 || _WIN64
			//system("PAUSE"); //DEBUG
#endif
		}
		else if (cMode == "getFrequency") {
			kASA::kASA kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
			}

			kASAObj.GetFrequencyK(contentFileIn, indexFile, indexFile + "_f.txt");
		}
		else if (cMode == "redundancy") {
			// Content file was created together with the index
			if (contentFileIn == "") {
				contentFileIn = indexFile + "_content.txt";
			}

			kASA::Shrink kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);
			uint32_t iIdxCounter = 1;
			ifstream content(contentFileIn);
			string sDummy = "";
			while (getline(content, sDummy)) {
				if (sDummy != "") {
					++iIdxCounter;
				}
			}

			ifstream fLibInfo(indexFile + "_info.txt");
			uint64_t iSizeOfLib = 0;
			fLibInfo >> iSizeOfLib;
			fLibInfo.close();

			stxxlFile libFile(indexFile, stxxlFile::RDONLY);
			unique_ptr<const contentVecType_32p> libVec(new const contentVecType_32p(&libFile, iSizeOfLib));
			const auto& iCutoffNumber = kASAObj.histogram(libVec,iIdxCounter);
			if (iCutoffNumber == 1) {
				cout << "OUT: 99% of the k-mers in your index have only one taxon. Using unique frequencies makes sense." << endl;
			}
			else {
				if (iCutoffNumber < 4) {
					cout << "OUT: 99% of the k-mers in your index have " << iCutoffNumber << " or less taxa. Using unique frequencies could make sense." << endl;
				}
				else {
					cout << "OUT: 99% of the k-mers in your index have " << iCutoffNumber << " or less taxa. You should consider looking at the non-unique frequencies as well." << endl;
				}
			}
		}
		else if (cMode == "trie") {
			kASA::kASA kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);

			ifstream fLibInfo(indexFile + "_info.txt");
			uint64_t iSizeOfLib = 0;
			fLibInfo >> iSizeOfLib;
			fLibInfo.close();

			auto start = std::chrono::high_resolution_clock::now();

			stxxlFile libFile(indexFile, stxxlFile::RDONLY);
			const contentVecType_32p libVec(&libFile, iSizeOfLib);

			//cout << stxxl::is_sorted(libVec.cbegin(), libVec.cend(), [](const packedBigPair& a, const packedBigPair& b) {return (a.first < b.first || (!(b.first < a.first) && a.second < b.second)); }) << endl;

			Trie T(static_cast<int8_t>(12), static_cast<int8_t>(iLowerK), iTrieDepth);
			T.SaveToStxxlVec(&libVec, indexFile);

			auto end = std::chrono::high_resolution_clock::now();
			cout << "OUT: Time: " << chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << endl;
		}
		else if (cMode == "half") {
			if (indexFile == sDBPathOut) {
				throw runtime_error("Paths and names of input and output are the same!");
			}
			kASA::Shrink kASAObj(sTempPath, iNumOfThreads, iHigherK, iLowerK, iNumOfCall, bVerbose, sStxxlMode);
			//kASAObj.GetFrequencyK(contentFileIn, databaseFile, sDBPathOut + "_f.txt");
			kASAObj.ShrinkLib(indexFile, sDBPathOut, kASA::Shrink::ShrinkingStrategy::TrieHalf, contentFileIn);
		}
		else if (cMode == "debug") {
			//kASA::Compare UnitTests(sTempPath, 1, 12, 9, 0, iNumOfBeasts);
			//UnitTests.testC2V_1();
			//UnitTests.testC2V_2();
			//UnitTests.testC2V_3();
			//UnitTests.testC2V_4();
			//UnitTests.testC2V_5();
			//UnitTests.testC2V_6();
			//UnitTests.testC2V_7();
			//UnitTests.testC2V_8();
			//UnitTests.testC2V_9();
		}
		else if (cMode == "test") {
			stxxlFile temp(indexFile, stxxl::file::RDONLY);
			ifstream fLibInfo(indexFile + "_info.txt");
			uint64_t iSizeOfLib = 0;
			fLibInfo >> iSizeOfLib;

			ifstream searchFile(vParameters[2]);
			string kMerString = "";
			vector<uint64_t> vSearchVec;
			while (getline(searchFile, kMerString)) {
				vSearchVec.push_back(kASA::kASA::aminoacidTokMer(kMerString));
			}
			uint32_t iSearchCounter = 0;

			const contentVecType_32p tempVec(&temp, iSizeOfLib);
			for (const auto& entry : tempVec) {
				const auto& kMer = entry.first;
				if (kMer < vSearchVec[iSearchCounter]) {
					continue;
				}
				else {
					if (!(vSearchVec[iSearchCounter] < kMer)) {
						cout << kASA::kASA::kMerToAminoacid(kMer, 12) << " " << entry.second << endl;
					}
					else {
						iSearchCounter++;
						if (iSearchCounter >= vSearchVec.size()) {
							break;
						}
						else {
							if (!(vSearchVec[iSearchCounter] < kMer)) {
								cout << kASA::kASA::kMerToAminoacid(kMer, 12) << " " << entry.second << endl;
							}
						}
					}
				}
			}
		}
		else if (cMode == "howmuchtaxids") {
			stxxlFile temp(indexFile, stxxl::file::RDONLY);
			ifstream fLibInfo(indexFile + "_info.txt");
			uint64_t iSizeOfLib = 0;
			fLibInfo >> iSizeOfLib;

			ofstream searchFile(sTempPath+"frequentkMers.txt");
			uint64_t iSeenkMer = 0;
			uint32_t iSearchCounter = 0;
			uint32_t iMagicNumber = 4;
			set<uint32_t> taxIDs;

			const contentVecType_32p tempVec(&temp, iSizeOfLib);
			for (const auto& entry : tempVec) {
				const auto& kMer = entry.first;
				if (kMer == iSeenkMer) {
					taxIDs.insert(entry.second);
					iSearchCounter++;
				}
				else {
					if (iSearchCounter >= iMagicNumber) {
						searchFile << kASA::kASA::kMerToAminoacid(kMer, 12);
						for (const auto& elem : taxIDs) {
							searchFile << " " << elem;
						}
						searchFile << endl;
					}
					iSearchCounter = 0;
					taxIDs.clear();
					iSeenkMer = kMer;
				}
			}
		}
		else if (cMode == "showVec") {
			stxxlFile temp(indexFile, stxxl::file::RDONLY);
			ifstream sizeFile(indexFile + "_info.txt");
			uint64_t iSize = 0;
			sizeFile >> iSize;
			bool index_t = false;
			sizeFile >> index_t;
			if (index_t) {
				kASA::kASA::showVec(index_t_p(&temp, iSize));
			}
			else {
				kASA::kASA::showVec(contentVecType_32p(&temp, iSize));
			}
		}
		else if (cMode == "transform") {
			stxxlFile temp(indexFile, stxxl::file::RDONLY);

			ifstream sizeFile(indexFile + "_info.txt"); // "_info.txt"
			uint64_t iSize = 0;
			sizeFile >> iSize;

			const contentVecType_32p tempV(&temp, iSize);

			Utilities::createFile(sDBPathOut);
			Utilities::createFile(sDBPathOut + "_2");

			stxxlFile tempOut(sDBPathOut, stxxl::file::RDWR);
			stxxl::VECTOR_GENERATOR<uint64_t, 4U, 4U, 2101248, stxxl::RC>::result tempPV(&tempOut, iSize);

			stxxlFile tempOut2(sDBPathOut+"_2", stxxl::file::RDWR);
			stxxl::VECTOR_GENERATOR<uint32_t, 4U, 4U, 2101248, stxxl::RC>::result tempPV2(&tempOut2, iSize);
			
			ofstream countsFile(sDBPathOut+"_counts.txt");

			auto t1It = tempV.cbegin();
			auto tOIt = tempPV.begin();
			auto tO2It = tempPV2.begin();
			uint64_t iCount = 0;
			uint64_t iSeen = 0;
			for (; t1It != tempV.cend();) {
				if (t1It->first == iSeen) {
					*tO2It = t1It->second;
					++t1It;
				}
				else {
					*tOIt = iSeen = t1It->first;
					countsFile << iCount << endl;
					*tO2It = t1It->second;
					++t1It; ++tOIt;
				}
				++tO2It;
				++iCount;
			}

			tempPV.resize(tOIt - tempPV.begin(), true);

			ofstream outSizeFile(sDBPathOut + "_info.txt");
			outSizeFile << tempPV.size() << endl << iSize;
			tempPV.export_files("_");
			tempPV2.export_files("_");
		}
	} catch (const exception& e) {
		cerr << "ERROR: " << e.what() << endl;
		return 1;
	}

#if (_WIN32 || _WIN64) && _DEBUG
	system("PAUSE");
#endif
	return 0;
}