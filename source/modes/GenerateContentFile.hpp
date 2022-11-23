/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*
*  Copyright (C) 2020 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
**************************************************************************/
#pragma once

#include "../kASA.hpp"

namespace kASA {

	class ContentFile : public kASA {
	public:
		ContentFile(const InputParameters& cParams) : kASA(cParams) {}
		ContentFile(const kASA& obj) : kASA(obj) {}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Generate temporary content file(s) in case not enough memory is available to hold the unordered_maps
		inline void generateTemporaryContentFile(const string& sTempCFilePath, const string& sTaxonomyPath, string sTaxonomicLevel, const bool& bTaxIdsAsStrings, const string& sAccToTaxFiles, pair<uint32_t, uint32_t> poolAndNames, unordered_map<string, bool>& vAccessions, unordered_map<string, uint32_t>& vEntriesWithoutAccNr, unordered_map<string, string>& vNamesFromFasta) {
			if (_bVerbose) {
				cout << "OUT: Generating temporary content file: " << sTempCFilePath << endl;
			}
			
			debugBarrier
			ofstream contentFile = Utilities::createFileAndGetIt(sTempCFilePath);

			size_t iIdentifiedCounter = 0;
			unordered_map<string, unordered_set<string>> taxWithAccNrs;
			unordered_map<string, string> taxToNames;
			bool bNotAllFound = true;

			if (_bVerbose) {
				cout << "OUT: Assigning taxids... " << endl;
			}

			if (sTaxonomicLevel == "lowest") {
				iIdentifiedCounter = 1;
				for (auto& entry : vAccessions) {
					auto elem = taxWithAccNrs.insert(make_pair(to_string(iIdentifiedCounter), unordered_set<string>()));
					if (elem.second) {
						elem.first->second.insert(entry.first);
						entry.second = true;

						taxToNames.insert(make_pair(to_string(iIdentifiedCounter), vNamesFromFasta.at(entry.first)));

						++iIdentifiedCounter;
					}
					else {
						throw runtime_error("AccTaxPair could not be added");
					}
				}
				bNotAllFound = false;
				debugBarrier
			}
			else {
				auto files = Utilities::gatherFilesFromPath(sAccToTaxFiles).first;

				////////////////////

				for (const auto& file : files) {
					string sDummy = "";
					size_t iIdxForAccession = 1, iIdxForTaxID = 2; // true = 4 columns, false = 2 columns

					bool isGzipped = file.second;
					ifstream acc2Tax;
					igzstream acc2TaxGZ;
					if (isGzipped) {
						acc2TaxGZ.open(file.first.c_str());
						getline(acc2TaxGZ, sDummy);
						const auto& columns = Utilities::split(sDummy, '\t');
						if (columns.size() == 2) {
							iIdxForAccession = 0;
							iIdxForTaxID = 1;
						}
						acc2TaxGZ.close();
						acc2TaxGZ.open(file.first.c_str());
					}
					else {
						acc2Tax.open(file.first);
						getline(acc2Tax, sDummy);
						const auto& columns = Utilities::split(sDummy, '\t');
						if (columns.size() == 2) {
							iIdxForAccession = 0;
							iIdxForTaxID = 1;
						}
						acc2Tax.seekg(0);
					}


					debugBarrier

					while (((isGzipped) ? getline(acc2TaxGZ, sDummy) : getline(acc2Tax, sDummy)) && bNotAllFound) {
						const auto& columns = Utilities::split(sDummy, '\t');
						auto res1 = vAccessions.find(columns[iIdxForAccession]);
						if (res1 != vAccessions.end()) {
							auto res2 = taxWithAccNrs.find(columns[iIdxForTaxID]);
							if (res2 != taxWithAccNrs.end()) {
								res2->second.insert(res1->first);
							}
							else {
								auto res3 = taxWithAccNrs.insert(make_pair(columns[iIdxForTaxID], unordered_set<string>()));
								if (res3.second) {
									res3.first->second.insert(res1->first);
								}
								else {
									throw runtime_error("AccTaxPair could not be added");
								}
							}
							res1->second = true;
							++iIdentifiedCounter;
							if (iIdentifiedCounter == vAccessions.size()) {
								bNotAllFound = false;
							}
						}
					}
					debugBarrier
				}
			}
			vNamesFromFasta.clear();

			debugBarrier
			////////////////////
			if (_bVerbose && bNotAllFound) {
				cout << "OUT: The accession numbers without a taxid will be written to " + _sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_accessionsWithoutTaxid.txt" << endl;
				ofstream noIDFile(_sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_accessionsWithoutTaxid.txt", ios::app);
				for (const auto& entry : vAccessions) {
					if (entry.second == false) {
						vEntriesWithoutAccNr.insert(make_pair(entry.first, 0));
						if (_bVerbose) {
							noIDFile << entry.first << endl;
						}
					}
				}
			}
			else {
				for (const auto& entry : vAccessions) {
					if (entry.second == false) {
						vEntriesWithoutAccNr.insert(make_pair(entry.first, 0));
					}
				}
			}
			vAccessions.clear();
			debugBarrier

			////////////////////
			if (_bVerbose) {
				cout << "OUT: Providing dummy taxid(s) if needed..." << endl;
			}

			uint32_t pool = poolAndNames.first;
			if (vEntriesWithoutAccNr.size() > pool) {
				throw runtime_error("Too many dummys, I can't handle this!");
			}
			for (auto& entry : vEntriesWithoutAccNr) {
				entry.second = pool;
				--pool;
			}

			debugBarrier
			////////////////////
			if (_bVerbose) {
				cout << "OUT: Fetching name(s)..." << endl;
			}

			//taxToNames
			if (taxToNames.empty()) {
				ifstream names(sTaxonomyPath + "names.dmp");
				string sDummy = "";
				while (getline(names, sDummy)) {
					const auto& line = Utilities::split(sDummy, '|');
					if (line[3] == "\tscientific name\t") {
						taxToNames.insert(make_pair(Utilities::rstrip(line[0]), Utilities::lstrip(Utilities::rstrip(line[1]))));
					}
				}
			}

			debugBarrier
			////////////////////
			if (_bVerbose) {
				cout << "OUT: Creating dictionary to map taxid(s) to your specified tax. level..." << endl;
			}

			unordered_map<string, pair<string, string>> taxToNodes;
			{
				ifstream nodes(sTaxonomyPath + "nodes.dmp");
				string sDummy = "";
				while (getline(nodes, sDummy)) {
					const auto& line = Utilities::split(sDummy, '|');
					taxToNodes.insert(make_pair(Utilities::rstrip(line[0]), make_pair(Utilities::lstrip(Utilities::rstrip(line[1])), Utilities::lstrip(Utilities::rstrip(line[2])))));
				}
			}

			debugBarrier
			////////////////////
			if (_bVerbose) {
				cout << "OUT: Linking to your specified tax. level..." << endl;
			}

			auto sortingFuncToContentFile = [&bTaxIdsAsStrings](const string& a, const string& b) { if (bTaxIdsAsStrings) { return a < b; } else { return stoull(a) < stoull(b); } };
			map<string, pair<unordered_set<string>, unordered_set<string>>, decltype(sortingFuncToContentFile)> taxToTaxWAccs(sortingFuncToContentFile);
			{
				std::transform(sTaxonomicLevel.begin(), sTaxonomicLevel.end(), sTaxonomicLevel.begin(), ::tolower);
				if (!(sTaxonomicLevel == "lowest" ||
					sTaxonomicLevel == "subspecies" ||
					sTaxonomicLevel == "species" ||
					sTaxonomicLevel == "genus" ||
					sTaxonomicLevel == "family" ||
					sTaxonomicLevel == "order" ||
					sTaxonomicLevel == "class" ||
					sTaxonomicLevel == "phylum" ||
					sTaxonomicLevel == "kingdom" ||
					sTaxonomicLevel == "superkingdom" ||
					sTaxonomicLevel == "domain"
					)) {
					cerr << "WARNING: No known tax. level specified. I'll just go with species..." << endl;
					sTaxonomicLevel = "species";
				}

				string upperTax = "";
				pair<string, string> entry;
				for (const auto& elem : taxWithAccNrs) {
					upperTax = elem.first;
					if (sTaxonomicLevel != "lowest") {
						auto res1 = taxToNodes.find(upperTax);
						if (res1 != taxToNodes.end()) {
							entry = res1->second;
						}
						else {
							entry = make_pair("1", "");
						}
						while (entry.second != sTaxonomicLevel && entry.first != "1") {
							upperTax = entry.first;
							entry = taxToNodes.at(upperTax);
						}
						if (entry.first == "1") {
							upperTax = elem.first;
						}


					}
					auto res2 = taxToTaxWAccs.find(upperTax);
					if (res2 != taxToTaxWAccs.end()) {
						res2->second.first.insert(elem.first);
						res2->second.second.insert(elem.second.begin(), elem.second.end());
					}
					else {
						auto it = taxToTaxWAccs.insert(make_pair(upperTax, make_pair(unordered_set<string>(), unordered_set<string>())));
						it.first->second.first.insert(elem.first);
						it.first->second.second.insert(elem.second.begin(), elem.second.end());
					}
				}
				taxWithAccNrs.clear();
			}

			debugBarrier
			////////////////////

			if (_bVerbose) {
				cout << "OUT: Writing to temporary content file..." << endl;
			}

			uint32_t iUnnamedCounter = 0;
			uint64_t iLineCounter = 1;
			for (const auto& elem : taxToTaxWAccs) {
				string taxa = "", accnrs = "";
				for (const auto& tax : elem.second.first) {
					taxa += tax + ";";
				}
				if (taxa != "") {
					taxa.erase(taxa.end() - 1);
				}
				for (const auto& accs : elem.second.second) {
					accnrs += accs + ";";
				}
				if (accnrs != "") {
					accnrs.erase(accnrs.end() - 1);
				}

				const auto res = taxToNames.find(elem.first);
				if (res != taxToNames.end()) {
					contentFile << Utilities::removeCharFromString(res->second, ',') << "\t" << elem.first << "\t" << taxa << "\t" << accnrs << ((bTaxIdsAsStrings) ? ("\t" + to_string(iLineCounter++)) : "") << "\n";
				}
				else {
					contentFile << "unnamed_" << iUnnamedCounter++ << "\t" << elem.first << "\t" << taxa << "\t" << accnrs << ((bTaxIdsAsStrings) ? ("\t" + to_string(iLineCounter++)) : "") << "\n";
				}
			}

			iUnnamedCounter = poolAndNames.second;
			for (const auto& elem : vEntriesWithoutAccNr) {
				contentFile << "EWAN_" << iUnnamedCounter++ << "\t" << elem.second << "\t" << elem.second << "\t" << elem.first << ((bTaxIdsAsStrings) ? ("\t" + to_string(iLineCounter++)) : "") << "\n";
			}

			vEntriesWithoutAccNr.clear();

			if (_bVerbose) {
				cout << "OUT: Finished generating temporary content file, either resuming with gathering accessions or writing to final content file... " << endl;
			}
			debugBarrier
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// C++ version of the python script that creates the content file
		inline void generateContentFile(const InputParameters& cParams, const string& sContentFile, pair<uint32_t, uint32_t> poolAndNames = make_pair(numeric_limits<uint32_t>::max() - 1, 0UL)) {
			try {
				if (cParams.sTaxLevel != "lowest") {
					if (!ifstream(cParams.sTaxonomyPath + "names.dmp") || !ifstream(cParams.sTaxonomyPath + "nodes.dmp")) {
						throw runtime_error("The taxonomy files couldn't be found");
					}
				}
				debugBarrier
				Utilities::checkIfFileCanBeCreated(sContentFile);

				auto files = Utilities::gatherFilesFromPath(cParams.sInput).first;

				if (_bVerbose) {
					cout << "OUT: Going through fasta(s), gathering accession number(s)..." << endl;
				}

				////////////////////
				unordered_map<string, bool> vAccessions;
				unordered_map<string, uint32_t> vEntriesWithoutAccNr;
				unordered_map<string, string> vNamesFromFasta;
				size_t iNumOfLegitAccessions = 0, iNumOfDummys = 0;

				uint32_t iTemporaryCounter = 0;
				uint64_t iMemoryAllocated = sizeof(unordered_map<string, bool>)
					+ sizeof(unordered_map<string, uint32_t>)
					+ 2 * sizeof(unordered_map<string, string>)
					+ sizeof(unordered_map<string, unordered_set<string>>); // all unordered_maps that are memory critical

				debugBarrier
				for (const auto& file : files) {
					ifstream fastaFile;
					igzstream fastaFileGZ;

					if (file.second) {
						fastaFileGZ.open(file.first.c_str());
						if (_bVerbose) {
							cout << "OUT: Current file: " << file.first.c_str() << endl;
						}
					}
					else {
						fastaFile.open(file.first);
						if (_bVerbose) {
							cout << "OUT: Current file: " << file.first << endl;
						}
					}

					debugBarrier
					string sDummy = "";
					while ((file.second) ? getline(fastaFileGZ, sDummy) : getline(fastaFile, sDummy)) {
						if (sDummy != "") {
							if (sDummy.front() == '>') {
								sDummy.erase(sDummy.begin());
								const auto& sNumbers = Utilities::split(Utilities::split(sDummy, ' ').at(0), '|');
								string acc = "";
								for (const auto& entry : sNumbers) {
									if (entry.find('.') != string::npos) {
										acc = entry;
										break;
									}
								}
								if (acc != "") {
									vAccessions.insert(make_pair(acc, false));
									iMemoryAllocated += sizeof(pair<string, bool>) + sizeof(char) * acc.length();
									if (cParams.sTaxLevel == "lowest") {
										vNamesFromFasta.insert(make_pair(acc, Utilities::replaceCharacter(sDummy, ',', ' ')));
										iMemoryAllocated += sizeof(pair<string, string>) + sizeof(char) * 2 * acc.length();
									}
								}
								else {
									vEntriesWithoutAccNr.insert(make_pair(sDummy, 0));
									iMemoryAllocated += sizeof(pair<string, uint32_t>) + sizeof(char) * sDummy.length();
								}
								iMemoryAllocated += sizeof(pair<string, unordered_set<string>>) + sizeof(char) * acc.length() + sizeof(char) * 9
									+ sizeof(pair<string, string>) + sizeof(char) * 30; // numbers are rough approximations of maximum string lengths for iIdentifiedCounter and length of species name
							}
							if (iMemoryAllocated > static_cast<uint64_t>(cParams.iMemorySizeAvail)) {
								iNumOfLegitAccessions += vAccessions.size();
								iNumOfDummys += vEntriesWithoutAccNr.size();
								debugBarrier
								generateTemporaryContentFile(_sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_" + Utilities::itostr(iTemporaryCounter) + ".txt", cParams.sTaxonomyPath, cParams.sTaxLevel, cParams.bTaxIdsAsStrings, cParams.sAccToTaxFiles, poolAndNames, vAccessions, vEntriesWithoutAccNr, vNamesFromFasta);
								++iTemporaryCounter;
								iMemoryAllocated = 0;
								debugBarrier
							}
						}
					}
				}
				debugBarrier

				iNumOfLegitAccessions += vAccessions.size();
				iNumOfDummys += vEntriesWithoutAccNr.size();

				if (iTemporaryCounter == 0) {
					if (_bVerbose) {
						cout << "OUT: " << iNumOfLegitAccessions << " legit accession number(s) found in " << files.size() << " file(s)." << endl
							<< "OUT: " << iNumOfDummys << " entr(y/ies) will get a dummy ID." << endl
							<< "OUT: Creating content file..." << endl;
					}
					debugBarrier
					generateTemporaryContentFile(sContentFile, cParams.sTaxonomyPath, cParams.sTaxLevel, cParams.bTaxIdsAsStrings, cParams.sAccToTaxFiles, poolAndNames, vAccessions, vEntriesWithoutAccNr, vNamesFromFasta);
				}
				else {
					if (vAccessions.size()) {
						debugBarrier
						generateTemporaryContentFile(_sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_" + Utilities::itostr(iTemporaryCounter) + ".txt", cParams.sTaxonomyPath, cParams.sTaxLevel, cParams.bTaxIdsAsStrings, cParams.sAccToTaxFiles, poolAndNames, vAccessions, vEntriesWithoutAccNr, vNamesFromFasta);
						++iTemporaryCounter;
					}

					if (_bVerbose) {
						cout << "OUT: " << iNumOfLegitAccessions << " legit accession number(s) found in " << files.size() << " file(s)." << endl
							<< "OUT: " << iNumOfDummys << " entr(y/ies) got a dummy ID." << endl
							<< "OUT: Merging temporary content files..." << endl;
					}

					debugBarrier
					string sTempFinal1 = _sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_final1.txt";
					string sTempFinal2 = _sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_final2.txt";
					mergeContentFiles(_sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_0.txt", _sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_1.txt", false, sTempFinal1);
					remove((_sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_0.txt").c_str());
					remove((_sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_1.txt").c_str());
					debugBarrier
					for (uint32_t iTemporaryCFilesIdx = 2; iTemporaryCFilesIdx < iTemporaryCounter; ++iTemporaryCFilesIdx) {
						string currentFile = _sTemporaryPath + "content_" + Utilities::itostr(_iNumOfCall) + "_" + Utilities::itostr(iTemporaryCFilesIdx) + ".txt";
						mergeContentFiles(sTempFinal1, currentFile, false, sTempFinal2);
						swap(sTempFinal1, sTempFinal2);
						remove(currentFile.c_str());
					}
					debugBarrier
					Utilities::moveFile(sTempFinal1, sContentFile);
					remove(sTempFinal2.c_str());

				}
				debugBarrier

			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Merge two content files
		inline pair<unordered_map<uint32_t, uint32_t>, unordered_map<uint32_t, uint32_t>> mergeContentFiles(const string& contentFile, const string& contentFile2, const bool& bMergeExistingIndices, const string& contentOut) {
			try {
				auto mergeEntries = [](const string& sTempLineSpecIDs, const string& sCurrLineSpecIDs, const string& sTempLineAccNrs, const string& sCurrLineAccNrs) {
					unordered_set<string> sSpecIDs, sAccNrs;
					// concatenate non-redundant species IDs
					const auto& newSpecIDs = Utilities::split(sTempLineSpecIDs, ';'); // tempLineContent[2]
					const auto& currentSpecIDs = Utilities::split(sCurrLineSpecIDs, ';'); // get<1>(entry->second)
					sSpecIDs.insert(newSpecIDs.cbegin(), newSpecIDs.cend());
					sSpecIDs.insert(currentSpecIDs.cbegin(), currentSpecIDs.cend());
					string sNewSpecIDs = "";
					for (const auto& elem : sSpecIDs) {
						sNewSpecIDs += elem + ";";
					}
					sNewSpecIDs.pop_back();

					// concatenate non-redundant accession numbers
					const auto& newAccNrs = Utilities::split(sTempLineAccNrs, ';'); // tempLineContent[3]
					const auto& currentAccNrs = Utilities::split(sCurrLineAccNrs, ';'); // get<2>(entry->second)
					sAccNrs.insert(newAccNrs.cbegin(), newAccNrs.cend());
					sAccNrs.insert(currentAccNrs.cbegin(), currentAccNrs.cend());
					string sNewAccNrs = "";
					for (const auto& elem : sAccNrs) {
						sNewAccNrs += elem + ";";
					}
					sNewAccNrs.pop_back();

					return make_pair(sNewSpecIDs, sNewAccNrs);
				};
				debugBarrier
				pair<uint32_t, uint32_t> pool(numeric_limits<uint32_t>::max(), 0); // first component counts downwards from the highest possible number to distinguish it from valid taxIDs, the second component gives the entry a name ID
				vector<string> vListOfDummys;
				unordered_map<uint32_t, uint32_t> mapOfOldDummyToNewDummy1;
				unordered_map<uint32_t, uint32_t> mapOfOldDummyToNewDummy2;

				ifstream fContent, fContent2;
				ofstream fContentOut;
				bool bTaxIdsAsStrings = false;
				//fContent.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
				fContent.open(contentFile);
				if (!fContent) {
					throw runtime_error("First content file couldn't be read!");
				}
				fContent2.open(contentFile2);
				if (!fContent2) {
					throw runtime_error("Second content file couldn't be read!");
				}
				fContentOut.open(contentOut);
				if (!fContentOut) {
					throw runtime_error("Resulting content file couldn't be opened for writing!");
				}

				uint64_t iLargestLineIndex = 1;

				// The following algorithm assumes that the content files are sorted regarding the taxID.
				string currLine1 = "", currLine2 = "";
				getline(fContent, currLine1);
				getline(fContent2, currLine2);
				if (currLine1 == "" || currLine2 == "") {
					cerr << "ERROR: One of the two content files has a leading empty line. This should not happen. Please look into both files and confirm that they are valid." << endl;
					throw runtime_error("Invalid content files!");
				}
				if ((Utilities::split(currLine1, '\t')).size() >= 5 && !bTaxIdsAsStrings) {
					bTaxIdsAsStrings = true;
				}
				auto compareFunc = [&bTaxIdsAsStrings](const string& a, const string& b) { if (bTaxIdsAsStrings) { return a < b; } else { return stoull(a) < stoull(b); } };

				debugBarrier

				while (fContent && fContent2) {
					const auto& currLine1Split = Utilities::split(currLine1, '\t');
					const auto& currLine2Split = Utilities::split(currLine2, '\t');

					if (currLine1Split[0].find("EWAN") != string::npos) {
						if (bMergeExistingIndices) {
							const auto& ewanTaxID = stoul(currLine1Split[1]);
							mapOfOldDummyToNewDummy1.insert(make_pair(ewanTaxID,pool.first)); // get current dummy taxID and map to new one
							pool.first--;
						}
						vListOfDummys.push_back(currLine1Split[3]);
						getline(fContent, currLine1);
						continue;
					}
					if (currLine2Split[0].find("EWAN") != string::npos) {
						if (bMergeExistingIndices) {
							const auto& ewanTaxID = stoul(currLine2Split[1]);
							mapOfOldDummyToNewDummy2.insert(make_pair(ewanTaxID, pool.first)); // get current dummy taxID and map to new one
							pool.first--;
						}
						vListOfDummys.push_back(currLine2Split[3]);
						getline(fContent2, currLine2);
						continue;
					}

					// <
					if (compareFunc(currLine1Split[1], currLine2Split[1])) {
						fContentOut << currLine1Split[0] << "\t" << currLine1Split[1] << "\t" << currLine1Split[2] << "\t" << currLine1Split[3] << ((bTaxIdsAsStrings) ? ("\t" + Utilities::itostr(iLargestLineIndex++)) : "") << "\n";
						getline(fContent, currLine1);
					}
					else {
						// ==
						if (!compareFunc(currLine2Split[1], currLine1Split[1])) {
							const auto& res = mergeEntries(currLine1Split[2], currLine2Split[2], currLine1Split[3], currLine2Split[3]);
							fContentOut << currLine2Split[0] << "\t" << currLine2Split[1] << "\t" << res.first << "\t" << res.second << ((bTaxIdsAsStrings) ? ("\t" + Utilities::itostr(iLargestLineIndex++)) : "") << "\n";
							getline(fContent, currLine1);
							getline(fContent2, currLine2);
						}
						// >
						else {
							fContentOut << currLine2Split[0] << "\t" << currLine2Split[1] << "\t" << currLine2Split[2] << "\t" << currLine2Split[3] << ((bTaxIdsAsStrings) ? ("\t" + Utilities::itostr(iLargestLineIndex++)) : "") << "\n";
							getline(fContent2, currLine2);
						}
					}
				}

				debugBarrier
				while (fContent) {
					const auto& currLine1Split = Utilities::split(currLine1, '\t');

					if (currLine1Split[0].find("EWAN") != string::npos) {
						if (bMergeExistingIndices) {
							const auto& ewanTaxID = stoul(currLine1Split[1]);
							mapOfOldDummyToNewDummy1.insert(make_pair(ewanTaxID, pool.first)); // get current dummy taxID and map to new one
							pool.first--;
						}
						vListOfDummys.push_back(currLine1Split[3]);
					}
					else {
						fContentOut << currLine1Split[0] << "\t" << currLine1Split[1] << "\t" << currLine1Split[2] << "\t" << currLine1Split[3] << ((bTaxIdsAsStrings) ? ("\t" + Utilities::itostr(iLargestLineIndex++)) : "") << "\n";
					}
					getline(fContent, currLine1);
				}
				debugBarrier
				while (fContent2) {
					const auto& currLine2Split = Utilities::split(currLine2, '\t');

					if (currLine2Split[0].find("EWAN") != string::npos) {
						if (bMergeExistingIndices) {
							const auto& ewanTaxID = stoul(currLine2Split[1]);
							mapOfOldDummyToNewDummy2.insert(make_pair(ewanTaxID, pool.first)); // get current dummy taxID and map to new one
							pool.first--;
						}
						vListOfDummys.push_back(currLine2Split[3]);
					}
					else {
						fContentOut << currLine2Split[0] << "\t" << currLine2Split[1] << "\t" << currLine2Split[2] << "\t" << currLine2Split[3] << ((bTaxIdsAsStrings) ? ("\t" + Utilities::itostr(iLargestLineIndex++)) : "") << "\n";
					}
					getline(fContent2, currLine2);
				}

				debugBarrier
				uint32_t iDummyID = numeric_limits<uint32_t>::max();
				for (const auto& entry : vListOfDummys) {
					const auto& IDStr = to_string(iDummyID--);
					const auto& nameStr = "EWAN_" + to_string(pool.second++);
					fContentOut << nameStr << "\t" << IDStr << "\t" << IDStr << "\t" << entry << ((bTaxIdsAsStrings) ? ("\t" + Utilities::itostr(iLargestLineIndex++)) : "") << "\n";
				}

				return make_pair(mapOfOldDummyToNewDummy1, mapOfOldDummyToNewDummy2);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Add new stuff to an existing content file
		inline pair<unordered_map<uint32_t, uint32_t>, unordered_map<uint32_t, uint32_t>> addToContentFile(const InputParameters& cParams, const string& sContentFileOut) {
			// Warning: taxonomic levels must not change in an update
			const string& temporaryContentFile = _sTemporaryPath + Utilities::itostr(_iNumOfCall) + "_tempContent.txt";
			const string& temporaryContentOutFile = _sTemporaryPath + Utilities::itostr(_iNumOfCall) + "_tempContentOut.txt";

			if (_bVerbose) {
				cout << "OUT: Generating temporary content file from new data..." << endl;
			}
			debugBarrier
			generateContentFile(cParams, temporaryContentFile);

			if (_bVerbose) {
				cout << "OUT: Merging content file from new data with the existing one..." << endl;
			}
			debugBarrier
			const auto& result = mergeContentFiles(cParams.contentFileIn, temporaryContentFile, true, temporaryContentOutFile);
			remove(temporaryContentFile.c_str());
			remove(sContentFileOut.c_str());
			Utilities::moveFile(temporaryContentOutFile, sContentFileOut);

			return result;
		}
	};
}