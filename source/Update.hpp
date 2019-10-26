/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*
*  Copyright (C) 2019 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
**************************************************************************/
#pragma once

#include "Read.hpp"
#include "Build.hpp"

namespace kASA {
	class Update : public Read, Build {
	public:
		Update(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const bool& bTranslated = false, const string& stxxl_mode = "") : Read(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, bTranslated, stxxl_mode), Build() {}

	public:


		void DeleteFromLib(const string& sLibFile, const string& fOutfile, const string& sDelNodesPath, const bool& bOverwrite, const uint64_t& iMemory) {
			try {
				// test if files exists
				if (!ifstream(sLibFile)) {
					throw runtime_error("Library file does not exist");
				}

				if (!ifstream(sDelNodesPath)) {
					throw runtime_error("delnodes.dmp not found");
				}

				ifstream deletedNodes(sDelNodesPath);
				string line = "";
				vector<uint32_t> vTaxIdsTBD;
				while (getline(deletedNodes, line)) {
					if (line != "") {
						vTaxIdsTBD.push_back(stoul((Utilities::split(line, '\t'))[0]));
					}
				}

				if (bOverwrite) {
					uint64_t iCounterOfHits = 0;
					ifstream infoFile(sLibFile + "_info.txt");
					uint64_t iLibSize;
					infoFile >> iLibSize;
					infoFile.close();
					stxxlFile libFile(sLibFile, stxxl::file::RDWR);
					contentVecType_32p vLib(&libFile, iLibSize);

					// Get an end iterator, that doesn't contain a value that should be deleted
					auto itEnd = vLib.end() - 1;
					while (find(vTaxIdsTBD.begin(), vTaxIdsTBD.end(), itEnd->second) != vTaxIdsTBD.end()) {
						--itEnd;
						++iCounterOfHits;
					}
					// iterate from front towards the end iterator, switching to-be-deleted values to the end and then cut them off later as well as changing any taxIds that need changing
					for (auto itFront = vLib.begin(); itFront != itEnd;) {
						//deletion-switching
						if (find(vTaxIdsTBD.begin(), vTaxIdsTBD.end(), itFront->second) != vTaxIdsTBD.end()) {
							++iCounterOfHits;
							const auto tempVal = *itEnd;
							*itEnd = *itFront;
							*itFront = tempVal;
							--itEnd;
							while (find(vTaxIdsTBD.begin(), vTaxIdsTBD.end(), itEnd->second) != vTaxIdsTBD.end() && itEnd != itFront) {
								--itEnd;
								++iCounterOfHits;
							}
							if (itEnd == itFront) {
								break;
							}
						}

						++itFront;
					}
					vLib.resize(vLib.size() - iCounterOfHits, true);

					stxxl::sort(vLib.begin(), vLib.end(), SCompareStructForSTXXLSort(), iMemory);
					vLib.export_files("_");

					ofstream infoFileOut(sLibFile + "_info.txt");
					infoFileOut << vLib.size();

					Trie T(static_cast<int8_t>(12), static_cast<int8_t>(_iMinK), 6);
					T.SaveToStxxlVec(&vLib, sLibFile);
				}
				else {
					ifstream infoFile(sLibFile + "_info.txt");
					uint64_t iLibSize;
					infoFile >> iLibSize;
					infoFile.close();
					stxxlFile libFile(sLibFile, stxxl::file::RDONLY);
					const contentVecType_32p vLib(&libFile, iLibSize);
					ofstream dummy;
					//dummy.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
					dummy.open(fOutfile);
					dummy.close();
					stxxlFile outFile(fOutfile, stxxl::file::RDWR);
					contentVecType_32p vOut(&outFile, iLibSize);

					// All entries with valid taxonomic ids will be written to the other library
					auto itOut = vOut.begin();
					for (auto itLib = vLib.cbegin(); itLib != vLib.cend(); ++itLib) {
						if (find(vTaxIdsTBD.begin(), vTaxIdsTBD.end(), itLib->second) == vTaxIdsTBD.end()) {
							*itOut++ = *itLib;
						}
					}
					vOut.resize(itOut - vOut.begin(), true);
					vOut.export_files("_");

					ofstream infoFileOut(fOutfile + "_info.txt");
					infoFileOut << vOut.size();

					Trie T(static_cast<int8_t>(12), static_cast<int8_t>(_iMinK), 6);
					T.SaveToStxxlVec(&vOut, fOutfile);
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Update an existing library with a fasta file
		void UpdateFromFasta(const string& contentFile, const string& sLibFile, const string& sDirectory, const string& fOutFile, const bool& bOverwrite, const uint64_t& iMemory, const float& fPercentageOfThrowAway) {
			try {
				// test if files exists
				if (!ifstream(contentFile) || !ifstream(sLibFile)) {
					throw runtime_error("One of the files does not exist");
				}

				// read new content file
				uint32_t iIdxCounter = 1;
				unordered_map<uint32_t, uint32_t> mIDsAsIdx; mIDsAsIdx[0] = 0;
				unordered_map<uint32_t, string> mIdxToName; mIdxToName[0] = "non_unique";
				unordered_map<string, uint32_t> mAccToID;
				ifstream content(contentFile);
				string sDummy = "";
				while (getline(content, sDummy)) {
					if (sDummy != "") {
						const auto& line = Utilities::split(sDummy, '\t');
						if (line.size() == 4) {
							mIDsAsIdx[stoul(line[1])] = iIdxCounter;
							mIdxToName[iIdxCounter] = line[0];
							++iIdxCounter;
							const auto& vAccessionNumbers = Utilities::split(line[3], ';');
							for (const auto& acc : vAccessionNumbers) {
								mAccToID.insert(make_pair(acc, stoul(line[1])));
							}
						}
					}
				}

				unique_ptr<uint64_t[]> arrFrequencies_(new uint64_t[iIdxCounter * 12]);
				for (uint64_t i = 0; i < iIdxCounter * 12; ++i) {
					arrFrequencies_[i] = 0;
				}

				/*uint32_t iEntryCounter = 0;
				ifstream freqFile(sLibFile + "_f.txt");
				while (getline(freqFile, sDummy)) {
					if (sDummy != "") {
						const auto& line = Utilities::split(sDummy, '\t');
						for (uint8_t k = 0; k < 12; ++k) {
							arrFrequencies[iEntryCounter * 12 + k] = stoull(line[k + 1]);
						}
						++iEntryCounter;
					}
				}
				freqFile.close();*/

				if (bOverwrite) {
					fstream fLibInfo(sLibFile + "_info.txt", ios::in);
					uint64_t iSizeOfLib = 0;
					fLibInfo >> iSizeOfLib;
					fLibInfo.close();
					unique_ptr<stxxlFile> stxxlLibFile(new stxxlFile(sLibFile, stxxl::file::RDWR));
					unique_ptr<contentVecType_32p> vLibIn(new contentVecType_32p(stxxlLibFile.get(), iSizeOfLib));

					// add kMers from fasta
					size_t overallCharsRead = 0;
					Build dummy;
					if (sDirectory.back() == '/') {
						auto filesAndSize = Utilities::gatherFilesFromPath(sDirectory);
						for (auto& fileName : filesAndSize.first) {
							if (_bVerbose) {
								cout << "OUT: Current file: " << fileName << endl;
							}
							ifstream fastaFile(fileName);
							fastaFile.seekg(0, fastaFile.end);
							const uint64_t& iFileLength = fastaFile.tellg();
							fastaFile.seekg(0, fastaFile.beg);
							readFasta(fastaFile, mAccToID, dummy, vLibIn, iFileLength, overallCharsRead, filesAndSize.second, fPercentageOfThrowAway);
						}
					}
					else {
						ifstream fastaFile;
						//fastaFile.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
						fastaFile.open(sDirectory);
						if (_bVerbose) {
							cout << "OUT: Current file: " << sDirectory << endl;
						}
						fastaFile.seekg(0, fastaFile.end);
						const uint64_t& iFileLength = fastaFile.tellg();
						fastaFile.seekg(0, fastaFile.beg);
						readFasta(fastaFile, mAccToID, dummy, vLibIn, iFileLength, overallCharsRead, iFileLength, fPercentageOfThrowAway);
					}


					stxxl::sort(vLibIn->begin(), vLibIn->end(), SCompareStructForSTXXLSort(), iMemory);

					auto funcForUnique = [](packedBigPair& a, packedBigPair& b) {
						if (a == b) {
							return true;
						}
						return false;
					};
					contentVecType_32p::iterator newEnd = std::unique(vLibIn->begin(), vLibIn->end(), funcForUnique);
					vLibIn->resize(newEnd - vLibIn->begin(), true);

					/*const auto& vFailed = verifyIntegrity(vLibIn, mIDsAsIdx, iIdxCounter);
					if (vFailed.size()) {
						cout << "Warning! The following disappeared:" << endl;
						for (const auto& entry : vFailed) {
							cout << mIdxToName[entry] << endl;
						}
					}*/

					ofstream fOutInfo(fOutFile + "_info.txt");
					fOutInfo << vLibIn->size();

					vLibIn->export_files("_");

					Trie T(static_cast<int8_t>(12), static_cast<int8_t>(_iMinK), 6);
					T.SaveToStxxlVec(vLibIn.get(), fOutFile);
				}
				else {
					// get Size of Lib and open it for reading
					fstream fLibInfo(sLibFile + "_info.txt", ios::in);
					uint64_t iSizeOfLib = 0;
					fLibInfo >> iSizeOfLib;
					fLibInfo.close();
					unique_ptr<stxxlFile> stxxlLibFile(new stxxlFile(sLibFile, stxxl::file::RDONLY));
					unique_ptr<const contentVecType_32p> vLibIn(new contentVecType_32p(stxxlLibFile.get(), iSizeOfLib));
					contentVecType_32p::bufreader_type vCBuff(*vLibIn);

					ofstream derp;
					//derp.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
					derp.open(fOutFile);
					derp.close();
					unique_ptr<stxxlFile> stxxlOutVec(new stxxlFile(fOutFile, stxxl::file::RDWR));
					unique_ptr<contentVecType_32p> vOutVec(new contentVecType_32p(stxxlOutVec.get(), 0));
					contentVecType_32p::bufwriter_type vNCBuff(*vOutVec);

					// add kMers from fasta
					derp.open(_sTemporaryPath + "_tempUpdate_" + to_string(_iNumOfCall));
					derp.close();
					unique_ptr<stxxlFile> stxxlTempVec(new stxxlFile(_sTemporaryPath + "_tempUpdate_" + to_string(_iNumOfCall), stxxl::file::RDWR));
					unique_ptr<contentVecType_32p> vTempVec(new contentVecType_32p(stxxlTempVec.get(), 0));
					
					size_t overallCharsRead = 0;
					Build dummy;
					if (sDirectory.back() == '/') {
						auto filesAndSize = Utilities::gatherFilesFromPath(sDirectory);
						for (auto& fileName : filesAndSize.first) {
							if (_bVerbose) {
								cout << "OUT: Current file: " << fileName << endl;
							}
							ifstream fastaFile(fileName);
							fastaFile.seekg(0, fastaFile.end);
							const uint64_t& iFileLength = fastaFile.tellg();
							fastaFile.seekg(0, fastaFile.beg);
							readFasta(fastaFile, mAccToID, dummy, vTempVec, iFileLength, overallCharsRead, filesAndSize.second, fPercentageOfThrowAway);
						}
					}
					else {
						ifstream fastaFile;
						//fastaFile.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
						fastaFile.open(sDirectory);
						if (_bVerbose) {
							cout << "OUT: Current file: " << sDirectory << endl;
						}
						fastaFile.seekg(0, fastaFile.end);
						const uint64_t& iFileLength = fastaFile.tellg();
						fastaFile.seekg(0, fastaFile.beg);
						readFasta(fastaFile, mAccToID, dummy, vTempVec, iFileLength, overallCharsRead, iFileLength, fPercentageOfThrowAway);
					}


					stxxl::sort(vTempVec->begin(), vTempVec->end(), SCompareStructForSTXXLSort(), iMemory);
					/*
					auto funcForUnique = [](packedBigPair& a, packedBigPair& b) {
						if (a == b) {
							return true;
						}
						return false;
					};
					contentVecType_32p::iterator newEnd = std::unique(vTempVec->begin(), vTempVec->end(), funcForUnique);
					vTempVec->resize(newEnd - vTempVec->begin(), true);*/ // will be done by merge
					// Merge existing db and new one
					const auto& vOutSize = merge(vNCBuff, vCBuff, iSizeOfLib, vTempVec->cbegin(), vTempVec->cend(), arrFrequencies_, mIDsAsIdx, true);
					vOutVec->resize(vOutSize, true);

					ofstream fOutInfo(fOutFile + "_info.txt");
					fOutInfo << vOutVec->size();

					ofstream outFile(fOutFile + "_f.txt");
					for (uint32_t j = 0; j < iIdxCounter; ++j) {
						outFile << Utilities::checkIfInMap(mIdxToName, j)->second << "\t";
						outFile << arrFrequencies_[j * 12];
						for (int32_t k = 1; k < 12; ++k) {
							outFile << "\t" << arrFrequencies_[j * 12 + k];
						}
						outFile << endl;
					}

					vOutVec->export_files("_");

					Trie T(static_cast<int8_t>(12), static_cast<int8_t>(_iMinK), 6);
					T.SaveToStxxlVec(vOutVec.get(), fOutFile);
				}
				remove((_sTemporaryPath + "_tempUpdate_" + to_string(_iNumOfCall)).c_str());
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

	};

}