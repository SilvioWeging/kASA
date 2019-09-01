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

#include "kASA.hpp"
#include "Build.hpp"
#include "Trie.hpp"

namespace kASA {
	class Read : public kASA {
		const bool _bInputAreAAs;

	public:
		Read(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const bool& bTranslated = false, const string& stxxl_mode = "") : kASA(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, stxxl_mode), _bInputAreAAs(bTranslated) {}


	protected:

		struct strTransfer {
			string name = "", overhang = "";
			size_t lengthOfDNA = 0;
			bool finished = true, addTail = false;
			uint8_t iExpectedInput = 0;
			pair<string, readIDType> lastLine;
			list<readIDType> vReadIDs;
			readIDType iCurrentReadID = 0;
			unordered_map<readIDType, uint64_t> mReadIDToArrayIdx;
			uint64_t iNumOfCharsRead = 0, iCurrentPercentage = 0, iNumOfAllCharsRead = 0, iCurrentOverallPercentage = 0, iNumOfNewReads = 0;
		};

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fastq as input for comparison
		template<typename T>
		inline uint64_t readFastq_partialSort(T& input, unordered_map<uint64_t, Utilities::rangeContainer>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const bool& bSpaced, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, unique_ptr<strTransfer>& transfer, const Trie& Tr) {
			uint64_t iSumOfkMers = 0;
			try {
				bool bNotFull = true;
				uint64_t iLocalReadID = 0;

				string sFalsekMerMarker = "";
				if (_bInputAreAAs) {
					for (int32_t i = 0; i < (_iMaxK - _iMinK); ++i) {
						sFalsekMerMarker += "^";
					}
				}
				else {
					for (int32_t i = 0; i < (_iMaxK - _iMinK) * 3; ++i) {
						sFalsekMerMarker += "X";
					}
				}
				////////////////////////////////////////////////////////
				uint8_t iExpectedInput = transfer->iExpectedInput; // 0 = name, 1 = sequence, 2 = + or quality
				string sName = transfer->name, sDNA = "", sOverhang = transfer->overhang;
				size_t iDNALength = transfer->lengthOfDNA, iQualityLength = 0;
				bool bNewRead = transfer->finished, bAddTail = transfer->addTail;

				size_t iLinesPerRun = 2;
				int32_t iProcID = 0;

				vector<pair<string, readIDType>> vLines(iLinesPerRun), vRCLines(iLinesPerRun);
				size_t iActualLines = 1;
				vLines[0] = transfer->lastLine;
				transfer->iNumOfNewReads = 0;

				while (input && bNotFull) {
					while (input && iActualLines < iLinesPerRun) {
						uint64_t iNumOfChars = 0;
						pair<std::string, bool> resultChunkPair = Utilities::getChunk(input, move(iNumOfChars));

						if (_bVerbose && iFileLength != 0) {
							transfer->iNumOfCharsRead += iNumOfChars;
							transfer->iNumOfAllCharsRead += iNumOfChars;
							double dPercentageOfInputRead = transfer->iNumOfCharsRead / double(iFileLength) * 100.;
							if (static_cast<uint64_t>(dPercentageOfInputRead) != transfer->iCurrentPercentage) {
								transfer->iCurrentPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
								cout << "OUT: Progress of current file " << transfer->iCurrentPercentage <<  "%" << endl;
							}
							if (iFileLength != overallFilesSize) {
								dPercentageOfInputRead = transfer->iNumOfAllCharsRead / double(overallFilesSize) * 100.;
								if (static_cast<uint64_t>(dPercentageOfInputRead) != transfer->iCurrentOverallPercentage) {
									transfer->iCurrentOverallPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
									cout << "OUT: Progress of all files " << transfer->iCurrentOverallPercentage << " %" << endl;
								}
							}
						}

						if (resultChunkPair.first.length() == 0) {
							continue;
						}

						if (resultChunkPair.first.front() == '+' && iExpectedInput != 2) {
							if (!resultChunkPair.second) {
								input.ignore(numeric_limits<streamsize>::max(), '\n'); // discard the rest of the +
							}
							iExpectedInput = 2;
							continue;
						}

						const uint8_t& iCurrentExpInput = iExpectedInput;
						switch (iCurrentExpInput) {

						case 0:
							sName = Utilities::lstrip(resultChunkPair.first, '@');
							if (!resultChunkPair.second) {
								input.ignore(numeric_limits<streamsize>::max(), '\n'); // discard the rest of the name
							}
							bNewRead = true;
							sOverhang = "";
							iDNALength = 0;
							iQualityLength = 0;

							//vReadNameAndLength.reserve(vReadNameAndLength.capacity() + 1);

							if (bReadIDsAreInteresting) {
								transfer->vReadIDs.push_back(transfer->iCurrentReadID);
								transfer->mReadIDToArrayIdx[transfer->iCurrentReadID] = iLocalReadID;

								iSoftMaxSize -= iAmountOfSpecies * sizeof(float) + sizeof(readIDType) + sizeof(pair<readIDType, uint64_t>);
							}
							transfer->iNumOfNewReads++;
							//cout << iSoftMaxSize << endl;

							iExpectedInput = 1;
							break;

						case 1:
							for (char& c : resultChunkPair.first) {
								if (c == '\t' || c == ' ') {
									throw runtime_error("Spaces or tabs inside read, please check your input. Error occured in:\n" + sName);
								}
								else {
									if (_bInputAreAAs) {
										if (c == '*') {
											c = '[';
										}
									}
									else {
										if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'a' && c != 'c' && c != 'g' && c != 't') {
											c = 'Z';
										}
									}
								}
							}

							sDNA = sOverhang + resultChunkPair.first;
							iDNALength += resultChunkPair.first.length();
							if (_bInputAreAAs) {
								if (sDNA.length() < size_t(_iHighestK)) {
									sOverhang = sDNA;
								}
								else {
									sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK);
								}
							}
							else {
								if (sDNA.length() < size_t(_iHighestK) * 3) {
									sOverhang = sDNA;
								}
								else {
									sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK * 3);
								}
							}

							if (bNewRead) {
								bNewRead = false;
								if (!_bInputAreAAs) {
									string tempDNA = reverseComplement(sDNA);
									while (tempDNA.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
										tempDNA += 'X';
									}
									vRCLines[iActualLines - 1] = (bSpaced) ? make_pair(tempDNA, transfer->iCurrentReadID) : make_pair(tempDNA + sFalsekMerMarker, transfer->iCurrentReadID);
								}
							}
							else {
								if (!_bInputAreAAs) {
									vRCLines[iActualLines - 1] = (bSpaced) ? make_pair(reverseComplement(sDNA), transfer->iCurrentReadID) : make_pair(reverseComplement(sDNA), transfer->iCurrentReadID);
								}
							}

							vLines[iActualLines] = (bSpaced) ? make_pair(sDNA, transfer->iCurrentReadID) : make_pair(sDNA, transfer->iCurrentReadID);
							++iActualLines;

							bAddTail = true;

							break;

						case 2:
							if (bAddTail) {
								// Pad very small reads
								if (_bInputAreAAs) {
									while (vLines[iActualLines - 1].first.length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
										vLines[iActualLines - 1].first += '^';
									}
								}
								else {
									while (vLines[iActualLines - 1].first.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
										vLines[iActualLines - 1].first += 'X';
									}
								}

								vLines[iActualLines - 1] = (bSpaced) ? vLines[iActualLines - 1] : make_pair(vLines[iActualLines - 1].first + sFalsekMerMarker, vLines[iActualLines - 1].second); // add false marker
								if (bReadIDsAreInteresting) {
									++(transfer->iCurrentReadID);
									++iLocalReadID;
									vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength + sFalsekMerMarker.length())));
									iSoftMaxSize -= sName.size() * sizeof(char) + sizeof(uint32_t);
								}
								bAddTail = false;
							}

							iQualityLength += resultChunkPair.first.length();


							if (iQualityLength == iDNALength) {
								iExpectedInput = 0;
							}
							else {
								if (iQualityLength > iDNALength) {
									throw runtime_error("Quality string length and DNA length don't match. Error occured in:\n" + sName);
								}
							}

							break;
						default:
							break;
						}
					}

					const auto& iNumOfkMers = convertLinesTokMers_(int32_t(iActualLines - 1), vLines, vRCLines, vOut, Tr, bSpaced);
					iSumOfkMers += iNumOfkMers;
					iProcID = (iProcID + 1) % _iNumOfThreads;
					//iSoftMaxSize -= iNumOfkMers * (sizeof(pair<uint32_t, uint32_t>) + sizeof(uint64_t));
					if (_iMinK <= 6) {
						iSoftMaxSize -= iNumOfkMers * (sizeof(pair<uint64_t, uint64_t>) + sizeof(Utilities::rangeContainer));
					}
					else {
						iSoftMaxSize -= iNumOfkMers * (sizeof(pair<uint32_t, uint32_t>) + sizeof(Utilities::rangeContainer));
					}
					//cout << iSoftMaxSize << endl;
					if (iSoftMaxSize <= 0) {
						bNotFull = false;
					}


					vLines[0] = vLines[iActualLines - 1];
					iActualLines = 1;
				}
				if (!input) {
					vRCLines[0] = make_pair("", 0);
					iSumOfkMers += convertLinesTokMers_(1, vLines, vRCLines, vOut, Tr, bSpaced);
					transfer->addTail = false; //signal that all is done
					transfer->lastLine = vLines[0];
				}
				else {
					transfer->lastLine = vLines[0];
					transfer->addTail = bAddTail;
					transfer->finished = bNewRead;
					transfer->iExpectedInput = iExpectedInput;
					transfer->lengthOfDNA = iDNALength;
					transfer->name = sName;
					transfer->overhang = sOverhang;
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
			////////////////////////////////////////////////////////
			return iSumOfkMers;
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fasta as input for comparison, no parallelization possible due to the unknown length of a contig/genome
		template<typename T>
		inline uint64_t readFasta_partialSort(T& input, unordered_map<uint64_t, Utilities::rangeContainer>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const bool& bSpaced, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, unique_ptr<strTransfer>& transfer, const Trie& Tr) {
			uint64_t iSumOfkMers = 0;
			try {
				bool bNotFull = true;
				uint64_t iLocalReadID = 0;

				/*for (int16_t i = 0; i < _iNumOfThreads; ++i) {
					vOut[i].reserve(static_cast<size_t>((iSoftMaxSize*1.1)/(sizeof(tuple<uint64_t, readIDType>) * (_iNumOfThreads))));
				}*/

				string sFalsekMerMarker = "";
				if (_bInputAreAAs) {
					for (int32_t i = 0; i < (_iMaxK - _iMinK); ++i) {
						sFalsekMerMarker += "^";
					}
				}
				else {
					for (int32_t i = 0; i < (_iMaxK - _iMinK) * 3; ++i) {
						sFalsekMerMarker += "X";
					}
				}
				////////////////////////////////////////////////////////
				uint8_t iExpectedInput = transfer->iExpectedInput; // 0 = name, 1 = sequence, 2 = signal end
				string sName = transfer->name, sDNA = "", sOverhang = transfer->overhang;
				size_t iDNALength = transfer->lengthOfDNA;
				bool bNewRead = transfer->finished, bAddTail = transfer->addTail;

				size_t iLinesPerRun = 2;
				int32_t iProcID = 0;

				vector<pair<string, readIDType>> vLines(iLinesPerRun), vRCLines(iLinesPerRun);
				size_t iActualLines = 1;
				vLines[0] = transfer->lastLine;

				while (input && bNotFull) {
					while (input && iActualLines < iLinesPerRun) {
						uint64_t iNumOfChars = 0;
						pair<std::string, bool> resultChunkPair = Utilities::getChunk(input, move(iNumOfChars));

						if (_bVerbose && iFileLength != 0) {
							transfer->iNumOfCharsRead += iNumOfChars;
							transfer->iNumOfAllCharsRead += iNumOfChars;
							double dPercentageOfInputRead = transfer->iNumOfCharsRead / double(iFileLength) * 100.;
							if (static_cast<uint64_t>(dPercentageOfInputRead) != transfer->iCurrentPercentage) {
								transfer->iCurrentPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
								cout << "OUT: Progress of current file " << transfer->iCurrentPercentage << " %" << endl;
							}
							if (iFileLength != overallFilesSize) {
								dPercentageOfInputRead = transfer->iNumOfAllCharsRead / double(overallFilesSize) * 100.;
								if (static_cast<uint64_t>(dPercentageOfInputRead) != transfer->iCurrentOverallPercentage) {
									transfer->iCurrentOverallPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
									cout << "OUT: Progress of all files " << transfer->iCurrentOverallPercentage << " %" << endl;
								}
							}
						}

						if (resultChunkPair.first.length() == 0) {
							continue;
						}

						if (resultChunkPair.first.front() == '>') {
							if (!resultChunkPair.second) {
								input.ignore(numeric_limits<streamsize>::max(), '\n'); // discard the rest of the name
							}
							iExpectedInput = 0;
						}

						const uint8_t& iCurrentExpInput = iExpectedInput; // In C++17, this can be written inside the switch brackets
						switch (iCurrentExpInput) {

						case 0:
							if (bAddTail) {
								// Pad very small reads
								if (_bInputAreAAs) {
									while (vLines[iActualLines - 1].first.length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
										vLines[iActualLines - 1].first += '^';
									}
								}
								else {
									while (vLines[iActualLines - 1].first.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
										vLines[iActualLines - 1].first += 'X';
									}
								}
								vLines[iActualLines - 1] = (bSpaced) ? vLines[iActualLines - 1] : make_pair(vLines[iActualLines - 1].first + sFalsekMerMarker, vLines[iActualLines - 1].second); // add false marker
								if (bReadIDsAreInteresting) {
									++(transfer->iCurrentReadID);
									++iLocalReadID;
									vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength + sFalsekMerMarker.length())));
									iSoftMaxSize -= sName.size() * sizeof(char) + sizeof(uint32_t);
								}
								bAddTail = false;
							}

							sName = Utilities::lstrip(resultChunkPair.first, '>');

							bNewRead = true;
							sOverhang = "";
							iDNALength = 0;

							//vReadNameAndLength.reserve(vReadNameAndLength.capacity() + 1);
							if (bReadIDsAreInteresting) {
								transfer->vReadIDs.push_back(transfer->iCurrentReadID);
								transfer->mReadIDToArrayIdx[transfer->iCurrentReadID] = iLocalReadID;

								iSoftMaxSize -= iAmountOfSpecies * sizeof(float) + sizeof(readIDType) + sizeof(pair<readIDType, uint64_t>);
							}

							iExpectedInput = 1;
							break;

						case 1:
							for (char& c : resultChunkPair.first) {
								if (c == '\t' || c == ' ') {
									throw runtime_error("Spaces or tabs inside read, please check your input. Error occured in:" + sName);
								}
								else {
									if (_bInputAreAAs) {
										if (c == '*') {
											c = '[';
										}
									}
									else {
										if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'a' && c != 'c' && c != 'g' && c != 't') {
											c = 'Z';
										}
									}
								}
							}

							sDNA = sOverhang + resultChunkPair.first;
							iDNALength += resultChunkPair.first.length();
							if (_bInputAreAAs) {
								if (sDNA.length() < size_t(_iHighestK)) {
									sOverhang = sDNA;
								}
								else {
									sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK);
								}
							}
							else {
								if (sDNA.length() < size_t(_iHighestK) * 3) {
									sOverhang = sDNA;
								}
								else {
									sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK * 3);
								}
							}

							if (bNewRead) {
								bNewRead = false;
								if (!_bInputAreAAs) {
									string tempDNA = reverseComplement(sDNA);
									while (tempDNA.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
										tempDNA += 'X';
									}
									vRCLines[iActualLines - 1] = (bSpaced) ? make_pair(tempDNA, transfer->iCurrentReadID) : make_pair(tempDNA + sFalsekMerMarker, transfer->iCurrentReadID);
								}
							}
							else {
								if (!_bInputAreAAs) {
									vRCLines[iActualLines - 1] = (bSpaced) ? make_pair(reverseComplement(sDNA), transfer->iCurrentReadID) : make_pair(reverseComplement(sDNA), transfer->iCurrentReadID);
								}
							}

							vLines[iActualLines] = (bSpaced) ? make_pair(sDNA, transfer->iCurrentReadID) : make_pair(sDNA, transfer->iCurrentReadID);
							++iActualLines;

							bAddTail = true;

							break;
						default:
							break;
						}
					}

					const auto& iNumOfkMers = convertLinesTokMers_(int32_t(iActualLines - 1), vLines, vRCLines, vOut, Tr, bSpaced);
					iSumOfkMers += iNumOfkMers;
					iProcID = (iProcID + 1) % _iNumOfThreads;
					if (_iMinK <= 6) {
						iSoftMaxSize -= iNumOfkMers * (sizeof(pair<uint64_t, uint64_t>) + sizeof(Utilities::rangeContainer));
					}
					else {
						iSoftMaxSize -= iNumOfkMers * (sizeof(pair<uint32_t, uint32_t>) + sizeof(Utilities::rangeContainer));
					}
					
					if (iSoftMaxSize <= 0) {
						bNotFull = false;
					}


					vLines[0] = vLines[iActualLines - 1];
					iActualLines = 1;
				}
				if (!input) {
					if (bAddTail) {
						// when finished, convert the last part

						// Pad very small reads
						if (_bInputAreAAs) {
							while (vLines[iActualLines - 1].first.length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
								vLines[iActualLines - 1].first += '^';
							}
						}
						else {
							while (vLines[iActualLines - 1].first.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
								vLines[iActualLines - 1].first += 'X';
							}
						}

						vLines[iActualLines - 1] = (bSpaced) ? vLines[iActualLines - 1] : make_pair(vLines[iActualLines - 1].first + sFalsekMerMarker, vLines[iActualLines - 1].second);
						vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength + sFalsekMerMarker.length())));
						iSoftMaxSize -= sName.size() * sizeof(char) + sizeof(uint32_t);
						bAddTail = false;
					}
					vRCLines[0] = make_pair("", 0);
					iSumOfkMers += convertLinesTokMers_(1, vLines, vRCLines, vOut, Tr, bSpaced);
					transfer->lastLine = vLines[0];
					transfer->addTail = bAddTail;
				}
				else {
					transfer->lastLine = vLines[0];
					transfer->addTail = bAddTail;
					transfer->finished = bNewRead;
					transfer->iExpectedInput = iExpectedInput;
					transfer->lengthOfDNA = iDNALength;
					transfer->name = sName;
					transfer->overhang = sOverhang;
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
			////////////////////////////////////////////////////////
			return iSumOfkMers;
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fasta and create a kMer-Vec, used in UpdateFromFasta(...)
		inline void readFasta(ifstream& input, unique_ptr<contentVecType_32p>& vOut, const float& fShrinkPercentage = 0.f, const bool& bSpaced = false, const unordered_map<string, uint32_t>& mAccToID = unordered_map<string, uint32_t>()) {
			try {
				if (!input.good()) {
					throw runtime_error("No input found!");
				}

				uint32_t iIdx = 0, iNextIdx = 0;

				bool bAccNumberFound = false;

				string sFalsekMerMarker = "";
				for (int32_t i = 0; i < (_iHighestK - _iLowestK) * 3; ++i) {
					sFalsekMerMarker += "X";
				}

				string sInitialString;
				while (input.good()) {
					getline(input, sInitialString);
					if (sInitialString == "") {
						continue;
					}
					if (sInitialString.front() == '>') {
						sInitialString.erase(sInitialString.begin());
						const auto& sNumbers = Utilities::split(Utilities::split(sInitialString, ' ').at(0), '|');
						string acc = "";
						for (const auto& entry : sNumbers) {
							if (entry.find('.') != string::npos) {
								acc = entry;
								break;
							}
						}
						auto isInMap = mAccToID.find(acc);
						if (isInMap != mAccToID.end()) {
							iIdx = isInMap->second;
							bAccNumberFound = true;
						}
						else {
							bAccNumberFound = false;
						}

						break;
					}
					else {
						throw runtime_error("No > found in input.");
					}
				}



				bool bAddFalseMarker = false, bAddFalseRCMarker = true, bNextIdx = false;

				const int32_t& iConcurrentLines = 20;
				string sOverhang = "", sDNA = "";


				// Either the file still contains dna or there is an Overhang left from a previous loop
				while (input.good() || sOverhang != "") {

					if (bNextIdx) {
						bAddFalseMarker = false;
						bAddFalseRCMarker = true;
					}

					bool bCreateOverhang = true;

					for (int32_t i = 0; i < iConcurrentLines; ++i) {
						char sTempCArr[100000];
						input.get(sTempCArr, 100000);
						auto iNumOfChars = input.gcount();
						bool bLineNotFinished = true;
						if (iNumOfChars != 99999 && input) {
							bLineNotFinished = false;
							input.get();
						}

						string sTempString(sTempCArr);
						if (sTempString == "") {
							if (input.fail() && !input.eof()) {
								input.clear();
								getline(input, sTempString);
							}
							else {
								bCreateOverhang = false;
								bAddFalseMarker = true;
								break;
							}
							--i;
							continue;
						}

						if (input) {

							//string sTempString = "";
							//if (getline(input, sTempString)) {

							// new entry in fasta
							if (sTempString.front() == '>') {
								sTempString.erase(sTempString.begin());
								while (bLineNotFinished) {
									// read the whole line
									input.get(sTempCArr, 100000);
									iNumOfChars = input.gcount();
									if (iNumOfChars != 99999 && input) {
										bLineNotFinished = false;
										input.get();
									}
									sTempString.append(sTempString);
								}
								const auto& sNumbers = Utilities::split(Utilities::split(sTempString, ' ').at(0), '|');
								string acc = "";
								for (const auto& entry : sNumbers) {
									if (entry.find('.') != string::npos) {
										acc = entry;
										break;
									}
								}
								auto isInMap = mAccToID.find(acc);
								if (isInMap != mAccToID.end()) {
									iNextIdx = isInMap->second;
									bAccNumberFound = true;
									bNextIdx = true;
								}
								else {
									bAccNumberFound = false;
								}

								bCreateOverhang = false;
								bAddFalseMarker = true;
								//bAddFalseRCMarker = true;
								break;
							}
							else {
								if (bAccNumberFound) {
									// append to string
									sDNA += sTempString;
								}
								else {
									--i;
								}
							}
						}
						else {
							// EOF
							bCreateOverhang = false;
							bAddFalseMarker = true;
							//bAddFalseRCMarker = false;
							break;
						}
					}

					sDNA = sOverhang + sDNA;

					if (sDNA.length() < uint32_t(_iHighestK * 3 + 2) && bCreateOverhang) { // e.g. 38 is the last frame for a kMer of Maxlength 12
						sOverhang = sDNA;
						sDNA = "";
						continue;
					}

					for (auto& chara : sDNA) {
						if (chara != 'A' && chara != 'C' && chara != 'G' && chara != 'T' && chara != 'a' && chara != 'c' && chara != 'g' && chara != 't') {
							chara = 'Z';
						}
					}
					string sRCDNA = reverseComplement(sDNA);

					if (!bSpaced) {
						if (bAddFalseRCMarker && !sRCDNA.empty()) {
							sRCDNA += sFalsekMerMarker;
							bAddFalseRCMarker = false;
						}
						if (bAddFalseMarker && !sDNA.empty()) {
							sDNA += sFalsekMerMarker;
							bAddFalseMarker = false;
						}
					}

					// not parallel due to issues...
					uint32_t iActualLines = 1;
					vector<string> vLines = { sDNA }, vRCLines = { sRCDNA };

					if (bSpaced) {
						//convertLinesToSpacedkMers(iActualLines, vLines, iIdx, false, vOut);
						//convertLinesToSpacedkMers(iActualLines, vRCLines, iIdx, true, vOut);
					}
					else {
						convertLinesTokMers(iActualLines, vLines, iIdx, vOut, fShrinkPercentage);
						convertLinesTokMers(iActualLines, vRCLines, iIdx, vOut, fShrinkPercentage);
					}

					if (bNextIdx) {
						iIdx = iNextIdx;
						bNextIdx = false;
					}

					if (bCreateOverhang) {
						sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK * 3);
					}
					else {
						sOverhang = "";
					}

					sDNA = "";
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// convert dna to k-mers in the building step
		inline void dnaTokMers(const string& sDNA, const uint32_t& iIdx, Build& vBricks, const float& fShrinkPercentage) {
			auto&& vResultingkMers = createStuff<contentVecType_32p>(nullptr, 1);
			const int32_t& iMaxKTimes3 = 3 * _iHighestK;

			for (int32_t i = 0; i < 1; ++i) {
				const string& sDna = sDNA;
				const int32_t& iMaxRange = int32_t(sDna.length()) - iMaxKTimes3 + 1;
				if (iMaxRange >= 1) {
					// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
					const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;

					vResultingkMers[i].resize(iMaxRange);
					int32_t ikMerCounter = 0;

					// Frameshifting, so that only one amino acid has to be computed
					// Compute initial frames
					uint64_t sAAFrames[3] = { 0, 0, 0 };
					uint32_t aDeletekMerCounter[3] = { 0, 0, 0 };
					for (int32_t j = 0; j < iNumFrames; ++j) {
						string sTempFrame = _sMaxKBlank;
						dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
						sAAFrames[j] = aminoacidTokMer(sTempFrame);
						// if a character is 'illegal' then don't save the kMer containing it
						const auto& iPosOfU = sTempFrame.find_last_of('U');
						if (iPosOfU == string::npos) {
							vResultingkMers[i][ikMerCounter + j] = returnStuff<contentVecType_32p>(nullptr, sAAFrames[j], iIdx);
						}
						else {
							vResultingkMers[i][ikMerCounter + j] = returnStuff<contentVecType_32p>(nullptr, 0, iIdx);
							aDeletekMerCounter[j] = uint32_t(iPosOfU);
						}

					}

					ikMerCounter += iNumFrames;

					if (iMaxRange > 3) {

						// If the Dna has an irregular (not mod 3) tail, a special case must be handled and an additional step in the loop is necessary thus the !(iMaxRange % 3)
						// To map the 2 onto 1, ! is applied twice
						const int32_t& iMaxRangeMod3 = iMaxRange % 3;
						const int32_t& iMaxRangeMod3Neg = !(!iMaxRangeMod3);

						for (int32_t j = 1; 3 * (j + iMaxRangeMod3Neg) < iMaxRange; ++j) {
							for (int32_t k = 0; k < 3; ++k) {
								int8_t sTempAA = ' ';
								dnaToAminoacid(sDna, k + iMaxKTimes3 + 3 * (j - 1), sTempAA);
								sAAFrames[k] = aminoacidTokMer(sAAFrames[k], sTempAA);
								//string sDEBUG = kMerToAminoacid(sAAFrames[k], 12);
								if (sTempAA != 'U' && aDeletekMerCounter[k] == 0) {
									vResultingkMers[i][ikMerCounter + k] = returnStuff<contentVecType_32p>(nullptr, sAAFrames[k], iIdx);
								}
								else {
									if (aDeletekMerCounter[k] != 0 && sTempAA != 'U') {
										--aDeletekMerCounter[k];
									}
									else {
										aDeletekMerCounter[k] = _iHighestK - 1;
									}
									vResultingkMers[i][ikMerCounter + k] = returnStuff<contentVecType_32p>(nullptr, 0, iIdx);
								}
							}
							ikMerCounter += 3;
						}

						// compute frames of said tail
						for (int32_t j = 0; j < iMaxRangeMod3; ++j) {
							int8_t sTempAA = ' ';
							dnaToAminoacid(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
							sAAFrames[j] = aminoacidTokMer(sAAFrames[j], sTempAA);
							if (sTempAA != 'U' && aDeletekMerCounter[j] == 0) {
								vResultingkMers[i][ikMerCounter + j] = returnStuff<contentVecType_32p>(nullptr, sAAFrames[j], iIdx);
							}
							else {
								if (aDeletekMerCounter[j] != 0 && sTempAA != 'U') {
									--aDeletekMerCounter[j];
								}
								else {
									aDeletekMerCounter[j] = _iHighestK - 1;
								}
								vResultingkMers[i][ikMerCounter + j] = returnStuff<contentVecType_32p>(nullptr, 0, iIdx);
							}
						}
					}
				}
			}

			/*uint64_t iResultingSize = 0;
			for (int32_t i = 0; i < 1; ++i) {
			iResultingSize += vResultingkMers[i].size();
			}*/

			// Write resulting kMers in stxxl vector but exclude those which are not useful
			double dStepSize = (fShrinkPercentage > 0.f) ? 100. / fShrinkPercentage : 0.;
			double dNextThrowOut = dStepSize;
			uint64_t iCounterOfThrowOut = 1;
			uint32_t iSizeCounterOfVec = 0;

			if (fShrinkPercentage > 0.f) {
				for (const auto& kMerResVec : vResultingkMers) {
					for (const auto& element : kMerResVec) {
						if (isNotZero<contentVecType_32p>(nullptr, element)) {
							if (iCounterOfThrowOut != static_cast<uint64_t>(dNextThrowOut)) {
								/*if (get<0>(element) == 37191153727705383) {
								cout << "derp" << endl;
								}*/
								/*bool bAdded = false;
								for (int i = 0; i < _iNumOfThreads - 1; ++i) {
								if (get<0>(element) < vDivisionArray[i]) {
								if (!(vBricks[i]->addToInt(element))) {

								vBricks[i]->IntToExt2();
								}
								bAdded = true;
								break;
								}
								}
								if (!bAdded) {*/
								//if (!(vBricks[_iNumOfThreads - 1]->addToInt(element))) {
								//	vBricks[_iNumOfThreads - 1]->IntToExt2();
								//}
								if (!(vBricks.addToInt(element))) {
									vBricks.IntToExt2();
								}
								//}
								++iSizeCounterOfVec;
							}
							else {
								dNextThrowOut += dStepSize;
							}
							++iCounterOfThrowOut;
						}
					}
				}
			}
			else {
				for (const auto& kMerResVec : vResultingkMers) {
					for (const auto& element : kMerResVec) {
						if (isNotZero<contentVecType_32p>(nullptr, element)) {
							if (!(vBricks.addToInt(element))) {
								vBricks.IntToExt2();
							}
							++iSizeCounterOfVec;
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// convert protein sequences to k-mers in the building step
		inline void proteinTokMers(const string& sAASequence, const uint32_t& iIdx, Build& vBricks, const float& fShrinkPercentage) {
			auto&& vResultingkMers = createStuff<contentVecType_32p>(nullptr, 1);
			const int32_t& iMaxRange = int32_t(sAASequence.length()) - _iHighestK + 1;
			if (iMaxRange > 0) {
				vResultingkMers[0].resize(iMaxRange);
				for (int32_t i = 0; i < iMaxRange; ++i) {
					const auto& kMer = aminoacidTokMer(sAASequence.cbegin() + i, sAASequence.cbegin() + i + 12);
					vResultingkMers[0][i] = returnStuff<contentVecType_32p>(nullptr, kMer, iIdx);
				}
			}

			// Write resulting kMers in stxxl vector but exclude those which are not useful
			double dStepSize = (fShrinkPercentage > 0.f) ? 100. / fShrinkPercentage : 0.;
			double dNextThrowOut = dStepSize;
			uint64_t iCounterOfThrowOut = 1;

			if (fShrinkPercentage > 0.f) {
				for (const auto& kMerResVec : vResultingkMers) {
					for (const auto& element : kMerResVec) {
						if (iCounterOfThrowOut != static_cast<uint64_t>(dNextThrowOut)) {
							if (!(vBricks.addToInt(element))) {
								vBricks.IntToExt2();
							}
						}
						else {
							dNextThrowOut += dStepSize;
						}
						++iCounterOfThrowOut;
					}
				}
			}
			else {
				for (const auto& kMerResVec : vResultingkMers) {
					for (const auto& element : kMerResVec) {
						if (!(vBricks.addToInt(element))) {
							vBricks.IntToExt2();
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fasta and create a kMer-Vec, used in BuildAll(...)
		inline void readFasta(ifstream& input, const unordered_map<string, uint32_t>& mAccToID, Build& vBricks, const uint64_t& iFileLength, size_t& overallCharsRead, const size_t& overallFilesSize, const float& fShrinkPercentage) {
			try {
				if (!input.good()) {
					throw runtime_error("No input found!");
				}

				uint32_t iIdx = 0, nextIdx = 0;

				bool bAccNumberFound = false;

				string sFalsekMerMarker = "";
				if (_bInputAreAAs) {
					for (int32_t i = 0; i < (_iHighestK - _iLowestK); ++i) {
						sFalsekMerMarker += "^";
					}
				}
				else {
					for (int32_t i = 0; i < (_iHighestK - _iLowestK) * 3; ++i) {
						sFalsekMerMarker += "X";
					}
				}

				string sInitialString;
				while (input.good()) {
					getline(input, sInitialString);
					if (sInitialString == "") {
						continue;
					}
					if (sInitialString.front() == '>') {
						sInitialString.erase(sInitialString.begin());
						const auto& sNumbers = Utilities::split(Utilities::split(sInitialString, ' ').at(0), '|');
						string acc = "";
						for (const auto& entry : sNumbers) {
							if (entry.find('.') != string::npos) {
								acc = entry;
								break;
							}
						}
						auto isInMap = mAccToID.find(acc);
						if (isInMap != mAccToID.end()) {
							iIdx = isInMap->second;
							bAccNumberFound = true;
						}
						else {
							// search if acc is a dummy
							isInMap = mAccToID.find(sInitialString);
							if (isInMap != mAccToID.end()) {
								iIdx = isInMap->second;
								bAccNumberFound = true;
							}
							else {
								bAccNumberFound = false;
							}
						}
						break;
					}
					else {
						throw runtime_error("No > found in input.");
					}
				}



				bool bAddFalseMarker = false,  bNextIdx = false; //bAddFalseRCMarker = true,

				const int32_t& iConcurrentLines = 20;
				string sOverhang = "", sDNA = "";

				uint64_t iNumOfCharsRead = 0, iCurrentPercentage = 0, iCurrentOverallPercentage = 0;

				// Either the file still contains dna or there is an Overhang left from a previous loop
				while (input.good() || sOverhang != "") {

					if (bNextIdx) {
						bAddFalseMarker = false;
						//bAddFalseRCMarker = true;
					}

					bool bCreateOverhang = true;

					for (int32_t i = 0; i < iConcurrentLines; ++i) {
						char sTempCArr[100000];
						input.get(sTempCArr, 100000);

						auto iNumOfChars = input.gcount();
						if (_bVerbose && iFileLength != 0) {
							iNumOfCharsRead += iNumOfChars;
							overallCharsRead += iNumOfChars;
							double dPercentageOfInputRead = iNumOfCharsRead / double(iFileLength) * 100.;
							if (static_cast<uint64_t>(dPercentageOfInputRead) != iCurrentPercentage) {
								iCurrentPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
								cout << "OUT: Progress of current file " << iCurrentPercentage << " %" << endl;
							}
							dPercentageOfInputRead = overallCharsRead / double(overallFilesSize) * 100.;
							if (static_cast<uint64_t>(dPercentageOfInputRead) != iCurrentOverallPercentage) {
								iCurrentOverallPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
								cout << "OUT: Progress of all files " << iCurrentOverallPercentage << " %" << endl;
							}
						}

						bool bLineNotFinished = true;
						if (iNumOfChars != 99999 && input) {
							bLineNotFinished = false;
							input.get();
						}
						string sTempString(sTempCArr);
						if (sTempString == "") {
							if (input.fail() && !input.eof()) {
								input.clear();
								getline(input, sTempString);
							}
							else {
								bCreateOverhang = false;
								bAddFalseMarker = true;
								break;
							}
							--i;
							continue;
						}
						if (input) {

							// new entry in fasta
							if (sTempString.front() == '>') {

								sTempString.erase(sTempString.begin());
								while (bLineNotFinished) {
									// read the whole line
									input.get(sTempCArr, 100000);
									iNumOfChars = input.gcount();
									if (iNumOfChars != 99999 && input) {
										bLineNotFinished = false;
										input.get();
									}
									sTempString.append(sTempString);
								}
								const auto& sNumbers = Utilities::split(Utilities::split(sTempString, ' ').at(0), '|');
								string acc = "";
								for (const auto& entry : sNumbers) {
									if (entry.find('.') != string::npos) {
										acc = entry;
										break;
									}
								}
								auto isInMap = mAccToID.find(acc);
								if (isInMap != mAccToID.end()) {
									nextIdx = isInMap->second;
									bAccNumberFound = true;
									bNextIdx = true;
								}
								else {
									// search if acc is a dummy
									isInMap = mAccToID.find(sTempString);
									if (isInMap != mAccToID.end()) {
										nextIdx = isInMap->second;
										bAccNumberFound = true;
										bNextIdx = true;
									}
									else {
										bAccNumberFound = false;
									}
								}

								bCreateOverhang = false;
								bAddFalseMarker = true;
								//bAddFalseRCMarker = true;
								break;
							}
							else {
								if (bAccNumberFound) {
									// append to string
									sDNA += sTempString;
								}
								else {
									--i;
								}
							}
						}
						else {
							if (bAccNumberFound) {
								// append to string
								sDNA += sTempString;
							}
							// EOF
							bCreateOverhang = false;
							bAddFalseMarker = true;
							//bAddFalseRCMarker = false;
							break;
						}
					}

					sDNA = sOverhang + sDNA;

					if (_bInputAreAAs) {
						if (sDNA.length() < size_t(_iHighestK + 2) && bCreateOverhang) { // e.g. 14 is the last frame for a kMer of Maxlength 12
							sOverhang = sDNA;
							sDNA = "";
							continue;
						}
					}
					else {
						if (sDNA.length() < size_t(_iHighestK * 3 + 2) && bCreateOverhang) { // e.g. 38 is the last frame for a kMer of Maxlength 12
							sOverhang = sDNA;
							sDNA = "";
							continue;
						}
					}

					for (auto& chara : sDNA) {
						if (chara == '\t' || chara == ' ') {
							throw runtime_error("Spaces or tabs inside reference, please check your input. Error occured in content entry Nr.: " + to_string(iIdx));
						}
						else {
							if (_bInputAreAAs) {
								if (chara == '*') {
									chara = '[';
								}
							}
							else {
								if (chara != 'A' && chara != 'C' && chara != 'G' && chara != 'T' && chara != 'a' && chara != 'c' && chara != 'g' && chara != 't') {
									chara = 'Z';
								}
							}
						}
					}
					//string sRCDNA = reverseComplement(sDNA);

					//if (bAddFalseRCMarker && !sRCDNA.empty()) {
					//	sRCDNA += sFalsekMerMarker;
					//	bAddFalseRCMarker = false;
					//}
					if (bAddFalseMarker && !sDNA.empty()) {
						sDNA += sFalsekMerMarker;
						bAddFalseMarker = false;
					}

					if (_bInputAreAAs) {
						proteinTokMers(sDNA, iIdx, vBricks, fShrinkPercentage);
					}
					else {
						dnaTokMers(sDNA, iIdx, vBricks, fShrinkPercentage);
						//dnaTokMers(sRCDNA, iIdx, vBricks, fShrinkPercentage);
					}

					if (bNextIdx) {
						iIdx = nextIdx;
						bNextIdx = false;
					}

					if (bCreateOverhang) {
						if (_bInputAreAAs) {
							sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK);
						}
						else {
							sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK * 3);
						}
					}
					else {
						sOverhang = "";
					}

					sDNA = "";
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Avoid copying code, instead use templates and just change the function for the respective instantiation
		template<typename vecType> inline vector<vector<tuple<uint64_t, uint32_t>>> createStuff(unique_ptr<vecType>&, const int32_t& iActualLines) {
			return vector<vector<tuple<uint64_t, uint32_t>>>(iActualLines, vector<tuple<uint64_t, uint32_t>>());
		}

		template<typename vecType> inline vector<vector<tuple<uint64_t, uint32_t>>> createStuff(vecType*, const int32_t& iActualLines) {
			return vector<vector<tuple<uint64_t, uint32_t>>>(iActualLines, vector<tuple<uint64_t, uint32_t>>());
		}

		inline vector<vector<tuple<uint64_t, readIDType>>> createStuff(unique_ptr<vector<tuple<uint64_t, readIDType>>>&, const int32_t& iActualLines) {
			return vector<vector<tuple<uint64_t, readIDType>>>(iActualLines, vector<tuple<uint64_t, readIDType>>());
		}


		template<typename vecType> inline tuple<uint64_t, uint32_t> returnStuff(unique_ptr<vecType>&, const uint64_t& kMer, const uint32_t& ID) {
			return tuple<uint64_t, uint32_t>(kMer, ID);
		}

		template<typename vecType> inline tuple<uint64_t, uint32_t> returnStuff(vecType*, const uint64_t& kMer, const uint32_t& ID) {
			return tuple<uint64_t, uint32_t>(kMer, ID);
		}

		inline tuple<uint64_t, uint64_t> returnStuff(unique_ptr<vector<tuple<uint64_t, readIDType>>>&, const uint64_t& kMer, const readIDType& ID) {
			return tuple<uint64_t, readIDType>(kMer, ID);
		}

		template<typename vecType, typename idType> inline bool isNotZero(unique_ptr<vecType>&, const tuple<uint64_t, idType>& elem) {
			return get<0>(elem) != 0;
		}

		template<typename vecType, typename idType> inline bool isNotZero(vecType*, const tuple<uint64_t, idType>& elem) {
			return get<0>(elem) != 0;
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// convert lines from fasta/fastq to kMers and save them in a (stxxl-)vector
		template<typename vecType, typename idType> inline void convertLinesTokMers(const int32_t& iActualLines, const vector<string>& vLines, const idType& iID, unique_ptr<vecType>& kMerVecOut, const float& fShrinkPercentage) {

			//vector<vector<tuple<uint64_t, uint32_t>>> vResultingkMers(iActualLines, vector<tuple<uint64_t, uint32_t>>());
			auto&& vResultingkMers = createStuff(kMerVecOut, iActualLines);
			const int32_t& iMaxKTimes3 = 3 * _iHighestK;

#pragma omp parallel for
			for (int32_t i = 0; i < iActualLines; ++i) {
				const string& sDna = vLines[i];
				const int32_t& iMaxRange = int32_t(sDna.length()) - iMaxKTimes3 + 1;
				if (iMaxRange >= 1) {
					// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
					const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;

					vResultingkMers[i].resize(iMaxRange);
					int32_t ikMerCounter = 0;

					// Frameshifting, so that only one amino acid has to be computed
					// Compute initial frames
					uint64_t sAAFrames[3] = { 0, 0, 0 };
					uint32_t aDeletekMerCounter[3] = { 0, 0, 0 };
					for (int32_t j = 0; j < iNumFrames; ++j) {
						string sTempFrame = _sMaxKBlank;
						dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
						sAAFrames[j] = aminoacidTokMer(sTempFrame);
						// if a character is 'illegal' then don't save the kMer containing it
						const auto& iPosOfU = sTempFrame.find_last_of('U');
						if (iPosOfU == string::npos) {
							vResultingkMers[i][ikMerCounter + j] = returnStuff(kMerVecOut, sAAFrames[j], iID);
						}
						else {
							vResultingkMers[i][ikMerCounter + j] = returnStuff(kMerVecOut, 0, iID);
							aDeletekMerCounter[j] = uint32_t(iPosOfU);
						}

					}

					ikMerCounter += iNumFrames;

					if (iMaxRange > 3) {

						// If the Dna has an irregular (not mod 3) tail, a special case must be handled and an additional step in the loop is necessary thus the !(iMaxRange % 3)
						// To map the 2 onto 1, ! is applied twice
						const int32_t& iMaxRangeMod3 = iMaxRange % 3;
						const int32_t& iMaxRangeMod3Neg = !(!iMaxRangeMod3);

						for (int32_t j = 1; 3 * (j + iMaxRangeMod3Neg) < iMaxRange; ++j) {
							for (int32_t k = 0; k < 3; ++k) {
								int8_t sTempAA = ' ';
								dnaToAminoacid(sDna, k + iMaxKTimes3 + 3 * (j - 1), sTempAA);
								sAAFrames[k] = aminoacidTokMer(sAAFrames[k], sTempAA);
								//string sDEBUG = kMerToAminoacid(sAAFrames[k], 12);
								if (sTempAA != 'U' && aDeletekMerCounter[k] == 0) {
									vResultingkMers[i][ikMerCounter + k] = returnStuff(kMerVecOut, sAAFrames[k], iID);
								}
								else {
									if (aDeletekMerCounter[k] != 0 && sTempAA != 'U') {
										--aDeletekMerCounter[k];
									}
									else {
										aDeletekMerCounter[k] = _iHighestK - 1;
									}
									vResultingkMers[i][ikMerCounter + k] = returnStuff(kMerVecOut, 0, iID);
								}
							}
							ikMerCounter += 3;
						}

						// compute frames of said tail
						for (int32_t j = 0; j < iMaxRangeMod3; ++j) {
							int8_t sTempAA = ' ';
							dnaToAminoacid(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
							sAAFrames[j] = aminoacidTokMer(sAAFrames[j], sTempAA);
							if (sTempAA != 'U' && aDeletekMerCounter[j] == 0) {
								vResultingkMers[i][ikMerCounter + j] = returnStuff(kMerVecOut, sAAFrames[j], iID);
							}
							else {
								if (aDeletekMerCounter[j] != 0 && sTempAA != 'U') {
									--aDeletekMerCounter[j];
								}
								else {
									aDeletekMerCounter[j] = _iHighestK - 1;
								}
								vResultingkMers[i][ikMerCounter + j] = returnStuff(kMerVecOut, 0, iID);
							}
						}
					}
				}
			}

			uint64_t iResultingSize = 0;
			for (int32_t i = 0; i < iActualLines; ++i) {
				iResultingSize += vResultingkMers[i].size();
			}

			// Write resulting kMers in stxxl vector but exclude those which are not useful
			double dStepSize = (fShrinkPercentage > 0.f) ? 100. / fShrinkPercentage : 0.;
			double dNextThrowOut = dStepSize;
			uint64_t iCounterOfThrowOut = 1;

			uint32_t iSizeCounterOfVec = 0;
			uint64_t iStart = kMerVecOut->size();
			if (fShrinkPercentage > 0.f) {
				kMerVecOut->resize(iStart + iResultingSize);
				for (const auto& kMerResVec : vResultingkMers) {
					for (const auto& element : kMerResVec) {
						if (isNotZero(kMerVecOut, element)) {
							if (iCounterOfThrowOut != static_cast<uint64_t>(dNextThrowOut)) {
								kMerVecOut->at(iStart + iSizeCounterOfVec++) = element;
							}
							else {
								dNextThrowOut += dStepSize;
							}
							++iCounterOfThrowOut;
						}
					}
				}
			}
			else {
				kMerVecOut->resize(iStart + iResultingSize);
				for (const auto& kMerResVec : vResultingkMers) {
					for (const auto& element : kMerResVec) {
						if (isNotZero(kMerVecOut, element)) {
							kMerVecOut->at(iStart + iSizeCounterOfVec++) = element;
						}
					}
				}
			}
			kMerVecOut->resize(iStart + iSizeCounterOfVec);
		}

		/////////////////////////////////////////////////
		// Not really necessary but whatever :D
		template <typename T>
		struct AAFrames {
			T first = 0, second = 0, third = 0; 
			T& operator[](const int32_t& i) {
				switch (i) {
				case 0:
					return first;
				case 1:
					return second;
				case 2:
					return third;
				}
				return first;
			}
		};
		/////////////////////////////////////////////////

		inline void convert_dnaTokMer(const string& sDna, const readIDType& iID, const int32_t& iMaxKTimes3, const Trie& T, unordered_map<uint64_t, Utilities::rangeContainer>& resultsVec, int64_t& iCurrentkMerCounter) {
			// This value gives the remaining length of the read which is the number of kMers created
			const int32_t& iMaxRange = int32_t(sDna.length() - iMaxKTimes3 + 1);
			if (iMaxRange >= 1) {
				// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
				const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;

				/*if (iSizeOfResultsBefore + iCurrentkMerCounter + iMaxRange >= iSoftMaxSize) {
				resultsVec[iProcIt].reserve(iSizeOfResultsBefore + iMaxRange * ((iProcIt + 1)*iNumOfResultsPerThread - i));
				}*/
				//resultsVec[iProcIt].resize(iSizeOfResultsBefore + iCurrentkMerCounter + iMaxRange);
				int32_t ikMerCounter = 0;

				// Frameshifting, so that only one amino acid has to be computed
				// Compute initial frames
				//uint64_t sAAFrames[3] = { 0, 0, 0 };
				//uint32_t aDeletekMerCounter[3] = { 0, 0, 0 };
				AAFrames<uint64_t> sAAFrames;
				//AAFrames<uint32_t> aDeletekMerCounter;
				for (int32_t j = 0; j < iNumFrames; ++j) {
					string sTempFrame = _sMaxKBlank;
					dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
					sAAFrames[j] = aminoacidTokMer(sTempFrame);

					//cout << kMerToAminoacid(sAAFrames[j], 12) << endl;
					const auto& range = T.GetIndexRange(sAAFrames[j] >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK));
					if (get<0>(range) != numeric_limits<uint64_t>::max()) {
						// if a character is 'illegal' then don't save the kMer containing it
						//const auto& iPosOfU = sTempFrame.find_last_of('U');
						//if (iPosOfU == string::npos) {
							auto& tempEntry = resultsVec[get<0>(range)];
							tempEntry.range = get<1>(range);
							//cout << kMerToAminoacid(sAAFrames[j], 12) << " " << kMerToAminoacid(static_cast<uint32_t>(sAAFrames[j] & 1073741823ULL), 12) << endl;
							if (_iMinK <= 6) {
								tempEntry.kMers_ST6.push_back(make_pair(sAAFrames[j], static_cast<uint32_t>(iID)));
							}
							else {
								tempEntry.kMers_GT6.push_back(make_pair(static_cast<uint32_t>(sAAFrames[j] & 1073741823ULL), static_cast<uint32_t>(iID)));
							}
							++iCurrentkMerCounter;
						//}
						//else {
							//aDeletekMerCounter[j] = uint32_t(iPosOfU);
						//}
					}
				}

				ikMerCounter += iNumFrames;

				if (iMaxRange > 3) {

					// If the Dna has an irregular (not mod 3) tail, a special case must be handled and an additional step in the loop is necessary thus the !(iMaxRange % 3)
					// To map the 2 onto 1, ! is applied twice
					const int32_t& iMaxRangeMod3 = iMaxRange % 3;
					const int32_t& iMaxRangeMod3Neg = !(!iMaxRangeMod3);

					for (int32_t j = 1; 3 * (j + iMaxRangeMod3Neg) < iMaxRange; ++j) {
						for (int32_t k = 0; k < 3; ++k) {
							int8_t sTempAA = ' ';
							dnaToAminoacid(sDna, k + iMaxKTimes3 + 3 * (j - 1), sTempAA);
							sAAFrames[k] = aminoacidTokMer(sAAFrames[k], sTempAA);

							//cout << kMerToAminoacid(sAAFrames[k], 12) << endl;

							const auto& range = T.GetIndexRange(sAAFrames[k] >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK));
							if (get<0>(range) != numeric_limits<uint64_t>::max()) {
								//string sDEBUG = kMerToAminoacid(sAAFrames[k], 12);
								//if (sTempAA != 'U' && aDeletekMerCounter[k] == 0) {
									auto& tempEntry = resultsVec[get<0>(range)];
									tempEntry.range = get<1>(range);
									//cout << kMerToAminoacid(static_cast<uint32_t>(sAAFrames[k] & 1073741823ULL), 12) << endl;
									if (_iMinK <= 6) {
										tempEntry.kMers_ST6.push_back(make_pair(sAAFrames[k], static_cast<uint32_t>(iID)));
									}
									else {
										tempEntry.kMers_GT6.push_back(make_pair(static_cast<uint32_t>(sAAFrames[k] & 1073741823ULL), static_cast<uint32_t>(iID)));
									}
									++iCurrentkMerCounter;
								//}
								//else {
								//	if (aDeletekMerCounter[k] != 0 && sTempAA != 'U') {
								//		--aDeletekMerCounter[k];
								//	}
								//	else {
								//		aDeletekMerCounter[k] = _iHighestK - 1;
								//	}
								//}
							}
						}
						ikMerCounter += 3;
					}

					// compute frames of said tail
					for (int32_t j = 0; j < iMaxRangeMod3; ++j) {
						int8_t sTempAA = ' ';
						dnaToAminoacid(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
						sAAFrames[j] = aminoacidTokMer(sAAFrames[j], sTempAA);

						//cout << kMerToAminoacid(sAAFrames[j], 12) << endl;

						const auto& range = T.GetIndexRange(sAAFrames[j] >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK));
						if (get<0>(range) != numeric_limits<uint64_t>::max()) {
							//if (sTempAA != 'U' && aDeletekMerCounter[j] == 0) {
								auto& tempEntry = resultsVec[get<0>(range)];
								tempEntry.range = get<1>(range);
								//cout << kMerToAminoacid(static_cast<uint32_t>(sAAFrames[j] & 1073741823ULL), 12) << endl;
								if (_iMinK <= 6) {
									tempEntry.kMers_ST6.push_back(make_pair(sAAFrames[j], static_cast<uint32_t>(iID)));
								}
								else {
									tempEntry.kMers_GT6.push_back(make_pair(static_cast<uint32_t>(sAAFrames[j] & 1073741823ULL), static_cast<uint32_t>(iID)));
								}
								++iCurrentkMerCounter;
							//}
							//else {
							//	if (aDeletekMerCounter[j] != 0 && sTempAA != 'U') {
							//		--aDeletekMerCounter[j];
							//	}
							//	else {
							//		aDeletekMerCounter[j] = _iHighestK - 1;
							//	}
							//}
						}
					}
				}
			}
		}

		inline void convert_alreadyTranslatedTokMers(const string& sProteinSequence, const readIDType& iReadID, const Trie& T, unordered_map<uint64_t, Utilities::rangeContainer>& resultsVec, int64_t& iCurrentkMerCounter) {
			const size_t& iLengthOfSequence = sProteinSequence.size();
			const int32_t& iMaxRange = int32_t(iLengthOfSequence - _iHighestK + 1);
			if (iMaxRange > 0) {
				//resultsVec[iProcIt].resize(iSizeOfResultsBefore + iCurrentkMerCounter + iMaxRange);
				for (; iCurrentkMerCounter < iMaxRange; ++iCurrentkMerCounter) {
					//resultsVec[iProcIt][iSizeOfResultsBefore + iCurrentkMerCounter] = make_tuple(aminoacidTokMer(sProteinSequence.cbegin() + iCurrentkMerCounter, sProteinSequence.cbegin() + iCurrentkMerCounter + 12), iReadID);
					
					const auto& kMer = aminoacidTokMer(sProteinSequence.cbegin() + iCurrentkMerCounter, sProteinSequence.cbegin() + iCurrentkMerCounter + 12);
					//cout << kMerToAminoacid(kMer, 12) << endl;
					const auto& range = T.GetIndexRange(kMer >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK));
					if (get<0>(range) != numeric_limits<uint64_t>::max()) {
						
						auto& tempEntry = resultsVec[get<0>(range)];
						tempEntry.range = get<1>(range);
						if (_iMinK <= 6) {
							tempEntry.kMers_ST6.push_back(make_pair(kMer, static_cast<uint32_t>(iReadID)));
						}
						else {
							tempEntry.kMers_GT6.push_back(make_pair(static_cast<uint32_t>(kMer & 1073741823ULL), static_cast<uint32_t>(iReadID)));
						}
						
					}
				}
			}
		}

		inline uint64_t convert_(unordered_map<uint64_t, Utilities::rangeContainer>& resultsVec, const vector<pair<string, readIDType>>& vLines, const int64_t& iActualLines, const Trie& T) {
			//int64_t iNumOfResultsPerThread = (iActualLines == 1) ? 0 : iActualLines / _iNumOfThreads;
			//int32_t iNumOfRemainingResults = (iActualLines == 1) ? 1 : iActualLines % _iNumOfThreads;

			const int32_t& iMaxKTimes3 = 3 * _iHighestK;
			int64_t iCurrentkMerCounter = 0;
			//uint64_t iSizeOfResultsBefore = resultsVec[iProcIt].size();

			/*int64_t iReserveSize = 0;

			for (int64_t i = 0; i < iActualLines; ++i) {
				int32_t iMaxRange = 0;
				if (_bInputAreAAs) {
					iMaxRange = int32_t(vLines[i].first.length() - _iHighestK + 1);
				}
				else {
					iMaxRange = int32_t(vLines[i].first.length() - iMaxKTimes3 + 1);
				}
				iReserveSize += iMaxRange;
			}

			if (iReserveSize > 0 && (iSizeOfResultsBefore + iCurrentkMerCounter + iReserveSize) > resultsVec[iProcIt].capacity()) {
				resultsVec[iProcIt].reserve(iSizeOfResultsBefore + iCurrentkMerCounter + iReserveSize);
			}*/

			//cout << iProcIt << " " << resultsVec[iProcIt].capacity() << " " << resultsVec[iProcIt].capacity() - resultsVec[iProcIt].size() << endl;

//#pragma omp parallel for
			//for (; iProcIt < _iNumOfThreads;) {
				
				//for (int64_t i = iProcIt * (iNumOfResultsPerThread); i < (iProcIt + 1)*iNumOfResultsPerThread; ++i) {
				for (int64_t i = 0; i < iActualLines; ++i) {
					if (vLines[i].first != "") {
						if (_bInputAreAAs) {
							convert_alreadyTranslatedTokMers(vLines[i].first, vLines[i].second, T, resultsVec, iCurrentkMerCounter);
						}
						else {
							convert_dnaTokMer(vLines[i].first, vLines[i].second, iMaxKTimes3, T, resultsVec, iCurrentkMerCounter);
						}
					}
				}

				// last part that could not be equally divided amongst the cpus
				//iSizeOfResultsBefore += iCurrentkMerCounter;
				//iCurrentkMerCounter = 0;
				/*if (iProcIt < iNumOfRemainingResults) {
					const string& sDna = vLines[iActualLines - iNumOfRemainingResults + iProcIt].first;
					const uint64_t& iID = vLines[iActualLines - iNumOfRemainingResults + iProcIt].second;
					convert_dnaTokMer(sDna, iID, iMaxKTimes3, resultsVec, iProcIt, iCurrentkMerCounter, iSizeOfResultsBefore);
					iSizeOfResultsBefore += iCurrentkMerCounter;
				}*/
				/*if (iSizeOfResultsBefore >= iSoftMaxSize) {
				resultsVec[iProcIt].reserve(iSizeOfResultsBefore);
				}*/
				//resultsVec[iProcIt].resize(iSizeOfResultsBefore);
			//}
				return iCurrentkMerCounter;
		}

		/*
		const string sSpacedMask = "XXOXXOOXOXXXXOOXXX";
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// convert lines from fasta/fastq to spaced-kMers and save them in a (stxxl-)vector
		template<typename vecType, typename idType> inline void convertLinesToSpacedkMers(const int32_t& iActualLines, const vector<string>& vLines, const idType& iID, const bool& bRC, unique_ptr<vecType>& kMerVecOut) {

			//vector<vector<tuple<uint64_t, uint32_t>>> vResultingkMers(iActualLines, vector<tuple<uint64_t, uint32_t>>());
			auto&& vResultingkMers = createStuff(kMerVecOut, iActualLines);
			const int32_t& iMaxKTimes3 = 3 * int32_t(sSpacedMask.size());

#pragma omp parallel for
			for (int32_t i = 0; i < iActualLines; ++i) {
				const string& sDna = vLines[i];
				const int32_t& iMaxRange = int32_t(sDna.length()) - iMaxKTimes3 + 1;
				if (iMaxRange >= 1) {
					// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
					const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;

					vResultingkMers[i].resize(iMaxRange);
					int32_t ikMerCounter = 0;

					// Frameshifting, so that only one amino acid has to be computed
					// Compute initial frames
					uint64_t sAAFrames[3] = { 0, 0, 0 };
					for (int32_t j = 0; j < iNumFrames; ++j) {
						string sTempFrame = sSpacedMask; // it'll be overwritten anyway, might as well use this
						dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
						string sNewTempFrame = "";
						if (bRC) {
							for (size_t l = 0, k = sSpacedMask.size() - 1; l < sSpacedMask.size(); ++l, --k) {
								if (sSpacedMask[k] == 'X') {
									sNewTempFrame += sTempFrame[l];
								}
							}
						}
						else {
							for (size_t l = 0; l < sSpacedMask.size(); ++l) {
								if (sSpacedMask[l] == 'X') {
									sNewTempFrame += sTempFrame[l];
								}
							}
						}
						sAAFrames[j] = aminoacidTokMer(sNewTempFrame);
						// if a character is 'illegal' then don't save the kMer containing it
						const auto& iPosOfU = sNewTempFrame.find_last_of('U');
						if (iPosOfU == string::npos) {
							vResultingkMers[i][ikMerCounter + j] = returnStuff(kMerVecOut, sAAFrames[j], iID);
						}
						else {
							vResultingkMers[i][ikMerCounter + j] = returnStuff(kMerVecOut, 0, iID);
						}

					}

					ikMerCounter += iNumFrames;

					if (iMaxRange > 3) {

						// If the Dna has an irregular (not mod 3) tail, a special case must be handled and an additional step in the loop is necessary thus the !(iMaxRange % 3)
						// To map the 2 onto 1, ! is applied twice
						const int32_t& iMaxRangeMod3 = iMaxRange % 3;
						const int32_t& iMaxRangeMod3Neg = !(!iMaxRangeMod3);

						for (int32_t j = 1; 3 * (j + iMaxRangeMod3Neg) < iMaxRange; ++j) {
							for (int32_t k = 0; k < 3; ++k) {
								string sTempFrame = sSpacedMask;
								dnaToAminoacid(sDna, iMaxKTimes3, k + 3 + 3 * (j - 1), &sTempFrame);
								string sNewTempFrame = "";
								if (bRC) {
									for (size_t l = 0, m = sSpacedMask.size() - 1; l < sSpacedMask.size(); ++l, --m) {
										if (sSpacedMask[m] == 'X') {
											sNewTempFrame += sTempFrame[l];
										}
									}
								}
								else {
									for (size_t l = 0; l < sSpacedMask.size(); ++l) {
										if (sSpacedMask[l] == 'X') {
											sNewTempFrame += sTempFrame[l];
										}
									}
								}
								sAAFrames[k] = aminoacidTokMer(sNewTempFrame);
								const auto& iPosOfU = sNewTempFrame.find_last_of('U');
								if (iPosOfU == string::npos) {
									vResultingkMers[i][ikMerCounter + k] = returnStuff(kMerVecOut, sAAFrames[k], iID);
								}
								else {
									vResultingkMers[i][ikMerCounter + k] = returnStuff(kMerVecOut, 0, iID);
								}
							}
							ikMerCounter += 3;
						}

						// compute frames of said tail
						for (int32_t j = 0; j < iMaxRangeMod3; ++j) {

							string sTempFrame = sSpacedMask;
							dnaToAminoacid(sDna, iMaxKTimes3, j + 3 + 3 * (iMaxRange / 3 - 1), &sTempFrame);
							string sNewTempFrame = "";
							if (bRC) {
								for (size_t m = 0, l = sSpacedMask.size() - 1; m < sSpacedMask.size(); ++m, --l) {
									if (sSpacedMask[l] == 'X') {
										sNewTempFrame += sTempFrame[m];
									}
								}
							}
							else {
								for (size_t m = 0; m < sSpacedMask.size(); ++m) {
									if (sSpacedMask[m] == 'X') {
										sNewTempFrame += sTempFrame[m];
									}
								}
							}
							sAAFrames[j] = aminoacidTokMer(sNewTempFrame);
							const auto& iPosOfU = sNewTempFrame.find_last_of('U');
							if (iPosOfU == string::npos) {
								vResultingkMers[i][ikMerCounter + j] = returnStuff(kMerVecOut, sAAFrames[j], iID);
							}
							else {
								vResultingkMers[i][ikMerCounter + j] = returnStuff(kMerVecOut, 0, iID);
							}
						}
					}
				}
			}

			uint64_t iResultingSize = 0;
			for (int32_t i = 0; i < iActualLines; ++i) {
				iResultingSize += vResultingkMers[i].size();
			}

			// Write resulting kMers in stxxl vector but exclude those which are not useful
			uint32_t iSizeCounterOfVec = 0;
			uint64_t iStart = kMerVecOut->size();
			kMerVecOut->resize(iStart + iResultingSize);
			for (const auto& kMerResVec : vResultingkMers) {
				for (const auto& element : kMerResVec) {
					if (isNotZero(kMerVecOut, element)) {
						kMerVecOut->at(iStart + iSizeCounterOfVec++) = element;
					}
				}
			}
			kMerVecOut->resize(iStart + iSizeCounterOfVec);
		}

		inline void convert_spaced_dnaToSKmer(const string& sDna, const readIDType& iID, const int32_t& iMaxKTimes3, unique_ptr<vector<tuple<uint64_t, readIDType>>[]>& resultsVec, const int32_t& iProcIt, int64_t& iCurrentkMerCounter, uint64_t& iSizeOfResultsBefore, const bool& bRC) {
			// This value gives the remaining length of the read which is the number of kMers created
			const int32_t& iMaxRange = int32_t(sDna.length() - iMaxKTimes3 + 1);
			if (iMaxRange >= 1) {
				// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
				const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;

				//if (iSizeOfResultsBefore + iCurrentkMerCounter + iMaxRange >= iSoftMaxSize) {
				//resultsVec[iProcIt].reserve(iSizeOfResultsBefore + iMaxRange * ((iProcIt + 1)*iNumOfResultsPerThread - i));
				//}
				resultsVec[iProcIt].resize(iSizeOfResultsBefore + iCurrentkMerCounter + iMaxRange);
				int32_t ikMerCounter = 0;

				// Frameshifting, so that only one amino acid has to be computed
				// Compute initial frames
				uint64_t sAAFrames[3] = { 0, 0, 0 };
				for (int32_t j = 0; j < iNumFrames; ++j) {
					string sTempFrame = sSpacedMask; // it'll be overwritten anyway, might as well use this
					dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
					string sNewTempFrame = "";
					if (bRC) {
						for (size_t i = 0, k = sSpacedMask.size() - 1; i < sSpacedMask.size(); ++i, --k) {
							if (sSpacedMask[k] == 'X') {
								sNewTempFrame += sTempFrame[i];
							}
						}
					}
					else {
						for (size_t i = 0; i < sSpacedMask.size(); ++i) {
							if (sSpacedMask[i] == 'X') {
								sNewTempFrame += sTempFrame[i];
							}
						}
					}
					sAAFrames[j] = aminoacidTokMer(sNewTempFrame);
					// if a character is 'illegal' then don't save the kMer containing it
					const auto& iPosOfU = sNewTempFrame.find_last_of('U');
					if (iPosOfU == string::npos) {
						resultsVec[iProcIt][iSizeOfResultsBefore + iCurrentkMerCounter++] = make_tuple(sAAFrames[j], iID);
					}

				}

				ikMerCounter += iNumFrames;

				if (iMaxRange > 3) {

					// If the Dna has an irregular (not mod 3) tail, a special case must be handled and an additional step in the loop is necessary thus the !(iMaxRange % 3)
					// To map the 2 onto 1, ! is applied twice
					const int32_t& iMaxRangeMod3 = iMaxRange % 3;
					const int32_t& iMaxRangeMod3Neg = !(!iMaxRangeMod3);

					for (int32_t j = 1; 3 * (j + iMaxRangeMod3Neg) < iMaxRange; ++j) {
						for (int32_t k = 0; k < 3; ++k) {
							string sTempFrame = sSpacedMask;
							dnaToAminoacid(sDna, iMaxKTimes3, k + 3 + 3 * (j - 1), &sTempFrame);
							string sNewTempFrame = "";
							if (bRC) {
								for (size_t i = 0, l = sSpacedMask.size() - 1; i < sSpacedMask.size(); ++i, --l) {
									if (sSpacedMask[l] == 'X') {
										sNewTempFrame += sTempFrame[i];
									}
								}
							}
							else {
								for (size_t i = 0; i < sSpacedMask.size(); ++i) {
									if (sSpacedMask[i] == 'X') {
										sNewTempFrame += sTempFrame[i];
									}
								}
							}
							sAAFrames[k] = aminoacidTokMer(sNewTempFrame);
							const auto& iPosOfU = sNewTempFrame.find_last_of('U');
							if (iPosOfU == string::npos) {
								resultsVec[iProcIt][iSizeOfResultsBefore + iCurrentkMerCounter++] = make_tuple(sAAFrames[k], iID);
							}
						}
						ikMerCounter += 3;
					}

					// compute frames of said tail
					for (int32_t j = 0; j < iMaxRangeMod3; ++j) {
						string sTempFrame = sSpacedMask;
						dnaToAminoacid(sDna, iMaxKTimes3, j + 3 + 3 * (iMaxRange / 3 - 1), &sTempFrame);
						string sNewTempFrame = "";
						if (bRC) {
							for (size_t i = 0, k = sSpacedMask.size() - 1; i < sSpacedMask.size(); ++i, --k) {
								if (sSpacedMask[k] == 'X') {
									sNewTempFrame += sTempFrame[i];
								}
							}
						}
						else {
							for (size_t i = 0; i < sSpacedMask.size(); ++i) {
								if (sSpacedMask[i] == 'X') {
									sNewTempFrame += sTempFrame[i];
								}
							}
						}
						sAAFrames[j] = aminoacidTokMer(sNewTempFrame);
						const auto& iPosOfU = sNewTempFrame.find_last_of('U');
						if (iPosOfU == string::npos) {
							resultsVec[iProcIt][iSizeOfResultsBefore + iCurrentkMerCounter++] = make_tuple(sAAFrames[j], iID);
						}
					}
				}
			}
		}

		inline void convert_spaced(unique_ptr<vector<tuple<uint64_t, readIDType>>[]>& resultsVec, const vector<pair<string, readIDType>>& vLines, const int64_t& iActualLines, const bool& bRC) {
			int64_t iNumOfResultsPerThread = (iActualLines == 1) ? 0 : iActualLines / _iNumOfThreads;
			int32_t iNumOfRemainingResults = (iActualLines == 1) ? 1 : iActualLines % _iNumOfThreads;

			const int32_t& iMaxKTimes3 = 3 * int32_t(sSpacedMask.size());

#pragma omp parallel for
			for (int32_t iProcIt = 0; iProcIt < _iNumOfThreads; ++iProcIt) {
				int64_t iCurrentkMerCounter = 0;
				uint64_t iSizeOfResultsBefore = resultsVec[iProcIt].size();
				for (int64_t i = iProcIt * (iNumOfResultsPerThread); i < (iProcIt + 1)*iNumOfResultsPerThread; ++i) {
					convert_spaced_dnaToSKmer(vLines[i].first, vLines[i].second, iMaxKTimes3, resultsVec, iProcIt, iCurrentkMerCounter, iSizeOfResultsBefore, bRC);
				}
				// last part that could not be equally divided amongst the cpus
				iSizeOfResultsBefore += iCurrentkMerCounter;
				iCurrentkMerCounter = 0;
				if (iProcIt < iNumOfRemainingResults) {
					const string& sDna = vLines[iActualLines - iNumOfRemainingResults + iProcIt].first;
					const readIDType& iID = vLines[iActualLines - iNumOfRemainingResults + iProcIt].second;
					convert_spaced_dnaToSKmer(sDna, iID, iMaxKTimes3, resultsVec, iProcIt, iCurrentkMerCounter, iSizeOfResultsBefore, bRC);
					iSizeOfResultsBefore += iCurrentkMerCounter;
				}
				//if (iSizeOfResultsBefore >= iSoftMaxSize) {
				//resultsVec[iProcIt].reserve(iSizeOfResultsBefore);
				//}
				resultsVec[iProcIt].resize(iSizeOfResultsBefore);
			}
		}
		*/

		inline uint64_t convertLinesTokMers_(const int32_t& iActualLines, const vector<pair<string, readIDType>>& vLines, const vector<pair<string, readIDType>>& vRCLines, unordered_map<uint64_t, Utilities::rangeContainer>& kMerVecOut, const Trie& T, const bool& bSpaced) {

			/*for (int32_t i = 0; i < _iNumOfThreads; ++i) {
				kMerVecOut[i].reserve(static_cast<uint64_t>(iSoftMaxSize));
			}*/

			//const uint64_t& iSoftSizeWithCopyBuffer = static_cast<uint64_t>(iSoftMaxSize);
			/*
			uint64_t iSumBefore = 0;
			for (int32_t i = 0; i < _iNumOfThreads; ++i) {
				iSumBefore += kMerVecOut[i].size();
			}*/

			uint64_t iNumOfkMers = 0;

			if (bSpaced) {
				//convert_spaced(kMerVecOut, vLines, iActualLines, false);
				//convert_spaced(kMerVecOut, vRCLines, iActualLines, true);
			}
			else {
				iNumOfkMers += convert_(kMerVecOut, vLines, iActualLines, T);
				if (!_bInputAreAAs) {
					iNumOfkMers += convert_(kMerVecOut, vRCLines, iActualLines, T);
				}
			}
			/*
			uint64_t outVal = 0;
			for (int32_t i = 0; i < _iNumOfThreads; ++i) {
				outVal += kMerVecOut[i].size();
			}*/

			return iNumOfkMers;//outVal - iSumBefore;
		}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public:

		void BuildAll(const string& fContentFile, const string& sDirectory, const string& fOutFile, const uint64_t& iMem, const float& fShrinkPercentage = 0.f) {
			try {
				// test if files exists
				if (!ifstream(fContentFile)) {
					throw runtime_error("Content file not found.");
				}

				/*unique_ptr<uint64_t[]> vDivisionArray(new uint64_t[_iNumOfThreads - 1]);
				const uint64_t& iStep = 0xF7BDEF7BDEF7BDE / _iNumOfThreads;
				for (int i = 0; i < _iNumOfThreads - 1; ++i) {
					vDivisionArray[i] = (i + 1) * iStep;
				}*/



				ofstream derp(fOutFile);
				if (derp.fail()) {
					throw runtime_error("File couldn't be created, maybe a wrong path was used?");
				}
				derp.close();
				//stxxlFile* stxxlOutFile = new stxxlFile(fOutFile, stxxl::file::RDWR);
				//contentVecType_32p* vOutVec = new contentVecType_32p(stxxlOutFile, 0);

				uint32_t iIdxCounter = 1;
				unordered_map<uint32_t, uint32_t> mIDsAsIdx; mIDsAsIdx[0] = 0;
				unordered_map<uint32_t, string> mIdxToName; mIdxToName[0] = "non_unique";
				unordered_map<string, uint32_t> mAccToID;
				ifstream content(fContentFile);
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

				Build brick(_sTemporaryPath, 0, iMem / (sizeof(packedBigPair)), iIdxCounter, mIDsAsIdx);
				size_t overallCharsRead = 0;
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
						readFasta(fastaFile, mAccToID, brick, iFileLength, overallCharsRead, filesAndSize.second, fShrinkPercentage);
					}
				}
				else {
					ifstream fastaFile;
					fastaFile.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
					fastaFile.open(sDirectory);
					if (_bVerbose) {
						cout << "OUT: Current file: " << sDirectory << endl;
					}
					fastaFile.seekg(0, fastaFile.end);
					const uint64_t& iFileLength = fastaFile.tellg();
					fastaFile.seekg(0, fastaFile.beg);
					readFasta(fastaFile, mAccToID, brick, iFileLength, overallCharsRead, iFileLength, fShrinkPercentage);
				}

				// Finalize
				brick.IntToExt2(true);
				const uint64_t& iSizeOfFinalIndex = brick.createDatabase(fOutFile, iIdxCounter, mIdxToName);

				// Create Trie out of final file
				stxxlFile libFile(fOutFile, stxxlFile::RDONLY);
				const contentVecType_32p libVec(&libFile, iSizeOfFinalIndex);
				Trie T(static_cast<int8_t>(12), static_cast<int8_t>(_iMinK), 6);
				T.SaveToStxxlVec(&libVec, fOutFile);
			}
			catch (invalid_argument&) {
				cerr << "ERROR: content file has not the right format. Taxids in the second row are required to be integers!" << endl;
				return;
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Create a library out of many fasta files given in a content file made by the python script to generate a database
		void CreateAll(const string& fContentFile, const string& sDirectory, const string& fOutFile, const bool& bSpaced, const uint64_t& iMemory, const float& fShrinkPercentage = 0.f) {
			// test if files exists
			if (!ifstream(fContentFile)) {
				throw runtime_error("The content file does not exist");
			}


			ofstream derp(fOutFile);
			derp.close();
			unique_ptr<stxxlFile> stxxlOutFile(new stxxlFile(fOutFile, stxxl::file::RDWR));
			unique_ptr<contentVecType_32p> vOutVec(new contentVecType_32p(stxxlOutFile.get(), 0));

			unordered_map<string, uint32_t> mAccToID;
			ifstream content(fContentFile);
			string sDummy = "";
			while (content.good()) {
				getline(content, sDummy);
				if (sDummy != "") {
					const auto& line = Utilities::split(sDummy, '\t');
					if (line.size() == 4) {
						const auto& vAccessionNumbers = Utilities::split(line[3], ';');
						for (const auto& acc : vAccessionNumbers) {
							mAccToID.insert(make_pair(acc, stoul(line[1])));
						}
					}
				}
			}

#if _WIN32 || _WIN64
			for (auto& fsPath : std::experimental::filesystem::directory_iterator(sDirectory)) {
				ifstream fastaFile(fsPath.path());
				readFasta(fastaFile, vOutVec, fShrinkPercentage, bSpaced, mAccToID);
			}
#else
#if __GNUC__
			DIR           *dirp;
			struct dirent *directory;

			dirp = opendir(sDirectory.c_str());
			if (dirp) {
				while ((directory = readdir(dirp)) != NULL)
				{
					const string fName(directory->d_name);
					if (fName != "." && fName != "..") {
						ifstream fastaFile(sDirectory + fName);
						readFasta(fastaFile, vOutVec, fShrinkPercentage, bSpaced, mAccToID);
					}
				}

				closedir(dirp);
			}

#else
		#error NO SUITABLE COMPILER FOUND
#endif
#endif



			//ifstream content(fContentFile);
			//string sDummy = "";
			//while (content.good()) {
			//	getline(content, sDummy);
			//	if (sDummy != "") {
			//		const auto& line = Utilities::split(sDummy, '\t');
			//		if (line.size() == 4) {
			//			ifstream fastaFile(line[3]);
			//			readFasta(fastaFile, iID, iMemory, vOutVec, false, NULL, true);
			//			++iID;
			//		}
			//	}
			//}

			stxxl::sort(vOutVec->begin(), vOutVec->end(), SCompareStructForSTXXLSort(), iMemory);

			/*auto funcForUnique = [](tuple<uint64_t, uint32_t>& a, tuple<uint64_t, uint32_t>& b) {
			if (get<0>(a) == get<0>(b)) {
			if (get<1>(a) != get<1>(b)) {
			get<1>(a) = 0;
			}
			return true;
			}
			return false;
			};
			contentVecType_32p::iterator newEnd = std::unique(vOutVec->begin(), vOutVec->end(), funcForUnique);
			vOutVec->resize(newEnd - vOutVec->begin(), true);
			*/

			auto funcForUnique = [](packedBigPair& a, packedBigPair& b) {
				if (a == b) {
					return true;
				}
				return false;
			};
			contentVecType_32p::iterator newEnd = std::unique(vOutVec->begin(), vOutVec->end(), funcForUnique);
			vOutVec->resize(newEnd - vOutVec->begin(), true);

			/*ofstream derpfile("F:/try/derpfile2.txt");
			for (const auto& elem : *vOutVec) {
			derpfile << get<0>(elem) << ", " << get<1>(elem) << endl;
			}*/

			ofstream fLibInfo(fOutFile + "_info.txt");
			fLibInfo << vOutVec->size();

			vOutVec->export_files("_");
		}
	};
}