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
		const bool bUnfunny;

	public:
		Read(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const bool& bTranslated = false, const string& stxxl_mode = "", const bool& bUnfunny = false, const bool& bSixFrames = false) : kASA(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, stxxl_mode, bSixFrames), _bInputAreAAs(bTranslated), bUnfunny(bUnfunny) {}

	protected:

		// test for garbage to correctly count matchable k-mers
		const uint64_t _tails[12] = { 0x1E, 0x3C0, 0x7800, 0xF0000, 0x1E00000, 0x3C000000, 0x780000000, 0xF00000000, 0x1E000000000, 0x3C0000000000, 0x7800000000000, 0xF000000000000 }; // ^, ^@, ^@@, ...

		inline void convert_dnaTokMer(const string& sDna, const readIDType& iID, const int32_t& iMaxKTimes3, const Trie&, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& resultsVec, int64_t& iCurrentkMerCounter, vector<uint64_t>& vCountGarbagekMerPerK) {
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

					auto kmerForSearchInTrie = sAAFrames[j];
					if (bUnfunny) {
						/*uint64_t tVal = 0;
						for (int32_t iLeftStart = 55, iShifted = 0; iLeftStart >= 5; iLeftStart -= 10, iShifted += 5) {
							tVal |= (sAAFrames[j] & (31ULL << iLeftStart)) << iShifted;
						}
						sAAFrames[j] = tVal;*/
						kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
					}


					//cout << kMerToAminoacid(sAAFrames[j], 12) << endl;uint64_t start = 0, range = 0;
					//uint64_t start = 0;
					//uint32_t range = 0;
					//T.GetIndexRange(kmerForSearchInTrie >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
					//if (start != numeric_limits<uint64_t>::max()) {

						for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
							vCountGarbagekMerPerK[kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal];
						}
						
						resultsVec.push_back(make_tuple(0, kmerForSearchInTrie, 0, static_cast<uint32_t>(iID)));

						++iCurrentkMerCounter;
					//}
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

							auto kmerForSearchInTrie = sAAFrames[k];
							if (bUnfunny) {
								/*uint64_t tVal = 0;
								for (int32_t iLeftStart = 55, iShifted = 0; iLeftStart >= 5; iLeftStart -= 10, iShifted += 5) {
									tVal |= (sAAFrames[k] & (31ULL << iLeftStart)) << iShifted;
								}
								sAAFrames[k] = tVal;*/
								kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
							}

							//uint64_t start = 0;
							//uint32_t range = 0;
							//T.GetIndexRange(kmerForSearchInTrie >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
							//if (start != numeric_limits<uint64_t>::max()) {

								for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
									vCountGarbagekMerPerK[kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal];
								}

								resultsVec.push_back(make_tuple(0, kmerForSearchInTrie, 0, static_cast<uint32_t>(iID)));
								++iCurrentkMerCounter;
							//}
						}
						ikMerCounter += 3;
					}

					// compute frames of said tail
					for (int32_t j = 0; j < iMaxRangeMod3; ++j) {
						int8_t sTempAA = ' ';
						dnaToAminoacid(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
						sAAFrames[j] = aminoacidTokMer(sAAFrames[j], sTempAA);

						//cout << kMerToAminoacid(sAAFrames[j], 12) << endl;
						auto kmerForSearchInTrie = sAAFrames[j];
						if (bUnfunny) {
							/*uint64_t tVal = 0;
							for (int32_t iLeftStart = 55, iShifted = 0; iLeftStart >= 5; iLeftStart -= 10, iShifted += 5) {
								tVal |= (sAAFrames[j] & (31ULL << iLeftStart)) << iShifted;
							}
							sAAFrames[j] = tVal;*/
							kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
						}

						//uint64_t start = 0;
						//uint32_t range = 0;
						//T.GetIndexRange(kmerForSearchInTrie >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
						//if (start != numeric_limits<uint64_t>::max()) {

							for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
								vCountGarbagekMerPerK[kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal];
							}

							resultsVec.push_back(make_tuple(0, kmerForSearchInTrie, 0, static_cast<uint32_t>(iID)));
							++iCurrentkMerCounter;
						//}
					}
				}
			}
		}

		inline void convert_alreadyTranslatedTokMers(const string& sProteinSequence, const readIDType& iReadID, const Trie&, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& resultsVec, int64_t& iCurrentkMerCounter) {
			const size_t& iLengthOfSequence = sProteinSequence.size();
			const int32_t& iMaxRange = int32_t(iLengthOfSequence - _iHighestK + 1);
			if (iMaxRange > 0) {
				//resultsVec[iProcIt].resize(iSizeOfResultsBefore + iCurrentkMerCounter + iMaxRange);
				for (; iCurrentkMerCounter < iMaxRange; ++iCurrentkMerCounter) {
					//resultsVec[iProcIt][iSizeOfResultsBefore + iCurrentkMerCounter] = make_tuple(aminoacidTokMer(sProteinSequence.cbegin() + iCurrentkMerCounter, sProteinSequence.cbegin() + iCurrentkMerCounter + 12), iReadID);

					auto kMer = aminoacidTokMer(sProteinSequence.cbegin() + iCurrentkMerCounter, sProteinSequence.cbegin() + iCurrentkMerCounter + 12);
					
					if (bUnfunny) {
						kMer = aminoAcidsToAminoAcid(kMer);;
					}
					
					//uint64_t start = 0;
					//uint32_t range = 0;
					//T.GetIndexRange(kMer >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
					//if (start != numeric_limits<uint64_t>::max()) {
						resultsVec.push_back(make_tuple(0, kMer, 0, iReadID));
					//}
				}
			}
		}

		inline uint64_t convertLinesTokMers_(const int32_t& iActualLines, const vector<pair<string, readIDType>>& vLines, const vector<pair<string, readIDType>>& vRCLines, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& kMerVecOut, const Trie& T, const bool&, vector<uint64_t>& vCountGarbagekMerPerK) {

			const int32_t& iMaxKTimes3 = 3 * _iHighestK;
			int64_t iNumOfkMers = 0;

			for (int64_t i = 0; i < iActualLines; ++i) {
				if (vLines[i].first != "" || vRCLines[i].first != "") {
					if (_bInputAreAAs) {
						convert_alreadyTranslatedTokMers(vLines[i].first, vLines[i].second, T, kMerVecOut, iNumOfkMers);
					}
					else {
						convert_dnaTokMer(vLines[i].first, vLines[i].second, iMaxKTimes3, T, kMerVecOut, iNumOfkMers, vCountGarbagekMerPerK);
						convert_dnaTokMer(vRCLines[i].first, vRCLines[i].second, iMaxKTimes3, T, kMerVecOut, iNumOfkMers, vCountGarbagekMerPerK);
					}
				}
			}

			return iNumOfkMers;
		}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
		inline uint64_t readFastq_partialSort(T& input, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const bool& bSpaced, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, unique_ptr<strTransfer>& transfer, const Trie& Tr, vector<uint64_t>& vCountGarbagekMerPerK) {
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
						pair<std::string, bool> resultChunkPair;
						Utilities::getChunk(input, move(resultChunkPair), move(iNumOfChars));

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
									vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));// + sFalsekMerMarker.length())));
									iSoftMaxSize -= sizeof(pair<string,uint32_t>) + sName.size() * sizeof(char) + sizeof(uint32_t);
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

					const auto& iNumOfkMers = convertLinesTokMers_(int32_t(iActualLines - 1), vLines, vRCLines, vOut, Tr, bSpaced, vCountGarbagekMerPerK);
					iSumOfkMers += iNumOfkMers;
					iProcID = (iProcID + 1) % _iNumOfThreads;
					//iSoftMaxSize -= iNumOfkMers * (sizeof(pair<uint32_t, uint32_t>) + sizeof(uint64_t));
					
					iSoftMaxSize -= iNumOfkMers * sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>);

					if (iSoftMaxSize <= 14399776 + 4 * iAmountOfSpecies) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 20 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
						bNotFull = false;
					}


					vLines[0] = vLines[iActualLines - 1];
					iActualLines = 1;
				}
				if (!input) {
					vRCLines[0] = make_pair("", 0);
					iSumOfkMers += convertLinesTokMers_(1, vLines, vRCLines, vOut, Tr, bSpaced, vCountGarbagekMerPerK);
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
		inline uint64_t readFasta_partialSort(T& input, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const bool& bSpaced, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, unique_ptr<strTransfer>& transfer, const Trie& Tr, vector<uint64_t>& vCountGarbagekMerPerK) {
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

				size_t iLinesPerRun = 2;//2;

				vector<pair<string, readIDType>> vLines(iLinesPerRun), vRCLines(iLinesPerRun);
				size_t iActualLines = 1;
				vLines[0] = transfer->lastLine;

				while (input && bNotFull) {
					while (input && iActualLines < iLinesPerRun) {
						uint64_t iNumOfChars = 0;
						pair<std::string, bool> resultChunkPair;
						Utilities::getChunk(input, move(resultChunkPair), move(iNumOfChars));

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
									vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));// + sFalsekMerMarker.length())));
									iSoftMaxSize -= sizeof(pair<string, uint32_t>) + sName.size() * sizeof(char) + sizeof(uint32_t);
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

					const auto& iNumOfkMers = convertLinesTokMers_(int32_t(iActualLines - 1), vLines, vRCLines, vOut, Tr, bSpaced, vCountGarbagekMerPerK);
					iSumOfkMers += iNumOfkMers;
					
					iSoftMaxSize -= iNumOfkMers * sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>);

					
					if (iSoftMaxSize <= 14399776 + 4 * iAmountOfSpecies) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 20 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
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
						vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));// + sFalsekMerMarker.length())));
						iSoftMaxSize -= sName.size() * sizeof(char) + sizeof(uint32_t);
						bAddTail = false;
					}
					vRCLines[0] = make_pair("", 0);
					iSumOfkMers += convertLinesTokMers_(1, vLines, vRCLines, vOut, Tr, bSpaced, vCountGarbagekMerPerK);
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

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// convert lines from fasta/fastq to kMers and save them in a (stxxl-)vector
		inline void dnaTokMers(const string& sDna, const uint32_t& iID, unique_ptr<contentVecType_32p>& kMerVecOut, Build& vBricks, const float& fShrinkPercentage) {
			vector<tuple<uint64_t, uint32_t>> vResultingkMers;
			const int32_t& iMaxKTimes3 = 3 * _iHighestK;

			const int32_t& iMaxRange = int32_t(sDna.length()) - iMaxKTimes3 + 1;
			if (iMaxRange >= 1) {
				// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
				const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;

				vResultingkMers.resize(iMaxRange);
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
						vResultingkMers[ikMerCounter + j] = make_tuple(sAAFrames[j], iID);
					}
					else {
						vResultingkMers[ikMerCounter + j] = make_tuple(0, iID);
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
								vResultingkMers[ikMerCounter + k] = make_tuple(sAAFrames[k], iID);
							}
							else {
								if (aDeletekMerCounter[k] != 0 && sTempAA != 'U') {
									--aDeletekMerCounter[k];
								}
								else {
									aDeletekMerCounter[k] = _iHighestK - 1;
								}
								vResultingkMers[ikMerCounter + k] = make_tuple(0, iID);
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
							vResultingkMers[ikMerCounter + j] = make_tuple(sAAFrames[j], iID);
						}
						else {
							if (aDeletekMerCounter[j] != 0 && sTempAA != 'U') {
								--aDeletekMerCounter[j];
							}
							else {
								aDeletekMerCounter[j] = _iHighestK - 1;
							}
							vResultingkMers[ikMerCounter + j] = make_tuple(0, iID);
						}
					}
				}
			}

			uint64_t iResultingSize = vResultingkMers.size();

			// Write resulting kMers in stxxl vector but exclude those which are not useful
			double dStepSize = (fShrinkPercentage > 0.f) ? 100. / fShrinkPercentage : 0.;
			double dNextThrowOut = dStepSize;
			uint64_t iCounterOfThrowOut = 1;

			uint32_t iSizeCounterOfVec = 0;
			uint64_t iStart = 0;
			if (kMerVecOut) {
				iStart = kMerVecOut->size();
			}
			if (fShrinkPercentage > 0.f) {
				if (kMerVecOut) {
					kMerVecOut->resize(iStart + iResultingSize);
				}
				for (auto& element : vResultingkMers) {
					if (get<0>(element) != 0) {
						if (iCounterOfThrowOut != static_cast<uint64_t>(dNextThrowOut)) {
							if (bUnfunny) {
								/*uint64_t tVal = 0;
								for (int32_t iLeftStart = 55, iShifted = 0; iLeftStart >= 5; iLeftStart -= 10, iShifted += 5) {
									tVal |= (get<0>(element) & (31ULL << iLeftStart)) << iShifted;
								}
								get<0>(element) = tVal;*/

								get<0>(element) = aminoAcidsToAminoAcid(get<0>(element));
							}

							if (kMerVecOut) {
								kMerVecOut->at(iStart + iSizeCounterOfVec++) = element;
							}
							else {
								if (!(vBricks.addToInt(element))) {
									vBricks.IntToExtPart();
								}
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
				if (kMerVecOut) {
					kMerVecOut->resize(iStart + iResultingSize);
				}
				for (auto& element : vResultingkMers) {
					if (get<0>(element) != 0) {
						if (bUnfunny) {
							/*uint64_t tVal = 0;
							for (int32_t iLeftStart = 55, iShifted = 0; iLeftStart >= 5; iLeftStart -= 10, iShifted += 5) {
								tVal |= (get<0>(element) & (31ULL << iLeftStart)) << iShifted;
							}
							get<0>(element) = tVal;*/
							
							get<0>(element) = aminoAcidsToAminoAcid(get<0>(element));
						}

						if (kMerVecOut) {
							kMerVecOut->at(iStart + iSizeCounterOfVec++) = element;
						}
						else {
							if (!(vBricks.addToInt(element))) {
								vBricks.IntToExtPart();
							}
						}
					}
				}
			}
			if (kMerVecOut) {
				kMerVecOut->resize(iStart + iSizeCounterOfVec);
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// convert protein sequences to k-mers in the building step
		inline void proteinTokMers(const string& sAASequence, const uint32_t& iIdx, unique_ptr<contentVecType_32p>& kMerVecOut, Build& vBricks, const float& fShrinkPercentage) {
			vector<tuple<uint64_t, uint32_t>> vResultingkMers;
			const int32_t& iMaxRange = int32_t(sAASequence.length()) - _iHighestK + 1;
			if (iMaxRange > 0) {
				vResultingkMers.resize(iMaxRange);
				for (int32_t i = 0; i < iMaxRange; ++i) {
					const auto& kMer = aminoacidTokMer(sAASequence.cbegin() + i, sAASequence.cbegin() + i + 12);
					vResultingkMers[i] = make_tuple(kMer, iIdx);
				}
			}

			// Write resulting kMers in stxxl vector but exclude those which are not useful
			double dStepSize = (fShrinkPercentage > 0.f) ? 100. / fShrinkPercentage : 0.;
			double dNextThrowOut = dStepSize;
			uint64_t iCounterOfThrowOut = 1;

			uint64_t iResultingSize = vResultingkMers.size();
			uint32_t iSizeCounterOfVec = 0;
			uint64_t iStart = 0;
			if (kMerVecOut) {
				iStart = kMerVecOut->size();
			}

			if (fShrinkPercentage > 0.f) {
				if (kMerVecOut) {
					kMerVecOut->resize(iStart + iResultingSize);
				}
				for (auto& element : vResultingkMers) {
					if (iCounterOfThrowOut != static_cast<uint64_t>(dNextThrowOut)) {
						if (bUnfunny) {
							get<0>(element) = aminoAcidsToAminoAcid(get<0>(element));
						}

						if (kMerVecOut) {
							kMerVecOut->at(iStart + iSizeCounterOfVec++) = element;
						}
						else {
							if (!(vBricks.addToInt(element))) {
								vBricks.IntToExtPart();
							}
						}
					}
					else {
						dNextThrowOut += dStepSize;
					}
					++iCounterOfThrowOut;
				}
			}
			else {
				if (kMerVecOut) {
					kMerVecOut->resize(iStart + iResultingSize);
				}
				for (auto& element : vResultingkMers) {
					if (bUnfunny) {
						get<0>(element) = aminoAcidsToAminoAcid(get<0>(element));
					}

					if (kMerVecOut) {
						kMerVecOut->at(iStart + iSizeCounterOfVec++) = element;
					}
					else {
						if (!(vBricks.addToInt(element))) {
							vBricks.IntToExtPart();
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fasta and create a kMer-Vec, used in BuildAll(...)
		inline void readFasta(ifstream& input, const unordered_map<string, uint32_t>& mAccToID, Build& vBricks, unique_ptr<contentVecType_32p>& vOut, const uint64_t& iFileLength, size_t& overallCharsRead, const size_t& overallFilesSize, const float& fShrinkPercentage) {
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



				bool bAddFalseMarker = false, bNextIdx = false, bAddFalseRCMarker = true;

				const int32_t& iConcurrentLines = 20;
				string sOverhang = "", sDNA = "";

				uint64_t iNumOfCharsRead = 0, iCurrentPercentage = 0, iCurrentOverallPercentage = 0;

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
								bAddFalseRCMarker = true;
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
							bAddFalseRCMarker = false;
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
					string sRCDNA = reverseComplement(sDNA);

					if (bAddFalseRCMarker && !sRCDNA.empty()) {
						sRCDNA += sFalsekMerMarker;
						bAddFalseRCMarker = false;
					}
					if (bAddFalseMarker && !sDNA.empty()) {
						sDNA += sFalsekMerMarker;
						bAddFalseMarker = false;
					}

					if (_bInputAreAAs) {
						proteinTokMers(sDNA, iIdx, vOut, vBricks, fShrinkPercentage);
					}
					else {
						dnaTokMers(sDNA, iIdx, vOut, vBricks, fShrinkPercentage);
						if (_bSixFrames) {
							dnaTokMers(sRCDNA, iIdx, vOut, vBricks, fShrinkPercentage);
						}
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



				Utilities::createFile(fOutFile);
				
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
						if (line.size() >= 4) {
							mIDsAsIdx[stoul(line[1])] = iIdxCounter;
							mIdxToName[iIdxCounter] = line[0];
							++iIdxCounter;
							const auto& vAccessionNumbers = Utilities::split(line[3], ';');
							for (const auto& acc : vAccessionNumbers) {
								mAccToID.insert(make_pair(acc, stoul(line[1])));
							}
						}
						else {
							throw runtime_error("Content file contains less than 4 columns, it may be damaged... The faulty line was: " + sDummy + "\n");
						}
					}
				}

				Build brick(_sTemporaryPath, _iNumOfCall, _iNumOfThreads, iMem / (sizeof(packedBigPair)), iIdxCounter, mIDsAsIdx);

				size_t overallCharsRead = 0;
				unique_ptr<contentVecType_32p> dummy;
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
						readFasta(fastaFile, mAccToID, brick, dummy, iFileLength, overallCharsRead, filesAndSize.second, fShrinkPercentage);
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
					readFasta(fastaFile, mAccToID, brick, dummy, iFileLength, overallCharsRead, iFileLength, fShrinkPercentage);
				}
				
				// Finalize
				brick.IntToExtPart();
				const uint64_t& iSizeOfFinalIndex = brick.mergeTemporaries(fOutFile);
				
				unique_ptr<uint64_t[]> arrFrequencies;
				arrFrequencies.reset(new uint64_t[iIdxCounter * 12]);
				for (uint64_t i = 0; i < iIdxCounter * 12; ++i) {
					arrFrequencies[i] = 0;
				}

				// Create Trie and frequencies out of final file
				stxxlFile* libFile = new stxxlFile(fOutFile, stxxlFile::RDONLY);
				const contentVecType_32p* libVec = new const contentVecType_32p(libFile, iSizeOfFinalIndex);
				Trie T(static_cast<int8_t>(12), static_cast<int8_t>(_iMinK), 6);
				T.SaveToStxxlVec(libVec, fOutFile, &arrFrequencies, mIDsAsIdx);

				// If taxaOnly index is desired, create it here:
				if (bUnfunny) {
					Utilities::createFile(fOutFile+"_taxOnly");
					stxxlFile* taxaOnlyFile = new stxxlFile(fOutFile+"_taxOnly", stxxlFile::RDWR);
					taxaOnly* taxaOnlyFileVec = new taxaOnly(taxaOnlyFile, iSizeOfFinalIndex);

					auto tOIt = taxaOnlyFileVec->begin();
					for (const auto& elem : *libVec) {
						*tOIt = static_cast<uint16_t>(mIDsAsIdx.find(elem.second)->second);
						++tOIt;
					}

					taxaOnlyFileVec->export_files("_");
					delete taxaOnlyFileVec;
					delete taxaOnlyFile;
					delete libVec;
					libFile->close_remove();
					delete libFile;

					remove(fOutFile.c_str());
					rename((fOutFile + "_taxOnly").c_str(), fOutFile.c_str());
				}

				ofstream outFile(fOutFile + "_f.txt");
				for (uint32_t j = 0; j < iIdxCounter; ++j) {
					outFile << Utilities::checkIfInMap(mIdxToName, j)->second << "\t";
					outFile << arrFrequencies[j * 12];
					for (int32_t k = 1; k < 12; ++k) {
						outFile << "\t" << arrFrequencies[j * 12 + k];
					}
					outFile << endl;
				}
			}
			catch (invalid_argument&) {
				cerr << "ERROR: content file has not the right format. Taxids in the second row are required to be integers!" << endl;
				return;
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
	};
}