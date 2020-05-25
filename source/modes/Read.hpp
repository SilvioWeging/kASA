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
#include "Build.hpp"
#include "Trie.hpp"
#include "../utils/WorkerThread.hpp"

namespace kASA {
	class Read : public kASA {
		const bool _bInputAreAAs;
		const bool bUnfunny;

	public:
		Read(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const bool& bTranslated = false, const string& stxxl_mode = "", const bool& bUnfunny = false, const bool& bSixFrames = false) : kASA(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, stxxl_mode, bSixFrames), _bInputAreAAs(bTranslated), bUnfunny(bUnfunny) {}

	protected:

		// test for garbage to correctly count matchable k-mers
		const uint64_t _tails[12] = { 0x1E, 0x3C0, 0x7800, 0xF0000, 0x1E00000, 0x3C000000, 0x780000000, 0xF00000000, 0x1E000000000, 0x3C0000000000, 0x7800000000000, 0xF000000000000 }; // ^, ^@, ^@@, ...

		/////////////////////////////////
		inline void convert_alreadyTranslatedTokMers(const string& sProteinSequence, const readIDType& iReadID, const int32_t& iNumberOfkMers, uint64_t& iPositionForOut, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& resultsVec) {
			const int32_t& iMaxRange = iNumberOfkMers;
			if (iMaxRange > 0) {
				for (int32_t iCurrentkMerCounter = 0; iCurrentkMerCounter < iMaxRange; ++iCurrentkMerCounter) {
					auto kMer = aminoacidTokMer(sProteinSequence.cbegin() + iCurrentkMerCounter, sProteinSequence.cbegin() + iCurrentkMerCounter + 12);
					
					if (bUnfunny) {
						kMer = aminoAcidsToAminoAcid(kMer);
					}

					get<1>(resultsVec[iPositionForOut]) = kMer;
					get<3>(resultsVec[iPositionForOut++]) = iReadID;
				}
			}
		}

		/////////////////////////////////
		inline void convert_dnaTokMer(const string& sDna, const readIDType& iID, const int32_t& iNumberOfkMers, uint64_t& iPositionForOut, const int32_t& iMaxKTimes3, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& resultsVec, const int32_t& iThreadID, unique_ptr<uint64_t[]>& vCountGarbagekMerPerK) {
			// This value gives the remaining length of the read which is the number of kMers created
			const int32_t& iMaxRange = iNumberOfkMers;
			const auto& garbageRangeCalc = iThreadID * _iNumOfK;
			if (iMaxRange >= 1) {
				// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
				const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;

				int32_t ikMerCounter = 0;

				// Frameshifting, so that only one amino acid has to be computed
				// Compute initial frames
				AAFrames<uint64_t> sAAFrames;

				for (int32_t j = 0; j < iNumFrames; ++j) {
					string sTempFrame = _sMaxKBlank;
					dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
					sAAFrames[j] = aminoacidTokMer(sTempFrame);

					auto kmerForSearchInTrie = sAAFrames[j];
					if (bUnfunny) {
						kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
					}

					for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
						vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal];
					}

					get<1>(resultsVec[iPositionForOut]) = kmerForSearchInTrie;
					get<3>(resultsVec[iPositionForOut]) = static_cast<uint32_t>(iID);
					++iPositionForOut;
				}

				ikMerCounter += iNumFrames;

				if (iMaxRange > 3) {

					// If the Dna has an irregular (not mod 3) tail, a special case must be handled and an additional step in the loop is necessary thus the !(iMaxRange % 3)
					// To map the 2 onto 1, ! is applied twice
					const int32_t& iMaxRangeMod3 = iMaxRange % 3;
					const int32_t& iMaxRangeMod3Neg = !(!iMaxRangeMod3);

					for (int32_t j = 1; 3 * (j + iMaxRangeMod3Neg) < iMaxRange; j++) {
						for (int32_t k = 0; k < 3; ++k) {
							int8_t sTempAA = ' ';
							dnaToAminoacid(sDna, k + iMaxKTimes3 + 3 * (j - 1), sTempAA);
							sAAFrames[k] = aminoacidTokMer(sAAFrames[k], sTempAA);

							auto kmerForSearchInTrie = sAAFrames[k];
							if (bUnfunny) {
								kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
							}

							for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
								vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal];
							}

							get<1>(resultsVec[iPositionForOut]) = kmerForSearchInTrie;
							get<3>(resultsVec[iPositionForOut]) = static_cast<uint32_t>(iID);
							++iPositionForOut;
						}
						ikMerCounter += 3;
					}

					// compute frames of said tail
					for (int32_t j = 0; j < iMaxRangeMod3; ++j) {
						int8_t sTempAA = ' ';
						dnaToAminoacid(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
						sAAFrames[j] = aminoacidTokMer(sAAFrames[j], sTempAA);

						auto kmerForSearchInTrie = sAAFrames[j];
						if (bUnfunny) {
							kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
						}

						for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
							vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal];
						}

						get<1>(resultsVec[iPositionForOut]) = kmerForSearchInTrie;
						get<3>(resultsVec[iPositionForOut]) = static_cast<uint32_t>(iID);
						++iPositionForOut;
					}
				}
			}
		}

		/////////////////////////////////
		inline void convertLinesTokMers_new(const size_t& start, const size_t& end, uint64_t iPositionForOut, const vector<tuple<string, readIDType, int32_t>>& vLines, const vector<tuple<string, readIDType, int32_t>>& vRCLines, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& kMerVecOut, const int32_t& iThreadID, unique_ptr<uint64_t[]>& vCountGarbagekMerPerK) {

			const int32_t& iMaxKTimes3 = 3 * _iHighestK;

			for (uint64_t i = start; i < end; ++i) {
				if (get<0>(vLines[i]) != "" || get<0>(vRCLines[i]) != "") {
					if (_bInputAreAAs) {
						convert_alreadyTranslatedTokMers(get<0>(vLines[i]), get<1>(vLines[i]), get<2>(vLines[i]), iPositionForOut, kMerVecOut);
					}
					else {
						convert_dnaTokMer(get<0>(vLines[i]), get<1>(vLines[i]), get<2>(vLines[i]), iPositionForOut, iMaxKTimes3, kMerVecOut, iThreadID, vCountGarbagekMerPerK);
						convert_dnaTokMer(get<0>(vRCLines[i]), get<1>(vRCLines[i]), get<2>(vRCLines[i]), iPositionForOut, iMaxKTimes3, kMerVecOut, iThreadID, vCountGarbagekMerPerK);
					}
				}
			}
		}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		struct strTransfer {
			string name = "", overhang = "";
			size_t lengthOfDNA = 0;
			bool finished = true, addTail = false, bNewRead = true;
			uint8_t iExpectedInput = 0;
			string lastLine;
			//list<readIDType> vReadIDs;
			//readIDType iCurrentReadID = 0;
			//unordered_map<readIDType, uint64_t> mReadIDToArrayIdx;
			uint64_t iNumOfCharsRead = 0, iCurrentPercentage = 0, iNumOfAllCharsRead = 0, iCurrentOverallPercentage = 0, iNumOfNewReads = 0;
		};

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fastq as input for comparison
		template<typename T>
		inline uint64_t readFastq_partialSort(T& input, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, unique_ptr<strTransfer>& transfer, vector<WorkerThread>& threadPool, vector<uint64_t>& vCountGarbagekMerPerK) {
			uint64_t iSumOfkMers = 0;
			try {
				bool bNotFull = true;
				readIDType iLocalReadID = 0;

				auto calculatekMerCount = [&,this](const string& str) {
					if (_bInputAreAAs) {
						return static_cast<int32_t>(str.length() - _iHighestK + 1);
					}
					else {
						return static_cast<int32_t>(str.length() - 3*_iHighestK + 1);
					}
				};

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
				bool bNewRead = transfer->bNewRead, bAddTail = transfer->addTail;

				vector<tuple<string, readIDType, int32_t>> vLines, vRCLines;
				if (transfer->lastLine != "") {
					auto lastLineKmerSize = calculatekMerCount(transfer->lastLine);
					vLines.emplace_back(transfer->lastLine, 0ul, lastLineKmerSize);
					vRCLines.emplace_back("", 0ul, 0);

					iSumOfkMers += lastLineKmerSize;
					iSoftMaxSize -= lastLineKmerSize * sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>);
					if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
						bNotFull = false;
					}
				}
				transfer->iNumOfNewReads = 0 + (transfer->addTail);
				transfer->finished = false;
				char bufferArrayForInput[10000]; //100000

				while (input && bNotFull) {
					uint64_t iNumOfChars = 0;
					pair<std::string, bool> resultChunkPair;
					Utilities::getChunk(input, move(resultChunkPair), bufferArrayForInput, move(iNumOfChars));

					if (_bVerbose && iFileLength != 0) {
						transfer->iNumOfCharsRead += iNumOfChars;
						transfer->iNumOfAllCharsRead += iNumOfChars;
						double dPercentageOfInputRead = transfer->iNumOfCharsRead / double(iFileLength) * 100.;
						if (static_cast<uint64_t>(dPercentageOfInputRead) != transfer->iCurrentPercentage) {
							transfer->iCurrentPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
							cout << "OUT: Progress of current file " << transfer->iCurrentPercentage << "%" << endl;
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

					const uint8_t iCurrentExpInput = iExpectedInput;
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
							//transfer->vReadIDs.push_back(transfer->iCurrentReadID);
							//transfer->mReadIDToArrayIdx[transfer->iCurrentReadID] = iLocalReadID;
							transfer->iNumOfNewReads++;
							iSoftMaxSize -= iAmountOfSpecies * sizeof(float);// + sizeof(readIDType) + sizeof(pair<readIDType, uint64_t>);
						}

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
								tempDNA += sFalsekMerMarker;
								vRCLines.emplace_back(tempDNA, iLocalReadID, calculatekMerCount(tempDNA));
							}
						}
						else {
							if (!_bInputAreAAs) {
								const auto& rc = reverseComplement(sDNA);
								vRCLines.emplace_back(rc, iLocalReadID, calculatekMerCount(rc));
							}
						}

						vLines.emplace_back(sDNA, iLocalReadID, calculatekMerCount(sDNA));

						bAddTail = true;

						break;

					case 2:
						if (bAddTail) {
							// Pad very small reads
							if (_bInputAreAAs) {
								while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
									get<0>(vLines.back()) += '^';
								}
							}
							else {
								while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
									get<0>(vLines.back()) += 'X';
								}
							}
							
							get<0>(vLines.back()) += sFalsekMerMarker; // add false marker
							get<2>(vLines.back()) += static_cast<int32_t>(sFalsekMerMarker.size());

							iSumOfkMers += sFalsekMerMarker.size();
							iSoftMaxSize -= sFalsekMerMarker.size() * sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>);

							if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
								bNotFull = false;
							}

							if (bReadIDsAreInteresting) {
								//++(transfer->iCurrentReadID);
								++iLocalReadID;
								vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));// + sFalsekMerMarker.length())));
								iSoftMaxSize -= sizeof(pair<string, uint32_t>) + sName.size() * sizeof(char) + sizeof(uint32_t);
								transfer->finished = true;
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
						// read completely parsed

						break;
					default:
						break;
					}

					if (bAddTail) {
						const auto& kMersForward = (vLines.size()) ? get<2>(vLines.back()) : 0;
						const auto& kMersBackwards = (vRCLines.size()) ? get<2>(vRCLines.back()) : 0;
						const auto& iNumOfkMers = kMersForward + kMersBackwards;
						iSumOfkMers += iNumOfkMers;

						iSoftMaxSize -= iNumOfkMers * sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>);

						if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
							bNotFull = false;
						}
					}
				}

				// convert it in parallel
				//auto startTIME = std::chrono::high_resolution_clock::now();
				try {
					vOut.resize(iSumOfkMers);
				}
				catch (const bad_alloc&) {
					cerr << "ERROR: Not enough memory available. Please try again with a lower number after - m" << endl;
					throw;
				}
				const auto& chunkSize = iSumOfkMers / _iNumOfThreads;
				size_t start = 0;
				uint64_t iCurrentkMerCount = 0, iTotalkMerCount = 0;
				int32_t iProcID = 0;
				unique_ptr<uint64_t[]> garbagekMersArray(new uint64_t[_iNumOfK*_iNumOfThreads]);
				memset(garbagekMersArray.get(), 0, _iNumOfK*_iNumOfThreads * sizeof(uint64_t));

				for (size_t iLineIdx = 0; iLineIdx < vLines.size(); ++iLineIdx) {
					if (iCurrentkMerCount >= chunkSize) {
						//call function with start, iLineIdx, iTotalkMerCount
						auto task = [&, start, iLineIdx, iTotalkMerCount, this](const int32_t& iTid) { convertLinesTokMers_new(start, iLineIdx, iTotalkMerCount, ref(vLines), ref(vRCLines), ref(vOut), iTid, garbagekMersArray); };
						threadPool[iProcID].pushTask(task);
						start = iLineIdx;
						iTotalkMerCount += iCurrentkMerCount;
						iCurrentkMerCount = 0;
						iProcID = (iProcID + 1) % _iNumOfThreads;
					}
					iCurrentkMerCount += get<2>(vLines[iLineIdx]) + get<2>(vRCLines[iLineIdx]);
				}
				// call function with start, vLines.size(), iTotalkMerCount
				auto task = [&, start, iTotalkMerCount, this](const int32_t& iTid) { convertLinesTokMers_new(start, vLines.size(), iTotalkMerCount, ref(vLines), ref(vRCLines), ref(vOut), iTid, garbagekMersArray); };
				threadPool[iProcID].pushTask(task);

				// start Threads, join threads
				for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
					threadPool[iThreadID].startThread();
				}
				for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
					threadPool[iThreadID].waitUntilFinished();
				}

				for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
					for (int32_t ik = 0; ik < _iNumOfK; ++ik) {
						vCountGarbagekMerPerK[ik] += garbagekMersArray[iThreadID*_iNumOfK + ik];
					}
				}
				//auto endTIME = std::chrono::high_resolution_clock::now();
				//cout << "Convert " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;

				if (!input) {
					transfer->addTail = false; //signal that all is done
				}
				else {
					transfer->lastLine = get<0>(vLines.back());
					transfer->addTail = bAddTail;
					transfer->bNewRead = bNewRead;
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
		inline uint64_t readFasta_partialSort(T& input, vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, unique_ptr<strTransfer>& transfer, vector<WorkerThread>& threadPool, vector<uint64_t>& vCountGarbagekMerPerK) {
			uint64_t iSumOfkMers = 0;
			try {
				bool bNotFull = true;
				readIDType iLocalReadID = 0;

				auto calculatekMerCount = [&, this](const string& str) {
					if (_bInputAreAAs) {
						return static_cast<int32_t>(str.length() - _iHighestK + 1);
					}
					else {
						return static_cast<int32_t>(str.length() - 3 * _iHighestK + 1);
					}
				};

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
				bool bNewRead = transfer->bNewRead, bAddTail = transfer->addTail;


				vector<tuple<string, readIDType, int32_t>> vLines, vRCLines;
				if (transfer->lastLine != "") {
					auto lastLineKmerSize = calculatekMerCount(transfer->lastLine);
					vLines.emplace_back(transfer->lastLine, 0ul, lastLineKmerSize);
					vRCLines.emplace_back("", 0ul, 0);

					iSumOfkMers += lastLineKmerSize;
					iSoftMaxSize -= lastLineKmerSize * sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>);
					if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
						bNotFull = false;
					}
				}
				transfer->iNumOfNewReads = 0 + (transfer->addTail);
				transfer->finished = false;

				char bufferArrayForInput[10000]; //100000

				while (input && bNotFull) {
					uint64_t iNumOfChars = 0;
					pair<std::string, bool> resultChunkPair;
					Utilities::getChunk(input, move(resultChunkPair), bufferArrayForInput, move(iNumOfChars));

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
								while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
									get<0>(vLines.back()) += '^';
								}
							}
							else {
								while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
									get<0>(vLines.back()) += 'X';
								}
							}

							get<0>(vLines.back()) += sFalsekMerMarker; // add false marker
							get<2>(vLines.back()) += static_cast<int32_t>(sFalsekMerMarker.size());

							iSumOfkMers += sFalsekMerMarker.size();
							iSoftMaxSize -= sFalsekMerMarker.size() * sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>);

							if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
								bNotFull = false;
							}

							if (bReadIDsAreInteresting) {
								//++(transfer->iCurrentReadID);
								transfer->finished = true;
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
							//transfer->vReadIDs.push_back(transfer->iCurrentReadID);
							//transfer->mReadIDToArrayIdx[transfer->iCurrentReadID] = iLocalReadID;
							transfer->iNumOfNewReads++;
							iSoftMaxSize -= iAmountOfSpecies * sizeof(float); //+ sizeof(readIDType) + sizeof(pair<readIDType, uint64_t>);
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
								tempDNA += sFalsekMerMarker;
								vRCLines.emplace_back(tempDNA, iLocalReadID, calculatekMerCount(tempDNA));
							}
						}
						else {
							if (!_bInputAreAAs) {
								const auto& rc = reverseComplement(sDNA);
								vRCLines.emplace_back(rc, iLocalReadID, calculatekMerCount(rc));
							}
						}

						vLines.emplace_back(sDNA, iLocalReadID, calculatekMerCount(sDNA));

						bAddTail = true;

						break;
					default:
						break;
					}

					if (bAddTail) {
						const auto& kMersForward = (vLines.size()) ? get<2>(vLines.back()) : 0;
						const auto& kMersBackwards = (vRCLines.size()) ? get<2>(vRCLines.back()) : 0;
						const auto& iNumOfkMers = kMersForward + kMersBackwards;
						iSumOfkMers += iNumOfkMers;

						iSoftMaxSize -= iNumOfkMers * sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>);


						if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
							bNotFull = false;
						}
					}

				}

				if (!input) {
					if (_bInputAreAAs) {
						while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
							get<0>(vLines.back()) += '^';
						}
					}
					else {
						while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
							get<0>(vLines.back()) += 'X';
						}
					}

					get<0>(vLines.back()) += sFalsekMerMarker; // add false marker
					get<2>(vLines.back()) += static_cast<int32_t>(sFalsekMerMarker.size());
					iSumOfkMers += get<2>(vLines.back());

					if (bReadIDsAreInteresting) {
						transfer->finished = true;
						vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));
					}
					transfer->addTail = false;
				}
				else {
					transfer->lastLine = get<0>(vLines.back());
					transfer->addTail = bAddTail;
					transfer->bNewRead = bNewRead;
					transfer->iExpectedInput = iExpectedInput;
					transfer->lengthOfDNA = iDNALength;
					transfer->name = sName;
					transfer->overhang = sOverhang;
				}

				// convert it in parallel
				try {
					vOut.resize(iSumOfkMers);
				}
				catch (const bad_alloc&) {
					cerr << "ERROR: Not enough memory available. Please try again with a lower number after - m" << endl;
					throw;
				}
				const auto& chunkSize = iSumOfkMers / _iNumOfThreads;
				size_t start = 0;
				uint64_t iCurrentkMerCount = 0, iTotalkMerCount = 0;
				int32_t iProcID = 0;
				unique_ptr<uint64_t[]> garbagekMersArray(new uint64_t[_iNumOfK*_iNumOfThreads]);
				memset(garbagekMersArray.get(), 0, _iNumOfK*_iNumOfThreads * sizeof(uint64_t));

				for (size_t iLineIdx = 0; iLineIdx < vLines.size(); ++iLineIdx) {
					if (iCurrentkMerCount >= chunkSize) {
						//call function with start, iLineIdx, iTotalkMerCount
						auto task = [&, start, iLineIdx, iTotalkMerCount, this](const int32_t& iTid) { convertLinesTokMers_new(start, iLineIdx, iTotalkMerCount, ref(vLines), ref(vRCLines), ref(vOut), iTid, garbagekMersArray); };
						threadPool[iProcID].pushTask(task);
						start = iLineIdx;
						iTotalkMerCount += iCurrentkMerCount;
						iCurrentkMerCount = 0;
						iProcID = (iProcID + 1) % _iNumOfThreads;
					}
					iCurrentkMerCount += get<2>(vLines[iLineIdx]) + get<2>(vRCLines[iLineIdx]);
				}
				// call function with start, vLines.size(), iTotalkMerCount
				auto task = [&, start, iTotalkMerCount, this](const int32_t& iTid) { convertLinesTokMers_new(start, vLines.size(), iTotalkMerCount, ref(vLines), ref(vRCLines), ref(vOut), iTid, garbagekMersArray); };
				threadPool[iProcID].pushTask(task);

				// start Threads, join threads
				for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
					threadPool[iThreadID].startThread();
				}
				for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
					threadPool[iThreadID].waitUntilFinished();
				}

				for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
					for (int32_t ik = 0; ik < _iNumOfK; ++ik) {
						vCountGarbagekMerPerK[ik] += garbagekMersArray[iThreadID*_iNumOfK + ik];
					}
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
		template<typename T>
		inline void readFasta(T& input, const unordered_map<string, uint32_t>& mAccToID, Build& vBricks, unique_ptr<contentVecType_32p>& vOut, const uint64_t& iFileLength, size_t& overallCharsRead, const size_t& overallFilesSize, const float& fShrinkPercentage) {
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

				auto filesAndSize = Utilities::gatherFilesFromPath(sDirectory);
				for (auto& fileName : filesAndSize.first) {
					if (_bVerbose) {
						cout << "OUT: Current file: " << fileName.first << endl;
					}

					// check if gzipped
					bool isGzipped = fileName.second;

					if (isGzipped) {
						if (_bVerbose) {
							cout << "OUT: File is gzipped, no progress output can be shown." << endl;
						}
						igzstream fastaFile_gz(fileName.first.c_str());
						readFasta(fastaFile_gz, mAccToID, brick, dummy, 0, overallCharsRead, filesAndSize.second, fShrinkPercentage);
					}
					else {
						ifstream fastaFile(fileName.first);
						fastaFile.seekg(0, fastaFile.end);
						const uint64_t& iFileLength = fastaFile.tellg();
						fastaFile.seekg(0, fastaFile.beg);
						readFasta(fastaFile, mAccToID, brick, dummy, iFileLength, overallCharsRead, filesAndSize.second, fShrinkPercentage);
					}
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
					//rename((fOutFile + "_taxOnly").c_str(), fOutFile.c_str());
					ifstream source(fOutFile + "_taxOnly", ios::binary);
					ofstream dest(fOutFile, ios::binary);

					istreambuf_iterator<char> begin_source(source);
					istreambuf_iterator<char> end_source;
					ostreambuf_iterator<char> begin_dest(dest);
					copy(begin_source, end_source, begin_dest);

					source.close();
					dest.close();
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