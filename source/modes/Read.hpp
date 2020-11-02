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
	template<class vecType, class elemType, class intType>
	class Read : public kASA {
		const bool _bInputAreAAs;
		const bool _bUnfunny;
	public:
		bool _bVisualize = false;
		vector<string> _translatedFramesForVisualization;

	public:
		Read(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const bool& bTranslated = false, const string& stxxl_mode = "", const bool& bSixFrames = false, const bool& bUnfunny = false) : kASA(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, stxxl_mode, bSixFrames), _bInputAreAAs(bTranslated), _bUnfunny(bUnfunny) {}
		Read(const kASA& obj, const bool& bTranslated = false, const bool& bUnfunny = false) : kASA(obj), _bInputAreAAs(bTranslated), _bUnfunny(bUnfunny) {}
	protected:

		// test for garbage to correctly count matchable k-mers
		//const uint64_t _tails[12] = { 0x1E, 0x3C0, 0x7800, 0xF0000, 0x1E00000, 0x3C000000, 0x780000000, 0xF00000000, 0x1E000000000, 0x3C0000000000, 0x7800000000000, 0xF000000000000 }; // ^, ^@, ^@@, ...

		/////////////////////////////////
		inline void convert_alreadyTranslatedTokMers(const string& sProteinSequence, const readIDType& iReadID, const int32_t& iNumberOfkMers, uint64_t& iPositionForOut, vector<tuple<uint64_t, intType, uint32_t, uint32_t>>& resultsVec) {
			const int32_t& iMaxRange = iNumberOfkMers;
			if (iMaxRange > 0) {
				for (int32_t iCurrentkMerCounter = 0; iCurrentkMerCounter < iMaxRange; ++iCurrentkMerCounter) {
					auto kMer = aminoacidTokMer<intType>(sProteinSequence.cbegin() + iCurrentkMerCounter, sProteinSequence.cbegin() + iCurrentkMerCounter + _iHighestK);
					
					if (_bUnfunny) {
						kMer = aminoAcidsToAminoAcid(kMer);
					}

					get<1>(resultsVec[iPositionForOut]) = kMer;
					get<3>(resultsVec[iPositionForOut++]) = iReadID;
				}
			}
		}

		/////////////////////////////////
		inline void convert_dnaTokMer(const string& sDna, const readIDType& iID, const int32_t& iNumberOfkMers, uint64_t& iPositionForOut, const int32_t& iMaxKTimes3, vector<tuple<uint64_t, intType, uint32_t, uint32_t>>& resultsVec) {
			// This value gives the remaining length of the read which is the number of kMers created
			const int32_t& iMaxRange = iNumberOfkMers;
			if (iMaxRange >= 1) {
				// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
				const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;
				if (_bVisualize && _translatedFramesForVisualization.size() == 0) {
					for (int32_t i = 0; i < iNumFrames; ++i) {
						_translatedFramesForVisualization.push_back("");
					}
				}

				int32_t ikMerCounter = 0;

				// Frameshifting, so that only one amino acid has to be computed
				// Compute initial frames
				AAFrames<intType> sAAFrames;

				for (int32_t j = 0; j < iNumFrames; ++j) {
					string sTempFrame = _sMaxKBlank;
					dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
					if (_bVisualize) {
						_translatedFramesForVisualization[j] += sTempFrame;
					}
					sAAFrames[j] = aminoacidTokMer<intType>(sTempFrame);

					auto kmerForSearchInTrie = sAAFrames[j];
					if (_bUnfunny) {
						kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
					}

					/*for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
						vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal]; //TODO: Slow
					}*/

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
							if (_bVisualize) {
								_translatedFramesForVisualization[k] += sTempAA;
							}
							if (_iMaxK <= 12) {
								sAAFrames[k] = aminoacidTokMer(static_cast<uint64_t>(sAAFrames[k]), sTempAA);
							}
							else {
								sAAFrames[k] = aminoacidTokMer(sAAFrames[k], sTempAA);
							}

							auto kmerForSearchInTrie = sAAFrames[k];
							if (_bUnfunny) {
								kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
							}

							/*for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
								vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal]; //TODO: Slow
							}*/

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
						if (_bVisualize) {
							_translatedFramesForVisualization[j] += sTempAA;
						}
						if (_iMaxK <= 12) {
							sAAFrames[j] = aminoacidTokMer(static_cast<uint64_t>(sAAFrames[j]), sTempAA);
						}
						else {
							sAAFrames[j] = aminoacidTokMer(sAAFrames[j], sTempAA);
						}

						auto kmerForSearchInTrie = sAAFrames[j];
						if (_bUnfunny) {
							kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
						}

						/*for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
							vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal]; //TODO: Slow
						}*/

						get<1>(resultsVec[iPositionForOut]) = kmerForSearchInTrie;
						get<3>(resultsVec[iPositionForOut]) = static_cast<uint32_t>(iID);
						++iPositionForOut;
					}
				}
			}
		}

		/////////////////////////////////
		inline void convertLinesTokMers_new(const vector<tuple<string, readIDType, int32_t>>& vLines, const vector<tuple<string, readIDType, int32_t>>& vRCLines, const size_t& start, const size_t& end, uint64_t iPositionForOut, vector<tuple<uint64_t, intType, uint32_t, uint32_t>>& kMerVecOut) {

			const int32_t& iMaxKTimes3 = 3 * _iHighestK;
			//const uint64_t iStartPositionForOut = iPositionForOut;
			for (uint64_t i = start; i < end; ++i) {
				if (get<0>(vLines[i]) != "" || get<0>(vRCLines[i]) != "") {
					if (_bInputAreAAs) {
						convert_alreadyTranslatedTokMers(get<0>(vLines[i]), get<1>(vLines[i]), get<2>(vLines[i]), iPositionForOut, kMerVecOut);
					}
					else {
						convert_dnaTokMer(get<0>(vLines[i]), get<1>(vLines[i]), get<2>(vLines[i]), iPositionForOut, iMaxKTimes3, kMerVecOut);
						if (_bSixFrames) {
							convert_dnaTokMer(get<0>(vRCLines[i]), get<1>(vRCLines[i]), get<2>(vRCLines[i]), iPositionForOut, iMaxKTimes3, kMerVecOut);
						}
					}
				}
			}

			
			//sort(kMerVecOut.begin() + iStartPositionForOut, kMerVecOut.begin() + iPositionForOut);
			
		}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public:
		struct strTransfer {
			string name = "", overhang = "", overhang2 = "";
			size_t lengthOfDNA = 0;
			bool finished = true, addTail = false, bNewRead = true;
			uint8_t iExpectedInput = 0;
			string lastLine = "", lastLine2 = "";
			//list<readIDType> vReadIDs;
			//readIDType iCurrentReadID = 0;
			//unordered_map<readIDType, uint64_t> mReadIDToArrayIdx;
			uint64_t iNumOfCharsRead = 0, iCurrentPercentage = 0, iNumOfAllCharsRead = 0, iCurrentOverallPercentage = 0, iNumOfNewReads = 0;
			//vector<uint64_t> vRangesOfOutVec;
		};

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fastq as input for comparison
		template<typename T>
		inline uint64_t readFastqa_partialSort(Utilities::FileReader<T>& input, Utilities::FileReader<T>& input2, vector<tuple<uint64_t, intType, uint32_t, uint32_t>>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, const bool& bIsFasta, unique_ptr<strTransfer>& transfer, vector<WorkerThread>& threadPool) {
			
			vector<tuple<string, readIDType, int32_t>> vLines, vRCLines, vLines2, vRCLines2;
			uint64_t iSumOfkMers = 0;
			const int64_t iAvailMemory = iSoftMaxSize / (1024 * 1024);
			try {
//#define TIME
#ifdef TIME
				auto startTIME = std::chrono::high_resolution_clock::now();
#endif

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
				string sName = transfer->name, sDNA = "", sOverhang = transfer->overhang, sDNA2 = "", sOverhang2 = transfer->overhang2;
				size_t iDNALength = transfer->lengthOfDNA, iQualityLength = 0;
				bool bNewRead = transfer->bNewRead, bAddTail = transfer->addTail;

				transfer->iNumOfNewReads = 0 + (transfer->addTail);
				transfer->finished = false;
				pair<std::string, bool> resultChunkPair, resultChunkPair2;
				bool bReadTooShort = false;

				while (!input.eof() && bNotFull) {
					uint64_t iNumOfChars = 0;
					resultChunkPair.first = "";
					resultChunkPair.second = false;
					input.getChunk(resultChunkPair, iNumOfChars);
					if (input2.notNull()) {
						resultChunkPair2.first = "";
						resultChunkPair2.second = false;
						input2.getChunk(resultChunkPair2, iNumOfChars);
					}

					if (iExpectedInput == 1 && (resultChunkPair.first.length() < static_cast<size_t>((_bInputAreAAs) ? (_iHighestK + 1) : (3 * _iHighestK + 1))) && resultChunkPair.second == false) {
						// currently read DNA string is too short and there is more DNA to be had, usually happens when the end of the buffer is read but the file still contains DNA
						// first, count number of characters read in the first try and then append more DNA
						transfer->iNumOfCharsRead += iNumOfChars;
						transfer->iNumOfAllCharsRead += iNumOfChars;
						input.getChunk(resultChunkPair, iNumOfChars); // new DNA will be appended
						if (input2.notNull()) {
							input2.getChunk(resultChunkPair2, iNumOfChars);
						}
					}

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
						if (input2.notNull()) {
							if (resultChunkPair2.first.front() != '+') {
								throw runtime_error("Paired-end files are out of sync!");
							}
						}
						if (!resultChunkPair.second) {
							input.ignore(iNumOfChars); // discard the rest of the +
							transfer->iNumOfCharsRead += iNumOfChars;
							transfer->iNumOfAllCharsRead += iNumOfChars;
							if (input2.notNull()) {
								input2.ignore(iNumOfChars);
							}
						}
						iExpectedInput = 2;
						continue;
					}

					if (bIsFasta && resultChunkPair.first.front() == '>') {
						if (bAddTail) {
							// Pad very small reads
							auto padding = [&](vector<tuple<string, readIDType, int32_t>>& vLines) {
								if (_bInputAreAAs) {
									while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
										get<0>(vLines.back()) += '^';
										get<2>(vLines.back())++;
									}
								}
								else {
									while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
										get<0>(vLines.back()) += 'X';
										get<2>(vLines.back())++;
									}
								}

								get<0>(vLines.back()) += sFalsekMerMarker; // add false marker
								get<2>(vLines.back()) += static_cast<int32_t>(sFalsekMerMarker.size());

								iSumOfkMers += sFalsekMerMarker.size();
								iSoftMaxSize -= sFalsekMerMarker.size() * sizeof(tuple<uint64_t, intType, uint32_t, uint32_t>);
							};
							padding(vLines);
							if (input2.notNull()) {
								padding(vLines2);
							}

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

						iExpectedInput = 0;
					}

					const uint8_t iCurrentExpInput = iExpectedInput;
					switch (iCurrentExpInput) {
						///////////////////////////////////////////////// case 0
					case 0:
					{
						if (bIsFasta) {
							sName = Utilities::lstrip(resultChunkPair.first, '>');
						}
						else {
							sName = Utilities::lstrip(resultChunkPair.first, '@');
						}

						while (!resultChunkPair.second) {
							resultChunkPair.first = "";
							input.getChunk(resultChunkPair, iNumOfChars); // discard the rest of the name
							transfer->iNumOfCharsRead += iNumOfChars;
							transfer->iNumOfAllCharsRead += iNumOfChars;
							sName += resultChunkPair.first;
							if (input2.notNull()) {
								resultChunkPair2.first = "";
								input2.getChunk(resultChunkPair2, iNumOfChars);
							}
						}
						bNewRead = true;
						sOverhang = "";
						sOverhang2 = "";
						iDNALength = 0;
						iQualityLength = 0;

						//vReadNameAndLength.reserve(vReadNameAndLength.capacity() + 1);

						transfer->iNumOfNewReads++;
						if (bReadIDsAreInteresting) {
							//transfer->vReadIDs.push_back(transfer->iCurrentReadID);
							//transfer->mReadIDToArrayIdx[transfer->iCurrentReadID] = iLocalReadID;
							iSoftMaxSize -= iAmountOfSpecies * sizeof(float);// + sizeof(readIDType) + sizeof(pair<readIDType, uint64_t>);
						}

						//cout << iSoftMaxSize << endl;

						iExpectedInput = 1;
						break;
					}
					///////////////////////////////////////////////// case 1
					case 1:
					{
						auto searchAndReplaceLetters = [&](pair<std::string, bool>& resultChunkPair) {
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
						};
						searchAndReplaceLetters(resultChunkPair);
						if (input2.notNull()) {
							searchAndReplaceLetters(resultChunkPair2);
						}


						auto funcReadTooShort = [&](string& sDNA, string& sOverhang, const pair<std::string, bool>& resultChunkPair) {
							sDNA = sOverhang + resultChunkPair.first;
							iDNALength += resultChunkPair.first.length();
							if (_bInputAreAAs) {
								if (sDNA.length() < size_t(_iHighestK)) {
									sOverhang = sDNA;
									bReadTooShort = true;
								}
								else {
									sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK);
									bReadTooShort = false;
								}
							}
							else {
								if (sDNA.length() < size_t(_iHighestK) * 3 && !resultChunkPair.second) {
									sOverhang = sDNA;
									bReadTooShort = true;
								}
								else {
									if (sDNA.length() < size_t(_iHighestK) * 3 && resultChunkPair.second) {
										sOverhang = "";
										bReadTooShort = false;
										while (sDNA.length() < size_t(_iHighestK) * 3) {
											sDNA += 'X';
										}
									}
									else {
										sOverhang = sDNA.substr(sDNA.length() + 1 - _iHighestK * 3);
										bReadTooShort = false;
									}
								}
							}
						};
						funcReadTooShort(sDNA, sOverhang, resultChunkPair);
						if (input2.notNull()) {
							funcReadTooShort(sDNA2, sOverhang2, resultChunkPair2);
						}

						auto emplaceBack = [&](const string& sDNA, vector<tuple<string, readIDType, int32_t>>& vLines, vector<tuple<string, readIDType, int32_t>>& vRCLines) {
							if (!bReadTooShort) {
								if (bNewRead) {
									bNewRead = false;
									if (!_bInputAreAAs && _bSixFrames) {
										string tempDNA = reverseComplement(sDNA);
										while (tempDNA.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
											tempDNA += 'X';
										}
										tempDNA += sFalsekMerMarker;
										vRCLines.emplace_back(tempDNA, iLocalReadID, calculatekMerCount(tempDNA));
										iSoftMaxSize -= tempDNA.length() * sizeof(char) + sizeof(readIDType) + sizeof(int32_t);
									}
								}
								else {
									if (!_bInputAreAAs && _bSixFrames) {
										const auto& rc = reverseComplement(sDNA);
										vRCLines.emplace_back(rc, iLocalReadID, calculatekMerCount(rc));
										iSoftMaxSize -= rc.length() * sizeof(char) + sizeof(readIDType) + sizeof(int32_t);
									}
								}

								vLines.emplace_back(sDNA, iLocalReadID, calculatekMerCount(sDNA));
								iSoftMaxSize -= sDNA.length() * sizeof(char) + sizeof(readIDType) + sizeof(int32_t);
								bAddTail = true;
							}
						};
						bool bNewReadTemp = bNewRead;
						emplaceBack(sDNA, vLines, vRCLines);
						if (input2.notNull()) {
							bNewRead = bNewReadTemp;
							emplaceBack(sDNA2, vLines2, vRCLines2);
						}

						break;
					}
					///////////////////////////////////////////////// case 2
					case 2:
					{
						if (vLines.size() == 0) {
							if (transfer->lastLine != "") {
								auto lastLineKmerSize = calculatekMerCount(transfer->lastLine);
								vLines.emplace_back(transfer->lastLine, 0ul, lastLineKmerSize);
								if (_bSixFrames) {
									vRCLines.emplace_back("", 0ul, 0);
								}

								if (input2.notNull()) {
									auto lastLineKmerSize2 = calculatekMerCount(transfer->lastLine2);
									vLines2.emplace_back(transfer->lastLine2, 0ul, lastLineKmerSize2);
									if (_bSixFrames) {
										vRCLines2.emplace_back("", 0ul, 0);
									}
									lastLineKmerSize += lastLineKmerSize2;
								}

								iSumOfkMers += lastLineKmerSize;
								iSoftMaxSize -= lastLineKmerSize * sizeof(tuple<uint64_t, intType, uint32_t, uint32_t>);
								if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
									bNotFull = false;
								}
							}
							else {
								throw runtime_error("No overhang line was found, something went wrong while parsing the input!");
							}
						}


						if (bAddTail) {
							// Pad very small reads
							auto padding = [&](vector<tuple<string, readIDType, int32_t>>& vLines) {
								if (_bInputAreAAs) {
									while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
										get<0>(vLines.back()) += '^';
										get<2>(vLines.back())++;
									}
								}
								else {
									while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
										get<0>(vLines.back()) += 'X';
										get<2>(vLines.back())++;
									}
								}

								get<0>(vLines.back()) += sFalsekMerMarker; // add false marker
								get<2>(vLines.back()) += static_cast<int32_t>(sFalsekMerMarker.size());
								iSoftMaxSize -= sFalsekMerMarker.length() * sizeof(char);

								iSumOfkMers += sFalsekMerMarker.size();
								iSoftMaxSize -= sFalsekMerMarker.size() * sizeof(tuple<uint64_t, intType, uint32_t, uint32_t>);
							};
							padding(vLines);
							if (input2.notNull()) {
								padding(vLines2);
							}

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
						if (input2.notNull()) {
							iQualityLength += resultChunkPair2.first.length();
						}
						while (!resultChunkPair.second || iQualityLength < iDNALength) {
							resultChunkPair.first = "";
							iNumOfChars = 0;
							input.getChunk(resultChunkPair, iNumOfChars); // get the whole quality string
							transfer->iNumOfCharsRead += iNumOfChars;
							transfer->iNumOfAllCharsRead += iNumOfChars;
							iQualityLength += resultChunkPair.first.length();
							if (input2.notNull()) {
								resultChunkPair2.first = "";
								input2.getChunk(resultChunkPair2, iNumOfChars);
								iQualityLength += resultChunkPair2.first.length();
							}
						}

						if (iQualityLength == iDNALength) {
							iExpectedInput = 0;
						}
						else {
							if (iQualityLength > iDNALength) {
								throw runtime_error("Quality string length and DNA length don't match. Difference in size: " + to_string(iQualityLength - iDNALength) + ". Error occured in:\n" + sName);
							}
							else {
								throw runtime_error("Quality string length and DNA length don't match. Difference in size: " + to_string(iDNALength - iQualityLength) + ". Error occured in:\n" + sName);
							}
						}
						// read completely parsed

						break;
					}
					default:
						break;
					}///////////////////////////////////////////////// end

					if (bAddTail) {
						const auto& kMersForward = (vLines.size()) ? get<2>(vLines.back()) : 0;
						const auto& kMersBackwards = (vRCLines.size()) ? get<2>(vRCLines.back()) : 0;
						const auto& iNumOfkMers = kMersForward + kMersBackwards;
						iSumOfkMers += iNumOfkMers;

						iSoftMaxSize -= iNumOfkMers * sizeof(tuple<uint64_t, intType, uint32_t, uint32_t>);

						if (input2.notNull()) {
							const auto& kMersForward2 = (vLines2.size()) ? get<2>(vLines2.back()) : 0;
							const auto& kMersBackwards2 = (vRCLines2.size()) ? get<2>(vRCLines2.back()) : 0;
							const auto& iNumOfkMers2 = kMersForward2 + kMersBackwards2;
							iSumOfkMers += iNumOfkMers2;

							iSoftMaxSize -= iNumOfkMers2 * sizeof(tuple<uint64_t, intType, uint32_t, uint32_t>);
						}

						if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
							bNotFull = false;
						}
					}
				}

				if (input.eof()) {
					if (bIsFasta) {
						auto handleEOF = [&](const string& sDNA, vector<tuple<string, readIDType, int32_t>>& vLines, vector<tuple<string, readIDType, int32_t>>& vRCLines) {
							if (bReadTooShort) {
								if (bNewRead) {
									bNewRead = false;
									if (!_bInputAreAAs && _bSixFrames) {
										string tempDNA = reverseComplement(sDNA);
										while (tempDNA.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
											tempDNA += 'X';
										}
										tempDNA += sFalsekMerMarker;
										vRCLines.emplace_back(tempDNA, iLocalReadID, calculatekMerCount(tempDNA));
										iSumOfkMers += get<2>(vRCLines.back());
									}
								}
								else {
									if (!_bInputAreAAs && _bSixFrames) {
										const auto& rc = reverseComplement(sDNA);
										vRCLines.emplace_back(rc, iLocalReadID, calculatekMerCount(rc));
										iSumOfkMers += get<2>(vRCLines.back());
									}
								}

								vLines.emplace_back(sDNA, iLocalReadID, calculatekMerCount(sDNA));
							}

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
						};
						handleEOF(sDNA, vLines, vRCLines);
						if (input2.notNull()) {
							handleEOF(sDNA2, vLines2, vRCLines2);
						}

						if (bReadIDsAreInteresting) {
							transfer->finished = true;
							vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));
						}
					}
					transfer->addTail = false; //signal that all is done
				}
				else {
					transfer->lastLine = sOverhang;
					transfer->lastLine2 = sOverhang2;
					transfer->addTail = bAddTail;
					transfer->bNewRead = bNewRead;
					transfer->iExpectedInput = iExpectedInput;
					transfer->lengthOfDNA = iDNALength;
					transfer->name = sName;
					transfer->overhang = sOverhang;
					transfer->overhang2 = sOverhang2;
				}

#ifdef TIME
				auto endTIME = std::chrono::high_resolution_clock::now();
				cout << "Input file load_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				/*uint64_t checkNumOfkMers = 0;
				for (const auto& entry : vLines) {
					checkNumOfkMers += get<2>(entry);
				}
				for (const auto& entry : vRCLines) {
					checkNumOfkMers += get<2>(entry);
				}
				if (checkNumOfkMers != iSumOfkMers) {
					cout << iSumOfkMers << " " << checkNumOfkMers << endl;
				}*/

				// in case of paired-end reads, merge both vectors containing the reads
				if (input2.notNull()) {
					vLines.insert(vLines.end(), vLines2.begin(), vLines2.end());
					vLines2.clear();
					vRCLines.insert(vRCLines.end(), vRCLines2.begin(), vRCLines2.end());
					vRCLines2.clear();
				}


				// convert it in parallel
				//auto startTIME = std::chrono::high_resolution_clock::now();
				try {
					// should the capacity be smaller than iSumOfkMers after the first allocation, bad_alloc will probably be thrown because there will likely be no larger contiguous chunk of memory available
					// this is mitigated by substracting 1% of the available memory after every time this is called (see Compare.hpp "reduce available memory")
					vOut.reserve(iSumOfkMers); 
					vOut.resize(iSumOfkMers);
				}
				catch (const bad_alloc&) {
					int64_t triedToAllocate = iSumOfkMers * sizeof(tuple<uint64_t, intType, uint32_t, uint32_t>) / (1024 * 1024);
					if (iAvailMemory < triedToAllocate) {
						cerr << "ERROR: Your system does not have enough contiguous memory available (which happens if the system is powered on over a long period of time). You might try to restart or use a lower number after -m. FYI: You tried to use " + to_string(triedToAllocate) + " MB." << endl;	
					} else {
						cerr << "ERROR: Not enough memory available. You tried to use " + to_string(triedToAllocate) + " MB. Please try again with a lower number after -m." << endl;
					}
					throw;
				}
				const auto& chunkSize = iSumOfkMers / _iNumOfThreads;
				size_t start = 0;
				uint64_t iCurrentkMerCount = 0, iTotalkMerCount = 0;
				int32_t iProcID = 0;

#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "Reserve memory_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				//vector<uint64_t> vRangesOfOutVec(_iNumOfThreads + 1);
				//cout << "Memory left: " << iSoftMaxSize << endl;
				for (size_t iLineIdx = 0; iLineIdx < vLines.size(); ++iLineIdx) {
					if (iCurrentkMerCount >= chunkSize) {
						//call function with start, iLineIdx, iTotalkMerCount
						auto task = [&, start, iLineIdx, iTotalkMerCount, this](const int32_t&) { convertLinesTokMers_new(vLines, vRCLines, start, iLineIdx, iTotalkMerCount, ref(vOut)); };
						threadPool[iProcID].pushTask(task);
						//vRangesOfOutVec[iProcID] = iTotalkMerCount;
						//cout << iTotalkMerCount << endl;
						start = iLineIdx;
						iTotalkMerCount += iCurrentkMerCount;
						iCurrentkMerCount = 0;
						++iProcID;
					}
					if (!_bSixFrames) {
						iCurrentkMerCount += uint64_t(get<2>(vLines[iLineIdx])); //+ get<2>(vRCLines[iLineIdx]);
					}
					else {
						iCurrentkMerCount += uint64_t(get<2>(vLines[iLineIdx])) + get<2>(vRCLines[iLineIdx]);
					}
				}
				// call function with start, vLines.size(), iTotalkMerCount
				auto task = [&, start, iTotalkMerCount, this](const int32_t&) { convertLinesTokMers_new(vLines, vRCLines, start, vLines.size(), iTotalkMerCount, ref(vOut)); };
				threadPool[iProcID].pushTask(task);
				//cout << iTotalkMerCount << endl;
				//vRangesOfOutVec[iProcID] = iTotalkMerCount;
				//vRangesOfOutVec[iProcID + 1] = vOut.size();
				//transfer->vRangesOfOutVec = vRangesOfOutVec;

				// start Threads, join threads
				for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
					threadPool[iThreadID].startThread();
				}
				for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
					threadPool[iThreadID].waitUntilFinished();
				}

#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "PAR Input translation_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				// merge sorted ranges. If number of threads is odd, one extra run must be done at the end
				// even part
				/*for (uint32_t j = 1; j < _iNumOfThreads; j *= 2) {
					for (uint32_t i = 0; i < _iNumOfThreads / j; i += 2) {
						//cout << i << " " << j << " " << vRangesOfOutVec[i * j] << " " << vRangesOfOutVec[(i + 1) * j] << " " << vRangesOfOutVec[(i + 2) * j] << endl;
						Utilities::my_inplace_merge(vOut.begin() + vRangesOfOutVec[i*j], vOut.begin() + vRangesOfOutVec[(i + 1)*j], vOut.begin() + vRangesOfOutVec[(i + 2)*j]);
					}
				}*/
				// odd part

				//cout << "is sorted? " << is_sorted(vOut.begin(), vOut.end()) << endl;

				//auto endTIME = std::chrono::high_resolution_clock::now();
				//cout << "Convert " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;

			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
			////////////////////////////////////////////////////////
			//this->vLines.clear();
			//this->vRCLines.clear();
			return iSumOfkMers;
		}

	protected:

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
		inline void dnaTokMers(const string& sDna, const uint32_t& iID, unique_ptr<vecType>& kMerVecOut, Build<vecType, elemType>& vBricks, const float& fShrinkPercentage) {
			vector<tuple<intType, uint32_t>> vResultingkMers;
			const int32_t& iMaxKTimes3 = 3 * _iHighestK;

			const int32_t& iMaxRange = int32_t(sDna.length()) - iMaxKTimes3 + 1;
			if (iMaxRange >= 1) {
				// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
				const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : iMaxRange;

				vResultingkMers.resize(iMaxRange);
				int32_t ikMerCounter = 0;

				// Frameshifting, so that only one amino acid has to be computed
				// Compute initial frames
				intType sAAFrames[3] = { 0, 0, 0 };
				uint32_t aDeletekMerCounter[3] = { 0, 0, 0 };
				for (int32_t j = 0; j < iNumFrames; ++j) {
					string sTempFrame = _sMaxKBlank;
					dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
					sAAFrames[j] = aminoacidTokMer<intType>(sTempFrame);
					// if a character is 'illegal' then don't save the kMer containing it
					const auto& iPosOfU = sTempFrame.find_last_of('_');
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
							if (sTempAA != '_' && aDeletekMerCounter[k] == 0) {
								vResultingkMers[ikMerCounter + k] = make_tuple(sAAFrames[k], iID);
							}
							else {
								if (aDeletekMerCounter[k] != 0 && sTempAA != '_') {
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
						if (sTempAA != '_' && aDeletekMerCounter[j] == 0) {
							vResultingkMers[ikMerCounter + j] = make_tuple(sAAFrames[j], iID);
						}
						else {
							if (aDeletekMerCounter[j] != 0 && sTempAA != '_') {
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
							if (_bUnfunny) {
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
						if (_bUnfunny) {
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
		inline void proteinTokMers(const string& sAASequence, const uint32_t& iIdx, unique_ptr<vecType>& kMerVecOut, Build<vecType, elemType>& vBricks, const float& fShrinkPercentage) {
			vector<tuple<intType, uint32_t>> vResultingkMers;
			const int32_t& iMaxRange = int32_t(sAASequence.length()) - _iHighestK + 1;
			if (iMaxRange > 0) {
				vResultingkMers.resize(iMaxRange);
				for (int32_t i = 0; i < iMaxRange; ++i) {
					const auto& kMer = aminoacidTokMer<intType>(sAASequence.cbegin() + i, sAASequence.cbegin() + i + _iHighestK);
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
						if (_bUnfunny) {
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
					if (_bUnfunny) {
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

	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fasta and create a kMer-Vec, used in BuildAll(...)
		template<typename T>
		inline void readFasta(T& input, const unordered_map<string, uint32_t>& mAccToID, Build<vecType, elemType>& vBricks, unique_ptr<vecType>& vOut, const uint64_t& iFileLength, size_t& overallCharsRead, const size_t& overallFilesSize, const float& fShrinkPercentage) {
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
						char sTempCArr[10000];
						input.get(sTempCArr, 10000);

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
						if (iNumOfChars != 9999 && input) {
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
									input.get(sTempCArr, 10000);
									iNumOfChars = input.gcount();
									if (iNumOfChars != 9999 && input) {
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
					string sRCDNA = (_bSixFrames) ? reverseComplement(sDNA) : "";

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
		// Build the index from fasta files
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
				bool bTaxIdsAsStrings = false;
				while (getline(content, sDummy)) {
					if (sDummy != "") {
						const auto& line = Utilities::split(sDummy, '\t');
						if (line.size() >= 5 && !bTaxIdsAsStrings) {
							bTaxIdsAsStrings = true;
						}
						if (line.size() >= 4) {
							mIdxToName[iIdxCounter] = line[0];
							const auto& vAccessionNumbers = Utilities::split(line[3], ';');
							if (bTaxIdsAsStrings) {
								mIDsAsIdx[stoul(line[4])] = iIdxCounter;
								for (const auto& acc : vAccessionNumbers) {
									mAccToID.insert(make_pair(acc, stoul(line[4])));
								}
							}
							else {
								mIDsAsIdx[stoul(line[1])] = iIdxCounter;
								for (const auto& acc : vAccessionNumbers) {
									mAccToID.insert(make_pair(acc, stoul(line[1])));
								}
							}
							++iIdxCounter;
						}
						else {
							throw runtime_error("Content file contains less than 4 columns, it may be damaged... The faulty line was: " + sDummy + "\n");
						}
					}
				}

				Build<vecType, elemType> brick(_sTemporaryPath, _iNumOfCall, _iNumOfThreads, iMem / (sizeof(elemType)), iIdxCounter);

				size_t overallCharsRead = 0;
				unique_ptr<vecType> dummy;

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

				if (iSizeOfFinalIndex == 0) {
					throw runtime_error("Index is empty, maybe the input were translated sequences and you forgot -z?");
				}

				unique_ptr<uint64_t[]> arrFrequencies;
				arrFrequencies.reset(new uint64_t[iIdxCounter * _iHighestK]);
				for (uint64_t i = 0; i < iIdxCounter * _iHighestK; ++i) {
					arrFrequencies[i] = 0;
				}

				// Create Trie and frequencies out of final file
				stxxlFile* libFile = new stxxlFile(fOutFile, stxxlFile::RDONLY);
				const vecType* libVec = new const vecType(libFile, iSizeOfFinalIndex);
				Trie<intType> T(static_cast<int8_t>(((_iMaxK > 12) ? HIGHESTPOSSIBLEK : 12)), static_cast<int8_t>(_iMinK), 6);
				T.SaveToStxxlVec(libVec, fOutFile, &arrFrequencies, _iHighestK, mIDsAsIdx);

				// If taxaOnly index is desired, create it here:
				if (_bUnfunny) {
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
					Utilities::copyFile(fOutFile + "_taxOnly", fOutFile);
				}

				ofstream outFile(fOutFile + "_f.txt");
				for (uint32_t j = 0; j < iIdxCounter; ++j) {
					outFile << Utilities::checkIfInMap(mIdxToName, j)->second << "\t";
					outFile << arrFrequencies[j * _iHighestK];
					for (int32_t k = 1; k < _iHighestK; ++k) {
						outFile << "\t" << arrFrequencies[j * _iHighestK + k];
					}
					outFile << endl;
				}
			}
			catch (invalid_argument&) {
				cerr << "ERROR: content file doesn't have the right format. If you'd like to use non-numeric taxids, please apply the --taxidasstr flag." << endl;
				return;
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Merge two existing indices into a new one
		void MergeTwoIndices(const string& index_1, const string& index_2, const string& index_out, const string& contentFile) {
			try {
				// test if files exists
				if (!ifstream(contentFile)) {
					throw runtime_error("Content file not found.");
				}

				// content file
				uint32_t iIdxCounter = 1;
				unordered_map<uint32_t, uint32_t> mIDsAsIdx; mIDsAsIdx[0] = 0;
				unordered_map<uint32_t, string> mIdxToName; mIdxToName[0] = "non_unique";
				ifstream content(contentFile);
				string sDummy = "";
				bool bTaxIdsAsStrings = false;
				while (getline(content, sDummy)) {
					if (sDummy != "") {
						const auto& line = Utilities::split(sDummy, '\t');
						if (line.size() >= 5 && !bTaxIdsAsStrings) {
							bTaxIdsAsStrings = true;
						}
						if (line.size() >= 4) {
							mIdxToName[iIdxCounter] = line[0];
							if (bTaxIdsAsStrings) {
								mIDsAsIdx[stoul(line[4])] = iIdxCounter;
							}
							else {
								mIDsAsIdx[stoul(line[1])] = iIdxCounter;
							}
							++iIdxCounter;
						}
						else {
							throw runtime_error("Content file contains less than 4 columns, it may be damaged...");
						}
					}
				}

				Build<vecType, elemType> brick;

				// indices
				ifstream fLibInfo(index_1 + "_info.txt");
				uint64_t iSizeOfFirstLib = 0, iSizeOfSecondLib = 0;
				fLibInfo >> iSizeOfFirstLib;
				fLibInfo.close();
				unique_ptr<stxxlFile> stxxlVecI1(new stxxlFile(index_1, stxxl::file::RDONLY));
				unique_ptr<vecType> vec1(new vecType(stxxlVecI1.get(), iSizeOfFirstLib));
				typename vecType::bufreader_type vec1Buff(*vec1);

				fLibInfo.open(index_2 + "_info.txt");
				fLibInfo >> iSizeOfSecondLib;
				fLibInfo.close();
				unique_ptr<stxxlFile> stxxlVecI2(new stxxlFile(index_2, stxxl::file::RDONLY));
				unique_ptr<vecType> vec2(new vecType(stxxlVecI2.get(), iSizeOfSecondLib));
				

				Utilities::createFile(index_out);
				unique_ptr<stxxlFile> stxxlVecOut(new stxxlFile(index_out, stxxl::file::RDWR));
				unique_ptr<vecType> vecOut(new vecType(stxxlVecOut.get(), iSizeOfFirstLib + iSizeOfSecondLib));
				typename vecType::bufwriter_type vecOutBuff(*vecOut);

				//merge
				unique_ptr<uint64_t[]> arrFrequencies;
				arrFrequencies.reset(new uint64_t[iIdxCounter * _iHighestK]);
				for (uint64_t i = 0; i < iIdxCounter * _iHighestK; ++i) {
					arrFrequencies[i] = 0;
				}

				const auto& vOutSize = brick.merge(vecOutBuff, vec1Buff, iSizeOfFirstLib, vec2->cbegin(), vec2->cend(), arrFrequencies, mIDsAsIdx, true);
				vecOut->resize(vOutSize, true);

				// save additional stuff
				vecOut->export_files("_");
				ofstream fOutInfo(index_out + "_info.txt");
				fOutInfo << vecOut->size();
				if (is_same<vecType, contentVecType_128>::value) {
					fOutInfo << endl << 128;
				}

				ofstream outFile(index_out + "_f.txt");
				for (uint32_t j = 0; j < iIdxCounter; ++j) {
					outFile << Utilities::checkIfInMap(mIdxToName, j)->second << "\t";
					outFile << arrFrequencies[j * _iHighestK];
					for (int32_t k = 1; k < _iHighestK; ++k) {
						outFile << "\t" << arrFrequencies[j * _iHighestK + k];
					}
					outFile << endl;
				}

				// trie
				Trie<intType> T(static_cast<int8_t>(((is_same<vecType, contentVecType_128>::value) ? HIGHESTPOSSIBLEK : 12)), static_cast<int8_t>(_iMinK), 6);
				T.SaveToStxxlVec(vecOut.get(), index_out);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
	};
}