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
		const bool _bUnfunny;

	public:
		bool _bVisualize = false;
		bool _bOnlyOneFrame = false;
		vector<string> _translatedFramesForVisualization;

	public:
		Read(const InputParameters& cParams) : kASA(cParams), _bUnfunny(cParams.bUnfunny) {}
		Read(const kASA& obj, const bool& bUnfunny = false) : kASA(obj), _bUnfunny(bUnfunny) {}
	protected:

		// test for garbage to correctly count matchable k-mers
		//const uint64_t _tails[12] = { 0x1E, 0x3C0, 0x7800, 0xF0000, 0x1E00000, 0x3C000000, 0x780000000, 0xF00000000, 0x1E000000000, 0x3C0000000000, 0x7800000000000, 0xF000000000000 }; // ^, ^@, ^@@, ...

		/////////////////////////////////
		inline uint64_t calculatekMerCount(const size_t& strLength) {
			if (this->_bProtein) {
				if (strLength > static_cast<size_t>(_iHighestK + 1)) {
					return strLength - _iHighestK + 1;
				}
			}
			else {
				if (_bOnlyOneFrame) {
					const auto& lengthDividedBy3 = strLength / 3;
					if (lengthDividedBy3 > static_cast<size_t>(_iHighestK + 1)) {
						return lengthDividedBy3 - _iHighestK + 1;
					}
				}
				else {
					if (strLength > static_cast<size_t>(3 * _iHighestK + 1)) {
						return strLength - 3 * _iHighestK + 1;
					}
				}
			}

			return 0;
		}

		/////////////////////////////////
		inline void convert_alreadyTranslatedTokMers(const string& sProteinSequence, const readIDType& iReadID, const uint64_t& iNumberOfkMers, uint64_t& iPositionForOut, uint32_t iPositionInString, InputType<intType>& resultsVec) {
			const uint64_t& iMaxRange = iNumberOfkMers;
			if (iMaxRange > 0) {
				for (uint64_t iCurrentkMerCounter = 0; iCurrentkMerCounter < iMaxRange; ++iCurrentkMerCounter) {
					auto kMer = aminoacidTokMer<intType>(sProteinSequence.cbegin() + iCurrentkMerCounter, sProteinSequence.cbegin() + iCurrentkMerCounter + _iHighestK);
					
					if (_bUnfunny) {
						kMer = aminoAcidsToAminoAcid(kMer);
					}

					resultsVec.setkMer(iPositionForOut, kMer);
					resultsVec.setReadID(iPositionForOut,iReadID);
					if (GlobalInputParameters.bPostProcess) {
						resultsVec.setFrame(iPositionForOut, 0);
						resultsVec.setPosition(iPositionForOut, iPositionInString);
						iPositionInString++;
					}

					iPositionForOut++;
				}
			}
		}

		/////////////////////////////////
		inline void convert_dnaTokMer(const string& sDna, const readIDType& iID, const uint64_t& iNumberOfkMers, uint64_t& iPositionForOut, const int32_t& iMaxKTimes3, const bool& bRC, uint32_t iPositionInString, InputType<intType>& resultsVec) {
			// This value gives the remaining length of the read which is the number of kMers created
			const uint64_t& iMaxRange = iNumberOfkMers;
			if (iMaxRange >= 1) {
				// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
				const int32_t& iNumFrames = (iMaxRange >= 3) ? 3 : static_cast<int32_t>(iMaxRange);
				if (_bVisualize && _translatedFramesForVisualization.size() == 0) {
					for (int32_t i = 0; i < iNumFrames; ++i) {
						_translatedFramesForVisualization.push_back("");
					}
				}

				//int32_t ikMerCounter = 0;

				// Frameshifting, so that only one amino acid has to be computed
				// Compute initial frames
				AAFrames<intType> sAAFrames;

				for (int32_t j = 0; j < iNumFrames; ++j) {
					string sTempFrame = _sMaxKBlank;
					//if (_bSpaced) {
						//dnaToAminoacidSpaced(sDna, iMaxKTimes3, j, &sTempFrame);
					//}
					//else {
						dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
					//}
					if (_bVisualize) {
						_translatedFramesForVisualization[j] += sTempFrame;
					}
					sAAFrames[j] = aminoacidTokMer<intType>(sTempFrame);
					//cout << kMerToAminoacid(sAAFrames[j], 25) << " " << sTempFrame << " " << sDna << endl;

					auto kmerForSearchInTrie = sAAFrames[j];
					if (_bUnfunny) {
						kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
					}

					/*for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
						vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal]; //TODO: Slow
					}*/


					resultsVec.setkMer(iPositionForOut, kmerForSearchInTrie);
					resultsVec.setReadID(iPositionForOut, static_cast<uint32_t>(iID));
					if (GlobalInputParameters.bPostProcess) {
						resultsVec.setFrame(iPositionForOut, bRC);
						resultsVec.setPosition(iPositionForOut, iPositionInString);
						iPositionInString++;
					}

					iPositionForOut++;
				}

				//ikMerCounter += iNumFrames;

				if (iMaxRange > 3) {

					// If the Dna has an irregular (not mod 3) tail, a special case must be handled and an additional step in the loop is necessary thus the !(iMaxRange % 3)
					// To map the 2 onto 1, ! is applied twice
					const int32_t& iMaxRangeMod3 = static_cast<int32_t>(iMaxRange % 3);
					const int32_t& iMaxRangeMod3Neg = !(!iMaxRangeMod3);

					for (uint64_t j = 1; 3 * (j + iMaxRangeMod3Neg) < iMaxRange; j++) {
						for (int32_t k = 0; k < 3; ++k) {
							int8_t sTempAA = ' ';
							//if (_bSpaced) {
							//	dnaToAminoacidSpaced(sDna, k + iMaxKTimes3 + 3 * (j - 1), sTempAA);
							//}
							//else {
								dnaToAminoacid(sDna, k + iMaxKTimes3 + 3 * (j - 1), sTempAA);
							//}
							if (_bVisualize) {
								_translatedFramesForVisualization[k] += sTempAA;
							}

							sAAFrames[k] = aminoacidTokMer(sAAFrames[k], sTempAA);

							auto kmerForSearchInTrie = sAAFrames[k];
							if (_bUnfunny) {
								kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
							}

							/*for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
								vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal]; //TODO: Slow
							}*/

							resultsVec.setkMer(iPositionForOut, kmerForSearchInTrie);
							resultsVec.setReadID(iPositionForOut, static_cast<uint32_t>(iID));
							if (GlobalInputParameters.bPostProcess) {
								resultsVec.setFrame(iPositionForOut, bRC);
								resultsVec.setPosition(iPositionForOut, iPositionInString);
								iPositionInString++;
							}

							iPositionForOut++;
						}
						//ikMerCounter += 3;
					}

					// compute frames of said tail
					for (int32_t j = 0; j < iMaxRangeMod3; ++j) {
						int8_t sTempAA = ' ';
						//if (_bSpaced) {
						//	dnaToAminoacidSpaced(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
						//}
						//else {
							dnaToAminoacid(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
						//}
						if (_bVisualize) {
							_translatedFramesForVisualization[j] += sTempAA;
						}

						sAAFrames[j] = aminoacidTokMer(sAAFrames[j], sTempAA);


						auto kmerForSearchInTrie = sAAFrames[j];
						if (_bUnfunny) {
							kmerForSearchInTrie = aminoAcidsToAminoAcid(kmerForSearchInTrie);
						}

						/*for (int32_t kVal = 12 - _iMaxK; kVal < _iMaxK - _iMinK; ++kVal) {
							vCountGarbagekMerPerK[garbageRangeCalc + kVal] += (kmerForSearchInTrie & _tails[kVal]) == _tails[kVal]; //TODO: Slow
						}*/

						resultsVec.setkMer(iPositionForOut, kmerForSearchInTrie);
						resultsVec.setReadID(iPositionForOut, static_cast<uint32_t>(iID));
						if (GlobalInputParameters.bPostProcess) {
							resultsVec.setFrame(iPositionForOut, bRC);
							resultsVec.setPosition(iPositionForOut, iPositionInString);
							iPositionInString++;
						}

						iPositionForOut++;
					}
				}
			}
		}

		/////////////////////////////////
		inline void convert_dnaTokMerOneFrame(const string& sDna, const readIDType& iID, const uint64_t& iNumberOfkMers, uint64_t& iPositionForOut, uint32_t iPositionInString, InputType<intType>& resultsVec) {
			// This value gives the remaining length of the read which is the number of k-mers created
			const uint64_t& iMaxRange = iNumberOfkMers;

			if (iMaxRange > 0) {
				// go through the dna, convert it to an aminoacid string and then to its coded k-mer representation
				if (_bVisualize && _translatedFramesForVisualization.size() == 0) {
					_translatedFramesForVisualization.push_back("");
				}

				string sAA(sDna.length() / 3, ' ');
				dnaToAminoacid(sDna, static_cast<uint64_t>(sDna.length()), 0, &sAA);
				sAA = Utilities::rstrip(sAA, ' ');

				if (_bVisualize) {
					_translatedFramesForVisualization[0] += sAA;
				}

				//cout << sDna.length() << " " << sDna.length() / 3 << " " << sAA.length() << " " << iMaxRange << endl;

				for (uint64_t iCurrentkMerCounter = 0; iCurrentkMerCounter < iMaxRange; ++iCurrentkMerCounter) {
					auto kMer = aminoacidTokMer<intType>(sAA.cbegin() + iCurrentkMerCounter, sAA.cbegin() + iCurrentkMerCounter + _iHighestK);

					if (_bUnfunny) {
						kMer = aminoAcidsToAminoAcid(kMer);
					}

					resultsVec.setkMer(iPositionForOut, kMer);
					resultsVec.setReadID(iPositionForOut, static_cast<uint32_t>(iID));
					if (GlobalInputParameters.bPostProcess) {
						resultsVec.setFrame(iPositionForOut, 0);
						resultsVec.setPosition(iPositionForOut, iPositionInString);
						iPositionInString++;
					}

					iPositionForOut++;
				}
			}
		}

		/////////////////////////////////
		inline void convertLinesTokMers_new(const vector<tuple<string, readIDType, uint64_t, uint32_t>>& vLines, const vector<tuple<string, readIDType, uint64_t, uint32_t>>& vRCLines, const size_t& start, const size_t& end, uint64_t iPositionForOut, InputType<intType>& kMerVecOut) {

			const int32_t& iMaxKTimes3 = 3 * _iHighestK;
			//const uint64_t iStartPositionForOut = iPositionForOut;
			for (uint64_t i = start; i < end; ++i) {
				if (get<0>(vLines[i]) != "" || get<0>(vRCLines[i]) != "") {
					if (this->_bProtein) {
						convert_alreadyTranslatedTokMers(get<0>(vLines[i]), get<1>(vLines[i]), get<2>(vLines[i]), iPositionForOut, get<3>(vLines[i]), kMerVecOut);
					}
					else {
						if (_bOnlyOneFrame) {
							convert_dnaTokMerOneFrame(get<0>(vLines[i]), get<1>(vLines[i]), get<2>(vLines[i]), iPositionForOut, get<3>(vLines[i]), kMerVecOut);
						}
						else {
							convert_dnaTokMer(get<0>(vLines[i]), get<1>(vLines[i]), get<2>(vLines[i]), iPositionForOut, iMaxKTimes3, false, get<3>(vLines[i]), kMerVecOut);
							if (_bSixFrames) {
								convert_dnaTokMer(get<0>(vRCLines[i]), get<1>(vRCLines[i]), get<2>(vRCLines[i]), iPositionForOut, iMaxKTimes3, true, get<3>(vLines[i]), kMerVecOut);
							}
						}
					}
				}
			}

			
			//sort(kMerVecOut.begin() + iStartPositionForOut, kMerVecOut.begin() + iPositionForOut);
			/*std::sort(kMerVecOut.begin() + iStartPositionForOut, kMerVecOut.begin() + iPositionForOut, [](const tuple<uint64_t, intType, uint32_t, uint32_t>& p1, const tuple<uint64_t, intType, uint32_t, uint32_t>& p2) {
					return get<1>(p1) < get<1>(p2);
					});*/
			
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public:
		void translateFileInOneFrame(const string& sFilename, const string& sOutput) {
			ifstream inFile(sFilename);
			ofstream outFile(sOutput);
			size_t iQualiLength = 0;
			int32_t iWhatNext = 0;
			while (inFile) {
				string sDummyString = "";
				getline(inFile, sDummyString);
				if (sDummyString != "") {
					switch (iWhatNext) {
					case 0:
						outFile << sDummyString << "\n";
						iWhatNext = 1;
						break;
					case 1:
						{
							for (char& c : sDummyString) {
								if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'a' && c != 'c' && c != 'g' && c != 't') {
									c = 'Z';
								}
							}
							string sAA(sDummyString.length() / 3, ' ');
							dnaToAminoacid(sDummyString, sDummyString.length(), 0, &sAA);
							sAA = Utilities::rstrip(sAA, ' ');
							iQualiLength = sAA.length();
							outFile << sAA << "\n";
							iWhatNext = 2;
						}
						break;
					case 2:
						outFile << sDummyString << "\n";
						iWhatNext = 3;
						break;
					case 3:
						{
							outFile << string(iQualiLength, 'I') << "\n";
							iWhatNext = 0;
						}
						break;
					};
				}
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		struct strTransfer {
			string name = string(""), overhang = string(""), overhang2 = string("");
			size_t lengthOfDNA = 0;
			bool finished = true, addTail = false, bNewRead = true, bStart = true;
			uint8_t iExpectedInput = 0;
			string lastLine = string(""), lastLine2 = string("");
			//list<readIDType> vReadIDs;
			//readIDType iCurrentReadID = 0;
			//unordered_map<readIDType, uint64_t> mReadIDToArrayIdx;
			uint64_t iNumOfCharsRead = 0, iNumOfAllCharsRead = 0, iCurrentOverallPercentage = 0, iNumOfNewReads = 0, iCurrentLineInInfoFile = 0, iCurrentLineInInfoFile2 = 0;
			double  iCurrentPercentage = 0.0;
			uint32_t iPositionInRead = 0;
			//vector<uint64_t> vRangesOfOutVec;
		};

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	private:
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline uint64_t calculateMemoryUsageFromkMerCount(const uint64_t& kMerCount, const size_t& sizeOfEachElementOfOutVec) {
			if (!this->_bProtein && _bSixFrames) {
				return kMerCount * sizeOfEachElementOfOutVec * 2;
			}
				
			return kMerCount * sizeOfEachElementOfOutVec;
		}
		

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template<typename FileType>
		inline void readFileAndGenerateInfos(Utilities::FileReader<FileType>& input, fstream& output, const bool& bIsFasta, const size_t& sizeOfEachElementOfOutVec) {
			// TODO optimize this, fastq version
			/*
			1,5,1 - Name line, 5 dna parts, 1 chunk
			3,100,5 - + line & quality line & name line, 100 dna parts, 5 chunks with this being the 5/5
			0,100,4 - no skip line, 100 dna parts of same read, 5 chunks with this being the 4/5
			0,100,3 ...
			0,100,2 ...
			0,95,1 - no skip line, 95 dna parts, 5 chunks with this being the last one
			3,31,1
			*/

			Utilities::BufferedWriter outFileStream(output, 524288ull);
			if (bIsFasta) {
				pair<std::string, bool> resultChunkPair("", false);
				uint64_t iDummy = 0, iNumOfReadChars = 0, iNumOfCurrentlyReadChars = 0;
				uint64_t iSkippedLines = 0, iImportantParts = 0, iChunkNumber = 0;
				vector<uint64_t> iSavedLineCountsForChunks;

				while (!input.eof()) {
					input.getChunk(resultChunkPair, iNumOfCurrentlyReadChars);
					if (!resultChunkPair.first.empty()) {

						if (resultChunkPair.first.front() == '>') {
							// save info from last read
							if (iChunkNumber == 1) {
								Utilities::itostr(iSkippedLines, outFileStream.getString());
								outFileStream += ",";
								Utilities::itostr(iImportantParts, outFileStream.getString());
								outFileStream += ",";
								Utilities::itostr(iChunkNumber, outFileStream.getString());
								outFileStream += "\n";
								iImportantParts = 0; //iChunkNumber = 1;
								iNumOfReadChars = 0;
								iSavedLineCountsForChunks.clear();
							}
							else {
								iSavedLineCountsForChunks.push_back(iImportantParts);
								while (iChunkNumber >= 1) {
									Utilities::itostr(iSkippedLines, outFileStream.getString());
									outFileStream += ",";
									Utilities::itostr(iSavedLineCountsForChunks[iSavedLineCountsForChunks.size() - iChunkNumber], outFileStream.getString());
									outFileStream += ",";
									Utilities::itostr(iChunkNumber, outFileStream.getString());
									outFileStream += "\n";
									iSkippedLines = 0;
									--iChunkNumber;
								}
								iImportantParts = 0; iChunkNumber = 1;
								iNumOfReadChars = 0;
								iSavedLineCountsForChunks.clear();
							}


							// read name line completely
							while (!resultChunkPair.second) {
								resultChunkPair.first = "";
								input.getChunk(resultChunkPair, iDummy);
							}
							iSkippedLines = 1;
						}
						else {
							iImportantParts++;
							iNumOfReadChars += iNumOfCurrentlyReadChars;
							uint64_t iCurrentlyPredictedMemoryUsage = calculateMemoryUsageFromkMerCount(calculatekMerCount(iNumOfReadChars), sizeOfEachElementOfOutVec);
							if (iCurrentlyPredictedMemoryUsage > uint64_t(100) * 1024 * 1024) {
								iChunkNumber++;
								iSavedLineCountsForChunks.push_back(iImportantParts);
								iImportantParts = 0;
								iNumOfReadChars = 0;
							}
						}

						resultChunkPair.first = "";
						resultChunkPair.second = false;
					}
					else {
						iImportantParts++;
					}
				}

				// save info from last read
				if (iChunkNumber == 1) {
					Utilities::itostr(iSkippedLines, outFileStream.getString());
					outFileStream += ",";
					Utilities::itostr(iImportantParts, outFileStream.getString());
					outFileStream += ",";
					Utilities::itostr(iChunkNumber, outFileStream.getString());
					outFileStream += "\n";
				}
				else {
					iSavedLineCountsForChunks.push_back(iImportantParts);
					while (iChunkNumber >= 1) {
						Utilities::itostr(iSkippedLines, outFileStream.getString());
						outFileStream += ",";
						Utilities::itostr(iSavedLineCountsForChunks[iSavedLineCountsForChunks.size() - iChunkNumber], outFileStream.getString());
						outFileStream += ",";
						Utilities::itostr(iChunkNumber, outFileStream.getString());
						outFileStream += "\n";
						iSkippedLines = 0;
						--iChunkNumber;
					}
				}

			}
			else { // fastq
				// todo

				pair<std::string, bool> resultChunkPair("", false);
				uint64_t iDummy = 0, iNumOfReadChars = 0, iNumOfCurrentlyReadChars = 0, iNumberOfCharactersInDNALines = 0, iNumberOfCharactersInQualityLines = 0;
				uint64_t iSkippedLines = 0, iImportantParts = 0, iChunkNumber = 1;
				uint64_t iCurrentlyPredictedMemoryUsage = 0;
				vector<uint64_t> iSavedLineCountsForChunks;
				int32_t iSwitchForLineType = 0; // 0 name line, 1 dna, 2 + string, 3 quality string

				while (!input.eof()) {
					input.getChunk(resultChunkPair, iNumOfCurrentlyReadChars);
					if (!resultChunkPair.first.empty()) {

						if (resultChunkPair.second) {
							iNumOfCurrentlyReadChars--; // newline character is of no use to us
						}

						if (resultChunkPair.first.front() == '+' && iSwitchForLineType == 1) {
							iSwitchForLineType = 2;
						}

						switch (iSwitchForLineType) {
						//////////////////////////////////////////////////////////////
						case 0: // new read, read name line completely
							while (!resultChunkPair.second) {
								resultChunkPair.first = "";
								input.getChunk(resultChunkPair, iDummy);
							}
							iSkippedLines++;
							iSwitchForLineType = 1;
							break;
						//////////////////////////////////////////////////////////////
						case 1: // read dna
							iImportantParts++;
							iNumOfReadChars += iNumOfCurrentlyReadChars;
							iNumberOfCharactersInDNALines += iNumOfCurrentlyReadChars;
							iCurrentlyPredictedMemoryUsage = calculateMemoryUsageFromkMerCount(calculatekMerCount(iNumOfReadChars), sizeOfEachElementOfOutVec);
							if (iCurrentlyPredictedMemoryUsage > uint64_t(100) * 1024 * 1024) {
								iChunkNumber++;
								iSavedLineCountsForChunks.push_back(iImportantParts);
								iImportantParts = 0;
								iNumOfReadChars = 0;
							}
							break;
						//////////////////////////////////////////////////////////////
						case 2: // + line
							// save stuff
							if (iChunkNumber == 1) {
								Utilities::itostr(iSkippedLines, outFileStream.getString());
								outFileStream += ",";
								Utilities::itostr(iImportantParts, outFileStream.getString());
								outFileStream += ",";
								Utilities::itostr(iChunkNumber, outFileStream.getString());
								outFileStream += "\n";
								iImportantParts = 0; //iChunkNumber = 1;
								iNumOfReadChars = 0;
								iSavedLineCountsForChunks.clear();
							}
							else {
								iSavedLineCountsForChunks.push_back(iImportantParts);
								while (iChunkNumber >= 1) {
									Utilities::itostr(iSkippedLines, outFileStream.getString());
									outFileStream += ",";
									Utilities::itostr(iSavedLineCountsForChunks[iSavedLineCountsForChunks.size() - iChunkNumber], outFileStream.getString());
									outFileStream += ",";
									Utilities::itostr(iChunkNumber, outFileStream.getString());
									outFileStream += "\n";
									iSkippedLines = 0;
									--iChunkNumber;
								}
								iImportantParts = 0; iChunkNumber = 1;
								iNumOfReadChars = 0;
								iSavedLineCountsForChunks.clear();
							}

							// fully read + line
							while (!resultChunkPair.second) {
								resultChunkPair.first = "";
								input.getChunk(resultChunkPair, iDummy);
							}
							iSkippedLines = 1;
							iSwitchForLineType = 3;
							break;
						//////////////////////////////////////////////////////////////
						case 3: // parse quality, must have equal number of characters as dna
							iNumberOfCharactersInQualityLines += iNumOfCurrentlyReadChars;
							
							iDummy = 0;
							while (!resultChunkPair.second) {
								resultChunkPair.first = "";
								input.getChunk(resultChunkPair, iDummy);
								iNumberOfCharactersInQualityLines += iDummy;
							}
							if (iDummy > 0) { // loop was used
								iNumberOfCharactersInQualityLines--;
							}


							if (iNumberOfCharactersInQualityLines == iNumberOfCharactersInDNALines) {
								iNumberOfCharactersInDNALines = 0;
								iNumberOfCharactersInQualityLines = 0;
								iSwitchForLineType = 0;
							}
							if (iNumberOfCharactersInQualityLines > iNumberOfCharactersInDNALines) {
								throw runtime_error("Quality string and DNA string do not have the same length!");
							}

							iSkippedLines++;
							break;
						default:
							break;
						};

						resultChunkPair.first = "";
						resultChunkPair.second = false;
					}
					else {
						iImportantParts++;
					}
				}


				// last part of the file is usually bullshit, so mark it as such
				Utilities::itostr(iSkippedLines, outFileStream.getString());
				outFileStream += ",";
				Utilities::itostr(0, outFileStream.getString());
				outFileStream += ",";
				Utilities::itostr(0, outFileStream.getString());
				outFileStream += "\n";
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void putReadIntoLocalMemory(const string& sFalsekMerMarker, const readIDType& iLocalReadID, const string& sRead, vector<tuple<string, readIDType, uint64_t, uint32_t>>& vLines, vector<tuple<string, readIDType, uint64_t, uint32_t>>& vRCLines, uint64_t& iSumOfkMers, int64_t& iSoftMaxSize, const size_t& sizeOfEachElementOfOutVec) {
			if (!this->_bProtein && _bSixFrames) {
				string sTempRead = reverseComplement(sRead);
				sTempRead += sFalsekMerMarker;
				auto kMerCount = calculatekMerCount(sTempRead.length());
				iSumOfkMers += kMerCount;
				vRCLines.emplace_back(sTempRead, iLocalReadID, kMerCount, 0);

				iSoftMaxSize -= kMerCount * static_cast<int64_t>(sizeOfEachElementOfOutVec); // size of each element in the resulting vector that is compared with the index
				iSoftMaxSize -= sTempRead.length() * sizeof(char) + sizeof(readIDType) + sizeof(uint64_t) + sizeof(uint32_t); // size of an element in vLines/vRCLines
			}
			string sTempRead = sRead + sFalsekMerMarker;
			auto kMerCount = calculatekMerCount(sTempRead.length());
			iSumOfkMers += kMerCount;
			vLines.emplace_back(sTempRead, iLocalReadID, kMerCount, 0);

			iSoftMaxSize -= kMerCount * static_cast<int64_t>(sizeOfEachElementOfOutVec); 
			iSoftMaxSize -= sTempRead.length() * sizeof(char) + sizeof(readIDType) + sizeof(uint64_t) + sizeof(uint32_t);
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void paddingOfSmallReads(string& sRead, const size_t& sFalsekMerMarkerLength) {
			// pad very small reads with bullshit
			if (sRead.length() > 0) {
				if (this->_bProtein) {
					while (sRead.length() + sFalsekMerMarkerLength < size_t(_iHighestK)) {
						sRead += '^';
					}
				}
				else {
					if (_bOnlyOneFrame) {
						while ((sRead.length() + sFalsekMerMarkerLength) / 3 < size_t(_iHighestK)) {
							sRead += 'X';
						}
					}
					else {
						while (sRead.length() + sFalsekMerMarkerLength < size_t(_iHighestK) * 3) {
							sRead += 'X';
						}
					}
				}
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void searchAndReplaceLettersOfRead(string& sRead) {
			for (char& c : sRead) {
				if (c == '\t' || c == ' ') {
					throw runtime_error("Spaces or tabs inside read, please check your input.");
				}
				else {
					if (this->_bProtein) {
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
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void generateOverhang(const string& sRead, string& sOverhang) {
			if (this->_bProtein) {
				if (sRead.length() < size_t(_iHighestK)) {
					sOverhang = sRead;
				}
				else {
					sOverhang = sRead.substr(sRead.length() + 1 - _iHighestK);
				}
			}
			else {
				if (sRead.length() < size_t(_iHighestK) * 3) {
					sOverhang = sRead;
				}
				else {
					sOverhang = sRead.substr(sRead.length() + 1 - _iHighestK * 3);
				}
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template<typename FileType>
		inline void processInput(const vector<uint64_t>& entriesOfLine, const string& sFalsekMerMarker, Utilities::FileReader<FileType>& input, uint64_t& iNumberOfCharsRead, string& sReadName, uint64_t& iLengthOfRead, readIDType& iLocalReadID, vector<tuple<string, readIDType, uint64_t, uint32_t>>& vLines, vector<tuple<string, readIDType, uint64_t, uint32_t>>& vRCLines, uint64_t& iSumOfkMers, string& sOverhang, int64_t& iSoftMaxSize, const size_t& sizeOfEachElementOfOutVec) {
			pair<std::string, bool> pairOfInput("", false);
			
			// skip lines
			for (uint64_t iSkipLineIdx = 0; iSkipLineIdx < entriesOfLine[0]; ++iSkipLineIdx) {
				pairOfInput.first = ""; pairOfInput.second = false;
				uint64_t iNumOfCharsOfThatPart = 0;
				while (!pairOfInput.second) {
					input.getChunk(pairOfInput, iNumOfCharsOfThatPart);
					iNumberOfCharsRead += iNumOfCharsOfThatPart;
				}
			}
			if (entriesOfLine[0]) {
				pairOfInput.first.erase(pairOfInput.first.begin()); // In Fastq and Fasta, the last line before DNA ist the name of the read! Remove first character tho
				sReadName += pairOfInput.first + " ";
			}

			// read parts of content
			pairOfInput.first = ""; pairOfInput.second = false;
			for (uint64_t iReadLineIdx = 0; iReadLineIdx < entriesOfLine[1]; ++iReadLineIdx) {
				uint64_t iNumOfCharsOfThatPart = 0;

				//pairOfInput.first = "";

				input.getChunk(pairOfInput, iNumOfCharsOfThatPart);

				/*if (pairOfInput.first.find(' ') != string::npos) {
					cout << "Number of accesses : " << input.getiNumberOfAccesses() << endl;
					cout << iReadLineIdx << " " << pairOfInput.first << endl;
				}*/

				iNumberOfCharsRead += iNumOfCharsOfThatPart;
				iLengthOfRead += iNumOfCharsOfThatPart;
			}

			// clean read of illegal letters
			searchAndReplaceLettersOfRead(pairOfInput.first);

			// if this part is the contiuation of a former read, attach the overhanging information of a prior run to the string
			if (sOverhang != "") {
				pairOfInput.first = sOverhang + pairOfInput.first;
			}

			// if read is too small, append some bs characters
			paddingOfSmallReads(pairOfInput.first, sFalsekMerMarker.length());

			// save read into temporary memory for later convert and sort operations
			putReadIntoLocalMemory(sFalsekMerMarker, iLocalReadID, pairOfInput.first, vLines, vRCLines, iSumOfkMers, iSoftMaxSize, sizeOfEachElementOfOutVec);
			
			// current line belongs to a new read, so increase read ID counter for next read
			if (entriesOfLine[2] == 1) {
				++iLocalReadID;
			}

			// if the next line of the info file belongs to the same read, generate overhang for next iteration here
			if (entriesOfLine[2] > 1) {
				generateOverhang(pairOfInput.first, sOverhang);
			}
			else {
				sOverhang = "";
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void convertAndSort(const int64_t& iAvailMemory, const uint64_t& iSumOfkMers, InputType<intType>& vOut, vector<WorkerThread>& threadPool, const vector<tuple<string, readIDType, uint64_t, uint32_t>>& vLines, const vector<tuple<string, readIDType, uint64_t, uint32_t>>& vRCLines) {
			// convert it in parallel
			//auto startTIME = std::chrono::high_resolution_clock::now();
			try {
				// should the capacity be smaller than iSumOfkMers after the first allocation, bad_alloc will probably be thrown because there will likely be no larger contiguous chunk of memory available
				// this is mitigated by substracting 1% of the available memory once after this is called (see Compare.hpp "reduce available memory")
				vOut.reserve(iSumOfkMers);
				vOut.resize(iSumOfkMers);
			}
			catch (const bad_alloc&) {
				int64_t triedToAllocate = iSumOfkMers * vOut.sizeOf() / (1024 * 1024);
				cerr << "ERROR: Your system does not have enough contiguous memory available (which happens if the system is powered on over a long period of time). You might try to restart or use a lower number after -m. FYI: You tried to use " + to_string(triedToAllocate) + " MB but had only " << to_string(iAvailMemory / (1024ull * 1024ull)) << " MB available." << endl;
				throw;
			}

			const auto& chunkSize = iSumOfkMers / threadPool.size();
			auto chunkSizeOverhead = iSumOfkMers % threadPool.size();
			size_t start = 0;
			uint64_t iCurrentkMerCount = 0, iTotalkMerCount = 0;
			int32_t iProcID = 0;

#ifdef TIME
			endTIME = std::chrono::high_resolution_clock::now();
			cout << "Reserve memory_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
			startTIME = std::chrono::high_resolution_clock::now();
#endif


			for (size_t iLineIdx = 0; iLineIdx < vLines.size(); ++iLineIdx) {
				if (iCurrentkMerCount >= chunkSize + chunkSizeOverhead) {
					//call function with start, iLineIdx, iTotalkMerCount
					auto task = [&, start, iLineIdx, iTotalkMerCount, this](const int32_t&) { convertLinesTokMers_new(vLines, vRCLines, start, iLineIdx, iTotalkMerCount, ref(vOut)); };
					threadPool[iProcID].pushTask(task);

					start = iLineIdx;
					iTotalkMerCount += iCurrentkMerCount;
					iCurrentkMerCount = 0;
					chunkSizeOverhead = 0;
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

			// start Threads, join threads
			for (size_t iThreadID = 0; iThreadID < threadPool.size(); ++iThreadID) {
				threadPool[iThreadID].startThread();
			}
			for (size_t iThreadID = 0; iThreadID < threadPool.size(); ++iThreadID) {
				threadPool[iThreadID].waitUntilFinished();
			}

#ifdef TIME
			endTIME = std::chrono::high_resolution_clock::now();
			cout << "PAR Input translation_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
			startTIME = std::chrono::high_resolution_clock::now();
#endif
		}

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a paired-end fasta/q files
		template<typename FileType>
		inline uint64_t readFastqa_pairedEnd(Utilities::FileReader<FileType>& input, Utilities::FileReader<FileType>& input2, const string& sInfoFileName, const string& sInfoFileName2, InputType<intType>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, const bool& bIsFasta, unique_ptr<strTransfer>& transfer, vector<WorkerThread>& threadPool) {
			
			// Total number of kmers expected to be in the output vector
			uint64_t iSumOfkMers = 0;
			
			try {
				// Given memory at start
				const int64_t iAvailMemory = iSoftMaxSize / (1024 * 1024);

				// local memory of parsed reads
				vector<tuple<string, readIDType, uint64_t, uint32_t>> vLines, vRCLines;

				// appendix to each read in order to correctly generate k-mers from them
				string sFalsekMerMarker = "";
				if (this->_bProtein) {
					for (int32_t i = 0; i < (_iHighestK - _iMinK); ++i) { // (_iHighestK - _iMaxK) + (_iMaxK - _iMinK)
						sFalsekMerMarker += "^";
					}
				}
				else {
					for (int32_t i = 0; i < (_iHighestK - _iMinK) * 3; ++i) {
						sFalsekMerMarker += "X";
					}
				}

				// files containing the information about the lines of the input(s)
				fstream fileWithLineInfosFirst;
				fstream fileWithLineInfosSecond; 
				if (transfer->bStart) {
					fileWithLineInfosFirst.open(sInfoFileName, ios::out);
					fileWithLineInfosSecond.open(sInfoFileName2, ios::out);
					// Generate files 
					readFileAndGenerateInfos<FileType>(input, fileWithLineInfosFirst, bIsFasta, vOut.sizeOf());
					readFileAndGenerateInfos<FileType>(input2, fileWithLineInfosSecond, bIsFasta, vOut.sizeOf());

					// set to beginning again
					fileWithLineInfosFirst.close();
					fileWithLineInfosFirst.open(sInfoFileName, ios::in);
					fileWithLineInfosSecond.close();
					fileWithLineInfosSecond.open(sInfoFileName2, ios::in);
					input.resetAndCleanBuffer();
					input2.resetAndCleanBuffer();
					transfer->bStart = false;
				}
				else {
					fileWithLineInfosFirst.open(sInfoFileName, ios::in);
					fileWithLineInfosSecond.open(sInfoFileName2, ios::in);
					transfer->bStart = false;
					// continue where you left off, loop through file until the correct line
					for (uint64_t iLineCounter = 0; iLineCounter < transfer->iCurrentLineInInfoFile; ++iLineCounter) {
						string sDummyString = "";
						getline(fileWithLineInfosFirst, sDummyString);
					}

					for (uint64_t iLineCounter = 0; iLineCounter < transfer->iCurrentLineInInfoFile2; ++iLineCounter) {
						string sDummyString = "";
						getline(fileWithLineInfosSecond, sDummyString);
					}
				}

				// read info files, files and calculate memory usage until there is no space left
				bool bNotFull = true;
				uint64_t iLineCountOfInfoFile = transfer->iCurrentLineInInfoFile, iLineCountOfInfoFile2 = transfer->iCurrentLineInInfoFile2;
				readIDType iLocalReadID = 0, iLocalReadID2 = 0;
				uint64_t iNumberOfCharsRead = 0;
				string sOverhang = transfer->overhang, sOverhang2 = transfer->overhang2;
				bool bFirstFileIsOkay = true, bSecondFileIsOkay = true;
				string sReadName = transfer->name;
				uint64_t iLengthOfRead = transfer->lengthOfDNA;
				transfer->addTail = true;
				while (bNotFull) {
					// Show current percentage
					if (_bVerbose && iFileLength != 0) {
						double dPercentageOfInputRead = transfer->iNumOfCharsRead / double(iFileLength) * 100.;
						if (static_cast<uint64_t>(dPercentageOfInputRead) != static_cast<uint64_t>(transfer->iCurrentPercentage)) {
							transfer->iCurrentPercentage = dPercentageOfInputRead;
							cout << "OUT: Progress of current file " << static_cast<uint64_t>(transfer->iCurrentPercentage) << "%" << endl;
						}
						if (iFileLength != overallFilesSize) {
							dPercentageOfInputRead = transfer->iNumOfAllCharsRead / double(overallFilesSize) * 100.;
							if (static_cast<uint64_t>(dPercentageOfInputRead) != transfer->iCurrentOverallPercentage) {
								transfer->iCurrentOverallPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
								cout << "OUT: Progress of all files " << transfer->iCurrentOverallPercentage << " %" << endl;
							}
						}
					}

					// Read lines of info files
					string lineOfFirstInfoFile = "";
					vector<uint64_t> entriesOfLineFirst({ 0,0,0 });
					if (iLocalReadID <= iLocalReadID2) {
						if (getline(fileWithLineInfosFirst, lineOfFirstInfoFile)) {
							// content of line: number of lines to skip, number of getChunk calls to parse from filebuffer, memory friendly chunks (1 or more)
							if (lineOfFirstInfoFile != "") {
								entriesOfLineFirst = Utilities::splitToUInt64s(lineOfFirstInfoFile, ',');
							}
							++iLineCountOfInfoFile;
						}
						else {
							bFirstFileIsOkay = false;
						}
					}

					string lineOfSecondInfoFile = "";
					vector<uint64_t> entriesOfLineSecond({ 0,0,0 });
					if (iLocalReadID >= iLocalReadID2) {
						if (getline(fileWithLineInfosSecond, lineOfSecondInfoFile)) {
							if (lineOfSecondInfoFile != "") {
								entriesOfLineSecond = Utilities::splitToUInt64s(lineOfSecondInfoFile, ',');
							}
							++iLineCountOfInfoFile2;
						}
						else {
							bSecondFileIsOkay = false;
						}
					}

					// if the next chunk would cost too much memory, both files have been read, or the readIDs are more than an integer, end here for now
					if ((iSoftMaxSize <= static_cast<int64_t>(100*1024*1024)) || (!bFirstFileIsOkay && !bSecondFileIsOkay) || (iLocalReadID == numeric_limits<readIDType>::max()) || (iLocalReadID2 == numeric_limits<readIDType>::max())) {
						bNotFull = false;

						// save iLineCountOfInfoFile as new starting point in next call of this function
						iLineCountOfInfoFile--;
						iLineCountOfInfoFile2--;

						break;
					}
				
					// First file
					if (bFirstFileIsOkay && entriesOfLineFirst[2] > 0) {
						processInput<FileType>(entriesOfLineFirst, sFalsekMerMarker, input, iNumberOfCharsRead, sReadName, iLengthOfRead, iLocalReadID, vLines, vRCLines, iSumOfkMers, sOverhang, iSoftMaxSize, vOut.sizeOf());
					}

					// Second file
					if (bSecondFileIsOkay && entriesOfLineSecond[2] > 0) {
						processInput<FileType>(entriesOfLineSecond, sFalsekMerMarker, input2, iNumberOfCharsRead, sReadName, iLengthOfRead, iLocalReadID2, vLines, vRCLines, iSumOfkMers, sOverhang2, iSoftMaxSize, vOut.sizeOf());
					}

					// See if read is finished in both files
					if ((iLocalReadID == iLocalReadID2) && (entriesOfLineFirst[2] == 1 && entriesOfLineSecond[2] == 1)) {
						transfer->finished = true;
						transfer->addTail = false;
						if (bReadIDsAreInteresting) {
							if (sReadName != "" && iLengthOfRead != 0) {
								iSoftMaxSize -= sizeof(pair<string, uint32_t>) + sReadName.size() * sizeof(char) + sizeof(uint32_t); // size for vReadNameAndLength
								vReadNameAndLength.push_back(make_pair(sReadName, uint32_t(iLengthOfRead)));
								transfer->name = "";
								transfer->lengthOfDNA = 0;
							}
						}
						sReadName = "";
						iLengthOfRead = 0;
					}
					else {
						if (entriesOfLineFirst[2] != 0 || entriesOfLineSecond[2] != 0) {
							transfer->finished = false;
							transfer->addTail = true;
							if (bReadIDsAreInteresting) {
								if (sReadName != "" && iLengthOfRead != 0) {
									transfer->name = sReadName;
									transfer->lengthOfDNA += iLengthOfRead;
								}
							}
						}
					}

					// Calculate memory usage when info about reads must be saved
					if (bReadIDsAreInteresting && transfer->finished) {
						if (GlobalInputParameters.bPostProcess) {
							iSoftMaxSize -= sizeof(float);
						}
						iSoftMaxSize -= iAmountOfSpecies * sizeof(float); // for matrix with results of comparison later
					}

					// Save more stuff
					transfer->iNumOfCharsRead += iNumberOfCharsRead;
					transfer->iNumOfAllCharsRead += iNumberOfCharsRead;
					iNumberOfCharsRead = 0;
				}


				// process data
				if (vLines.size()) {
					convertAndSort(iAvailMemory, iSumOfkMers, vOut, threadPool, vLines, vRCLines);
				}

				// save current state
				if (bFirstFileIsOkay || bSecondFileIsOkay) {
					transfer->overhang = sOverhang;
					transfer->overhang2 = sOverhang2;
					transfer->iCurrentLineInInfoFile = iLineCountOfInfoFile;
					transfer->iCurrentLineInInfoFile2 = iLineCountOfInfoFile2;
				}

				transfer->iNumOfNewReads = max(iLocalReadID, iLocalReadID2) + transfer->addTail; // save total number of reads

				// if there is still data, that's fine but if there is not and the eof flag hasn't been set yet, it should be set
				if (bNotFull || (!bFirstFileIsOkay && !bSecondFileIsOkay)) {
					pair<std::string, bool> pairOfInput("", false);
					uint64_t iNumOfCharsOfThatChunk = 0;
					while (!input.eof()) {
						pairOfInput.first = "";
						input.getChunk(pairOfInput, iNumOfCharsOfThatChunk);
					}
					while (!input2.eof()) {
						pairOfInput.first = "";
						input2.getChunk(pairOfInput, iNumOfCharsOfThatChunk);
					}
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}

			return iSumOfkMers;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read one fasta/q file
		template<typename FileType>
		inline uint64_t readFastqa_singleEnd(Utilities::FileReader<FileType>& input, const string& sInfoFileName, InputType<intType>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, const bool& bIsFasta, unique_ptr<strTransfer>& transfer, vector<WorkerThread>& threadPool) {

			// Total number of kmers expected to be in the output vector
			uint64_t iSumOfkMers = 0;

			try {

				// Given memory at start
				const int64_t iAvailMemory = iSoftMaxSize / (1024 * 1024);

				// local memory of parsed reads
				vector<tuple<string, readIDType, uint64_t, uint32_t>> vLines, vRCLines;

				// appendix to each read in order to correctly generate k-mers from them
				string sFalsekMerMarker = "";
				if (this->_bProtein) {
					for (int32_t i = 0; i < (_iHighestK - _iMinK); ++i) { // (_iHighestK - _iMaxK) + (_iMaxK - _iMinK)
						sFalsekMerMarker += "^";
					}
				}
				else {
					for (int32_t i = 0; i < (_iHighestK - _iMinK) * 3; ++i) {
						sFalsekMerMarker += "X";
					}
				}

				// files containing the information about the lines of the input(s)
				fstream fileWithLineInfosFirst;
				fstream fileWithLineInfosSecond;
				if (transfer->bStart) {
					fileWithLineInfosFirst.open(sInfoFileName, ios::out);
					// Generate files 
					readFileAndGenerateInfos<FileType>(input, fileWithLineInfosFirst, bIsFasta, vOut.sizeOf());

					// set to beginning again
					fileWithLineInfosFirst.close();
					fileWithLineInfosFirst.open(sInfoFileName, ios::in);
					input.resetAndCleanBuffer();
					transfer->bStart = false;
				}
				else {
					fileWithLineInfosFirst.open(sInfoFileName, ios::in);
					transfer->bStart = false;
					// continue where you left off, loop through file until the correct line
					for (uint64_t iLineCounter = 0; iLineCounter < transfer->iCurrentLineInInfoFile; ++iLineCounter) {
						string sDummyString = "";
						getline(fileWithLineInfosFirst, sDummyString);
					}
				}

				// read info files, files and calculate memory usage until there is no space left
				bool bNotFull = true;
				uint64_t iLineCountOfInfoFile = transfer->iCurrentLineInInfoFile;
				readIDType iLocalReadID = 0;
				uint64_t iNumberOfCharsRead = 0;
				string sOverhang = transfer->overhang;
				bool bFirstFileIsOkay = true;
				string sReadName = transfer->name;
				uint64_t iLengthOfRead = transfer->lengthOfDNA;
				transfer->addTail = true;
				while (bNotFull) {
					// Show current percentage
					if (_bVerbose && iFileLength != 0) {
						double dPercentageOfInputRead = transfer->iNumOfCharsRead / double(iFileLength) * 100.;
						if (static_cast<uint64_t>(dPercentageOfInputRead) != static_cast<uint64_t>(transfer->iCurrentPercentage)) {
							transfer->iCurrentPercentage = dPercentageOfInputRead;
							cout << "OUT: Progress of current file " << static_cast<uint64_t>(transfer->iCurrentPercentage) << "%" << endl;
						}
						if (iFileLength != overallFilesSize) {
							dPercentageOfInputRead = transfer->iNumOfAllCharsRead / double(overallFilesSize) * 100.;
							if (static_cast<uint64_t>(dPercentageOfInputRead) != transfer->iCurrentOverallPercentage) {
								transfer->iCurrentOverallPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
								cout << "OUT: Progress of all files " << transfer->iCurrentOverallPercentage << " %" << endl;
							}
						}
					}

					// Read lines of info files
					string lineOfFirstInfoFile = "";
					vector<uint64_t> entriesOfLineFirst({ 0,0,0 });

					if (getline(fileWithLineInfosFirst, lineOfFirstInfoFile)) {
						// content of line: number of lines to skip, number of getChunk calls to parse from filebuffer, memory friendly chunks (1 or more)
						if (lineOfFirstInfoFile != "") {
							entriesOfLineFirst = Utilities::splitToUInt64s(lineOfFirstInfoFile, ',');
						}
						++iLineCountOfInfoFile;
					}
					else {
						bFirstFileIsOkay = false;
						}

					// if the next chunk would cost too much memory, both files have been read, or the readIDs are more than an integer, end here for now
					if ((iSoftMaxSize <= static_cast<int64_t>(100 * 1024 * 1024)) || !bFirstFileIsOkay || (iLocalReadID == numeric_limits<readIDType>::max())) {
						bNotFull = false;

						// save iLineCountOfInfoFile as new starting point in next call of this function
						iLineCountOfInfoFile--;

						break;
					}

					// Process file
					if (bFirstFileIsOkay && entriesOfLineFirst[2] > 0) {
						processInput<FileType>(entriesOfLineFirst, sFalsekMerMarker, input, iNumberOfCharsRead, sReadName, iLengthOfRead, iLocalReadID, vLines, vRCLines, iSumOfkMers, sOverhang, iSoftMaxSize, vOut.sizeOf());
					}

					// See if read is finished
					if (entriesOfLineFirst[2] == 1) {
						transfer->finished = true;
						transfer->addTail = false;
						if (bReadIDsAreInteresting) {
							if (sReadName != "" && iLengthOfRead != 0) {
								iSoftMaxSize -= sizeof(pair<string, uint32_t>) + sReadName.size() * sizeof(char) + sizeof(uint32_t); // size for vReadNameAndLength
								vReadNameAndLength.push_back(make_pair(sReadName, uint32_t(iLengthOfRead)));
								transfer->name = "";
								transfer->lengthOfDNA = 0;
							}
						}
						sReadName = "";
						iLengthOfRead = 0;
					}
					else {
						if (entriesOfLineFirst[2] != 0) {
							transfer->finished = false;
							transfer->addTail = true;
							if (bReadIDsAreInteresting) {
								if (sReadName != "" && iLengthOfRead != 0) {
									transfer->name = sReadName;
									transfer->lengthOfDNA += iLengthOfRead;
								}
							}
						}
					}

					// Calculate memory usage when info about reads must be saved
					if (bReadIDsAreInteresting && transfer->finished) {
						if (GlobalInputParameters.bPostProcess) {
							iSoftMaxSize -= sizeof(float);
						}
						iSoftMaxSize -= iAmountOfSpecies * sizeof(float); // for matrix with results of comparison later
					}

					// Save more stuff
					transfer->iNumOfCharsRead += iNumberOfCharsRead;
					transfer->iNumOfAllCharsRead += iNumberOfCharsRead;
					iNumberOfCharsRead = 0;
				}


				// process data
				if (vLines.size()) {
					convertAndSort(iAvailMemory, iSumOfkMers, vOut, threadPool, vLines, vRCLines);
				}

				// save current state
				if (bFirstFileIsOkay) {
					transfer->overhang = sOverhang;
					transfer->iCurrentLineInInfoFile = iLineCountOfInfoFile;
				}

				transfer->iNumOfNewReads = iLocalReadID + transfer->addTail; // save total number of reads

				// if there is still data, that's fine but if there is not and the eof flag hasn't been set yet, it should be set
				if (bNotFull || !bFirstFileIsOkay) {
					pair<std::string, bool> pairOfInput("", false);
					uint64_t iNumOfCharsOfThatChunk = 0;
					while (!input.eof()) {
						pairOfInput.first = "";
						input.getChunk(pairOfInput, iNumOfCharsOfThatChunk);
					}
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}

			return iSumOfkMers;
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fastq as input for comparison - old version, not working anymore!
		template<typename T>
		inline uint64_t readFastqa_partialSort(Utilities::FileReader<T>& input, Utilities::FileReader<T>& input2, InputType<intType>& vOut, list<pair<string, uint32_t>>& vReadNameAndLength, int64_t iSoftMaxSize, const uint64_t& iAmountOfSpecies, const uint64_t& iFileLength, const size_t& overallFilesSize, const bool& bReadIDsAreInteresting, const bool& bIsFasta, unique_ptr<strTransfer>& transfer, vector<WorkerThread>& threadPool) {
			
			vector<tuple<string, readIDType, uint64_t, uint32_t>> vLines, vRCLines, vLines2, vRCLines2;

			uint32_t iStartingPositionInRead = transfer->iPositionInRead;
			uint64_t iSumOfkMers = 0;
			const int64_t iAvailMemory = iSoftMaxSize / (1024 * 1024);
			try {
//#define TIME
#ifdef TIME
				auto startTIME = std::chrono::high_resolution_clock::now();
#endif

				bool bNotFull = true;
				readIDType iLocalReadID = 0;

				string sFalsekMerMarker = "";
				if (this->_bProtein) {
					for (int32_t i = 0; i < (_iHighestK - _iMinK); ++i) { // (_iHighestK - _iMaxK) + (_iMaxK - _iMinK)
						sFalsekMerMarker += "^";
					}
				}
				else {
					for (int32_t i = 0; i < (_iHighestK - _iMinK) * 3; ++i) {
						sFalsekMerMarker += "X";
					}
				}

				////////////////////////////////////////////////////////
				uint8_t iExpectedInput = transfer->iExpectedInput; // 0 = name, 1 = sequence, 2 = + or quality
				string sName = transfer->name, sDNA = "", sOverhang = transfer->overhang, sDNA2 = "", sOverhang2 = transfer->overhang2;
				size_t iDNALength = transfer->lengthOfDNA, iQualityLength = 0;
				bool bNewRead = transfer->bNewRead, bAddTail = transfer->addTail;

				if (sOverhang == "" && sName != "") {
					if (bReadIDsAreInteresting) {
						Utilities::checkOverFlow(iLocalReadID, __FILE__, __LINE__);
						++iLocalReadID;
						vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));
						iSoftMaxSize -= sizeof(pair<string, uint32_t>) + sName.size() * sizeof(char) + sizeof(uint32_t);
						if (GlobalInputParameters.bPostProcess) {
							iSoftMaxSize -= sizeof(float);
						}
					}
				}

				transfer->iNumOfNewReads = 0 + (transfer->addTail);
				transfer->finished = false;
				pair<std::string, bool> resultChunkPair("", false), resultChunkPair2("", false);
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

					if (resultChunkPair.first.length() == 0) {
						if (input2.notNull()) {
							if (resultChunkPair2.first.length() == 0) {
								continue;
							}
							else {

							}
						}
						else {
							continue;
						}
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
						}
						if (input2.notNull()) {
							if (!resultChunkPair2.second) {
								input2.ignore(iNumOfChars);
							}
						}
						iExpectedInput = 2;
						continue;
					}

					if (iExpectedInput == 1 && (resultChunkPair.first.length() < static_cast<size_t>((this->_bProtein) ? (_iHighestK + 1) : (3 * _iHighestK + 1))) && resultChunkPair.second == false) {
						// currently read DNA string is too short and there is more DNA to be had, usually happens when the end of the buffer is read but the file still contains DNA
						// first, count number of characters read in the first try and then append more DNA
						transfer->iNumOfCharsRead += iNumOfChars;
						transfer->iNumOfAllCharsRead += iNumOfChars;
						input.getChunk(resultChunkPair, iNumOfChars); // new DNA will be appended
						if (input2.notNull()) {
							if (resultChunkPair2.second == false) {
								input2.getChunk(resultChunkPair2, iNumOfChars); // new DNA will be appended
							}
						}
					}
					if (input2.notNull()) {
						// vice versa for other pair
						if (iExpectedInput == 1 && (resultChunkPair2.first.length() < static_cast<size_t>((this->_bProtein) ? (_iHighestK + 1) : (3 * _iHighestK + 1))) && resultChunkPair2.second == false) {
							transfer->iNumOfCharsRead += iNumOfChars;
							transfer->iNumOfAllCharsRead += iNumOfChars;
							input2.getChunk(resultChunkPair2, iNumOfChars); // new DNA will be appended

							if (resultChunkPair.second == false) {
								input.getChunk(resultChunkPair, iNumOfChars); // new DNA will be appended
							}
						}
					}

					// if the paired-ends are not of similar lengths, it can happen that one read was read just fine, but the other wasn't.
					if (input2.notNull()) {
						if (iExpectedInput == 1) {
							if (resultChunkPair.second != resultChunkPair2.second) {
								if (!resultChunkPair.second) {
									transfer->iNumOfCharsRead += iNumOfChars;
									transfer->iNumOfAllCharsRead += iNumOfChars;
									input.getChunk(resultChunkPair, iNumOfChars); // new DNA will be appended
								}
								else {
									input2.getChunk(resultChunkPair2, iNumOfChars); // new DNA will be appended
								}
							}
						}
					}

					if (_bVerbose && iFileLength != 0) {
						transfer->iNumOfCharsRead += iNumOfChars;
						transfer->iNumOfAllCharsRead += iNumOfChars;
						double dPercentageOfInputRead = transfer->iNumOfCharsRead / double(iFileLength) * 100.;
						if (static_cast<uint64_t>(dPercentageOfInputRead) != static_cast<uint64_t>(transfer->iCurrentPercentage)) {
							transfer->iCurrentPercentage = dPercentageOfInputRead;
							cout << "OUT: Progress of current file " << static_cast<uint64_t>(transfer->iCurrentPercentage) << "%" << endl;
						}
						if (iFileLength != overallFilesSize) {
							dPercentageOfInputRead = transfer->iNumOfAllCharsRead / double(overallFilesSize) * 100.;
							if (static_cast<uint64_t>(dPercentageOfInputRead) != transfer->iCurrentOverallPercentage) {
								transfer->iCurrentOverallPercentage = static_cast<uint64_t>(dPercentageOfInputRead);
								cout << "OUT: Progress of all files " << transfer->iCurrentOverallPercentage << " %" << endl;
							}
						}
					}

					if (bIsFasta && resultChunkPair.first.front() == '>') {
						if (bAddTail) {
							// Pad very small reads
							auto padding = [&](vector<tuple<string, readIDType, uint64_t, uint32_t>>& vLines) {
								if (vLines.size()) {
									if (this->_bProtein) {
										while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
											get<0>(vLines.back()) += '^';
										}
									}
									else {
										if (_bOnlyOneFrame) {
											while ((get<0>(vLines.back()).length() + sFalsekMerMarker.length()) / 3 < size_t(_iHighestK)) {
												get<0>(vLines.back()) += 'X';
											}
										}
										else {
											while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
												get<0>(vLines.back()) += 'X';
											}
										}
									}

									get<0>(vLines.back()) += sFalsekMerMarker; // add false marker

									iSumOfkMers -= get<2>(vLines.back());
									iSoftMaxSize += get<2>(vLines.back()) * vOut.sizeOf();

									get<2>(vLines.back()) = calculatekMerCount(get<0>(vLines.back()).length());

									iSumOfkMers += get<2>(vLines.back());
									iSoftMaxSize -= get<2>(vLines.back()) * vOut.sizeOf();
								}
							};
							padding(vLines);
							if (input2.notNull()) {
								padding(vLines2);
							}

							if (iSoftMaxSize <= static_cast<int64_t>(14399760 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 + 4 bytes of memory
								bNotFull = false;
							}

							if (bReadIDsAreInteresting) {
								//++(transfer->iCurrentReadID);
								transfer->finished = true;
								Utilities::checkOverFlow(iLocalReadID, __FILE__, __LINE__);
								++iLocalReadID;
								vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));// + sFalsekMerMarker.length())));
								iSoftMaxSize -= sizeof(pair<string, uint32_t>) + sName.size() * sizeof(char) + sizeof(uint32_t);
								if (GlobalInputParameters.bPostProcess) {
									iSoftMaxSize -= sizeof(float);
								}
							}
							iStartingPositionInRead = 0;
							bAddTail = false;
						}

						iExpectedInput = 0;
					}

					const uint8_t iCurrentExpInput = iExpectedInput;
					switch (iCurrentExpInput) {
						///////////////////////////////////////////////// case 0
					case 0: // New read
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
						}
						if (input2.notNull()) {
							while (!resultChunkPair2.second) {
								resultChunkPair2.first = "";
								input2.getChunk(resultChunkPair2, iNumOfChars); // discard the rest of the name
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
					case 1: // DNA / Protein sequence
					{
						auto searchAndReplaceLetters = [&](pair<std::string, bool>& resultChunkPair) {
							for (char& c : resultChunkPair.first) {
								if (c == '\t' || c == ' ') {
									throw runtime_error("Spaces or tabs inside read, please check your input. Error occured in:\n" + sName);
								}
								else {
									if (this->_bProtein) {
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
							if (this->_bProtein) {
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
									if (sDNA.length() < size_t(_iHighestK) * 3 && resultChunkPair.second) { // sequence was fully read
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

						auto emplaceBack = [&](const string& sDNA, vector<tuple<string, readIDType, uint64_t, uint32_t>>& vLines, vector<tuple<string, readIDType, uint64_t, uint32_t>>& vRCLines) {
							if (!bReadTooShort) {
								if (bNewRead) {
									bNewRead = false;
									if (!this->_bProtein && _bSixFrames) {
										string tempDNA = reverseComplement(sDNA);
										while (tempDNA.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
											tempDNA += 'X';
										}
										tempDNA += sFalsekMerMarker;
										vRCLines.emplace_back(tempDNA, iLocalReadID, calculatekMerCount(tempDNA.length()), iStartingPositionInRead);
										iSoftMaxSize -= tempDNA.length() * sizeof(char) + sizeof(readIDType) + sizeof(int32_t);
									}
								}
								else {
									if (!this->_bProtein && _bSixFrames) {
										const auto& rc = reverseComplement(sDNA);
										vRCLines.emplace_back(rc, iLocalReadID, calculatekMerCount(rc.length()), iStartingPositionInRead);
										iSoftMaxSize -= rc.length() * sizeof(char) + sizeof(readIDType) + sizeof(int32_t);
									}
								}

								vLines.emplace_back(sDNA, iLocalReadID, calculatekMerCount(sDNA.length()), iStartingPositionInRead);
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
						iStartingPositionInRead += calculatekMerCount(sDNA.length());

						break;
					}
					///////////////////////////////////////////////// case 2
					case 2: // Score
					{
						if (vLines.size() == 0) {
							if (transfer->lastLine != "") {
								auto lastLineKmerSize = calculatekMerCount(transfer->lastLine.length());
								vLines.emplace_back(transfer->lastLine, 0ul, lastLineKmerSize, iStartingPositionInRead);
								if (_bSixFrames) {
									vRCLines.emplace_back("", 0ul, 0, 0);
								}

								if (input2.notNull()) {
									auto lastLineKmerSize2 = calculatekMerCount(transfer->lastLine2.length());
									vLines2.emplace_back(transfer->lastLine2, 0ul, lastLineKmerSize2, iStartingPositionInRead);
									if (_bSixFrames) {
										vRCLines2.emplace_back("", 0ul, 0, 0);
									}
									lastLineKmerSize += lastLineKmerSize2;
								}

								iSumOfkMers += lastLineKmerSize;
								iSoftMaxSize -= lastLineKmerSize * vOut.sizeOf();
								if (iSoftMaxSize <= static_cast<int64_t>(14399760 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 +4 bytes of memory
									bNotFull = false;
								}
							}
							//else {

								//throw runtime_error("No overhang line was found, something went wrong while parsing the input!");
							//}
						}


						if (bAddTail ) { //&& transfer->lastLine != ""
							// Pad very small reads
							auto padding = [&](vector<tuple<string, readIDType, uint64_t, uint32_t>>& vLines) {
								if (vLines.size()) {
									if (this->_bProtein) {
										while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
											get<0>(vLines.back()) += '^';
										}
									}
									else {
										if (_bOnlyOneFrame) {
											while ((get<0>(vLines.back()).length() + sFalsekMerMarker.length()) / 3 < size_t(_iHighestK)) {
												get<0>(vLines.back()) += 'X';
											}
										}
										else {
											while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
												get<0>(vLines.back()) += 'X';
											}
										}
									}

									get<0>(vLines.back()) += sFalsekMerMarker; // add false marker

									iSumOfkMers -= get<2>(vLines.back());
									iSoftMaxSize += get<2>(vLines.back()) * vOut.sizeOf();

									get<2>(vLines.back()) = calculatekMerCount(get<0>(vLines.back()).length());

									iSumOfkMers += get<2>(vLines.back());
									iSoftMaxSize -= get<2>(vLines.back()) * vOut.sizeOf();
								}
							};
							padding(vLines);
							if (input2.notNull()) {
								padding(vLines2);
							}

							if (iSoftMaxSize <= static_cast<int64_t>(14399760 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 + 4 bytes of memory
								bNotFull = false;
							}

							if (bReadIDsAreInteresting) {
								//++(transfer->iCurrentReadID);
								Utilities::checkOverFlow(iLocalReadID, __FILE__, __LINE__);
								++iLocalReadID;
								vReadNameAndLength.push_back(make_pair(sName, uint32_t(iDNALength)));// + sFalsekMerMarker.length())));
								iSoftMaxSize -= sizeof(pair<string, uint32_t>) + sName.size() * sizeof(char) + sizeof(uint32_t);
								if (GlobalInputParameters.bPostProcess) {
									iSoftMaxSize -= sizeof(float);
								}
								transfer->finished = true;
							}
							iStartingPositionInRead = 0;
							bAddTail = false;
						}

						iQualityLength += resultChunkPair.first.length();
						if (input2.notNull()) {
							iQualityLength += resultChunkPair2.first.length();
						}
						while (iQualityLength < iDNALength) {
							if (!resultChunkPair.second) {
								resultChunkPair.first = "";
								iNumOfChars = 0;
								input.getChunk(resultChunkPair, iNumOfChars); // get the whole quality string
								transfer->iNumOfCharsRead += iNumOfChars;
								transfer->iNumOfAllCharsRead += iNumOfChars;
								iQualityLength += resultChunkPair.first.length();
							}
							if (input2.notNull() && !resultChunkPair2.second) {
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

						iSoftMaxSize -= iNumOfkMers * vOut.sizeOf();

						if (input2.notNull()) {
							const auto& kMersForward2 = (vLines2.size()) ? get<2>(vLines2.back()) : 0;
							const auto& kMersBackwards2 = (vRCLines2.size()) ? get<2>(vRCLines2.back()) : 0;
							const auto& iNumOfkMers2 = kMersForward2 + kMersBackwards2;
							iSumOfkMers += iNumOfkMers2;

							iSoftMaxSize -= iNumOfkMers2 * vOut.sizeOf();
						}

						if (iSoftMaxSize <= static_cast<int64_t>(14399756 + 4 * iAmountOfSpecies)) { // next chunk would at most need (100033 - 12 * 3 + 1) * 6 * 24 + 4 * iAmountOfSpecies + 40 + 4 bytes of memory
							bNotFull = false;
						}
					}
				}

				if (input.eof()) {
					if (bIsFasta) {
						auto handleEOF = [&](const string& sDNA, vector<tuple<string, readIDType, uint64_t, uint32_t>>& vLines, vector<tuple<string, readIDType, uint64_t, uint32_t>>& vRCLines) {
							if (bReadTooShort) {
								if (bNewRead) {
									bNewRead = false;
									if (!this->_bProtein && _bSixFrames) {
										string tempDNA = reverseComplement(sDNA);
										while (tempDNA.length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
											tempDNA += 'X';
										}
										tempDNA += sFalsekMerMarker;
										vRCLines.emplace_back(tempDNA, iLocalReadID, calculatekMerCount(tempDNA.length()), iStartingPositionInRead);
										iSumOfkMers += get<2>(vRCLines.back());
									}
								}
								else {
									if (!this->_bProtein && _bSixFrames) {
										const auto& rc = reverseComplement(sDNA);
										vRCLines.emplace_back(rc, iLocalReadID, calculatekMerCount(rc.length()), iStartingPositionInRead);
										iSumOfkMers += get<2>(vRCLines.back());
									}
								}

								vLines.emplace_back(sDNA, iLocalReadID, calculatekMerCount(sDNA.length()), iStartingPositionInRead);
							}

							if (this->_bProtein) {
								while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK)) {
									get<0>(vLines.back()) += '^';
								}
							}
							else {
								if (_bOnlyOneFrame) {
									while ((get<0>(vLines.back()).length() + sFalsekMerMarker.length()) / 3 < size_t(_iHighestK)) {
										get<0>(vLines.back()) += 'X';
									}
								}
								else {
									while (get<0>(vLines.back()).length() + sFalsekMerMarker.length() < size_t(_iHighestK) * 3) {
										get<0>(vLines.back()) += 'X';
									}
								}
							}

							get<0>(vLines.back()) += sFalsekMerMarker; // add false marker

							iSumOfkMers -= get<2>(vLines.back());
							iSoftMaxSize += get<2>(vLines.back()) * vOut.sizeOf();

							get<2>(vLines.back()) = calculatekMerCount(get<0>(vLines.back()).length());

							iSumOfkMers += get<2>(vLines.back());
							iSoftMaxSize -= get<2>(vLines.back()) * vOut.sizeOf();
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
					transfer->iPositionInRead = 0;
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
					if (iExpectedInput != 0) {
						transfer->iPositionInRead = iStartingPositionInRead;
					}

					if (sOverhang.length() == 0 && !input.eof()) {
						cerr << "Info: Overhanging line is empty!" << endl;
						if (bReadIDsAreInteresting) {
							iLocalReadID--;
							vReadNameAndLength.pop_back();
						}
					}
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
					vLines2.shrink_to_fit();
					vRCLines.insert(vRCLines.end(), vRCLines2.begin(), vRCLines2.end());
					vRCLines2.clear();
					vRCLines2.shrink_to_fit();
				}


				// convert it in parallel
				//auto startTIME = std::chrono::high_resolution_clock::now();
				try {
					// should the capacity be smaller than iSumOfkMers after the first allocation, bad_alloc will probably be thrown because there will likely be no larger contiguous chunk of memory available
					// this is mitigated by substracting 1% of the available memory once after this is called (see Compare.hpp "reduce available memory")
					vOut.reserve(iSumOfkMers); 
					vOut.resize(iSumOfkMers);
				}
				catch (const bad_alloc&) {
					int64_t triedToAllocate = iSumOfkMers * vOut.sizeOf() / (1024 * 1024);
					cerr << "ERROR: Your system does not have enough contiguous memory available (which happens if the system is powered on over a long period of time). You might try to restart or use a lower number after -m. FYI: You tried to use " + to_string(triedToAllocate) + " MB but had only " << to_string(iAvailMemory/(1024ull*1024ull)) << " MB available." << endl;	
					throw;
				}

				const auto& chunkSize = iSumOfkMers / threadPool.size();
				auto chunkSizeOverhead = iSumOfkMers % threadPool.size();
				size_t start = 0;
				uint64_t iCurrentkMerCount = 0, iTotalkMerCount = 0;
				int32_t iProcID = 0;

#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "Reserve memory_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				//vector<uint64_t> vRangesOfOutVec(threadPool.size() + 1);
				//cout << "Memory left: " << iSoftMaxSize << endl;
				for (size_t iLineIdx = 0; iLineIdx < vLines.size(); ++iLineIdx) {
					if (iCurrentkMerCount >= chunkSize + chunkSizeOverhead) {
						//call function with start, iLineIdx, iTotalkMerCount
						auto task = [&, start, iLineIdx, iTotalkMerCount, this](const int32_t&) { convertLinesTokMers_new(vLines, vRCLines, start, iLineIdx, iTotalkMerCount, ref(vOut)); };
						threadPool[iProcID].pushTask(task);
						//vRangesOfOutVec[iProcID] = iTotalkMerCount;
						//cout << iTotalkMerCount << endl;
						start = iLineIdx;
						iTotalkMerCount += iCurrentkMerCount;
						iCurrentkMerCount = 0;
						chunkSizeOverhead = 0;
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
				for (size_t iThreadID = 0; iThreadID < threadPool.size(); ++iThreadID) {
					threadPool[iThreadID].startThread();
				}
				for (size_t iThreadID = 0; iThreadID < threadPool.size(); ++iThreadID) {
					threadPool[iThreadID].waitUntilFinished();
				}

#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "PAR Input translation_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				// merge sorted ranges. If number of threads is odd, one extra run must be done at the end
				// even part
				/*for (uint32_t j = 1; j < threadPool.size(); j *= 2) {
					for (uint32_t i = 0; i < threadPool.size() / j; i += 2) {
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
		inline void dnaTokMers(const string& sDna, const uint32_t& iID, Build<vecType, elemType>& vBricks, const float& fShrinkPercentage) {
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
					//if (_bSpaced) {
					//	dnaToAminoacidSpaced(sDna, iMaxKTimes3, j, &sTempFrame);
					//}
					//else {
						dnaToAminoacid(sDna, iMaxKTimes3, j, &sTempFrame);
					//}
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
							//if (_bSpaced) {
							//	dnaToAminoacidSpaced(sDna, k + iMaxKTimes3 + 3 * (j - 1), sTempAA);
							//}
							//else {
								dnaToAminoacid(sDna, k + iMaxKTimes3 + 3 * (j - 1), sTempAA);
							//}
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
						//if (_bSpaced) {
						//	dnaToAminoacidSpaced(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
						//}
						//else {
							dnaToAminoacid(sDna, j + iMaxKTimes3 + 3 * (iMaxRange / 3 - 1), sTempAA);
						//}
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

			// Write resulting kMers in stxxl vector but exclude those which are not useful
			double dStepSize = (fShrinkPercentage > 0.f) ? 100. / fShrinkPercentage : 0.;
			double dNextThrowOut = dStepSize;
			uint64_t iCounterOfThrowOut = 1;

			if (fShrinkPercentage > 0.f) {
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

							if (!(vBricks.addToInt(element))) {
								vBricks.IntToExtPart();
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

						if (!(vBricks.addToInt(element))) {
							vBricks.IntToExtPart();
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// convert lines from fasta/fastq to kMers and save them in a (stxxl-)vector, use only one frame
		inline void dnaTokMers_oneFrame(const string& sDna, const uint32_t& iID, Build<vecType, elemType>& vBricks, const float& fShrinkPercentage) {
			vector<tuple<intType, uint32_t>> vResultingkMers;
			const int32_t& iMaxKTimes3 = 3 * _iHighestK;

			const int32_t& iMaxRange = int32_t(sDna.length()) - iMaxKTimes3 + 1;
			if (iMaxRange > 0) {
				// go through the dna, convert it framewise to an aminoacid kMer and then to its coded representation
				const int32_t& iNumFrames = 1;

				vResultingkMers.resize(iMaxRange);
				int32_t ikMerCounter = 0;

				// Frameshifting, so that only one amino acid has to be computed
				// Compute initial frames
				intType sAAFrames[1] = { 0 };
				uint32_t aDeletekMerCounter[1] = { 0 };
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

					for (int32_t j = 1; 3 * j < iMaxRange; ++j) {
						for (int32_t k = 0; k < 1; ++k) {
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
						ikMerCounter++;
					}
				}
			}

			// Write resulting kMers in stxxl vector but exclude those which are not useful
			double dStepSize = (fShrinkPercentage > 0.f) ? 100. / fShrinkPercentage : 0.;
			double dNextThrowOut = dStepSize;
			uint64_t iCounterOfThrowOut = 1;

			if (fShrinkPercentage > 0.f) {
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

							if (!(vBricks.addToInt(element))) {
								vBricks.IntToExtPart();
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

						if (!(vBricks.addToInt(element))) {
							vBricks.IntToExtPart();
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// convert protein sequences to k-mers in the building step
		inline void proteinTokMers(const string& sAASequence, const uint32_t& iIdx, Build<vecType, elemType>& vBricks, const float& fShrinkPercentage) {
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

			if (fShrinkPercentage > 0.f) {
				for (auto& element : vResultingkMers) {
					if (iCounterOfThrowOut != static_cast<uint64_t>(dNextThrowOut)) {
						if (_bUnfunny) {
							get<0>(element) = aminoAcidsToAminoAcid(get<0>(element));
						}

						if (!(vBricks.addToInt(element))) {
							vBricks.IntToExtPart();
						}
					}
					else {
						dNextThrowOut += dStepSize;
					}
					++iCounterOfThrowOut;
				}
			}
			else {
				for (auto& element : vResultingkMers) {
					if (_bUnfunny) {
						get<0>(element) = aminoAcidsToAminoAcid(get<0>(element));
					}

					if (!(vBricks.addToInt(element))) {
						vBricks.IntToExtPart();
					}
				}
			}
		}

	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fasta and create a kMer-Vec, used in BuildAll(...)
		template<typename T>
		inline void readFasta(T& input, const unordered_map<string, uint32_t>& mAccToID, Build<vecType, elemType>& vBricks, const uint64_t& iFileLength, size_t& overallCharsRead, const size_t& overallFilesSize, const float& fShrinkPercentage) {
			try {
				if (!input.good()) {
					throw runtime_error("No input found!");
				}

				Utilities::FileReader<T> inputReader; 
				inputReader.setFile(&input);
				pair<std::string, bool> resultChunkPair;
				uint64_t iNumOfChars = 0;

				uint64_t iNumOfCharsRead = 0, iCurrentPercentage = 0, iCurrentOverallPercentage = 0;

				uint32_t iIdx = 0, nextIdx = 0;

				bool bAccNumberFound = false;

				string sFalsekMerMarker = "";
				if (this->_bProtein) {
					for (int32_t i = 0; i < (_iHighestK - _iLowestK); ++i) {
						sFalsekMerMarker += "^";
					}
				}
				else {
					for (int32_t i = 0; i < (_iHighestK - _iLowestK) * 3; ++i) {
						sFalsekMerMarker += "X";
					}
				}

				while (!inputReader.eof()) {
					resultChunkPair.first = "";
					resultChunkPair.second = false;
					inputReader.getChunk(resultChunkPair, iNumOfChars);
					iNumOfCharsRead += iNumOfChars;

					if (resultChunkPair.first == "") {
						continue;
					}
					if (resultChunkPair.first.front() == '>') {
						resultChunkPair.first.erase(resultChunkPair.first.begin());
						const auto& sNumbers = Utilities::split(Utilities::split(resultChunkPair.first, ' ').at(0), '|');
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
							isInMap = mAccToID.find(resultChunkPair.first);
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

				// Either the file still contains dna or there is an Overhang left from a previous loop
				while (!inputReader.eof() || sOverhang != "") {

					if (bNextIdx) {
						bAddFalseMarker = false;
						bAddFalseRCMarker = true;
					}

					bool bCreateOverhang = true;

					for (int32_t i = 0; i < iConcurrentLines; ++i) {
						resultChunkPair.first = "";
						resultChunkPair.second = false;
						inputReader.getChunk(resultChunkPair, iNumOfChars);

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

						string sTempString(resultChunkPair.first);
						if (sTempString == "") {
							if (!inputReader.eof()) {
								break;
							}
							else {
								bCreateOverhang = false;
								bAddFalseMarker = true;
								break;
							}
							--i;
							continue;
						}
						if (!inputReader.eof()) {

							// new entry in fasta
							if (sTempString.front() == '>') {

								sTempString.erase(sTempString.begin());
								while (resultChunkPair.second == false) {
									// read the whole line
									resultChunkPair.first = "";
									resultChunkPair.second = false;
									inputReader.getChunk(resultChunkPair, iNumOfChars);
									iNumOfCharsRead += iNumOfChars;
									overallCharsRead += iNumOfChars;
									sTempString.append(resultChunkPair.first);
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

					if (this->_bProtein) {
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
							if (this->_bProtein) {
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

					if (this->_bProtein) {
						proteinTokMers(sDNA, iIdx, vBricks, fShrinkPercentage);
					}
					else {
						if (_bOnlyOneFrame) {
							dnaTokMers_oneFrame(sDNA, iIdx, vBricks, fShrinkPercentage);
						}
						else {
							dnaTokMers(sDNA, iIdx, vBricks, fShrinkPercentage);
							if (_bSixFrames) {
								dnaTokMers(sRCDNA, iIdx, vBricks, fShrinkPercentage);
							}
						}
					}

					if (bNextIdx) {
						iIdx = nextIdx;
						bNextIdx = false;
					}

					if (bCreateOverhang) {
						if (this->_bProtein) {
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

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Go through file and map accs to taxids
		inline void accToTaxST(unordered_map<string, tuple<uint32_t, string>>& mAccTo_ID_DNA_NumOfkMers, const string& content) {
			debugBarrier
			ifstream contentFile(content);
			size_t iNumOfEntries = mAccTo_ID_DNA_NumOfkMers.size(), iCounterOfFoundElements = 0;
			string item, sDummy;
			while (getline(contentFile, sDummy) && iCounterOfFoundElements < iNumOfEntries) {
				stringstream ss(sDummy);
				getline(ss, item, '\t');
				getline(ss, item, '\t');
				const string tID = move(item);
				getline(ss, item, '\t');
				while (getline(ss, item, ';')) {
					auto res = mAccTo_ID_DNA_NumOfkMers.find(item);
					if (res != mAccTo_ID_DNA_NumOfkMers.end()) {
						get<0>(res->second) = stoul(tID);
						iCounterOfFoundElements++;
					}
				}
			}
			debugBarrier
		}

		inline void findAccInLine(const string& sLine, unordered_map<string, tuple<uint32_t, string>>& mAccTo_ID_DNA_NumOfkMers) {
			string item;
			stringstream ss(sLine);
			getline(ss, item, '\t');
			getline(ss, item, '\t');
			const string tID = move(item);
			getline(ss, item, '\t');
			while (getline(ss, item, ';')) {
				auto res = mAccTo_ID_DNA_NumOfkMers.find(item);
				if (res != mAccTo_ID_DNA_NumOfkMers.end()) {
					get<0>(res->second) = stoul(tID);
				}
			}
		}

		inline void accToTaxMT(unordered_map<string, tuple<uint32_t, string>>& mAccTo_ID_DNA_NumOfkMers, const string& content, WorkerQueueWithIDs& Q) {
			debugBarrier
			ifstream contentFile(content);
			Utilities::FileReader<ifstream> fileReader;
			fileReader.setFile(&contentFile);
			pair<string, bool> pChunk("",false);
			uint64_t iDummyNumOfChars = 0;
			string sLine = "", item;
			int32_t iCounterOfConcurrentLinesBeingProcessed = 0;
			while (!fileReader.eof()) {
				if (iCounterOfConcurrentLinesBeingProcessed == 2*_iNumOfThreads) {
					Q.waitUntilFinished();
					iCounterOfConcurrentLinesBeingProcessed = 0;
				}
				fileReader.getChunk(pChunk, iDummyNumOfChars);
				while (!pChunk.second) {
					sLine += pChunk.first; // maybe this is not a good idea if the line is very large...
					pChunk.first = "";
					fileReader.getChunk(pChunk, iDummyNumOfChars);
				}
				sLine += pChunk.first;
				Q.pushTask([&mAccTo_ID_DNA_NumOfkMers, this, sLine](const int32_t&) { findAccInLine(sLine, mAccTo_ID_DNA_NumOfkMers); });
				sLine = "";
				pChunk.first = "";
				iCounterOfConcurrentLinesBeingProcessed++;
			}
			Q.waitUntilFinished();
			debugBarrier
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Translate and convert for readFastaAlternativeMode
		inline void translateAndConvert(const unordered_map<string, tuple<uint32_t, string>>& mAccTo_ID_DNA_NumOfkMers, const string& sAccThatHasBeenPushedAlready, Build<vecType, elemType>& vBricks, const float& fShrinkPercentage) {
			debugBarrier
			string sFalsekMerMarker = "";
			if (this->_bProtein) {
				for (int32_t i = 0; i < (_iHighestK - _iLowestK); ++i) {
					sFalsekMerMarker += "^";
				}
			}
			else {
				for (int32_t i = 0; i < (_iHighestK - _iLowestK) * 3; ++i) {
					sFalsekMerMarker += "X";
				}
			}
			
			for (auto mapIt = mAccTo_ID_DNA_NumOfkMers.cbegin(); mapIt != mAccTo_ID_DNA_NumOfkMers.cend(); mapIt++) {
				if (get<0>(mapIt->second)) { // only translate entries with valid taxids (0 is a signal for non-valid) 
					if (this->_bProtein) {
						proteinTokMers(get<1>(mapIt->second), get<0>(mapIt->second), vBricks, fShrinkPercentage);
					}
					else {
						if (_bOnlyOneFrame) {
							dnaTokMers_oneFrame(get<1>(mapIt->second), get<0>(mapIt->second), vBricks, fShrinkPercentage);
						}
						else {
							dnaTokMers(get<1>(mapIt->second), get<0>(mapIt->second), vBricks, fShrinkPercentage);
							if (_bSixFrames) {
								string sRCDNA = (_bSixFrames) ? reverseComplement(get<1>(mapIt->second)) : "";
								if (mapIt->first != sAccThatHasBeenPushedAlready) {
									sRCDNA += sFalsekMerMarker;
								}
								dnaTokMers(sRCDNA, get<0>(mapIt->second), vBricks, fShrinkPercentage);
							}
						}
					}
				}
			}
			debugBarrier
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read a fasta and create a kMer-Vec, used in BuildAll(...)
		template<typename T>
		inline void readFastaAlternativeMode(Utilities::FileReader<T>& input, const string& contentFile, Build<vecType, elemType>& vBricks, const uint64_t& iFileLength, size_t& overallCharsRead, const size_t& overallFilesSize, const float& fShrinkPercentage, const int64_t& iMemoryAvail) {
			try {
				debugBarrier

				WorkerQueueWithIDs Q(_iNumOfThreads);
				pair<std::string, bool> resultChunkPair;

				const float fShrinkFactor = 1.0f - fShrinkPercentage / 100.0f;

				unordered_map<string, tuple<uint32_t, string>> mAccTo_ID_DNA_NumOfkMers;
				unordered_map<string, tuple<uint32_t, string>>::iterator currElement;
				string sAccThatHasBeenPushedAlready = "";
				int64_t iMemoryUsed = sizeof(unordered_map<string, tuple<uint32_t, string>>);

				bool bNewAcc = true;

				string sOverhang = "", sDNA = "", currAcc = "", newAcc = "";

				uint64_t iNumOfCharsRead = 0, iNumOfChars = 0, iCurrentPercentage = 0, iCurrentOverallPercentage = 0, iRequiredSizeOfInternal = 0;

				// add those to correctly parse smaller k-mers from larger ones
				string sFalsekMerMarker = "";
				if (this->_bProtein) {
					for (int32_t i = 0; i < (_iHighestK - _iLowestK); ++i) {
						sFalsekMerMarker += "^";
					}
				}
				else {
					for (int32_t i = 0; i < (_iHighestK - _iLowestK) * 3; ++i) {
						sFalsekMerMarker += "X";
					}
				}
				debugBarrier
				// Read data from file, skip irrelevant lines
				while (!input.eof()) {

					for (int32_t i = 0; i < 1; ++i) {
						resultChunkPair.first = "";
						resultChunkPair.second = false;
						input.getChunk(resultChunkPair, iNumOfChars);
						debugBarrier
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

						string sTempString(resultChunkPair.first);
						if (sTempString == "") {
							if (input.eof()) {
								break;
							}
							--i;
							continue;
						}
						if (!input.eof()) {
							debugBarrier
							// new entry in fasta
							if (sTempString.front() == '>') {

								sTempString.erase(sTempString.begin());
								while (resultChunkPair.second == false) {
									// read the whole line
									resultChunkPair.first = "";
									resultChunkPair.second = false;
									input.getChunk(resultChunkPair, iNumOfChars);
									iNumOfCharsRead += iNumOfChars;
									overallCharsRead += iNumOfChars;
									sTempString.append(resultChunkPair.first);
								}
								const auto& sNumbers = Utilities::split(Utilities::split(sTempString, ' ').at(0), '|');
								string acc = "";
								for (const auto& entry : sNumbers) {
									if (entry.find('.') != string::npos) {
										acc = entry;
										break;
									}
								}
								
								newAcc = acc;
								bNewAcc = true;

								break;
							}
							else {
								sDNA = sTempString;
							}
						}
						else {
							// EOF
							if (currAcc == "") {
								throw runtime_error("No > found in input.");
							}
							debugBarrier
							sDNA = sTempString;
							break;
						}
					}
					debugBarrier
					// check for irregular letters
					for (auto& chara : sDNA) {
						if (chara == '\t' || chara == ' ') {
							throw runtime_error("Spaces or tabs inside reference, please check your input. Error occured in content entry: " + currAcc);
						}
						else {
							if (this->_bProtein) {
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
					debugBarrier

					// new acc has been found, save the old one
					if (bNewAcc) {
						if (currAcc != "") {
							get<1>(currElement->second) += sFalsekMerMarker;
							iRequiredSizeOfInternal += _iHighestK - 1;

							iMemoryUsed += sizeof(char) * sFalsekMerMarker.length() + sizeof(elemType) * (_iHighestK - 1);
						}
						debugBarrier
						currAcc = newAcc;
						currElement = mAccTo_ID_DNA_NumOfkMers.insert(make_pair(newAcc, make_tuple(0, ""))).first;
						iMemoryUsed += sizeof(char) * newAcc.length() + static_cast<int64_t>(1.5*sizeof(pair< string, tuple<uint32_t, string>>));

						bNewAcc = false;
					}
					else {
						if (input.eof()) {
							get<1>(currElement->second) += sFalsekMerMarker;
							iRequiredSizeOfInternal += _iHighestK - 1;
							debugBarrier
							// save to external
							// the entry in sAccThatHasBeenPushedAlready shall receive no false kMers for RC (since it had them already), check if it is still in map
							accToTaxMT(mAccTo_ID_DNA_NumOfkMers, contentFile, Q);
							debugBarrier
							vBricks.setInternalSize(iRequiredSizeOfInternal);
							debugBarrier
							translateAndConvert(mAccTo_ID_DNA_NumOfkMers, sAccThatHasBeenPushedAlready, vBricks, fShrinkPercentage);
							debugBarrier
							vBricks.IntToExtPart();
							debugBarrier
						}
						else {

							uint64_t iNumOfkMers = 0;
							if (get<1>(currElement->second) == "") {
								iNumOfkMers = (_bSixFrames) ? 2 * calculatekMerCount(sDNA.length()) : calculatekMerCount(sDNA.length());
							}
							else {
								iNumOfkMers = (_bSixFrames) ? 2 * (calculatekMerCount(sDNA.length()) + _iHighestK - 1) : (calculatekMerCount(sDNA.length()) + _iHighestK - 1);
							}
							debugBarrier
							int64_t iAdditionalMemoryUsage = sizeof(char) * sDNA.length() + sizeof(elemType) * static_cast<size_t>(fShrinkFactor * iNumOfkMers); // depending on RC or not
							if (iMemoryUsed + iAdditionalMemoryUsage > iMemoryAvail) {
								// save to external
								// the entry in sAccThatHasBeenPushedAlready shall receive no false kMers for RC (since it had them already), check if it is still in map
								debugBarrier
								accToTaxMT(mAccTo_ID_DNA_NumOfkMers, contentFile, Q);
								debugBarrier
								vBricks.setInternalSize(iRequiredSizeOfInternal);
								debugBarrier
								translateAndConvert(mAccTo_ID_DNA_NumOfkMers, sAccThatHasBeenPushedAlready, vBricks, fShrinkPercentage);
								debugBarrier
								vBricks.IntToExtPart();
								debugBarrier
								// create overhang to new dna from old dna which will now be converted
								string sTempDNA = get<1>(currElement->second);
								sOverhang = "";
								if (sTempDNA.length()) {
									if (this->_bProtein) {
										if (sTempDNA.length() + 1 > static_cast<size_t>(_iHighestK)) {
											sOverhang = sTempDNA.substr(sTempDNA.length() + 1 - _iHighestK);
										}
										else {
											sOverhang = sTempDNA;
										}
									}
									else {
										if (sTempDNA.length() + 1 > static_cast<size_t>(_iHighestK * 3)) {
											sOverhang = sTempDNA.substr(sTempDNA.length() + 1 - _iHighestK * 3);
										}
										else {
											sOverhang = sTempDNA;
										}
									}
								}

								sAccThatHasBeenPushedAlready = currElement->first;
								auto tempElement = *currElement;
								mAccTo_ID_DNA_NumOfkMers.clear();
								currElement = mAccTo_ID_DNA_NumOfkMers.insert(tempElement).first;
								get<1>(currElement->second) = sOverhang;
								const auto& kmerCount = (_bSixFrames) ? 2 * calculatekMerCount(sOverhang.length()) : calculatekMerCount(sOverhang.length());
								iRequiredSizeOfInternal = static_cast<size_t>(fShrinkFactor * kmerCount);
								debugBarrier
								iMemoryUsed = sizeof(unordered_map<string, tuple<uint32_t, string>>) + sizeof(char) * sOverhang.length() + sizeof(elemType) * static_cast<size_t>(fShrinkFactor * kmerCount);
							}

							get<1>(currElement->second) += sDNA;
							iRequiredSizeOfInternal += static_cast<size_t>(fShrinkFactor * iNumOfkMers);

							iMemoryUsed += iAdditionalMemoryUsage;
							debugBarrier
						}
					}
					debugBarrier
					sDNA = "";
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}



		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Build the index from fasta files
		tuple< unordered_map<uint32_t, uint32_t>, uint64_t >  BuildAll(const InputParameters& cParams, const string& fOutFile, const bool& bContinue = false, const bool& bUpdate = false) {
			try {

				int64_t iMem = static_cast<uint64_t>(cParams.iMemorySizeAvail * 0.9 - GIGABYTEASBYTES);

				// test if files exists
				if (!ifstream(cParams.contentFileIn)) {
					throw runtime_error("Content file not found.");
				}

				/*unique_ptr<uint64_t[]> vDivisionArray(new uint64_t[_iNumOfThreads - 1]);
				const uint64_t& iStep = 0xF7BDEF7BDEF7BDE / _iNumOfThreads;
				for (int i = 0; i < _iNumOfThreads - 1; ++i) {
					vDivisionArray[i] = (i + 1) * iStep;
				}*/

				debugBarrier
				Utilities::checkIfFileCanBeCreated(fOutFile);
				
				//stxxlFile* stxxlOutFile = new stxxlFile(fOutFile, stxxl::file::RDWR);
				//contentVecType_32p* vOutVec = new contentVecType_32p(stxxlOutFile, 0);

				if (_bVerbose) {
					cout << "OUT: Reading content file... " << endl;
				}

				int64_t iMemGiven = iMem;
				uint32_t iIdxCounter = 1;
				unordered_map<uint32_t, uint32_t> mIDsAsIdx; mIDsAsIdx[0] = 0;
				//unordered_map<uint32_t, string> mIdxToName; mIdxToName[0] = "non_unique";
				unordered_map<string, uint32_t> mAccToID;
				ifstream content(cParams.contentFileIn);
				string sDummy = "";

				bool bTaxIdsAsStrings = false, bAlternativeMode = false;
				
				while (getline(content, sDummy)) {
					if (iMem < iMemGiven / 2 && !bAlternativeMode) {
						cerr << "WARNING: Your content file is quite large and keeping the hash tables generated from it in memory consumes more than half of the memory given. \
kASA will therefore use an alternative but slower way of generating the index. You can Ctrl+C now if you don't want this and instead reduce the number of entries in the content file or abstract to a higher taxonomic rank and generate it anew. \
Sorry!" << endl;
						bAlternativeMode = true;
					}

					if (sDummy != "") {
						const auto& line = Utilities::split(sDummy, '\t');
						if (line.size() >= 5 && !bTaxIdsAsStrings) {
							bTaxIdsAsStrings = true;
						}
						if (line.size() >= 4) {
							if (bTaxIdsAsStrings) {
								mIDsAsIdx[stoul(line[4])] = iIdxCounter;
								iMem -= sizeof(pair<uint32_t, uint32_t>);
								if (!bAlternativeMode && !bContinue) {
									const auto& vAccessionNumbers = Utilities::split(line[3], ';');
									for (const auto& acc : vAccessionNumbers) {
										mAccToID.insert(make_pair(acc, stoul(line[4])));
										iMem -= static_cast<int64_t>(1.5*sizeof(pair< string, uint32_t>)) + sizeof(char) * acc.length();
									}
								}
							}
							else {
								mIDsAsIdx[stoul(line[1])] = iIdxCounter;
								iMem -= sizeof(pair<uint32_t, uint32_t>);
								if (!bAlternativeMode && !bContinue) {
									const auto& vAccessionNumbers = Utilities::split(line[3], ';');
									for (const auto& acc : vAccessionNumbers) {
										mAccToID.insert(make_pair(acc, stoul(line[1])));
										iMem -= static_cast<int64_t>(1.5*sizeof(pair< string, uint32_t>)) + sizeof(char) * acc.length();
									}
								}
							}
							++iIdxCounter;
						}
						else {
							throw runtime_error("Content file contains less than 4 columns, it may be damaged... The faulty line was: " + sDummy + "\n");
						}
					}
				}
				content.close();
				content.clear();

				if (bAlternativeMode) {
					iMemGiven -= Utilities::calculateSizeInByteOfUnorderedMap(mIDsAsIdx);
					mAccToID.clear();
				}
				
				debugBarrier

				if (_bVerbose && !bAlternativeMode) {
					cout << "OUT: Creating internal k-mer container... " << endl;
				}

				Build<vecType, elemType> brick;

				if (bAlternativeMode || bContinue) {
					brick.setMemberVariables(_sTemporaryPath, _iNumOfCall, _iNumOfThreads, 0);
				}
				else {
					brick.setMemberVariables(_sTemporaryPath, _iNumOfCall, _iNumOfThreads, static_cast<size_t>(iMem) / sizeof(elemType));
				}

				debugBarrier

				size_t overallCharsRead = 0;

				auto filesAndSize = Utilities::gatherFilesFromPath(cParams.sInput);
				if (!bContinue) {
					for (auto& fileName : filesAndSize.first) {
						if (_bVerbose) {
							cout << "OUT: Current file: " << fileName.first << endl;
						}

						// check if gzipped
						bool isGzipped = fileName.second;

						debugBarrier

							if (isGzipped) {
								if (_bVerbose) {
									cout << "OUT: File is gzipped, no progress output can be shown." << endl;
								}
								unique_ptr<igzstream> fastaFile_gz(new igzstream(fileName.first.c_str()));
								
								// Determine if protein or DNA sequence
								if (!this->_bProtein) {
									const string& sFirstSequence = Utilities::getFirstSequenceOfFile(*fastaFile_gz);
									this->detectAlphabet(sFirstSequence);
									fastaFile_gz.reset(new igzstream(fileName.first.c_str()));
								}
								
								if (bAlternativeMode) {
									Utilities::FileReader<igzstream> fastaFileReader;
									fastaFileReader.setFile(fastaFile_gz.get());
									readFastaAlternativeMode(fastaFileReader, cParams.contentFileIn, brick, 0, overallCharsRead, filesAndSize.second, cParams.fPercentageOfThrowAway, iMemGiven);
								}
								else {
									readFasta(*fastaFile_gz, mAccToID, brick, 0, overallCharsRead, filesAndSize.second, cParams.fPercentageOfThrowAway);
								}
							}
							else {
								ifstream fastaFile(fileName.first);

								// Determine if protein or DNA sequence
								if (!this->_bProtein) {
									const string& sFirstSequence = Utilities::getFirstSequenceOfFile(fastaFile);
									this->detectAlphabet(sFirstSequence);
								}

								const uint64_t& iFileLength = Utilities::getSizeOfFile(fileName.first);

								if (bAlternativeMode) {
									Utilities::FileReader<ifstream> fastaFileReader;
									fastaFileReader.setFile(&fastaFile);
									readFastaAlternativeMode(fastaFileReader, cParams.contentFileIn, brick, iFileLength, overallCharsRead, filesAndSize.second, cParams.fPercentageOfThrowAway, iMemGiven);
								}
								else {
									readFasta(fastaFile, mAccToID, brick, iFileLength, overallCharsRead, filesAndSize.second, cParams.fPercentageOfThrowAway);
								}
							}
						debugBarrier
					}
				}

				if (_bVerbose) {
					cout << "OUT: Merging temporary index files... " << endl;
				}
				
				debugBarrier
				// Finalize
				if (!bAlternativeMode && !bContinue) {
					brick.IntToExtPart();
				}

				if (bContinue) {
					auto vecOfSizes = brick.getVectorSizesVec();
					auto filesInTemp = Utilities::gatherFilesAndSizesFromPath(_sTemporaryPath);
					sort(filesInTemp.begin(), filesInTemp.end(), [](const pair<string, size_t>& a, const pair<string, size_t>& b) { return a.first < b.first; });
					brick.setNumOfContainers(filesInTemp.size());
					for (const auto& entry : filesInTemp) {
						vecOfSizes->push_back(entry.second / sizeof(elemType));
					}
				}
				const uint64_t& iSizeOfFinalIndex = (GlobalInputParameters.bIGotSpace) ? brick.mergeTemporariesWithEnoughSpace(fOutFile) : brick.mergeTemporaries(fOutFile);

				debugBarrier

				if (iSizeOfFinalIndex == 0) {
					throw runtime_error("Index is empty, are all input files okay?");
				}

				if (!bUpdate) {

					if (_bVerbose) {
						cout << "OUT: Creating trie... " << endl;
					}

					debugBarrier

					// Create Trie and frequencies out of final file
					unique_ptr<stxxlFile> libFile(new stxxlFile(fOutFile, stxxlFile::RDONLY));
					unique_ptr<const vecType> libVec(new const vecType(libFile.get(), iSizeOfFinalIndex));
					Trie<intType> T(static_cast<int8_t>(((_iMaxK > 12) ? HIGHESTPOSSIBLEK : 12)), static_cast<int8_t>(_iMinK), 6);
					T.SaveToStxxlVec(libVec.get(), fOutFile);
					debugBarrier

					// If taxaOnly index is desired, create it here:
					if (_bUnfunny) {
						Utilities::checkIfFileCanBeCreated(fOutFile + "_taxOnly");
						stxxlFile* taxaOnlyFile = new stxxlFile(fOutFile + "_taxOnly", stxxlFile::RDWR);
						taxaOnly* taxaOnlyFileVec = new taxaOnly(taxaOnlyFile, iSizeOfFinalIndex);

						auto tOIt = taxaOnlyFileVec->begin();
						for (const auto& elem : *libVec) {
							*tOIt = static_cast<uint16_t>(mIDsAsIdx.find(elem.second)->second);
							++tOIt;
						}

						taxaOnlyFileVec->export_files("_");
						delete taxaOnlyFileVec;
						delete taxaOnlyFile;

						remove(fOutFile.c_str());
						Utilities::copyFile(fOutFile + "_taxOnly", fOutFile);
					}

					debugBarrier

					if (_bVerbose) {
						cout << "OUT: Creating frequency file... " << endl;
					}

					libVec.reset();
					libFile.reset();
					// create frequency file
					this->GetFrequencyK<vecType>(cParams, cParams.contentFileIn, fOutFile, mIDsAsIdx);
				}

				debugBarrier
				return make_tuple(mIDsAsIdx, iSizeOfFinalIndex);
			}
			catch (invalid_argument&) {
				cerr << "ERROR: content file doesn't have the right format. If you'd like to use non-numeric taxids, please apply the --taxidasstr flag." << endl;
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
			return make_tuple(unordered_map<uint32_t, uint32_t>(), 0);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Merge two existing indices into a new one
		void MergeTwoIndices(const InputParameters& cParams, const pair<unordered_map<uint32_t, uint32_t>, unordered_map<uint32_t, uint32_t>>& mapsForDummys = make_pair(unordered_map<uint32_t, uint32_t>(), unordered_map<uint32_t, uint32_t>())) {
			try {

				// test if files exists
				if (!ifstream(cParams.contentFileAfterUpdate)) {
					throw runtime_error("Content file not found.");
				}

				Build<vecType, elemType> brick;

				// indices
				ifstream fLibInfo(cParams.firstOldIndex + "_info.txt");
				uint64_t iSizeOfFirstLib = 0, iSizeOfSecondLib = 0;
				fLibInfo >> iSizeOfFirstLib;
				fLibInfo.close();
				unique_ptr<stxxlFile> stxxlVecI1(new stxxlFile(cParams.firstOldIndex, stxxl::file::RDONLY));
				unique_ptr<vecType> vec1(new vecType(stxxlVecI1.get(), iSizeOfFirstLib));
				typename vecType::bufreader_type vec1Buff(*vec1);

				fLibInfo.open(cParams.secondOldIndex + "_info.txt");
				fLibInfo >> iSizeOfSecondLib;
				fLibInfo.close();
				unique_ptr<stxxlFile> stxxlVecI2(new stxxlFile(cParams.secondOldIndex, stxxl::file::RDONLY));
				unique_ptr<vecType> vec2(new vecType(stxxlVecI2.get(), iSizeOfSecondLib));
				

				Utilities::checkIfFileCanBeCreated(cParams.sDBPathOut);
				unique_ptr<stxxlFile> stxxlVecOut(new stxxlFile(cParams.sDBPathOut, stxxl::file::RDWR));
				unique_ptr<vecType> vecOut(new vecType(stxxlVecOut.get(), iSizeOfFirstLib + iSizeOfSecondLib));
				typename vecType::bufwriter_type vecOutBuff(*vecOut);

				if (_bVerbose) {
					cout << "OUT: Merging files... " << endl;
				}

				//merge
				const auto& vOutSize = brick.merge(vecOutBuff, vec1Buff, iSizeOfFirstLib, vec2->cbegin(), vec2->cend(), mapsForDummys);
				vecOut->resize(vOutSize, true);
				vecOutBuff.finish();

				// save additional stuff
				vecOut->export_files("_");
				
				if (_bVerbose) {
					cout << "OUT: Creating trie... " << endl;
				}

				// trie
				Trie<intType> T(static_cast<int8_t>(((is_same<vecType, contentVecType_128>::value) ? HIGHESTPOSSIBLEK : 12)), static_cast<int8_t>(_iMinK), 6);
				T.SaveToStxxlVec(vecOut.get(), cParams.sDBPathOut);

				vecOut.reset();
				stxxlVecOut.reset();

				if (_bVerbose) {
					cout << "OUT: Creating frequency file... " << endl;
				}

				this->GetFrequencyK<vecType>(cParams, cParams.contentFileAfterUpdate, cParams.sDBPathOut);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
	};
}