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
#include "Read.hpp"
#include "Trie.hpp"
#include "../utils/WorkerThread.hpp"
#include "../utils/ParallelQuicksort.hpp"

#include "../utils/dToStr.h"
#include "../utils/iToStr.hpp"

namespace kASA {
	class Compare : public Read {

		const bool _bTranslated = false;
		bool bUnfunny = false;
		mutex m_exceptionLock;
		exception_ptr someThingWentWrong;

	public:
		Compare(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const int32_t& iNumOfBeasts, const bool& bVerbose = false, const bool& bProtein = false, const string& stxxl_mode = "", const bool& bUnfunny = false, const bool& bSixFrames = false) : Read(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, bProtein, stxxl_mode, bUnfunny, bSixFrames), _bTranslated(bProtein), bUnfunny(bUnfunny), iNumOfBeasts(iNumOfBeasts) {}

		// for output
		int32_t iNumOfBeasts = 3;

		enum OutputFormat
		{
			Kraken,
			Json,
			JsonL,
			tsv
		} format;

	private:

		string outputFormatFileEnding() {
			if (format == OutputFormat::Kraken) {
				return ".ktsv";
			}
			if (format == OutputFormat::Json) {
				return ".json";
			}
			if (format == OutputFormat::JsonL) {
				return ".jsonl";
			}
			if (format == OutputFormat::tsv) {
				return ".tsv";
			}
			return ".rtt";
		}

		///////////////////////////////////////////////////////
		const float arrWeightingFactors[12] = { 1, 121.f / 144.f , 100.f / 144.f , 81.f / 144.f , 64.f / 144.f , 49.f / 144.f , 36.f / 144.f , 25.f / 144.f , 16.f / 144.f , 9.f / 144.f , 4.f / 144.f , 1.f / 144.f };
		//const float arrWeightingFactors[12] = { 1, 1331.f / 1728.f , 1000.f / 1728.f , 729.f / 1728.f , 512.f / 1728.f , 343.f / 1728.f , 216.f / 1728.f , 125.f / 1728.f , 64.f / 1728.f , 27.f / 1728.f , 8.f / 1728.f , 1.f / 1728.f };
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void markTaxIDs(const uint64_t& codedTaxIDs, Utilities::sBitArray& vMemoryOfTaxIDs_k, const unordered_map<uint32_t, uint32_t>& mTaxToIdx, unique_ptr<const contentVecType_32p>&) {
			try {
				vMemoryOfTaxIDs_k.set(Utilities::checkIfInMap(mTaxToIdx, static_cast<uint32_t>(codedTaxIDs))->second);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
		inline void markTaxIDs(const uint64_t& codedTaxIDs, Utilities::sBitArray& vMemoryOfTaxIDs_k, const unordered_map<uint32_t, uint32_t>& mTaxToIdx, const vector<packedBigPair>*) {
			try {
				vMemoryOfTaxIDs_k.set(Utilities::checkIfInMap(mTaxToIdx, static_cast<uint32_t>(codedTaxIDs))->second);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
		inline void markTaxIDs(const uint64_t& codedTaxIDs, Utilities::sBitArray& vMemoryOfTaxIDs_k, const unordered_map<uint32_t, uint32_t>&, const vector<packedPair>*) {
			try {
				vMemoryOfTaxIDs_k.set(codedTaxIDs);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
		inline void markTaxIDs(const uint64_t& codedTaxIDs, Utilities::sBitArray& vMemoryOfTaxIDs_k, const unordered_map<uint32_t, uint32_t>&, unique_ptr<const index_t_p>&) {
			try {
				/*while (codedTaxIDs != 0ul) {
					const uint8_t& numOfBits = codedTaxIDs & 31;
					codedTaxIDs >>= 5;
					vMemoryOfTaxIDs_k.set(static_cast<uint16_t>(codedTaxIDs & Utilities::_bitMasks[numOfBits]));
					codedTaxIDs >>= numOfBits;
				}*/
				vMemoryOfTaxIDs_k.set(codedTaxIDs);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}


		///////////////////////////////////////////////////////
		inline unique_ptr<const contentVecType_32p>& getVec(const unique_ptr<unique_ptr<const contentVecType_32p>[]>& vec, const int& iThreadID) {
			return vec[iThreadID];
		}
		//inline const unique_ptr<const contentVecType_32p[]>& getVec(const unique_ptr<const contentVecType_32p[]>& vec, const int&) {
		//	return vec;
		//}
		inline unique_ptr<const index_t_p>& getVec(const unique_ptr<unique_ptr<const index_t_p>[]>& vec, const int& iThreadID) {
			return vec[iThreadID];
		}
		inline const vector<packedBigPair>* getVec(const vector<packedBigPair>* vec, const int&) {
			return vec;
		}
		inline const vector<packedPair>* getVec(const vector<packedPair>* vec, const int&) {
			return vec;
		}

		inline void addToMatchedReadID(vector<uint64_t>& vReadIDs, uint64_t& position, const uint64_t& value) {
			if (vReadIDs.size() <= position) {
				vReadIDs.resize(position + 100);
			}

			vReadIDs[position] = value;
			position++;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Compare as many as #Number-of-processors vectors with an index lying on a HDD/SSD and note all similarities for any k. 
		// To minimize hard disk access, the order is as follows: Get kMer from RAM Vec -> Search in Prefix-Trie -> Get range of possible hit -> binary search in that range -> note if hit
		template <typename vecType>
		inline void compareWithDatabase(const int32_t& iThreadID, const vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& vIn, const uint64_t& vInStart, const uint64_t& vInEnd, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, unique_ptr<float[]>& vReadIDtoGenID, const uint32_t& iSpecIDRange, const unordered_map<uint32_t, uint32_t>& mTaxToIdx) {//, const unordered_map<readIDType, uint64_t>& mReadIDToArrayIdx) {

			try {

				//const uint64_t iReadIDTODEBUG = 0;//424734

				unique_ptr<vector<uint64_t>[]> vReadIDs(new vector<uint64_t>[_iNumOfK]);
				vector<uint64_t> vPositions(_iNumOfK);
				//unique_ptr<uint32_t[]> vReadIDs_(new uint32_t[_iNumOfK * (mReadIDToArrayIdx.size() + 1)]);
				vector<uint64_t> vMemoryOfSeenkMers(_iNumOfK);
				//vector<uint32_t> vMemoryCounterOnly(_iNumOfK, 0);
				vector<Utilities::sBitArray> vMemoryOfTaxIDs(_iNumOfK, Utilities::sBitArray(iSpecIDRange));
				const int32_t& ikDifferenceTop = _iHighestK - _iMaxK;
				const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
				auto seenResultIt = libBeginIt;
				auto rangeBeginIt = libBeginIt, rangeEndIt = libBeginIt;

				tuple<uint64_t, uint64_t> iSeenInput = make_pair(0, 0);
				int16_t ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1);
				int32_t shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);
				const auto shiftVal = [&shift](const uint64_t& val) { return val >> shift; };
				//uint64_t iSeenShiftedkMer = 0;

				//uint32_t iTEST = 0;
				//bool bTEST = false;

				auto vInIndex = vInStart;
				while (vInIndex < vInEnd) {
					// determining the range once is better than checking everytime if it still the same
					uint64_t iSeenRange = get<0>(vIn[vInIndex]);
					uint32_t iRangeLength = get<2>(vIn[vInIndex]);
					uint64_t iInputRangeIdx = vInIndex;

					if (iSeenRange == numeric_limits<uint64_t>::max()) {
						++vInIndex;
						continue;
					}

					while (vInIndex < vInEnd) {
						const uint64_t& currentRange = get<0>(vIn[vInIndex]);
						if (currentRange != iSeenRange && currentRange != numeric_limits<uint64_t>::max()) {
							break;
						}
						else {
							++vInIndex;
						}
					}

					// reset stuff
					for (int32_t i = 0; i < _iNumOfK; ++i) {
						//vReadIDs[i].clear();
						vPositions[i] = 0;
						vMemoryOfSeenkMers[i] = 0;
						vMemoryOfTaxIDs[i].clear();
					}
					iSeenInput = make_pair(0, 0);

					// set range in index
					rangeBeginIt = seenResultIt = libBeginIt + iSeenRange;
					rangeEndIt = libBeginIt + iSeenRange + iRangeLength;

					//auto iStartIdx = iInputRangeIdx;

					bool bDetermineBeginForMatching = true;
					bool bInputIterated = true;

					for (; iInputRangeIdx < vInIndex; ++iInputRangeIdx) {
						// ignore missmatches
						if (get<0>(vIn[iInputRangeIdx]) == numeric_limits<uint64_t>::max()) {
							continue;
						}

						// reset local stuff
						ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1);
						shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);

						// get kmer
						const pair<uint64_t, uint32_t>& iCurrentkMer = make_pair(get<1>(vIn[iInputRangeIdx]), get<3>(vIn[iInputRangeIdx]));
						auto iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;
						bInputIterated = true;
						//if (bTEST) {
							//cout << "input: " << kMerToAminoacid(iCurrentkMer.first, 12) << endl;
						//}

						//if (iCurrentkMer.second == iReadIDTODEBUG) {
						//	cout << iTEST++ << " " << get<0>(vIn[iInputRangeIdx]) << " " << get<2>(vIn[iInputRangeIdx]) << " H:" << vReadIDtoGenID[iReadIDTODEBUG * iSpecIDRange + 1] << " M:" << vReadIDtoGenID[iReadIDTODEBUG * iSpecIDRange + 2] << endl; //424570
							//if (iTEST == 5) {
								/*cout << "input:" << endl;
								for (uint64_t i = iStartIdx; i < vInIndex; ++i) {
									cout << kMerToAminoacid(get<1>(vIn[i]), 12) << endl;
								}
								cout << endl << "index: " << endl;
								for (auto i = rangeBeginIt; i < rangeEndIt + 1; ++i) {
									cout << kMerToAminoacid(i->first, 12) << " " << i->second << " " << (i - libBeginIt) << endl;
								}
								cout << endl;*/
								//	bTEST = true;
								//}
								//else {
								//	bTEST = false;
								//}
							//}

							// determine first occurence inside index to save matching time
						if ((get<0>(iSeenInput) != get<0>(iCurrentkMer)) && (shiftVal(seenResultIt->first) != iCurrentkMerShifted) && bDetermineBeginForMatching) {
							if (shiftVal(rangeBeginIt->first) == iCurrentkMerShifted) {
								seenResultIt = rangeBeginIt;
							}
							else {
								if (shiftVal(rangeEndIt->first) == iCurrentkMerShifted) {
									// end matches but we need the first occurence in the database
									uint64_t iTemp = 1;
									while (shiftVal((rangeEndIt - iTemp)->first) == iCurrentkMerShifted) {
										++iTemp;
									}
									seenResultIt = rangeEndIt - (iTemp - 1);
								}
								else {
									if (iCurrentkMerShifted < shiftVal(rangeBeginIt->first) || iCurrentkMerShifted > shiftVal(rangeEndIt->first)) {
										// kmer not inside index
										continue;
									}
									else {
										// binary search
										seenResultIt = lower_bound(rangeBeginIt, rangeEndIt + 1, iCurrentkMerShifted, [&shift](const decltype(*libBeginIt)& a, const uint64_t& val) { return (a.first >> shift) < val; });
									}
								}
							}
						}
						bDetermineBeginForMatching = false;

						// Now for the real part: Trying to match the kmer with those from the index

						// If the ending is ^ it's not going to hit anyway, might as well stop here
						if ((iCurrentkMerShifted & 31) == 30) {
							continue;
						}


						// Count duplicates that matched too
						if ((get<0>(iSeenInput) == get<0>(iCurrentkMer)) || (seenResultIt == rangeEndIt + 1)) {
							for (int32_t ik = _iNumOfK - 1; ik > -1; --ik) {
								const int32_t& shift_ = 5 * (_iHighestK - _aOfK[ik]);
								const auto& iCurrentkMerShifted_ = get<0>(iCurrentkMer) >> shift_;
								if (iCurrentkMerShifted_ == vMemoryOfSeenkMers[ik]) {
									//vReadIDs[ik].push_back(mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
									//addToMatchedReadID(vReadIDs[ik], vPositions[ik], mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
									addToMatchedReadID(vReadIDs[ik], vPositions[ik], get<1>(iCurrentkMer));
								}
							}
							continue;
						}
						else {
							iSeenInput = iCurrentkMer;
						}



						bool bBreakOut = false;
						while (seenResultIt != rangeEndIt + 1 && !bBreakOut) {

							const tuple<uint64_t, uint32_t>& iCurrentLib = make_tuple(seenResultIt->first, seenResultIt->second);
							//if (bTEST) {
							//	cout << "lib: " << kMerToAminoacid(get<0>(iCurrentLib), 12) << endl;
							//}

							ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1);
							for (; ikLengthCounter >= 0; --ikLengthCounter) {

								shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);
								iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;

								// No matching of stupid stuff!
								/*if ((iCurrentkMerShifted & 31) == 30) {
									bBreakOut = true;
									break;
								}*/


								const auto& iCurrentLibkMerShifted = get<0>(iCurrentLib) >> shift;

								if (iCurrentkMerShifted < iCurrentLibkMerShifted) {
									if (bInputIterated) {
										for (int32_t ik = ikLengthCounter; ik > -1; --ik) {
											const int32_t& shift_ = 5 * (_iHighestK - _aOfK[ik]);
											const auto& iCurrentkMerShifted_ = get<0>(iCurrentkMer) >> shift_;
											if (iCurrentkMerShifted_ == vMemoryOfSeenkMers[ik]) {
												//vReadIDs[ik].push_back(mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
												//addToMatchedReadID(vReadIDs[ik], vPositions[ik], mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
												addToMatchedReadID(vReadIDs[ik], vPositions[ik], get<1>(iCurrentkMer));
											}
											else {
												break;
											}
										}
									}
									bBreakOut = true;
									break;
								}
								else {
									if (!(iCurrentLibkMerShifted < iCurrentkMerShifted)) {

										// If the ending is ^, it's not a valid match, might as well stop here
										if ((iCurrentkMerShifted & 31) == 30) {
											bBreakOut = true;
											break;
										}

										// Delayed scoring: First gather everything and then score it instead of scoring everytime you encounter a new read or tax id.
										if (iCurrentkMerShifted == vMemoryOfSeenkMers[ikLengthCounter]) {
											// We've seen that already, just add it. 
											markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));
											if (bInputIterated) {
												//vReadIDs[ikLengthCounter].push_back(mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
												//addToMatchedReadID(vReadIDs[ikLengthCounter], vPositions[ikLengthCounter], mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
												addToMatchedReadID(vReadIDs[ikLengthCounter], vPositions[ikLengthCounter], get<1>(iCurrentkMer));
											}
										}
										else {
											// For this k, the kmer is different. Save the gathered information.
											//auto valBefore = vReadIDtoGenID[iReadIDTODEBUG * iSpecIDRange + 1];

											const auto& numOfEntries = vMemoryOfTaxIDs[ikLengthCounter].numOfEntries();
											auto it = vMemoryOfTaxIDs[ikLengthCounter].begin();
											it.SetNumOfEntries(numOfEntries);
											for (; it != vMemoryOfTaxIDs[ikLengthCounter].end() && numOfEntries != 0; ++it) {
												const auto& tempIndex = iSpecIDRange * _iNumOfK * iThreadID + iSpecIDRange * ikLengthCounter + (*it);
												const auto& numOfHits = vPositions[ikLengthCounter];//vReadIDs[ikLengthCounter].size();
												//#pragma omp atomic
												vCount[tempIndex] += double(numOfHits) / numOfEntries;

												if (numOfEntries == 1) {
													//#pragma omp atomic
													vCountUnique[tempIndex] += numOfHits;
												}

												const auto& entry = *it;
												const auto& weight = arrWeightingFactors[ikDifferenceTop + ikLengthCounter];
												const auto& score = weight * (1.f / numOfEntries);

												/*for (const auto& readID : vReadIDs[ikLengthCounter]) {
													vReadIDtoGenID[readID * iSpecIDRange + entry] += score;
												}*/
												for (uint64_t pos = 0; pos < numOfHits; ++pos) {
													vReadIDtoGenID[vReadIDs[ikLengthCounter][pos] * iSpecIDRange + entry] += score;
												}
											}

											//auto valAfter = vReadIDtoGenID[iReadIDTODEBUG * iSpecIDRange + 1];
											//if (valBefore != valAfter) {
											//	cout << ikLengthCounter << " " << vReadIDtoGenID[iReadIDTODEBUG * iSpecIDRange + 1] << " " << numOfEntries << " " << vReadIDs[ikLengthCounter].size() << endl;
												/*if (vReadIDs[ikLengthCounter].size() > 2) {
													for (const auto& entry : vReadIDs[ikLengthCounter]) {
														cout << entry << endl;
													}
												}*/
												//}

											//vReadIDs[ikLengthCounter].clear();
											//vReadIDs[ikLengthCounter].push_back(mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
											vPositions[ikLengthCounter] = 0;
											//addToMatchedReadID(vReadIDs[ikLengthCounter], vPositions[ikLengthCounter], mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
											addToMatchedReadID(vReadIDs[ikLengthCounter], vPositions[ikLengthCounter], get<1>(iCurrentkMer));

											vMemoryOfTaxIDs[ikLengthCounter].clear();
											markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));

											vMemoryOfSeenkMers[ikLengthCounter] = iCurrentkMerShifted;
										}

									}
									else {

										// index kmer is smaller than input kmer. Iterate linearly and score everything that also matched before (albeit with a smaller k)
										uint64_t iTempCounter = 1;
										while (seenResultIt + iTempCounter != rangeEndIt + 1) {
											const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + iTempCounter)->first);
											if (iCurrentkMerShifted > (iNextLibSuffix >> shift)) {
												int16_t iUntilK = static_cast<int16_t>(_iNumOfK - 1);
												for (; iUntilK > -1; --iUntilK) {
													if (vMemoryOfSeenkMers[iUntilK] == (iNextLibSuffix >> 5 * (_iHighestK - _aOfK[iUntilK]))) {
														markTaxIDs((seenResultIt + iTempCounter)->second, vMemoryOfTaxIDs[iUntilK], mTaxToIdx, getVec(vLib, iThreadID));
													}
													else {
														break;
													}
												}
												if (iUntilK < static_cast<int16_t>(_iNumOfK - 1)) {
													++iTempCounter;
												}
												else {
													break;
												}
											}
											else {
												break;
											}
										}
										seenResultIt += iTempCounter;

										//bInputIterated = false;
										//seenResultIt = lower_bound(seenResultIt, rangeEndIt + 1, iCurrentkMerShifted, [&shift, this](const decltype(*libBeginIt)& a, const uint64_t& val) { return (_iMinK > 6) ? ((a.first & 1073741823ULL) >> shift) < val : (a.first >> shift) < val; });
										break;
									}
								}
							}
							// loop through to find other hits in the library (index is redundant so equal kmers with different tax ids are possible)
							if (ikLengthCounter == -1) {
								++seenResultIt;
							}

							bInputIterated = false;
						}
					}

					// look through the rest of the range for possible tax id matches

					uint64_t iTempCounter = 0;
					while (seenResultIt + iTempCounter != rangeEndIt + 1) {
						const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + iTempCounter)->first);
						//if (bTEST) {
						//	cout << "libEnd: " << kMerToAminoacid(iNextLibSuffix, 12) << endl;
						//}
						int16_t iUntilK = static_cast<int16_t>(_iNumOfK - 1);
						for (; iUntilK > -1; --iUntilK) {
							if (vMemoryOfSeenkMers[iUntilK] == (iNextLibSuffix >> 5 * (_iHighestK - _aOfK[iUntilK]))) {
								markTaxIDs((seenResultIt + iTempCounter)->second, vMemoryOfTaxIDs[iUntilK], mTaxToIdx, getVec(vLib, iThreadID));
							}
							else {
								break;
							}
						}
						if (iUntilK < static_cast<int16_t>(_iNumOfK - 1)) {
							++iTempCounter;
						}
						else {
							break;
						}
					}

					// save everything from this range
					for (int16_t ikLC = static_cast<int16_t>(_iNumOfK - 1); ikLC >= 0; --ikLC) {

						//auto valBefore = vReadIDtoGenID[iReadIDTODEBUG * iSpecIDRange + 1];

						const auto& numOfEntries = vMemoryOfTaxIDs[ikLC].numOfEntries();
						auto it = vMemoryOfTaxIDs[ikLC].begin();
						it.SetNumOfEntries(numOfEntries);
						for (; it != vMemoryOfTaxIDs[ikLC].end() && numOfEntries != 0; ++it) {
							const auto& tempIndex = iSpecIDRange * _iNumOfK * iThreadID + iSpecIDRange * ikLC + (*it);
							const auto& numOfHits = vPositions[ikLC];//vReadIDs[ikLC].size();
							//#pragma omp atomic
							vCount[tempIndex] += double(numOfHits) / numOfEntries;

							if (numOfEntries == 1) {
								//#pragma omp atomic
								vCountUnique[tempIndex] += numOfHits;
							}

							const auto& entry = *it;
							const auto& weight = arrWeightingFactors[ikDifferenceTop + ikLC];
							const auto& score = weight * (1.f / numOfEntries);
							for (uint64_t pos = 0; pos < numOfHits; ++pos) {
								vReadIDtoGenID[vReadIDs[ikLC][pos] * iSpecIDRange + entry] += score;
							}

							/*
							for (const auto& readID : vReadIDs[ikLC]) {
								vReadIDtoGenID[readID * iSpecIDRange + entry] += score;
							}
							*/
						}

						//auto valAfter = vReadIDtoGenID[iReadIDTODEBUG * iSpecIDRange + 1];
						//if (valBefore != valAfter) {
						//	cout << "eor " << ikLC << " " << vReadIDtoGenID[iReadIDTODEBUG * iSpecIDRange + 1] << " " << numOfEntries << " " << vReadIDs[ikLC].size() << endl;
							/*if (vReadIDs[ikLC].size() > 2) {
								for (const auto& entry : vReadIDs[ikLC]) {
									cout << entry << endl;
								}
							}*/
							//}
					}
				}

				// Don't forget the last saved part
				/*for (int16_t ikL = static_cast<int16_t>(_iNumOfK - 1); ikL >= 0; --ikL) {
					const auto& numOfEntries = vMemoryOfTaxIDs[ikL].numOfEntries();
					auto it = vMemoryOfTaxIDs[ikL].begin();
					it.SetNumOfEntries(numOfEntries);
					for (; it != vMemoryOfTaxIDs[ikL].end() && numOfEntries != 0; ++it) {
						const auto& tempIndex = iSpecIDRange * _iNumOfK * iThreadID + iSpecIDRange * ikL + (*it);
						const auto& numOfHits = vReadIDs[ikL].size();
						vCount[tempIndex] += float(numOfHits) / numOfEntries;

						if (numOfEntries == 1) {
							vCountUnique[tempIndex] += numOfHits;
						}

						const auto& entry = *it;
						const auto& weight = arrWeightingFactors[ikDifferenceTop + ikL];
						const auto& score = weight * (1.f / numOfEntries);

						for (const auto& readID : vReadIDs[ikL]) {
							vReadIDtoGenID[readID * iSpecIDRange + entry] += score;
						}
					}
				}*/
			}
			catch (...) {
				if (someThingWentWrong == nullptr) {
					m_exceptionLock.lock();
					someThingWentWrong = current_exception();
					m_exceptionLock.unlock();
				}
				return;
			}

		}

		//#define DEBUGOUT
#ifdef DEBUGOUT
		int64_t amountArr[5] = { 0,0,0,0,0 };
		int64_t amountTaxArr[5] = { 0,0,0,0,0 };
#endif

		template <typename vecType>
		inline void createProfile(const int32_t iThreadID, const vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>>& vIn, const uint64_t& vInStart, const uint64_t& vInEnd, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, const uint32_t& iSpecIDRange, const unordered_map<uint32_t, uint32_t>& mTaxToIdx) {
			try {
				vector<uint64_t> vMemoryOfSeenkMers(_iNumOfK);
				vector<uint32_t> vMemoryCounterOnly(_iNumOfK, 0);
				vector<Utilities::sBitArray> vMemoryOfTaxIDs(_iNumOfK, Utilities::sBitArray(iSpecIDRange));
				const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
				auto seenResultIt = libBeginIt;
				auto rangeBeginIt = libBeginIt, rangeEndIt = libBeginIt;

				tuple<uint64_t, uint64_t> iSeenInput = make_pair(0, 0);
				int16_t ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1);
				int32_t shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);
				const auto shiftVal = [&shift](const uint64_t& val) { return val >> shift; };

				auto vInIndex = vInStart;
				while (vInIndex < vInEnd) {
					// determining the range once is better than checking everytime if it still the same
					uint64_t iSeenRange = get<0>(vIn[vInIndex]);
					uint32_t iRangeLength = get<2>(vIn[vInIndex]);
					uint64_t iInputRangeIdx = vInIndex;

					if (iSeenRange == numeric_limits<uint64_t>::max()) {
						++vInIndex;
						continue;
					}

					while (vInIndex < vInEnd) {
						const uint64_t& currentRange = get<0>(vIn[vInIndex]);
						if (currentRange != iSeenRange && currentRange != numeric_limits<uint64_t>::max()) {
							break;
						}
						else {
							++vInIndex;
						}
					}

					// reset stuff
					for (int32_t i = 0; i < _iNumOfK; ++i) {
						vMemoryCounterOnly[i] = 0;
						vMemoryOfSeenkMers[i] = 0;
						vMemoryOfTaxIDs[i].clear();
					}
					iSeenInput = make_pair(0, 0);

					// set range in index
					rangeBeginIt = seenResultIt = libBeginIt + iSeenRange;
					rangeEndIt = libBeginIt + iSeenRange + iRangeLength;

#ifdef DEBUGOUT
					cout << "Range: " << endl;
					for (auto rangeIt = rangeBeginIt; rangeIt != rangeEndIt + 1; ++rangeIt) {
						cout << kMerToAminoacid(rangeIt->first, 12) << " " << rangeIt->second << endl;
					}
					cout << endl;
#endif
					bool bDetermineBeginForMatching = true;
					bool bInputIterated = true;

					for (; iInputRangeIdx < vInIndex; ++iInputRangeIdx) {
						// ignore missmatches
						if (get<0>(vIn[iInputRangeIdx]) == numeric_limits<uint64_t>::max()) {
							continue;
						}


						//bool bMatched = false;


						// reset local stuff
						ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1);
						shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);

						// get kmer
						const pair<uint64_t, uint32_t>& iCurrentkMer = make_pair(get<1>(vIn[iInputRangeIdx]), get<3>(vIn[iInputRangeIdx]));
						auto iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;
						bInputIterated = true;
#ifdef DEBUGOUT
						cout << "Input: " << kMerToAminoacid(get<0>(iCurrentkMer), 12) << " " << get<1>(iCurrentkMer) << endl;
#endif
						// determine first occurence inside index to save matching time
						if ((get<0>(iSeenInput) != get<0>(iCurrentkMer)) && (shiftVal(seenResultIt->first) != iCurrentkMerShifted) && bDetermineBeginForMatching) {
							if (shiftVal(rangeBeginIt->first) == iCurrentkMerShifted) {
								seenResultIt = rangeBeginIt;
							}
							else {
								if (shiftVal(rangeEndIt->first) == iCurrentkMerShifted) {
									// end matches but we need the first occurence in the database
									uint64_t iTemp = 1;
									while (shiftVal((rangeEndIt - iTemp)->first) == iCurrentkMerShifted) {
										++iTemp;
									}
									seenResultIt = rangeEndIt - (iTemp - 1);
								}
								else {
									if (iCurrentkMerShifted < shiftVal(rangeBeginIt->first) || iCurrentkMerShifted > shiftVal(rangeEndIt->first)) {
										// kmer not inside index
										continue;
									}
									else {
										// binary search
										seenResultIt = lower_bound(rangeBeginIt, rangeEndIt + 1, iCurrentkMerShifted, [&shift](const decltype(*libBeginIt)& a, const uint64_t& val) { return (a.first >> shift) < val; });
									}
								}
							}
						}
						bDetermineBeginForMatching = false;

						// Now for the real part: Trying to match the kmer with those from the index

						// If the ending is ^ it's not going to hit anyway, might as well stop here
						if ((iCurrentkMerShifted & 31) == 30) {
							continue;
						}


						// Count duplicates that matched too
						if ((get<0>(iSeenInput) == get<0>(iCurrentkMer)) || (seenResultIt == rangeEndIt + 1)) {
							for (int32_t ik = _iNumOfK - 1; ik > -1; --ik) {
								const int32_t& shift_ = 5 * (_iHighestK - _aOfK[ik]);
								const auto& iCurrentkMerShifted_ = get<0>(iCurrentkMer) >> shift_;
								if (iCurrentkMerShifted_ == vMemoryOfSeenkMers[ik]) {
#ifdef DEBUGOUT
									if (ik == _iNumOfK - 1) {
										amountArr[0]++;
									}
#endif
									++vMemoryCounterOnly[ik];
									//bMatched = true;
								}
							}
							continue;
						}
						else {
							iSeenInput = iCurrentkMer;
						}


						/*if (vInIndex - iInputRangeIdx > 1) {
							cout << vInIndex - iInputRangeIdx << endl;
						}*/

						bool bBreakOut = false;
						while (seenResultIt != rangeEndIt + 1 && !bBreakOut) {

							const tuple<uint64_t, uint32_t>& iCurrentLib = make_tuple(seenResultIt->first, seenResultIt->second);

							ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1);
							for (; ikLengthCounter >= 0; --ikLengthCounter) {

								shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);
								iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;

								const auto& iCurrentLibkMerShifted = get<0>(iCurrentLib) >> shift;
#ifdef DEBUGOUT
								cout << kMerToAminoacid(get<0>(iCurrentkMer), 12) << " " << kMerToAminoacid(get<0>(iCurrentLib), 12) << " " << kMerToAminoacid(iCurrentkMerShifted, 12) << " " << kMerToAminoacid(iCurrentLibkMerShifted, 12) << endl;
#endif
								if (iCurrentkMerShifted < iCurrentLibkMerShifted) {
									if (bInputIterated) {
										for (int32_t ik = ikLengthCounter; ik > -1; --ik) {
											const int32_t& shift_ = 5 * (_iHighestK - _aOfK[ik]);
											const auto& iCurrentkMerShifted_ = get<0>(iCurrentkMer) >> shift_;
											if (iCurrentkMerShifted_ == vMemoryOfSeenkMers[ik]) {
#ifdef DEBUGOUT
												if (ik == _iNumOfK - 1) {
													amountArr[1]++;
												}
#endif
												++vMemoryCounterOnly[ik];
												//bMatched = true;
											}
											else {
												break;
											}
										}
									}
									bBreakOut = true;
									break;
								}
								else {
									if (!(iCurrentLibkMerShifted < iCurrentkMerShifted)) {

										// If the ending is ^, it's not a valid match, might as well stop here
										if ((iCurrentkMerShifted & 31) == 30) {
											bBreakOut = true;
											break;
										}

										// Delayed scoring: First gather everything and then score it instead of scoring everytime you encounter a new read or tax id.
										if (iCurrentkMerShifted == vMemoryOfSeenkMers[ikLengthCounter]) {
											// We've seen that already, just add it. 
											markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));

#ifdef DEBUGOUT
											if (ikLengthCounter == _iNumOfK - 1) {
												amountTaxArr[0]++;
											}
#endif
											if (bInputIterated) {
#ifdef DEBUGOUT
												if (ikLengthCounter == _iNumOfK - 1) {
													amountArr[2]++;
												}
#endif
												++vMemoryCounterOnly[ikLengthCounter];
												//bMatched = true;
											}
										}
										else {
											// For this k, the kmer is different. Save the gathered information.
											const auto& numOfEntries = vMemoryOfTaxIDs[ikLengthCounter].numOfEntries();
											auto it = vMemoryOfTaxIDs[ikLengthCounter].begin();
											it.SetNumOfEntries(numOfEntries);
											for (; it != vMemoryOfTaxIDs[ikLengthCounter].end() && numOfEntries != 0; ++it) {
												const auto& tempIndex = iSpecIDRange * _iNumOfK * iThreadID + iSpecIDRange * ikLengthCounter + (*it);
												const auto& numOfHits = vMemoryCounterOnly[ikLengthCounter];
												vCount[tempIndex] += double(numOfHits) / numOfEntries;

												if (numOfEntries == 1) {
													vCountUnique[tempIndex] += numOfHits;
												}
											}

											vMemoryCounterOnly[ikLengthCounter] = 1;
											//bMatched = true;
#ifdef DEBUGOUT
											if (ikLengthCounter == _iNumOfK - 1) {
												amountArr[3]++;
											}
#endif

											vMemoryOfTaxIDs[ikLengthCounter].clear();
											markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));
#ifdef DEBUGOUT
											if (ikLengthCounter == _iNumOfK - 1) {
												amountTaxArr[1]++;
											}
#endif
											vMemoryOfSeenkMers[ikLengthCounter] = iCurrentkMerShifted;
										}

									}
									else {

										// index kmer is smaller than input kmer. Iterate linearly and score everything that also matched before (albeit with a smaller k)
										uint64_t iTempCounter = 1;
										while (seenResultIt + iTempCounter != rangeEndIt + 1) {
											const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + iTempCounter)->first);
#ifdef DEBUGOUT
											cout << kMerToAminoacid(iNextLibSuffix, 12) << endl;
#endif							
											if (iCurrentkMerShifted > (iNextLibSuffix >> shift)) {
												int16_t iUntilK = static_cast<int16_t>(_iNumOfK - 1);
												for (; iUntilK > -1; --iUntilK) {
													if (vMemoryOfSeenkMers[iUntilK] == (iNextLibSuffix >> 5 * (_iHighestK - _aOfK[iUntilK]))) {
														markTaxIDs((seenResultIt + iTempCounter)->second, vMemoryOfTaxIDs[iUntilK], mTaxToIdx, getVec(vLib, iThreadID));
#ifdef DEBUGOUT
														if (iUntilK == _iNumOfK - 1) {
															amountTaxArr[2]++;
														}
#endif
													}
													else {
														break;
													}
												}
												if (iUntilK < static_cast<int16_t>(_iNumOfK - 1)) {
													++iTempCounter;
												}
												else {
													break;
												}
											}
											else {
												break;
											}
										}
										seenResultIt += iTempCounter;

										//bInputIterated = false;
										break;
									}
								}
							}
							// loop through to find other hits in the library (index is redundant so equal kmers with different tax ids are possible)
							if (ikLengthCounter == -1) {
								++seenResultIt;
								/*const uint64_t& iCurrentLibSuffix = static_cast<uint64_t>(seenResultIt->first);
								while (seenResultIt + 1 != rangeEndIt + 1) {
									const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + 1)->first);
									const auto& iNextLibIdx = (seenResultIt + 1)->second;
#ifdef DEBUGOUT
									cout << kMerToAminoacid(iNextLibSuffix, 12) << endl;
#endif
									if ((iCurrentLibSuffix >> 5 * (_iHighestK - _aOfK[0])) == (iNextLibSuffix >> 5 * (_iHighestK - _aOfK[0]))) {
										for (int16_t ikLengthCounter_ = static_cast<int16_t>(_iNumOfK - 1); ikLengthCounter_ > -1; --ikLengthCounter_) {
											markTaxIDs(iNextLibIdx, vMemoryOfTaxIDs[ikLengthCounter_], mTaxToIdx, getVec(vLib, iThreadID)); // to identify multiple hits
											if (ikLengthCounter_ == _iNumOfK - 1) {
												amountTaxArr[3]++;
											}
										}

										++seenResultIt;
									}
									else {
										++seenResultIt;
										break;
									}
								}
								bBreakOut = true;
								break;
								//bInputIterated = false;*/
							}

							bInputIterated = false;
						}

#ifdef DEBUGOUT
						//if (bMatched == false) {
						//	cout << kMerToAminoacid(iCurrentkMer.first, 12) << endl;
						//}
#endif
					}

					// look through the rest of the range for possible tax id matches
					uint64_t iTempCounter = 0;
					while (seenResultIt + iTempCounter != rangeEndIt + 1) {
						const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + iTempCounter)->first);
#ifdef DEBUGOUT
						cout << kMerToAminoacid(iNextLibSuffix, 12) << endl;
#endif		
						int16_t iUntilK = static_cast<int16_t>(_iNumOfK - 1);
						for (; iUntilK > -1; --iUntilK) {
							if (vMemoryOfSeenkMers[iUntilK] == (iNextLibSuffix >> 5 * (_iHighestK - _aOfK[iUntilK]))) {
								markTaxIDs((seenResultIt + iTempCounter)->second, vMemoryOfTaxIDs[iUntilK], mTaxToIdx, getVec(vLib, iThreadID));
#ifdef DEBUGOUT
								if (iUntilK == _iNumOfK - 1) {
									amountTaxArr[4]++;
								}
#endif
							}
							else {
								break;
							}
						}
						if (iUntilK < static_cast<int16_t>(_iNumOfK - 1)) {
							++iTempCounter;
						}
						else {
							break;
						}
					}

					// save everything from this range
					for (int16_t ikLC = static_cast<int16_t>(_iNumOfK - 1); ikLC >= 0; --ikLC) {

						const auto& numOfEntries = vMemoryOfTaxIDs[ikLC].numOfEntries();
						auto it = vMemoryOfTaxIDs[ikLC].begin();
						it.SetNumOfEntries(numOfEntries);
						for (; it != vMemoryOfTaxIDs[ikLC].end() && numOfEntries != 0; ++it) {
							const auto& tempIndex = iSpecIDRange * _iNumOfK * iThreadID + iSpecIDRange * ikLC + (*it);
							const auto& numOfHits = vMemoryCounterOnly[ikLC];
							vCount[tempIndex] += double(numOfHits) / numOfEntries;

							if (numOfEntries == 1) {
								vCountUnique[tempIndex] += numOfHits;
							}
						}
					}
				}

			}
			catch (...) {
				if (someThingWentWrong == nullptr) {
					m_exceptionLock.lock();
					someThingWentWrong = current_exception();
					m_exceptionLock.unlock();
				}
				return;
			}
		}

		//////////////////////////////////////////////////
		// For the sloppy index

		inline const vector<uint16_t>* getVec(const vector<uint16_t>* vec, const int&) {
			return vec;
		}
		inline unique_ptr<const taxaOnly>& getVec(const unique_ptr<unique_ptr<const taxaOnly>[]>& vec, const int& iThreadID) {
			return vec[iThreadID];
		}

		template <typename vecType>
		inline void compareWithDatabase_sloppy(const int32_t& iThreadID, const vector<pair<uint64_t, Utilities::rangeContainer>>& vIn, const uint64_t& vInStart, const uint64_t& vInEnd, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, unique_ptr<float[]>& vReadIDtoGenID, const uint32_t& iSpecIDRange, const unordered_map<readIDType, uint64_t>& mReadIDToArrayIdx) {

			try {

				//auto start = std::chrono::high_resolution_clock::now();
				for (auto vInIndex = vInStart; vInIndex < vInEnd; ++vInIndex) {

					const uint64_t& ivInSize = vIn[vInIndex].second.kMers_ST6.size();

					const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
					auto seenResultIt = libBeginIt;

					const auto rangeBeginIt = libBeginIt + vIn[vInIndex].first, rangeEndIt = libBeginIt + static_cast<uint64_t>(vIn[vInIndex].first) + vIn[vInIndex].second.range;
					seenResultIt = rangeBeginIt;

					while (seenResultIt != rangeEndIt + 1) {
						const auto& iCurrentLib = *seenResultIt;

						vCount[iSpecIDRange*iThreadID + iCurrentLib] += static_cast<double>(ivInSize) / (vIn[vInIndex].second.range + 1);


						if (vIn[vInIndex].second.range == 0) {
							vCountUnique[iSpecIDRange*iThreadID + iCurrentLib] += ivInSize;
						}

						const auto& score = arrWeightingFactors[0] * (1.f / ivInSize);
						for (const auto& elem : vIn[vInIndex].second.kMers_ST6) {
							const auto& arrayIdx = mReadIDToArrayIdx.find(get<1>(elem))->second * iSpecIDRange + iCurrentLib;
							//#pragma omp atomic
							vReadIDtoGenID[arrayIdx] += score;

						}

						++seenResultIt;
					}

					/*auto end = std::chrono::high_resolution_clock::now();
					criticalSectionLock.lock();
					cout << "OUT: Time: " << iThreadID << " " <<  chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << endl;
					criticalSectionLock.unlock();*/
				}
			}
			catch (...) {
				if (someThingWentWrong == nullptr) {
					m_exceptionLock.lock();
					someThingWentWrong = current_exception();
					m_exceptionLock.unlock();
				}
				return;
			}

		}

		template <typename vecType>
		inline void createProfile_sloppy(const int32_t iThreadID, const vector<pair<uint64_t, Utilities::rangeContainer>>& vIn, const uint64_t& vInStart, const uint64_t& vInEnd, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, const uint32_t& iSpecIDRange) {

			try {
				for (auto vInIndex = vInStart; vInIndex < vInEnd; ++vInIndex) {

					const uint64_t& ivInSize = vIn[vInIndex].second.kMers_ST6.size();

					const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
					auto seenResultIt = libBeginIt;

					const auto rangeBeginIt = libBeginIt + vIn[vInIndex].first, rangeEndIt = libBeginIt + static_cast<uint64_t>(vIn[vInIndex].first) + vIn[vInIndex].second.range;
					seenResultIt = rangeBeginIt;

					while (seenResultIt != rangeEndIt + 1) {
						const auto& iCurrentLib = *seenResultIt;
						//#pragma omp atomic
						vCount[iSpecIDRange*iThreadID + iCurrentLib] += static_cast<double>(ivInSize) / (vIn[vInIndex].second.range + 1);


						if (vIn[vInIndex].second.range == 0) {
							//#pragma omp atomic
							vCountUnique[iSpecIDRange*iThreadID + iCurrentLib] += ivInSize;
						}
						++seenResultIt;
					}
				}
			}
			catch (...) {
				if (someThingWentWrong == nullptr) {
					m_exceptionLock.lock();
					someThingWentWrong = current_exception();
					m_exceptionLock.unlock();
				}
				return;
			}

		}


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void scoringFunc(const unique_ptr<float[]>& vReadIDtoGenID, const uint64_t& iStart, const uint64_t& iEnd, const uint64_t& iRealReadIDStart, list<pair<string, readIDType>>& vReadNameAndLength, const unique_ptr<uint64_t[]>& mFrequencies, const vector<uint32_t>& mIdxToTax, const vector<string>& mOrganisms, const float& fThreshold, ofstream&& fOut) {
			try {

				const size_t& iAmountOfSpecies = mIdxToTax.size();
				vector<tuple<size_t, float, double>> resultVec(iAmountOfSpecies);

				for (uint64_t readIdx = iStart; readIdx < iEnd; ++readIdx) {
					auto currentReadLengthAndName = vReadNameAndLength.front();

					float bestScore = 0.f;
					for (int32_t i = _iMinK; i <= _iMaxK; ++i) {
						bestScore += (currentReadLengthAndName.second - i * 3 + 1)*arrWeightingFactors[_iHighestK - i];
					}


					float kMerScore = 0.f;
					double relativeScore = 0.0;
					int16_t iCountOfHits = 0;
					for (size_t iSpecIdx = 1; iSpecIdx < iAmountOfSpecies; ++iSpecIdx) {
						if (vReadIDtoGenID[readIdx*iAmountOfSpecies + iSpecIdx] > 0.f) {
							kMerScore = vReadIDtoGenID[readIdx*iAmountOfSpecies + iSpecIdx];

							if (_bTranslated) {
								relativeScore = kMerScore / (1.0 + log2(mFrequencies[iSpecIdx] * double(currentReadLengthAndName.second - _iHighestK + 1)));
							}
							else {
								relativeScore = kMerScore / (1.0 + log2(mFrequencies[iSpecIdx] * double(currentReadLengthAndName.second - _iHighestK * 3 + 1)));
							}
							if (relativeScore >= fThreshold) {
								get<0>(resultVec[iSpecIdx]) = iSpecIdx;
								get<1>(resultVec[iSpecIdx]) = kMerScore;
								get<2>(resultVec[iSpecIdx]) = relativeScore;
								++iCountOfHits;
							}
						}
					}

					if (iCountOfHits == 0) {
						switch (format) {
						case OutputFormat::tsv:
							fOut << iRealReadIDStart + readIdx << "\t" << currentReadLengthAndName.first << "\t-\t-\t-" << "\n";
							break;

						case OutputFormat::Json:
							if (iRealReadIDStart + readIdx == 0) {
								fOut << "{" << "\n";
							}
							else {
								fOut << "," << "\n" << "{" << "\n";
							}
							fOut << "\t\"Read number\": " << iRealReadIDStart + readIdx << ",\n\t\"Specifier from input file\": \"" + currentReadLengthAndName.first + "\",\n\t\"Top hits\": [\n\t],\n\t\"Further hits\": [\n\t]\n}";
							break;

						case OutputFormat::JsonL:
							fOut << "{ \"Read number\": " << iRealReadIDStart + readIdx << "," << " \"Specifier from input file\": \"" + currentReadLengthAndName.first + "\", " << "\"Top hits\": [], \"Further hits\": [] }\n";
							break;

						case OutputFormat::Kraken:
							fOut << "U\t" << currentReadLengthAndName.first << "\t0\t" << currentReadLengthAndName.second << "\tA:00\n";
							break;
						};

					}
					else {

						//partial_sort(resultVec.begin(), resultVec.begin() + iCountOfHits, resultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<2>(a) > get<2>(b); });
						//sort(std::execution::par_unseq, resultVec.begin(), resultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<2>(a) > get<2>(b); });
# if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
						__gnu_parallel::sort(resultVec.begin(), resultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<2>(a) > get<2>(b); }, __gnu_parallel::balanced_quicksort_tag(_iNumOfThreads));
#else					
#if __has_include(<execution>)
						sort(std::execution::par_unseq, resultVec.begin(), resultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<2>(a) > get<2>(b); });
#else
						sort(resultVec.begin(), resultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<2>(a) > get<2>(b); });
#endif
#endif

						auto maxValue = get<1>(*(max_element(resultVec.begin(), resultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<1>(a) < get<1>(b); })));
						int16_t iTopHitCounter = 1;
						for (int16_t i = 1; i < iCountOfHits && i < iNumOfBeasts; ++i) {
							if ((get<1>(resultVec[i])) / (maxValue) > 0.8f) { // determine "close enough" with normalization relative to the highest score
								++iTopHitCounter;
							}
							else {
								break;
							}
						}

						//cout << iRealReadIDStart + readIdx << endl;
						string sOut = "", sOut2 = "", sOut3 = "";
						sOut.reserve(1000);
						auto it = resultVec.begin();
						float iValueBefore = 0;

						switch (format) {
						case OutputFormat::tsv:
							// tsv
							// += is faster than () + ()
							Utilities::itostr(iRealReadIDStart + readIdx, sOut);
							sOut += "\t";
							sOut += currentReadLengthAndName.first + "\t";

							for (int16_t j = 0, i = 0; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
								sOut += ";";
								sOut2 += mOrganisms[get<0>(*it)];
								sOut2 += ";";


								dtoa_milo(get<2>(*it), sOut3);
								sOut3 += ",";
								dtoa_milo(get<1>(*it), sOut3);
								sOut3 += ";";

								if (iValueBefore != get<1>(*it)) {
									iValueBefore = get<1>(*it);
									++j;
								}
							}
							if (sOut.back() == ';') {
								sOut.pop_back();
							}
							if (sOut2.back() == ';') {
								sOut2.pop_back();
							}
							if (sOut3.back() == ';') {
								sOut3.pop_back();
							}
							if (sOut2.length()) {
								fOut << sOut + "\t" + sOut2 + "\t" + sOut3 + "\t" << (bestScore - get<1>(resultVec[0])) / bestScore << "\n";
							}
							break;

						case OutputFormat::Json:
							// json
					 
							if (iRealReadIDStart + readIdx == 0) {
								sOut += "{\n";
							}
							else {
								sOut += ",\n{\n";
							}

							sOut += "\t\"Read number\": ";
							Utilities::itostr(iRealReadIDStart + readIdx, sOut);
							sOut += ",\n\t\"Specifier from input file\": \"";
							sOut += currentReadLengthAndName.first;
							sOut += "\",\n\t\"Top hits\": [\n";

							for (int16_t i = 0; i < iTopHitCounter; ++i, ++it) {
								if (i == 0) {
									sOut += "\t{\n";
								}
								else {
									sOut += ",\n\t{\n";
								}

								sOut += "\t\t\"tax ID\": \"";
								Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
								sOut += "\",\n";
								sOut += "\t\t\"Name\": \"";
								sOut += mOrganisms[get<0>(*it)];
								sOut += "\",\n";
								sOut += "\t\t\"k-mer Score\": ";
								dtoa_milo(get<1>(*it), sOut);
								sOut += ",\n";
								sOut += "\t\t\"Relative Score\": ";
								dtoa_milo(get<2>(*it), sOut);
								sOut += ",\n";
								sOut += "\t\t\"Error\": ";
								dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut);
								sOut += "\n\t}";
							}

							sOut += "\n\t],\n\t\"Further hits\": [\n";

							for (int16_t j = iTopHitCounter, i = iTopHitCounter; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								if (j == iTopHitCounter) {
									sOut += "\t{\n";
								}
								else {
									sOut += ",\n\t{\n";
								}

								sOut += "\t\t\"tax ID\": \"";
								Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
								sOut += "\",\n";
								sOut += "\t\t\"Name\": \"";
								sOut += mOrganisms[get<0>(*it)];
								sOut += "\",\n";
								sOut += "\t\t\"k-mer Score\": ";
								dtoa_milo(get<1>(*it), sOut);
								sOut += ",\n";
								sOut += "\t\t\"Relative Score\": ";
								dtoa_milo(get<2>(*it), sOut);
								sOut += ",\n";
								sOut += "\t\t\"Error\": ";
								dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut);
								sOut += "\n\t}";

								if (iValueBefore != get<1>(*it)) {
									iValueBefore = get<1>(*it);
									++j;
								}
							}

							sOut += "\n\t]\n}";
							fOut << sOut;
							break;

						case OutputFormat::JsonL:
							// json lines
							sOut += "{ \"Read number\": ";
							Utilities::itostr(iRealReadIDStart + readIdx, sOut);
							sOut += ", \"Specifier from input file\": \"";
							sOut += currentReadLengthAndName.first;
							sOut += "\", \"Top hits\": [";

							for (int16_t i = 0; i < iTopHitCounter; ++i, ++it) {
								if (i == 0) {
									sOut += "{";
								}
								else {
									sOut += ",{";
								}

								sOut += " \"tax ID\": \"";
								Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
								sOut += "\",";
								sOut += " \"Name\": \"";
								sOut += mOrganisms[get<0>(*it)];
								sOut += "\",";
								sOut += " \"k-mer Score\": ";
								dtoa_milo(get<1>(*it), sOut);
								sOut += ",";
								sOut += " \"Relative Score\": ";
								dtoa_milo(get<2>(*it), sOut);
								sOut += ",";
								sOut += " \"Error\": ";
								dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut);
								sOut += "}";
							}

							sOut += "], \"Further hits\": [";

							for (int16_t j = iTopHitCounter, i = iTopHitCounter; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								if (j == iTopHitCounter) {
									sOut += "{";
								}
								else {
									sOut += ", {";
								}

								sOut += " \"tax ID\": \"";
								Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
								sOut += "\",";
								sOut += " \"Name\": \"";
								sOut += mOrganisms[get<0>(*it)];
								sOut += "\",";
								sOut += " \"k-mer Score\": "; 
								dtoa_milo(get<1>(*it), sOut);
								sOut += ",";
								sOut += " \"Relative Score\": "; 
								dtoa_milo(get<2>(*it), sOut);
								sOut += ",";
								sOut += " \"Error\": ";
								dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut);
								sOut += "}";

								if (iValueBefore != get<1>(*it)) {
									iValueBefore = get<1>(*it);
									++j;
								}
							}

							sOut += "] }\n";
							fOut << sOut;
							break;

						case OutputFormat::Kraken:
							// Kraken
							sOut += "C\t";
							sOut += currentReadLengthAndName.first;
							sOut += "\t";
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += "\t";
							Utilities::itostr(currentReadLengthAndName.second, sOut);
							sOut += "\t";

							for (int16_t i = 0; i < iTopHitCounter; ++i, ++it) {
								Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
								sOut += ":";
								dtoa_milo(get<1>(*it), sOut);
								sOut += " ";
							}

							for (int16_t j = iTopHitCounter, i = iTopHitCounter; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
								sOut += ":";
								dtoa_milo(get<1>(*it), sOut);
								sOut += " ";
								if (iValueBefore != get<1>(*it)) {
									iValueBefore = get<1>(*it);
									++j;
								}
							}

							fOut << sOut << "\n";
							break;
						};


						for (auto& entry : resultVec) {
							get<0>(entry) = 0;
							get<1>(entry) = 0.f;
							get<2>(entry) = 0.0;
						}
					}

					vReadNameAndLength.pop_front();
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl;
				throw;
			}
		}


		inline void scoringFunc(vector<tuple<readIDType, float, double>>&& vTempResultVec, const uint64_t& iReadNum, const pair<string, uint32_t>& vReadNameAndLength, const unique_ptr<uint64_t[]>& mFrequencies, const vector<uint32_t>& mIdxToTax, const vector<string>& mOrganisms, const float& fThreshold, ofstream&& fOut) {
			try {
				float bestScore = 0.f;
				for (int32_t i = _iMinK; i <= _iMaxK; ++i) {
					bestScore += (vReadNameAndLength.second - i * 3 + 1)*arrWeightingFactors[_iHighestK - i];
					//cout << (vReadNameAndLength.second - i * 3 + 1)*arrWeightingFactors[_iHighestK - i] << " " << vReadNameAndLength.second << " " << i << " " << vReadNameAndLength.second - i * 3 + 1 <<" " << arrWeightingFactors[_iHighestK - i] << endl;
				}

				for (auto it = vTempResultVec.begin(); it != vTempResultVec.end();) {
					double relativeScore = 0.0;
					if (_bTranslated) {
						relativeScore = double(get<1>(*it)) / (1.0 + log2(mFrequencies[get<0>(*it)] * double(vReadNameAndLength.second - _iHighestK + 1)));
					}
					else {
						relativeScore = double(get<1>(*it)) / (1.0 + log2(mFrequencies[get<0>(*it)] * double(vReadNameAndLength.second - _iHighestK * 3 + 1)));
					}
					if (relativeScore >= fThreshold) {
						get<2>(*it) = relativeScore;
						++it;
					}
					else {
						it = vTempResultVec.erase(it);
					}
				}

				if (vTempResultVec.size() == 0) {
					switch (format) {
					case OutputFormat::tsv:
						fOut << iReadNum << "\t" << vReadNameAndLength.first << "\t-\t-\t-" << "\n";
						break;

					case OutputFormat::Json:
						if (iReadNum == 0) {
							fOut << "{" << "\n";
						}
						else {
							fOut << "," << "\n" << "{" << "\n";
						}
						fOut << "\t\"Read number\": " << iReadNum << ",\n\t\"Specifier from input file\": \"" + vReadNameAndLength.first + "\",\n\t\"Top hits\": [\n\t],\n\t\"Further hits\": [\n\t]\n}";
						break;

					case OutputFormat::JsonL:
						fOut << "{ \"Read number\": " << iReadNum << "," << " \"Specifier from input file\": \"" + vReadNameAndLength.first + "\", " << "\"Top hits\": [], \"Further hits\": [] }\n";
						break;

					case OutputFormat::Kraken:
						fOut << "U\t" << vReadNameAndLength.first << "\t0\t" << vReadNameAndLength.second << "\tA:00\n";
						break;
					};
				}
				else {

# if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
					__gnu_parallel::sort(vTempResultVec.begin(), vTempResultVec.end(), [](const tuple<uint64_t, float, double>& p1, const tuple<uint64_t, float, double>& p2) { return get<2>(p1) > get<2>(p2); }, __gnu_parallel::balanced_quicksort_tag(_iNumOfThreads));
#else					
#if __has_include(<execution>)
					sort(std::execution::par_unseq, vTempResultVec.begin(), vTempResultVec.end(), [](const tuple<uint64_t, float, double>& p1, const tuple<uint64_t, float, double>& p2) { return get<2>(p1) > get<2>(p2); });
#else
					sort(vTempResultVec.begin(), vTempResultVec.end(), [](const tuple<uint64_t, float, double>& p1, const tuple<uint64_t, float, double>& p2) { return get<2>(p1) > get<2>(p2); });
#endif
#endif

					auto maxValue = get<1>(*(max_element(vTempResultVec.begin(), vTempResultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<1>(a) < get<1>(b); })));
					int32_t iTopHitCounter = 1;
					for (int16_t i = 1; i < static_cast<int16_t>(vTempResultVec.size()) && i < iNumOfBeasts; ++i) {
						if ((get<1>(vTempResultVec[i])) / (maxValue) > 0.8f) {
							++iTopHitCounter;
						}
						else {
							break;
						}
					}
					string sOut = "", sOut2 = "", sOut3 = "";
					sOut.reserve(1000);
					auto it = vTempResultVec.begin();
					float iValueBefore = 0;

					switch (format) {
					case OutputFormat::tsv:
						// tsv

						// += is faster than () + ()
						Utilities::itostr(iReadNum, sOut);
						sOut += "\t";
						sOut += vReadNameAndLength.first;
						sOut += "\t";

						for (int16_t j = 0; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += ";";
							sOut2 += mOrganisms[get<0>(*it)];
							sOut2 += ";";

							dtoa_milo(get<2>(*it), sOut3);
							sOut3 += ",";
							dtoa_milo(get<1>(*it), sOut3);
							sOut3 += ";";

							if (iValueBefore != get<1>(*it)) {
								iValueBefore = get<1>(*it);
								++j;
							}
						}
						if (sOut.back() == ';') {
							sOut.pop_back();
						}
						if (sOut2.back() == ';') {
							sOut2.pop_back();
						}
						if (sOut3.back() == ';') {
							sOut3.pop_back();
						}
						if (sOut2.length()) {
							fOut << sOut + "\t" + sOut2 + "\t" + sOut3 + "\t" << (bestScore - get<1>(vTempResultVec[0])) / bestScore << "\n";
						}
						break;

					case OutputFormat::Json:
						if (iReadNum == 0) {
							sOut += "{\n";
						}
						else {
							sOut += ",\n{\n";
						}

						sOut += "\t\"Read number\": ";
						Utilities::itostr(iReadNum, sOut);
						sOut += ",\n";
						sOut += "\t\"Specifier from input file\": \"";
						sOut += vReadNameAndLength.first;
						sOut += "\",\n\t\"Top hits\": [\n";

						for (int16_t i = 0; i < iTopHitCounter; ++i, ++it) {
							if (i == 0) {
								sOut += "\t{\n";
							}
							else {
								sOut += ",\n\t{\n";
							}

							sOut += "\t\t\"tax ID\": \"";
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += "\",\n";
							sOut += "\t\t\"Name\": \"";
							sOut += mOrganisms[get<0>(*it)];
							sOut += "\",\n";
							sOut += "\t\t\"k-mer Score\": ";
							dtoa_milo(get<1>(*it), sOut);
							sOut += ",\n";
							sOut += "\t\t\"Relative Score\": ";
							dtoa_milo(get<2>(*it), sOut);
							sOut += ",\n";
							sOut += "\t\t\"Error\": ";
							dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut);
							sOut += "\n\t}";
						}

						sOut += "\n\t],\n\t\"Further hits\": [\n";

						for (int32_t j = iTopHitCounter; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
							if (j == iTopHitCounter) {
								sOut += "\t{\n";
							}
							else {
								sOut += ",\n\t{\n";
							}

							sOut += "\t\t\"tax ID\": \"";
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += "\",\n";
							sOut += "\t\t\"Name\": \"";
							sOut += mOrganisms[get<0>(*it)];
							sOut += "\",\n";
							sOut += "\t\t\"k-mer Score\": ";
							dtoa_milo(get<1>(*it), sOut);
							sOut += ",\n";
							sOut += "\t\t\"Relative Score\": ";
							dtoa_milo(get<2>(*it), sOut);
							sOut += ",\n";
							sOut += "\t\t\"Error\": ";
							dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut);
							sOut += "\n\t}";

							if (iValueBefore != get<1>(*it)) {
								iValueBefore = get<1>(*it);
								++j;
							}
						}

						sOut += "\n\t]\n}";
						fOut << sOut;
						break;

					case OutputFormat::JsonL:
						// json lines
						sOut += "{ \"Read number\": ";
						Utilities::itostr(iReadNum, sOut);
						sOut += ",";
						sOut += " \"Specifier from input file\": \"";
						sOut += vReadNameAndLength.first;
						sOut += "\", \"Top hits\": [";

						for (int16_t i = 0; i < iTopHitCounter; ++i, ++it) {
							if (i == 0) {
								sOut += "{";
							}
							else {
								sOut += ",{";
							}

							sOut += " \"tax ID\": \"";
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += "\",";
							sOut += " \"Name\": \"";
							sOut += mOrganisms[get<0>(*it)];
							sOut += "\",";
							sOut += " \"k-mer Score\": ";
							dtoa_milo(get<1>(*it), sOut);
							sOut += ",";
							sOut += " \"Relative Score\": ";
							dtoa_milo(get<2>(*it), sOut);
							sOut += ",";
							sOut += " \"Error\": ";
							dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut);
							sOut += "}";
						}

						sOut += "], \"Further hits\": [";

						for (int32_t j = iTopHitCounter; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
							if (j == iTopHitCounter) {
								sOut += "{";
							}
							else {
								sOut += ", {";
							}

							sOut += " \"tax ID\": \"";
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += "\",";
							sOut += " \"Name\": \"";
							sOut += mOrganisms[get<0>(*it)];
							sOut += "\",";
							sOut += " \"k-mer Score\": ";
							dtoa_milo(get<1>(*it), sOut);
							sOut += ",";
							sOut += " \"Relative Score\": ";
							dtoa_milo(get<2>(*it), sOut);
							sOut += ",";
							sOut += " \"Error\": ";
							dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut);
							sOut += "}";

							if (iValueBefore != get<1>(*it)) {
								iValueBefore = get<1>(*it);
								++j;
							}
						}

						sOut += "] }\n";
						fOut << sOut;
						break;

					case OutputFormat::Kraken:
						sOut += "C\t";
						sOut += vReadNameAndLength.first;
						sOut += "\t";
						Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
						sOut += "\t";
						Utilities::itostr(vReadNameAndLength.second, sOut);
						sOut += "\t";

						for (int16_t i = 0; i < iTopHitCounter; ++i, ++it) {
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += ":";
							dtoa_milo(get<1>(*it), sOut);
							sOut += " ";
						}

						for (int32_t j = iTopHitCounter; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += ":";
							dtoa_milo(get<1>(*it), sOut);
							sOut += " ";
							if (iValueBefore != get<1>(*it)) {
								iValueBefore = get<1>(*it);
								++j;
							}
						}

						fOut << sOut << "\n";
						break;
					};
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl;
				throw;
			}
		}

	public:
		/////////////////////////////////////////////////////////////////////////////////
		void CompareWithLib_partialSort(const string& contentFile, const string& sLibFile, const string& fInFile, const string& fOutFile, const string& fTableFile, const uint8_t& iTrieDepth, const int64_t& iMemory, const bool& bSpaced, bool bRAM, const bool& bUnique, const uint8_t& iPrefixCheckMode, const float& fThreshold) {

			try {
				// test if files exists
				if (!ifstream(contentFile) || !ifstream(sLibFile) || !ifstream(sLibFile + "_f.txt") || !ifstream(sLibFile + "_trie.txt")) {
					throw runtime_error("One of the files does not exist");
				}
//#define TIME
#ifdef TIME
				auto startTIME = std::chrono::high_resolution_clock::now();
#endif

				// get names of idxes
				ifstream fContent(contentFile);
				uint32_t iAmountOfSpecies = 1;
				string sTempLine = "";
				vector<string> mOrganisms;
				unordered_map<uint32_t, uint32_t> mTaxToIdx;
				vector<uint32_t> mIdxToTax;
				mTaxToIdx[0] = 0;
				mOrganisms.push_back("non_unique");
				mIdxToTax.push_back(0);
				while (getline(fContent, sTempLine)) {
					if (sTempLine != "") {
						const auto& tempLineContent = Utilities::split(sTempLine, '\t');
						mOrganisms.push_back(Utilities::removeCharFromString(tempLineContent[0], ','));
						mIdxToTax.push_back(stoul(tempLineContent[1]));
						mTaxToIdx[stoul(tempLineContent[1])] = iAmountOfSpecies++;
					}
				}
				fContent.close();

				// get frequencies
				ifstream fFrequencies(sLibFile + "_f.txt");
				unique_ptr<uint64_t[]> mFrequencies(new uint64_t[iAmountOfSpecies]);
				uint32_t iCounterForFreqs = 0;
				while (getline(fFrequencies, sTempLine)) {
					if (sTempLine != "") {
						const auto& vLine = Utilities::split(sTempLine, '\t');
						mFrequencies[iCounterForFreqs++] = stoull(vLine[1]);
					}
				}


				// get Size of Lib and open it for reading
				ifstream fLibInfo(sLibFile + "_info.txt");
				uint64_t iSizeOfLib = 0;
				fLibInfo >> iSizeOfLib;
				bool bPartitioned = false;
				uint32_t iTempVecType = 0;
				fLibInfo >> iTempVecType;

				if (iTempVecType == 3) {
					bPartitioned = true;
				}
				fLibInfo.close();

				// Assert, that if the index was shrunken via trieHalf MinK can only be larger than 6
				if (bPartitioned && _iMinK <= 6) {
					cerr << "ERROR: k can only be larger than 6 if your index was shrunken via strategy 2. Setting it to 7..." << endl;
					_iMinK = 7;
					_iNumOfK = _iMaxK - _iMinK + 1;
				}
				if (bUnfunny) {
					_iMinK = 6;
					_iMaxK = 6;
					_iNumOfK = _iMaxK - _iMinK + 1;
				}

				// Create threadpool(s), in stxxl mode we can only have synced parallelism (as in no two threads should not access the same vector instance)
				vector<WorkerThread> workerThreadPool(_iNumOfThreads);
				for (int32_t i = 0; i < _iNumOfThreads; ++i) {
					workerThreadPool[i].setID(i);
				}

				// load index
				unique_ptr<stxxlFile> stxxlLibFile(new stxxlFile(sLibFile, stxxl::file::RDONLY));
				unique_ptr<unique_ptr<const contentVecType_32p>[]> vLib;
				unique_ptr<unique_ptr<const index_t_p>[]> vLibParted_p;
				vector<packedBigPair> vLib_RAM_Full;
				vector<packedPair> vLib_RAM_Half;
				vector<uint16_t> vLib_RAM_taxaOnly;

				unique_ptr<unique_ptr<const taxaOnly>[]> vLib_taxaOnly;

				int64_t iBytesUsedBySTXXLVectors = 0;
				if (bPartitioned) {
					vLibParted_p.reset(new unique_ptr<const index_t_p>[_iNumOfThreads]);
					//stxxlLibFile.reset(new unique_ptr<stxxlFile>[_iNumOfThreads]);
					for (int32_t i = 0; i < _iNumOfThreads; ++i) {
						//stxxlLibFile[i].reset(new stxxlFile(sLibFile + "_"+ to_string(i), stxxl::file::RDONLY));
						vLibParted_p[i].reset(new const index_t_p(stxxlLibFile.get(), iSizeOfLib));
					}
					iBytesUsedBySTXXLVectors = _iNumOfThreads * index_t_p::block_size * index_t_p::page_size * (vLibParted_p[0])->numpages();
				}
				else {
					if (bUnfunny) {
						vLib_taxaOnly.reset(new unique_ptr<const taxaOnly>[_iNumOfThreads]);
						for (int32_t i = 0; i < _iNumOfThreads; ++i) {
							vLib_taxaOnly[i].reset(new const taxaOnly(stxxlLibFile.get(), iSizeOfLib));
						}
						iBytesUsedBySTXXLVectors = _iNumOfThreads * taxaOnly::block_size * taxaOnly::page_size * (vLib_taxaOnly[0])->numpages();
					}
					else {
						vLib.reset(new unique_ptr<const contentVecType_32p>[_iNumOfThreads]);
						for (int32_t i = 0; i < _iNumOfThreads; ++i) {
							vLib[i].reset(new const contentVecType_32p(stxxlLibFile.get(), iSizeOfLib));
						}
						iBytesUsedBySTXXLVectors = _iNumOfThreads * contentVecType_32p::block_size * contentVecType_32p::page_size * (vLib[0])->numpages();
					}
				}

				// load Trie
				const uint8_t& iTD = iTrieDepth;
				Trie T(static_cast<int8_t>(_iMaxK), static_cast<int8_t>(_iMinK), iTD, _iNumOfThreads, iPrefixCheckMode);
				T.LoadFromStxxlVec(sLibFile);
				T.SetForIsInTrie((_iMinK < 6) ? static_cast<uint8_t>(_iMinK) : static_cast<uint8_t>(6));

				if (_bVerbose) {
					T.GetIfVecIsUsed();
				}

#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "Trie load " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				// This holds the hits for each organism
				// Initialize this here to calculate the average memory usage
				const uint64_t& iMult = _iNumOfThreads * _iNumOfK * uint64_t(iAmountOfSpecies);
				unique_ptr<double[]> vCount_all(new double[iMult]);
				unique_ptr<uint64_t[]> vCount_unique(new uint64_t[iMult]);

				for (uint64_t i = 0; i < iMult; ++i) {
					vCount_all[i] = 0.;
					vCount_unique[i] = 0;
				}

				// If the user wishes to load the index into RAM, check if there is enough space.
				// First calculate the average memory usage: 
				//	Size of trie + size of counts_all/_unique + index size + memory for tax ids + chunksize from input + overhead buffer for various allocations 
				int64_t iMemoryUsageOnAverage = T.GetSize() + iMult * sizeof(uint64_t) + iBytesUsedBySTXXLVectors + _iNumOfThreads * Utilities::sBitArray(iAmountOfSpecies).sizeInBytes() + 14399756 + 4 * iAmountOfSpecies + 1024ull * 1024ull * 1024ull;
				
				// Set memory boundaries
				int64_t iSoftMaxMemoryUsage = 0;
				if (iMemory > iMemoryUsageOnAverage) {
					iSoftMaxMemoryUsage = iMemory - iMemoryUsageOnAverage;
				}
				else {
					cerr << "ERROR: Not enough memory given, try to download more RAM. Adding 1GB. May lead to bad_alloc errors..." << endl;
					//iSoftMaxMemoryUsage = 1024ull * 1024ull * 1ull / (sizeof(tuple<uint64_t, uint64_t>)*_iNumOfThreads);
					iSoftMaxMemoryUsage = static_cast<int64_t>(1024ull * 1024ull * 1024ull); // 1024ull * 1024ull * 1024ull
				}

				if (bRAM) {
					try {
						if (bPartitioned || (_iMinK > 6 && iAmountOfSpecies <= 65535 && !bUnfunny)) {
							if (iSoftMaxMemoryUsage + iBytesUsedBySTXXLVectors - static_cast<int64_t>(iSizeOfLib * sizeof(packedPair)) < 0) {
								cerr << "ERROR: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
								bRAM = false;
							}
							else {
								vLib_RAM_Half.reserve(iSizeOfLib);

								if (bPartitioned) {
									stxxl::vector_bufreader<index_t_p::const_iterator> bufferedReader(vLibParted_p[0]->cbegin(), vLibParted_p[0]->cend(), 0);
									for (; !bufferedReader.empty(); ++bufferedReader) {
										vLib_RAM_Half.push_back(*bufferedReader);
									}

									for (int32_t i = 0; i < _iNumOfThreads; ++i) {
										vLibParted_p[i].reset();
									}
									vLibParted_p.reset();
								}
								else {
									stxxl::vector_bufreader<contentVecType_32p::const_iterator> bufferedReader(vLib[0]->cbegin(), vLib[0]->cend(), 0);
									for (; !bufferedReader.empty(); ++bufferedReader) {
										vLib_RAM_Half.push_back(packedPair(static_cast<uint32_t>(bufferedReader->first & 1073741823ULL), static_cast<uint16_t>(mTaxToIdx[bufferedReader->second])));
									}

									for (int32_t i = 0; i < _iNumOfThreads; ++i) {
										vLib[i].reset();
									}
									vLib.reset();
									bPartitioned = true;
								}

								iSoftMaxMemoryUsage = iSoftMaxMemoryUsage + iBytesUsedBySTXXLVectors - iSizeOfLib * sizeof(packedPair);
							}
						}
						else {
							if (bUnfunny) {
								if (iSoftMaxMemoryUsage + iBytesUsedBySTXXLVectors - static_cast<int64_t>(iSizeOfLib * sizeof(uint16_t)) < 0) {
									cerr << "ERROR: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
									bRAM = false;
								}
								else {
									vLib_RAM_taxaOnly.reserve(iSizeOfLib);
									stxxl::vector_bufreader<taxaOnly::const_iterator> bufferedReader(vLib_taxaOnly[0]->cbegin(), vLib_taxaOnly[0]->cend(), 0);
									for (; !bufferedReader.empty(); ++bufferedReader) {
										vLib_RAM_taxaOnly.push_back(*bufferedReader);
									}

									for (int32_t i = 0; i < _iNumOfThreads; ++i) {
										vLib_taxaOnly[i].reset();
									}
									vLib_taxaOnly.reset();
									iSoftMaxMemoryUsage = iSoftMaxMemoryUsage + iBytesUsedBySTXXLVectors - iSizeOfLib * sizeof(uint16_t);
								}
							}
							else {
								if (iSoftMaxMemoryUsage + iBytesUsedBySTXXLVectors - static_cast<int64_t>(iSizeOfLib * sizeof(packedBigPair)) < 0) {
									cerr << "ERROR: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
									bRAM = false;
								}
								else {
									vLib_RAM_Full.reserve(iSizeOfLib);
									stxxl::vector_bufreader<contentVecType_32p::const_iterator> bufferedReader(vLib[0]->cbegin(), vLib[0]->cend(), 0);
									for (; !bufferedReader.empty(); ++bufferedReader) {
										vLib_RAM_Full.push_back(*bufferedReader);
									}

									for (int32_t i = 0; i < _iNumOfThreads; ++i) {
										vLib[i].reset();
									}
									vLib.reset();
									iSoftMaxMemoryUsage = iSoftMaxMemoryUsage + iBytesUsedBySTXXLVectors - iSizeOfLib * sizeof(packedBigPair);
								}
							}
						}
					}
					catch (const bad_alloc&) {
						cerr << "ERROR: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
						bRAM = false;
						vLib_RAM_Half.clear();
						vLib_RAM_taxaOnly.clear();
						vLib_RAM_Full.clear();
					}
				}

#ifdef TIME
				auto endTIME = std::chrono::high_resolution_clock::now();
				cout << "Initialization " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif



				uint64_t iTimeFastq = 0, iTimeCompare = 0, iNumOfReads = 0, iNumOfReadsOld = 0, iNumOfReadsSum = 0;
				vector<tuple<uint64_t, uint64_t, uint32_t, uint32_t>> vInputVec;
				/*while (true) {
					try {
						vInputVec.reserve(iSoftMaxMemoryUsage / sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>));
						//vInputVec.resize(iSoftMaxMemoryUsage / sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>));
						break;
					}
					catch (const bad_alloc&) {
						//vInputVec.clear();
						vInputVec.shrink_to_fit();
						if (iSoftMaxMemoryUsage > 1073741824ull) { // 1073741824 == 1024^3 == 1GB
							iSoftMaxMemoryUsage -= 1073741824ull;
						}
						else {
							cerr << "ERROR: Not enough memory available. Please try again with a lower number after - m" << endl;
							throw;
						}
					}
				}*/

				const bool& bReadIDsAreInteresting = fOutFile != "";

				unique_ptr<float[]> vReadIDtoTaxID; // Array of #species times #reads with scores as values
				auto getVecOfScored = [&iAmountOfSpecies](const unique_ptr<float[]>& vReadIDtoGenID, const uint64_t& iIdx) {
					try {
						vector<tuple<uint32_t, float, double>> output;
						for (uint32_t i = 1; i < iAmountOfSpecies; ++i) {
							if (vReadIDtoGenID[iIdx*iAmountOfSpecies + i] > 0.f) {
								output.push_back(make_tuple(i, vReadIDtoGenID[iIdx*iAmountOfSpecies + i], 0.0));
							}
						}
						return output;
					}
					catch (...) {
						cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
					}
				};

				list<pair<string, uint32_t>> vReadNameAndLength;

				uint64_t iNumberOfkMersInInput = 0;
				vector<uint64_t> vNumberOfGarbagekMersPerK(_iNumOfK, 0);

				auto pathAndSize = Utilities::gatherFilesFromPath(fInFile);
				const auto& vInputFiles = pathAndSize.first;
				size_t overallFileSize = pathAndSize.second, allFilesProgress = 0, charsReadOverall = 0;

#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "Reserve memory " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				// allow multiple input files
				for (const auto& inFile : vInputFiles) {

					if (_bVerbose) {
						cout << "OUT: Current file: " << inFile.first << endl;
					}


					string fileName = "";
					if (vInputFiles.size() > 1) { // get file name without path and ending
						const auto& vRawNameSplit = Utilities::split(inFile.first.substr(fInFile.size(), inFile.first.size()), '.');
						if (vRawNameSplit.size() == 1) {
							fileName = vRawNameSplit[0];
						}
						else {
							for (size_t rawC = 0; rawC < vRawNameSplit.size() - 1; ++rawC) {
								fileName += vRawNameSplit[rawC] + ".";
							}
							fileName.pop_back();
						}
					}

					// see if input is gziped or not and if it's a fasta or fastq file
					bool isGzipped = inFile.second; //= (inFile[inFile.length() - 3] == '.' && inFile[inFile.length() - 2] == 'g' && inFile[inFile.length() - 1] == 'z');

					bool bIsGood = false, bIsFasta = false;

					unique_ptr<ifstream> fast_q_a_File;
					unique_ptr<igzstream> fast_q_a_File_gz;
					uint64_t iFileLength = 0;
					vector<char> inFilebuffer(16777216);

					if (isGzipped) {
						if (_bVerbose) {
							cout << "OUT: File is gzipped, no progress output can be shown." << endl;
						}

						fast_q_a_File_gz.reset(new igzstream());
						fast_q_a_File_gz->rdbuf()->pubsetbuf(&inFilebuffer[0], 16777216);
						fast_q_a_File_gz->open(inFile.first.c_str());
						bIsGood = fast_q_a_File_gz->good();
						char cLessOrAt = static_cast<char>(fast_q_a_File_gz->peek());
						if (cLessOrAt == '>') {
							bIsFasta = true;
						}
						else {
							if (cLessOrAt == '@') {
								bIsFasta = false;
							}
							else {
								throw runtime_error("Input does not start with @ or >.");
							}
						}

						fast_q_a_File_gz->exceptions(std::ios_base::badbit);
					}
					else {
						fast_q_a_File.reset(new ifstream());
						fast_q_a_File->rdbuf()->pubsetbuf(&inFilebuffer[0], 16777216);
						fast_q_a_File->open(inFile.first);
						bIsGood = fast_q_a_File->good();

						char cLessOrAt = static_cast<char>(fast_q_a_File->peek());
						if (cLessOrAt == '>') {
							bIsFasta = true;
						}
						else {
							if (cLessOrAt == '@') {
								bIsFasta = false;
							}
							else {
								throw runtime_error("Input does not start with @ or >.");
							}
						}

						fast_q_a_File->seekg(0, fast_q_a_File->end);
						iFileLength = fast_q_a_File->tellg();
						fast_q_a_File->seekg(0, fast_q_a_File->beg);
					}

					ofstream fOut;
					// a larger buffer works better for SSDs or HPCCs
					vector<char> outFilebuffer(16777216);
					fOut.rdbuf()->pubsetbuf(&outFilebuffer[0], 16777216);
					//fOut.exceptions(std::ifstream::failbit | std::ifstream::badbit);
					if (bReadIDsAreInteresting) {
						fOut.open((vInputFiles.size() > 1) ? fOutFile + fileName + outputFormatFileEnding() : fOutFile); // in case of multiple input files, specify only beginning of the output and the rest will be appended
						if (fOut) {
							if (format == OutputFormat::tsv) {
								fOut << "#Read number\tSpecifier from input file\tMatched taxa\tNames\tScores{relative,k-mer}\tError" << "\n";
							}
							else {
								if (format == OutputFormat::Json) {
									fOut << "[" << "\n";
								}
							}
						}
						else {
							throw runtime_error("Readwise output file could not be created!");
						}
					}
					// Test if profile file can be written (it would suck to see your calculation gone to waste because it is only tested at the end now wouldn't it?!)
					if (fTableFile != "") {
						ofstream tableFileStream((vInputFiles.size() > 1) ? fTableFile + fileName + ".csv" : fTableFile);
						if (!tableFileStream) {
							throw runtime_error("Profile file couldn't be opened for writing!");
						}
					}

#ifdef TIME
					endTIME = std::chrono::high_resolution_clock::now();
					cout << "Prepare input " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
					startTIME = std::chrono::high_resolution_clock::now();
#endif

					unique_ptr<strTransfer> transferBetweenRuns(new strTransfer);
					transferBetweenRuns->iCurrentOverallPercentage = allFilesProgress;
					transferBetweenRuns->iNumOfAllCharsRead = charsReadOverall;
					vector<tuple<readIDType, float, double>> vSavedScores;
					//readIDType iReadIDofSavedScores = 0;

					// read input
					while (bIsGood) {
						ios::sync_with_stdio(false);
						std::cin.tie(nullptr);
						auto start = std::chrono::high_resolution_clock::now();
						if (bIsFasta) {

							//iSoftMaxMemoryUsage = 19199040;

							if (isGzipped) {
								iNumberOfkMersInInput += readFasta_partialSort(*fast_q_a_File_gz, vInputVec, vReadNameAndLength, iSoftMaxMemoryUsage, iAmountOfSpecies, iFileLength, iFileLength, bReadIDsAreInteresting, transferBetweenRuns, workerThreadPool, vNumberOfGarbagekMersPerK);
							}
							else {
								iNumberOfkMersInInput += readFasta_partialSort(*fast_q_a_File, vInputVec, vReadNameAndLength, iSoftMaxMemoryUsage, iAmountOfSpecies, iFileLength, overallFileSize, bReadIDsAreInteresting, transferBetweenRuns, workerThreadPool, vNumberOfGarbagekMersPerK);
							}

							//iNumOfReads = (transferBetweenRuns->iNumOfNewReads > 0) ? transferBetweenRuns->iNumOfNewReads : 1; // vReadIDs.size();
							iNumOfReads = transferBetweenRuns->iNumOfNewReads;
						}
						else {

							//iSoftMaxMemoryUsage = 19199040;

							if (isGzipped) {
								iNumberOfkMersInInput += readFastq_partialSort(*fast_q_a_File_gz, vInputVec, vReadNameAndLength, iSoftMaxMemoryUsage, iAmountOfSpecies, iFileLength, iFileLength, bReadIDsAreInteresting, transferBetweenRuns, workerThreadPool, vNumberOfGarbagekMersPerK);
							}
							else {
								iNumberOfkMersInInput += readFastq_partialSort(*fast_q_a_File, vInputVec, vReadNameAndLength, iSoftMaxMemoryUsage, iAmountOfSpecies, iFileLength, overallFileSize, bReadIDsAreInteresting, transferBetweenRuns, workerThreadPool, vNumberOfGarbagekMersPerK);
							}

							//iNumOfReads = (transferBetweenRuns->iNumOfNewReads > 0) ? transferBetweenRuns->iNumOfNewReads : 1; // vReadIDs.size();
							iNumOfReads = transferBetweenRuns->iNumOfNewReads;
						}

#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "Input " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						// sort prefixes for each range in parallel
# if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
						// Some context:
						// The gnu parallel quicksort implementation only uses two cores whereas the gnu parallel merge uses all but isn't in-place. This significantly worsened performance in Linux environments.
						// Futhermore, the gcc 9.* compiler uses Threadblocks for its stl conform parallel implementation as of now (2020). This is unacceptable for kASA.
						// Therefore I set out to find a better one and found the preliminary implementation for C++17 inside Visual Studio which was published under the Apache Software License 2.0.
						// This is compatible with kASA and as near to the current implementation inside Visual Studio 2017 as can be without risking copyright infringement towards Microsoft.
						// It can be found inside the Utilities.hpp header.
						Utilities::parallelQuicksort(vInputVec.begin(), vInputVec.end(), [](const tuple<uint64_t, uint64_t, uint32_t, uint32_t>& p1, const tuple<uint64_t, uint64_t, uint32_t, uint32_t>& p2) {
							return get<1>(p1) < get<1>(p2);
						}, _iNumOfThreads);

						if (bPartitioned) {
							__gnu_parallel::for_each(vInputVec.begin(), vInputVec.end(), [this, &T](tuple<uint64_t, uint64_t, uint32_t, uint32_t>& a) {
								uint64_t start = 0;
								uint32_t range = 0;
								T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
								if (start != numeric_limits<uint64_t>::max()) {
									get<0>(a) = start;
									get<1>(a) = get<1>(a) & 1073741823ULL;
									get<2>(a) = range;
								}
								else {
									get<0>(a) = start;
								}
							});
						}
						else {
							__gnu_parallel::for_each(vInputVec.begin(), vInputVec.end(), [this, &T](tuple<uint64_t, uint64_t, uint32_t, uint32_t>& a) {
								uint64_t start = 0;
								uint32_t range = 0;
								T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
								if (start != numeric_limits<uint64_t>::max()) {
									get<0>(a) = start;
									get<2>(a) = range;
								}
								else {
									get<0>(a) = start;
								}
							});
						}
#else					
#if __has_include(<execution>)

						sort(std::execution::par_unseq, vInputVec.begin(), vInputVec.end(), [](const tuple<uint64_t, uint64_t, uint32_t, uint32_t>& p1, const tuple<uint64_t, uint64_t, uint32_t, uint32_t>& p2) {
							return get<1>(p1) < get<1>(p2);
						});

						if (bPartitioned) {
							for_each(std::execution::par_unseq, vInputVec.begin(), vInputVec.end(), [this, &T](tuple<uint64_t, uint64_t, uint32_t, uint32_t>& a) {
								uint64_t start = 0;
								uint32_t range = 0;
								T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
								if (start != numeric_limits<uint64_t>::max()) {
									get<0>(a) = start;
									get<1>(a) = get<1>(a) & 1073741823ULL;
									get<2>(a) = range;
								}
								else {
									get<0>(a) = start;
								}
							});
						}
						else {
							for_each(std::execution::par_unseq, vInputVec.begin(), vInputVec.end(), [this, &T](tuple<uint64_t, uint64_t, uint32_t, uint32_t>& a) {
								uint64_t start = 0;
								uint32_t range = 0;
								T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
								if (start != numeric_limits<uint64_t>::max()) {
									get<0>(a) = start;
									get<2>(a) = range;
								}
								else {
									get<0>(a) = start;
								}
							});
						}
#else
						Utilities::parallelQuicksort(vInputVec.begin(), vInputVec.end(), [](const tuple<uint64_t, uint64_t, uint32_t, uint32_t>& p1, const tuple<uint64_t, uint64_t, uint32_t, uint32_t>& p2) {
							return get<1>(p1) < get<1>(p2);
						}, _iNumOfThreads);

						if (bPartitioned) {
							for_each(vInputVec.begin(), vInputVec.end(), [this, &T](tuple<uint64_t, uint64_t, uint32_t, uint32_t>& a) {
								uint64_t start = 0;
								uint32_t range = 0;
								T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
								if (start != numeric_limits<uint64_t>::max()) {
									get<0>(a) = start;
									get<1>(a) = get<1>(a) & 1073741823ULL;
									get<2>(a) = range;
								}
								else {
									get<0>(a) = start;
								}
							});
						}
						else {
							for_each(vInputVec.begin(), vInputVec.end(), [this, &T](tuple<uint64_t, uint64_t, uint32_t, uint32_t>& a) {
								uint64_t start = 0;
								uint32_t range = 0;
								T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((_iMinK > 6) ? 6 : _iMinK), move(start), move(range));
								if (start != numeric_limits<uint64_t>::max()) {
									get<0>(a) = start;
									get<2>(a) = range;
								}
								else {
									get<0>(a) = start;
								}
							});
						}
#endif
#endif


#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "Elimination " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						if (bUnique) {
							auto newEnd = unique(vInputVec.begin(), vInputVec.end());
							vInputVec.resize(newEnd - vInputVec.begin());
						}

						auto end = std::chrono::high_resolution_clock::now();
						iTimeFastq += chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();


						if (bReadIDsAreInteresting) {
							// This holds the mapping read ID -> Genus ID
							if (iNumOfReads <= iNumOfReadsOld) {
								memset(vReadIDtoTaxID.get(), 0, iNumOfReads*iAmountOfSpecies * sizeof(float));
								/*for (uint64_t i = 0; i < iNumOfReads*iAmountOfSpecies; ++i) {
									vReadIDtoTaxID[i] = 0.f;
								}*/
							}
							else {
								vReadIDtoTaxID.reset();
								vReadIDtoTaxID.reset(new float[iNumOfReads*iAmountOfSpecies]);

								memset(vReadIDtoTaxID.get(), 0, iNumOfReads*iAmountOfSpecies * sizeof(float));
								/*for (uint64_t i = 0; i < iNumOfReads*iAmountOfSpecies; ++i) {
									vReadIDtoTaxID[i] = 0.f;
								}*/
							}
						}
						else {
							iNumOfReadsSum += transferBetweenRuns->iNumOfNewReads;
						}

						function<void(const int32_t&, const uint64_t&, const uint64_t&)> foo;

						// now compare with index
						start = std::chrono::high_resolution_clock::now();

						if (bReadIDsAreInteresting) {
							if (bRAM) {
								if (bPartitioned) {
									foo = bind(&Compare::compareWithDatabase< vector<packedPair>*>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, &vLib_RAM_Half, ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(mTaxToIdx));//, ref(transferBetweenRuns->mReadIDToArrayIdx));
								}
								else {
									if (bUnfunny) {
										//foo = bind(&Compare::compareWithDatabase_sloppy< vector<uint16_t>*>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, &vLib_RAM_taxaOnly, ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(transferBetweenRuns->mReadIDToArrayIdx));
									}
									else {
										foo = bind(&Compare::compareWithDatabase< vector<packedBigPair>*>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, &vLib_RAM_Full, ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(mTaxToIdx));//, ref(transferBetweenRuns->mReadIDToArrayIdx));
									}
								}
							}
							else {
								if (bPartitioned) {
									foo = bind(&Compare::compareWithDatabase< unique_ptr<unique_ptr<const index_t_p>[]>>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(vLibParted_p), ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(mTaxToIdx));//, ref(transferBetweenRuns->mReadIDToArrayIdx));
								}
								else {
									if (bUnfunny) {
										//foo = bind(&Compare::compareWithDatabase_sloppy<unique_ptr<unique_ptr<const taxaOnly>[]>>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(vLib_taxaOnly), ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(transferBetweenRuns->mReadIDToArrayIdx));
									}
									else {
										foo = bind(&Compare::compareWithDatabase<unique_ptr<unique_ptr<const contentVecType_32p>[]>>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(vLib), ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(mTaxToIdx));//, ref(transferBetweenRuns->mReadIDToArrayIdx));
									}
								}
							}
						}
						else {
							if (bRAM) {
								if (bPartitioned) {
									foo = bind(&Compare::createProfile<vector<packedPair>*>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, &vLib_RAM_Half, ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies), ref(mTaxToIdx));
								}
								else {
									if (bUnfunny) {
										//foo = bind(&Compare::createProfile_sloppy<vector<uint16_t>*>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, &vLib_RAM_taxaOnly, ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies));
									}
									else {
										foo = bind(&Compare::createProfile< vector<packedBigPair>* >, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, &vLib_RAM_Full, ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies), ref(mTaxToIdx));
									}
								}
							}
							else {
								if (bPartitioned) {
									foo = bind(&Compare::createProfile< unique_ptr<unique_ptr<const index_t_p>[]> >, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(vLibParted_p), ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies), ref(mTaxToIdx));
								}
								else {
									if (bUnfunny) {
										//foo = bind(&Compare::createProfile_sloppy< unique_ptr<unique_ptr<const taxaOnly>[]> >, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(vLib_taxaOnly), ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies));
									}
									else {
										foo = bind(&Compare::createProfile< unique_ptr<unique_ptr<const contentVecType_32p>[]> >, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(vLib), ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies), ref(mTaxToIdx));
									}
								}

							}
						}

						// because the stxxl is not threadsafe (as in two threads cannot access the same location on the drive), 
						//   we need to separate the work and make sure that ranges are disjoint for every thread
						size_t iChunkSize = vInputVec.size() / _iNumOfThreads;
						size_t iStart = 0, iEnd = iChunkSize + vInputVec.size() % _iNumOfThreads;
						for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
							uint64_t iSeenRange = get<0>(vInputVec[iEnd - 1]);
							while (iEnd < vInputVec.size()) {
								if (get<0>(vInputVec[iEnd]) == iSeenRange) {
									iEnd++;
								}
								else {
									break;
								}
							}
							auto task = bind(foo, placeholders::_1, iStart, iEnd);
							workerThreadPool[iThreadID].pushTask(task);
							iStart = iEnd;
							iEnd += iChunkSize;
							if (iEnd > vInputVec.size()) {
								iEnd = vInputVec.size();
							}
						}

#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "Divide for Compare " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
							workerThreadPool[iThreadID].startThread();
						}
						for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
							workerThreadPool[iThreadID].waitUntilFinished();
						}

						if (someThingWentWrong) {
							rethrow_exception(someThingWentWrong);
						}

#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "Compare " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						end = std::chrono::high_resolution_clock::now();
						iTimeCompare += chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

						////////////////////////////////////////////////////////////////////////////////////////////////
						// save results
						if (bReadIDsAreInteresting) {

							uint64_t iReadIDStart = 0;
							// check if there is still some unfinished read which is now complete
							if (vSavedScores.size() && transferBetweenRuns->finished) {
								//if (transferBetweenRuns->lastLine.second != iReadIDofSavedScores) {
									//auto lastScoreVec = getVecOfScored(vReadIDtoTaxID, Utilities::checkIfInMap(transferBetweenRuns->mReadIDToArrayIdx, iReadIDofSavedScores)->second);
								auto lastScoreVec = getVecOfScored(vReadIDtoTaxID, 0);
								vSavedScores.insert(vSavedScores.end(), lastScoreVec.cbegin(), lastScoreVec.cend());
								lastScoreVec.clear();
								sort(vSavedScores.begin(), vSavedScores.end(), [](const tuple<uint32_t, float, double>& p1, const tuple<uint32_t, float, double>& p2) { return get<0>(p1) < get<0>(p2); });
								auto seen = vSavedScores[0];
								for (auto it = vSavedScores.begin() + 1; it != vSavedScores.end(); ++it) {
									if (get<0>(*it) != get<0>(seen)) {
										lastScoreVec.push_back(seen);
										seen = *it;
									}
									else {
										get<1>(seen) += get<1>(*it);
									}
								}
								lastScoreVec.push_back(seen);
								vSavedScores.swap(lastScoreVec);
								scoringFunc(move(vSavedScores), (iReadIDStart++) + iNumOfReadsSum, vReadNameAndLength.front(), mFrequencies, mIdxToTax, mOrganisms, fThreshold, move(fOut));
								vSavedScores.clear();

								//transferBetweenRuns->vReadIDs.erase(find(transferBetweenRuns->vReadIDs.begin(), transferBetweenRuns->vReadIDs.end(), iReadIDofSavedScores));
								//transferBetweenRuns->mReadIDToArrayIdx.erase(iReadIDofSavedScores);
								vReadNameAndLength.pop_front();
								//}
							}

							if (transferBetweenRuns->addTail) {
								// last read is not yet finished

								// save the score of the not yet finished
								//auto resultOfUnfinished = getVecOfScored(vReadIDtoTaxID, Utilities::checkIfInMap(transferBetweenRuns->mReadIDToArrayIdx, transferBetweenRuns->lastLine.second)->second);
								auto resultOfUnfinished = getVecOfScored(vReadIDtoTaxID, iNumOfReads - 1);
								if (resultOfUnfinished.size()) {
									vSavedScores.insert(vSavedScores.end(), resultOfUnfinished.cbegin(), resultOfUnfinished.cend());
									resultOfUnfinished.clear();
									sort(vSavedScores.begin(), vSavedScores.end(), [](const tuple<uint32_t, float, double>& p1, const tuple<uint32_t, float, double>& p2) { return get<0>(p1) < get<0>(p2); });
									auto seen = vSavedScores[0];
									for (auto it = vSavedScores.begin() + 1; it != vSavedScores.end(); ++it) {
										if (get<0>(*it) != get<0>(seen)) {
											resultOfUnfinished.push_back(seen);
											seen = *it;
										}
										else {
											get<1>(seen) += get<1>(*it);
										}
									}
									resultOfUnfinished.push_back(seen);
									vSavedScores.swap(resultOfUnfinished);
								}

								// save the finished ones
								scoringFunc(vReadIDtoTaxID, iReadIDStart, iNumOfReads - 1, iNumOfReadsSum, vReadNameAndLength, mFrequencies, mIdxToTax, mOrganisms, fThreshold, move(fOut));

								//for (; i < iNumOfReads - 1; ++i) {
									//auto tempID = transferBetweenRuns->vReadIDs.front();
									//transferBetweenRuns->vReadIDs.pop_front();
									//transferBetweenRuns->mReadIDToArrayIdx.erase(tempID);
									//vReadNameAndLength.pop_front();
								//}


								//iReadIDofSavedScores = transferBetweenRuns->iCurrentReadID;
								//transferBetweenRuns->mReadIDToArrayIdx[iReadIDofSavedScores] = 0;

								iNumOfReadsSum += iNumOfReads - 1;
								iNumOfReadsOld = (iNumOfReads - 1 < iNumOfReadsOld) ? iNumOfReadsOld : iNumOfReads - 1;
							}
							else {
								// reads finished

								// write down the saved one
								//uint64_t IReadIDStart = 0;
								//if (vSavedScores.size()) {
								//	auto resultOfFinished = getVecOfScored(vReadIDtoTaxID, 0);
								//	vSavedScores.insert(vSavedScores.end(), resultOfFinished.cbegin(), resultOfFinished.cend());
								//	resultOfFinished.clear();
								//	sort(vSavedScores.begin(), vSavedScores.end(), [](const tuple<uint32_t, float, double>& p1, const tuple<uint32_t, float, double>& p2) { return get<0>(p1) < get<0>(p2); });
								//	auto seen = vSavedScores[0];
								//	for (auto it = vSavedScores.begin() + 1; it != vSavedScores.end(); ++it) {
								//		if (get<0>(*it) != get<0>(seen)) {
								//			resultOfFinished.push_back(seen);
								//			seen = *it;
								//		}
								//		else {
								//			get<1>(seen) += get<1>(*it);
								//		}
								//	}
								//	resultOfFinished.push_back(seen);
								//	vSavedScores.swap(resultOfFinished);

								//	scoringFunc(move(vSavedScores), iNumOfReadsSum, vReadNameAndLength.front(), mFrequencies, mIdxToTax, mOrganisms, move(fOut));
								//	IReadIDStart = 1;
								//	//transferBetweenRuns->vReadIDs.erase(find(transferBetweenRuns->vReadIDs.begin(), transferBetweenRuns->vReadIDs.end(), iReadIDofSavedScores));
								//	//transferBetweenRuns->mReadIDToArrayIdx.erase(iReadIDofSavedScores);
								//	vReadNameAndLength.pop_front();
								//	vSavedScores.clear();
								//}

								// and now the regular ones
								scoringFunc(vReadIDtoTaxID, iReadIDStart, iNumOfReads, iNumOfReadsSum, vReadNameAndLength, mFrequencies, mIdxToTax, mOrganisms, fThreshold, move(fOut));
								//for (; i < iNumOfReads; ++i) {
								//	auto tempID = transferBetweenRuns->vReadIDs.front();
								//	transferBetweenRuns->vReadIDs.pop_front();
								//	transferBetweenRuns->mReadIDToArrayIdx.erase(tempID);
									//vReadNameAndLength.pop_front();
								//}

								iNumOfReadsSum += iNumOfReads;
								iNumOfReadsOld = iNumOfReads;
								vReadNameAndLength.clear();
							}
						}


						vInputVec.clear();

						if (isGzipped) {
							bIsGood = fast_q_a_File_gz->good();
						}
						else {
							bIsGood = fast_q_a_File->good();
						}
						charsReadOverall = transferBetweenRuns->iNumOfAllCharsRead;

#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "Saving results " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						//iterate until no dna is left
					}

					///////////////////////////////////////////////////////////////////////////////////////////////

					// if json is the output format for readToTaxa, end it with a ]
					if (bReadIDsAreInteresting && format == OutputFormat::Json) {
						fOut << "\n" << "]";
					}
					fOut.flush(); // empty the buffer to avoid memory leak

					// sum up parallel results
					for (int32_t iThreadID = 1; iThreadID < _iNumOfThreads; ++iThreadID) {
						const uint64_t& iStepsize = iThreadID * _iNumOfK * iAmountOfSpecies;
						for (uint64_t iIdx = 0; iIdx < _iNumOfK*iAmountOfSpecies; ++iIdx) {
							vCount_all[iIdx] += vCount_all[iStepsize + iIdx];
							vCount_unique[iIdx] += vCount_unique[iStepsize + iIdx];
						}
					}
#ifdef DEBUGOUT
					for (int i = 0; i < 4; ++i) {
						cout << amountArr[i] << endl;
					}
					cout << "tax:" << endl;
					for (int i = 0; i < 5; ++i) {
						cout << amountTaxArr[i] << endl;
					}
#endif
					// get profiling results
					vector<uint64_t> vSumOfUniquekMers(_iNumOfK);
					vector<double> vSumOfNonUniques(_iNumOfK);
					vector<tuple<string, vector<pair<double, uint64_t>>, uint32_t>> vOut(iAmountOfSpecies, tuple<string, vector<pair<double, uint64_t>>, uint32_t>("", vector<pair<double, uint64_t>>(_iNumOfK), 0));
					for (uint32_t iSpecIdx = 1; iSpecIdx < iAmountOfSpecies; ++iSpecIdx) {
						vector<pair<double, uint64_t>> vTemp(_iNumOfK);
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							const uint64_t& iTempScore = vCount_unique[iSpecIdx + ikMerlength * iAmountOfSpecies];
							vSumOfUniquekMers[ikMerlength] += iTempScore;
							vSumOfNonUniques[ikMerlength] += vCount_all[iSpecIdx + ikMerlength * iAmountOfSpecies];
							vTemp[ikMerlength] = make_pair(vCount_all[iSpecIdx + ikMerlength * iAmountOfSpecies], iTempScore);
						}
						vOut[iSpecIdx] = make_tuple(Utilities::replaceCharacter(mOrganisms[iSpecIdx], ',', ' '), vTemp, mIdxToTax[iSpecIdx]);
					}
					sort(vOut.begin(), vOut.end(), [](const tuple<string, vector<pair<double, uint64_t>>, uint32_t>& a, const tuple<string, vector<pair<double, uint64_t>>, uint32_t>& b) {
						for (uint64_t i = 0; i < get<1>(a).size(); ++i) {
							if (get<1>(a)[i].second == get<1>(b)[i].second) {
								continue;
							}
							else {
								return get<1>(a)[i].second > get<1>(b)[i].second;
							}
						}
						return false;
					});

					// save to file(s)
					//auto orgBuf = cout.rdbuf();
					double allSumOfIdentified = 0;
					if (fTableFile != "") {
						//cout.rdbuf(tableFileStream.rdbuf());
						ofstream tableFileStream((vInputFiles.size() > 1) ? fTableFile + fileName + ".csv" : fTableFile);
						//tableFileStream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
						if (!tableFileStream) {
							throw runtime_error("Profile file couldn't be opened for writing!");
						}

						/*if (bHumanReadable) {
							// short version: taxID,Name,Unique Percentage of highest k,Non-unique Percentage of highest k\n
							bool bBreakOut = false;
							double dSumOfIdentified = 0;
							ostringstream sOutStr;
							tableFileStream << "#tax ID,Name,Unique rel. freq. in %,Non-unique rel. freq. in %,Overall rel. freq. in %\n";
							for (const auto& entry : vOut) {
								if (get<1>(entry)[_iNumOfK - 1].first > 0 && !bBreakOut) {
									sOutStr << get<2>(entry) << "," << get<0>(entry);
									if (get<1>(entry)[0].second == 0) {
										sOutStr << "," << 0.0;
										bBreakOut = true;
									}
									else {
										sOutStr << "," << static_cast<double>(get<1>(entry)[0].second) / vSumOfUniquekMers[0] * 100.0;
									}
									if (get<1>(entry)[0].first == 0) {
										sOutStr << "," << 0.0;
									}
									else {
										sOutStr << "," << static_cast<double>(get<1>(entry)[0].first) / vSumOfNonUniques[0] * 100.0;
									}
									dSumOfIdentified += get<1>(entry)[0].first;
									sOutStr << "," << static_cast<double>(get<1>(entry)[0].first) / (iNumberOfkMersInInput - vNumberOfGarbagekMersPerK[0]) * 100.;
									sOutStr << "\n";
								}
							}

							// last entry
							tableFileStream << "0,not identified,"
											<< "0.0,0.0,"
											<< ((static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[0]) - dSumOfIdentified) / (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[0]))) * 100.;
							tableFileStream << "\n" << sOutStr.str();
							allSumOfIdentified = dSumOfIdentified;
						}
						else {*/
						ostringstream sOutStr;
						// long version: taxID,Name,Unique Counts,Unique rel. freq. x in 0.x,Non-unique Counts,Non-unique rel. freq. x in 0.x\n
						tableFileStream << "#taxID,Name";
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Unique counts k=" << _iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Unique rel. freq. k=" << _iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Non-unique counts k=" << _iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Non-unique rel. freq. k=" << _iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Overall rel. freq. k=" << _iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Overall unique rel. freq. k=" << _iMaxK - ikMerlength;
						}
						tableFileStream << "\n";

						vector<double> vSumOfIdentified(_iNumOfK, 0), vSumOfUniqueIdentified(_iNumOfK, 0);
						for (const auto& entry : vOut) {
							if (get<1>(entry)[_iNumOfK - 1].first > 0) {
								// unique count
								sOutStr << get<2>(entry) << "," << get<0>(entry);
								for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
									sOutStr << "," << get<1>(entry)[ikMerlength].second;
								}
								// unique rel freq
								for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
									if (get<1>(entry)[ikMerlength].second == 0) {
										sOutStr << "," << 0.0;
									}
									else {
										sOutStr << "," << static_cast<double>(get<1>(entry)[ikMerlength].second) / vSumOfUniquekMers[ikMerlength];
									}
								}
								// non-unique count
								for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
									sOutStr << "," << get<1>(entry)[ikMerlength].first;
								}
								// non-unique rel freq
								for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
									if (bSpaced) {
										sOutStr << "," << static_cast<double>(get<1>(entry)[ikMerlength].first) / vSumOfNonUniques[ikMerlength];
									}
									else {
										//sOutStr << "," << static_cast<double>(get<1>(entry)[ikMerlength]) / (iNumberOfkMersInInput - aNonUniqueHits[ikMerlength] - (_iMaxK - _iMinK - ikMerlength) * 6 * iNumOfReadsSum);
										sOutStr << "," << static_cast<double>(get<1>(entry)[ikMerlength].first) / vSumOfNonUniques[ikMerlength];
									}
								}
								// Overall rel freq
								for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
									vSumOfIdentified[ikMerlength] += get<1>(entry)[ikMerlength].first;
									sOutStr << "," << static_cast<double>(get<1>(entry)[ikMerlength].first) / (iNumberOfkMersInInput - vNumberOfGarbagekMersPerK[ikMerlength]);
								}

								// Overall unique rel freq
								for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
									vSumOfUniqueIdentified[ikMerlength] += get<1>(entry)[ikMerlength].second;
									sOutStr << "," << static_cast<double>(get<1>(entry)[ikMerlength].second) / (iNumberOfkMersInInput - vNumberOfGarbagekMersPerK[ikMerlength]);
								}
								sOutStr << "\n";
							}
						}
						// last entry
						tableFileStream << "0,not identified";
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK * 4; ++ikMerlength) {
							// all counts for unique and non-unique relate to the identified number of counts so no value other than 0 can be written here
							tableFileStream << "," << 0.0;
						}

						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							//cout << iNumberOfkMersInInput << " " << vNumberOfGarbagekMersPerK[ikMerlength] << " " << vSumOfIdentified[ikMerlength] << endl;
							tableFileStream << "," << (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]) - static_cast<double>(vSumOfIdentified[ikMerlength])) / (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]));
						}
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							tableFileStream << "," << (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]) - static_cast<double>(vSumOfUniqueIdentified[ikMerlength])) / (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]));
						}
						tableFileStream << "\n" << sOutStr.str();

						allSumOfIdentified = vSumOfIdentified[0];
						//}
					}
					else {
						for (const auto& entry : vOut) {
							allSumOfIdentified += get<1>(entry)[0].first;
						}
					}

#ifdef TIME
					endTIME = std::chrono::high_resolution_clock::now();
					cout << "Save profile " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
					startTIME = std::chrono::high_resolution_clock::now();
#endif

					/*if (fTableFile != "") {
						cout.rdbuf(orgBuf);
					}*/
					if (_bVerbose) {
						cout << "OUT: Number of k-mers in input: " << iNumberOfkMersInInput << " of which " << allSumOfIdentified / (iNumberOfkMersInInput - vNumberOfGarbagekMersPerK[0]) * 100. << " % were identified." << endl;
						cout << "OUT: Number of uniques:";
						for (int32_t j = 0; j < _iNumOfK; ++j) {
							cout << " " << vSumOfUniquekMers[j];
						}
						cout << endl;
						cout << "OUT: Time fastq: " << iTimeFastq << " ns" << endl;
						cout << "OUT: Time compare: " << iTimeCompare << " ns" << endl;
					}

					iNumOfReadsSum = 0;
					iNumberOfkMersInInput = 0;
					iTimeFastq = 0;
					iTimeCompare = 0;
					vReadNameAndLength.clear();
					fill(vNumberOfGarbagekMersPerK.begin(), vNumberOfGarbagekMersPerK.end(), 0);

					for (uint64_t i = 0; i < iMult; ++i) {
						vCount_all[i] = 0.;
						vCount_unique[i] = 0;
					}

					allFilesProgress = transferBetweenRuns->iCurrentOverallPercentage;
				}
#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "End " << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
		/*
		void testC2V_1();
		void testC2V_2();
		void testC2V_3();
		void testC2V_4();
		void testC2V_5();
		void testC2V_6();
		void testC2V_7();
		void testC2V_8();
		void testC2V_9();*/
	};
}