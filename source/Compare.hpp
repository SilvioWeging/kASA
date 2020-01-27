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
#include "Trie.hpp"
#include "WorkerThread.hpp"

namespace kASA {
	class Compare : public Read {

		const bool _bTranslated = false;
		bool bUnfunny = false;
		mutex m_exceptionLock;
		exception_ptr someThingWentWrong;

	public:
		Compare(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const int32_t& iNumOfBeasts, const bool& bVerbose = false, const bool& bProtein = false, const string& stxxl_mode = "", const bool& bUnfunny = false) : Read(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, bProtein, stxxl_mode, bUnfunny), _bTranslated(bProtein), bUnfunny(bUnfunny), iNumOfBeasts(iNumOfBeasts) {}

		// for output
		bool bHumanReadable = false;
		int32_t iNumOfBeasts = 3;

	private:

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
		inline void markTaxIDs(const uint64_t& codedTaxIDs, Utilities::sBitArray& vMemoryOfTaxIDs_k, const unordered_map<uint32_t, uint32_t>& , const vector<packedPair>*) {
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


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Compare as many as #Number-of-processors vectors with an index lying on a HDD/SSD and note all similarities for any k. 
		// To minimize hard disk access, the order is as follows: Get kMer from RAM Vec -> Search in Prefix-Trie -> Get range of possible hit -> binary search in that range -> note if hit
		template <typename vecType>
		inline void compareWithDatabase(const int32_t& iThreadID, const pair<uint64_t, Utilities::rangeContainer>* vIn, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, unique_ptr<float[]>& vReadIDtoGenID, const uint32_t& iSpecIDRange, const unordered_map<uint32_t, uint32_t>& mTaxToIdx, const unordered_map<readIDType, uint64_t>& mReadIDToArrayIdx) {

			try {

				unique_ptr<vector<uint64_t>[]> vReadIDs(new vector<uint64_t>[_iNumOfK]);
				//unique_ptr<uint32_t[]> vReadIDs_(new uint32_t[_iNumOfK * (mReadIDToArrayIdx.size() + 1)]);
				vector<uint64_t> vMemoryOfSeenkMers(_iNumOfK);
				//vector<uint32_t> vMemoryCounterOnly(_iNumOfK, 0);
				vector<Utilities::sBitArray> vMemoryOfTaxIDs(_iNumOfK, Utilities::sBitArray(iSpecIDRange));


				const uint64_t& ivInSize = (_iMinK <= 6) ? vIn->second.kMers_ST6.size() : vIn->second.kMers_GT6.size();

				uint64_t iIdxIn = 0;

				/*for (int32_t iK = 0; iK < _iNumOfK; ++iK) {
					vReadIDs[iK].clear();
					vReadIDs[iK].reserve(100);
					vMemoryOfSeenkMers[iK] = 0;
					//vMemoryCounterOnly[iK] = 0;
					vMemoryOfTaxIDs[iK].clear();
				}*/

				tuple<uint64_t, uint64_t> iSeenInput = make_pair(0, 0);

				const int32_t& ikDifferenceTop = _iHighestK - _iMaxK;

				const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
				auto seenResultIt = libBeginIt;
				bool bInputIterated = false;
				//uint64_t iIdxLib = 0;
				//uint64_t iBitMask = 1152921503533105152ULL;
				//uint64_t iLastRangeBegin = 0;


				while (iIdxIn < ivInSize) {

					const pair<uint64_t, uint32_t>& iCurrentkMer = (_iMinK <= 6) ? static_cast<pair<uint64_t, uint32_t>>(vIn->second.kMers_ST6[iIdxIn]) : static_cast<pair<uint64_t, uint32_t>>(vIn->second.kMers_GT6[iIdxIn]);
					
					int16_t ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1);
					int32_t shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);
					auto iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;

					// If the ending is ^ it's not going to hit anyway, might as well stop here
					if ((iCurrentkMerShifted & 31) == 30) {
						++iIdxIn;
						bInputIterated = true;
						continue;
					}

					////////////////////////////////////////////////////////////////////////

					// Count duplicates too
					if (get<0>(iSeenInput) == get<0>(iCurrentkMer) && bInputIterated) {
						for (int32_t ik = _iNumOfK - 1; ik >= 0; --ik) {
							const int32_t& shift_ = 5 * (_iHighestK - _aOfK[ik]);
							const auto& iCurrentkMerShifted_ = get<0>(iCurrentkMer) >> shift_;
							if (iCurrentkMerShifted_ == vMemoryOfSeenkMers[ik]) {
								//++vMemoryCounterOnly[ikLengthCounter];
								vReadIDs[ik].push_back(mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
								//vReadIDs_[mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second * _iNumOfK + ikLengthCounter]++;
							}
						}


						++iIdxIn;
						bInputIterated = true;
						continue;
					}
					else {
						iSeenInput = iCurrentkMer;
					}

					const auto shiftVal = [&shift, this](const uint64_t& val) { return (_iMinK > 6) ? (val & 1073741823ULL) >> shift : (val >> shift); };
					const auto rangeBeginIt = libBeginIt + vIn->first, rangeEndIt = libBeginIt + static_cast<uint64_t>(vIn->first) + vIn->second.range;


					if (shiftVal(rangeBeginIt->first) == iCurrentkMerShifted) {
						seenResultIt = rangeBeginIt;
					}
					else {
						if (shiftVal(rangeEndIt->first) == iCurrentkMerShifted) {
							// we need the first occurence in the database
							uint64_t iTemp = 1;
							while (shiftVal((rangeEndIt - iTemp)->first) == iCurrentkMerShifted) {
								++iTemp;
							}
							seenResultIt = rangeEndIt - (iTemp - 1);
						}
						else {
							if (iCurrentkMerShifted < shiftVal(rangeBeginIt->first) || iCurrentkMerShifted > shiftVal(rangeEndIt->first)) {
								//cout << kMerToAminoacid(get<0>(iCurrentkMer),12) << " " << kMerToAminoacid(iCurrentkMerShifted,12) << " " << kMerToAminoacid(rangeBeginIt->first, 12) << " " << kMerToAminoacid(shiftVal(rangeBeginIt->first), 12) << " " << kMerToAminoacid(shiftVal(rangeEndIt->first), 12) << endl;
								++iIdxIn;
								bInputIterated = true;
								continue;
							}
							else {
								seenResultIt = lower_bound(rangeBeginIt, rangeEndIt + 1, iCurrentkMerShifted, [&shift, this](const decltype(*libBeginIt)& a, const uint64_t& val) { return (_iMinK > 6) ? ((a.first & 1073741823ULL) >> shift) < val : (a.first >> shift) < val; });
							}
						}
					}

					bool bBreakOut = false;
					while (seenResultIt != rangeEndIt + 1 && !bBreakOut) {

						const tuple<uint64_t, uint32_t>& iCurrentLib = (_iMinK > 6) ? static_cast<tuple<uint64_t, uint32_t>>(make_tuple(seenResultIt->first & 1073741823ULL, seenResultIt->second)) : tuple<uint64_t, uint32_t>(make_tuple(seenResultIt->first, seenResultIt->second));

						for (; ikLengthCounter >= 0; --ikLengthCounter) {

							shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);
							iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;

							// No matching of stupid stuff!
							if ((iCurrentkMerShifted & 31) == 30) {
								bBreakOut = true;
								break;
							}


							const auto& iCurrentLibkMerShifted = get<0>(iCurrentLib) >> shift;

							if (iCurrentkMerShifted < iCurrentLibkMerShifted) {
								bBreakOut = true;
								break;
							}
							else {
								if (!(iCurrentLibkMerShifted < iCurrentkMerShifted)) {

									if (iCurrentLibkMerShifted == vMemoryOfSeenkMers[ikLengthCounter]) {

										//++vMemoryCounterOnly[ikLengthCounter];
										markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));
										vReadIDs[ikLengthCounter].push_back(mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
										//vReadIDs_[mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second * _iNumOfK + ikLengthCounter]++;
									}
									else {
										const auto& numOfEntries = vMemoryOfTaxIDs[ikLengthCounter].numOfEntries();
										auto it = vMemoryOfTaxIDs[ikLengthCounter].begin();
										it.SetNumOfEntries(numOfEntries);
										for (; it != vMemoryOfTaxIDs[ikLengthCounter].end() && numOfEntries != 0; ++it) {
											const auto& tempIndex = iSpecIDRange * _iNumOfK*iThreadID + (*it)*_iNumOfK + ikLengthCounter;
											const auto& numOfHits = vReadIDs[ikLengthCounter].size();
											//#pragma omp atomic
											vCount[tempIndex] += double(numOfHits) / numOfEntries;

											if (numOfEntries == 1) {
												//#pragma omp atomic
												vCountUnique[tempIndex] += numOfHits;
											}

											const auto& entry = *it;
											const auto& weight = arrWeightingFactors[ikDifferenceTop + ikLengthCounter];

											for (const auto& readID : vReadIDs[ikLengthCounter]) {
												//const auto& readIDScore = 1.f / vReadIDs_[readID * _iNumOfK + ikLengthCounter];
												//if (readIDScore) {
												const auto& score = weight * (1.f / numOfEntries);//(1.f + logf(static_cast<float>(numOfEntries)));
												const auto& arrayIdx = readID * iSpecIDRange + entry;
												//#pragma omp atomic
												vReadIDtoGenID[arrayIdx] += score;

												//}
											}
										}

										//for (const auto& readID : vReadIDs[ikLengthCounter]) {
										//	vReadIDs_[readID * _iNumOfK + ikLengthCounter] = 0;
										//}

										vReadIDs[ikLengthCounter].clear();
										vReadIDs[ikLengthCounter].push_back(mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second);
										//vReadIDs_[mReadIDToArrayIdx.find(get<1>(iCurrentkMer))->second * _iNumOfK + ikLengthCounter] = 1;

										vMemoryOfTaxIDs[ikLengthCounter].clear();
										markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));

										//vMemoryCounterOnly[ikLengthCounter] = 1;

										vMemoryOfSeenkMers[ikLengthCounter] = iCurrentLibkMerShifted;
									}

								}
								else {
									uint64_t iTempCounter = 1;
									const uint64_t& iCurrentLibSuffix = seenResultIt->first;
									while (seenResultIt + iTempCounter != rangeEndIt + 1) {
										const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + iTempCounter)->first);
										int16_t iUntilK = static_cast<int16_t>(_iNumOfK - 1);
										for (; iUntilK > ikLengthCounter; --iUntilK) {
											if ((iCurrentLibSuffix >> 5 * iUntilK) == (iNextLibSuffix >> 5 * iUntilK)) {
												markTaxIDs((seenResultIt + iTempCounter)->second, vMemoryOfTaxIDs[iUntilK], mTaxToIdx, getVec(vLib, iThreadID));
											}
											else {
												break;
											}
										}
										if (iCurrentLibSuffix == iNextLibSuffix) {
											++iTempCounter;
										}
										else {
											break;
										}
									}
									//iIdxLib += iTempCounter;
									seenResultIt += iTempCounter;
									bInputIterated = false;
									break;
								}
							}
						}
						// loop through to find other hits in the library
						if (ikLengthCounter == -1) {
							uint64_t iTempCounter = 1;
							const uint64_t& iCurrentLibSuffix = static_cast<uint64_t>(seenResultIt->first);
							while (seenResultIt + iTempCounter != rangeEndIt + 1) {
								const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + iTempCounter)->first);
								const auto& iNextLibIdx = (seenResultIt + iTempCounter)->second;
								if (iCurrentLibSuffix == iNextLibSuffix) {
									for (int16_t ikLengthCounter_ = static_cast<int16_t>(_iNumOfK - 1); ikLengthCounter_ > ikLengthCounter; --ikLengthCounter_) {
										markTaxIDs(iNextLibIdx, vMemoryOfTaxIDs[ikLengthCounter_], mTaxToIdx, getVec(vLib, iThreadID)); // to identify multiple hits 
									}

									++iTempCounter;
								}
								else {
									bBreakOut = true;
									break;
								}
							}

							seenResultIt += iTempCounter;
						}
					}
					++iIdxIn;
					bInputIterated = true;
				}


				// Don't forget the last saved part
				for (int16_t ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1); ikLengthCounter >= 0; --ikLengthCounter) {
					const auto& numOfEntries = vMemoryOfTaxIDs[ikLengthCounter].numOfEntries();
					auto it = vMemoryOfTaxIDs[ikLengthCounter].begin();
					it.SetNumOfEntries(numOfEntries);
					for (; it != vMemoryOfTaxIDs[ikLengthCounter].end() && numOfEntries != 0; ++it) {
						const auto& tempIndex = (*it)*_iNumOfK + ikLengthCounter;
						const auto& numOfHits = vReadIDs[ikLengthCounter].size();
						//#pragma omp atomic
						vCount[tempIndex] += float(numOfHits) / numOfEntries;

						if (numOfEntries == 1) {
							//#pragma omp atomic
							vCountUnique[tempIndex] += numOfHits;
						}

						const auto& entry = *it;
						const auto& weight = arrWeightingFactors[ikDifferenceTop + ikLengthCounter];

						for (const auto& readID : vReadIDs[ikLengthCounter]) {
							//const auto& readIDScore = 1.f / vReadIDs_[readID * _iNumOfK + ikLengthCounter];
							//if (readIDScore) {
							const auto& score = weight * (1.f / numOfEntries);//(1.f + logf(static_cast<float>(numOfEntries)));
							const auto& arrayIdx = readID * iSpecIDRange + entry;
							//#pragma omp atomic
							vReadIDtoGenID[arrayIdx] += score;

							//}
						}
					}
					//for (const auto& readID : vReadIDs[ikLengthCounter]) {
					//	vReadIDs_[readID * _iNumOfK + ikLengthCounter] = 0;
					//}
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
		inline void createProfile(const int32_t iThreadID, const pair<uint64_t, Utilities::rangeContainer>* vIn, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, const uint32_t& iSpecIDRange, const unordered_map<uint32_t, uint32_t>& mTaxToIdx) {

			try {

				vector<uint64_t> vMemoryOfSeenkMers(_iNumOfK);
				vector<uint32_t> vMemoryCounterOnly(_iNumOfK, 0);
				vector<Utilities::sBitArray> vMemoryOfTaxIDs(_iNumOfK, Utilities::sBitArray(iSpecIDRange));


				const uint64_t& ivInSize = (_iMinK <= 6) ? vIn->second.kMers_ST6.size() : vIn->second.kMers_GT6.size();

				uint64_t iIdxIn = 0;

				/*for (int32_t iK = 0; iK < _iNumOfK; ++iK) {
					vMemoryOfSeenkMers[iK] = 0;
					vMemoryCounterOnly[iK] = 0;
					vMemoryOfTaxIDs[iK].clear();
				}*/

				tuple<uint64_t, uint64_t> iSeenInput = make_pair(0, 0);

				const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
				auto seenResultIt = libBeginIt;

				bool bInputIterated = false;
				//uint64_t iIdxLib = 0;
				//uint64_t iLastRangeBegin = 0;

				while (iIdxIn < ivInSize) {

					const pair<uint64_t, uint64_t>& iCurrentkMer = (_iMinK <= 6) ? static_cast<pair<uint64_t, uint64_t>>(vIn->second.kMers_ST6[iIdxIn]) : static_cast<pair<uint64_t, uint64_t>>(vIn->second.kMers_GT6[iIdxIn]);

					int16_t ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1);
					int32_t shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);
					auto iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;

					// If the ending is ^ it's not going to hit anyway, might as well stop here
					if ((iCurrentkMerShifted & 31) == 30) {
						++iIdxIn;
						bInputIterated = true;
						continue;
					}

					////////////////////////////////////////////////////////////////////////

					// Count duplicates too
					if (get<0>(iSeenInput) == get<0>(iCurrentkMer) && bInputIterated) {
						for (int32_t ik = _iNumOfK - 1; ik >= 0; --ik) {
							const int32_t& shift_ = 5 * (_iHighestK - _aOfK[ik]);
							const auto& iCurrkMerShifted_ = get<0>(iCurrentkMer) >> shift_;
							if (iCurrkMerShifted_ == vMemoryOfSeenkMers[ik]) {
								++vMemoryCounterOnly[ik];
							}
						}

						++iIdxIn;
						bInputIterated = true;
						continue;
					}
					else {
						iSeenInput = iCurrentkMer;
					}


					const auto shiftVal = [&shift, this](const uint64_t& val) { return (_iMinK > 6) ? (val & 1073741823ULL) >> shift : (val >> shift); };
					const auto rangeBeginIt = libBeginIt + vIn->first, rangeEndIt = libBeginIt + static_cast<uint64_t>(vIn->first) + vIn->second.range;

					if (shiftVal(rangeBeginIt->first) == iCurrentkMerShifted) {
						seenResultIt = rangeBeginIt;
					}
					else {
						if (shiftVal(rangeBeginIt->first) == iCurrentkMerShifted) {
							// we need the first occurence in the database
							uint64_t iTemp = 1;
							while (shiftVal((rangeEndIt - iTemp)->first) == iCurrentkMerShifted) {
								++iTemp;
							}
							seenResultIt = rangeEndIt - (iTemp - 1);
						}
						else {
							if (iCurrentkMerShifted < shiftVal(rangeBeginIt->first) || iCurrentkMerShifted > shiftVal(rangeEndIt->first)) {
								++iIdxIn;
								bInputIterated = true;
								continue;
							}
							else {
								seenResultIt = lower_bound(rangeBeginIt, rangeEndIt + 1, iCurrentkMerShifted, [&shift, this](const decltype(*libBeginIt)& a, const uint64_t& val) { return (_iMinK > 6) ? ((a.first & 1073741823ULL) >> shift) < val : (a.first >> shift) < val; });
							}
						}
					}

					bool bBreakOut = false;
					while (seenResultIt != rangeEndIt + 1 && !bBreakOut) {
						const tuple<uint64_t, uint32_t>& iCurrentLib = (_iMinK > 6) ? static_cast<tuple<uint64_t, uint32_t>>(make_tuple(seenResultIt->first & 1073741823ULL, seenResultIt->second)) : tuple<uint64_t, uint32_t>(make_tuple(seenResultIt->first, seenResultIt->second));
						//cout << get<0>(iCurrentLib) << " " << get<1>(iCurrentLib) << endl;

						for (; ikLengthCounter >= 0; --ikLengthCounter) {

							shift = 5 * (_iHighestK - _aOfK[ikLengthCounter]);
							iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;
							
							// No matching of stupid stuff!
							if ((iCurrentkMerShifted & 31) == 30) {
								bBreakOut = true;
								break;
							}
							
							const auto& iCurrentLibkMerShifted = get<0>(iCurrentLib) >> shift;

							if (iCurrentkMerShifted < iCurrentLibkMerShifted) {
								bBreakOut = true;
								break;
							}
							else {
								if (!(iCurrentLibkMerShifted < iCurrentkMerShifted)) {
									/*if (!bInputIterated) {
										if (iCurrentLibkMerShifted == vMemoryOfSeenkMers[ikLengthCounter] && get<1>(iCurrentLib) == vMemoryOfSeenTaxIDs[ikLengthCounter]) {
											continue;
										}
										else {
											vMemoryOfSeenTaxIDs[ikLengthCounter] = get<1>(iCurrentLib);
										}
									}*/

									if (iCurrentLibkMerShifted == vMemoryOfSeenkMers[ikLengthCounter]) {
										++vMemoryCounterOnly[ikLengthCounter];
										markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));
									}
									else {
										const auto& numOfEntries = vMemoryOfTaxIDs[ikLengthCounter].numOfEntries();
										auto it = vMemoryOfTaxIDs[ikLengthCounter].begin();
										it.SetNumOfEntries(numOfEntries);
										for (; it != vMemoryOfTaxIDs[ikLengthCounter].end() && numOfEntries != 0; ++it) {
											const auto& tempIndex = iSpecIDRange * _iNumOfK*iThreadID + (*it)*_iNumOfK + ikLengthCounter;

											//#pragma omp atomic
											vCount[tempIndex] += double(vMemoryCounterOnly[ikLengthCounter]) / numOfEntries;


											if (numOfEntries == 1) {
												//#pragma omp atomic
												vCountUnique[tempIndex] += vMemoryCounterOnly[ikLengthCounter];
											}
										}

										vMemoryOfTaxIDs[ikLengthCounter].clear();
										markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));

										vMemoryCounterOnly[ikLengthCounter] = 1;

										vMemoryOfSeenkMers[ikLengthCounter] = iCurrentLibkMerShifted;
									}
								}
								else {
									uint64_t iTempCounter = 1;
									const uint64_t& iCurrentLibSuffix = static_cast<uint64_t>(seenResultIt->first);
									while (seenResultIt + iTempCounter != rangeEndIt + 1) {
										const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + iTempCounter)->first);
										int16_t iUntilK = static_cast<int16_t>(_iNumOfK - 1);
										for (; iUntilK > ikLengthCounter; --iUntilK) {
											if ((iCurrentLibSuffix >> 5 * iUntilK) == (iNextLibSuffix >> 5 * iUntilK)) {
												markTaxIDs((seenResultIt + iTempCounter)->second, vMemoryOfTaxIDs[iUntilK], mTaxToIdx, getVec(vLib, iThreadID));
											}
											else {
												break;
											}
										}
										if (iCurrentLibSuffix == iNextLibSuffix) {
											++iTempCounter;
										}
										else {
											break;
										}
									}
									//iIdxLib += iTempCounter;
									seenResultIt += iTempCounter;
									bInputIterated = false;
									break;
								}
							}
						}
						// loop through to find other hits in the library
						if (ikLengthCounter == -1) {
							uint64_t iTempCounter = 1;
							const uint64_t& iCurrentLibSuffix = static_cast<uint64_t>(seenResultIt->first);
							while (seenResultIt + iTempCounter != rangeEndIt + 1) {
								const uint64_t& iNextLibSuffix = static_cast<uint64_t>((seenResultIt + iTempCounter)->first);
								const auto& iNextLibIdx = (seenResultIt + iTempCounter)->second;
								if (iCurrentLibSuffix == iNextLibSuffix) {
									for (int16_t ikLengthCounter_ = static_cast<int16_t>(_iNumOfK - 1); ikLengthCounter_ > ikLengthCounter; --ikLengthCounter_) {
										markTaxIDs(iNextLibIdx, vMemoryOfTaxIDs[ikLengthCounter_], mTaxToIdx, getVec(vLib, iThreadID)); // to identify multiple hits  
									}

									++iTempCounter;
								}
								else {
									bBreakOut = true;
									break;
								}
							}

							seenResultIt += iTempCounter;
						}
					}
					++iIdxIn;
					bInputIterated = true;
				}


				// Don't forget the last saved part
				for (int16_t ikLengthCounter = static_cast<int16_t>(_iNumOfK - 1); ikLengthCounter >= 0; --ikLengthCounter) {
					const auto& numOfEntries = vMemoryOfTaxIDs[ikLengthCounter].numOfEntries();
					auto it = vMemoryOfTaxIDs[ikLengthCounter].begin();
					it.SetNumOfEntries(numOfEntries);
					for (; it != vMemoryOfTaxIDs[ikLengthCounter].end() && numOfEntries != 0; ++it) {
						const auto& tempIndex = (*it)*_iNumOfK + ikLengthCounter;
						//#pragma omp atomic
						vCount[tempIndex] += float(vMemoryCounterOnly[ikLengthCounter]) / numOfEntries;

						if (numOfEntries == 1) {
							//#pragma omp atomic
							vCountUnique[tempIndex] += vMemoryCounterOnly[ikLengthCounter];
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
		inline void compareWithDatabase_sloppy(const int32_t& iThreadID, const pair<uint64_t, Utilities::rangeContainer>* vIn, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, unique_ptr<float[]>& vReadIDtoGenID, const uint32_t& iSpecIDRange, const unordered_map<readIDType, uint64_t>& mReadIDToArrayIdx) {

			try {

				//auto start = std::chrono::high_resolution_clock::now();

				const uint64_t& ivInSize = vIn->second.kMers_ST6.size();

				const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
				auto seenResultIt = libBeginIt;

				const auto rangeBeginIt = libBeginIt + vIn->first, rangeEndIt = libBeginIt + static_cast<uint64_t>(vIn->first) + vIn->second.range;
				seenResultIt = rangeBeginIt;

				while (seenResultIt != rangeEndIt + 1) {
					const auto& iCurrentLib = *seenResultIt;

					vCount[iCurrentLib] += static_cast<double>(ivInSize) / (vIn->second.range + 1);


					if (vIn->second.range == 0) {
						vCountUnique[iCurrentLib] += ivInSize;
					}

					const auto& score = arrWeightingFactors[0] * (1.f / ivInSize);
					for (const auto& elem : vIn->second.kMers_ST6) {
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
		inline void createProfile_sloppy(const int32_t iThreadID, const pair<uint64_t, Utilities::rangeContainer>* vIn, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique) {

			try {


				const uint64_t& ivInSize = vIn->second.kMers_ST6.size();

				const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
				auto seenResultIt = libBeginIt;

				const auto rangeBeginIt = libBeginIt + vIn->first, rangeEndIt = libBeginIt + static_cast<uint64_t>(vIn->first) + vIn->second.range;
				seenResultIt = rangeBeginIt;

				while (seenResultIt != rangeEndIt + 1) {
					const auto& iCurrentLib = *seenResultIt;
					//#pragma omp atomic
					vCount[iCurrentLib] += static_cast<double>(ivInSize) / (vIn->second.range + 1);


					if (vIn->second.range == 0) {
						//#pragma omp atomic
						vCountUnique[iCurrentLib] += ivInSize;
					}
					++seenResultIt;
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
		inline void scoringFunc(const unique_ptr<float[]>& vReadIDtoGenID, const uint64_t& iStart, const uint64_t& iEnd, const uint64_t& iRealReadIDStart, list<pair<string, readIDType>>& vReadNameAndLength, const unique_ptr<uint64_t[]>& mFrequencies, const vector<uint32_t>& mIdxToTax, const vector<string>& mOrganisms, ofstream&& fOut) {
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
							get<0>(resultVec[iSpecIdx]) = iSpecIdx;
							get<1>(resultVec[iSpecIdx]) = kMerScore; 
							get<2>(resultVec[iSpecIdx]) = relativeScore;
							++iCountOfHits;
						}
					}
					
					if (iCountOfHits == 0) {
						if (bHumanReadable) {
							fOut << iRealReadIDStart + readIdx << "\t" << currentReadLengthAndName.first << "\t-\t-\t-" << "\n";
						}
						else {
							if (readIdx == 0) {
								fOut << "{" << "\n";
							}
							else {
								fOut << "," << "\n" << "{" << "\n";
							}
							fOut << "\t\"Read number\": " << readIdx << "," << "\n" << "\t\"Specifier from input file\": \"" + currentReadLengthAndName.first + "\"," << "\n" << "\t\"Matched taxa\": [" << "\n" << "\t]" << "\n" << "}";
						}
					}
					else {

						partial_sort(resultVec.begin(), resultVec.begin() + iNumOfBeasts, resultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<2>(a) > get<2>(b); });

						//cout << iRealReadIDStart + readIdx << endl;

						if (bHumanReadable) {
							string sOut = "", sOut2 = "", sOut3 = "";

							sOut += to_string(iRealReadIDStart + readIdx) + "\t" + currentReadLengthAndName.first + "\t";
							auto it = resultVec.begin();
							float iValueBefore = 0;
							for (int16_t j = 0, i = 0; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								sOut += to_string(mIdxToTax[get<0>(*it)]) + ";";
								ostringstream e_value;
								e_value.precision(5);
								e_value << std::scientific << get<2>(*it) << "," << std::defaultfloat << get<1>(*it);

								sOut2 += mOrganisms[get<0>(*it)] + ";";
								sOut3 += e_value.str() + ";";

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
						}
						else {
							// json
							ostringstream sOutStr;
							if (iRealReadIDStart + readIdx == 0) {
								sOutStr << "{" << "\n";
							}
							else {
								sOutStr << "," << "\n" << "{" << "\n";
							}

							sOutStr << "\t\"Read number\": " << iRealReadIDStart + readIdx << ",\n" << "\t\"Specifier from input file\": \"" + currentReadLengthAndName.first + "\",\n" << "\t\"Matched taxa\": [\n";
							auto it = resultVec.begin();
							float iValueBefore = 0;
							for (int16_t j = 0, i = 0; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								if (j == 0) {
									sOutStr << "\t{\n";
								}
								else {
									sOutStr << ",\n\t{\n";
								}

								sOutStr << "\t\t\"tax ID\": \"" << mIdxToTax[get<0>(*it)] << "\",\n"
									<< "\t\t\"Name\": \"" << mOrganisms[get<0>(*it)] << "\",\n"
									<< "\t\t\"k-mer Score\": " << std::defaultfloat << get<1>(*it) << ",\n"
									<< "\t\t\"Relative Score\": " << std::scientific << get<2>(*it) << ",\n"
									<< "\t\t\"Error\": " << std::defaultfloat << (bestScore - get<1>(*it)) / bestScore << "\n"
									<< "\t}";

								if (iValueBefore != get<1>(*it)) {
									iValueBefore = get<1>(*it);
									++j;
								}
							}

							sOutStr << "\n\t]\n}";
							fOut << sOutStr.str();
						}

						fill(resultVec.begin(), resultVec.begin() + iCountOfHits, make_tuple(0ULL, 0.f, 0.0));
					}

					vReadNameAndLength.pop_front();
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl;
				throw;
			}
		}


		inline void scoringFunc(vector<tuple<readIDType, float, double>>&& vTempResultVec, const uint64_t& iReadNum, const pair<string, uint32_t>& vReadNameAndLength, const unique_ptr<uint64_t[]>& mFrequencies, const vector<uint32_t>& mIdxToTax, const vector<string>& mOrganisms, ofstream&& fOut) {
			try {
				float bestScore = 0.f;
				for (int32_t i = _iMinK; i <= _iMaxK; ++i) {
					bestScore += (vReadNameAndLength.second - i * 3 + 1)*arrWeightingFactors[_iHighestK - i];
					//cout << (vReadNameAndLength.second - i * 3 + 1)*arrWeightingFactors[_iHighestK - i] << " " << vReadNameAndLength.second << " " << i << " " << vReadNameAndLength.second - i * 3 + 1 <<" " << arrWeightingFactors[_iHighestK - i] << endl;
				}

				for (auto it = vTempResultVec.begin(); it != vTempResultVec.end();  ++it) {
					if (_bTranslated) {
						get<2>(*it) = double(get<1>(*it)) / (1.0 + log2(mFrequencies[get<0>(*it)] * double(vReadNameAndLength.second - _iHighestK + 1)));
					}
					else {
						get<2>(*it) = double(get<1>(*it)) / (1.0 + log2(mFrequencies[get<0>(*it)] * double(vReadNameAndLength.second - _iHighestK * 3 + 1)));
					}
				}

				partial_sort(vTempResultVec.begin(), vTempResultVec.begin() + iNumOfBeasts, vTempResultVec.end(), [](const tuple<uint64_t, float, double>& p1, const tuple<uint64_t, float, double>& p2) { return get<2>(p1) > get<2>(p2); });

				if (bHumanReadable) {
					string sOut = "", sOut2 = "", sOut3 = "";

					sOut += to_string(iReadNum) + "\t" + vReadNameAndLength.first + "\t";
					auto it = vTempResultVec.begin();
					float iValueBefore = 0;
					for (int16_t j = 0; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
						sOut += to_string(mIdxToTax[get<0>(*it)]) + ";";
						ostringstream e_value;
						e_value.precision(5);
						e_value << std::scientific << get<2>(*it) << "," << std::defaultfloat << get<1>(*it);

						sOut2 += mOrganisms[get<0>(*it)] + ";";
						sOut3 += e_value.str() + ";";

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
				}
				else {
					// json
					ostringstream sOutStr;
					if (iReadNum == 0) {
						sOutStr << "{" << "\n";
					}
					else {
						sOutStr << "," << "\n" << "{" << "\n";
					}

					sOutStr << "\t\"Read number\": " << iReadNum << ",\n" << "\t\"Specifier from input file\": \"" + vReadNameAndLength.first + "\",\n" << "\t\"Matched taxa\": [\n";
					auto it = vTempResultVec.begin();
					float iValueBefore = 0;
					for (int16_t j = 0; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
						if (j == 0) {
							sOutStr << "\t{\n";
						}
						else {
							sOutStr << ",\n\t{\n";
						}
						 
						 sOutStr << "\t\t\"tax ID\": \"" << mIdxToTax[get<0>(*it)] << "\",\n"
							 << "\t\t\"Name\": \"" << mOrganisms[get<0>(*it)] << "\",\n"
							 << "\t\t\"k-mer Score\": " << std::defaultfloat << get<1>(*it) << ",\n" 
							 << "\t\t\"Relative Score\": " << std::scientific << get<2>(*it) << ",\n"
							 << "\t\t\"Error\": " << std::defaultfloat << (bestScore - get<1>(*it))/bestScore << "\n"
							 << "\t}";
						
						if (iValueBefore != get<1>(*it)) {
							iValueBefore = get<1>(*it);
							++j;
						}
					}

					sOutStr << "\n\t]\n}";
					fOut << sOutStr.str();
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl;
				throw;
			}
		}

	public:
		/////////////////////////////////////////////////////////////////////////////////
		void CompareWithLib_partialSort(const string& contentFile, const string& sLibFile, const string& fInFile, const string& fOutFile, const string& fTableFile, const uint8_t& iTrieDepth, const uint64_t& iMemory, const bool& bSpaced, bool bRAM, const bool& bUnique) {
			
			try {
				// test if files exists
				if (!ifstream(contentFile) || !ifstream(sLibFile) || !ifstream(sLibFile + "_f.txt") || !ifstream(sLibFile + "_trie.txt")) {
					throw runtime_error("One of the files does not exist");
				}

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
						mOrganisms.push_back(tempLineContent[0]);
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

				//unique_ptr<unique_ptr<stxxlFile>[]> stxxlLibFile;
				unique_ptr<stxxlFile> stxxlLibFile(new stxxlFile(sLibFile, stxxl::file::RDONLY));
				unique_ptr<unique_ptr<const contentVecType_32p>[]> vLib;
				unique_ptr<unique_ptr<const index_t_p>[]> vLibParted_p;
				vector<packedBigPair> vLib_RAM_Full;
				vector<packedPair> vLib_RAM_Half;
				vector<uint16_t> vLib_RAM_taxaOnly;

				unique_ptr<unique_ptr<const taxaOnly>[]> vLib_taxaOnly;

				//vector<tuple<uint64_t, uint32_t>> vLibInRAM;

				uint64_t iBytesUsedByVectors = 0;
				if (bPartitioned) {
					vLibParted_p.reset(new unique_ptr<const index_t_p>[_iNumOfThreads]);
					//stxxlLibFile.reset(new unique_ptr<stxxlFile>[_iNumOfThreads]);
					for (int32_t i = 0; i < _iNumOfThreads; ++i) {
						//stxxlLibFile[i].reset(new stxxlFile(sLibFile + "_"+ to_string(i), stxxl::file::RDONLY));
						vLibParted_p[i].reset(new const index_t_p(stxxlLibFile.get(), iSizeOfLib));
					}
					iBytesUsedByVectors = _iNumOfThreads * index_t_p::block_size * index_t_p::page_size * (vLibParted_p[0])->numpages();
				}
				else {
					if (bUnfunny) {
						vLib_taxaOnly.reset(new unique_ptr<const taxaOnly>[_iNumOfThreads]);
						for (int32_t i = 0; i < _iNumOfThreads; ++i) {
							vLib_taxaOnly[i].reset(new const taxaOnly(stxxlLibFile.get(), iSizeOfLib));
						}
						iBytesUsedByVectors = _iNumOfThreads * taxaOnly::block_size * taxaOnly::page_size * (vLib_taxaOnly[0])->numpages();
					}
					else {
						vLib.reset(new unique_ptr<const contentVecType_32p>[_iNumOfThreads]);
						for (int32_t i = 0; i < _iNumOfThreads; ++i) {
							vLib[i].reset(new const contentVecType_32p(stxxlLibFile.get(), iSizeOfLib));
						}
						iBytesUsedByVectors = _iNumOfThreads * contentVecType_32p::block_size * contentVecType_32p::page_size * (vLib[0])->numpages();
					}
				}

				if (bRAM) {
					try {
						if (bPartitioned) {
							if ((iSizeOfLib * sizeof(packedPair) + 2048ULL * 1024 * 1024) >= iMemory) {
								cerr << "ERROR: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
								bRAM = false;
							}
							else {
								vLib_RAM_Half.reserve(iSizeOfLib);
								stxxl::vector_bufreader<index_t_p::const_iterator> bufferedReader(vLibParted_p[0]->cbegin(), vLibParted_p[0]->cend(), 0);
								for (; !bufferedReader.empty(); ++bufferedReader) {
									vLib_RAM_Half.push_back(*bufferedReader);
								}

								for (int32_t i = 0; i < _iNumOfThreads; ++i) {
									vLibParted_p[i].reset();
								}
								vLibParted_p.reset();
								iBytesUsedByVectors = iSizeOfLib * sizeof(packedPair);
							}
						}
						else {
							if (bUnfunny) {
								if ((iSizeOfLib * sizeof(uint16_t) + 2048ULL * 1024 * 1024) >= iMemory) {
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
									iBytesUsedByVectors = iSizeOfLib * sizeof(uint16_t);
								}
							}
							else {
								if ((iSizeOfLib * sizeof(packedBigPair) + 2048ULL * 1024 * 1024) >= iMemory) {
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
									iBytesUsedByVectors = iSizeOfLib * sizeof(packedBigPair);
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

				// Create threadpool(s), in stxxl mode we can only have synced parallelism (as in no two threads should not access the same vector instance)
				vector<WorkerThread> workerThreadPool(_iNumOfThreads);
				

				// load Trie
				const uint8_t& iTD = iTrieDepth;
				Trie T(static_cast<int8_t>(_iMaxK), static_cast<int8_t>(_iMinK), iTD, _iNumOfThreads, (_iMinK >= 6) && (iMemory - (_iNumOfK * uint64_t(iAmountOfSpecies) * 8 + iBytesUsedByVectors + _iNumOfThreads * Utilities::sBitArray(iAmountOfSpecies).sizeInBytes())) > (kASA::kASA::aminoacidTokMer("]^^^^^")*sizeof(packedBigPair)+1024ULL*1024ULL*1024ULL));
				T.LoadFromStxxlVec(sLibFile);
				T.SetForIsInTrie( (_iMinK < 6) ? static_cast<uint8_t>(_iMinK) : static_cast<uint8_t>(6));

				if (_bVerbose) {
					T.GetIfVecIsUsed();
				}

				// This holds the hits for each organism

				const uint64_t& iMult = _iNumOfThreads * _iNumOfK * uint64_t(iAmountOfSpecies);
				unique_ptr<double[]> vCount_all(new double[iMult]);
				unique_ptr<uint64_t[]> vCount_unique(new uint64_t[iMult]);

				for (uint64_t i = 0; i < iMult; ++i) {
					vCount_all[i] = 0.;
					vCount_unique[i] = 0;
				}

				uint64_t iTimeFastq = 0, iTimeCompare = 0, iNumOfReads = 0, iNumOfReadsOld = 0, iNumOfReadsSum = 0, iSoftMaxSizeOfInputVecs = 0;
				
				// Set memory boundaries
				if (iMemory > T.GetSize() + iMult * 8 + iBytesUsedByVectors + _iNumOfThreads * Utilities::sBitArray(iAmountOfSpecies).sizeInBytes()) {
					iSoftMaxSizeOfInputVecs = static_cast<uint64_t>(0.9 * (iMemory - T.GetSize() - iMult * sizeof(double) - iMult * sizeof(uint64_t) - iBytesUsedByVectors - _iNumOfThreads * Utilities::sBitArray(iAmountOfSpecies).sizeInBytes()));
				}
				else {
					cerr << "ERROR: Not enough Memory given, try to download more RAM. Setting to 1GB..." << endl;
					//iSoftMaxSizeOfInputVecs = 1024ull * 1024ull * 1ull / (sizeof(tuple<uint64_t, uint64_t>)*_iNumOfThreads);
					iSoftMaxSizeOfInputVecs = static_cast<uint64_t>(0.9 * 1024ull * 1024ull * 1024ull); // 1024ull * 1024ull * 1024ull
				}

				unordered_map<uint64_t, Utilities::rangeContainer> vInputMap;

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
				vector<string> vInputFiles = pathAndSize.first;
				size_t overallFileSize = pathAndSize.second, allFilesProgress = 0, charsReadOverall = 0;

				
				// allow multiple input files
				for (const auto& inFile : vInputFiles) {

					if (_bVerbose) {
						cout << "OUT: Current file: " << inFile << endl;
					}

					
					string fileName = "";
					if (vInputFiles.size() > 1) { // get file name without path and ending
						const auto& vRawNameSplit = Utilities::split(inFile.substr(fInFile.size(), inFile.size()), '.');
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
					bool isGzipped = (inFile[inFile.length() - 3] == '.' && inFile[inFile.length() - 2] == 'g' && inFile[inFile.length() - 1] == 'z');
					bool bIsGood = false, bIsFasta = false;

					unique_ptr<ifstream> fast_q_a_File;
					unique_ptr<igzstream> fast_q_a_File_gz;
					uint64_t iFileLength = 0;
					
					if (isGzipped) {
						if (_bVerbose) {
							cout << "OUT: File is gzipped, no progress output can be shown." << endl;
						}

						fast_q_a_File_gz.reset(new igzstream(inFile.c_str()));
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
						fast_q_a_File.reset(new ifstream(inFile));
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
					vector<char> mybuffer(16777216);
					fOut.rdbuf()->pubsetbuf(&mybuffer[0], 16777216);
					//fOut.exceptions(std::ifstream::failbit | std::ifstream::badbit);
					if (bReadIDsAreInteresting) {
						fOut.open((vInputFiles.size() > 1) ? fOutFile + fileName + ((bHumanReadable) ? ".rtt" : ".json") : fOutFile); // in case of multiple input files, specify only beginning of the output and the rest will be appended
						if (fOut) {
							if (bHumanReadable) {
								fOut << "#Read number\tSpecifier from input file\tMatched taxa\tNames\tScores{relative,k-mer}\tError" << "\n";
							}
							else {
								fOut << "[" << "\n";
							}
						}
						else {
							throw runtime_error("Readwise output file could not be created!");
						}
					}
					
					
					unique_ptr<strTransfer> transferBetweenRuns(new strTransfer);
					transferBetweenRuns->iCurrentOverallPercentage = allFilesProgress;
					transferBetweenRuns->iNumOfAllCharsRead = charsReadOverall;
					vector<tuple<readIDType, float, double>> vSavedScores;
					readIDType iReadIDofSavedScores = 0;

					// read input
					while (bIsGood) {
						auto start = std::chrono::high_resolution_clock::now();
						if (bIsFasta) {
							if (isGzipped) {
								iNumberOfkMersInInput += readFasta_partialSort(*fast_q_a_File_gz, vInputMap, vReadNameAndLength, iSoftMaxSizeOfInputVecs, iAmountOfSpecies, bSpaced, iFileLength, iFileLength, bReadIDsAreInteresting, transferBetweenRuns, T, vNumberOfGarbagekMersPerK);
							}
							else {
								iNumberOfkMersInInput += readFasta_partialSort(*fast_q_a_File, vInputMap, vReadNameAndLength, iSoftMaxSizeOfInputVecs, iAmountOfSpecies, bSpaced, iFileLength, overallFileSize, bReadIDsAreInteresting, transferBetweenRuns, T, vNumberOfGarbagekMersPerK);
							}

							iNumOfReads = transferBetweenRuns->vReadIDs.size();
						}
						else {
							if (isGzipped) {
								iNumberOfkMersInInput += readFastq_partialSort(*fast_q_a_File_gz, vInputMap, vReadNameAndLength, iSoftMaxSizeOfInputVecs, iAmountOfSpecies, bSpaced, iFileLength, iFileLength, bReadIDsAreInteresting, transferBetweenRuns, T, vNumberOfGarbagekMersPerK);
							}
							else {
								iNumberOfkMersInInput += readFastq_partialSort(*fast_q_a_File, vInputMap, vReadNameAndLength, iSoftMaxSizeOfInputVecs, iAmountOfSpecies, bSpaced, iFileLength, overallFileSize, bReadIDsAreInteresting, transferBetweenRuns, T, vNumberOfGarbagekMersPerK);
							}

							iNumOfReads = transferBetweenRuns->vReadIDs.size();
						}

						// sort suffixes for each range in parallel
						vector<pair<uint64_t, Utilities::rangeContainer>> vInputVec;

						//vInputVec.reserve(vInputMap.size());
						for (auto it = vInputMap.begin(); it != vInputMap.end();) {
							vInputVec.push_back(*it);
							vInputMap.erase(it++);
						}
						vInputMap.clear();
# if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
						__gnu_parallel::sort(vInputVec.begin(), vInputVec.end(), [](const pair<uint64_t, Utilities::rangeContainer>& p1, const pair<uint64_t, Utilities::rangeContainer>& p2) { return p1.first < p2.first; }, __gnu_parallel::balanced_quicksort_tag());
#else					
#if __has_include(<execution>)
						sort(std::execution::par_unseq, vInputVec.begin(), vInputVec.end(), [](const pair<uint64_t, Utilities::rangeContainer>& p1, const pair<uint64_t, Utilities::rangeContainer>& p2) { return p1.first < p2.first; });
#else
						sort(vInputVec.begin(), vInputVec.end(), [](const pair<uint64_t, Utilities::rangeContainer>& p1, const pair<uint64_t, Utilities::rangeContainer>& p2) { return p1.first < p2.first; });
#endif
#endif
						// sort inside each vector (in parallel)
						auto sortingFunction = [this, &vInputVec, &bUnique](const int32_t iThreadID) {
							size_t iParallelCount = 0;

							for (auto it = vInputVec.begin(); it != vInputVec.end(); ++it, ++iParallelCount) {
								if (static_cast<int32_t>(iParallelCount%_iNumOfThreads) != iThreadID) {
									continue;
								}
								if (_iMinK <= 6) {
									sort(it->second.kMers_ST6.begin(), it->second.kMers_ST6.end(), [](const pair<uint64_t, readIDType>& a, const pair<uint64_t, readIDType>& b) { return a < b; });
									if (bUnique) {
										const auto& newIt = unique(it->second.kMers_ST6.begin(), it->second.kMers_ST6.end(), [](const pair<uint64_t, readIDType>& a, const pair<uint64_t, readIDType>& b) { return a == b; });
										it->second.kMers_ST6.resize(newIt - it->second.kMers_ST6.begin());
									}
								}
								else {
									sort(it->second.kMers_GT6.begin(), it->second.kMers_GT6.end(), [](const pair<uint32_t, readIDType>& a, const pair<uint32_t, readIDType>& b) { return a < b; });
									if (bUnique) {
										const auto& newIt = unique(it->second.kMers_GT6.begin(), it->second.kMers_GT6.end(), [](const pair<uint32_t, readIDType>& a, const pair<uint32_t, readIDType>& b) { return a == b; });
										it->second.kMers_GT6.resize(newIt - it->second.kMers_GT6.begin());
									}
								}
							}
						};

						// enqueue and run
						for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
							workerThreadPool[iThreadID].pushTask(bind(sortingFunction, iThreadID));
							workerThreadPool[iThreadID].startThread();
						}

						for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
							workerThreadPool[iThreadID].waitUntilFinished();
						}

						auto end = std::chrono::high_resolution_clock::now();
						iTimeFastq += chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();


						if (bReadIDsAreInteresting) {
							// This holds the mapping read ID -> Genus ID
							if (iNumOfReads <= iNumOfReadsOld) {
								for (uint64_t i = 0; i < iNumOfReads*iAmountOfSpecies; ++i) {
									vReadIDtoTaxID[i] = 0.f;
								}
							}
							else {
								vReadIDtoTaxID.reset();
								vReadIDtoTaxID.reset(new float[iNumOfReads*iAmountOfSpecies]);

								for (uint64_t i = 0; i < iNumOfReads*iAmountOfSpecies; ++i) {
									vReadIDtoTaxID[i] = 0.f;
								}
							}
						}
						else {
							iNumOfReadsSum += transferBetweenRuns->iNumOfNewReads;
						}

						function<void(const int32_t&, const pair<uint64_t, Utilities::rangeContainer>*)> foo;

						// now compare with index
						start = std::chrono::high_resolution_clock::now();

						if (bReadIDsAreInteresting) {
							if (bRAM) {
								if (bPartitioned) {
									foo = bind(&Compare::compareWithDatabase< vector<packedPair>*>, this, placeholders::_1, placeholders::_2, &vLib_RAM_Half, ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(mTaxToIdx), ref(transferBetweenRuns->mReadIDToArrayIdx));
								}
								else {
									if (bUnfunny) {
										foo = bind(&Compare::compareWithDatabase_sloppy< vector<uint16_t>*>, this, placeholders::_1, placeholders::_2, &vLib_RAM_taxaOnly, ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(transferBetweenRuns->mReadIDToArrayIdx));
									}
									else {
										foo = bind(&Compare::compareWithDatabase< vector<packedBigPair>*>, this, placeholders::_1, placeholders::_2, &vLib_RAM_Full, ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(mTaxToIdx), ref(transferBetweenRuns->mReadIDToArrayIdx));
									}
								}
							}
							else {
								if (bPartitioned) {
									foo = bind(&Compare::compareWithDatabase< unique_ptr<unique_ptr<const index_t_p>[]>>, this, placeholders::_1, placeholders::_2, ref(vLibParted_p), ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(mTaxToIdx), ref(transferBetweenRuns->mReadIDToArrayIdx));
								}
								else {
									if (bUnfunny) {
										foo = bind(&Compare::compareWithDatabase_sloppy<unique_ptr<unique_ptr<const taxaOnly>[]>>, this, placeholders::_1, placeholders::_2, ref(vLib_taxaOnly), ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(transferBetweenRuns->mReadIDToArrayIdx));
									}
									else {
										foo = bind(&Compare::compareWithDatabase<unique_ptr<unique_ptr<const contentVecType_32p>[]>>, this, placeholders::_1, placeholders::_2, ref(vLib), ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(mTaxToIdx), ref(transferBetweenRuns->mReadIDToArrayIdx));
									}
								}
							}
						}
						else {
							if (bRAM) {
								if (bPartitioned) {
									foo = bind(&Compare::createProfile<vector<packedPair>*>, this, placeholders::_1, placeholders::_2, &vLib_RAM_Half, ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies), ref(mTaxToIdx));
								}
								else {
									if (bUnfunny) {
										foo = bind(&Compare::createProfile_sloppy<vector<uint16_t>*>, this, placeholders::_1, placeholders::_2, &vLib_RAM_taxaOnly, ref(vCount_all), ref(vCount_unique));
									}
									else {
										foo = bind(&Compare::createProfile< vector<packedBigPair>* >, this, placeholders::_1, placeholders::_2, &vLib_RAM_Full, ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies), ref(mTaxToIdx));
									}
								}
							}
							else {
								if (bPartitioned) {
									foo = bind(&Compare::createProfile< unique_ptr<unique_ptr<const index_t_p>[]> >, this, placeholders::_1, placeholders::_2, ref(vLibParted_p), ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies), ref(mTaxToIdx));
								}
								else {
									if (bUnfunny) {
										foo = bind(&Compare::createProfile_sloppy< unique_ptr<unique_ptr<const taxaOnly>[]> >, this, placeholders::_1, placeholders::_2, ref(vLib_taxaOnly), ref(vCount_all), ref(vCount_unique));
									}
									else {
										foo = bind(&Compare::createProfile< unique_ptr<unique_ptr<const contentVecType_32p>[]> >, this, placeholders::_1, placeholders::_2, ref(vLib), ref(vCount_all), ref(vCount_unique), ref(iAmountOfSpecies), ref(mTaxToIdx));
									}
								}

							}
						}

						// because the stxxl is not threadsafe (as in two threads cannot access different locations on the drive), we need to separate the work
						uint64_t iSumOfRanges = 0;
						for (const auto& entry : vInputVec) {
							iSumOfRanges += (entry.second.range + 1)*( (_iMinK >= 6) ? entry.second.kMers_GT6.size() : entry.second.kMers_ST6.size() );
						}
						size_t iDiv = iSumOfRanges / _iNumOfThreads;
						auto inputIt = vInputVec.cbegin();

						for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {

							workerThreadPool[iThreadID].setNumberOfTasks(vInputVec.size() / static_cast<size_t>(_iNumOfThreads) + 1);

							uint64_t iCurrentSum = 0;
							while (iCurrentSum < iDiv && inputIt != vInputVec.cend()) {
								auto task = bind(foo, iThreadID, &(*inputIt));
								workerThreadPool[iThreadID].pushTask(task);
								iCurrentSum += (inputIt->second.range + 1)*((_iMinK >= 6) ? inputIt->second.kMers_GT6.size() : inputIt->second.kMers_ST6.size());
								++inputIt;
							}
						}

						for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
							workerThreadPool[iThreadID].startThread();
						}
						for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
							workerThreadPool[iThreadID].waitUntilFinished();
						}
					
						if (someThingWentWrong) {
							rethrow_exception(someThingWentWrong);
						}

						end = std::chrono::high_resolution_clock::now();
						iTimeCompare += chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
						
						////////////////////////////////////////////////////////////////////////////////////////////////
						// save results
						if (bReadIDsAreInteresting) {
							if (transferBetweenRuns->addTail) {
								// last read is not yet finished
								
								uint64_t i = 0;
								// check if there is still some unfinished read which is now complete
								if (vSavedScores.size()) {
									if (transferBetweenRuns->lastLine.second != iReadIDofSavedScores) {
										auto lastScoreVec = getVecOfScored(vReadIDtoTaxID, Utilities::checkIfInMap(transferBetweenRuns->mReadIDToArrayIdx, iReadIDofSavedScores)->second);
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
										scoringFunc(move(vSavedScores), (i++) + iNumOfReadsSum, vReadNameAndLength.front(), mFrequencies, mIdxToTax, mOrganisms, move(fOut));
										vSavedScores.clear();

										transferBetweenRuns->vReadIDs.erase(find(transferBetweenRuns->vReadIDs.begin(), transferBetweenRuns->vReadIDs.end(), iReadIDofSavedScores));
										transferBetweenRuns->mReadIDToArrayIdx.erase(iReadIDofSavedScores);
										vReadNameAndLength.pop_front();
									}
								}
								
								// save the score of the not yet finished
								auto resultOfUnfinished = getVecOfScored(vReadIDtoTaxID, Utilities::checkIfInMap(transferBetweenRuns->mReadIDToArrayIdx, transferBetweenRuns->lastLine.second)->second);
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
								scoringFunc(vReadIDtoTaxID, i, iNumOfReads - 1, iNumOfReadsSum, vReadNameAndLength, mFrequencies, mIdxToTax, mOrganisms, move(fOut));

								for (; i < iNumOfReads - 1; ++i) {
									auto tempID = transferBetweenRuns->vReadIDs.front();
									transferBetweenRuns->vReadIDs.pop_front();
									transferBetweenRuns->mReadIDToArrayIdx.erase(tempID);
									//vReadNameAndLength.pop_front();
								}
								

								iReadIDofSavedScores = transferBetweenRuns->iCurrentReadID;
								transferBetweenRuns->mReadIDToArrayIdx[iReadIDofSavedScores] = 0;

								iNumOfReadsSum += iNumOfReads - 1;
								iNumOfReadsOld = (iNumOfReads - 1 < iNumOfReadsOld) ? iNumOfReadsOld : iNumOfReads - 1;
							}
							else {
								// reads finished
								
								// write down the saved one
								uint64_t i = 0;
								if (vSavedScores.size()) {
									auto resultOfFinished = getVecOfScored(vReadIDtoTaxID, 0);
									vSavedScores.insert(vSavedScores.end(), resultOfFinished.cbegin(), resultOfFinished.cend());
									resultOfFinished.clear();
									sort(vSavedScores.begin(), vSavedScores.end(), [](const tuple<uint32_t, float, double>& p1, const tuple<uint32_t, float, double>& p2) { return get<0>(p1) < get<0>(p2); });
									auto seen = vSavedScores[0];
									for (auto it = vSavedScores.begin() + 1; it != vSavedScores.end(); ++it) {
										if (get<0>(*it) != get<0>(seen)) {
											resultOfFinished.push_back(seen);
											seen = *it;
										}
										else {
											get<1>(seen) += get<1>(*it);
										}
									}
									resultOfFinished.push_back(seen);
									vSavedScores.swap(resultOfFinished);

									scoringFunc(move(vSavedScores), iNumOfReadsSum, vReadNameAndLength.front(), mFrequencies, mIdxToTax, mOrganisms, move(fOut));
									i = 1;
									transferBetweenRuns->vReadIDs.erase(find(transferBetweenRuns->vReadIDs.begin(), transferBetweenRuns->vReadIDs.end(), iReadIDofSavedScores));
									transferBetweenRuns->mReadIDToArrayIdx.erase(iReadIDofSavedScores);
									vReadNameAndLength.pop_front();
									vSavedScores.clear();
								}
								
								// and now the regular ones
								scoringFunc(vReadIDtoTaxID, i, iNumOfReads, iNumOfReadsSum, vReadNameAndLength, mFrequencies, mIdxToTax, mOrganisms, move(fOut));
								for (; i < iNumOfReads; ++i) {
									auto tempID = transferBetweenRuns->vReadIDs.front();
									transferBetweenRuns->vReadIDs.pop_front();
									transferBetweenRuns->mReadIDToArrayIdx.erase(tempID);
									//vReadNameAndLength.pop_front();
								}
								
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

						//iterate until no dna is left
					}

					///////////////////////////////////////////////////////////////////////////////////////////////

					// if json is the output format for readToTaxa, end it with a ]
					if (bReadIDsAreInteresting && !bHumanReadable) {
						fOut << "\n" << "]";
					}
					fOut.flush(); // empty the buffer to avoid memory leak

					// sum up parallel results
					for (int32_t iThreadID = 1; iThreadID < _iNumOfThreads; ++iThreadID) {
						const uint64_t& iStepsize = iThreadID * _iNumOfK*iAmountOfSpecies;
						for (uint64_t iIdx = 0; iIdx < _iNumOfK*iAmountOfSpecies; ++iIdx) {
							vCount_all[iIdx] += vCount_all[iStepsize + iIdx];
							vCount_unique[iIdx] += vCount_unique[iStepsize + iIdx];
						}
					}

					// get profiling results
					vector<uint64_t> vSumOfUniquekMers(_iNumOfK);
					vector<double> vSumOfNonUniques(_iNumOfK);
					vector<tuple<string, vector<pair<double, uint64_t>>, uint32_t>> vOut(iAmountOfSpecies, tuple<string, vector<pair<double, uint64_t>>, uint32_t>("", vector<pair<double, uint64_t>>(_iNumOfK), 0));
					for (uint32_t iSpecIdx = 1; iSpecIdx < iAmountOfSpecies; ++iSpecIdx) {
						vector<pair<double, uint64_t>> vTemp(_iNumOfK);
						for (int32_t ikMerlength = 0; ikMerlength < _iNumOfK; ++ikMerlength) {
							const uint64_t& iTempScore = vCount_unique[iSpecIdx*uint64_t(_iNumOfK) + ikMerlength];
							vSumOfUniquekMers[ikMerlength] += iTempScore;
							vSumOfNonUniques[ikMerlength] += vCount_all[iSpecIdx*uint64_t(_iNumOfK) + ikMerlength];
							vTemp[ikMerlength] = make_pair(vCount_all[iSpecIdx*uint64_t(_iNumOfK) + ikMerlength], iTempScore);
						}
						vOut[iSpecIdx] = make_tuple(mOrganisms[iSpecIdx], vTemp, mIdxToTax[iSpecIdx]);
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
					ofstream tableFileStream;
					//tableFileStream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
					tableFileStream.open((vInputFiles.size() > 1) ? fTableFile + fileName + ".csv" : fTableFile);
					//auto orgBuf = cout.rdbuf();
					if (fTableFile != "") {
						//cout.rdbuf(tableFileStream.rdbuf());

						if (bHumanReadable) {
							// short version: taxID,Name,Unique Percentage of highest k,Non-unique Percentage of highest k\n
							bool bBreakOut = false;
							double iSumOfIdentified = 0;
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
									iSumOfIdentified += get<1>(entry)[0].first;
									sOutStr << "," << static_cast<double>(get<1>(entry)[0].first) / (iNumberOfkMersInInput - vNumberOfGarbagekMersPerK[0]) * 100.;
									sOutStr << "\n";
								}
							}

							// last entry
							tableFileStream << "0,not identified,"
											<< "0.0,0.0,"
											<< ((static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[0]) - iSumOfIdentified) / (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[0]))) * 100.;
							tableFileStream << "\n" << sOutStr.str();
						}
						else {
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
							tableFileStream << "\n";

							vector<double> iSumOfIdentified(_iNumOfK, 0);
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
										iSumOfIdentified[ikMerlength] += get<1>(entry)[ikMerlength].first;
										sOutStr << "," << static_cast<double>(get<1>(entry)[ikMerlength].first) / (iNumberOfkMersInInput*(ikMerlength + 1) - vNumberOfGarbagekMersPerK[ikMerlength]);
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
								//cout << iNumberOfkMersInInput * (ikMerlength + 1) << " " << vNumberOfGarbagekMersPerK[ikMerlength] << " " << iSumOfIdentified[ikMerlength] << endl;
								tableFileStream << "," << (static_cast<double>(iNumberOfkMersInInput*(ikMerlength + 1)) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]) - static_cast<double>(iSumOfIdentified[ikMerlength])) / (static_cast<double>(iNumberOfkMersInInput*(ikMerlength + 1)) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]));
							}
							tableFileStream << "\n" << sOutStr.str();
						}
					}
					/*if (fTableFile != "") {
						cout.rdbuf(orgBuf);
					}*/
					if (_bVerbose) {
						cout << "OUT: Number of k-mers in input: " << iNumberOfkMersInInput << " of which " << vSumOfNonUniques[0]/iNumberOfkMersInInput * 100. << " % were identified." << endl;
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