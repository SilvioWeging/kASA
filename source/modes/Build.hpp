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
#include "../utils/ParallelQuicksort.hpp"

namespace kASA {
	template<class vecType, class elemType>
	class Build {
	private:
		unique_ptr<vector<elemType>> vInternal;
		//uint64_t iInternalCounter = 0;
		vector<uint64_t> vVectorSizes;

		string _sTempPath = string("");
		//uint64_t _iConstSize = 0;
		int32_t _iCounterOfContainers = 0;

		//size_t _iAmountOfSpace = 1024000;
		size_t _iSoftSize = 0;

#if __GNUC__ || defined(__llvm__)
		int32_t _iNumOfThreads_ = 1;
#endif

		//unique_ptr<uint64_t[]> arrFrequencies;
		//unordered_map<uint32_t, uint32_t> _mContent;

	public:

		Build(const string& path, const int32_t& iNumOfCall, const int32_t& iNumOfThreads, const size_t& iSoftLimit, const uint64_t&) {
			//ofstream derp;
			//derp.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
			//derp.open(_sTempPath + to_string(_iFlagOfContainerIdx));
			//derp.close();
			//derp.open(_sTempPath + to_string(!_iFlagOfContainerIdx));

			//_iAmountOfSpace = iSoftLimit;

			/*arrFrequencies.reset(new uint64_t[iNumOfTaxa* _iHighestK]);
			for (uint64_t i = 0; i < iNumOfTaxa * _iHighestK; ++i) {
				arrFrequencies[i] = 0;
			}*/
			setMemberVariables(path, iNumOfCall, iNumOfThreads, iSoftLimit);
		}

		Build() {}

		/////////////////////////////////////////////////////////////////////////////////////////
	public:

		/////////////////////////////////////////////////////////////////////////////////////////
		// a posteriori constructor
		inline void setMemberVariables(const string& path, const int32_t& iNumOfCall, const int32_t& iNumOfThreads, const size_t& iSoftLimit) {
			_sTempPath = path + "_temp_" + to_string(iNumOfCall) + "_";

			try {
				vInternal.reset(new vector<elemType>());
				vInternal->reserve(iSoftLimit + 1);
			}
			catch (const bad_alloc&) {
				cerr << "ERROR: Not enough memory provided. Quitting now..." << endl;
				throw;
			}

			_iSoftSize = iSoftLimit;

#if __GNUC__ || defined(__llvm__)
			_iNumOfThreads_ = iNumOfThreads;
#endif
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// For the alternative mode, where the size of the internal vector might change
		inline void setInternalSize(const size_t& iSoftLimit) {
			try {
				vInternal.reset(new vector<elemType>());
				vInternal->reserve(iSoftLimit + 1);
			}
			catch (const bad_alloc&) {
				cerr << "ERROR: Not enough memory provided. Quitting now..." << endl;
				throw;
			}

			_iSoftSize = iSoftLimit;
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// If someone wants to use this directly
		inline vector<elemType>* getInternal() {
			return &vInternal;
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// If the merging should be continued, the sizes must be recalculated
		inline vector<uint64_t>* getVectorSizesVec() {
			return &vVectorSizes;
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// Set num of containers 
		inline void setNumOfContainers(const size_t& numberOfContainers) {
			_iCounterOfContainers = static_cast<int32_t>(numberOfContainers);
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// add elements to internal vector until it is full
		inline bool addToInt(const tuple<uint64_t,uint32_t>& elem) {
			vInternal->push_back(elemType(get<0>(elem), get<1>(elem)));
			if (vInternal->size() >= _iSoftSize) {
				return false;
			}
			return true;
		}

		inline bool addToInt(const tuple<uint128_t, uint32_t>& elem) {
			vInternal->push_back(elemType(get<0>(elem), get<1>(elem)));
			if (vInternal->size() >= _iSoftSize) {
				return false;
			}
			return true;
		}

		inline void visualize(const vecType* vC, typename vecType::const_iterator& endIt) {
			ofstream derp(_sTempPath + "derp.txt");
			for (auto it = vC->cbegin(); it != endIt; ++it) {
				derp << it->first << ", " << it->second << endl;
			}
		}

		inline void visualize(typename vecType::const_iterator& bg, typename vecType::const_iterator& end, ofstream& outFile) {
			for (auto it = bg; it != end; ++it) {
				if (outFile.is_open()) {
					outFile << it->first << ", " << it->second << endl;
				}
				else {
					cout << it->first << ", " << it->second << endl;
				}
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// It's like merging two sorted decks of cards into one
		template<typename W, typename R, typename I>
		inline uint64_t merge(W& vNCIt, R& it, const uint64_t& vecSize, I&& vIntItC, const I&& vIntEndItC, const pair<unordered_map<uint32_t, uint32_t>, unordered_map<uint32_t, uint32_t>>& mapsForDummys) {
			try {
				elemType tSeenInt; //eliminate duplicates
				bool bIndexIntChanged = false;

				uint64_t ivCSizeCounter = 0, iElemCounter = 0;

				while (ivCSizeCounter < vecSize && vIntItC != vIntEndItC) {

					//cout << it->first << "," << get<1>(*it) << "," << bIndexCChanged << " " << get<0>(*vIntItC) << "," << get<1>(*vIntItC) << "," << bIndexIntChanged << endl;

					if (*vIntItC == tSeenInt && bIndexIntChanged) {
						++vIntItC;
						continue;
					}
					else {
						tSeenInt = *vIntItC;
						bIndexIntChanged = false;
					}

					if (it->first < vIntItC->first) {
						elemType pair = *it;
						if (mapsForDummys.first.size()) {
							auto res = mapsForDummys.first.find(pair.second);
							if (res != mapsForDummys.first.end()) {
								pair.second = res->second;
							}
						}
						vNCIt << pair;
						++it;
						++ivCSizeCounter;
						++iElemCounter;
						bIndexIntChanged = false;
					}
					else {
						if (it->first == vIntItC->first) {
							if (it->second < vIntItC->second) {
								elemType pair = *it;
								if (mapsForDummys.first.size()) {
									auto res = mapsForDummys.first.find(pair.second);
									if (res != mapsForDummys.first.end()) {
										pair.second = res->second;
									}
								}
								vNCIt << pair;
								++it;
								++ivCSizeCounter;
								++iElemCounter;
								bIndexIntChanged = false;
							}
							else {
								if (it->second > vIntItC->second) {
									elemType pair = *vIntItC;
									if (mapsForDummys.second.size()) {
										auto res = mapsForDummys.second.find(pair.second);
										if (res != mapsForDummys.second.end()) {
											pair.second = res->second;
										}
									}
									vNCIt << pair;
									++vIntItC;
									++iElemCounter;
									bIndexIntChanged = true;
								}
								else {
									elemType pair = *vIntItC;
									if (mapsForDummys.second.size()) {
										auto res = mapsForDummys.second.find(pair.second);
										if (res != mapsForDummys.second.end()) {
											pair.second = res->second;
										}
									}
									vNCIt << pair;
									++it;
									++ivCSizeCounter;
									++vIntItC;
									++iElemCounter;
									bIndexIntChanged = true;
								}
							}
						}
						else {
							elemType pair = *vIntItC;
							if (mapsForDummys.second.size()) {
								auto res = mapsForDummys.second.find(pair.second);
								if (res != mapsForDummys.second.end()) {
									pair.second = res->second;
								}
							}
							vNCIt << pair;
							++vIntItC;
							++iElemCounter;
							bIndexIntChanged = true;
						}
					}
				}

				// Send the rest into the other vector respectively
				//while (it != vC->cend()) {
				while (ivCSizeCounter < vecSize) {
					elemType pair = *it;
					if (mapsForDummys.first.size()) {
						auto res = mapsForDummys.first.find(pair.second);
						if (res != mapsForDummys.first.end()) {
							pair.second = res->second;
						}
					}
					vNCIt << pair;
					++it;
					++iElemCounter;
					++ivCSizeCounter;
				}

				while (vIntItC != vIntEndItC) {
					if (*vIntItC == tSeenInt && bIndexIntChanged) {
						++vIntItC;
						continue;
					}
					else {
						tSeenInt = *vIntItC;
						bIndexIntChanged = false;
					}

					elemType pair = *vIntItC;
					if (mapsForDummys.second.size()) {
						auto res = mapsForDummys.second.find(pair.second);
						if (res != mapsForDummys.second.end()) {
							pair.second = res->second;
						}
					}
					vNCIt << pair;
					++vIntItC;
					++iElemCounter;
				}
				return iElemCounter;
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
			return 0;
		}


		/////////////////////////////////////////////////////////////////////////////////////////
		// if the internal vector is full enough, save it to an external one
		inline void IntToExtPart() {
			try {

# if __GNUC__ || defined(__llvm__)
				Utilities::parallelQuicksort(vInternal->begin(), vInternal->end(), [](const elemType& a, const elemType& b) { return a < b; }, _iNumOfThreads_);
#else
#if __has_include(<execution>)
				sort(std::execution::par_unseq, vInternal->begin(), vInternal->end());
#else
				sort(vInternal->begin(), vInternal->end());
#endif
#endif
				

				// squeeze it a little
				auto newEnd = std::unique(vInternal->begin(), vInternal->end(), [](elemType& a, elemType& b) {
					if (a == b) {
						return true;
					}
					return false;
				});
				vInternal->resize(newEnd - vInternal->begin());
				// if sorting would be faster, this would make sense since you'd want to maximize the amount of data inside each temporary file but alas, this takes too long
				//if (float(newEnd - vInternal->begin())/_iSoftSize < 0.95 && !bLastCall) {
				//	return;
				//}
				// else:

				// save to external vector
				Utilities::checkIfFileCanBeCreated(_sTempPath + to_string(_iCounterOfContainers));

				unique_ptr<stxxlFile> fNCFile(new stxxlFile(_sTempPath + to_string(_iCounterOfContainers), stxxl::file::RDWR));
				unique_ptr<vecType> vNC(new vecType(fNCFile.get(),vInternal->size()));
				vVectorSizes.push_back(vInternal->size());
				typename vecType::bufwriter_type vNCIt(*vNC);

				auto vIntEndItC = vInternal->cend();
				for (auto it = vInternal->cbegin(); it != vIntEndItC; ++it) {
					vNCIt << *(it);
				}

				vNCIt.finish();

				vNC->export_files("_");
				_iCounterOfContainers++;

				vInternal->clear();
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// Merge as many external vectors as possible in one call until only one remains
		inline uint64_t mergeTemporaries(const string& sOutPath) {
			try {
				vInternal.reset(); //not needed anymore

				auto remainingFiles = _iCounterOfContainers;
				int32_t iCurrentIdx = 0;
				auto comparisonFunction = [](const pair<typename vecType::const_iterator, uint32_t>& it1, const pair<typename vecType::const_iterator, uint32_t>& it2) { return *(it1.first) < *(it2.first); };
				priority_queue<pair<typename vecType::const_iterator,uint32_t>, vector<pair<typename vecType::const_iterator, uint32_t>>, decltype(comparisonFunction)> vCurrentIterators(comparisonFunction);

				/*auto print_queue = [](priority_queue<pair<typename vecType::const_iterator, uint32_t>, vector<pair<typename vecType::const_iterator, uint32_t>>, decltype(comparisonFunction)> Q) {
					while (!Q.empty()) {
						cout << kASA::kMerToAminoacid(Q.top().first->first, 12) << endl;
						Q.pop();
					}
					cout << endl;
				};*/

				size_t iNumOfVecsThatCanBeOpened = _iSoftSize * sizeof(elemType) / (vecType::block_size * vecType::page_size * 4); // It could be, that the memory needed to open the files is larger than the memory provided
				//cout << "iNumOfVecsThatCanBeOpened: " << iNumOfVecsThatCanBeOpened << endl;
				iNumOfVecsThatCanBeOpened = (iNumOfVecsThatCanBeOpened > FOPEN_MAX - 1) ? FOPEN_MAX - 1 : iNumOfVecsThatCanBeOpened;

				auto maxNrOfFiles = remainingFiles;
				while (remainingFiles != 1) {
					// first gather as many files as permitted
					vector<unique_ptr<stxxlFile>> currentFiles(FOPEN_MAX - 1);
					vector<unique_ptr<vecType>> currentVecs(FOPEN_MAX - 1);

					uint64_t iSumSize = 0;
					for (int32_t i = 0; i < static_cast<int32_t>(iNumOfVecsThatCanBeOpened) && iCurrentIdx < maxNrOfFiles; ++iCurrentIdx, ++i) {
						currentFiles[i].reset(new stxxlFile(_sTempPath + to_string(iCurrentIdx), stxxl::file::RDWR));
						currentVecs[i].reset(new vecType(currentFiles[i].get(), vVectorSizes[iCurrentIdx]));
						vCurrentIterators.push(make_pair(currentVecs[i]->cend() - 1, i));
						iSumSize += vVectorSizes[iCurrentIdx];
						remainingFiles--;
					}

					//print_queue(vCurrentIterators);

					// create merge-file
					Utilities::checkIfFileCanBeCreated(_sTempPath + to_string(_iCounterOfContainers));

					unique_ptr<stxxlFile> fNCFile(new stxxlFile(_sTempPath + to_string(_iCounterOfContainers), stxxl::file::RDWR));
					//unique_ptr<vecType> vNC(new vecType(fNCFile.get(), iSumSize));
					size_t sizeOfAGigabyte = 1024 * 1024 * 1024 / sizeof(elemType);
					unique_ptr<vecType> vNC(new vecType(fNCFile.get(), 0));
					vNC->reserve(sizeOfAGigabyte);
					auto addElement = [&](const uint64_t& iCurrentSize, const elemType& elem) {
						if (iCurrentSize < vNC->size()) {
							vNC->push_back(elem);
						}
						else {
							vNC->reserve(vNC->size() + sizeOfAGigabyte);
							vNC->push_back(elem);
						}
					};
					//typename vecType::bufwriter_type vNCIt(*vNC);

					elemType pSeen;
					uint64_t iSizeOfNewVec = 0;
					while (!vCurrentIterators.empty()) {
						// take smallest value and check if that has been seen before
						auto it = vCurrentIterators.top();

						if (!(*(it.first) == pSeen)) {
							// if not: write it to the output
							pSeen = *(it.first);
							//vNCIt << pSeen;
							addElement(iSizeOfNewVec, pSeen);
							++iSizeOfNewVec;
						}
						
						// The end of that vector has been reached: take this iterator (and temporary vector) out of consideration
						if (it.first == currentVecs[it.second]->cbegin()) {
							currentVecs[it.second].reset();

							currentFiles[it.second]->close_remove();
							currentFiles[it.second].reset();

							vCurrentIterators.pop();
						}
						else {
							it.first--;
							currentVecs[it.second]->resize(currentVecs[it.second]->size() - 1, true);
							
							vCurrentIterators.pop();
							vCurrentIterators.push(it);
						}

						//print_queue(vCurrentIterators);
					}
					//vNCIt.finish();
					vVectorSizes.push_back(iSizeOfNewVec);
					vNC->resize(iSizeOfNewVec, true);
					std::reverse(vNC->begin(), vNC->end());

					vNC->export_files("_");
					++_iCounterOfContainers;
					remainingFiles++;
					++maxNrOfFiles;
				}

				// move the last temporary file to the given path
				remove(sOutPath.c_str());
				// attempt to move, if that fails, copy it
				Utilities::moveFile(_sTempPath + to_string(_iCounterOfContainers - 1), sOutPath);
				remove((_sTempPath + to_string(_iCounterOfContainers - 1)).c_str());
				ofstream fLibInfo(sOutPath + "_info.txt");
				fLibInfo << vVectorSizes[_iCounterOfContainers - 1];
				if (is_same<vecType, contentVecType_128>::value) {
					fLibInfo << endl << 128;
				}

				return vVectorSizes[_iCounterOfContainers - 1];
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// BuildAll can be found in Read.hpp
	};
}