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

namespace kASA {
	class Build {
	private:
		unique_ptr<vector<packedBigPair>> vInternal;
		vector<packedBigPair>::iterator vInternalIt;
		vector<uint64_t> vVectorSizes;

		string _sTempPath = "";
		uint64_t _iConstSize = 0;
		int16_t _iCounterOfContainers = 0;

		size_t _iAmountOfSpace = 1024000, _iSoftSize = 0;

		unique_ptr<uint64_t[]> arrFrequencies;
		unordered_map<uint32_t, uint32_t> _mContent;

	public:

		Build(const string& path, const int32_t& iNumOfCall, const size_t& iSoftLimit, const uint64_t& iNumOfTaxa, const unordered_map<uint32_t, uint32_t>& mContent) : _iSoftSize(iSoftLimit), _mContent(mContent) {
			_sTempPath = path + "_temp_" + to_string(iNumOfCall) + "_";
			//ofstream derp;
			//derp.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
			//derp.open(_sTempPath + to_string(_iFlagOfContainerIdx));
			//derp.close();
			//derp.open(_sTempPath + to_string(!_iFlagOfContainerIdx));

			vInternal.reset(new vector<packedBigPair>(iSoftLimit + 1));
			vInternalIt = vInternal->begin();

			_iAmountOfSpace = iSoftLimit;

			arrFrequencies.reset(new uint64_t[iNumOfTaxa*12]);
			for (uint64_t i = 0; i < iNumOfTaxa * 12; ++i) {
				arrFrequencies[i] = 0;
			}
		}

		Build() {}

		/////////////////////////////////////////////////////////////////////////////////////////
	public:
		inline bool addToInt(const tuple<uint64_t,uint32_t>& elem) {
			vInternalIt->first = get<0>(elem);
			vInternalIt->second = get<1>(elem);
			vInternalIt++;
			if (size_t(vInternalIt - vInternal->begin()) >= _iSoftSize) {
				return false;
			}
			return true;
		}

		inline void visualize(const contentVecType_32p* vC, contentVecType_32p::const_iterator& endIt) {
			ofstream derp(_sTempPath + "derp.txt");
			for (auto it = vC->cbegin(); it != endIt; ++it) {
				derp << it->first << ", " << it->second << endl;
			}
		}

		inline void visualize(contentVecType_32p::const_iterator& bg, contentVecType_32p::const_iterator& end, ofstream& outFile) {
			for (auto it = bg; it != end; ++it) {
				if (outFile.is_open()) {
					outFile << it->first << ", " << it->second << endl;
				}
				else {
					cout << it->first << ", " << it->second << endl;
				}
			}
		}

	protected:
		/////////////////////////////////////////////////////////////////////////////////////////
		// It's like merging two sorted decks of cards into one
		template<typename W, typename R, typename I>
		inline uint64_t merge(W& vNCIt, R& it, const uint64_t& vecSize, I&& vIntItC, const I&& vIntEndItC, unique_ptr<uint64_t[]>& freqArray, const unordered_map<uint32_t, uint32_t>& mContent, const bool& bLastCall = false) {
			try {
				packedBigPair tSeenInt; //eliminate duplicates
				bool bIndexIntChanged = false;

				auto countFreqs = [&mContent,&freqArray,&bLastCall](const packedBigPair& pair) {
					if (bLastCall) {
						const auto& idx = Utilities::checkIfInMap(mContent, pair.second)->second * 12;
						for (uint8_t k = 0; k < 12; ++k) {
							if (((pair.first >> 5 * k) & 31) != 30) {
								freqArray[idx + k]++;
							}
						}
					}
				};

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
						//*(vNCIt++) = *it;
						vNCIt << *(it);
						countFreqs(*it);
						++it;
						++ivCSizeCounter;
						++iElemCounter;
						bIndexIntChanged = false;
					}
					else {
						if (it->first == vIntItC->first) {
							if (it->second < vIntItC->second) {
								//*(vNCIt++) = *it;
								vNCIt << *(it);
								countFreqs(*it);
								++it;
								++ivCSizeCounter;
								++iElemCounter;
								bIndexIntChanged = false;
							}
							else {
								if (it->second > vIntItC->second) {
									//*(vNCIt++) = *vIntItC;
									vNCIt << *(vIntItC);
									countFreqs(*vIntItC);
									++vIntItC;
									++iElemCounter;
									bIndexIntChanged = true;
								}
								else {
									//*(vNCIt++) = *vIntItC;
									vNCIt << *(vIntItC);
									countFreqs(*vIntItC);
									++it;
									++ivCSizeCounter;
									++vIntItC;
									++iElemCounter;
									bIndexIntChanged = true;
								}
							}
						}
						else {
							//*(vNCIt++) = *vIntItC;
							vNCIt << *(vIntItC);
							countFreqs(*vIntItC); 
							++vIntItC;
							++iElemCounter;
							bIndexIntChanged = true;
						}
					}
				}

				// Send the rest into the other vector respectively
				//while (it != vC->cend()) {
				while (ivCSizeCounter < vecSize) {
					//*(vNCIt++) = *it;
					vNCIt << *(it);
					countFreqs(*it);
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

					//*(vNCIt++) = *vIntItC;
					vNCIt << *(vIntItC);
					countFreqs(*vIntItC);
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

	public:

		/////////////////////////////////////////////////////////////////////////////////////////
		// if the internal vector is full enough, save it to an external one
		inline void IntToExtPart() {
			try {

#if __has_include(<execution>)
				sort(std::execution::par_unseq, vInternal->begin(), vInternalIt);
#else
# if __GNUC__
				__gnu_parallel::sort(vInternal->begin(), vInternalIt, __gnu_parallel::balanced_quicksort_tag()); 
#else
				sort(vInternal->begin(), vInternalIt);
#endif
#endif
				

				// squeeze it a little
				auto newEnd = std::unique(vInternal->begin(), vInternalIt, [](packedBigPair& a, packedBigPair& b) {
					if (a == b) {
						return true;
					}
					return false;
				});
				vInternalIt = newEnd;
				// if sorting would be faster, this would make sense since you'd want to maximize the amount of data inside each temporary file but alas, this takes too long
				//if (float(newEnd - vInternal->begin())/_iSoftSize < 0.95 && !bLastCall) {
				//	return;
				//}
				// else:

				// save to external vector
				Utilities::createFile(_sTempPath + to_string(_iCounterOfContainers));

				unique_ptr<stxxlFile> fNCFile(new stxxlFile(_sTempPath + to_string(_iCounterOfContainers), stxxl::file::RDWR));
				unique_ptr<contentVecType_32p> vNC(new contentVecType_32p(fNCFile.get(), vInternalIt - vInternal->cbegin()));
				vVectorSizes.push_back(vInternalIt - vInternal->cbegin());
				contentVecType_32p::bufwriter_type vNCIt(*vNC);

				auto vIntEndItC = static_cast<vector<packedBigPair>::const_iterator>(vInternalIt);
				for (auto it = vInternal->cbegin(); it != vIntEndItC; ++it) {
					vNCIt << *(it);
				}

				vNCIt.finish();

				vNC->export_files("_");
				_iCounterOfContainers++;

				vInternalIt = vInternal->begin();
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
				int16_t iCurrentIdx = 0;
				vector<contentVecType_32p::const_iterator> vCurrentIterators;
				

				while (remainingFiles != 1) {
					// first gather as many files as permitted
					vector<unique_ptr<stxxlFile>> currentFiles(FOPEN_MAX - 1);
					vector<unique_ptr<contentVecType_32p>> currentVecs(FOPEN_MAX - 1);
					auto maxNrOfFiles = remainingFiles;
					uint64_t iSumSize = 0;
					for (int16_t i=0; i < (FOPEN_MAX - 1) && iCurrentIdx < maxNrOfFiles; ++iCurrentIdx, ++i) {
						currentFiles[i].reset(new stxxlFile(_sTempPath + to_string(iCurrentIdx), stxxl::file::RDONLY));
						currentVecs[i].reset(new contentVecType_32p(currentFiles[i].get(), vVectorSizes[iCurrentIdx]));
						vCurrentIterators.push_back(currentVecs[i]->cbegin());
						iSumSize += vVectorSizes[iCurrentIdx];
						remainingFiles--;
					}

					// create merge-file
					Utilities::createFile(_sTempPath + to_string(_iCounterOfContainers));

					unique_ptr<stxxlFile> fNCFile(new stxxlFile(_sTempPath + to_string(_iCounterOfContainers), stxxl::file::RDWR));
					unique_ptr<contentVecType_32p> vNC(new contentVecType_32p(fNCFile.get(), iSumSize));
					contentVecType_32p::bufwriter_type vNCIt(*vNC);

					packedBigPair pSeen;
					uint64_t iSizeOfNewVec = 0;
					while (!vCurrentIterators.empty()) {
						// take smallest value and check if that has been seen before
						auto it = std::min_element(vCurrentIterators.begin(), vCurrentIterators.end(), [](const contentVecType_32p::const_iterator& it1, const contentVecType_32p::const_iterator& it2) { return *it1 < *it2; });
						auto& it2 = *it; 
						
						if (!(*it2 == pSeen)) {
							// if not: write it to the output
							pSeen = *it2;
							vNCIt << pSeen;
							++iSizeOfNewVec;
						}
						it2++;
						
						// The end of that vector has been reached: take this iterator (and temporary vector) out of consideration
						if (it2 == currentVecs[distance(vCurrentIterators.begin(), it)]->cend()) {
							const auto& whichOne = distance(vCurrentIterators.begin(), it);
							currentVecs[whichOne].reset();
							currentVecs.erase(currentVecs.begin() + whichOne);
							currentFiles[whichOne]->close_remove();
							currentFiles[whichOne].reset();
							currentFiles.erase(currentFiles.begin() + whichOne);
							vCurrentIterators.erase(it); 
							cout << whichOne;
							//for (size_t i = 0; i < vCurrentIterators.size(); ++i) {
							//	cout << " " << currentVecs[i]->size() << endl;
							//}
						}
					}
					vNCIt.finish();
					vVectorSizes.push_back(iSizeOfNewVec);
					vNC->resize(iSizeOfNewVec, true);

					vNC->export_files("_");
					++_iCounterOfContainers;
					remainingFiles++;
				}

				// move the last temporary file to the given path
				remove(sOutPath.c_str());
				rename((_sTempPath + to_string(_iCounterOfContainers - 1)).c_str(), sOutPath.c_str());
				ofstream fLibInfo(sOutPath + "_info.txt");
				fLibInfo << vVectorSizes[_iCounterOfContainers - 1];

				return vVectorSizes[_iCounterOfContainers - 1];
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}


	};
}