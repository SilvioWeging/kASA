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

		string _sTempPath = "";
		uint64_t _iConstSize = 0;
		int8_t _iFlagOfContainerIdx = 0;

		size_t _iAmountOfSpace = 1024000, _iSoftSize = 0;

		unique_ptr<uint64_t[]> arrFrequencies;
		unordered_map<uint32_t, uint32_t> _mContent;

	public:

		Build(const string& path, const int32_t& iNumOfCall, const size_t& iSoftLimit, const uint64_t& iNumOfTaxa, const unordered_map<uint32_t, uint32_t>& mContent) : _iSoftSize(iSoftLimit), _mContent(mContent) {
			_sTempPath = path + "_temp_" + to_string(iNumOfCall) + "_";
			ofstream derp;
			//derp.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
			derp.open(_sTempPath + to_string(_iFlagOfContainerIdx));
			derp.close();
			derp.open(_sTempPath + to_string(!_iFlagOfContainerIdx));

			vInternal.reset(new vector<packedBigPair>(iSoftLimit + 1));
			vInternalIt = vInternal->begin();

			_iAmountOfSpace = iSoftLimit;

			arrFrequencies.reset(new uint64_t[iNumOfTaxa*12]);
			for (uint64_t i = 0; i < iNumOfTaxa * 12; ++i) {
				arrFrequencies[i] = 0;
			}
		}
	protected:
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

		inline void IntToExt2(const bool& bLastCall = false) {
			try {
				if (!ifstream(_sTempPath + to_string(_iFlagOfContainerIdx))) {
					ofstream derp;
					//derp.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
					derp.open(_sTempPath + to_string(_iFlagOfContainerIdx));
				}
				unique_ptr<stxxlFile> fNCFile(new stxxlFile(_sTempPath + to_string(_iFlagOfContainerIdx), stxxl::file::RDWR));
				unique_ptr<contentVecType_32p> vNC(new contentVecType_32p(fNCFile.get(), _iConstSize + (vInternalIt - vInternal->cbegin())));
				//auto vNCIt = vNC->begin();
				contentVecType_32p::bufwriter_type vNCIt(*vNC);
#if __has_include(<execution>)
				sort(std::execution::par_unseq, vInternal->begin(), vInternalIt);
#else
				sort(vInternal->begin(), vInternalIt);
#endif

				packedBigPair tSeenInt; //eliminate duplicates

				if (_iConstSize == 0) {
					auto vIntEndItC = static_cast<vector<packedBigPair>::const_iterator>(vInternalIt);
					//uint32_t iDebugCounter = 0;
					for (auto it = vInternal->cbegin(); it != vIntEndItC;) {
						/*if (get<0>(*it) == 37191153727705383 ) {
							cout << "derp" << endl;
						}*/
						if (*it == tSeenInt) {
							++it;
							continue;
						}
						else {
							tSeenInt = *it;
						}

						//*(vNCIt++) = *(it++);
						vNCIt << *(it);
						if (bLastCall) {
							const auto& idx = Utilities::checkIfInMap(_mContent, it->second)->second * 12;
							for (uint8_t k = 0; k < 12; ++k) {
								if (((it->first >> 5 * k) & 31) != 30) {
									arrFrequencies[idx + k]++;
								}
							}
						}
						++it;
						//++iDebugCounter;
						++_iConstSize;
					}
				}
				else {
					if (!ifstream(_sTempPath + to_string(!_iFlagOfContainerIdx))) {
						ofstream derp;
						//derp.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
						derp.open(_sTempPath + to_string(!_iFlagOfContainerIdx));
					}
					unique_ptr<stxxlFile> fCFile(new stxxlFile(_sTempPath + to_string(!_iFlagOfContainerIdx), stxxl::file::RDONLY));
					unique_ptr<const contentVecType_32p> vC(new const contentVecType_32p(fCFile.get(), _iConstSize));

					contentVecType_32p::bufreader_type it(*vC);

					_iConstSize = merge(vNCIt, it, vC->size(), vInternal->cbegin(), static_cast<vector<packedBigPair>::const_iterator>(vInternalIt), arrFrequencies, _mContent, bLastCall);
				}


				vNCIt.finish();
				//_iConstSize = vNCIt - vNC->begin();

				vNC->resize(_iConstSize + 1, true);
				vNC->export_files("_");
				_iFlagOfContainerIdx = !_iFlagOfContainerIdx;

				//cout << stxxl::is_sorted(vNC->cbegin(), vNC->cend()) << endl;

				vInternal->clear();
				vInternal->resize(_iAmountOfSpace);
				vInternalIt = vInternal->begin();
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////
	public:
		
		inline uint64_t createDatabase(const string& sOutPath, const uint64_t& iNumOfTaxa, const unordered_map<uint32_t, string>& mIdxToName) {
			try {
				remove(sOutPath.c_str()); // "Override"
				remove((_sTempPath + to_string(_iFlagOfContainerIdx)).c_str()); // delete the other temporary
				rename((_sTempPath + to_string(!_iFlagOfContainerIdx)).c_str(), sOutPath.c_str()); // just move the file
				ofstream fLibInfo(sOutPath + "_info.txt");
				fLibInfo << _iConstSize;

				ofstream outFile(sOutPath + "_f.txt");
				for (uint32_t j = 0; j < iNumOfTaxa; ++j) {
					outFile << Utilities::checkIfInMap(mIdxToName, j)->second << "\t";
					outFile << arrFrequencies[j * 12];
					for (int32_t k = 1; k < 12; ++k) {
						outFile << "\t" << arrFrequencies[j * 12 + k];
					}
					outFile << endl;
				}
				
				return _iConstSize;
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
	};
}