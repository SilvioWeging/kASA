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

#include "Trie.hpp"



namespace kASA {
	class Shrink : public kASA {

	public:
		Shrink(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const string& stxxl_mode = "") : kASA(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, stxxl_mode) {}

		enum ShrinkingStrategy
		{
			EveryNth,
			TrieHalf,
			Overrepresented
		};

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// count occurrence of taxIDs per kMer and locate 99% of kMers cutoff. Most kMers have only one taxon.
		inline uint32_t histogram(const unique_ptr<const contentVecType_32p>& vLib, const uint32_t& iNumOfTaxIDs) {
			uint64_t iSeenkMer = (vLib->at(0)).first;
			vector<uint32_t> taxIDs(iNumOfTaxIDs + 1);
			uint32_t iCounter = 0;
			uint64_t iUniquekMerCounter = 0;

			for (const auto& entry : *vLib) {
				const auto& kMer = entry.first;
				if (kMer == iSeenkMer) {
					iCounter++;
				}
				else {
					taxIDs[iCounter]++;
					iCounter = 1;
					iSeenkMer = kMer;
					++iUniquekMerCounter;
				}
			}
			taxIDs[iCounter]++;

			double percentage = 0.0;
			for (uint32_t i = 1; i < iNumOfTaxIDs; ++i) {
				if (percentage >= 0.99) {
					return i - 1;
				}
				percentage += static_cast<double>(taxIDs[i]) / iUniquekMerCounter;
			}
			return 4; // Some magic number
		}

	private:

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Put first six letters in Trie, the others in index and restrict number of taxa to 2^16 - 1
		inline void putHalfInTrie(const unique_ptr<const contentVecType_32p>& vLib, unique_ptr<index_t_p>& vOut, const unordered_map<uint32_t, uint32_t>& mIDsAsIdx, const string& sOutfile, const bool& bTrieIsNotThere) {
			try {
				unique_ptr<stxxlFile> trieFile;
				unique_ptr<trieVector> trieVec;
				if (bTrieIsNotThere) {
					ofstream dummyFile(sOutfile + "_trie");
					dummyFile.close();
					trieFile.reset(new stxxlFile(sOutfile + "_trie", stxxl::file::RDWR));
					trieVec.reset(new trieVector(trieFile.get(), 0));
				}
				int64_t iCount = 0;
				uint64_t iHigher6LettersBitmask = 31;
				for (uint8_t i = 1; i < 6; ++i) {
					iHigher6LettersBitmask |= 31ULL << (5 * i);
				}
				iHigher6LettersBitmask <<= 30;

				uint64_t iLower6LettersBitMask = 1073741823ULL, currentShortMer = (vLib->cbegin())->first & iHigher6LettersBitmask;
				auto outIt = vOut->begin();

				for (auto libIt = vLib->cbegin(); libIt != vLib->cend(); ++libIt) {
					// no need for if (outIt != vOut->end) since the output should at most be of the same size as the original one

					const auto& entry = *libIt;

					outIt->second = static_cast<uint16_t>(Utilities::checkIfInMap(mIDsAsIdx, entry.second)->second);
					outIt->first = static_cast<uint32_t>(entry.first & iLower6LettersBitMask);
					++outIt;
					++iCount;

					//cout << kMerToAminoacid(get<0>(entry), 12) << " " << iCount << endl;

					const uint64_t& tempMer = entry.first & iHigher6LettersBitmask;
					if (tempMer != currentShortMer && bTrieIsNotThere) {
						trieVec->push_back(packedBigPairTrie(static_cast<uint32_t>(currentShortMer >> 30), iCount - 1));
						currentShortMer = tempMer;
						iCount = 1;
					}
				}

				if (bTrieIsNotThere) {
					if (iCount != 1) {
						trieVec->push_back(packedBigPairTrie(uint32_t(currentShortMer >> 30), iCount - 1));
					}
					else {
						trieVec->push_back(packedBigPairTrie(uint32_t(currentShortMer >> 30), 1));
					}
					ofstream sizeFile(sOutfile + "_trie.txt");
					sizeFile << trieVec->size();
					trieVec->export_files("_");
				}

				vOut->resize(outIt - vOut->begin(), true);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}

		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Do the same as in "putHalfInTrie" but delete overrepresented kMers first
		/*inline void deleteOverrepresentedAndputHalfInTrie(const unique_ptr<const contentVecType_32p>& vLib, unique_ptr<index_t_p>& vOut, const unordered_map<uint32_t, uint32_t>& mIDsAsIdx, const string& sOutfile, unique_ptr<uint64_t[]>& arrFrequencies) {
			try {
				// Gather Histogram for deletion of useless k-mers
				//const auto& howManyTaxIDsAreOkay = 1000;//mIDsAsIdx.size();

				//TODO
				throw runtime_error("Not implemented yet");
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
		*/
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Count the frequencies as a helper function
		inline void countFreqs(const unordered_map<uint32_t, uint32_t>& mContent, unique_ptr<uint64_t[]>& freqArray, const packedBigPair& pair) {
			try {
				const auto& idx = Utilities::checkIfInMap(mContent, pair.second)->second * 12;
				for (uint8_t k = 0; k < 12; ++k) {
					if (((pair.first >> 5 * k) & 31) != 30) {
						freqArray[idx + k]++;
					}
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Throw out kMers for every ID where fN is the percentage of how many are thrown out of the respective ID (e.g. fN = 60 -> throw out 60%, iCount[ID] = 1200 -> 480 stay, 720 are thrown away)
		inline void deleteEveryNth(const unique_ptr<const contentVecType_32p>& vLib, unique_ptr<contentVecType_32p>& vOut, const float& fN, const unordered_map<uint32_t, uint32_t>& mContent, unique_ptr<uint64_t[]>&freqArray) {
			try {
				vector<uint64_t> vSteps(mContent.size(), 1);
				const double& dStepSize = 100. / fN;
				vector<double> vNextThrowOutIdx(mContent.size(), dStepSize);


				auto itOut = vOut->begin();
				uint64_t iCounter = 0;
				for (const auto& entry : *vLib) {
					const auto& idx = Utilities::checkIfInMap(mContent, entry.second)->second;
					if (vSteps[idx] != static_cast<uint64_t>(vNextThrowOutIdx[idx])) {
						*itOut = entry;
						countFreqs(mContent, freqArray, entry);
						itOut++;
						++iCounter;
					}
					else {
						vNextThrowOutIdx[idx] += dStepSize;
					}
					++vSteps[idx];
				}

				vOut->resize(iCounter, true);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
		

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Shrink the library by throwing all redundant kMers away
		void ShrinkLib(const string& sLibFile, const string& fOutFile, const ShrinkingStrategy& eDeleteStrategy, const string& sContentFile, const float& iPercOfThrownAway = 50.f) {
			try {
				// test if files exists
				if (!ifstream(sLibFile)) {
					throw runtime_error("The library file could not be found");
				}
				if (!ifstream(sContentFile)) {
					throw runtime_error("Content file not found");
				}

				uint32_t iIdxCounter = 1;
				unordered_map<uint32_t, uint32_t> mIDsAsIdx; mIDsAsIdx[0] = 0;
				unordered_map<uint32_t, string> mIdxToName; mIdxToName[0] = "non_unique";
				ifstream content(sContentFile);
				string sDummy = "";
				while (getline(content, sDummy)) {
					if (sDummy != "") {
						const auto& line = Utilities::split(sDummy, '\t');
						if (line.size() == 4) {
							mIDsAsIdx[stoul(line[1])] = iIdxCounter;
							mIdxToName[iIdxCounter] = line[0];
							++iIdxCounter;
						}
					}
				}

				

				unique_ptr<uint64_t[]> arrFrequencies(new uint64_t[iIdxCounter * 12]);
				for (uint64_t i = 0; i < iIdxCounter * 12; ++i) {
					arrFrequencies[i] = 0;
				}

				// get lib
				fstream fLibInfo(sLibFile + "_info.txt", ios::in);
				uint64_t iSizeOfLib = 0;
				fLibInfo >> iSizeOfLib;
				fLibInfo.close();
				unique_ptr<stxxlFile> stxxlLibFile(new stxxlFile(sLibFile, stxxl::file::RDONLY));
				unique_ptr<unique_ptr<const contentVecType_32p>[]> vLibIn(new unique_ptr<const contentVecType_32p>[_iNumOfThreads]);
				for (int32_t i = 0; i < _iNumOfThreads; ++i) {
					vLibIn[i].reset(new const contentVecType_32p(stxxlLibFile.get(), iSizeOfLib));
				}

				// create reduced vec
				ofstream derp;
				derp.exceptions(std::ifstream::failbit | std::ifstream::badbit);
				derp.open(fOutFile);
				derp.close();


				unique_ptr<stxxlFile> stxxlOutFile;
				unique_ptr<contentVecType_32p> vOutVec;

				unique_ptr<stxxlFile> stxxlFileP;
				unique_ptr<index_t_p> vOutPVec;

				if (eDeleteStrategy == TrieHalf || eDeleteStrategy == Overrepresented) {
					stxxlFileP.reset(new stxxlFile(fOutFile, stxxl::file::RDWR));
					vOutPVec.reset(new index_t_p(stxxlFileP.get(), iSizeOfLib));
				}
				else {
					stxxlOutFile.reset(new stxxlFile(fOutFile, stxxl::file::RDWR));
					vOutVec.reset(new contentVecType_32p(stxxlOutFile.get(), iSizeOfLib));
				}

				// copy and rename in case of half
				ifstream  oldFreqFile;
				ofstream  newFreqFile;

				// reduce
				switch (eDeleteStrategy) {
				case EveryNth:
				{
					deleteEveryNth(vLibIn[0], vOutVec, fabsf(iPercOfThrownAway), mIDsAsIdx, arrFrequencies);

					newFreqFile.open(fOutFile + "_f.txt");
					for (uint32_t j = 0; j < iIdxCounter; ++j) {
						newFreqFile << Utilities::checkIfInMap(mIdxToName, j)->second << "\t";
						newFreqFile << arrFrequencies[j * 12];
						for (int32_t k = 1; k < 12; ++k) {
							newFreqFile << "\t" << arrFrequencies[j * 12 + k];
						}
						newFreqFile << endl;
					}

					fLibInfo.open(fOutFile + "_info.txt", ios::out);
					fLibInfo << vOutVec->size();
					vOutVec->export_files("_");

					Trie T(static_cast<int8_t>(12), static_cast<int8_t>(_iMinK), 6);
					T.SaveToStxxlVec(vOutVec.get(), fOutFile);

					break;
				}
				case TrieHalf:
				{
					assert(iIdxCounter <= 65535);

					bool bTrieIsNotThere = false;
					if (!ifstream(sLibFile + "_trie.txt")) {
						bTrieIsNotThere = true;
					}

					putHalfInTrie(vLibIn[0], vOutPVec, mIDsAsIdx, fOutFile, bTrieIsNotThere);

					fLibInfo.open(fOutFile + "_info.txt", ios::out);
					fLibInfo << vOutPVec->size() << endl << 3;
					vOutPVec->export_files("_");
					fLibInfo.close();

					if (!bTrieIsNotThere) { // Trie is already there so copy it
						oldFreqFile.open(sLibFile + "_trie.txt", std::ios::binary);
						newFreqFile.open(fOutFile + "_trie.txt", std::ios::binary);
						newFreqFile << oldFreqFile.rdbuf();
						oldFreqFile.close();
						newFreqFile.close();

						oldFreqFile.open(sLibFile + "_trie", std::ios::binary);
						newFreqFile.open(fOutFile + "_trie", std::ios::binary);
						newFreqFile << oldFreqFile.rdbuf();
						oldFreqFile.close();
						newFreqFile.close();
					}


					oldFreqFile.open(sLibFile + "_f.txt", std::ios::binary);
					newFreqFile.open(fOutFile + "_f.txt", std::ios::binary);

					newFreqFile << oldFreqFile.rdbuf();

					break;
				}
				case Overrepresented:
				{
					//deleteOverrepresentedAndputHalfInTrie(vLibIn[0], vOutPVec, mIDsAsIdx, fOutFile, arrFrequencies);

					fLibInfo.open(fOutFile + "_info.txt", ios::out);
					fLibInfo << vOutPVec->size() << endl << 3;
					vOutPVec->export_files("_");
					fLibInfo.close();

					newFreqFile.open(fOutFile + "_f.txt");
					for (uint32_t j = 0; j < iIdxCounter; ++j) {
						newFreqFile << Utilities::checkIfInMap(mIdxToName, j)->second << "\t";
						newFreqFile << arrFrequencies[j * 12];
						for (int32_t k = 1; k < 12; ++k) {
							newFreqFile << "\t" << arrFrequencies[j * 12 + k];
						}
						newFreqFile << endl;
					}

					break;
				}
				};

			
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

	};

}