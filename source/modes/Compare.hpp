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

	enum OutputFormat
	{
		Kraken,
		Json,
		JsonL,
		tsv
	};

	template<class contentVecType, class elemType, class intType>
	class Compare : public Read<contentVecType, elemType, intType> {

		typedef Read<contentVecType, elemType, intType> Base;

		bool bUnfunny = false, bCoverage = false;
		mutex m_exceptionLock;
		exception_ptr someThingWentWrong;

		vector<pair<string, uint32_t>> _matchedkMers;

	public:
		Compare(const InputParameters& cParams) : Read<contentVecType, elemType, intType>(cParams), bUnfunny(cParams.bUnfunny), bCoverage(cParams.bCoverage), iNumOfBeasts(cParams.iNumOfBeasts) {}

		// for output
		int32_t iNumOfBeasts = 3;

		OutputFormat format = OutputFormat::Json;

		struct ReadIndex {
			shared_ptr<stxxlFile> stxxlLibFile;
			unique_ptr<unique_ptr<const contentVecType>[]> vLib;
			unique_ptr<unique_ptr<const index_t_p>[]> vLibParted_p;
			unique_ptr<unique_ptr<const taxaOnly>[]> vLib_taxaOnly;
			vector<elemType> vLib_RAM_Full;
			vector<packedPair> vLib_RAM_Half;
			vector<uint16_t> vLib_RAM_taxaOnly;

			shared_ptr<Trie<intType>> trieForVector;


			bool bRAM = false;
			bool bUnfunny = false;
			bool bPartitioned = false;
			bool bIdentifyMultiple = false;
			int32_t iNumOfThreads = 1;
			int64_t iBytesUsedBySTXXLVectors = 0;
			uint64_t iNumOfSpecies = 0;
			uint64_t iSizeOfLib = 0;

			shared_ptr< vector<string> > mOrganisms;
			shared_ptr< unordered_map<uint32_t, uint32_t> > mTaxToIdx;
			shared_ptr< vector<uint32_t> > mIdxToTax;
			vector<vector<uint64_t>> mFrequencies;

			ReadIndex() {}

			ReadIndex(const ReadIndex& existingIndexInstance) {
				stxxlLibFile = existingIndexInstance.stxxlLibFile;
				trieForVector = existingIndexInstance.trieForVector;
				mOrganisms = existingIndexInstance.mOrganisms;
				mTaxToIdx = existingIndexInstance.mTaxToIdx;
				mIdxToTax = existingIndexInstance.mIdxToTax;
				mFrequencies = existingIndexInstance.mFrequencies;
				bPartitioned = existingIndexInstance.bPartitioned;
				bIdentifyMultiple = existingIndexInstance.bIdentifyMultiple;
				iNumOfSpecies = existingIndexInstance.iNumOfSpecies;
				iSizeOfLib = existingIndexInstance.iSizeOfLib;
			}

			void set(const bool& ram, const bool& unfunny, const int32_t& threads) {
				bRAM = ram;
				bUnfunny = unfunny;
				iNumOfThreads = threads;
			}

			void setFile(const string& sLibFile) {
				if (!ifstream(sLibFile)) {
					throw runtime_error("The index file cannot be found!");
				}
				ifstream fLibInfo(sLibFile + "_info.txt");
				fLibInfo >> iSizeOfLib;
				uint32_t iTempVecType = 0;
				fLibInfo >> iTempVecType;
				if (iTempVecType == 3) {
					bPartitioned = true;
				}
				fLibInfo.close();
				stxxlLibFile.reset(new stxxlFile(sLibFile, stxxl::file::RDONLY));
			}

			inline void loadContentAndFrequencyFiles(const int32_t& iMaxK, const int32_t& iMinK, const string& sLibFile, const string& sContentFile, int64_t& iAvailableMemory) {

				if (!ifstream(sContentFile) || !ifstream(sLibFile + "_f.txt")) {
					throw runtime_error("The content file or the frequency file cannot be found!");
				}

				this->mOrganisms.reset(new vector<string>);
				this->mTaxToIdx.reset(new unordered_map<uint32_t, uint32_t>);
				this->mIdxToTax.reset(new vector<uint32_t>);

				ifstream fContent(sContentFile);
				uint32_t iNumberOfSpecies = 1;
				string sTempLine = "";
				bool bTaxIdsAsStrings = false;
				this->mTaxToIdx->insert(make_pair(0, 0));
				this->mOrganisms->push_back("non_unique");
				this->mIdxToTax->push_back(0);
				while (getline(fContent, sTempLine)) {
					if (sTempLine != "") {
						const auto& tempLineContent = Utilities::split(sTempLine, '\t');
						if (tempLineContent.size() >= 5 && !bTaxIdsAsStrings) {
							bTaxIdsAsStrings = true;
						}
						if (tempLineContent.size() >= 4) {
							this->mOrganisms->push_back(Utilities::removeCharFromString(tempLineContent[0], ','));
							iAvailableMemory -= this->mOrganisms->back().size() * sizeof(char);
							if (bTaxIdsAsStrings) {
								this->mIdxToTax->push_back(stoul(tempLineContent[4]));
								this->mTaxToIdx->insert(make_pair(stoul(tempLineContent[4]), iNumberOfSpecies++));
							}
							else {
								this->mIdxToTax->push_back(stoul(tempLineContent[1]));
								this->mTaxToIdx->insert(make_pair(stoul(tempLineContent[1]), iNumberOfSpecies++));
							}
						}
						else {
							throw runtime_error("Content file contains less than 4 columns, it may be damaged... The faulty line was: " + sTempLine + "\n");
						}
					}
				}
				fContent.close();

				this->iNumOfSpecies = iNumberOfSpecies;

				// hash tables do cost some memory!
				iAvailableMemory -= Utilities::calculateSizeInByteOfUnorderedMap(*mTaxToIdx);
				iAvailableMemory -= this->mIdxToTax->size() * sizeof(uint32_t);
				iAvailableMemory -= iNumberOfSpecies * sizeof(uint64_t) * (iMaxK - iMinK); // frequencies

				if (iAvailableMemory < 0) {
					cerr << "WARNING: Loading the metadata (content file, frequency file) consumes more memory than given. kASA will try to continue with 1GB nevertheless..." << endl;
					iAvailableMemory = GIGABYTEASBYTES;
				}

				// get frequencies
				ifstream fFrequencies(sLibFile + "_f.txt");
				this->mFrequencies.resize(iNumberOfSpecies);
				uint32_t iCounterForFreqs = 0;
				while (getline(fFrequencies, sTempLine)) {
					if (sTempLine != "") {
						const auto& vLine = Utilities::split(sTempLine, '\t');
						const auto& numOfK = vLine.size() - 1;
						mFrequencies[iCounterForFreqs].resize(numOfK);
						for (uint8_t i = iMaxK, j = 0; i >= iMinK; --i, ++j) {
							mFrequencies[iCounterForFreqs][j] = stoull(vLine[1 + numOfK - i]);
						}
						iCounterForFreqs++;
					}
				}
			}

			inline void loadIndex(const int32_t& iMaxk, const int32_t& iMink, int64_t& iAvailableMemory) {
				if (bRAM) {
					try {
						if (bPartitioned || (iMink > 6 && mTaxToIdx->size() <= 65535 && !bUnfunny && iMaxk <= 12)) {
							if (iAvailableMemory - static_cast<int64_t>(iSizeOfLib * sizeof(packedPair)) < 0) {
								cerr << "WARNING: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
								bRAM = false;
							}
							else {
								vLib_RAM_Half.reserve(iSizeOfLib);

								if (bPartitioned) {
									vLibParted_p.reset(new unique_ptr<const index_t_p>[1]);
									vLibParted_p[0].reset(new const index_t_p(stxxlLibFile.get(), iSizeOfLib));

									stxxl::vector_bufreader<index_t_p::const_iterator> bufferedReader(vLibParted_p[0]->cbegin(), vLibParted_p[0]->cend(), 0);
									for (; !bufferedReader.empty(); ++bufferedReader) {
										vLib_RAM_Half.push_back(*bufferedReader);
									}

									vLibParted_p[0].reset();
									vLibParted_p.reset();
									stxxlLibFile.reset();
								}
								else {
									vLib.reset(new unique_ptr<const contentVecType>[1]);
									vLib[0].reset(new const contentVecType(stxxlLibFile.get(), iSizeOfLib));

									stxxl::vector_bufreader<typename contentVecType::const_iterator> bufferedReader(vLib[0]->cbegin(), vLib[0]->cend(), 0);
									for (; !bufferedReader.empty(); ++bufferedReader) {
										vLib_RAM_Half.push_back(packedPair(static_cast<uint32_t>(bufferedReader->first & 1073741823ULL), static_cast<uint16_t>(mTaxToIdx->at(bufferedReader->second))));
									}

									vLib[0].reset();
									vLib.reset();
									stxxlLibFile.reset();
									bPartitioned = true;
								}

								iAvailableMemory -= iSizeOfLib * sizeof(packedPair);
							}
						}
						else {
							if (bUnfunny) {
								if (iAvailableMemory - static_cast<int64_t>(iSizeOfLib * sizeof(uint16_t)) < 0) {
									cerr << "WARNING: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
									bRAM = false;
								}
								else {
									vLib_taxaOnly.reset(new unique_ptr<const taxaOnly>[1]);
									vLib_taxaOnly[0].reset(new const taxaOnly(stxxlLibFile.get(), iSizeOfLib));

									vLib_RAM_taxaOnly.reserve(iSizeOfLib);
									stxxl::vector_bufreader<taxaOnly::const_iterator> bufferedReader(vLib_taxaOnly[0]->cbegin(), vLib_taxaOnly[0]->cend(), 0);
									for (; !bufferedReader.empty(); ++bufferedReader) {
										vLib_RAM_taxaOnly.push_back(*bufferedReader);
									}

									vLib_taxaOnly[0].reset();
									vLib_taxaOnly.reset();
									stxxlLibFile.reset();
									iAvailableMemory -= iSizeOfLib * sizeof(uint16_t);
								}
							}
							else {
								if (iAvailableMemory - static_cast<int64_t>(iSizeOfLib * sizeof(elemType)) < 0) {
									cerr << "WARNING: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
									bRAM = false;
								}
								else {
									vLib.reset(new unique_ptr<const contentVecType>[1]);
									vLib[0].reset(new const contentVecType(stxxlLibFile.get(), iSizeOfLib));

									vLib_RAM_Full.reserve(iSizeOfLib);
									stxxl::vector_bufreader<typename contentVecType::const_iterator> bufferedReader(vLib[0]->cbegin(), vLib[0]->cend(), 0);
									for (; !bufferedReader.empty(); ++bufferedReader) {
										vLib_RAM_Full.push_back(*bufferedReader);
									}

									vLib[0].reset();
									vLib.reset();
									stxxlLibFile.reset();
									iAvailableMemory -= iSizeOfLib * sizeof(elemType);
								}
							}
						}
					}
					catch (const bad_alloc&) {
						cerr << "WARNING: Not enough RAM available to load index into it. Resuming with secondary memory approach..." << endl;
						bRAM = false;
						vLib_RAM_Half.clear();
						vLib_RAM_taxaOnly.clear();
						vLib_RAM_Full.clear();
					}
				}
				else {
					try {
						if (bPartitioned) {
							iBytesUsedBySTXXLVectors = uint64_t(iNumOfThreads) * index_t_p::block_size * index_t_p::page_size * 4; // 4 = number of cache pages

							if (iAvailableMemory - iBytesUsedBySTXXLVectors < 0) {
								cerr << "WARNING: The overhead from creating so many threads is consuming more memory than given. Please consider using less threads or switch to RAM mode, else kASA may crash!" << endl;
							}

							vLibParted_p.reset(new unique_ptr<const index_t_p>[iNumOfThreads]);
							//stxxlLibFile.reset(new unique_ptr<stxxlFile>[iNumOfThreads]);
							for (int32_t i = 0; i < iNumOfThreads; ++i) {
								//stxxlLibFile[i].reset(new stxxlFile(sLibFile + "_"+ to_string(i), stxxl::file::RDONLY));
								vLibParted_p[i].reset(new const index_t_p(stxxlLibFile.get(), iSizeOfLib));
							}

						}
						else {
							if (bUnfunny) {
								iBytesUsedBySTXXLVectors = uint64_t(iNumOfThreads) * taxaOnly::block_size * taxaOnly::page_size * 4;

								if (iAvailableMemory - iBytesUsedBySTXXLVectors < 0) {
									cerr << "WARNING: The overhead from creating so many threads is consuming more memory than given. Please consider using less threads or switch to RAM mode, else kASA may crash!" << endl;
								}

								vLib_taxaOnly.reset(new unique_ptr<const taxaOnly>[iNumOfThreads]);
								for (int32_t i = 0; i < iNumOfThreads; ++i) {
									vLib_taxaOnly[i].reset(new const taxaOnly(stxxlLibFile.get(), iSizeOfLib));
								}

							}
							else {
								iBytesUsedBySTXXLVectors = uint64_t(iNumOfThreads) * contentVecType::block_size * contentVecType::page_size * 4;

								if (iAvailableMemory - iBytesUsedBySTXXLVectors < 0) {
									cerr << "WARNING: The overhead from creating so many threads is consuming more memory than given. Please consider using less threads or switch to RAM mode, else kASA may crash!" << endl;
								}

								vLib.reset(new unique_ptr<const contentVecType>[iNumOfThreads]);
								for (int32_t i = 0; i < iNumOfThreads; ++i) {
									vLib[i].reset(new const contentVecType(stxxlLibFile.get(), iSizeOfLib));
								}
							}
						}
					}
					catch (...) {
						throw;
					}

					iAvailableMemory -= iBytesUsedBySTXXLVectors;
				}
			}

			inline void loadTrie(const string& sLibFile, const int32_t& iMaxk, const int32_t& iMink) {
				if (!ifstream(sLibFile + "_trie.txt") || !ifstream(sLibFile + "_trie")) {
					throw runtime_error("The trie file cannot be found!");
				}
				trieForVector.reset(new Trie<intType>(iMaxk, iMink, 6));
				trieForVector->LoadFromStxxlVec(sLibFile);
				trieForVector->SetForIsInTrie((iMink < 6) ? static_cast<uint8_t>(iMink) : static_cast<uint8_t>(6));
			}

			inline unique_ptr<unique_ptr<const contentVecType>[]>* getContentVecType() {
				return &vLib;
			}

			inline unique_ptr<unique_ptr<const index_t_p>[]>* getIndexTP() {
				return &vLibParted_p;
			}

			inline unique_ptr<unique_ptr<const taxaOnly>[]>* getTaxaOnly() {
				return &vLib_taxaOnly;
			}

			inline vector<elemType>* getLibRAMFull() {
				return &vLib_RAM_Full;
			}

			inline vector<packedPair>* getLibRAMHalf() {
				return &vLib_RAM_Half;
			}

			inline vector<uint16_t>* getLibRAMTaxaOnly() {
				return &vLib_RAM_taxaOnly;
			}

		} index;

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

		public:
		inline string getOutputFormatFileEnding() {
			return outputFormatFileEnding();
		}

		private:

		///////////////////////////////////////////////////////
		//const float arrWeightingFactors[12] = { 1.f , 121.f / 144.f , 100.f / 144.f , 81.f / 144.f , 64.f / 144.f , 49.f / 144.f , 36.f / 144.f , 25.f / 144.f , 16.f / 144.f , 9.f / 144.f , 4.f / 144.f , 1.f / 144.f };
		const float arrWeightingFactors[HIGHESTPOSSIBLEK] = { 1.f , 576.f / 625.f , 529.f / 625.f , 484.f / 625.f , 441.f / 625.f , 400.f / 625.f , 361.f / 625.f , 324.f / 625.f , 289.f / 625.f , 256.f / 625.f , 225.f / 625.f , 196.f / 625.f , 169.f / 625.f , 144.f / 625.f , 121.f / 625.f , 100.f / 625.f , 81.f / 625.f , 64.f / 625.f , 49.f / 625.f , 36.f / 625.f , 25.f / 625.f , 16.f / 625.f , 9.f / 625.f , 4.f / 625.f , 1.f / 625.f };
		//const float arrWeightingFactors[HIGHESTPOSSIBLEK] = { 1.f, 0.84934656, 0.71639296, 0.59969536, 0.49787136, 0.4096, 0.33362176, 0.26873856, 0.21381376, 0.16777216, 0.1296, 0.09834496, 0.07311616, 0.05308416, 0.03748096, 0.0256, 0.01679616, 0.01048576, 0.00614656, 0.00331776, 0.0016, 0.00065536, 0.00020736, 4.096E-05, 2.56E-06 };
		//const float arrWeightingFactors[12] = { 1, 1331.f / 1728.f , 1000.f / 1728.f , 729.f / 1728.f , 512.f / 1728.f , 343.f / 1728.f , 216.f / 1728.f , 125.f / 1728.f , 64.f / 1728.f , 27.f / 1728.f , 8.f / 1728.f , 1.f / 1728.f };
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inline void markTaxIDs(const uint64_t& codedTaxIDs, Utilities::sBitArray& vMemoryOfTaxIDs_k, const unordered_map<uint32_t, uint32_t>& mTaxToIdx, unique_ptr<const contentVecType>&) {
			try {
				vMemoryOfTaxIDs_k.set(Utilities::checkIfInMap(mTaxToIdx, static_cast<uint32_t>(codedTaxIDs))->second);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}
		inline void markTaxIDs(const uint64_t& codedTaxIDs, Utilities::sBitArray& vMemoryOfTaxIDs_k, const unordered_map<uint32_t, uint32_t>& mTaxToIdx, const vector<packedLargePair>*) {
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
		inline unique_ptr<const contentVecType_128>& getVec(const unique_ptr<unique_ptr<const contentVecType_128>[]>& vec, const int& iThreadID) {
			return vec[iThreadID];
		}
		//inline const unique_ptr<const contentVecType_32p[]>& getVec(const unique_ptr<const contentVecType_32p[]>& vec, const int&) {
		//	return vec;
		//}
		inline unique_ptr<const index_t_p>& getVec(const unique_ptr<unique_ptr<const index_t_p>[]>& vec, const int& iThreadID) {
			return vec[iThreadID];
		}
		inline const vector<packedLargePair>* getVec(const vector<packedLargePair>* vec, const int&) {
			return vec;
		}
		inline const vector<packedBigPair>* getVec(const vector<packedBigPair>* vec, const int&) {
			return vec;
		}
		inline const vector<packedPair>* getVec(const vector<packedPair>* vec, const int&) {
			return vec;
		}

		const vector<string> vMasks{ "1111111011111111111111111", "1111111011011111111111111", "1111111011010111111111111", "1111111011010101111111111", "1111111011010100111111111", "1111111011010100101111111", "1111111011010100101011111", "1111111011010100101011011", "1111111011010100101011001"  }; // 13, 14, 15, 16, 17, 19, 21, 24, 25

		inline uint8_t compareTwoKmers(const intType& inKmer, const intType& idxKmer, const int32_t&) {
			if (inKmer < idxKmer) {
				return 0;
			}
			else {
				if (!(idxKmer < inKmer)) {
					return 1;
				}
				else {
					return 2;
				}
			}

			return 3;
		}

		inline uint8_t compareTwoKmersSpaced(const intType& inKmer, const intType& idxKmer, const int32_t& k) {
			
			if (vMasks[Base::_iWhichMask][k] == '0') {
				//if (inKmer != idxKmer) {
					//cout << k << " " << Base::kMerToAminoacid(inKmer, 25) << " " << Base::kMerToAminoacid(idxKmer, 25) << endl;
				//}

				return 1;
			}
			return compareTwoKmers(inKmer, idxKmer, 0);
		}

		inline void scoreMatchNonAVX(Utilities::Non_contiguousArray& vReadIDtoGenID, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, unique_ptr<uint64_t[]>& vCountTotal, function<void(unique_ptr<uint64_t[]>&, const uint64_t&)>& countTotal, const Utilities::sBitArray& taxIDs, const uint64_t& numOfEntries, const uint64_t& iPartialTempIndex, const vector<uint64_t>& vReadIDs, const uint64_t& numOfHits, const double& counts, const float& score) {
			for (auto it = taxIDs.begin(); it != taxIDs.end(); ++it) {
				const auto& tempIndex = iPartialTempIndex + (*it);

				vCount[tempIndex] += counts;
				
				countTotal(vCountTotal, tempIndex);

				if (numOfEntries == 1) {
					vCountUnique[tempIndex] += numOfHits;
				}

				for (auto readIDIt = vReadIDs.cbegin(); readIDIt != vReadIDs.cbegin() + numOfHits; readIDIt++) {
					vReadIDtoGenID[*readIDIt][*it] += score;
				}
			}
		}

		inline void scoreMatchAVX(Utilities::Non_contiguousArray& vReadIDtoGenID, float* tmpStore, array<pair<uint64_t, uint64_t>,8>& aOfIndices, double* tmpStore2, std::array<uint64_t,4>& aOfTempIdx, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountTotal, function<void(unique_ptr<uint64_t[]>&, const uint64_t&)>& countTotal, const Utilities::sBitArray& taxIDs, const uint64_t& iPartialTempIndex, const vector<uint64_t>& vReadIDs, const uint64_t& numOfHits, const double& counts, const float& score) {
			__m256 avx_score = _mm256_broadcast_ss(&score);
			__m256d avx_counts = _mm256_broadcast_sd(&counts);
			uint32_t iBlockCount = 0, iIdxCount = 0;
			for (auto it = taxIDs.begin(); it != taxIDs.end(); ++it) {
				const auto& tempIndex = iPartialTempIndex + (*it);
				aOfTempIdx[iIdxCount] = tempIndex;
				tmpStore2[iIdxCount] = vCount[tempIndex];
				iIdxCount++;

				if (iIdxCount == 4) {
					//vCount[tempIndex] += counts;
					//__m256d avx_taxIDs = _mm256_setr_pd(vCount[aOfTempIdx[0]], vCount[aOfTempIdx[1]], vCount[aOfTempIdx[2]], vCount[aOfTempIdx[3]]);
					__m256d avx_taxIDs = _mm256_load_pd(tmpStore2);
					__m256d avx_result = _mm256_add_pd(avx_taxIDs, avx_counts);
					_mm256_store_pd(tmpStore2, avx_result);

					for (uint32_t iResCount = 0; iResCount < 4; iResCount++) {
						vCount[aOfTempIdx[iResCount]] = tmpStore2[iResCount];
					}
					iIdxCount = 0;
				}

				countTotal(vCountTotal, tempIndex);

				for (auto readIDIt = vReadIDs.cbegin(); iBlockCount < 8 && readIDIt != vReadIDs.cbegin() + numOfHits; iBlockCount++, readIDIt++) {
					aOfIndices[iBlockCount].first = *readIDIt;
					aOfIndices[iBlockCount].second = *it;
					tmpStore[iBlockCount] = vReadIDtoGenID[*readIDIt][*it];
				}
				if (iBlockCount == 8) {
					//__m256 avx_readIDs = _mm256_setr_ps(vReadIDtoGenID[aOfIndices[0].first][aOfIndices[0].second], vReadIDtoGenID[aOfIndices[1].first][aOfIndices[1].second], vReadIDtoGenID[aOfIndices[2].first][aOfIndices[2].second], vReadIDtoGenID[aOfIndices[3].first][aOfIndices[3].second], vReadIDtoGenID[aOfIndices[4].first][aOfIndices[4].second], vReadIDtoGenID[aOfIndices[5].first][aOfIndices[5].second], vReadIDtoGenID[aOfIndices[6].first][aOfIndices[6].second], vReadIDtoGenID[aOfIndices[7].first][aOfIndices[7].second]);
					__m256 avx_readIDs = _mm256_load_ps(tmpStore);
					__m256 avx_result = _mm256_add_ps(avx_readIDs, avx_score);
					_mm256_store_ps(tmpStore, avx_result);

					for (uint32_t iResCount = 0; iResCount < 8; iResCount++) {
						vReadIDtoGenID[aOfIndices[iResCount].first][aOfIndices[iResCount].second] = tmpStore[iResCount];
					}

					iBlockCount = 0;
				}
			}

			if (iIdxCount > 0) {
				__m256d avx_taxIDs = _mm256_load_pd(tmpStore2);
				__m256d avx_result = _mm256_add_pd(avx_taxIDs, avx_counts);
				_mm256_store_pd(tmpStore2, avx_result);

				for (uint32_t iResCount = 0; iResCount < iIdxCount; iResCount++) {
					vCount[aOfTempIdx[iResCount]] = tmpStore2[iResCount];
				}
			}

			if (iBlockCount > 0) {
				__m256 avx_readIDs = _mm256_load_ps(tmpStore);
				__m256 avx_result = _mm256_add_ps(avx_readIDs, avx_score);
				_mm256_store_ps(tmpStore, avx_result);

				for (uint32_t iResCount = 0; iResCount < iBlockCount; iResCount++) {
					vReadIDtoGenID[aOfIndices[iResCount].first][aOfIndices[iResCount].second] = tmpStore[iResCount];
				}
			}
		}

		inline void scoreMatchNonAVX_countsOnly(unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, unique_ptr<uint64_t[]>& vCountTotal, function<void(unique_ptr<uint64_t[]>&, const uint64_t&)>& countTotal, const Utilities::sBitArray& taxIDs, const uint64_t& numOfEntries, const uint64_t& numOfHits, const uint64_t& iPartialTempIndex, const double& counts) {
			for (auto it = taxIDs.begin(); it != taxIDs.end(); ++it) {
				const auto& tempIndex = iPartialTempIndex + (*it);

				vCount[tempIndex] += counts;
				
				countTotal(vCountTotal,tempIndex);
				
				if (numOfEntries == 1) {
					vCountUnique[tempIndex] += numOfHits;
				}
			}
		}

		inline void scoreMatchAVX_countsOnly(double* tmpStore2, array<uint64_t, 4>& aOfTempIdx, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountTotal, function<void(unique_ptr<uint64_t[]>&, const uint64_t&)>& countTotal, const Utilities::sBitArray& taxIDs, const uint64_t& iPartialTempIndex, const double& counts) {
			__m256d avx_counts = _mm256_broadcast_sd(&counts);
			uint32_t iIdxCount = 0;
			for (auto it = taxIDs.begin(); it != taxIDs.end(); ++it) {
				const auto& tempIndex = iPartialTempIndex + (*it);
				aOfTempIdx[iIdxCount] = tempIndex;
				iIdxCount++;

				if (iIdxCount == 4) {
					//vCount[tempIndex] += counts;
					__m256d avx_taxIDs = _mm256_setr_pd(vCount[aOfTempIdx[0]], vCount[aOfTempIdx[1]], vCount[aOfTempIdx[2]], vCount[aOfTempIdx[3]]);
					__m256d avx_result = _mm256_add_pd(avx_taxIDs, avx_counts);
					_mm256_store_pd(tmpStore2, avx_result);

					for (uint32_t iResCount = 0; iResCount < 4; iResCount++) {
						vCount[aOfTempIdx[iResCount]] = tmpStore2[iResCount];
					}
					iIdxCount = 0;
				}

				countTotal(vCountTotal, tempIndex);
			}

			if (iIdxCount > 0) {
				__m256d avx_taxIDs = _mm256_setr_pd(vCount[aOfTempIdx[0]], vCount[aOfTempIdx[1]], vCount[aOfTempIdx[2]], vCount[aOfTempIdx[3]]);
				__m256d avx_result = _mm256_add_pd(avx_taxIDs, avx_counts);
				_mm256_store_pd(tmpStore2, avx_result);

				for (uint32_t iResCount = 0; iResCount < iIdxCount; iResCount++) {
					vCount[aOfTempIdx[iResCount]] = tmpStore2[iResCount];
				}
			}
		}


		inline void scoreMatchForReadIDsAndTaxIDs(Utilities::Non_contiguousArray& vReadIDtoGenID, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, unique_ptr<uint64_t[]>& vCountTotal, function<void(unique_ptr<uint64_t[]>&, const uint64_t&)>& countTotal, float* tmpStore, array<pair<uint64_t,uint64_t>,8>& aOfIndices, double* tmpStore2, array<uint64_t, 4>& aOfTempIdx, const Utilities::sBitArray& taxIDs, const uint64_t& numOfEntries, const uint64_t& iPartialTempIndex, const vector<uint64_t>& vReadIDs, const uint64_t& numOfHits, const double& counts, const float& score) {
			if (vReadIDtoGenID.size()) {
#if defined(__AVX__) || defined(_INCLUDED_IMM)
				if (numOfEntries > 3) {
					scoreMatchAVX(vReadIDtoGenID, tmpStore, aOfIndices, tmpStore2, aOfTempIdx, vCount, vCountTotal, countTotal, taxIDs, iPartialTempIndex, vReadIDs, numOfHits, counts, score);
				}
				else {
					scoreMatchNonAVX(vReadIDtoGenID, vCount, vCountUnique, vCountTotal, countTotal, taxIDs, numOfEntries, iPartialTempIndex, vReadIDs, numOfHits, counts, score);
				}
#else
				scoreMatchNonAVX(vReadIDtoGenID, vCount, vCountUnique, vCountTotal, countTotal, taxIDs, numOfEntries, iPartialTempIndex, vReadIDs, numOfHits, counts, score);
#endif // __AVX__
			}
			else {
#if defined(__AVX__) || defined(_INCLUDED_IMM)
				if (numOfEntries > 3) {
					scoreMatchAVX_countsOnly(tmpStore2, aOfTempIdx, vCount, vCountTotal, countTotal, taxIDs, iPartialTempIndex, counts);
				}
				else {
					scoreMatchNonAVX_countsOnly(vCount, vCountUnique, vCountTotal, countTotal, taxIDs, numOfEntries, numOfHits, iPartialTempIndex, counts);
				}
#else
				scoreMatchNonAVX_countsOnly(vCount, vCountUnique, vCountTotal, countTotal, taxIDs, numOfEntries, numOfHits, iPartialTempIndex, counts);
#endif // __AVX__
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Compare as many as #Number-of-processors vectors with an index lying on a HDD/SSD and note all similarities for any k. 
		// To minimize hard disk access, the order is as follows: Get kMer from RAM Vec -> Search in Prefix-Trie -> Get range of possible hit -> binary search in that range -> note if hit
		template <typename vecType>
		inline void compareWithDatabase(const int32_t& iThreadID, InputType<intType>& vIn, const uint64_t& vInStart, const uint64_t& vInEnd, const vecType& vLib, unique_ptr<double[]>& vCount, unique_ptr<uint64_t[]>& vCountUnique, unique_ptr<uint64_t[]>& vCountTotal, Utilities::Non_contiguousArray& vReadIDtoGenID, const uint64_t& iSpecIDRange, const unordered_map<uint32_t, uint32_t>& mTaxToIdx) {//, const unordered_map<readIDType, uint64_t>& mReadIDToArrayIdx) {

			try {
				debugBarrier

				//double dAverageLoopCount = 0.;
				//size_t iNumOfLoops = 0;
				// In case of profiling only, this is the default
				function<void(vector<uint64_t>&, uint64_t&, const uint64_t&)> addToMatchedReadID = [](vector<uint64_t>&, uint64_t& position, const uint64_t&) {
					position++;
				};
				function<void(unique_ptr<uint64_t[]>&, const uint64_t&)> countTotal = [](unique_ptr<uint64_t[]>& vCountTotal, const uint64_t& idx) { vCountTotal[idx] += 1; };
				if (!this->bCoverage) {
					countTotal = [](unique_ptr<uint64_t[]>&, const uint64_t&) {};
				}

				alignas(32) float tmpStore[8] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };
				alignas(64) double tmpStore2[4] = { 0.f, 0.f, 0.f, 0.f};
				array<pair<uint64_t, uint64_t>,8> aOfIndices = { make_pair(0ull,0ull), make_pair(0ull,0ull), make_pair(0ull,0ull), make_pair(0ull,0ull), make_pair(0ull,0ull), make_pair(0ull,0ull), make_pair(0ull,0ull), make_pair(0ull,0ull) };
				array<uint64_t, 4> aOfTempIdx = {0ull, 0ull, 0ull, 0ull};
				/////////////////////////////////////// Spaced try
				function<uint8_t(const uint64_t&, const uint64_t&, const int32_t&)> compare;
				if (Base::_bSpaced) {
					compare = bind(&Compare::compareTwoKmersSpaced, this, placeholders::_1, placeholders::_2, placeholders::_3);
				}
				else {
					compare = bind(&Compare::compareTwoKmers, this, placeholders::_1, placeholders::_2, 0);
				}
				///////////////////////////////////////

				/////// Statistics
//#define STATISTICS
#ifdef STATISTICS
				int64_t iterationVar = 0;
				vector<double> statsVec_numOfHits(3, 0.0); // lowest, highest, average for numOfHits ...
				vector<double> statsVec_numOfEntries(3, 0.0); // ... for numOfEntries ...
				vector<double> statsVec_numOfBoth(3, 0.0); // ... and for both together
#endif

				// In case read IDs are relevant
				if (vReadIDtoGenID.size()) {

					addToMatchedReadID = [](vector<uint64_t>& vReadIDs, uint64_t& position, const uint64_t& value) {
						if (vReadIDs.size() <= position) {
							vReadIDs.resize(position + 100);
						}

						vReadIDs[position] = value;
						position++;
					};

				}
				debugBarrier
				unique_ptr<vector<uint64_t>[]> vReadIDs(new vector<uint64_t>[Base::_iNumOfK]);
				vector<uint64_t> vPositions(Base::_iNumOfK);
				vector<intType> vMemoryOfSeenkMers(Base::_iNumOfK);
				vector<Utilities::sBitArray> vMemoryOfTaxIDs(Base::_iNumOfK, Utilities::sBitArray(iSpecIDRange));
				const int32_t& ikDifferenceTop = HIGHESTPOSSIBLEK - Base::_iMaxK; // Highest possible k and length of arrWeightingFactors
				const auto libBeginIt = getVec(vLib, iThreadID)->cbegin();
				auto seenResultIt = libBeginIt;
				auto rangeBeginIt = libBeginIt, rangeEndIt = libBeginIt;

				tuple<intType, uint64_t> iSeenInput = make_pair(0, 0);
				int16_t ikLengthCounter = static_cast<int16_t>(Base::_iNumOfK - 1);
				int32_t shift = 5 * (Base::_iHighestK - Base::_aOfK[ikLengthCounter]);
				const auto shiftVal = [&shift](const intType& val) { return val >> shift; };

				auto vInIndex = vInStart;
				while (vInIndex < vInEnd) {
					// determining the range once is better than checking everytime if its still the same
					uint64_t iSeenRange = vIn.getRangeStart(vInIndex);
					uint32_t iRangeLength = vIn.getRangeEnd(vInIndex);
					uint64_t iInputRangeIdx = vInIndex;
					debugBarrier
					if (iSeenRange == numeric_limits<uint64_t>::max()) {
						++vInIndex;
						continue;
					}

					while (vInIndex < vInEnd) {
						const uint64_t& currentRange = vIn.getRangeStart(vInIndex);
						if (currentRange != iSeenRange && currentRange != numeric_limits<uint64_t>::max()) {
							break;
						}
						else {
							++vInIndex;
						}
					}

					// reset stuff
					for (int32_t i = 0; i < Base::_iNumOfK; ++i) {
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
					debugBarrier
					for (; iInputRangeIdx < vInIndex; ++iInputRangeIdx) {
						// ignore missmatches
						if (vIn.getRangeStart(iInputRangeIdx) == numeric_limits<uint64_t>::max()) {
							continue;
						}

						// reset local stuff
						ikLengthCounter = static_cast<int16_t>(Base::_iNumOfK - 1);
						shift = 5 * (Base::_iHighestK - Base::_aOfK[ikLengthCounter]);

						// get kmer
						const pair<intType, uint32_t>& iCurrentkMer = make_pair(vIn.getkMer(iInputRangeIdx), vIn.getReadID(iInputRangeIdx));
						auto iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;
						bInputIterated = true;

						//cout << kASA::kMerToAminoacid(iCurrentkMer.first, 25) << " " << kASA::kMerToAminoacid(shiftVal(iCurrentkMer.first), 25) << " " << kASA::kMerToAminoacid(iCurrentkMerShifted, 25) << " " << kASA::kMerToAminoacid(shiftVal(rangeBeginIt->first), 25) << " " << kASA::kMerToAminoacid(iCurrentkMer.first, 25) << endl;
						debugBarrier
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
										debugBarrier
										seenResultIt = lower_bound(rangeBeginIt, rangeEndIt + 1, iCurrentkMerShifted, [&shift](const decltype(*libBeginIt)& a, const decltype(iCurrentkMerShifted)& val) { return (a.first >> shift) < val; });
										//cout << get<0>(iCurrentkMer) << " " << iCurrentkMerShifted << " " << (rangeBeginIt - libBeginIt) << " " << (rangeEndIt - libBeginIt) << " " << (seenResultIt - libBeginIt) << " " << seenResultIt - rangeBeginIt << endl;
									}
								}
							}
						}
						bDetermineBeginForMatching = false;
						debugBarrier

						// Now for the real part: Trying to match the kmer with those from the index

						// If the ending is ^ it's not going to hit anyway, might as well stop here
						if ((iCurrentkMerShifted & 31) == 30) {
							continue;
						}

						// Count duplicates that matched too
						if ((get<0>(iSeenInput) == get<0>(iCurrentkMer)) || (seenResultIt == rangeEndIt + 1)) {
							for (int32_t ik = Base::_iNumOfK - 1; ik > -1; --ik) {
								const int32_t& shift_ = 5 * (Base::_iHighestK - Base::_aOfK[ik]);
								const auto& iCurrentkMerShifted_ = get<0>(iCurrentkMer) >> shift_;
								if (compare(iCurrentkMerShifted_, vMemoryOfSeenkMers[ik], Base::_aOfK[ik]) == 1) { //  iCurrentkMerShifted_ == vMemoryOfSeenkMers[ik]
									addToMatchedReadID(vReadIDs[ik], vPositions[ik], get<1>(iCurrentkMer));
									if (GlobalInputParameters.bPostProcess) {
										vIn.setMatchLength(iInputRangeIdx, Base::_aOfK[ik]);
									}
								}
							}
							continue;
						}
						else {
							iSeenInput = iCurrentkMer;
						}


						debugBarrier
						bool bBreakOut = false;
						while (seenResultIt != rangeEndIt + 1 && !bBreakOut) {
							const tuple<intType, uint32_t>& iCurrentLib = make_tuple(seenResultIt->first, seenResultIt->second);

							ikLengthCounter = static_cast<int16_t>(Base::_iNumOfK - 1);
							for (; ikLengthCounter >= 0; --ikLengthCounter) {

								shift = 5 * (Base::_iHighestK - Base::_aOfK[ikLengthCounter]);
								iCurrentkMerShifted = get<0>(iCurrentkMer) >> shift;

								const auto& iCurrentLibkMerShifted = get<0>(iCurrentLib) >> shift;

								//cout << Base::kMerToAminoacid(iCurrentkMerShifted, 25) << " " << Base::kMerToAminoacid(iCurrentLibkMerShifted, 25) << endl;

								const uint8_t iResultFromComparison = compare(iCurrentkMerShifted, iCurrentLibkMerShifted, Base::_aOfK[ikLengthCounter]);
								if (iResultFromComparison == 0) {
									if (bInputIterated) {
										for (int32_t ik = ikLengthCounter; ik > -1; --ik) {
											const int32_t& shift_ = 5 * (Base::_iHighestK - Base::_aOfK[ik]);
											const auto& iCurrentkMerShifted_ = get<0>(iCurrentkMer) >> shift_;
											if (compare(iCurrentkMerShifted_, vMemoryOfSeenkMers[ik], Base::_aOfK[ik]) == 1) { //  iCurrentkMerShifted_ == vMemoryOfSeenkMers[ik]
												addToMatchedReadID(vReadIDs[ik], vPositions[ik], get<1>(iCurrentkMer));
												if (GlobalInputParameters.bPostProcess) {
													vIn.setMatchLength(iInputRangeIdx, Base::_aOfK[ik]);
												}
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
									if (iResultFromComparison == 1) {
										// If the ending is ^, it's not a valid match, might as well stop here
										if ((iCurrentkMerShifted & 31) == 30) {
											bBreakOut = true;
											break;
										}

										if (this->_bVisualize) {
											_matchedkMers.push_back(make_pair(Base::kMerToAminoacid(iCurrentLibkMerShifted, Base::_aOfK[ikLengthCounter]), get<1>(iCurrentLib)));
										}

										// Delayed scoring: First gather everything and then score it instead of scoring everytime you encounter a new read or tax id.
										if (iCurrentkMerShifted == vMemoryOfSeenkMers[ikLengthCounter]) {
											// We've seen that already, just add it. 
											markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));
											if (bInputIterated) {
												addToMatchedReadID(vReadIDs[ikLengthCounter], vPositions[ikLengthCounter], get<1>(iCurrentkMer));
												if (GlobalInputParameters.bPostProcess) {
													vIn.setMatchLength(iInputRangeIdx, Base::_aOfK[ikLengthCounter]);
												}
											}
										}
										else {
											debugBarrier
												// For this k, the kmer is different. Save the gathered information.
											const auto& numOfEntries = vMemoryOfTaxIDs[ikLengthCounter].numOfEntries();
											const auto& numOfHits = vPositions[ikLengthCounter];
											const uint64_t& iPartialTempIndex = uint64_t(iSpecIDRange) * Base::_iNumOfK * (iThreadID)+iSpecIDRange * uint64_t(ikLengthCounter);
											const auto& weight = arrWeightingFactors[ikDifferenceTop + ikLengthCounter];
											const auto& score = weight * (1.f / numOfEntries);
											const auto& counts = double(numOfHits) / numOfEntries;
											

#ifdef STATISTICS
											statsVec_numOfHits[0] = min(statsVec_numOfHits[0], static_cast<double>(numOfHits));
											statsVec_numOfHits[1] = max(statsVec_numOfHits[1], static_cast<double>(numOfHits));
											statsVec_numOfHits[2] = statsVec_numOfHits[2] + (numOfHits - statsVec_numOfHits[2]) / (iterationVar + 1);

											statsVec_numOfEntries[0] = min(statsVec_numOfEntries[0], static_cast<double>(numOfEntries));
											statsVec_numOfEntries[1] = max(statsVec_numOfEntries[1], static_cast<double>(numOfEntries));
											statsVec_numOfEntries[2] = statsVec_numOfEntries[2] + (numOfEntries - statsVec_numOfEntries[2]) / (iterationVar + 1);

											statsVec_numOfBoth[0] = min(statsVec_numOfBoth[0], static_cast<double>(numOfHits * numOfEntries));
											statsVec_numOfBoth[1] = max(statsVec_numOfBoth[1], static_cast<double>(numOfHits * numOfEntries));
											statsVec_numOfBoth[2] = statsVec_numOfBoth[2] + (numOfHits * numOfEntries - statsVec_numOfBoth[2]) / (iterationVar + 1);

											iterationVar++;
#endif

											scoreMatchForReadIDsAndTaxIDs(vReadIDtoGenID, vCount, vCountUnique, vCountTotal, countTotal, tmpStore, aOfIndices, tmpStore2, aOfTempIdx, vMemoryOfTaxIDs[ikLengthCounter], numOfEntries, iPartialTempIndex, vReadIDs[ikLengthCounter], numOfHits, counts, score);

											vPositions[ikLengthCounter] = 0;
											addToMatchedReadID(vReadIDs[ikLengthCounter], vPositions[ikLengthCounter], get<1>(iCurrentkMer));
											if (GlobalInputParameters.bPostProcess) { vIn.setMatchLength(iInputRangeIdx, Base::_aOfK[ikLengthCounter]); }

											vMemoryOfTaxIDs[ikLengthCounter].clear();
											markTaxIDs(get<1>(iCurrentLib), vMemoryOfTaxIDs[ikLengthCounter], mTaxToIdx, getVec(vLib, iThreadID));

											vMemoryOfSeenkMers[ikLengthCounter] = iCurrentkMerShifted;
											debugBarrier
										}
									}
									else {
										// index kmer is smaller than input kmer. Iterate linearly and score everything that also matched before (albeit with a smaller k)
										uint64_t iTempCounter = 1;
										while (seenResultIt + iTempCounter != rangeEndIt + 1) {
											const intType& iNextLibSuffix = static_cast<intType>((seenResultIt + iTempCounter)->first);
											//cout << kASA::kMerToAminoacid(iCurrentkMerShifted, 12) << " " << kASA::kMerToAminoacid(iNextLibSuffix, 12) << endl;
											if (iCurrentkMerShifted > (iNextLibSuffix >> shift)) {
												int16_t iUntilK = static_cast<int16_t>(Base::_iNumOfK - 1);
												for (; iUntilK > -1; --iUntilK) {
													if (compare(vMemoryOfSeenkMers[iUntilK], (iNextLibSuffix >> 5 * (Base::_iHighestK - Base::_aOfK[iUntilK])), Base::_aOfK[iUntilK]) == 1) { //  vMemoryOfSeenkMers[iUntilK] == (iNextLibSuffix >> 5 * (Base::_iHighestK - Base::_aOfK[iUntilK]))
														//if (this->_bVisualize) {
															//_matchedkMers.push_back(make_pair(Base::kMerToAminoacid((iNextLibSuffix >> 5 * (Base::_iHighestK - Base::_aOfK[iUntilK])), Base::_aOfK[iUntilK]), (seenResultIt + iTempCounter)->second));
														//}
														markTaxIDs((seenResultIt + iTempCounter)->second, vMemoryOfTaxIDs[iUntilK], mTaxToIdx, getVec(vLib, iThreadID));
													}
													else {
														break;
													}
												}
												if (iUntilK < static_cast<int16_t>(Base::_iNumOfK - 1)) {
													++iTempCounter;
												}
												else {
													auto nextValuablekMer = lower_bound(seenResultIt + iTempCounter, rangeEndIt + 1, iCurrentkMerShifted, [&shift](const decltype(*libBeginIt)& a, const decltype(iCurrentkMerShifted)& val) { return (a.first >> shift) < val; });
													//cout << kASA::kMerToAminoacid((nextValuablekMer - 1)->first, 12) << " " << kASA::kMerToAminoacid(nextValuablekMer->first, 12) << " " << nextValuablekMer - seenResultIt << endl;
													iTempCounter = nextValuablekMer - seenResultIt;
													break;
												}
											}
											else {
												break;
											}
										}
										seenResultIt += iTempCounter;
										debugBarrier
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
					debugBarrier
					uint64_t iTempCounter = 0;
					while (seenResultIt + iTempCounter != rangeEndIt + 1) {
						const intType& iNextLibSuffix = static_cast<intType>((seenResultIt + iTempCounter)->first);
						int16_t iUntilK = static_cast<int16_t>(Base::_iNumOfK - 1);
						for (; iUntilK > -1; --iUntilK) {
							if (compare(vMemoryOfSeenkMers[iUntilK], (iNextLibSuffix >> 5 * (Base::_iHighestK - Base::_aOfK[iUntilK])), Base::_aOfK[iUntilK]) == 1) { //   vMemoryOfSeenkMers[iUntilK] == (iNextLibSuffix >> 5 * (Base::_iHighestK - Base::_aOfK[iUntilK]))
								//if (this->_bVisualize) {
									//_matchedkMers.push_back(make_pair(Base::kMerToAminoacid((iNextLibSuffix >> 5 * (Base::_iHighestK - Base::_aOfK[iUntilK])), Base::_aOfK[iUntilK]), (seenResultIt + iTempCounter)->second));
								//}
								markTaxIDs((seenResultIt + iTempCounter)->second, vMemoryOfTaxIDs[iUntilK], mTaxToIdx, getVec(vLib, iThreadID));
							}
							else {
								break;
							}
						}
						if (iUntilK < static_cast<int16_t>(Base::_iNumOfK - 1)) {
							++iTempCounter;
						}
						else {
							break;
						}
					}

					// save everything from this range
					debugBarrier
					for (int16_t ikLC = static_cast<int16_t>(Base::_iNumOfK - 1); ikLC >= 0; --ikLC) {

						const auto& numOfEntries = vMemoryOfTaxIDs[ikLC].numOfEntries();
						const auto& numOfHits = vPositions[ikLC];
						const uint64_t& iPartialTempIndex = uint64_t(iSpecIDRange) * Base::_iNumOfK * (iThreadID) + uint64_t(iSpecIDRange) * ikLC;
						const auto& weight = arrWeightingFactors[ikDifferenceTop + ikLC];
						const auto& score = weight * (1.f / numOfEntries);
						const auto& counts = double(numOfHits) / numOfEntries;
						scoreMatchForReadIDsAndTaxIDs(vReadIDtoGenID, vCount, vCountUnique, vCountTotal, countTotal, tmpStore, aOfIndices, tmpStore2, aOfTempIdx, vMemoryOfTaxIDs[ikLC], numOfEntries, iPartialTempIndex, vReadIDs[ikLC], numOfHits, counts, score);
					}
					debugBarrier
				}
				debugBarrier
					//m_exceptionLock.lock();
					//cout << "Loops: " << dAverageLoopCount / iNumOfLoops  << " " << dAverageLoopCount << " " << iNumOfLoops << " " << vInEnd - vInStart << endl;
					//m_exceptionLock.unlock();

#ifdef STATISTICS
					m_exceptionLock.lock();
					cout << "STATISTICS:" << endl;
					cout << "NumOfHits:\n Min: " << statsVec_numOfHits[0] << "\n Max: " << statsVec_numOfHits[1] << "\n Avg: " << statsVec_numOfHits[2] << endl;
					cout << "NumOfEntries:\n Min: " << statsVec_numOfEntries[0] << "\n Max: " << statsVec_numOfEntries[1] << "\n Avg: " << statsVec_numOfEntries[2] << endl;
					cout << "NumOfBoth:\n Min: " << statsVec_numOfBoth[0] << "\n Max: " << statsVec_numOfBoth[1] << "\n Avg: " << statsVec_numOfBoth[2] << endl;
					cout << "Iterations: " << iterationVar << endl;
					m_exceptionLock.unlock();
#endif // STATISTICS

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
		// Sort input (in parallel)
		inline void sortInputAndCheckInvalidkMers_sta(InputType<intType>& vInputVec, const bool& bPartitioned, const Trie<intType>& T, const int32_t& iLocalNumOfThreads) {

			if (iLocalNumOfThreads == 1) {
				std::sort(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [](const typename InputType<intType>::staTuple& p1, const typename InputType<intType>::staTuple& p2) {
					return get<1>(p1) < get<1>(p2);
					});
				

				if (bPartitioned) {
					for_each(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [this, &T](typename InputType<intType>::staTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
					for_each(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [this, &T](typename InputType<intType>::staTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						uint64_t kMerForTrie = 0;
						if (is_same<intType, uint128_t>::value) {
							kMerForTrie = static_cast<uint64_t>((static_cast<uint128_t>(get<1>(a)) >> 95));
						}
						else {
							kMerForTrie = (get<1>(a) >> 30);
						}
						//cout << kASA::kMerToAminoacid(get<1>(a), 25) << " " << kASA::kMerToAminoacid(kMerForTrie, 12) << endl;
						T.GetIndexRange(kMerForTrie, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
						if (start != numeric_limits<uint64_t>::max()) {
							get<0>(a) = start;
							get<2>(a) = range;
						}
						else {
							get<0>(a) = start;
						}
						});
				}
			} 
			else {


# if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
				// Some context:
				// The gnu parallel quicksort implementation only uses two cores whereas the gnu parallel merge uses all but isn't in-place. This significantly worsened performance in Linux environments.
				// Futhermore, the gcc 9.* compiler uses Threadblocks for its stl conform parallel implementation as of now (2020). This is unacceptable for kASA.
				// Therefore I set out to find a better one and found the preliminary implementation for C++17 inside Visual Studio which was published under the Apache Software License 2.0.
				// This is compatible with kASA and as near to the current implementation inside Visual Studio 2017 as can be without risking copyright infringement towards Microsoft.
				// It can be found inside the ParallelQuicksort.hpp header.
				Utilities::parallelQuicksort(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [](const typename InputType<intType>::staTuple& p1, const typename InputType<intType>::staTuple& p2) {
					return get<1>(p1) < get<1>(p2);
					}, iLocalNumOfThreads);
				
				if (bPartitioned) {
					__gnu_parallel::for_each(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [this, &T](typename InputType<intType>::staTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
					__gnu_parallel::for_each(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [this, &T](typename InputType<intType>::staTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						uint64_t kMerForTrie = 0;
						if (is_same<intType, uint128_t>::value) {
							kMerForTrie = static_cast<uint64_t>((static_cast<uint128_t>(get<1>(a)) >> 95));
						}
						else {
							kMerForTrie = (get<1>(a) >> 30);
						}
						T.GetIndexRange(kMerForTrie, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
#if __has_include(<execution>) && !defined(__llvm__)


			std::sort(std::execution::par_unseq, vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [](const typename InputType<intType>::staTuple& p1, const typename InputType<intType>::staTuple& p2) {
					return get<1>(p1) < get<1>(p2);
			});


				if (bPartitioned) {
					for_each(std::execution::par_unseq, vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [this, &T](typename InputType<intType>::staTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
					for_each(std::execution::par_unseq, vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [this, &T](typename InputType<intType>::staTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						uint64_t kMerForTrie = 0;
						if (is_same<intType, uint128_t>::value) {
							kMerForTrie = static_cast<uint64_t>((static_cast<uint128_t>(get<1>(a)) >> 95));
						}
						else {
							kMerForTrie = (get<1>(a) >> 30);
						}
						T.GetIndexRange(kMerForTrie, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
				
				Utilities::parallelQuicksort(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [](const typename InputType<intType>::staTuple& p1, const typename InputType<intType>::staTuple& p2) {
					return get<1>(p1) < get<1>(p2);
					}, iLocalNumOfThreads);
				
				if (bPartitioned) {
					for_each(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [this, &T](typename InputType<intType>::staTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
					for_each(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), [this, &T](typename InputType<intType>::staTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						uint64_t kMerForTrie = 0;
						if (is_same<intType, uint128_t>::value) {
							kMerForTrie = static_cast<uint64_t>((static_cast<uint128_t>(get<1>(a)) >> 95));
						}
						else {
							kMerForTrie = (get<1>(a) >> 30);
						}
						T.GetIndexRange(kMerForTrie, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
			}
		}

		inline void sortInputAndCheckInvalidkMers_pp(InputType<intType>& vInputVec, const bool& bPartitioned, const Trie<intType>& T, const int32_t& iLocalNumOfThreads) {

			if (iLocalNumOfThreads == 1) {
				std::sort(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [](const typename InputType<intType>::ppTuple& p1, const typename InputType<intType>::ppTuple& p2) {
					return get<1>(p1) < get<1>(p2);
					});


				if (bPartitioned) {
					for_each(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [this, &T](typename InputType<intType>::ppTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
					for_each(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [this, &T](typename InputType<intType>::ppTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						uint64_t kMerForTrie = 0;
						if (is_same<intType, uint128_t>::value) {
							kMerForTrie = static_cast<uint64_t>((static_cast<uint128_t>(get<1>(a)) >> 95));
						}
						else {
							kMerForTrie = (get<1>(a) >> 30);
						}
						//cout << kASA::kMerToAminoacid(get<1>(a), 25) << " " << kASA::kMerToAminoacid(kMerForTrie, 12) << endl;
						T.GetIndexRange(kMerForTrie, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
						if (start != numeric_limits<uint64_t>::max()) {
							get<0>(a) = start;
							get<2>(a) = range;
						}
						else {
							get<0>(a) = start;
						}
						});
				}
			}
			else {


# if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
				// Some context:
				// The gnu parallel quicksort implementation only uses two cores whereas the gnu parallel merge uses all but isn't in-place. This significantly worsened performance in Linux environments.
				// Futhermore, the gcc 9.* compiler uses Threadblocks for its stl conform parallel implementation as of now (2020). This is unacceptable for kASA.
				// Therefore I set out to find a better one and found the preliminary implementation for C++17 inside Visual Studio which was published under the Apache Software License 2.0.
				// This is compatible with kASA and as near to the current implementation inside Visual Studio 2017 as can be without risking copyright infringement towards Microsoft.
				// It can be found inside the ParallelQuicksort.hpp header.
				Utilities::parallelQuicksort(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [](const typename InputType<intType>::ppTuple& p1, const typename InputType<intType>::ppTuple& p2) {
					return get<1>(p1) < get<1>(p2);
					}, iLocalNumOfThreads);

				if (bPartitioned) {
					__gnu_parallel::for_each(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [this, &T](typename InputType<intType>::ppTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
					__gnu_parallel::for_each(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [this, &T](typename InputType<intType>::ppTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						uint64_t kMerForTrie = 0;
						if (is_same<intType, uint128_t>::value) {
							kMerForTrie = static_cast<uint64_t>((static_cast<uint128_t>(get<1>(a)) >> 95));
						}
						else {
							kMerForTrie = (get<1>(a) >> 30);
						}
						T.GetIndexRange(kMerForTrie, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
#if __has_include(<execution>) && !defined(__llvm__)


				std::sort(std::execution::par_unseq, vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [](const typename InputType<intType>::ppTuple& p1, const typename InputType<intType>::ppTuple& p2) {
					return get<1>(p1) < get<1>(p2);
					});


				if (bPartitioned) {
					for_each(std::execution::par_unseq, vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [this, &T](typename InputType<intType>::ppTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
					for_each(std::execution::par_unseq, vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [this, &T](typename InputType<intType>::ppTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						uint64_t kMerForTrie = 0;
						if (is_same<intType, uint128_t>::value) {
							kMerForTrie = static_cast<uint64_t>((static_cast<uint128_t>(get<1>(a)) >> 95));
						}
						else {
							kMerForTrie = (get<1>(a) >> 30);
						}
						T.GetIndexRange(kMerForTrie, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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

				Utilities::parallelQuicksort(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [](const typename InputType<intType>::ppTuple& p1, const typename InputType<intType>::ppTuple& p2) {
					return get<1>(p1) < get<1>(p2);
					}, iLocalNumOfThreads);

				if (bPartitioned) {
					for_each(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [this, &T](typename InputType<intType>::ppTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						T.GetIndexRange(get<1>(a) >> 30, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
					for_each(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), [this, &T](typename InputType<intType>::ppTuple& a) {
						uint64_t start = 0;
						uint32_t range = 0;
						uint64_t kMerForTrie = 0;
						if (is_same<intType, uint128_t>::value) {
							kMerForTrie = static_cast<uint64_t>((static_cast<uint128_t>(get<1>(a)) >> 95));
						}
						else {
							kMerForTrie = (get<1>(a) >> 30);
						}
						T.GetIndexRange(kMerForTrie, static_cast<int8_t>((Base::_iMinK > 6) ? 6 : Base::_iMinK), move(start), move(range));
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
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// this score is needed for the error score
		inline float calculateBestScore(const size_t& iReadLength) {
			float bestScore = 0.f;

			if (Base::_bProtein) {
				for (int32_t i = Base::_iMinK; i <= Base::_iMaxK; ++i) {
					bestScore += (iReadLength - i + 1) * arrWeightingFactors[HIGHESTPOSSIBLEK - i];
				}
			}
			else {
				switch (Base::_iNumOfFrames) {
				case 1:
					for (int32_t i = Base::_iMinK; i <= Base::_iMaxK; ++i) {
						bestScore += (iReadLength/3 - i + 1) * arrWeightingFactors[HIGHESTPOSSIBLEK - i];
					}
					break;
				case 3:
					for (int32_t i = Base::_iMinK; i <= Base::_iMaxK; ++i) {
						bestScore += (iReadLength - i * 3 + 1) * arrWeightingFactors[HIGHESTPOSSIBLEK - i];

					}
					break;
				case 6:
					for (int32_t i = Base::_iMinK; i <= Base::_iMaxK; ++i) {
						bestScore += 2 * (iReadLength - i * 3 + 1) * arrWeightingFactors[HIGHESTPOSSIBLEK - i];
					}
					break;
				}
			}
			return bestScore;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Score all of them
		inline void scoringFunc(const Utilities::Non_contiguousArray& vReadIDtoGenID, const uint64_t& iStart, const uint64_t& iEnd, const uint64_t& iRealReadIDStart, list<pair<string, readIDType>>& vReadNameAndLength, const uint64_t& iNumberOfSpecies, const vector<vector<uint64_t>>& mFrequencies, const vector<uint32_t>& mIdxToTax, const vector<string>& mOrganisms, const float& fThreshold, vector<uint64_t>& vReadIDsToBeFiltered, const vector<float>& vCoherence, ofstream&& fOut) {
			try {
				debugBarrier
				vector<tuple<size_t, float, double>> resultVec(iNumberOfSpecies);
				int64_t iOldCountOfHits = 0;
				Utilities::BufferedWriter outStreamer(fOut, 104857600ull);

				for (uint64_t readIdx = iStart; readIdx < iEnd; ++readIdx) {
					auto currentReadLengthAndName = vReadNameAndLength.front();

					float bestScore = calculateBestScore(currentReadLengthAndName.second);

					debugBarrier
					float kMerScore = 0.f;
					double relativeScore = 0.0;
					int64_t iCountOfHits = 0;
					for (size_t iSpecIdx = 1; iSpecIdx < iNumberOfSpecies; ++iSpecIdx) {
						if (vReadIDtoGenID[readIdx][iSpecIdx] > 0.f) {
							debugBarrier
							kMerScore = vReadIDtoGenID[readIdx][iSpecIdx];
							debugBarrier
							if (Base::_bProtein) {
								relativeScore = kMerScore / (1.0 + log2(mFrequencies[iSpecIdx][0] * double(currentReadLengthAndName.second - Base::_iHighestK + 1)));
							}
							else {
								relativeScore = kMerScore / (1.0 + log2(mFrequencies[iSpecIdx][0] * double(currentReadLengthAndName.second - Base::_iHighestK * 3 + 1)));
							}
							debugBarrier
								
							if (relativeScore >= fThreshold) {
								get<0>(resultVec[iCountOfHits]) = iSpecIdx;
								get<1>(resultVec[iCountOfHits]) = kMerScore;
								get<2>(resultVec[iCountOfHits]) = relativeScore;
								++iCountOfHits;
							}
							debugBarrier
						}
					}

					debugBarrier
					if (iCountOfHits == 0) {
						switch (format) {
						case OutputFormat::tsv:
							Utilities::itostr(iRealReadIDStart + readIdx, outStreamer.getString());
							outStreamer += "\t";
							outStreamer += currentReadLengthAndName.first;
							outStreamer += "\t-\t-\t-\t-";
							if (GlobalInputParameters.bPostProcess) {
								outStreamer += "\t-";
							}
							outStreamer += "\n";
							break;

						case OutputFormat::Json:
							if (iRealReadIDStart + readIdx == 0) {
								outStreamer += "{\n";
							}
							else {
								outStreamer += ",\n{\n";
							}
							outStreamer += "\t\"Read number\": ";
							Utilities::itostr(iRealReadIDStart + readIdx, outStreamer.getString());
							outStreamer += ",\n\t\"Specifier from input file\": \"";
							outStreamer += currentReadLengthAndName.first;
							outStreamer += "\",\n\t\"Length\": ";
							Utilities::itostr(currentReadLengthAndName.second, outStreamer.getString());
							outStreamer += ",\n\t\"Top hits\": [\n\t],\n\t\"Further hits\": [\n\t]\n}";
							break;

						case OutputFormat::JsonL:
							outStreamer += "{ \"Read number\": ";
							Utilities::itostr(iRealReadIDStart + readIdx, outStreamer.getString());
							outStreamer += ", \"Specifier from input file\": \"";
							outStreamer += currentReadLengthAndName.first;
							outStreamer += "\", \"Length\": ";
							Utilities::itostr(currentReadLengthAndName.second, outStreamer.getString());
							outStreamer += ", \"Top hits\": [], \"Further hits\": [] }\n";
							break;

						case OutputFormat::Kraken:
							outStreamer += "U\t";
							outStreamer += currentReadLengthAndName.first;
							outStreamer += "\t0\t";
							outStreamer += currentReadLengthAndName.second;
							outStreamer += "\tA:00\n";
							break;
						};

					}
					else {
						debugBarrier
						if (iCountOfHits < iOldCountOfHits) {
							fill(resultVec.begin() + iCountOfHits, resultVec.begin() + iOldCountOfHits, tuple<size_t, float, double>(0, 0.f, 0.0));
						}
						iOldCountOfHits = iCountOfHits;


						sort(resultVec.begin(), resultVec.begin() + iCountOfHits, [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<2>(a) > get<2>(b); });

						debugBarrier
						auto maxValue = get<1>(*(max_element(resultVec.begin(), resultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<1>(a) < get<1>(b); })));
						int64_t iTopHitCounter = 1;
						for (int64_t i = 1; i < iCountOfHits && i < iNumOfBeasts; ++i) {
							if ((get<1>(resultVec[i])) / (maxValue) > 0.8f) { // determine "close enough" with normalization relative to the highest score
								++iTopHitCounter;
							}
							else {
								break;
							}
						}

						// Calculate error and if Error is too high, filter it out later
						if (GlobalInputParameters.bFilter) {
							if ((bestScore - double(maxValue)) / bestScore < GlobalInputParameters.fErrorThreshold) {
								vReadIDsToBeFiltered.push_back(iRealReadIDStart + readIdx);
							}
							else {
								if (GlobalInputParameters.bPostProcess) {
									if (vCoherence[readIdx] >= GlobalInputParameters.fCoherenceThreshold) {
										vReadIDsToBeFiltered.push_back(iRealReadIDStart + readIdx);
									}
								}
							}
						}

						string sOut = "", sOut2 = "", sOut3 = "", sOut4 = "";
						auto it = resultVec.begin();
						float iValueBefore = 0;

						debugBarrier
						switch (format) {
						case OutputFormat::tsv:
							// tsv
							// += is faster than () + ()
							Utilities::itostr(iRealReadIDStart + readIdx, sOut);
							sOut += "\t";
							sOut += currentReadLengthAndName.first + "\t";

							for (int64_t j = 0, i = 0; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
								sOut += ";";
								sOut2 += mOrganisms[get<0>(*it)];
								sOut2 += ";";

								dtoa_milo(get<2>(*it), sOut3);
								sOut3 += ",";
								dtoa_milo(get<1>(*it), sOut3);
								sOut3 += ";";

								dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut4);
								sOut4 += ";";

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
							if (sOut4.back() == ';') {
								sOut4.pop_back();
							}
							if (sOut2.length()) {
								outStreamer += sOut;
								outStreamer += "\t";
								outStreamer += sOut2;
								outStreamer += "\t";
								outStreamer += sOut3;
								outStreamer += "\t";
								outStreamer += sOut4;
								if (GlobalInputParameters.bPostProcess) {
									outStreamer += "\t";
									dtoa_milo(vCoherence[readIdx], outStreamer.getString());
								}
								outStreamer += "\n";
							}
							debugBarrier
							break;

						case OutputFormat::Json:
							// json
					 
							if (iRealReadIDStart + readIdx == 0) {
								outStreamer += "{\n";
							}
							else {
								outStreamer += ",\n{\n";
							}

							outStreamer += "\t\"Read number\": ";
							Utilities::itostr(iRealReadIDStart + readIdx, outStreamer.getString());
							outStreamer += ",\n\t\"Specifier from input file\": \"";
							outStreamer += currentReadLengthAndName.first;
							outStreamer += "\",\n\t\"Length\": ";
							Utilities::itostr(currentReadLengthAndName.second, outStreamer.getString());
							outStreamer += ",\n\t\"Top hits\": [\n";

							for (int64_t i = 0; i < iTopHitCounter; ++i, ++it) {
								if (i == 0) {
									outStreamer += "\t{\n";
								}
								else {
									outStreamer += ",\n\t{\n";
								}

								outStreamer += "\t\t\"tax ID\": \"";
								Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
								outStreamer += "\",\n";
								outStreamer += "\t\t\"Name\": \"";
								outStreamer += mOrganisms[get<0>(*it)];
								outStreamer += "\",\n";
								outStreamer += "\t\t\"k-mer Score\": ";
								dtoa_milo(get<1>(*it), outStreamer.getString());
								outStreamer += ",\n";
								outStreamer += "\t\t\"Relative Score\": ";
								dtoa_milo(get<2>(*it), outStreamer.getString());
								outStreamer += ",\n";
								outStreamer += "\t\t\"Error\": ";
								dtoa_milo((bestScore - get<1>(*it)) / bestScore, outStreamer.getString());
								if (GlobalInputParameters.bPostProcess) {
									outStreamer += ",\n";
									outStreamer += "\t\t\"Coherence\": ";
									dtoa_milo(vCoherence[readIdx], outStreamer.getString());
								}
								outStreamer += "\n\t}";
							}

							outStreamer += "\n\t],\n\t\"Further hits\": [\n";

							for (int64_t j = iTopHitCounter, i = iTopHitCounter; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								if (j == iTopHitCounter) {
									outStreamer += "\t{\n";
								}
								else {
									outStreamer += ",\n\t{\n";
								}

								outStreamer += "\t\t\"tax ID\": \"";
								Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
								outStreamer += "\",\n";
								outStreamer += "\t\t\"Name\": \"";
								outStreamer += mOrganisms[get<0>(*it)];
								outStreamer += "\",\n";
								outStreamer += "\t\t\"k-mer Score\": ";
								dtoa_milo(get<1>(*it), outStreamer.getString());
								outStreamer += ",\n";
								outStreamer += "\t\t\"Relative Score\": ";
								dtoa_milo(get<2>(*it), outStreamer.getString());
								outStreamer += ",\n";
								outStreamer += "\t\t\"Error\": ";
								dtoa_milo((bestScore - get<1>(*it)) / bestScore, outStreamer.getString());
								if (GlobalInputParameters.bPostProcess) {
									outStreamer += ",\n";
									outStreamer += "\t\t\"Coherence\": ";
									dtoa_milo(vCoherence[readIdx], outStreamer.getString());
								}
								outStreamer += "\n\t}";

								if (iValueBefore != get<1>(*it)) {
									iValueBefore = get<1>(*it);
									++j;
								}
							}

							debugBarrier
							outStreamer += "\n\t]\n}";
							break;

						case OutputFormat::JsonL:
							// json lines
							outStreamer += "{ \"Read number\": ";
							Utilities::itostr(iRealReadIDStart + readIdx, outStreamer.getString());
							outStreamer += ", \"Specifier from input file\": \"";
							outStreamer += currentReadLengthAndName.first;
							outStreamer += "\", \"Length\": ";
							Utilities::itostr(currentReadLengthAndName.second, outStreamer.getString());
							outStreamer += ", \"Top hits\": [";

							for (int64_t i = 0; i < iTopHitCounter; ++i, ++it) {
								if (i == 0) {
									outStreamer += "{";
								}
								else {
									outStreamer += ",{";
								}

								outStreamer += " \"tax ID\": \"";
								Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
								outStreamer += "\",";
								outStreamer += " \"Name\": \"";
								outStreamer += mOrganisms[get<0>(*it)];
								outStreamer += "\",";
								outStreamer += " \"k-mer Score\": ";
								dtoa_milo(get<1>(*it), outStreamer.getString());
								outStreamer += ",";
								outStreamer += " \"Relative Score\": ";
								dtoa_milo(get<2>(*it), outStreamer.getString());
								outStreamer += ",";
								outStreamer += " \"Error\": ";
								dtoa_milo((bestScore - get<1>(*it)) / bestScore, outStreamer.getString());
								if (GlobalInputParameters.bPostProcess) {
									outStreamer += ",";
									outStreamer += "\"Coherence\": ";
									dtoa_milo(vCoherence[readIdx], outStreamer.getString());
								}
								outStreamer += "}";
							}

							outStreamer += "], \"Further hits\": [";

							for (int64_t j = iTopHitCounter, i = iTopHitCounter; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								if (j == iTopHitCounter) {
									outStreamer += "{";
								}
								else {
									outStreamer += ", {";
								}

								outStreamer += " \"tax ID\": \"";
								Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
								outStreamer += "\",";
								outStreamer += " \"Name\": \"";
								outStreamer += mOrganisms[get<0>(*it)];
								outStreamer += "\",";
								outStreamer += " \"k-mer Score\": "; 
								dtoa_milo(get<1>(*it), outStreamer.getString());
								outStreamer += ",";
								outStreamer += " \"Relative Score\": "; 
								dtoa_milo(get<2>(*it), outStreamer.getString());
								outStreamer += ",";
								outStreamer += " \"Error\": ";
								dtoa_milo((bestScore - get<1>(*it)) / bestScore, outStreamer.getString());
								if (GlobalInputParameters.bPostProcess) {
									outStreamer += ",";
									outStreamer += "\"Coherence\": ";
									dtoa_milo(vCoherence[readIdx], outStreamer.getString());
								}
								outStreamer += "}";

								if (iValueBefore != get<1>(*it)) {
									iValueBefore = get<1>(*it);
									++j;
								}
							}

							debugBarrier
							outStreamer += "] }\n";
							break;

						case OutputFormat::Kraken:
							// Kraken
							outStreamer += "C\t";
							outStreamer += currentReadLengthAndName.first;
							outStreamer += "\t";
							Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
							outStreamer += "\t";
							Utilities::itostr(currentReadLengthAndName.second, outStreamer.getString());
							outStreamer += "\t";

							for (int64_t i = 0; i < iTopHitCounter; ++i, ++it) {
								Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
								outStreamer += ":";
								dtoa_milo(get<1>(*it), outStreamer.getString());
								outStreamer += " ";
							}

							for (int64_t j = iTopHitCounter, i = iTopHitCounter; i < iCountOfHits && j < iNumOfBeasts; ++it, ++i) {
								Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
								outStreamer += ":";
								dtoa_milo(get<1>(*it), outStreamer.getString());
								outStreamer += " ";
								if (iValueBefore != get<1>(*it)) {
									iValueBefore = get<1>(*it);
									++j;
								}
							}

							debugBarrier
							outStreamer += "\n";
							break;
						};

						//fill(resultVec.begin(), resultVec.begin() + iCountOfHits, tuple<size_t, float, double>(0, 0.f, 0.0));
						/*for (auto& entry : resultVec) {
							get<0>(entry) = 0;
							get<1>(entry) = 0.f;
							get<2>(entry) = 0.0;
						}*/
					}

					vReadNameAndLength.pop_front();
				}
				debugBarrier
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl;
				throw;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Score one of them
		inline void scoringFunc(vector<tuple<readIDType, float, double>>&& vTempResultVec, const uint64_t& iReadNum, const pair<string, uint32_t>& vReadNameAndLength, const vector<vector<uint64_t>>& mFrequencies, const vector<uint32_t>& mIdxToTax, const vector<string>& mOrganisms, const float& fThreshold, vector<uint64_t>& vReadIDsToBeFiltered, const float& vCoherence, ofstream&& fOut) {
			try {
				debugBarrier
				Utilities::BufferedWriter outStreamer(fOut, 524288ull);
				float bestScore = calculateBestScore(vReadNameAndLength.second);
				debugBarrier

				for (auto it = vTempResultVec.begin(); it != vTempResultVec.end();) {
					double relativeScore = 0.0;
					if (Base::_bProtein) {
						relativeScore = double(get<1>(*it)) / (1.0 + log2(mFrequencies[get<0>(*it)][0] * double(vReadNameAndLength.second - Base::_iHighestK + 1)));
					}
					else {
						relativeScore = double(get<1>(*it)) / (1.0 + log2(mFrequencies[get<0>(*it)][0] * double(vReadNameAndLength.second - Base::_iHighestK * 3 + 1)));
					}
					if (relativeScore >= fThreshold) {
						get<2>(*it) = relativeScore;
						++it;
					}
					else {
						it = vTempResultVec.erase(it);
					}
				}
				debugBarrier

				if (vTempResultVec.size() == 0) {
					switch (format) {
					case OutputFormat::tsv:
						Utilities::itostr(iReadNum, outStreamer.getString());
						outStreamer += "\t";
						outStreamer += vReadNameAndLength.first;
						outStreamer += "\t-\t-\t-\t-";
						if (GlobalInputParameters.bPostProcess) {
							outStreamer += "\t-";
						}
						outStreamer += "\n";
						break;

					case OutputFormat::Json:
						if (iReadNum == 0) {
							outStreamer += "{\n";
						}
						else {
							outStreamer += ",\n{\n";
						}
						outStreamer += "\t\"Read number\": ";
						Utilities::itostr(iReadNum, outStreamer.getString());
						outStreamer += ",\n\t\"Specifier from input file\": \"";
						outStreamer += vReadNameAndLength.first;
						outStreamer += "\",\n\t\"Length\": ";
						Utilities::itostr(vReadNameAndLength.second, outStreamer.getString());
						outStreamer += ",\n\t\"Top hits\": [\n\t],\n\t\"Further hits\": [\n\t]\n}";
						break;

					case OutputFormat::JsonL:
						outStreamer += "{ \"Read number\": ";
						Utilities::itostr(iReadNum, outStreamer.getString());
						outStreamer += ",";
						outStreamer += " \"Specifier from input file\": \"";
						outStreamer += vReadNameAndLength.first;
						outStreamer += "\", \"Length\": ";
						Utilities::itostr(vReadNameAndLength.second, outStreamer.getString());
						outStreamer += ", \"Top hits\": [], \"Further hits\": [] }\n";
						break;

					case OutputFormat::Kraken:
						outStreamer += "U\t";
						outStreamer += vReadNameAndLength.first;
						outStreamer += "\t0\t";
						outStreamer += vReadNameAndLength.second;
						outStreamer += "\tA:00\n";
						break;
					};
					debugBarrier
				}
				else {
					debugBarrier
					sort(vTempResultVec.begin(), vTempResultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<2>(a) > get<2>(b); });

					auto maxValue = get<1>(*(max_element(vTempResultVec.begin(), vTempResultVec.end(), [](const tuple<size_t, float, double>& a, const tuple<size_t, float, double>& b) {return get<1>(a) < get<1>(b); })));
					int32_t iTopHitCounter = 1;
					for (int64_t i = 1; i < static_cast<int64_t>(vTempResultVec.size()) && i < iNumOfBeasts; ++i) {
						if ((get<1>(vTempResultVec[i])) / (maxValue) > 0.8f) {
							++iTopHitCounter;
						}
						else {
							break;
						}
					}

					// Calculate error and if Error is too high, filter it out later
					if (GlobalInputParameters.bFilter) {
						if ((bestScore - double(maxValue)) / bestScore < GlobalInputParameters.fErrorThreshold) {
							vReadIDsToBeFiltered.push_back(iReadNum);
						}
						else {
							if (GlobalInputParameters.bPostProcess) {
								if (vCoherence >= GlobalInputParameters.fCoherenceThreshold) {
									vReadIDsToBeFiltered.push_back(iReadNum);
								}
							}
						}
					}

					debugBarrier
					string sOut = "", sOut2 = "", sOut3 = "", sOut4 = "";
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

						for (int64_t j = 0; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
							Utilities::itostr(mIdxToTax[get<0>(*it)], sOut);
							sOut += ";";
							sOut2 += mOrganisms[get<0>(*it)];
							sOut2 += ";";

							dtoa_milo(get<2>(*it), sOut3);
							sOut3 += ",";
							dtoa_milo(get<1>(*it), sOut3);
							sOut3 += ";";

							dtoa_milo((bestScore - get<1>(*it)) / bestScore, sOut4);
							sOut4 += ";";

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
						if (sOut4.back() == ';') {
							sOut4.pop_back();
						}
						if (sOut2.length()) {
							outStreamer += sOut;
							outStreamer += "\t";
							outStreamer += sOut2;
							outStreamer += "\t";
							outStreamer += sOut3;
							outStreamer += "\t";
							outStreamer += sOut4;
							if (GlobalInputParameters.bPostProcess) {
								outStreamer += "\t";
								dtoa_milo(vCoherence, outStreamer.getString());
							}
							outStreamer += "\n";
						}
						debugBarrier
						break;

					case OutputFormat::Json:
						if (iReadNum == 0) {
							outStreamer += "{\n";
						}
						else {
							outStreamer += ",\n{\n";
						}

						outStreamer += "\t\"Read number\": ";
						Utilities::itostr(iReadNum, outStreamer.getString());
						outStreamer += ",\n";
						outStreamer += "\t\"Specifier from input file\": \"";
						outStreamer += vReadNameAndLength.first;
						outStreamer += "\",\n\t\"Length\": ";
						Utilities::itostr(vReadNameAndLength.second, outStreamer.getString());
						outStreamer += ",\n\t\"Top hits\": [\n";

						for (int64_t i = 0; i < iTopHitCounter; ++i, ++it) {
							if (i == 0) {
								outStreamer += "\t{\n";
							}
							else {
								outStreamer += ",\n\t{\n";
							}

							outStreamer += "\t\t\"tax ID\": \"";
							Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
							outStreamer += "\",\n";
							outStreamer += "\t\t\"Name\": \"";
							outStreamer += mOrganisms[get<0>(*it)];
							outStreamer += "\",\n";
							outStreamer += "\t\t\"k-mer Score\": ";
							dtoa_milo(get<1>(*it), outStreamer.getString());
							outStreamer += ",\n";
							outStreamer += "\t\t\"Relative Score\": ";
							dtoa_milo(get<2>(*it), outStreamer.getString());
							outStreamer += ",\n";
							outStreamer += "\t\t\"Error\": ";
							dtoa_milo((bestScore - get<1>(*it)) / bestScore, outStreamer.getString());
							if (GlobalInputParameters.bPostProcess) {
								outStreamer += ",\n";
								outStreamer += "\t\t\"Coherence\": ";
								dtoa_milo(vCoherence, outStreamer.getString());
							}
							outStreamer += "\n\t}";
						}

						outStreamer += "\n\t],\n\t\"Further hits\": [\n";

						for (int32_t j = iTopHitCounter; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
							if (j == iTopHitCounter) {
								outStreamer += "\t{\n";
							}
							else {
								outStreamer += ",\n\t{\n";
							}

							outStreamer += "\t\t\"tax ID\": \"";
							Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
							outStreamer += "\",\n";
							outStreamer += "\t\t\"Name\": \"";
							outStreamer += mOrganisms[get<0>(*it)];
							outStreamer += "\",\n";
							outStreamer += "\t\t\"k-mer Score\": ";
							dtoa_milo(get<1>(*it), outStreamer.getString());
							outStreamer += ",\n";
							outStreamer += "\t\t\"Relative Score\": ";
							dtoa_milo(get<2>(*it), outStreamer.getString());
							outStreamer += ",\n";
							outStreamer += "\t\t\"Error\": ";
							dtoa_milo((bestScore - get<1>(*it)) / bestScore, outStreamer.getString());
							if (GlobalInputParameters.bPostProcess) {
								outStreamer += ",\n";
								outStreamer += "\t\t\"Coherence\": ";
								dtoa_milo(vCoherence, outStreamer.getString());
							}
							outStreamer += "\n\t}";

							if (iValueBefore != get<1>(*it)) {
								iValueBefore = get<1>(*it);
								++j;
							}
						}
						debugBarrier
						outStreamer += "\n\t]\n}";
						break;

					case OutputFormat::JsonL:
						// json lines
						outStreamer += "{ \"Read number\": ";
						Utilities::itostr(iReadNum, outStreamer.getString());
						outStreamer += ",";
						outStreamer += " \"Specifier from input file\": \"";
						outStreamer += vReadNameAndLength.first;
						outStreamer += "\", \"Length\": ";
						Utilities::itostr(vReadNameAndLength.second, outStreamer.getString());
						outStreamer += ", \"Top hits\": [";

						for (int64_t i = 0; i < iTopHitCounter; ++i, ++it) {
							if (i == 0) {
								outStreamer += "{";
							}
							else {
								outStreamer += ",{";
							}

							outStreamer += " \"tax ID\": \"";
							Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
							outStreamer += "\",";
							outStreamer += " \"Name\": \"";
							outStreamer += mOrganisms[get<0>(*it)];
							outStreamer += "\",";
							outStreamer += " \"k-mer Score\": ";
							dtoa_milo(get<1>(*it), outStreamer.getString());
							outStreamer += ",";
							outStreamer += " \"Relative Score\": ";
							dtoa_milo(get<2>(*it), outStreamer.getString());
							outStreamer += ",";
							outStreamer += " \"Error\": ";
							dtoa_milo((bestScore - get<1>(*it)) / bestScore, outStreamer.getString());
							if (GlobalInputParameters.bPostProcess) {
								outStreamer += ",";
								outStreamer += "\"Coherence\": ";
								dtoa_milo(vCoherence, outStreamer.getString());
							}
							outStreamer += "}";
						}

						outStreamer += "], \"Further hits\": [";

						for (int32_t j = iTopHitCounter; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
							if (j == iTopHitCounter) {
								outStreamer += "{";
							}
							else {
								outStreamer += ", {";
							}

							outStreamer += " \"tax ID\": \"";
							Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
							outStreamer += "\",";
							outStreamer += " \"Name\": \"";
							outStreamer += mOrganisms[get<0>(*it)];
							outStreamer += "\",";
							outStreamer += " \"k-mer Score\": ";
							dtoa_milo(get<1>(*it), outStreamer.getString());
							outStreamer += ",";
							outStreamer += " \"Relative Score\": ";
							dtoa_milo(get<2>(*it), outStreamer.getString());
							outStreamer += ",";
							outStreamer += " \"Error\": ";
							dtoa_milo((bestScore - get<1>(*it)) / bestScore, outStreamer.getString());
							if (GlobalInputParameters.bPostProcess) {
								outStreamer += ",";
								outStreamer += "\"Coherence\": ";
								dtoa_milo(vCoherence, outStreamer.getString());
							}
							outStreamer += "}";

							if (iValueBefore != get<1>(*it)) {
								iValueBefore = get<1>(*it);
								++j;
							}
						}
						debugBarrier
						outStreamer += "] }\n";
						break;

					case OutputFormat::Kraken:
						outStreamer += "C\t";
						outStreamer += vReadNameAndLength.first;
						outStreamer += "\t";
						Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
						outStreamer += "\t";
						Utilities::itostr(vReadNameAndLength.second, outStreamer.getString());
						outStreamer += "\t";

						for (int64_t i = 0; i < iTopHitCounter; ++i, ++it) {
							Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
							outStreamer += ":";
							dtoa_milo(get<1>(*it), outStreamer.getString());
							outStreamer += " ";
						}

						for (int32_t j = iTopHitCounter; it != vTempResultVec.end() && j < iNumOfBeasts; ++it) {
							Utilities::itostr(mIdxToTax[get<0>(*it)], outStreamer.getString());
							outStreamer += ":";
							dtoa_milo(get<1>(*it), outStreamer.getString());
							outStreamer += " ";
							if (iValueBefore != get<1>(*it)) {
								iValueBefore = get<1>(*it);
								++j;
							}
						}
						debugBarrier
						outStreamer += "\n";
						break;
					};
				}
				debugBarrier
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl;
				throw;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Calculate relative score and filter out the ones crossing the threshold
		inline void scoringFuncWithFilterNoOutput(const Utilities::Non_contiguousArray& vReadIDtoGenID, const uint64_t& iStart, const uint64_t& iEnd, const uint64_t& iRealReadIDStart, list<pair<string, readIDType>>& vReadNameAndLength, const uint64_t& iNumberOfSpecies, vector<uint64_t>& vReadIDsToBeFiltered, const vector<float>& vCoherence) {
			debugBarrier
			for (uint64_t readIdx = iStart; readIdx < iEnd; ++readIdx) {
				auto currentReadLengthAndName = vReadNameAndLength.front();

				float bestScore = calculateBestScore(currentReadLengthAndName.second);

				for (size_t iSpecIdx = 1; iSpecIdx < iNumberOfSpecies; ++iSpecIdx) {
					if (vReadIDtoGenID[readIdx][iSpecIdx] > 0.f) {
						// Calculate error and if Error is too high, filter it out
						if ((bestScore - double(vReadIDtoGenID[readIdx][iSpecIdx])) / bestScore < GlobalInputParameters.fErrorThreshold) {
							vReadIDsToBeFiltered.push_back(iRealReadIDStart + readIdx);
							break; // it suffices, that one species was hit "good enough"
						}
						else {
							if (GlobalInputParameters.bPostProcess) {
								if (vCoherence[readIdx] >= GlobalInputParameters.fCoherenceThreshold) {
									vReadIDsToBeFiltered.push_back(iRealReadIDStart + readIdx);
									break;
								}
							}
						}
					}
				}
			}
			debugBarrier
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Calculate error and filter out the ones below the threshold
		inline void scoringFuncWithFilterNoOutputOne(vector<tuple<readIDType, float, double>>&& vTempResultVec, const uint64_t& iReadNum, const pair<string, uint32_t>& vReadNameAndLength, vector<uint64_t>& vReadIDsToBeFiltered, const float& vCoherence) {
			debugBarrier
			float bestScore = calculateBestScore(vReadNameAndLength.second);
			
			for (auto it = vTempResultVec.begin(); it != vTempResultVec.end();) {
				// Calculate error and if Error is too high, filter it out
				if ((bestScore - double(get<1>(*it))) / bestScore < GlobalInputParameters.fErrorThreshold) {
					vReadIDsToBeFiltered.push_back(iReadNum);
					break;
				}
				else {
					if (GlobalInputParameters.bPostProcess) {
						if (vCoherence >= GlobalInputParameters.fCoherenceThreshold) {
							vReadIDsToBeFiltered.push_back(iReadNum);
						}
					}
				}
			}
			debugBarrier
		}

		////////////////////////////////////////////////////////////////////////////////////////////////
		// save results from read to tax ID
		inline void saveResults(const bool transferBetweenRuns_finished, const bool transferBetweenRuns_addTail, vector<tuple<readIDType, float, double>>& vSavedScores, const Utilities::Non_contiguousArray& vReadIDtoTaxID, list<pair<string, uint32_t>> vReadNameAndLength, const uint64_t& iNumberOfSpecies, const vector<uint32_t>& mIdxToTax, const vector<string>& mOrganisms, const vector<vector<uint64_t>>& mFrequencies, uint64_t& iNumOfReadsSum, uint64_t& iNumOfReadsOld, uint64_t iNumOfReads, const float& fThreshold, vector<uint64_t>& vReadIDsToBeFiltered, const vector<float> vCoherence, ofstream& fOut) {

			auto getVecOfScored = [&iNumberOfSpecies](const Utilities::Non_contiguousArray& vReadIDtoGenID, const uint64_t& iIdx) {
				try {
					vector<tuple<uint32_t, float, double>> output;
					for (uint32_t i = 1; i < iNumberOfSpecies; ++i) {
						if (vReadIDtoGenID[iIdx][i] > 0.f) {
							output.push_back(make_tuple(i, vReadIDtoGenID[iIdx][i], 0.0));
						}
					}
					return output;
				}
				catch (...) {
					cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
				}
			};
			debugBarrier

			uint64_t iReadIDStart = 0;
			// check if there is still some unfinished read which is now complete
			if (vSavedScores.size() && transferBetweenRuns_finished) {
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
				debugBarrier

				if (fOut.good()) {
					scoringFunc(move(vSavedScores), iReadIDStart + iNumOfReadsSum, vReadNameAndLength.front(), mFrequencies, mIdxToTax, mOrganisms, fThreshold, vReadIDsToBeFiltered, vCoherence[iReadIDStart], move(fOut));
					iReadIDStart++;
				}
				else {
					if (GlobalInputParameters.bFilter) {
						scoringFuncWithFilterNoOutputOne(move(vSavedScores), iReadIDStart + iNumOfReadsSum, vReadNameAndLength.front(), vReadIDsToBeFiltered, vCoherence[iReadIDStart]);
						iReadIDStart++;
					}
					else {
						cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw runtime_error("Cannot write to identification output!");
					}
				}
				
				vSavedScores.clear();

				vReadNameAndLength.pop_front();
			}
			debugBarrier

			if (transferBetweenRuns_addTail) {
				// last read is not yet finished
				// save the score of the not yet finished
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

				debugBarrier
				// save the finished ones
				if (fOut.good()) {
					scoringFunc(vReadIDtoTaxID, iReadIDStart, iNumOfReads - 1, iNumOfReadsSum, vReadNameAndLength, iNumberOfSpecies, mFrequencies, mIdxToTax, mOrganisms, fThreshold, vReadIDsToBeFiltered, vCoherence, move(fOut));
				}
				else {
					if (GlobalInputParameters.bFilter) {
						scoringFuncWithFilterNoOutput(vReadIDtoTaxID, iReadIDStart, iNumOfReads - 1, iNumOfReadsSum, vReadNameAndLength, iNumberOfSpecies, vReadIDsToBeFiltered, vCoherence);
					}
					else {
						cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw runtime_error("Cannot write to identification output!");
					}
				}

				iNumOfReadsSum += iNumOfReads - 1;
				iNumOfReadsOld = (iNumOfReads - 1 < iNumOfReadsOld) ? iNumOfReadsOld : iNumOfReads - 1;
			}
			else {
				debugBarrier
				// reads finished
				if (fOut.good()) {
					scoringFunc(vReadIDtoTaxID, iReadIDStart, iNumOfReads, iNumOfReadsSum, vReadNameAndLength, iNumberOfSpecies, mFrequencies, mIdxToTax, mOrganisms, fThreshold, vReadIDsToBeFiltered, vCoherence, move(fOut));
				}
				else {
					if (GlobalInputParameters.bFilter) {
						scoringFuncWithFilterNoOutput(vReadIDtoTaxID, iReadIDStart, iNumOfReads, iNumOfReadsSum, vReadNameAndLength, iNumberOfSpecies, vReadIDsToBeFiltered, vCoherence);
					}
					else {
						cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw runtime_error("Cannot write to identification output!");
					}
				}

				iNumOfReadsSum += iNumOfReads;
				iNumOfReadsOld = iNumOfReads;
			}
			debugBarrier
		}

		template <typename inFileType, typename outFileType>
		inline void filter(const vector<uint64_t>& vReadIDsToBeFiltered, unique_ptr<inFileType>& inFile1, unique_ptr<inFileType>& inFile2, const bool& bIsFasta, const bool& bIsPairedEnd) {
			unique_ptr<outFileType> fOutClean_1;
			unique_ptr<outFileType> fOutClean_2;
			unique_ptr<outFileType> fOutContaminated_1;
			unique_ptr<outFileType> fOutContaminated_2;
			string fastqOrA = (bIsFasta) ? "a" : "q";
			string gzipOrNot = (GlobalInputParameters.bGzipOut) ? ".gz" : "";

			// open files for writing
			if (bIsPairedEnd) {
				if (GlobalInputParameters.sFilteredCleanOut != "_") {
					fOutClean_1.reset(new outFileType((GlobalInputParameters.sFilteredCleanOut + "_1.fast" + fastqOrA + gzipOrNot).c_str()));
					fOutClean_2.reset(new outFileType((GlobalInputParameters.sFilteredCleanOut + "_2.fast" + fastqOrA + gzipOrNot).c_str()));
					if (fOutClean_1->fail() || fOutClean_2->fail()) {
						throw runtime_error("Filtered output files could not be opened for writing, did you use the correct path?");
					}
				}
				if (GlobalInputParameters.sFilteredContaminantsOut != "_") {
					fOutContaminated_1.reset(new outFileType((GlobalInputParameters.sFilteredContaminantsOut + "_1.fast" + fastqOrA + gzipOrNot).c_str()));
					fOutContaminated_2.reset(new outFileType((GlobalInputParameters.sFilteredContaminantsOut + "_2.fast" + fastqOrA + gzipOrNot).c_str()));
					if (fOutContaminated_1->fail() || fOutContaminated_2->fail()) {
						throw runtime_error("Contaminants output files could not be opened for writing, did you use the correct path?");
					}
				}
			}
			else {
				if (GlobalInputParameters.sFilteredCleanOut != "_") {
					fOutClean_1.reset(new outFileType((GlobalInputParameters.sFilteredCleanOut + ".fast" + fastqOrA + gzipOrNot).c_str()));
					if (fOutClean_1->fail()) {
						throw runtime_error("Filtered output files could not be opened for writing, did you use the correct path?");
					}
				}
				if (GlobalInputParameters.sFilteredContaminantsOut != "_") {
					fOutContaminated_1.reset(new outFileType((GlobalInputParameters.sFilteredContaminantsOut + ".fast" + fastqOrA + gzipOrNot).c_str()));
					if (fOutContaminated_1->fail()) {
						throw runtime_error("Contaminants output files could not be opened for writing, did you use the correct path?");
					}
				}
			}

			// if no contaminants were found, copy the content
			if (vReadIDsToBeFiltered.empty() && GlobalInputParameters.sFilteredCleanOut != "_") {
				*fOutClean_1 << inFile1->rdbuf();
				if (bIsPairedEnd) {
					*fOutClean_2 << inFile2->rdbuf();
				}
				return;
			}

			// go through input
			uint64_t iReadID = 0, iFilterVecIdx = 0;
			bool bCleanOrContaminated = true;
			if (bIsFasta) {
				string sDummyString_1 = "", sDummyString_2 = "";
				while (getline(*inFile1, sDummyString_1)) {
					if (bIsPairedEnd) {
						getline(*inFile2, sDummyString_2);
					}
					if (sDummyString_1 == "") {
						continue;
					}

					if (sDummyString_1.front() == '>') {
						if (iReadID == vReadIDsToBeFiltered[iFilterVecIdx]) {
							// contaminated
							bCleanOrContaminated = false;
							iReadID++;
							iFilterVecIdx++;
						}
						else {
							// Clean
							bCleanOrContaminated = true;
							iReadID++;
						}
					}

					// write to respective output
					if (bCleanOrContaminated) {
						if (GlobalInputParameters.sFilteredCleanOut != "_") {
							(*fOutClean_1) << sDummyString_1 << endl;
							if (bIsPairedEnd) {
								(*fOutClean_2) << sDummyString_2 << endl;
							}
						}
					}
					else {
						if (GlobalInputParameters.sFilteredContaminantsOut != "_") {
							(*fOutContaminated_1) << sDummyString_1 << endl;
							if (bIsPairedEnd) {
								(*fOutContaminated_2) << sDummyString_2 << endl;
							}
						}
					}
				}
			}
			else { // fastq
				while (inFile1->good()) {
					array<string, 4> fastqContent, fastqContent2;
					getline(*inFile1, fastqContent[0]); // specifier
					getline(*inFile1, fastqContent[1]); // dna/protein sequence
					getline(*inFile1, fastqContent[2]); // +
					getline(*inFile1, fastqContent[3]); // quality

					if (bIsPairedEnd) {
						getline(*inFile2, fastqContent2[0]);
						getline(*inFile2, fastqContent2[1]);
						getline(*inFile2, fastqContent2[2]);
						getline(*inFile2, fastqContent2[3]);
					}

					if (fastqContent[0] == "") {
						continue;
					}

					if (iReadID == vReadIDsToBeFiltered[iFilterVecIdx]) {
						// contaminated
						bCleanOrContaminated = false;
						iReadID++;
						iFilterVecIdx++;
					}
					else {
						// Clean
						bCleanOrContaminated = true;
						iReadID++;
					}

					// write to respective output
					if (bCleanOrContaminated) {
						if (GlobalInputParameters.sFilteredCleanOut != "_") {
							for (int i = 0; i < 4; ++i) {
								(*fOutClean_1) << fastqContent[i] << endl;
							}

							if (bIsPairedEnd) {
								for (int i = 0; i < 4; ++i) {
									(*fOutClean_2) << fastqContent2[i] << endl;
								}
							}
						}
					}
					else {
						if (GlobalInputParameters.sFilteredContaminantsOut != "_") {
							for (int i = 0; i < 4; ++i) {
								(*fOutContaminated_1) << fastqContent[i] << endl;
							}

							if (bIsPairedEnd) {
								for (int i = 0; i < 4; ++i) {
									(*fOutContaminated_2) << fastqContent2[i] << endl;
								}
							}
						}
					}
				}
			}
		}


		inline void postProcess(InputType<intType>& vIn, vector<float>& vCoherenceScores, const int32_t& iLocalNumOfThreads) {
			// sort
			auto sortfun = [](const typename InputType<intType>::ppTuple& p1, const typename InputType<intType>::ppTuple& p2) {
				return (get<3>(p1) < get<3>(p2))
					|| (get<3>(p1) == get<3>(p2)
						&& ((get<5>(p1) >> 5) < (get<5>(p2) >> 5)))
					|| (get<3>(p1) == get<3>(p2)
						&& ((get<5>(p1) >> 5) == (get<5>(p2) >> 5))
						&& (get<4>(p1) < get<4>(p2)));};

			if (iLocalNumOfThreads == 1) {
				std::sort(vIn.getPPVec()->begin(), vIn.getPPVec()->end(), sortfun);
			}
# if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
			Utilities::parallelQuicksort(vIn.getPPVec()->begin(), vIn.getPPVec()->end(), sortfun, iLocalNumOfThreads);
#else
#if __has_include(<execution>) && !defined(__llvm__)
			std::sort(std::execution::par_unseq, vIn.getPPVec()->begin(), vIn.getPPVec()->end(), sortfun);
#else
			Utilities::parallelQuicksort(vIn.getPPVec()->begin(), vIn.getPPVec()->end(), sortfun, iLocalNumOfThreads);
#endif
#endif
			// Search for first match
			readIDType iReadID = 0;
			uint64_t iIdxOfIn = 0;
			uint32_t iLastMatchPosAndLength = 0;
			uint32_t iCurrentOverlap = 0;
			uint32_t iCountOfMax = 0;
			while (iIdxOfIn < vIn.size()) {
				const auto& iMatchLength = vIn.getLength(iIdxOfIn);
				if (iMatchLength != 0) {
					iReadID = vIn.getReadID(iIdxOfIn);
					iLastMatchPosAndLength = vIn.getPosition(iIdxOfIn) + iMatchLength;
					++iIdxOfIn;
					break;
				}
				else {
					++iIdxOfIn;
				}
			}

			auto setClusterLength = [](float& iCurrentClusterLength, const float& iNewClusterLength) {
				iCurrentClusterLength = max(iCurrentClusterLength, iNewClusterLength);
			};
			auto determineMaxAndCount = [&iCountOfMax,&iCurrentOverlap](const uint32_t& nextOverlap) {
				if (nextOverlap > iCurrentOverlap) {
					iCurrentOverlap = nextOverlap;
					iCountOfMax = 1;
				} else {
					if (nextOverlap == iCurrentOverlap) {
						iCountOfMax++;
					}
				}
			};

			// go through by read ID
			for (; iReadID < vCoherenceScores.size() && iIdxOfIn < vIn.size(); ++iReadID) {
				// and by frame (forward or backward)
				for (int8_t iForwardOrBackward = 0; iForwardOrBackward < 1 + Base::_bSixFrames;) {
					// start cluster

					const auto& iMatchLength = vIn.getLength(iIdxOfIn);
					if (iMatchLength != 0) {
						// if position + length > next pos -> elongate cluster
						if (vIn.getPosition(iIdxOfIn) <= iLastMatchPosAndLength) {
							if (vIn.getPosition(iIdxOfIn) + iMatchLength < iLastMatchPosAndLength) {
								determineMaxAndCount(static_cast<uint32_t>(iMatchLength));
							}
							else {
								auto overlap = static_cast<int32_t>(iLastMatchPosAndLength) - static_cast<int32_t>(vIn.getPosition(iIdxOfIn));
								determineMaxAndCount(overlap);
							}
							
						}
						else {
							// else add last cluster to score for that read ID
							setClusterLength(vCoherenceScores[iReadID], static_cast<float>(iCurrentOverlap) + 1.0f - 1.0f/iCountOfMax);
							iCurrentOverlap = 0;
						}
						iLastMatchPosAndLength = vIn.getPosition(iIdxOfIn) + iMatchLength;
					}

					++iIdxOfIn;
					// Check if vIn was completey processed
					if (iIdxOfIn == vIn.size()) {
						setClusterLength(vCoherenceScores[iReadID], static_cast<float>(iCurrentOverlap) + 1.0f - 1.0f / iCountOfMax);
						break;
					}

					// Check if read is finished
					if (vIn.getReadID(iIdxOfIn) != iReadID) {
						// save score and reset everything
						setClusterLength(vCoherenceScores[iReadID], static_cast<float>(iCurrentOverlap) + 1.0f - 1.0f / iCountOfMax);
						iLastMatchPosAndLength = numeric_limits<uint32_t>::max();
						iCurrentOverlap = 0;
						iCountOfMax = 0;
						break; // get new readID
					}

					// Check if RC is now active
					if (vIn.getFrame(iIdxOfIn) != iForwardOrBackward) {
						// save score and reset everything
						setClusterLength(vCoherenceScores[iReadID], static_cast<float>(iCurrentOverlap) + 1.0f - 1.0f / iCountOfMax);
						iCurrentOverlap = 0;
						iCountOfMax = 0;
						++iForwardOrBackward;

						while (iIdxOfIn < vIn.size()) {
							const auto& iMatchLength = vIn.getLength(iIdxOfIn);
							if (iMatchLength != 0) {
								iLastMatchPosAndLength = vIn.getPosition(iIdxOfIn) + iMatchLength;
								++iIdxOfIn;
								break;
							}
							else {
								++iIdxOfIn;
							}
						}
					}
				}
			}
		}

	public:

		/////////////////////////////////////////////////////////////////////////////////
		void CompareWithLib_partialSort(const string& fInFile, const string& fOutFile, const string& fTableFile, const int64_t& iMemory, bool bRAM, const bool& bUnique, const float& fThreshold, int32_t iLocalNumOfThreads = 0) {
			try {
				if (iLocalNumOfThreads == 0) {
					iLocalNumOfThreads = Base::_iNumOfThreads;
				}

//#define TIME
#ifdef TIME
				auto startTIME = std::chrono::high_resolution_clock::now();
#endif

				bool bShrinkAvailMemoryAfterFirstStepOnce = true;

				debugBarrier

				// Assert, that if the index was shrunken via trieHalf MinK can only be larger than 6
				if (index.bPartitioned && Base::_iMinK <= 6 && Base::_iMaxK <= 12) {
					cerr << "WARNING: k can only be larger than 6 and smaller than 12 if your index was shrunken via strategy 2. Setting it right..." << endl;
					if (Base::_iMinK <= 6) {
						Base::_iMinK = 7;
					}
					if (Base::_iMaxK > 12) {
						Base::_iMaxK = 12;
					}
					Base::_iNumOfK = Base::_iMaxK - Base::_iMinK + 1;
				}
				if (bUnfunny) {
					Base::_iMinK = 6;
					Base::_iMaxK = 6;
					Base::_iNumOfK = Base::_iMaxK - Base::_iMinK + 1;
				}

				// Create threadpool(s), in stxxl mode we can only have synced parallelism (as in no two threads should not access the same vector instance)
				vector<WorkerThread> workerThreadPool(iLocalNumOfThreads);
				for (int32_t i = 0; i < iLocalNumOfThreads; ++i) {
					workerThreadPool[i].setID(i);
				}


				// load index and get trie
				unique_ptr<ReadIndex> localIndex;
				if (index.bIdentifyMultiple) {
					localIndex.reset(new ReadIndex(index));
					localIndex->set(false, bUnfunny, iLocalNumOfThreads);
				}

				debugBarrier


#ifdef TIME
				auto endTIME = std::chrono::high_resolution_clock::now();
				cout << "Index load_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				// are we using a trie or a different method?
				if (Base::_bVerbose) {
					index.trieForVector->GetIfVecIsUsed();
				}

#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "Trie load_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif
				debugBarrier

				// If the user wishes to load the index into RAM, check if there is enough space.
				// First calculate the average memory usage: 
				//	Size of trie + size of counts_all/_unique + index size + memory for tax ids + chunksize from input + overhead buffer for various allocations
				const uint64_t& iMult = uint64_t(iLocalNumOfThreads) * Base::_iNumOfK * uint64_t(index.iNumOfSpecies);
				int64_t iMemoryUsageOnAverage = ((index.bIdentifyMultiple) ? (GIGABYTEASBYTES / iLocalNumOfThreads) : GIGABYTEASBYTES) + iMult * (sizeof(uint64_t) + sizeof(double) + sizeof(uint64_t)) + uint64_t(iLocalNumOfThreads) * Utilities::sBitArray(index.iNumOfSpecies).sizeInBytes() + 14399756 + 4 * index.iNumOfSpecies;

				// Set memory boundaries
				int64_t iSoftMaxMemoryAvailable = 0;
				if (iMemory > iMemoryUsageOnAverage) {
					iSoftMaxMemoryAvailable = iMemory - iMemoryUsageOnAverage;
				}
				else {
					cerr << "WARNING: Not enough memory given, trying to use the amount given but may lead to bad_alloc errors..." << endl;
					//iSoftMaxMemoryUsage = 1024ull * 1024ull * 1ull / (sizeof(tuple<uint64_t, uint64_t>)*iLocalNumOfThreads);
					iSoftMaxMemoryAvailable = iMemory; // 1024ull * 1024ull * 1024ull
				}

				bRAM = index.bRAM; // if not enough primary memory is available, the variable will be set to false inside "index"
				if (index.bIdentifyMultiple && !bRAM) {
					localIndex->loadIndex(Base::_iMaxK, Base::_iMinK, iSoftMaxMemoryAvailable);
					if (iSoftMaxMemoryAvailable < 0)  {
						cerr << "WARNING: Not enough memory given, trying with 1GB but may lead to bad_alloc errors..." << endl;
						iSoftMaxMemoryAvailable = static_cast<int64_t>(GIGABYTEASBYTES); // 1024ull * 1024ull * 1024ull
					}
				}
				

				// This holds the hits for each organism
				// Initialize this here to calculate the average memory usage
				
				unique_ptr<double[]> vCount_all(new double[iMult]);
				unique_ptr<uint64_t[]> vCount_total;
				unique_ptr<uint64_t[]> vCount_unique(new uint64_t[iMult]);
				if (this->bCoverage) {
					vCount_total.reset(new uint64_t[iMult]);
					fill_n(vCount_total.get(), iMult, 0);
				}

				fill_n(vCount_all.get(), iMult, 0.);
				fill_n(vCount_unique.get(), iMult, 0);

				debugBarrier

#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "Load Index from Memory_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
				startTIME = std::chrono::high_resolution_clock::now();
#endif

				
				uint64_t iTimeFastq = 0, iTimeCompare = 0, iNumOfReads = 0, iNumOfReadsOld = 0, iNumOfReadsSum = 0;
				InputType<intType> vInputVec(GlobalInputParameters.bPostProcess);

				// Memory management
				m_exceptionLock.lock();
				while (true) {
					try {
						vInputVec.reserve((iSoftMaxMemoryAvailable) / vInputVec.sizeOf());
						vInputVec.shrink_to_fit(); // just to see, if this is possible. If it is, we can reset it to zero here as we will work with reserve later anyway
						//vInputVec.resize(iSoftMaxMemoryUsage / sizeof(tuple<uint64_t, uint64_t, uint32_t, uint32_t>));
						break;
					}
					catch (const bad_alloc&) {
						cout << "OUT: Not enough memory available, reducing usage by 1GB..." << endl;
						//vInputVec.clear();
						vInputVec.shrink_to_fit();
						if ((iSoftMaxMemoryAvailable + 524288000) > static_cast<int64_t>(GIGABYTEASBYTES)) { // 1073741824 == 1024^3 == 1GB
							iSoftMaxMemoryAvailable -= GIGABYTEASBYTES;
						}
						else {
							cerr << "ERROR: Your system has not enough contiguous memory available. Please restart the system!" << endl;
							m_exceptionLock.unlock();
							throw;
						}
					}
				}
				m_exceptionLock.unlock();
				debugBarrier
				const bool& bReadIDsAreInteresting = (fOutFile != "" || GlobalInputParameters.bFilter);

				//Utilities::Vector2D<float> vReadIDtoTaxID; // Array of #species times #reads with scores as values
				Utilities::Non_contiguousArray vReadIDtoTaxID;

				list<pair<string, uint32_t>> vReadNameAndLength;

				vector<uint64_t> vReadIDsToBeFiltered; // in case filtering is enabled

				uint64_t iNumberOfkMersInInput = 0;

				pair<vector<pair<string, int8_t>>, size_t> pathAndSize;
				
				const bool& bIsPairedEnd = fInFile.find_first_of(char(0xB)) != string::npos;
				vector<string> vPairedEndFiles;
				if (bIsPairedEnd) {
					vPairedEndFiles = Utilities::split(fInFile, char(0xB));
					pathAndSize = Utilities::gatherFilesFromPath(vPairedEndFiles[0]);
				}
				else {
					pathAndSize = Utilities::gatherFilesFromPath(fInFile);
				}
				const auto& vInputFiles = pathAndSize.first;
				
				
				size_t overallFileSize = pathAndSize.second, allFilesProgress = 0, charsReadOverall = 0;

				debugBarrier
				// allow multiple input files
				for (const auto& inFile : vInputFiles) {

					if (Base::_bVerbose) {
						if (bIsPairedEnd) {
							cout << "OUT: Paired files: " << vPairedEndFiles[0] << " and " << vPairedEndFiles[1] << endl;
						}
						else {
							cout << "OUT: Current file: " << inFile.first << endl;
						}
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
					Utilities::FileReader<ifstream> fileReaderObject;
					Utilities::FileReader<igzstream> fileReaderObject_gz;
					// paired
					unique_ptr<ifstream> fast_q_a_File2;
					unique_ptr<igzstream> fast_q_a_File_gz2;
					Utilities::FileReader<ifstream> fileReaderObject2;
					Utilities::FileReader<igzstream> fileReaderObject_gz2;

					uint64_t iFileLength = 0;
					vector<char> inFilebuffer(2097152);
					vector<char> inFilebuffer2(2097152);

					if (isGzipped) {
						if (Base::_bVerbose) {
							cout << "OUT: File is gzipped, no progress output can be shown." << endl;
						}

						if (bIsPairedEnd) {
							fast_q_a_File_gz2.reset(new igzstream());
							fast_q_a_File_gz2->rdbuf()->pubsetbuf(&inFilebuffer2[0], 2097152);
							fast_q_a_File_gz2->open(vPairedEndFiles[1].c_str());
							fileReaderObject_gz2.setFile(fast_q_a_File_gz2.get());
							fast_q_a_File_gz2->exceptions(std::ios_base::badbit);
						}

						fast_q_a_File_gz.reset(new igzstream());
						fast_q_a_File_gz->rdbuf()->pubsetbuf(&inFilebuffer[0], 2097152);
						fast_q_a_File_gz->open(inFile.first.c_str());
						bIsGood = fast_q_a_File_gz->good();
						
						// Determine if protein or DNA sequence
						if (!Base::_bProtein) {
							const string& sFirstSequence = Utilities::getFirstSequenceOfFile(*fast_q_a_File_gz);
							Base::detectAlphabet(sFirstSequence);
							// gzstream does not allow seek, therefore the file must be resetted the hard way
							fast_q_a_File_gz.reset(new igzstream());
							fast_q_a_File_gz->rdbuf()->pubsetbuf(&inFilebuffer[0], 2097152);
							fast_q_a_File_gz->open(inFile.first.c_str());
						}

						// Determine if fasta or fastq
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
						fileReaderObject_gz.setFile(fast_q_a_File_gz.get());
						fast_q_a_File_gz->exceptions(std::ios_base::badbit);
					}
					else {

						if (bIsPairedEnd) {
							fast_q_a_File2.reset(new ifstream());
							fast_q_a_File2->rdbuf()->pubsetbuf(&inFilebuffer2[0], 2097152);
							fast_q_a_File2->open(vPairedEndFiles[1].c_str());
							fileReaderObject2.setFile(fast_q_a_File2.get());
						}

						fast_q_a_File.reset(new ifstream());
						fast_q_a_File->rdbuf()->pubsetbuf(&inFilebuffer[0], 2097152);
						fast_q_a_File->open(inFile.first);
						bIsGood = fast_q_a_File->good();

						// Determine if protein or DNA sequence
						// Determine if protein or DNA sequence
						if (!Base::_bProtein) {
							const string& sFirstSequence = Utilities::getFirstSequenceOfFile(*fast_q_a_File);
							Base::detectAlphabet(sFirstSequence);
						}

						// Determine if fasta or fastq
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

						fileReaderObject.setFile(fast_q_a_File.get());

						iFileLength = Utilities::getSizeOfFile(inFile.first);

					}
					debugBarrier
					// output file
					ofstream fOut;
					// a larger buffer works better for SSDs or HPCCs
					vector<char> outFilebuffer(104857600); // 104857600 = 100MB, 2097152 = 2MB
					fOut.rdbuf()->pubsetbuf(&outFilebuffer[0], 104857600);
					//fOut.exceptions(std::ifstream::failbit | std::ifstream::badbit);
					if (bReadIDsAreInteresting) {
						if (fOutFile != "") {
							fOut.open((vInputFiles.size() > 1) ? fOutFile + fileName + outputFormatFileEnding() : fOutFile); // in case of multiple input files, specify only beginning of the output and the rest will be appended
							if (fOut) {
								if (format == OutputFormat::tsv) {
									if (GlobalInputParameters.bPostProcess) {
										fOut << "#Read number\tSpecifier from input file\tMatched taxa\tNames\tScores{relative,k-mer}\tError\tCoherence" << "\n";
									}
									else {
										fOut << "#Read number\tSpecifier from input file\tMatched taxa\tNames\tScores{relative,k-mer}\tError" << "\n";
									}
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
					}

					// Thread that handles output while other stuff happens
					unique_ptr<thread> tOutputThread;

					// Test if profile file can be written (it would suck to see your calculation gone to waste because it is only tested at the end now wouldn't it?!)
					if (fTableFile != "") {
						ofstream tableFileStream((vInputFiles.size() > 1) ? fTableFile + fileName + ".csv" : fTableFile);
						if (!tableFileStream) {
							throw runtime_error("Profile file couldn't be opened for writing!");
						}
					}

#ifdef TIME
					endTIME = std::chrono::high_resolution_clock::now();
					cout << "Prepare input_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
					startTIME = std::chrono::high_resolution_clock::now();
#endif
					debugBarrier
					unique_ptr<typename Base::strTransfer> transferBetweenRuns(new typename Base::strTransfer);
					transferBetweenRuns->iCurrentOverallPercentage = allFilesProgress;
					transferBetweenRuns->iNumOfAllCharsRead = charsReadOverall;
					uint64_t iLastProgress = 0;
					auto timeWholePartStart = std::chrono::high_resolution_clock::now();
					vector<tuple<readIDType, float, double>> vSavedScores;
					//readIDType iReadIDofSavedScores = 0;

					// read input
					while (bIsGood) {
						ios::sync_with_stdio(false);
						std::cin.tie(nullptr);
						auto timeInputStart = std::chrono::high_resolution_clock::now();
						

						//iSoftMaxMemoryUsage = 19199040;

						if (isGzipped) {
							iNumberOfkMersInInput += Base::readFastqa_partialSort(fileReaderObject_gz, fileReaderObject_gz2, vInputVec, vReadNameAndLength, iSoftMaxMemoryAvailable, index.iNumOfSpecies, iFileLength, iFileLength, bReadIDsAreInteresting, bIsFasta, transferBetweenRuns, workerThreadPool);

						}
						else {
							iNumberOfkMersInInput += Base::readFastqa_partialSort(fileReaderObject, fileReaderObject2, vInputVec, vReadNameAndLength, iSoftMaxMemoryAvailable, index.iNumOfSpecies, iFileLength, overallFileSize, bReadIDsAreInteresting, bIsFasta, transferBetweenRuns, workerThreadPool);
						}
	
						// reduce available memory
						// context: should the input vector increase in size after the first call, bad_alloc will probably be thrown because there will likely be no larger contiguous chunk of memory available
						// this is mitigated by substracting 0.1% of the available memory after every time this is called
						if (iSoftMaxMemoryAvailable - static_cast<int64_t>(iSoftMaxMemoryAvailable * 0.001) > 0 && bShrinkAvailMemoryAfterFirstStepOnce) {
							iSoftMaxMemoryAvailable -= static_cast<int64_t>(iSoftMaxMemoryAvailable * 0.001);
							bShrinkAvailMemoryAfterFirstStepOnce = false;
						}
						debugBarrier
						//iNumOfReads = (transferBetweenRuns->iNumOfNewReads > 0) ? transferBetweenRuns->iNumOfNewReads : 1; // vReadIDs.size();
						iNumOfReads = transferBetweenRuns->iNumOfNewReads;

						/*for (int i = 0; i < 5; i++) {
							cout << Base::kMerToAminoacid(get<1>(vInputVec[i]), 25) << endl;
						}*/


#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "Input_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						// sort prefixes for each range in parallel
						if (GlobalInputParameters.bPostProcess) {
							sortInputAndCheckInvalidkMers_pp(vInputVec, index.bPartitioned, *(index.trieForVector), iLocalNumOfThreads);
						}
						else {
							sortInputAndCheckInvalidkMers_sta(vInputVec, index.bPartitioned, *(index.trieForVector), iLocalNumOfThreads);
						}

						/*for (int i = 0; i < 5; i++) {
							cout << Base::kMerToAminoacid(get<1>(vInputVec[i]), 25) << endl;
						}*/
						debugBarrier

#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "PAR Sort and Elimination_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						if (bUnique) {
							if (GlobalInputParameters.bPostProcess) {
								auto ppFun = [&](typename InputType<intType>::ppTuple& a, typename InputType<intType>::ppTuple& b) { return (get<1>(a) == get<1>(b)) && (get<3>(a) == get<3>(b)); };
								auto newEnd = unique(vInputVec.getPPVec()->begin(), vInputVec.getPPVec()->end(), ppFun);
								vInputVec.resize(newEnd - vInputVec.getPPVec()->begin());
							}
							else {
								auto staFun = [&](typename InputType<intType>::staTuple& a, typename InputType<intType>::staTuple& b) { return (get<1>(a) == get<1>(b)) && (get<3>(a) == get<3>(b)); };
								auto newEnd = unique(vInputVec.getStaVec()->begin(), vInputVec.getStaVec()->end(), staFun);
								vInputVec.resize(newEnd - vInputVec.getStaVec()->begin());
							}
						}
						debugBarrier
						auto timeInputEnd = std::chrono::high_resolution_clock::now();
						iTimeFastq += chrono::duration_cast<std::chrono::nanoseconds>(timeInputEnd - timeInputStart).count();


						if (bReadIDsAreInteresting) {
							// join previous thread should it still be running
							if (tOutputThread) {
								tOutputThread->join();
							}
							debugBarrier
								// This holds the mapping read ID -> Genus ID

							vReadIDtoTaxID.generate(iNumOfReads, index.iNumOfSpecies);

						}
						else {
							iNumOfReadsSum += transferBetweenRuns->iNumOfNewReads - transferBetweenRuns->addTail;
						}
						debugBarrier

						
						auto timeWholePartEnd = std::chrono::high_resolution_clock::now();
						if (Base::_bVerbose && !isGzipped) {
							if (iLastProgress != 0) {
								cout << "OUT: Estimated remaining time needed for this file: " << static_cast<double>(iFileLength - transferBetweenRuns->iNumOfCharsRead) / (transferBetweenRuns->iNumOfCharsRead - iLastProgress) * static_cast<double>(chrono::duration_cast<std::chrono::seconds>(timeWholePartEnd - timeWholePartStart).count()) << "s" << endl;
								timeWholePartStart = std::chrono::high_resolution_clock::now();
							}
						}

						iLastProgress = transferBetweenRuns->iNumOfCharsRead;


						debugBarrier
						function<void(const int32_t&, const uint64_t&, const uint64_t&)> foo;

						// now compare with index
						auto timeCompareStart = std::chrono::high_resolution_clock::now();


						if (bRAM) {
							if (index.bPartitioned) {
								foo = bind(&Compare::compareWithDatabase< vector<packedPair>*>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, index.getLibRAMHalf(), ref(vCount_all), ref(vCount_unique), ref(vCount_total), ref(vReadIDtoTaxID), ref(index.iNumOfSpecies), ref(*index.mTaxToIdx));//, ref(transferBetweenRuns->mReadIDToArrayIdx));
							}
							else {
								if (bUnfunny) {
									//foo = bind(&Compare::compareWithDatabase_sloppy< vector<uint16_t>*>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, &vLib_RAM_taxaOnly, ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(transferBetweenRuns->mReadIDToArrayIdx));
								}
								else {
									foo = bind(&Compare::compareWithDatabase< vector<elemType>*>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, index.getLibRAMFull(), ref(vCount_all), ref(vCount_unique), ref(vCount_total), ref(vReadIDtoTaxID), ref(index.iNumOfSpecies), ref(*index.mTaxToIdx));//, ref(transferBetweenRuns->mReadIDToArrayIdx));
								}
							}
						}
						else {
							if (index.bPartitioned) {
								unique_ptr<unique_ptr<const index_t_p>[]>* indexPtr = nullptr;
								if (index.bIdentifyMultiple) {
									indexPtr = localIndex->getIndexTP();
								}
								else {
									indexPtr = index.getIndexTP();
								}
								foo = bind(&Compare::compareWithDatabase< unique_ptr<unique_ptr<const index_t_p>[]>>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(*indexPtr), ref(vCount_all), ref(vCount_unique), ref(vCount_total), ref(vReadIDtoTaxID), ref(index.iNumOfSpecies), ref(*index.mTaxToIdx));//, ref(transferBetweenRuns->mReadIDToArrayIdx));
							}
							else {
								if (bUnfunny) {
									//foo = bind(&Compare::compareWithDatabase_sloppy<unique_ptr<unique_ptr<const taxaOnly>[]>>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(vLib_taxaOnly), ref(vCount_all), ref(vCount_unique), ref(vReadIDtoTaxID), ref(iAmountOfSpecies), ref(transferBetweenRuns->mReadIDToArrayIdx));
								}
								else {
									unique_ptr<unique_ptr<const contentVecType>[]>* indexPtr = nullptr;
									if (index.bIdentifyMultiple) {
										indexPtr = localIndex->getContentVecType();
									}
									else {
										indexPtr = index.getContentVecType();
									}
									foo = bind(&Compare::compareWithDatabase<unique_ptr<unique_ptr<const contentVecType>[]>>, this, placeholders::_1, ref(vInputVec), placeholders::_2, placeholders::_3, ref(*indexPtr), ref(vCount_all), ref(vCount_unique), ref(vCount_total), ref(vReadIDtoTaxID), ref(index.iNumOfSpecies), ref(*index.mTaxToIdx));//, ref(transferBetweenRuns->mReadIDToArrayIdx));
								}
							}
						}


						// because the stxxl is not threadsafe (as in two threads cannot access the same location on the drive), 
						//   we need to separate the work and make sure that ranges are disjoint for every thread
						size_t iChunkSize = vInputVec.size() / iLocalNumOfThreads;
						size_t iStart = 0, iEnd = iChunkSize + vInputVec.size() % iLocalNumOfThreads;
						for (int32_t iThreadID = 0; iThreadID < iLocalNumOfThreads; ++iThreadID) {
							uint64_t iSeenRange = vInputVec.getRangeStart(iEnd - 1);
							while (iEnd < vInputVec.size()) {
								if (vInputVec.getRangeStart(iEnd) == iSeenRange) {
									iEnd++;
								}
								else {
									break;
								}
							}
							//cout << iStart << " " << iEnd << endl;
							auto task = bind(foo, placeholders::_1, iStart, iEnd);
							workerThreadPool[iThreadID].pushTask(task);
							iStart = iEnd;
							iEnd += iChunkSize;
							if (iEnd > vInputVec.size()) {
								iEnd = vInputVec.size();
							}
						}
						debugBarrier
						/*
						size_t iStart = 0, iEnd = transferBetweenRuns->vRangesOfOutVec[1];
						for (int32_t iThreadID = 0; iThreadID < iLocalNumOfThreads; ++iThreadID) {

							//cout << iStart << " " << iEnd << endl;
							auto task = bind(foo, placeholders::_1, iStart, iEnd);
							workerThreadPool[iThreadID].pushTask(task);
							iStart = iEnd;
							if (iLocalNumOfThreads != 1) {
								iEnd = transferBetweenRuns->vRangesOfOutVec[iThreadID + 2];
							}
						}*/

						
#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "Divide for Compare_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						for (int32_t iThreadID = 0; iThreadID < iLocalNumOfThreads; ++iThreadID) {
							workerThreadPool[iThreadID].startThread();
						}
						for (int32_t iThreadID = 0; iThreadID < iLocalNumOfThreads; ++iThreadID) {
							workerThreadPool[iThreadID].waitUntilFinished();
						}
						debugBarrier
						if (someThingWentWrong) {
							rethrow_exception(someThingWentWrong);
						}

						// Post processing to calculate coherence
						vector<float> vCoherenceScores;
						if (GlobalInputParameters.bPostProcess) {
							vCoherenceScores.resize(iNumOfReads);
							postProcess(vInputVec, vCoherenceScores, iLocalNumOfThreads);
						}
						

#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "PAR Compare_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif

						if (this->_bVisualize) {
							vector<tuple<string,size_t,size_t,uint32_t>> vResultingLinesToPrint;
							vector<pair<uint32_t, uint64_t>> vResultScores;
							for (const auto& entry : this->_translatedFramesForVisualization) {
								cout << entry << endl;
								for (const auto& elem : this->_matchedkMers) {
									string sTempString = Utilities::lstrip(elem.first, '@');
									//cout << sTempString << endl;
									auto pos = entry.find(sTempString);
									if (pos != string::npos) {
										string visualize_outString = "";
										for (size_t i = 0; i < pos; ++i) {
											visualize_outString += " ";
										}
										visualize_outString += sTempString + "," + to_string(elem.second);
										for (size_t i = visualize_outString.length() ; i < entry.length(); ++i) {
											visualize_outString += " ";
										}
										vResultingLinesToPrint.push_back(make_tuple(visualize_outString, pos, sTempString.length(), elem.second));
									}
								}
								debugBarrier
								if (!vResultingLinesToPrint.empty()) {
									sort(vResultingLinesToPrint.begin(), vResultingLinesToPrint.end(), [&](const tuple<string, size_t, size_t, uint32_t>& a, const tuple<string, size_t, size_t, uint32_t>& b) { return make_tuple(get<1>(a), get<2>(a), get<3>(a)) < make_tuple(get<1>(b), get<2>(b), get<3>(b)); });

									for (const auto& printableLines : vResultingLinesToPrint) {
										cout << get<0>(printableLines) << endl;
									}

									sort(vResultingLinesToPrint.begin(), vResultingLinesToPrint.end(), [&](const tuple<string, size_t, size_t, uint32_t>& a, const tuple<string, size_t, size_t, uint32_t>& b) { return get<3>(a) < get<3>(b); });

									uint32_t iSeentaxID = get<3>(vResultingLinesToPrint[0]);
									uint64_t iScore = 0;
									for (const auto& elem : vResultingLinesToPrint) {
										if (iSeentaxID == get<3>(elem)) {
											iScore += get<2>(elem);
										}
										else {
											vResultScores.push_back(make_pair(iSeentaxID, iScore));
											iScore = get<2>(elem);
											iSeentaxID = get<3>(elem);
										}
									}
									vResultScores.push_back(make_pair(iSeentaxID, iScore));
									sort(vResultScores.begin(), vResultScores.end(), [](const pair<uint32_t, uint64_t>& a, const pair<uint32_t, uint64_t>& b) {return a.second > b.second; });
									cout << "Scores: " << endl;
									for (const auto& elem : vResultScores) {
										cout << elem.first << " " << elem.second << endl;
									}
									cout << endl;

									vResultScores.clear();

									vResultingLinesToPrint.clear();
								}
							}
						}
						debugBarrier

						///////////////////////////
						if (bReadIDsAreInteresting) {
							tOutputThread.reset(new thread(&Compare::saveResults, this, transferBetweenRuns->finished, transferBetweenRuns->addTail, ref(vSavedScores), ref(vReadIDtoTaxID), vReadNameAndLength, ref(index.iNumOfSpecies), ref(*index.mIdxToTax), ref(*index.mOrganisms), ref(index.mFrequencies), ref(iNumOfReadsSum), ref(iNumOfReadsOld), iNumOfReads, ref(fThreshold), ref(vReadIDsToBeFiltered), vCoherenceScores, ref(fOut)));
						}

						vCoherenceScores.clear();
						vReadNameAndLength.clear();
						vInputVec.clear();

						if (isGzipped) {
							bIsGood = fast_q_a_File_gz->good();
						}
						else {
							bIsGood = fast_q_a_File->good();
						}
						charsReadOverall = transferBetweenRuns->iNumOfAllCharsRead;

						if (bIsGood && this->bCoverage) {
							cerr << "WARNING: File could not be processed completely in one run. This messes up the coverage values so please increase given RAM or use a smaller file!" << endl;
						}

#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
						cout << "Saving results_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
						startTIME = std::chrono::high_resolution_clock::now();
#endif
						auto timeCompareEnd = std::chrono::high_resolution_clock::now();
						iTimeCompare += chrono::duration_cast<std::chrono::nanoseconds>(timeCompareEnd - timeCompareStart).count();
						//iterate until no dna is left
					}
					debugBarrier
					///////////////////////////////////////////////////////////////////////////////////////////////
					if (tOutputThread) {
						tOutputThread->join();
					}

					// if json is the output format for readToTaxa, end it with a ]
					if (bReadIDsAreInteresting && format == OutputFormat::Json && fOutFile != "") {
						fOut << "\n" << "]";
					}
					if (fOutFile != "") {
						fOut.flush(); // empty the buffer to avoid memory leak
					}
					debugBarrier
					// sum up parallel results
						for (int32_t iThreadID = 1; iThreadID < iLocalNumOfThreads; ++iThreadID) {
							const uint64_t& iStepsize = uint64_t(iThreadID) * Base::_iNumOfK * index.iNumOfSpecies;
							for (uint64_t iIdx = 0; iIdx < uint64_t(Base::_iNumOfK) * index.iNumOfSpecies; ++iIdx) {
								vCount_all[iIdx] += vCount_all[iStepsize + iIdx];
								vCount_unique[iIdx] += vCount_unique[iStepsize + iIdx];
								if (this->bCoverage) {
									vCount_total[iIdx] += vCount_total[iStepsize + iIdx];
								}
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
					debugBarrier
						// get profiling results
						vector<uint64_t> vSumOfUniquekMers(Base::_iNumOfK);
					vector<double> vSumOfNonUniques(Base::_iNumOfK);
					vector<tuple<string, vector<tuple<double, uint64_t, uint64_t>>, uint32_t, uint64_t>> vOut(index.iNumOfSpecies, tuple<string, vector<tuple<double, uint64_t, uint64_t>>, uint32_t, uint64_t>("", vector<tuple<double, uint64_t, uint64_t>>(Base::_iNumOfK), 0, 0));
					for (uint32_t iSpecIdx = 1; iSpecIdx < index.iNumOfSpecies; ++iSpecIdx) {
						vector<tuple<double, uint64_t, uint64_t>> vTemp(Base::_iNumOfK);
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							const uint64_t& iTempScore = vCount_unique[iSpecIdx + uint64_t(ikMerlength) * index.iNumOfSpecies];
							vSumOfUniquekMers[ikMerlength] += iTempScore;
							vSumOfNonUniques[ikMerlength] += vCount_all[iSpecIdx + uint64_t(ikMerlength) * index.iNumOfSpecies];
							if (this->bCoverage) {
								vTemp[ikMerlength] = make_tuple(vCount_all[iSpecIdx + uint64_t(ikMerlength) * index.iNumOfSpecies], iTempScore, vCount_total[iSpecIdx + uint64_t(ikMerlength) * index.iNumOfSpecies]);
							}
							else {
								vTemp[ikMerlength] = make_tuple(vCount_all[iSpecIdx + uint64_t(ikMerlength) * index.iNumOfSpecies], iTempScore, 0);
							}
						}
						vOut[iSpecIdx] = make_tuple(Utilities::replaceCharacter(index.mOrganisms->at(iSpecIdx), ',', ' '), vTemp, index.mIdxToTax->at(iSpecIdx), iSpecIdx);
					}
					sort(vOut.begin(), vOut.end(), [](const tuple<string, vector<tuple<double, uint64_t, uint64_t>>, uint32_t, uint64_t>& a, const tuple<string, vector<tuple<double, uint64_t, uint64_t>>, uint32_t, uint64_t>& b) {
						for (uint64_t i = 0; i < get<1>(a).size(); ++i) {
							if (get<1>(get<1>(a)[i]) == get<1>(get<1>(b)[i])) {
								continue;
							}
							else {
								return get<1>(get<1>(a)[i]) > get<1>(get<1>(b)[i]);
							}
						}
						return false;
						});

					debugBarrier
						// Because of the padding at the end to get useful k-mers between max_k and min_k (see sFalsekMerMarker), overall relative frequency is distorted
						// therefore we need to calculate how many there were to get the number right, how many could have possibly matched
						vector<uint64_t> vNumberOfGarbagekMersPerK(Base::_iNumOfK, 0);
					const uint64_t iFrameMultiplier = (this->_bOnlyOneFrame) ? 1 : ((Base::_bSixFrames) ? 6 : 3);
					for (int32_t i = Base::_iMaxK - Base::_iMinK, j = 0; i > 0; --i, ++j) {
						vNumberOfGarbagekMersPerK[j] = iNumOfReadsSum * iFrameMultiplier * i; //Number of reads * 2(forward,reverse complement) * 3(frames) * dist(max_k,min_k)
					}

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
								if (get<1>(entry)[Base::_iNumOfK - 1].first > 0 && !bBreakOut) {
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
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Unique counts k=" << Base::_iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Unique rel. freq. k=" << Base::_iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Non-unique counts k=" << Base::_iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Non-unique rel. freq. k=" << Base::_iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Overall rel. freq. k=" << Base::_iMaxK - ikMerlength;
						}
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							tableFileStream << "," << "Overall unique rel. freq. k=" << Base::_iMaxK - ikMerlength;
						}
						if (this->bCoverage) {
							for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
								tableFileStream << "," << "Special Counts k=" << Base::_iMaxK - ikMerlength;
							}
							for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
								tableFileStream << "," << "Genome Coverage k=" << Base::_iMaxK - ikMerlength;
							}
						}
						tableFileStream << "\n";

						vector<double> vSumOfIdentified(Base::_iNumOfK, 0), vSumOfUniqueIdentified(Base::_iNumOfK, 0);
						uint64_t iSpecIdx = 1;
						for (const auto& entry : vOut) {
							if (get<0>(get<1>(entry)[Base::_iNumOfK - 1]) > 0) {
								// unique count
								sOutStr << get<2>(entry) << "," << get<0>(entry);
								for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
									sOutStr << "," << get<1>(get<1>(entry)[ikMerlength]);
								}
								// unique rel freq
								for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
									if (get<1>(get<1>(entry)[ikMerlength]) == 0) {
										sOutStr << "," << 0.0;
									}
									else {
										sOutStr << "," << static_cast<double>(get<1>(get<1>(entry)[ikMerlength])) / vSumOfUniquekMers[ikMerlength];
									}
								}
								// non-unique count
								for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
									sOutStr << "," << get<0>(get<1>(entry)[ikMerlength]);
								}
								// non-unique rel freq
								for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
									if (get<0>(get<1>(entry)[ikMerlength]) == 0) {
										sOutStr << "," << 0.0;
									}
									else {
										//sOutStr << "," << static_cast<double>(get<1>(entry)[ikMerlength]) / (iNumberOfkMersInInput - aNonUniqueHits[ikMerlength] - (Base::_iMaxK - Base::_iMinK - ikMerlength) * 6 * iNumOfReadsSum);
										sOutStr << "," << static_cast<double>(get<0>(get<1>(entry)[ikMerlength])) / vSumOfNonUniques[ikMerlength];
									}
								}
								// Overall rel freq
								for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
									vSumOfIdentified[ikMerlength] += get<0>(get<1>(entry)[ikMerlength]);
									sOutStr << "," << static_cast<double>(get<0>(get<1>(entry)[ikMerlength])) / (iNumberOfkMersInInput - vNumberOfGarbagekMersPerK[ikMerlength]);
								}

								// Overall unique rel freq
								for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
									vSumOfUniqueIdentified[ikMerlength] += get<1>(get<1>(entry)[ikMerlength]);
									sOutStr << "," << static_cast<double>(get<1>(get<1>(entry)[ikMerlength])) / (iNumberOfkMersInInput - vNumberOfGarbagekMersPerK[ikMerlength]);
								}
								if (this->bCoverage) {
									// Special count
									for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
										sOutStr << "," << get<2>(get<1>(entry)[ikMerlength]);
									}

									// Genome coverage
									for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
										sOutStr << "," << static_cast<double>(get<2>(get<1>(entry)[ikMerlength])) / index.mFrequencies[get<3>(entry)][ikMerlength];
									}
								}
								sOutStr << "\n";
							}
							iSpecIdx++;
						}
						// last entry
						tableFileStream << "0,not identified";
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK * 4; ++ikMerlength) {
							// all counts for unique and non-unique relate to the identified number of counts so no value other than 0 can be written here
							tableFileStream << "," << 0.0;
						}

						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							//cout << iNumberOfkMersInInput << " " << vNumberOfGarbagekMersPerK[ikMerlength] << " " << vSumOfIdentified[ikMerlength] << endl;
							tableFileStream << "," << (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]) - static_cast<double>(vSumOfIdentified[ikMerlength])) / (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]));
						}
						for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK; ++ikMerlength) {
							tableFileStream << "," << (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]) - static_cast<double>(vSumOfUniqueIdentified[ikMerlength])) / (static_cast<double>(iNumberOfkMersInInput) - static_cast<double>(vNumberOfGarbagekMersPerK[ikMerlength]));
						}
						if (this->bCoverage) {
							for (int32_t ikMerlength = 0; ikMerlength < Base::_iNumOfK * 2; ++ikMerlength) {
								tableFileStream << "," << 0.0;
							}
						}
						tableFileStream << "\n" << sOutStr.str();

						allSumOfIdentified = vSumOfIdentified[Base::_iNumOfK - 1];
						//}
					}
					else {
						for (const auto& entry : vOut) {
							allSumOfIdentified += get<0>(get<1>(entry)[Base::_iNumOfK - 1]);
						}
							}
					debugBarrier

#ifdef TIME
						endTIME = std::chrono::high_resolution_clock::now();
					cout << "Save profile_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
					startTIME = std::chrono::high_resolution_clock::now();
#endif

					/*if (fTableFile != "") {
						cout.rdbuf(orgBuf);
					}*/
					if (Base::_bVerbose) {
						cout << "OUT: Number of k-mers in input: " << iNumberOfkMersInInput << " of which " << allSumOfIdentified / iNumberOfkMersInInput * 100. << " % were identified." << endl;
						cout << "OUT: Number of uniques:";
						for (int32_t j = 0; j < Base::_iNumOfK; ++j) {
							cout << " " << vSumOfUniquekMers[j];
						}
						cout << endl;
						cout << "OUT: Time fastq: " << iTimeFastq << " ns" << endl;
						cout << "OUT: Time compare: " << iTimeCompare << " ns" << endl;
					}


					// filtering
					if (GlobalInputParameters.bFilter) {
						if (isGzipped) {
							// reset input files
							fast_q_a_File_gz.reset(new igzstream());
							fast_q_a_File_gz->rdbuf()->pubsetbuf(&inFilebuffer[0], 2097152);
							fast_q_a_File_gz->open(inFile.first.c_str());

							if (bIsPairedEnd) {
								fast_q_a_File_gz2.reset(new igzstream());
								fast_q_a_File_gz2->rdbuf()->pubsetbuf(&inFilebuffer2[0], 2097152);
								fast_q_a_File_gz2->open(vPairedEndFiles[1].c_str());
							}

							if (GlobalInputParameters.bGzipOut) {
								filter<igzstream,ogzstream>(vReadIDsToBeFiltered, fast_q_a_File_gz, fast_q_a_File_gz2, bIsFasta, bIsPairedEnd);
							}
							else {
								filter<igzstream, ofstream>(vReadIDsToBeFiltered, fast_q_a_File_gz, fast_q_a_File_gz2, bIsFasta, bIsPairedEnd);
							}
						}
						else {
							fast_q_a_File.reset(new ifstream());
							fast_q_a_File->rdbuf()->pubsetbuf(&inFilebuffer[0], 2097152);
							fast_q_a_File->open(inFile.first.c_str());

							if (bIsPairedEnd) {
								fast_q_a_File2.reset(new ifstream());
								fast_q_a_File2->rdbuf()->pubsetbuf(&inFilebuffer2[0], 2097152);
								fast_q_a_File2->open(vPairedEndFiles[1].c_str());
							}

							if (GlobalInputParameters.bGzipOut) {
								filter<ifstream, ogzstream>(vReadIDsToBeFiltered, fast_q_a_File, fast_q_a_File2, bIsFasta, bIsPairedEnd);
							}
							else {
								filter<ifstream, ofstream>(vReadIDsToBeFiltered, fast_q_a_File, fast_q_a_File2, bIsFasta, bIsPairedEnd);
							}
						}
						vReadIDsToBeFiltered.clear();
					}

					iNumOfReadsSum = 0;
					iNumberOfkMersInInput = 0;
					iTimeFastq = 0;
					iTimeCompare = 0;
					vReadNameAndLength.clear();


					fill_n(vCount_all.get(), iMult, 0.);
					fill_n(vCount_unique.get(), iMult, 0);
					if (this->bCoverage) {
						fill_n(vCount_total.get(), iMult, 0);
					}

					allFilesProgress = transferBetweenRuns->iCurrentOverallPercentage;
				}
				debugBarrier
#ifdef TIME
				endTIME = std::chrono::high_resolution_clock::now();
				cout << "End_" << chrono::duration_cast<std::chrono::nanoseconds>(endTIME - startTIME).count() << endl;
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