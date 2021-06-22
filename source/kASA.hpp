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

#include "MetaHeader.h"

#include "utils/Utilities.hpp"
#include "utils/iToStr.hpp"


namespace kASA {


	//////////////////////////////////////////////
	// typedefs can be found in "Metaheader.h"
	class kASA {

	public: //protected:

		struct SCompareStructForSTXXLSort {
			bool operator() (const packedBigPair& a, const packedBigPair& b) const {
				return a < b;
			}

			packedBigPair min_value() const { packedBigPair t; t.first = numeric_limits<uint64_t>::min(); t.second = numeric_limits<uint32_t>::min(); return t; }
			packedBigPair max_value() const { packedBigPair t; t.first = numeric_limits<uint64_t>::max(); t.second = numeric_limits<uint32_t>::max(); return t; }
		};


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Variables
	public:
		string _sTemporaryPath = string("");

		int32_t _iNumOfThreads;
		int32_t _iNumOfCall;

		uint32_t _iMaxMemUsePerThread;

		int32_t _iMaxK, _iMinK, _iNumOfK, _iHighestK = 12, _iLowestK = 1;
		unique_ptr<int32_t[]> _aOfK;
		string _sMaxKBlank = string("");

		static int8_t _sAminoAcids_bs[], _sAminoAcids_aas[];
		const static int8_t _sAminoAcids_cs[], _sAminoAcids_un[];
		const int8_t _aRevComp[6] = { 'T','G','A','C','X','Z' };
		//const string _spacedMasks[21] = { "1111111111111*111111111111111111111*11111111111111111111111111111111111111111", "1111111111111*11111111111*111111111*1111111111111111111*111111111111111111111", "1111111111111*1111111*111*11111*111*1111111111111111111*111111111111111111111", "1111111111111*1111111*111*11111*111*11111*111111*111111*111111111111111111111", "1111111*11111*1111111*111*11111*111*11111*111111*111111*11111*111111111111111", "1111111*11111*1111111*111*11*11*111*11*11*111111*111111*11111*111111111111111", "1111111*111*1*1111111*111*11*11*111*11*11*11*111*111111*11111*111111111111111", "1111111*111*1*1111*11*111*11*11*111*11*11*11*111*111*11*11111*111111111111111", "1111*11*111*1*1111*11*111*11*11*111*11*11*11*111*111*11*1*111*111111111111111", "1111*111*111**1*111**1*11*111111111*111*111**1*111**1*11*111111111*111*111*11", "1111*11*111*1*1111*11*111*11*11*111*11*11*11*111*111*11*1*11**111111*11111111", "1111*11*111*1*111**11*111*11*1**111*11*11*11*111*111*11*1*11**111111*11111111", "1111*11*111*1**11**11*111*11*1**111*11*11*11*111*111*11*1*11**1111*1*11111111", "1111*11*111*1**11**11*111*11*1**11**11*11*11*111*111*11*1*11**11*1*1*11111111", "1111*11**11*1**11**11*111*11*1**11**11*11*11*11**111*11*1*11**11*1*1*11111111", "1111*11**11*1**11**11*111*11*1**11**11*11*11*11**11**11*1*11**11*1*1*1111*111", "1111*11***1*1**11**11*111*11*1**11**11*11*11*11**11**11*1*11**11*1*1**111*111", "1111*11***1*1**11**11**11*11*1**11**11*11*11*11**11**11***11**11*1*1**111*111", "1111*11***1*1**11***1**11*11*1**11**11*11*11*11**11**1****11**11*1*1**111*111" };
		int16_t _iWhichMask = 0;

		const bool _bVerbose = false;
		bool _bSixFrames = false;
		bool _bProtein = false;
		bool _bSpaced = false;
		
		// Functions

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Converts a string of dna to a string of aminoacids by triplets;
		// sDNA is a reference to the dna string, length is the input length, pointer position is the starting position in the dna string
		inline void dnaToAminoacid(const string& sDna, const int32_t& iLength, const int32_t& iPointerPosition, string* outString) {
			// map A,C,T,G,X,Z to 0,1,2,3,4,5 by & with 1110 = 14 and then shift to the position in 000 000 000, this gives the position in the array where the respective amino acid is written
			// small letters are also correctly converted
			const uint32_t& iDnaSize = iLength / 3;

			for (uint32_t i = 0, j = 0; i < iDnaSize; ++i, j += 3) {
				const int32_t& iIndex = ((sDna[iPointerPosition + j] & 14) << 5) | ((sDna[iPointerPosition + j + 1] & 14) << 2) | ((sDna[iPointerPosition + j + 2] & 14) >> 1);
				(*outString)[i] = _sAminoAcids_bs[iIndex];
			}

		}

		inline void dnaToAminoacid(const string& sDna, const int32_t& iPointerPosition, int8_t& outString) {
			// map A,C,T,G,X,Z to 0,1,2,3,4,5 by & with 1110 = 14 and then shift to the position in 000 000 000, this gives the position in the array where the respective amino acid is written
			// small letters are also correctly converted

			const int32_t& iIndex = ((sDna[iPointerPosition] & 14) << 5) | ((sDna[iPointerPosition + 1] & 14) << 2) | ((sDna[iPointerPosition + 2] & 14) >> 1);
			outString = _sAminoAcids_bs[iIndex];
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/*inline char getLetterOfMask(const string& sDNA, const uint32_t& iIdx) {
			const uint32_t& spacedMaskIdx = iIdx % (3* HIGHESTPOSSIBLEK);
			if (_spacedMasks[_iWhichMask][spacedMaskIdx] == '1') {
				if (_spacedMasks[_iWhichMask][spacedMaskIdx +1] == '1') {
					if (_spacedMasks[_iWhichMask][spacedMaskIdx +2] == '1') {
						const int32_t& iIndex = ((sDNA[iIdx] & 14) << 5) | ((sDNA[iIdx + 1] & 14) << 2) | ((sDNA[iIdx + 2] & 14) >> 1);
						return _sAminoAcids_bs[iIndex]; /// 111
					}
					else {
						return 'Z'; // 11*
					}
				}
				else {
					if (_spacedMasks[_iWhichMask][spacedMaskIdx + 2] == '1') {
						return '\\'; // 1*1
					}
					else {
						return 'X'; // 1**
					}
				}
			}
			else {
				if (_spacedMasks[_iWhichMask][spacedMaskIdx + 1] == '1') {
					if (_spacedMasks[_iWhichMask][spacedMaskIdx + 2] == '1') {
						return 'O'; // *11
					}
					else {
						return 'U'; // *1*
					}
				}
				else {
					if (_spacedMasks[_iWhichMask][spacedMaskIdx + 2] == '1') {
						return 'J'; // **1
					}
					else {
						return 'B'; // ***
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Converts a string of dna to a string of aminoacids by triplets with a spaced mask
		inline void dnaToAminoacidSpaced(const string& sDna, const int32_t& iLength, const int32_t& iPointerPosition, string* outString) {
			const uint32_t& iDnaSize = iLength / 3;

			for (uint32_t i = 0, j = 0; i < iDnaSize; ++i, j += 3) {
				(*outString)[i] = getLetterOfMask(sDna, iPointerPosition + j);
			}
		}

		inline void dnaToAminoacidSpaced(const string& sDna, const int32_t& iPointerPosition, int8_t& outString) {
			outString = getLetterOfMask(sDna, iPointerPosition);
		}*/
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// translates with and sets up the lookup table for the sloppy mode
		inline uint64_t aminoAcidsToAminoAcid(const uint64_t& kmer) {
			uint64_t iOut = 0;
			for (int i = 0, j = 0; i < 12; i += 2, ++j) {
				const int32_t& iShift = 5 * (10 - i);
				// the &31 is necessary because the letters in the lookpup table are in the column of printable letters of the ascii table
				// 1023 gets two letters which are joined via the lookup table
				// the last shift ensures, that the 6 letter word is on the left side of the integer
				iOut |= static_cast<uint64_t>(_sAminoAcids_aas[(kmer & (1023ULL << iShift)) >> iShift] & 31) << (55 - j * 5);
			}
			return iOut;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// detect which alphabet and set variable accordingly
		inline void detectAlphabet(const string& s) {
			regex self_regex_dna("^[ACGTURYKMSWBDHVN-]+$",  std::regex_constants::icase);
			if (std::regex_search(s, self_regex_dna)) {
				if (_bVerbose) {
					cout << "OUT: DNA sequences detected." << endl;
				}
				_bProtein = false;
			}
			else {
				regex self_regex_p("^[ABCDEFGHIJKLMNOPQRSTUVWXYZ*-]+$", std::regex_constants::icase);
				if (!std::regex_search(s, self_regex_p)) {
					cerr << "ERROR: The sequence is neither recognized as protein nor DNA. It will be treated as protein sequence but it may fail... Sequence was: " << s << endl;
				}
				else {
					if (_bVerbose) {
						cout << "OUT: Protein sequences detected." << endl;
					}
				}

				_bProtein = true;
				_bSixFrames = false;
			}
		}

		static inline void setAAToAATable(const string& file) {
			string sDummy = "";
			ifstream mappings(file); // C:/Users/Silvio/Desktop/out.csv
			if (!mappings) {
				throw runtime_error("Mapping file couldn't be opened for reading!");
			}
			while (getline(mappings, sDummy)) {
				if (sDummy != "") {
					auto spl = Utilities::split(Utilities::removeSpaceAndEndline(sDummy), ',');
					const uint32_t& index = ((static_cast<uint32_t>(spl[0][0]) & 31) << 5) | (static_cast<uint32_t>(spl[0][1]) & 31);
					_sAminoAcids_aas[index] = static_cast<uint32_t>(spl[1][0]) & 31;
				}
			}
#ifdef _DEBUG
			ofstream dummyOut("I:/DA/lookupJ.txt");
			for (int i = 0; i < 900; ++i) {
				if (static_cast<char>(_sAminoAcids_aas[i] + 64) == '\\') {
					dummyOut << "'\\', ";
				}
				else {
					dummyOut << "'" << static_cast<char>(_sAminoAcids_aas[i] + 64) << "', ";
				}
			}
#endif
		}
	protected:

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Creates the reverse complement of a DNA string
		inline string reverseComplement(const string& sDNA) {
			const int32_t& iDNALength = int32_t(sDNA.length());
			string outString(iDNALength, ' ');
			for (int32_t i = iDNALength - 1, j = 0; i >= 0; --i, ++j) {
				outString[i] = _aRevComp[(sDNA[j] >> 1) & 7];
			}
			return outString;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Creates default config for stxxl, only called once
		inline void createConfig(const string& tempFileName, const int32_t& iNumCall, const string& mode = "", const bool& bDelete = true) {
#if (_WIN32 || _WIN64) && !__APPLE__
			string IOCall = (mode != "") ? mode : "wincall autogrow";
#endif
#if (__GNUC__ || defined(__llvm__)) && !_MSC_VER
			string IOCall = (mode != "") ? mode : "syscall unlink autogrow"; 
#endif
			if (bDelete) {
				IOCall += " delete";
			}


			int32_t localNumCall = iNumCall;
			while (ifstream(_sTemporaryPath + tempFileName + to_string(localNumCall) + ".tmp")) {
				++localNumCall;
			}

			stxxl::config* cfg = stxxl::config::get_instance();
			if (bDelete) {
				stxxl::disk_config disk(_sTemporaryPath + tempFileName + to_string(localNumCall) + ".tmp", 1024 * 1024, IOCall);
				disk.direct = stxxl::disk_config::DIRECT_ON;
				cfg->add_disk(disk);
			}
			else {
				stxxl::disk_config disk(_sTemporaryPath + tempFileName + to_string(localNumCall) + ".tmp", 1024 * 1024, IOCall);
				disk.direct = stxxl::disk_config::DIRECT_ON;
				cfg->add_disk(disk);
			}
		}


	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Standard Constructor
		kASA() : _bVerbose(false), _bSixFrames(false) {
			_sTemporaryPath = "";
			_iNumOfThreads = 1;
			_iNumOfCall = 1;
			_iMaxK = 1;
			_iMinK = 1;
			_iNumOfK = 1;

			_aOfK.reset(new int32_t[_iNumOfK]);
			for (int32_t i = 0; i < _iNumOfK; ++i) {
				_aOfK[i] = _iMaxK - i;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Constructor with paths e.g. "derp/"
		kASA(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHighestK, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const string& stxxl_mode = "", const bool& bSixFrames = false) : _sTemporaryPath(tmpPath), _iNumOfThreads(iNumOfProcs), _iNumOfCall(iNumOfCall), _iMaxK(iHigherK), _iMinK(iLowerK), _iNumOfK(_iMaxK - _iMinK + 1), _bVerbose(bVerbose), _bSixFrames(bSixFrames) {
#ifdef ENVIRONMENT32
			_iMaxMemUsePerThread = (1024 / _iNumOfThreads) * 1024 * 1024;
#else
			_iMaxMemUsePerThread = (3072 / iNumOfProcs) * 1024 * 1024;
#endif

#if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
			omp_set_num_threads(_iNumOfThreads);
			omp_set_nested(1);
			omp_set_dynamic(_iNumOfThreads);
#endif
			_iHighestK = iHighestK;
			_iMaxK = (iHigherK <= _iHighestK && iHigherK >= iLowerK) ? iHigherK : _iHighestK;
			_iMinK = (iLowerK < _iLowestK) ? _iLowestK : _iMinK;
			_iNumOfK = _iMaxK - _iMinK + 1;

			_sMaxKBlank = string(_iHighestK, ' ');

			_aOfK.reset(new int32_t[_iNumOfK]);
			for (int32_t i = 0; i < _iNumOfK; ++i) {
				_aOfK[i] = _iMaxK - i;
			}

			createConfig("stxxl_temp_", iNumOfCall, stxxl_mode);
		}
		
		kASA(const kASA& obj) : _bVerbose(obj._bVerbose) {
			_iMaxMemUsePerThread = obj._iMaxMemUsePerThread;
			_iHighestK = obj._iHighestK;
			_iMaxK = obj._iMaxK;
			_iMinK = obj._iMinK;
			_iNumOfK = obj._iNumOfK;
			_aOfK.reset(new int32_t[_iNumOfK]);
			for (int32_t i = 0; i < _iNumOfK; ++i) {
				_aOfK[i] = _iMaxK - i;
			}

			//config already created
			_iNumOfCall = obj._iNumOfCall;
			_sTemporaryPath = obj._sTemporaryPath;
			_iNumOfThreads = obj._iNumOfThreads;
			_sMaxKBlank = obj._sMaxKBlank;
			_iWhichMask = obj._iWhichMask;
			_bSixFrames = obj._bSixFrames;
			_bProtein = obj._bProtein;
			_bSpaced = obj._bSpaced;
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Converts one aminoacid string to a coded kMer, using 5 bits per amino acid letter (max k = 12 => 60 bits used)
		template<class intType>
		static inline intType aminoacidTokMer(const string& aminoacid) {
			intType ikMer = 0;

			const int32_t& iLengthAA = int32_t(aminoacid.length());
			for (int32_t i = 0; i < iLengthAA - 1; ++i) {
				uint32_t iConvertedChar = aminoacid[i] & 31;
				ikMer |= iConvertedChar;
				ikMer <<= 5;
			}
			uint32_t iConvertedChar = aminoacid[iLengthAA - 1] & 31;
			ikMer |= iConvertedChar;

			return ikMer;
		}

		template<class intType>
		inline intType aminoacidTokMer(const string::const_iterator& begin, const string::const_iterator& end) {
			intType ikMer = 0;

			const int32_t& iLengthAA = int32_t(end - begin);
			for (int32_t i = 0; i < iLengthAA - 1; ++i) {
				uint32_t iConvertedChar = *(begin + i) & 31;
				ikMer |= iConvertedChar;
				ikMer <<= 5;
			}
			uint32_t iConvertedChar = *(begin + iLengthAA - 1) & 31;
			ikMer |= iConvertedChar;

			return ikMer;
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Concatenates one aminoacid (int_8_t) with a coded kMer
		inline uint64_t aminoacidTokMer(const uint64_t& kMer, const int8_t& aa) {
			/// shift away the first letter and gather only the relevant 60 bits that remain, then copy the 5 relevant bits of the input onto it
			uint64_t ikMer = (0 | (kMer << 5)) & 0xFFFFFFFFFFFFFFF;
			ikMer |= aa & 31;
			return ikMer;
		}

		inline uint128_t aminoacidTokMer(const uint128_t& kMer, const int8_t& aa) {
			/// shift away the first letter and gather only the relevant bits that remain, then copy the 5 relevant bits of the input onto it
			uint128_t ikMer = (uint128_t(0) | (kMer << 5)) & uint128_t(0x1FFFFFFFFFFFFFFF, numeric_limits<uint64_t>::max());
			ikMer |= aa & 31;
			return ikMer;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Reverts a coded kMer to the respective aminoacid string
		static inline string kMerToAminoacid(const int64_t& kMer, const int32_t& iMaxK) {
			string aacid = "";
			int32_t counter = iMaxK - 1;
			const int8_t& bitmask = 31;
			while (counter >= 0) {
				int8_t&& toBeChar = int8_t(kMer >> counter * 5);
				toBeChar &= bitmask;
				toBeChar |= 64;
				aacid += toBeChar;
				counter--;
			}

			return aacid;
		}

		static inline string kMerToAminoacid(const uint128_t& kMer, const int32_t& iMaxK) {
			string aacid = "";
			int32_t counter = iMaxK - 1;
			const int8_t& bitmask = 31;
			while (counter >= 0) {
				uint8_t&& toBeChar = uint8_t(kMer >> counter * 5);
				toBeChar &= bitmask;
				toBeChar |= 64;
				aacid += toBeChar;
				counter--;
			}

			return aacid;
		}


		template<typename T> inline static void showVec(const T& in, const uint64_t& iStartIdx = 0) {
			uint32_t iCounter = 0;
			string sDummy = "", lookup = "";
			for (auto elem = in.begin() + iStartIdx; elem != in.end(); ++elem) {
				if (iCounter == 20) {
					iCounter = 0;
					if (lookup == "") {
						cin >> sDummy;
					}
					if (sDummy == "q" || sDummy == "Q") {
						return;
					}
					if (sDummy == "l") {
						cin >> lookup;
					}
					if (sDummy == "e") {
						elem = in.end() - 20;
					}
				}
				if (lookup != "") {
					if (kMerToAminoacid(elem->first, (is_same<T,contentVecType_128>::value) ? HIGHESTPOSSIBLEK : 12) == lookup) {
						cout << elem->first << " " << kMerToAminoacid(elem->first, (is_same<T, contentVecType_128>::value) ? HIGHESTPOSSIBLEK : 12) << " " << elem->second << endl;
						lookup = "";
					}
				}
				else {
					cout << elem->first << " " << kMerToAminoacid(elem->first, (is_same<T, contentVecType_128>::value) ? HIGHESTPOSSIBLEK : 12) << " " << elem->second << endl;
					++iCounter;
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Counts the number of kMers for each taxon (frequency). 
		// Usually this is counted while creating the index. If for some reason the frequency file got lost, it can be recreated with this.
		template<class vecType>
		inline void GetFrequencyK(const string& contentFile, const string& sLibFile, const int64_t& iMemoryAvail, const unordered_map<uint32_t, uint32_t>& mGivenContentMap = unordered_map<uint32_t, uint32_t>()) {
			try {
				// test if files exists
				if (!ifstream(contentFile) || !ifstream(sLibFile)) {
					throw runtime_error("OUT: One of the files does not exist");
				}

				int64_t iMemoryUsage = 0;
				const int32_t& iMaxNumK = _iHighestK - _iLowestK + 1;

				ifstream fContent(contentFile);
				uint32_t iIdxCounter = 1;
				const unordered_map<uint32_t, uint32_t>* mContentThatWillBeUsed = nullptr;
				unordered_map<uint32_t, uint32_t> mContent; mContent[0] = 0;
				if (mGivenContentMap.empty()) {
					while (fContent.good()) {
						string dummy = "";
						getline(fContent, dummy);
						if (dummy != "") {
							const auto& line = Utilities::split(dummy, '\t');
							mContent[stoul(line[1])] = iIdxCounter;
							++iIdxCounter;
						}
					}
					mContentThatWillBeUsed = &mContent;
				}
				else {
					iIdxCounter = static_cast<uint32_t>(mGivenContentMap.size());
					mContentThatWillBeUsed = &mGivenContentMap;
				}
				fContent.close();
				fContent.clear();

				iMemoryUsage += Utilities::calculateSizeInByteOfUnorderedMap(*mContentThatWillBeUsed);

				// Determine how many threads can be created
				int32_t iLocalNumOfThreads = 0;
				while (iMemoryUsage < iMemoryAvail && iLocalNumOfThreads < _iNumOfThreads) {
					iMemoryUsage += iIdxCounter * iMaxNumK * sizeof(uint64_t);
					iMemoryUsage += vecType::block_size * vecType::page_size * 4;
					iLocalNumOfThreads++;
				}
				if (iLocalNumOfThreads == 1 && iMemoryUsage > iMemoryAvail) {
					cerr << "OUT: WARNING! Due to the large content file, creating the frequency file will consume more memory than given. kASA may crash or slow down..." << endl;
				}
				

				ifstream fInfo(sLibFile + "_info.txt");
				uint64_t iSizeOfVec = 0;
				fInfo >> iSizeOfVec;
				stxxlFile libfile(sLibFile, stxxl::file::RDONLY);
				unique_ptr<unique_ptr<const vecType>[]> libvec(new unique_ptr<const vecType>[iLocalNumOfThreads]);
				for (int32_t i = 0; i < iLocalNumOfThreads; ++i) {
					libvec[i].reset(new const vecType(&libfile, iSizeOfVec));
				}
				
				unique_ptr<uint64_t[]> aFrequencyArray(new uint64_t[size_t(iLocalNumOfThreads) * iIdxCounter * iMaxNumK]);
				memset(aFrequencyArray.get(), 0, size_t(iLocalNumOfThreads) * iIdxCounter * iMaxNumK * sizeof(uint64_t));

				// count up
				vector<thread> workerThreads;
				uint64_t iPartialSize = iSizeOfVec / iLocalNumOfThreads;
				uint64_t iPartialRemainder = iSizeOfVec % iLocalNumOfThreads;
				auto func = [&](const int32_t iThreadID) {
					for (uint64_t i = iThreadID * iPartialSize; i < (iThreadID + 1) * iPartialSize + (((iThreadID + 1) * iPartialSize + iPartialRemainder == iSizeOfVec) ? iPartialRemainder : 0); ++i) {
						const auto& entry = (libvec[iThreadID])->at(i);
						const auto& idx = (iThreadID * iIdxCounter + Utilities::checkIfInMap(*mContentThatWillBeUsed, entry.second)->second) * uint64_t(iMaxNumK);
						for (int32_t k = 0; k < iMaxNumK; ++k) {
							const int32_t& shift = 5 * k;
							if (((entry.first >> shift) & 31) != 30) {
								++(aFrequencyArray[idx + k]);
							}
						}
					}
				};

				for (int32_t iThreadID = 0; iThreadID < iLocalNumOfThreads; ++iThreadID) {
					workerThreads.push_back(thread(func, iThreadID));
				}
				for (int32_t iThreadID = 0; iThreadID < iLocalNumOfThreads; ++iThreadID) {
					if (workerThreads[iThreadID].joinable()) {
						workerThreads[iThreadID].join();
					}
				}

				// Gather everything
				for (int32_t i = 1; i < iLocalNumOfThreads; ++i) {
					for (uint32_t j = 0; j < iIdxCounter; ++j) {
						for (int32_t k = 0; k < iMaxNumK; ++k) {
							aFrequencyArray[j * iMaxNumK + k] += aFrequencyArray[(i * iIdxCounter + j) * iMaxNumK + k];
						}
					}
				}

				// Write to file
				fContent.open(contentFile);
				ofstream outFile(sLibFile + "_f.txt");
				if (!outFile) {
					throw runtime_error("Frequency file couldn't be opened for writing!");
				}
				// First line
				outFile << "non_unique";
				for (int32_t k = 0; k < iMaxNumK; ++k) {
					outFile << "\t" << aFrequencyArray[k];
				}
				outFile << endl;

				// rest
				string sDummy = "";
				for (uint32_t j = 1; j < iIdxCounter && getline(fContent, sDummy); ++j) {
					const auto& line = Utilities::split(sDummy, '\t');
					outFile << line[0] << "\t";
					outFile << aFrequencyArray[j * iMaxNumK];
					for (int32_t k = 1; k < iMaxNumK; ++k) {
						outFile << "\t" << aFrequencyArray[j * iMaxNumK + k];
					}
					outFile << endl;
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Use a different codon table
		inline void setCodonTable(const string& codonFile, const string& iID) {
			try {
				ifstream ncbiFile(codonFile);

				string line = "";

				bool bFound = false;
				while (getline(ncbiFile, line)) {
					if (line.find("  id " + iID + " ,") != string::npos) {
						bFound = true;
						break;
					}
				}
				
				if (bFound) {
					string aminoAcids = "", Base1 = "", Base2 = "", Base3 = "", dummy = "";
					getline(ncbiFile, aminoAcids);
					getline(ncbiFile, dummy);
					getline(ncbiFile, Base1);
					getline(ncbiFile, Base2);
					getline(ncbiFile, Base3);

					size_t iPosAA = aminoAcids.find_first_of('"') + 1, iPosB = Base1.find_first_of("TGCA"); // first match of any of the bases suffices
					for (; iPosB < Base1.size(); ++iPosB, ++iPosAA) {
						//shift to correct position just like in dnaToAminoAcid
						_sAminoAcids_bs[((Base1[iPosB] & 14) << 5) | ((Base2[iPosB] & 14) << 2) | ((Base3[iPosB] & 14) >> 1)] = (aminoAcids[iPosAA] == '*') ? '[' : int8_t(aminoAcids[iPosAA]);
					}

				}
				else {
					cerr << "WARNING: codon table not found in file. Using built-in." << endl;
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

	};

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Lookup table for Amino Acids, used in dnaToAminoAcid, can be modified via file
	int8_t kASA::_sAminoAcids_bs[] = {
		'K', 'N', 'N', 'K', '^', '_', ' ', ' ',  // AAA	AAC	AAT	AAG	AAX	AAZ
		'T', 'T', 'T', 'T', '^', '_', ' ', ' ',	 // ACA ACC ACT ACG	ACX ACZ
		'I', 'I', 'I', 'M', '^', '_', ' ', ' ',	 // ATA ATC ATT ATG ATX ATZ
		'R', 'S', 'S', 'R', '^', '_', ' ', ' ',  // AGA AGC AGT AGG AGX AGZ
		'^', '^', '^', '^', '^', '_', ' ', ' ',  // AXA AXC AXT AXG AXX AXZ
		'_', '_', '_', '_', '_', '_', ' ', ' ',  // AZA AZC AZT AZG AZX AZZ
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',  // Empty spaces because of the gap towards the next number
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		'Q', 'H', 'H', 'Q', '^', '_', ' ', ' ',  // CAA CAC CAT ...
		'P', 'P', 'P', 'P', '^', '_', ' ', ' ',
		'L', 'L', 'L', 'L', '^', '_', ' ', ' ',
		'R', 'R', 'R', 'R', '^', '_', ' ', ' ',
		'^', '^', '^', '^', '^', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		'[', 'Y', 'Y', '[', '^', '_', ' ', ' ',
		'S', 'S', 'S', 'S', '^', '_', ' ', ' ',
		'L', 'F', 'F', 'L', '^', '_', ' ', ' ',
		']', 'C', 'C', 'W', '^', '_', ' ', ' ',
		'^', '^', '^', '^', '^', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		'E', 'D', 'D', 'E', '^', '_', ' ', ' ',
		'A', 'A', 'A', 'A', '^', '_', ' ', ' ',
		'V', 'V', 'V', 'V', '^', '_', ' ', ' ',
		'G', 'G', 'G', 'G', '^', '_', ' ', ' ',
		'^', '^', '^', '^', '^', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		'^', '^', '^', '^', '^', '_', ' ', ' ',
		'^', '^', '^', '^', '^', '_', ' ', ' ',
		'^', '^', '^', '^', '^', '_', ' ', ' ',
		'^', '^', '^', '^', '^', '_', ' ', ' ',
		'^', '^', '^', '^', '^', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_', ' ', ' ',
		'_', '_', '_', '_', '_', '_' };

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Lookup table for Amino Acids to AA conversion, used in aminoAcidsToAminoAcid, can be modified via file
	int8_t kASA::_sAminoAcids_aas[900] = {
		'@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', '@', 
		'G', '\\', '[', 'P', 'I', 'L', '[', ']', 'B', 'D', 'M', 'X', 'T', 'X', 'Z', 'Z', 'W', 'U', 'C', 'Y', ']', 'U', 'D', 'W', 'X', 'J', 'S', 'S', 'W', '^', '@', '@', 'K', 
		'O', 'C', 'C', 'B', 'G', 'G', 'O', 'B', 'V', 'L', 'I', 'K', 'I', 'B', 'N', 'W', 'F', 'X', 'A', 'Q', 'D', '\\', 'S', 'Q', ']', 'A', 'C', 'U', '^', '@', '@', 'B', 'L', 
		'S', 'Z', '[', 'L', 'G', 'I', 'U', '[', 'H', 'D', 'W', '\\', ']', 'U', 'D', 'P', 'A', 'J', ']', 'H', 'S', 'G', 'V', 'V', 'C', 'R', 'Z', '^', '@', '@', 'V', 'I', '\\', 
		'Z', 'M', 'Q', 'Y', 'I', 'S', 'G', 'I', 'J', '[', 'F', 'Y', '[', 'J', 'J', 'C', 'Y', 'U', 'J', 'A', 'F', 'D', 'K', 'L', 'B', 'B', '^', '@', '@', 'F', 'G', 'R', 'Z', 'A', 
		'C', 'O', 'J', 'V', 'N', 'H', 'P', 'N', 'X', 'N', 'C', 'L', 'Q', 'K', 'V', 'X', 'K', 'B', 'O', 'N', 'W', 'L', 'S', 'D', '^', '@', '@', 'Y', 'J', 'K', 'I', 'Q', 'X', 'I', 
		'J', 'M', 'G', '\\', '[', 'M', 'V', 'W', 'M', 'A', 'P', 'F', 'V', 'A', 'G', 'Z', 'B', 'Z', 'D', 'S', '\\', 'M', '^', '@', '@', 'V', 'Z', 'M', 'J', '\\', 'X', 'F', 'T', 
		'V', 'E', 'W', 'C', 'U', 'R', '[', 'Z', 'U', 'H', 'S', 'I', 'W', 'F', 'C', 'N', '\\', 'N', 'V', 'W', 'F', '^', '@', '@', 'X', 'W', 'B', 'B', 'R', 'U', 'V', 'O', 'U', '\\',
		'R', 'Y', 'S', 'Z', 'Q', 'C', 'G', 'L', 'M', 'W', 'Y', 'P', 'Z', 'F', 'G', 'U', 'D', 'S', 'V', '^', '@', '@', 'V', 'A', 'U', 'S', 'R', 'L', 'B', 'G', 'N', 'I', 'F', '\\', 
		'F', 'P', 'M', 'K', 'C', 'F', 'B', 'X', 'U', 'Y', 'D', 'K', 'V', 'W', 'O', 'N', 'N', '^', '@', '@', 'Z', 'U', 'S', 'O', 'I', 'Z', 'J', 'Q', 'J', 'O', 'Z', 'X', 'A', 'X',
		'R', 'C', 'G', '[', '[', 'H', 'P', 'Z', 'N', 'Z', 'D', 'H', 'J', '\\', 'T', '^', '@', '@', 'S', 'W', 'G', 'Z', 'A', 'X', 'H', 'D', 'H', 'Y', 'D', 'Z', 'E', 'K', 'H', 'H',
		'Q', 'H', '\\', 'L', 'O', 'Y', 'S', 'V', 'I', 'X', 'G', ']', 'R', '^', '@', '@', 'Y', 'Z', 'H', 'T', '\\', 'C', '[', 'L', 'D', 'I', 'U', 'G', 'S', '\\', 'V', 'I', 'S', '[',
		'I', 'X', 'E', 'G', '\\', 'A', 'D', 'X', 'R', 'I', 'Y', '^', '@', '@', 'A', 'I', 'O', 'W', 'P', 'A', 'R', 'U', 'I', 'H', 'H', 'S', 'V', ']', 'D', '\\', 'U', 'U', 'T', 'K', 
		'M', 'N', 'J', 'T', 'J', '[', 'A', 'W', 'I', '^', '@', '@', 'P', 'M', 'G', 'Z', 'N', 'X', 'F', '[', 'Q', 'D', 'Y', 'Y', 'N', 'K', 'R', 'H', 'Q', 'O', 'T', 'C', 'Z', 'M', 'Z',
		'I', 'Z', 'X', 'W', 'D', '[', '^', '@', '@', 'A', 'Q', 'X', 'P', 'I', 'F', 'T', 'H', 'H', 'Q', 'V', '[', 'P', 'M', 'U', 'X', 'K', ']', 'E', 'U', 'E', 'R', 'O', 'K', 'J', '\\', 
		'I', 'A', 'E', '^', '@', '@', 'Z', 'S', 'G', 'A', 'L', 'X', 'L', 'I', 'Q', 'O', 'H', '\\', 'H', 'G', 'F', 'B', ']', 'U', 'H', 'J', 'Z', 'J', 'O', 'F', 'Q', ']', 'A', 'H', 'E', 
		'^', '@', '@', 'B', 'J', 'W', 'P', 'N', 'E', 'U', 'V', 'I', ']', 'C', 'N', 'E', 'Y', 'I', 'J', 'O', 'E', 'W', 'R', 'Y', 'G', 'K', 'F', 'C', 'K', 'A', 'Y', 'Q', '^', '@', '@', 'G',
		'\\', 'M', 'G', 'N', 'K', 'Z', 'F', 'I', 'J', 'N', 'G', 'E', 'Y', 'P', 'Z', 'U', 'I', 'C', 'N', 'Q', 'Q', 'R', 'K', 'W', 'U', 'R', 'X', 'T', '^', '@', '@', 'V', 'W', 'G', 'I', 'W',
		'B', 'S', 'R', 'H', 'R', 'J', 'K', 'T', 'X', 'N', 'J', 'X', 'U', 'F', 'F', ']', 'R', 'J', 'C', 'Z', 'G', 'F', ']', 'G', '^', '@', '@', 'G', 'F', 'X', '[', 'H', 'Y', 'S', 'T', '\\',
		'Q', 'F', 'W', 'B', 'J', 'S', 'H', 'W', 'U', ']', 'S', 'K', 'C', 'U', 'A', 'N', 'A', 'U', 'V', 'J', '^', '@', '@', 'T', 'T', 'F', 'M', 'X', 'F', 'A', 'Q', 'Y', 'G', 'N', 'L', 'A',
		'\\', 'M', 'E', ']', 'N', 'B', 'A', 'Q', 'Y', 'T', 'E', 'O', 'X', 'V', 'C', 'J', '^', '@', '@', 'E', 'Q', 'O', ']', 'H', 'N', 'S', '\\', 'P', 'Y', 'J', 'Q', 'D', 'A', 'L', 'E', 
		'V', 'S', 'R', 'M', 'N', 'U', 'Q', 'A', 'B', 'P', 'T', 'P', 'F', '^', '@', '@', 'R', '[', 'D', '[', 'Y', 'M', 'C', 'Q', '\\', 'L', 'Q', '[', 'T', 'N', 'H', 'B', 'N', 'B', 'M', 
		'L', 'P', 'E', 'Y', 'X', 'J', 'W', 'C', 'E', 'C', '^', '@', '@', 'N', '[', 'V', '[', 'X', 'N', 'R', 'B', 'P', 'V', 'H', 'W', 'O', 'Y', 'T', 'A', 'P', 'M', 'F', 'K', 'A', 'A', 
		'E', 'S', 'D', ']', 'S', 'E', 'H', '^', '@', '@', 'Y', 'O', 'Q', 'R', 'V', 'M', 'O', 'L', 'Q', 'K', 'P', 'C', 'M', 'Y', '[', 'M', 'L', 'S', 'H', 'O', 'M', '\\', 'E', 'E', 'V',
		'K', '[', 'L', 'O', '^', '@', '@', 'T', 'Q', 'T', 'T', '[', 'Y', 'O', 'Q', '[', 'Y', 'F', 'V', 'W', 'S', 'W', 'O', 'K', 'P', 'R', 'P', 'D', '\\', 'T', 'K', 'T', ']', 'M', 'T', 'K',
		'^', '@', '@', 'W', 'K', ']', '\\', 'B', 'E', 'O', 'R', 'M', ']', 'K', 'P', '[', 'F', 'L', 'L', 'L', 'L', 'O', 'E', 'D', 'B', 'E', 'R', 'D', 'K', 'P', '\\', 'B', '^', '@', '@', 'B', 'M', 'R'
	};


}//namespace