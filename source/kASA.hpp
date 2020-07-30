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
	protected:
		string _sTemporaryPath = "";

		int32_t _iNumOfThreads;
		int32_t _iNumOfCall;

		uint32_t _iMaxMemUsePerThread;

		int32_t _iMaxK, _iMinK, _iNumOfK, _iHighestK = 12, _iLowestK = 1;
		unique_ptr<int32_t[]> _aOfK;
		const string _sMaxKBlank;

		static int8_t _sAminoAcids_bs[], _sAminoAcids_aas[];
		const static int8_t _sAminoAcids_cs[], _sAminoAcids_un[];
		const int8_t _aRevComp[6] = { 'T','G','A','C','X','Z' };

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

	public:
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
#if (_WIN32 || _WIN64) && !defined(__llvm__)
			string IOCall = (mode != "") ? mode : "wincall autogrow";
#endif
#if __GNUC__ || defined(__llvm__)
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
		const bool _bVerbose;
		const bool _bSixFrames;

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
		kASA(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const string& stxxl_mode = "", const bool& bSixFrames = false) : _sTemporaryPath(tmpPath), _iNumOfThreads(iNumOfProcs), _iNumOfCall(iNumOfCall), _iMaxK(iHigherK), _iMinK(iLowerK), _iNumOfK(_iMaxK - _iMinK + 1), _sMaxKBlank(12, ' '), _bVerbose(bVerbose), _bSixFrames(bSixFrames) {
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

			_iMaxK = (iHigherK <= _iHighestK && iHigherK >= iLowerK) ? iHigherK : _iHighestK;
			_iMinK = (iLowerK < _iLowestK) ? _iLowestK : _iMinK;
			//_iNumOfK = _iMaxK - _iMinK + 1;

			_aOfK.reset(new int32_t[_iNumOfK]);
			for (int32_t i = 0; i < _iNumOfK; ++i) {
				_aOfK[i] = _iMaxK - i;
			}

			createConfig("stxxl_temp_", iNumOfCall, stxxl_mode);
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Converts one aminoacid string to a coded kMer, using 5 bits per amino acid letter (max k = 12 => 60 bits used)
		static inline uint64_t aminoacidTokMer(const string& aminoacid) {
			uint64_t ikMer = 0;

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

		inline uint64_t aminoacidTokMer(const string::const_iterator& begin, const string::const_iterator& end) {
			uint64_t ikMer = 0;

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
				}
				if (lookup != "") {
					if (kMerToAminoacid(elem->first, 12) == lookup) {
						cout << elem->first << " " << kMerToAminoacid(elem->first, 12) << " " << elem->second << endl;
						lookup = "";
					}
				}
				else {
					cout << elem->first << " " << kMerToAminoacid(elem->first, 12) << " " << elem->second << endl;
					++iCounter;
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Counts the number of kMers for each taxon (frequency). 
		// Usually this is counted while creating the index. If for some reason the frequency file got lost, it can be recreated with this.
		inline void GetFrequencyK(const string& contentFile, const string& sLibFile, const string& sOutFile) {
			// test if files exists
			if (!ifstream(contentFile) || !ifstream(sLibFile)) {
				throw runtime_error("OUT: One of the files does not exist");
			}

			ifstream fContent(contentFile);
			uint32_t iIdxCounter = 1;
			unordered_map<uint32_t, uint32_t> mContent; mContent[0] = 0;
			unordered_map<uint32_t, string> mIdxToName; mIdxToName[0] = "non_unique";
			while (fContent.good()) {
				string dummy = "";
				getline(fContent, dummy);
				if (dummy != "") {
					const auto& line = Utilities::split(dummy, '\t');
					mContent[stoul(line[1])] = iIdxCounter;
					mIdxToName[iIdxCounter] = line[0];
					++iIdxCounter;
				}
			}


			ifstream fInfo(sLibFile + "_info.txt");
			uint64_t iSizeOfVec = 0;
			fInfo >> iSizeOfVec;
			stxxlFile libfile(sLibFile, stxxl::file::RDONLY);
			unique_ptr<unique_ptr<const contentVecType_32p>[]> libvec(new unique_ptr<const contentVecType_32p>[_iNumOfThreads]);
			for (int32_t i = 0; i < _iNumOfThreads; ++i) {
				libvec[i].reset(new const contentVecType_32p(&libfile, iSizeOfVec));
			}
			const int32_t& iMaxNumK = _iHighestK - _iLowestK + 1;
			unique_ptr<uint64_t[]> aFrequencyArray(new uint64_t[_iNumOfThreads * iIdxCounter * iMaxNumK]);
			for (uint64_t j = 0; j < _iNumOfThreads * iIdxCounter * iMaxNumK; ++j) {
				aFrequencyArray[j] = 0;
			}

			// count up
			vector<thread> workerThreads;
			uint64_t iPartialSize = iSizeOfVec / _iNumOfThreads;
			uint64_t iPartialRemainder = iSizeOfVec % _iNumOfThreads;
			auto func = [&](const int32_t iThreadID) {
				for (uint64_t i = iThreadID * iPartialSize; i < (iThreadID + 1) * iPartialSize + (((iThreadID + 1) * iPartialSize + iPartialRemainder == iSizeOfVec) ? iPartialRemainder : 0); ++i) {
					const auto& entry = (libvec[iThreadID])->at(i);
					const auto& idx = (iThreadID * iIdxCounter + Utilities::checkIfInMap(mContent, entry.second)->second) * uint64_t(iMaxNumK);
					for (int32_t k = 0; k < iMaxNumK; ++k) {
						const int32_t& shift = 5 * k;
						if (((entry.first >> shift) & 31) != 30) {
							++(aFrequencyArray[idx + k]);
						}
					}
				}
			};

			for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
				workerThreads.push_back(thread(func, iThreadID));
			}
			for (int32_t iThreadID = 0; iThreadID < _iNumOfThreads; ++iThreadID) {
				if (workerThreads[iThreadID].joinable()) {
					workerThreads[iThreadID].join();
				}
			}

			// Gather everything
			for (int32_t i = 1; i < _iNumOfThreads; ++i) {
				for (uint32_t j = 0; j < iIdxCounter; ++j) {
					for (int32_t k = 0; k < iMaxNumK; ++k) {
						aFrequencyArray[j*iMaxNumK + k] += aFrequencyArray[(i*iIdxCounter + j)*iMaxNumK + k];
					}
				}
			}

			// Write to file
			ofstream outFile(sOutFile);
			if (!outFile) {
				throw runtime_error("Frequency file couldn't be opened for writing!");
			}
			for (uint32_t j = 0; j < iIdxCounter; ++j) {
				outFile << Utilities::checkIfInMap(mIdxToName, j)->second << "\t";
				outFile << aFrequencyArray[j*iMaxNumK];
				for (int32_t k = 1; k < iMaxNumK; ++k) {
					outFile << "\t" << aFrequencyArray[j*iMaxNumK + k];
				}
				outFile << endl;
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
					cerr << "ERROR: codon table not found in file. Using built-in." << endl;
				}
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// C++ version of the python script that creates the content file
		inline void generateContentFile(const string& sTaxonomyPath, const string& sAccToTaxFiles, const string& sInput, const string& sOutput, string sTaxonomicLevel, pair<uint32_t,uint32_t> poolAndNames = make_pair(numeric_limits<uint32_t>::max() - 1,0UL)) {
			try {
				if (sTaxonomicLevel != "lowest") {
					if (!ifstream(sTaxonomyPath + "names.dmp") || !ifstream(sTaxonomyPath + "nodes.dmp")) {
						throw runtime_error("The taxonomy files couldn't be found");
					}
				}

				
				Utilities::createFile(sOutput);
				

				auto files = Utilities::gatherFilesFromPath(sInput).first;

				if (_bVerbose) {
					cout << "OUT: Going through fasta(s), gathering accession number(s)..." << endl;
				}

				////////////////////
				unordered_map<string, bool> vAccessions;
				unordered_map<string, uint32_t> vEntriesWithoutAccNr;
				unordered_map<string, string> vNamesFromFasta;
				for (const auto& file : files) {
					ifstream fastaFile;
					igzstream fastaFileGZ;

					if (file.second) {
						fastaFileGZ.open(file.first.c_str());
					}
					else {
						fastaFile.open(file.first);
					}

					string sDummy = "";
					while ((file.second) ? getline(fastaFileGZ, sDummy) : getline(fastaFile, sDummy)) {
						if (sDummy != "") {
							if (sDummy.front() == '>') {
								sDummy.erase(sDummy.begin());
								const auto& sNumbers = Utilities::split(Utilities::split(sDummy, ' ').at(0), '|');
								string acc = "";
								for (const auto& entry : sNumbers) {
									if (entry.find('.') != string::npos) {
										acc = entry;
										break;
									}
								}
								if (acc != "") {
									vAccessions.insert(make_pair(acc, false));
									if (sTaxonomicLevel == "lowest") {
										vNamesFromFasta.insert(make_pair(acc,Utilities::replaceCharacter(sDummy, ',', ' ')));
									}
								}
								else {
									vEntriesWithoutAccNr.insert(make_pair(sDummy, 0));
								}
							}
						}
					}
				}

				if (_bVerbose) {
					cout << "OUT: " << vAccessions.size() << " legit accession number(s) found in " << files.size() << " file(s)." << endl
						<< "OUT: " << vEntriesWithoutAccNr.size() << " entr(y/ies) will get a dummy ID." << endl
						<< "OUT: Translating to taxid(s)..." << endl;
				}

				files.clear();

				size_t iIdentifiedCounter = 0;
				unordered_map<string, unordered_set<string>> taxWithAccNrs;
				unordered_map<string, string> taxToNames;
				bool bNotAllFound = true;

				if (sTaxonomicLevel == "lowest") {
					for (auto& entry : vAccessions) {
						auto elem = taxWithAccNrs.insert(make_pair(to_string(iIdentifiedCounter), unordered_set<string>()));
						if (elem.second) {
							elem.first->second.insert(entry.first);
							entry.second = true;

							taxToNames.insert(make_pair(to_string(iIdentifiedCounter), vNamesFromFasta.at(entry.first)));

							++iIdentifiedCounter;
						}
						else {
							throw runtime_error("AccTaxPair could not be added");
						}
					}
					bNotAllFound = false;
				}
				else {
					files = Utilities::gatherFilesFromPath(sAccToTaxFiles).first;

					////////////////////
					
					for (const auto& file : files) {
						bool isGzipped = file.second;
						ifstream acc2Tax;
						igzstream acc2TaxGZ;
						if (isGzipped) {
							acc2TaxGZ.open(file.first.c_str());
						}
						else {
							acc2Tax.open(file.first);
						}

						string sDummy = "";
						while (((isGzipped) ? getline(acc2TaxGZ, sDummy) : getline(acc2Tax, sDummy)) && bNotAllFound) {
							const auto& columns = Utilities::split(sDummy, '\t');
							auto res1 = vAccessions.find(columns[1]);
							if (res1 != vAccessions.end()) {
								auto res2 = taxWithAccNrs.find(columns[2]);
								if (res2 != taxWithAccNrs.end()) {
									res2->second.insert(res1->first);
								}
								else {
									auto res3 = taxWithAccNrs.insert(make_pair(columns[2], unordered_set<string>()));
									if (res3.second) {
										res3.first->second.insert(res1->first);
									}
									else {
										throw runtime_error("AccTaxPair could not be added");
									}
								}
								res1->second = true;
								++iIdentifiedCounter;
								if (iIdentifiedCounter == vAccessions.size()) {
									bNotAllFound = false;
								}
							}
						}
					}
				}

				////////////////////
				if (_bVerbose && bNotAllFound) {
					cout << "OUT: The following accession numbers have no taxid:" << endl;
				}

				for (const auto& entry : vAccessions) {
					if (entry.second == false) {
						vEntriesWithoutAccNr.insert(make_pair(entry.first, 0));
						if (_bVerbose) {
							cout << entry.first << endl;
						}
					}
				}
				vAccessions.clear();

				////////////////////
				if (_bVerbose) {
					cout << "OUT: Providing dummy taxid(s) if needed..." << endl;
				}

				uint32_t pool = poolAndNames.first;
				if (vEntriesWithoutAccNr.size() > pool) {
					throw runtime_error("Too many dummys, I can't handle this!");
				}
				for (auto& entry : vEntriesWithoutAccNr) {
					entry.second = pool;
					--pool;
				}

				////////////////////
				if (_bVerbose) {
					cout << "OUT: Fetching name(s)..." << endl;
				}

				//taxToNames
				if (taxToNames.empty())
				{
					ifstream names(sTaxonomyPath + "names.dmp");
					string sDummy = "";
					while (getline(names, sDummy)) {
						const auto& line = Utilities::split(sDummy, '|');
						if (line[3] == "\tscientific name\t") {
							taxToNames.insert(make_pair(Utilities::rstrip(line[0]), Utilities::lstrip(Utilities::rstrip(line[1]))));
						}
					}
				}

				////////////////////
				if (_bVerbose) {
					cout << "OUT: Creating dictionary to map taxid(s) to your specified tax. level..." << endl;
				}

				unordered_map<string, pair<string, string>> taxToNodes;
				{
					ifstream nodes(sTaxonomyPath + "nodes.dmp");
					string sDummy = "";
					while (getline(nodes, sDummy)) {
						const auto& line = Utilities::split(sDummy, '|');
						taxToNodes.insert(make_pair(Utilities::rstrip(line[0]), make_pair(Utilities::lstrip(Utilities::rstrip(line[1])), Utilities::lstrip(Utilities::rstrip(line[2])))));
					}
				}

				////////////////////
				if (_bVerbose) {
					cout << "OUT: Linking to your specified tax. level..." << endl;
				}

				unordered_map<string, pair<unordered_set<string>, unordered_set<string>>> taxToTaxWAccs;
				{
					std::transform(sTaxonomicLevel.begin(), sTaxonomicLevel.end(), sTaxonomicLevel.begin(), ::tolower);
					if (!(sTaxonomicLevel == "lowest" ||
						sTaxonomicLevel == "subspecies" ||
						sTaxonomicLevel == "species" ||
						sTaxonomicLevel == "genus" ||
						sTaxonomicLevel == "family" ||
						sTaxonomicLevel == "order" ||
						sTaxonomicLevel == "class" ||
						sTaxonomicLevel == "phylum" ||
						sTaxonomicLevel == "kingdom" ||
						sTaxonomicLevel == "superkingdom" ||
						sTaxonomicLevel == "domain"
						)) {
						cerr << "ERROR: No known tax. level specified. I'll just go with species..." << endl;
						sTaxonomicLevel = "species";
					}

					string upperTax = "";
					pair<string, string> entry;
					for (const auto& elem : taxWithAccNrs) {
						upperTax = elem.first;
						if (sTaxonomicLevel != "lowest") {
							auto res1 = taxToNodes.find(upperTax);
							if (res1 != taxToNodes.end()) {
								entry = res1->second;
							}
							else {
								entry = make_pair("1", "");
							}
							while (entry.second != sTaxonomicLevel && entry.first != "1") {
								upperTax = entry.first;
								entry = taxToNodes.at(upperTax);
							}
							if (entry.first == "1") {
								upperTax = elem.first;
							}


						}
						auto res2 = taxToTaxWAccs.find(upperTax);
						if (res2 != taxToTaxWAccs.end()) {
							res2->second.first.insert(elem.first);
							res2->second.second.insert(elem.second.begin(), elem.second.end());
						}
						else {
							auto it = taxToTaxWAccs.insert(make_pair(upperTax, make_pair(unordered_set<string>(), unordered_set<string>())));
							it.first->second.first.insert(elem.first);
							it.first->second.second.insert(elem.second.begin(), elem.second.end());
						}
					}
					taxWithAccNrs.clear();
				}

				////////////////////

				if (_bVerbose) {
					cout << "OUT: Creating content file..." << endl;
				}

				ofstream contentFile(sOutput);
				if (!contentFile) {
					throw runtime_error("Content file couldn't be opened for writing!");
				}
				uint32_t iUnnamedCounter = 0;
				for (const auto& elem : taxToTaxWAccs) {
					string taxa = "", accnrs = "";
					for (const auto& tax : elem.second.first) {
						taxa += tax + ";";
					}
					if (taxa != "") {
						taxa.erase(taxa.end() - 1);
					}
					for (const auto& accs : elem.second.second) {
						accnrs += accs + ";";
					}
					if (accnrs != "") {
						accnrs.erase(accnrs.end() - 1);
					}

					const auto res = taxToNames.find(elem.first);
					if (res != taxToNames.end()) {
						contentFile << Utilities::removeCharFromString(res->second, ',') << "\t" << elem.first << "\t" << taxa << "\t" << accnrs << endl;
					}
					else {
						contentFile << "unnamed_" << iUnnamedCounter++ << "\t" << elem.first << "\t" << taxa << "\t" << accnrs << endl;
					}
				}

				iUnnamedCounter = poolAndNames.second;
				for (const auto& elem : vEntriesWithoutAccNr) {
					contentFile << "EWAN_" << iUnnamedCounter++ << "\t" << elem.second << "\t" << elem.second << "\t" << elem.first << endl;
				}

			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Add new stuff to the existing content file
		inline void addToContentFile(const string& sTaxonomyPath, const string& sAccToTaxFiles, const string& sInput, const string& sTaxonomicLevel, const string& contentFile) {
			try {
				// Warning: taxonomic levels must not change in an update
				if (!ifstream(sTaxonomyPath + "names.dmp") || !ifstream(sTaxonomyPath + "nodes.dmp")) {
					throw runtime_error("The taxonomy files couldn't be found");
				}

				// TODO: Deal with dummys
				auto func = [](const string& sTempLineSpecIDs, const string& sCurrLineSpecIDs, const string& sTempLineAccNrs, const string& sCurrLineAccNrs) {
					unordered_set<string> sSpecIDs, sAccNrs;
					// concatenate non-redundant species IDs
					const auto& newSpecIDs = Utilities::split(sTempLineSpecIDs, ';'); // tempLineContent[2]
					const auto& currentSpecIDs = Utilities::split(sCurrLineSpecIDs, ';'); // get<1>(entry->second)
					sSpecIDs.insert(newSpecIDs.cbegin(), newSpecIDs.cend());
					sSpecIDs.insert(currentSpecIDs.cbegin(), currentSpecIDs.cend());
					string sNewSpecIDs = "";
					for (const auto& elem : sSpecIDs) {
						sNewSpecIDs += elem + ";";
					}
					sNewSpecIDs.pop_back();

					// concatenate non-redundant accession numbers
					const auto& newAccNrs = Utilities::split(sTempLineAccNrs, ';'); // tempLineContent[3]
					const auto& currentAccNrs = Utilities::split(sCurrLineAccNrs, ';'); // get<2>(entry->second)
					sAccNrs.insert(newAccNrs.cbegin(), newAccNrs.cend());
					sAccNrs.insert(currentAccNrs.cbegin(), currentAccNrs.cend());
					string sNewAccNrs = "";
					for (const auto& elem : sAccNrs) {
						sNewAccNrs += elem + ";";
					}
					sNewAccNrs.pop_back();

					return make_pair(sNewSpecIDs, sNewAccNrs);
				};


				if (_bVerbose) {
					cout << "OUT: reading and checking existing content file..." << endl;
				}

				pair<uint32_t, uint32_t> pool(0, 0);

				ifstream fContent;
				//fContent.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
				fContent.open(contentFile);
				//uint32_t iAmountOfSpecies = 1;
				string sTempLine = "";
				unordered_map<string, tuple<string, string, string>> mOrganisms;
				vector<string> vDuplicateDummys;
				while (getline(fContent, sTempLine)) {
					if (sTempLine != "") {
						const auto& tempLineContent = Utilities::split(sTempLine, '\t');
						bool bDummy = false;
						// check for dummys to get the current counter
						if (tempLineContent[0].find("EWAN") != string::npos) {
							bDummy = true;
							const auto& ewanNumber = stoul(Utilities::split(tempLineContent[0], '_')[1]);
							pool.second = (ewanNumber > pool.second) ? ewanNumber : pool.second; // get highest dummy name counter
							const auto& ewanTaxID = stoul(tempLineContent[1]);
							pool.first = (ewanTaxID < pool.first) ? ewanTaxID : pool.first; // get lowest dummy taxID
						}
						auto entry = mOrganisms.find(tempLineContent[1]);
						if (entry != mOrganisms.end()) {
							if (bDummy) {
								vDuplicateDummys.push_back(tempLineContent[3]);
							}
							else {
								// content file is corrupted
								if (_bVerbose) {
									cout << "OUT: Content file is corrupted, duplicate entries" << tempLineContent[0] << " and " << get<0>(entry->second) << " were found. Merging them now..." << endl;
								}

								const auto& res = func(tempLineContent[2], get<1>(entry->second), tempLineContent[3], get<2>(entry->second));

								entry->second = make_tuple(get<0>(entry->second), res.first, res.second);
							}
						}
						else {
							mOrganisms[tempLineContent[1]] = make_tuple(tempLineContent[0], tempLineContent[2], tempLineContent[3]); // [taxID] -> (Name, species ID, Acc Nrs.) 
						}
					}
				}
				fContent.close();

				for (const auto& entry : vDuplicateDummys) {
					const auto& IDStr = to_string(--pool.first);
					const auto& nameStr = "EWAN_" + to_string(++pool.second);
					mOrganisms[IDStr] = make_tuple(nameStr, IDStr, entry);
				}
				if (pool.first == 0) {
					pool.first = numeric_limits<uint32_t>::max() - 1;
				}

				if (_bVerbose) {
					cout << "OUT: generating temporary content file..." << endl;
				}

				generateContentFile(sTaxonomyPath, sAccToTaxFiles, sInput, _sTemporaryPath + "tempContent.txt", sTaxonomicLevel, pool);

				if (_bVerbose) {
					cout << "OUT: merging content files..." << endl;
				}

				ifstream tempCFile;
				//tempCFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
				tempCFile.open(_sTemporaryPath + "tempContent.txt");
				// same as above
				while (getline(tempCFile, sTempLine)) {
					if (sTempLine != "") {
						const auto& tempLineContent = Utilities::split(sTempLine, '\t');
						auto entry = mOrganisms.find(tempLineContent[1]);
						if (entry != mOrganisms.end()) {
							const auto& res = func(tempLineContent[2], get<1>(entry->second), tempLineContent[3], get<2>(entry->second));

							entry->second = make_tuple(tempLineContent[0], res.first, res.second); // overwrite old name with new one
						}
						else {
							mOrganisms[tempLineContent[1]] = make_tuple(tempLineContent[0], tempLineContent[2], tempLineContent[3]); // [taxID] -> (Name, species ID, Acc Nrs.) 
						}
					}
				}
				tempCFile.close();

				if (_bVerbose) {
					cout << "OUT: writing new content file..." << endl;
				}

				// write content file
				ofstream fContentOS(contentFile);
				if (!fContentOS) {
					throw runtime_error("Content file couldn't be opened for writing!");
				}
				for (const auto& entry : mOrganisms) {
					fContentOS << get<0>(entry.second) << "\t" << entry.first << "\t" << get<1>(entry.second) << "\t" << get<2>(entry.second) << endl;
				}

				remove((_sTemporaryPath + "tempContent.txt").c_str());
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