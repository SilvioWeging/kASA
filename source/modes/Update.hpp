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
#include "Build.hpp"

namespace kASA {
	template<class vecType, class elemType, class intType>
	class Update : public Read<vecType, elemType, intType>, Build<vecType, elemType> {

		typedef Read<vecType, elemType, intType> Base;

	public:
		Update(const string& tmpPath, const int32_t& iNumOfProcs, const int32_t& iHigherK, const int32_t& iLowerK, const int32_t& iNumOfCall, const bool& bVerbose = false, const string& stxxl_mode = "", const bool& bSixFrames = false) : Read<vecType, elemType, intType>(tmpPath, iNumOfProcs, iHigherK, iLowerK, iNumOfCall, bVerbose, stxxl_mode, bSixFrames), Build<vecType, elemType>() {}
		Update(const kASA& obj, const bool& bUnfunny = false) : Read<vecType, elemType, intType>(obj, bUnfunny), Build<vecType, elemType>() {}

	public:


		void DeleteFromLib(const string& contentFile, const string& sLibFile, const string& fOutfile, const string& sDelNodesPath, const bool& bOverwrite, int64_t iMemory) {
			try {
				// test if files exists
				if (!ifstream(sLibFile)) {
					throw runtime_error("Library file does not exist");
				}

				if (!ifstream(sDelNodesPath)) {
					throw runtime_error("delnodes.dmp not found");
				}

				ifstream deletedNodes(sDelNodesPath);
				string line = "";
				unordered_set<uint32_t> vTaxIdsTBD;
				while (getline(deletedNodes, line)) {
					if (line != "") {
						vTaxIdsTBD.insert(stoul((Utilities::split(line, '\t'))[0]));
					}
				}

				// Create new vector
				string sTempFile = Base::_sTemporaryPath + "_update_temp_" + to_string(Base::_iNumOfCall);
				{ // need a scope so that the stxxl releases the vector before copying in the case of overwrite == true
					ifstream infoFile(sLibFile + "_info.txt");
					uint64_t iLibSize;
					infoFile >> iLibSize;
					infoFile.close();
					stxxlFile libFile(sLibFile, stxxl::file::RDONLY);
					const vecType vLib(&libFile, iLibSize);
					ofstream dummy;
					//dummy.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
					dummy.open((bOverwrite) ? sTempFile : fOutfile);
					dummy.close();
					stxxlFile outFile((bOverwrite) ? sTempFile : fOutfile, stxxl::file::RDWR);
					vecType vOut(&outFile, iLibSize);

					// All entries with valid taxonomic ids will be written to the other library
					auto itOut = vOut.begin();
					for (auto itLib = vLib.cbegin(); itLib != vLib.cend(); ++itLib) {
						if (vTaxIdsTBD.find(itLib->second) == vTaxIdsTBD.end()) {
							*itOut++ = *itLib;
						}
					}
					vOut.resize(itOut - vOut.begin(), true);
					vOut.export_files("_");

					ofstream infoFileOut(fOutfile + "_info.txt");
					infoFileOut << vOut.size();

					Trie<intType> T(static_cast<int8_t>(((is_same<vecType, contentVecType_128>::value) ? HIGHESTPOSSIBLEK : 12)), static_cast<int8_t>(Base::_iMinK), 6);
					T.SaveToStxxlVec(&vOut, fOutfile);
				}
				if (bOverwrite) {
					remove(sLibFile.c_str());
					Utilities::copyFile(sTempFile, fOutfile);
				}

				if (Base::_bVerbose) {
					cout << "OUT: Creating frequency file... " << endl;
				}

				Base::template GetFrequencyK<vecType>(contentFile, fOutfile, iMemory);
			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Update an existing library with a fasta file
		void UpdateFromFasta(const string& contentFile, const string& sLibFile, const string& sDirectory, string fOutFile, const bool& bOverwrite, const uint64_t& iMemory, const float& fPercentageOfThrowAway, const pair<unordered_map<uint32_t, uint32_t>, unordered_map<uint32_t, uint32_t>>& mapsForDummys) {
			try {

				// test if files exists
				if (!ifstream(contentFile) || !ifstream(sLibFile)) {
					throw runtime_error("One of the files does not exist");
				}

				// Create new index
				string sMergedIndexFileName = Base::_sTemporaryPath + "_tempUpdate_out_" + to_string(Base::_iNumOfCall);
				string sTempIndexFileName = Base::_sTemporaryPath + "_tempUpdate_" + to_string(Base::_iNumOfCall);
				unordered_map<uint32_t, uint32_t> mIDsAsIdx;
				{
					auto buildResults = this->BuildAll(contentFile, sDirectory, sTempIndexFileName, iMemory, fPercentageOfThrowAway, false, true);

					mIDsAsIdx = get<0>(buildResults);
					const uint64_t& iSizeOfTempIndex = get<1>(buildResults);

					// final index file
					Utilities::checkIfFileCanBeCreated(sMergedIndexFileName);
					unique_ptr<stxxlFile> stxxlOutVec(new stxxlFile(sMergedIndexFileName, stxxl::file::RDWR));
					unique_ptr<vecType> vOutVec(new vecType(stxxlOutVec.get(), 0));
					typename vecType::bufwriter_type vNCBuff(*vOutVec);

					{ // new scope for merging so that the temporary can be removed afterwards
						// get Size of Lib and open it for reading
						fstream fLibInfo(sLibFile + "_info.txt", ios::in);
						uint64_t iSizeOfLib = 0;
						fLibInfo >> iSizeOfLib;
						fLibInfo.close();
						unique_ptr<stxxlFile> stxxlLibFile(new stxxlFile(sLibFile, stxxl::file::RDONLY));
						unique_ptr<const vecType> vLibIn(new vecType(stxxlLibFile.get(), iSizeOfLib));
						typename vecType::bufreader_type vCBuff(*vLibIn);

						// just created temporary index from new data
						unique_ptr<stxxlFile> stxxlTempVec(new stxxlFile(sTempIndexFileName, stxxl::file::RDONLY));
						unique_ptr<vecType> vTempVec(new vecType(stxxlTempVec.get(), iSizeOfTempIndex));

						// Merge existing index and new one
						const auto& vOutSize = this->merge(vNCBuff, vCBuff, iSizeOfLib, vTempVec->cbegin(), vTempVec->cend(), mapsForDummys);
						vOutVec->resize(vOutSize, true);
					}
					vNCBuff.finish();

					remove(sTempIndexFileName.c_str());
					remove((sTempIndexFileName + "_info.txt").c_str());

					ofstream fOutInfo(fOutFile + "_info.txt");
					fOutInfo << vOutVec->size();
					if (is_same<vecType, contentVecType_128>::value) {
						fOutInfo << endl << 128;
					}

					if (Base::_bVerbose) {
						cout << "OUT: Creating trie... " << endl;
					}

					vOutVec->export_files("_");

					Trie<intType> T(static_cast<int8_t>(((is_same<vecType, contentVecType_128>::value) ? HIGHESTPOSSIBLEK : 12)), static_cast<int8_t>(Base::_iMinK), 6);
					T.SaveToStxxlVec(vOutVec.get(), fOutFile);
				}

				// "overwrite" old index
				if (bOverwrite) {
					remove(sLibFile.c_str());
				}
				Utilities::moveFile(sMergedIndexFileName, fOutFile);

				if (Base::_bVerbose) {
					cout << "OUT: Creating frequency file... " << endl;
				}

				Base::template GetFrequencyK<vecType>(contentFile, fOutFile, iMemory, mIDsAsIdx);

			}
			catch (...) {
				cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
			}
		}

	};

}