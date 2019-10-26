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

#include "MetaHeader.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace TrieSnG { // Statics and globals used inside this header
	const uint8_t _iNumOfAminoAcids = 31; // convenience, better than looking through the list everytime

	///DEBUG
	//const uint8_t LISTOFAA[31] = { '@','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','[','\\', ']', '^' };

	///

	//static uint64_t iNumOfNodes[5] = { 0,0,0,0,0 };
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Leaf5 {
	// Sizeof(Leaf): 372 Byte
	uint64_t _vRanges1[TrieSnG::_iNumOfAminoAcids];
	uint32_t _vRanges2[TrieSnG::_iNumOfAminoAcids];


	Leaf5() {
		for (uint8_t i = 0; i < TrieSnG::_iNumOfAminoAcids; ++i) {
			_vRanges1[i] = numeric_limits<uint64_t>::max();
			_vRanges2[i] = 0ul;
		}
	}

	inline uint64_t GetFirst() {
		for (uint8_t i = 0; i < TrieSnG::_iNumOfAminoAcids; ++i) {
			if (_vRanges1[i] != numeric_limits<uint64_t>::max()) {
				return _vRanges1[i];
			}
		}
		return numeric_limits<uint64_t>::max();
	}

	inline uint64_t GetLast() {
		for (int8_t i = static_cast<int8_t>(TrieSnG::_iNumOfAminoAcids) - 1; i >= 0; --i) {
			if (_vRanges1[i] != numeric_limits<uint64_t>::max()) {
				return _vRanges1[i] + _vRanges2[i];
			}
		}
		return 0;
	}

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Node from a Graph(Trie) with a list of its Edges and current level
class Node {
private:

	// Sizeof(Node): 256 Byte
	void* _Edges[TrieSnG::_iNumOfAminoAcids];

	int8_t _iLevelDiff;

public:

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	inline void IncreaseIndex(const uint32_t& kMer, const uint64_t& iStartIdx, const uint32_t& iEndIdx, uint64_t& iSizeOfTrie) {
		if (_iLevelDiff - 2 > 0) {
			const uint8_t& tempVal = (uint8_t)((kMer >> (_iLevelDiff - 1) * 5) & 31); // read only the relevant 5 bytes (and start with level 1)
			//cout << bitset<15>(tempVal) << endl;
			//cout << TrieSnG::LISTOFAA[tempVal] << " " << kASA::kASA::kMerToAminoacid(kMer, 12) << endl;

			Node* edge = static_cast<Node*>(_Edges[tempVal]);
			if (edge == nullptr) {
				//iNumOfNodes[_iLevel] += 1;
				_Edges[tempVal] = new Node(kMer, _iLevelDiff - 1, iStartIdx, iEndIdx, iSizeOfTrie);
				iSizeOfTrie += 256;
			}
			else {
				edge->IncreaseIndex(kMer, iStartIdx, iEndIdx, iSizeOfTrie);
			}
		}
		else {
			const uint8_t& tempVal = (uint8_t)((kMer >> (_iLevelDiff - 1) * 5) & 31);
			const uint8_t& lastVal = (uint8_t)((kMer >> (_iLevelDiff - 2) * 5) & 31);

			//cout << TrieSnG::LISTOFAA[tempVal] << " " << TrieSnG::LISTOFAA[lastVal]  << " " << kASA::kASA::kMerToAminoacid(kMer,12) << endl;

			Leaf5* edge = static_cast<Leaf5*>(_Edges[tempVal]);
			if (edge == nullptr) {
				_Edges[tempVal] = new Leaf5();
				iSizeOfTrie += sizeof(Leaf5);
				static_cast<Leaf5*>(_Edges[tempVal])->_vRanges1[lastVal] = iStartIdx;
				static_cast<Leaf5*>(_Edges[tempVal])->_vRanges2[lastVal] = iEndIdx;
			}
			else {
				edge->_vRanges1[lastVal] = iStartIdx;
				edge->_vRanges2[lastVal] = iEndIdx;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	explicit Node(const int8_t& iMaxLevel) : _Edges(), _iLevelDiff(iMaxLevel) {
		// for root node
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Node(const uint32_t& kMer, const int8_t& iLevel, const uint64_t& iStartIdx, const uint32_t& iEndIdx, uint64_t& iSizeOfTrie) : _Edges(), _iLevelDiff(iLevel) {
		IncreaseIndex(kMer, iStartIdx, iEndIdx, iSizeOfTrie);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	~Node() {
		if (_iLevelDiff - 2 > 0) {
			for (uint8_t i = 0; i < TrieSnG::_iNumOfAminoAcids; ++i) {
				delete static_cast<Node*>(_Edges[i]);
			}
		}
		else {
			for (uint8_t i = 0; i < TrieSnG::_iNumOfAminoAcids; ++i) {
				delete static_cast<Leaf5*>(_Edges[i]);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	inline std::tuple<uint64_t, uint32_t> GetRange(const uint64_t& kMer, const int8_t& iCurrentK, const int8_t& iMinK, const int8_t& iMaxLevel) {
		/// First handle the usual case
		if (_iLevelDiff != iMaxLevel - iCurrentK) {
			if (_iLevelDiff - 2 > 0) {
				const uint8_t& tempVal = (uint8_t)((kMer >> (_iLevelDiff - 1) * 5) & 31); // read only the relevant 5 bytes (and start with level 1)
				//cout << TrieSnG::LISTOFAA[tempVal] << endl;
				Node* edge = static_cast<Node*>(_Edges[tempVal]);
				if (edge == nullptr) {
					return std::tuple<uint64_t, uint32_t>(numeric_limits<uint64_t>::max(),0);
				}
				else {
					return edge->GetRange(kMer, iCurrentK, iMinK, iMaxLevel);
				}
			}
			else {
				const uint8_t& tempVal = (uint8_t)((kMer >> (_iLevelDiff - 1) * 5) & 31);
				//cout << TrieSnG::LISTOFAA[tempVal] << endl;
				Leaf5* edge = static_cast<Leaf5*>(_Edges[tempVal]);
				if (edge == nullptr) {
					return std::tuple<uint64_t, uint32_t>(numeric_limits<uint64_t>::max(), 0);
				}
				else {
					const uint8_t& lastVal = (uint8_t)((kMer >> (_iLevelDiff - 2) * 5) & 31);
					//cout << TrieSnG::LISTOFAA[lastVal] << endl;
					return make_tuple(edge->_vRanges1[lastVal], edge->_vRanges2[lastVal]);
				}
			}
		}
		else {
			std::tuple<uint64_t, uint32_t> outTuple(numeric_limits<uint64_t>::max(),0);
			get<0>(outTuple) = GetRangeAfterMinK(iMaxLevel, true);
			get<1>(outTuple) = static_cast<uint32_t>(GetRangeAfterMinK(iMaxLevel, false) - get<0>(outTuple));
			return outTuple;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	inline bool isInTrie(const uint32_t& kMer, const int8_t& iMinK, const int8_t& iMaxLevel) const {
		/// First handle the usual case
		if (_iLevelDiff != iMaxLevel - iMinK) {
			if (_iLevelDiff - 2 > 0) {
				const uint8_t& tempVal = (uint8_t)((kMer >> (5 + _iLevelDiff - iMaxLevel) * 5) & 31); // read only the relevant 5 bytes (and start with level 1)
				
				//cout << TrieSnG::LISTOFAA[tempVal] << endl;
				
				Node* edge = static_cast<Node*>(_Edges[tempVal]);
				if (edge == nullptr) {
					return false;
				}
				else {
					return edge->isInTrie(kMer, iMinK, iMaxLevel);
				}
			}
			else {
				const uint8_t& tempVal = (uint8_t)((kMer >> (5 + _iLevelDiff - iMaxLevel) * 5) & 31);

				//cout << TrieSnG::LISTOFAA[tempVal] << endl;

				Leaf5* edge = static_cast<Leaf5*>(_Edges[tempVal]);
				if (edge == nullptr) {
					return false;
				}
				else {
					const uint8_t& lastVal = (uint8_t)((kMer >> (5 + _iLevelDiff - iMaxLevel - 1) * 5) & 31);

					//cout << TrieSnG::LISTOFAA[lastVal] << endl;

					return (edge->_vRanges1[lastVal] != numeric_limits<uint64_t>::max());
				}
			}
		}
		else {
			return true;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	inline uint64_t GetRangeAfterMinK(const int8_t& iMaxLevel, const bool& bStart) {
		if (_iLevelDiff - 2 > 0) {
			if (bStart) {
				for (uint8_t i = 0; i < TrieSnG::_iNumOfAminoAcids; ++i) {
					if (_Edges[i]) {
						return static_cast<Node*>(_Edges[i])->GetRangeAfterMinK(iMaxLevel, true);
					}
				}
				return 0;
			}
			else {
				for (int8_t i = static_cast<int8_t>(TrieSnG::_iNumOfAminoAcids) - 1; i >= 0; --i) {
					if (_Edges[i]) {
						return static_cast<Node*>(_Edges[i])->GetRangeAfterMinK(iMaxLevel, false);
					}
				}
				return numeric_limits<uint64_t>::max();
			}
			
		}
		else {
			if (bStart) {
				for (uint8_t i = 0; i < TrieSnG::_iNumOfAminoAcids; ++i) {
					if (_Edges[i]) {
						return static_cast<Leaf5*>(_Edges[i])->GetFirst();
					}
				}
				return 0;
			}
			else {
				for (int8_t i = static_cast<int8_t>(TrieSnG::_iNumOfAminoAcids) - 1; i >= 0; --i) {
					if (_Edges[i]) {
						return static_cast<Leaf5*>(_Edges[i])->GetLast();
					}
				}
			}
			return numeric_limits<uint64_t>::max();
		}
	}


};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Prefix-tree for a logarithmic access of the first few letters
class Trie {
private:
	
	int8_t _iMaxK, _iMinK, _iMaxLevel, _ikForIsInTrie;
	unique_ptr<uint64_t[]> _iTempKMer;
	uint64_t _Bitmask = 31;
	unique_ptr<std::tuple<uint64_t, uint32_t>[]> _tempTuple;

	uint32_t _iMaxRange = 0;

	// this counts the size of the trie in bytes to substract it later from the available memory in RAM
	uint64_t iSizeOfTrie = 0;

	unique_ptr<Node> _root;

	vector<packedBigPair> _vAllPrefixes;

public:
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Trie(const int8_t& iMaxK, const int8_t& iMinK, const int8_t& iMaxLevel, const int32_t& iThreads = 1, const bool& bShitLoadOfRAM = false) : _iMaxK(iMaxK), _iMinK(iMinK), _iMaxLevel(iMaxLevel), _ikForIsInTrie(6) {
		for (uint8_t i = 1; i < iMaxLevel; ++i) {
			_Bitmask |= 31ULL << (5 * i);
		}
		_Bitmask <<= (iMaxK - iMaxLevel) * 5;
		_iTempKMer.reset(new uint64_t[iThreads]);
		_tempTuple.reset(new std::tuple<uint64_t, uint32_t>[iThreads]);
		for (int32_t i = 0; i < iThreads; ++i) {
			_iTempKMer[i] = 0;
		}

		if (bShitLoadOfRAM) {
			_vAllPrefixes.resize(kASA::kASA::aminoacidTokMer("]^^^^^"), packedBigPair(0,0));
			iSizeOfTrie = _vAllPrefixes.size() * sizeof(packedBigPair);
		}
	}

	inline void SetForIsInTrie(const uint8_t& val) {
		_ikForIsInTrie = val;
	}

	inline uint64_t GetSize() {
		return iSizeOfTrie;
	}

	inline uint32_t GetMaxRange() const {
		return _iMaxRange;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Save to file
	template<typename T>
	inline void Save(const T* vKMerVec, const string& savePath) {
		try {
			fstream file(savePath, ios::out);
			//ostream console(cout.rdbuf());
			//_root->Save(file, _iMaxLevel);
			//uint64_t itmpCount = 0;
			uint64_t iCount = 1;
			uint64_t iKnownShortMer = get<0>(vKMerVec->at(0)) & _Bitmask;
			stxxl::vector_bufreader<typename T::const_iterator> bufferedReader(vKMerVec->cbegin() + 1, vKMerVec->cend(), 1);
			for (; !bufferedReader.empty(); ++bufferedReader) {
				const auto& entry = get<0>(*bufferedReader);
				const uint64_t& iTemp = entry & _Bitmask;
				if (iTemp != iKnownShortMer) {
					file << (iKnownShortMer >> (_iMaxK - _iMaxLevel) * 5) << ' ' << iCount << endl;
					iKnownShortMer = iTemp;
					//cout << itmpCount << " " << itmpCount + iCount << endl;
					//itmpCount += iCount;
					iCount = 1;
				}
				else {
					++iCount;
				}
			}
			file << (iKnownShortMer >> (_iMaxK - _iMaxLevel) * 5) << ' ' << iCount << endl;
			//cout << itmpCount << " " << itmpCount + iCount << endl;
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}
	////////////////
	template<typename T>
	inline void SaveToStxxlVec(const T* vKMerVec, const string& savePath, unique_ptr<uint64_t[]>* arrOfFreqs = nullptr, const unordered_map<uint32_t, uint32_t>& mContent = unordered_map<uint32_t, uint32_t>()) { //passing a non-const reference with default value is a hassle
		try {
			ofstream dummyFile(savePath + "_trie");
			dummyFile.close();
			stxxlFile trieFile(savePath+"_trie", stxxl::file::RDWR);
			trieVector trieVec(&trieFile, 0);
			uint64_t iCount = 1;
			uint64_t iKnownShortMer = vKMerVec->at(0).first & _Bitmask;
			stxxl::vector_bufreader<typename T::const_iterator> bufferedReader(vKMerVec->cbegin() + 1, vKMerVec->cend(), 0);
			for (; !bufferedReader.empty(); ++bufferedReader) {
				const auto& entry = *bufferedReader;
				if (!mContent.empty()) {
					const auto& idx = Utilities::checkIfInMap(mContent, entry.second)->second * 12;
					for (uint8_t k = 0; k < 12; ++k) {
						if (((entry.first >> 5 * k) & 31) != 30) {
							(*arrOfFreqs)[idx + k]++;
						}
					}
				}
				const uint64_t& iTemp = entry.first & _Bitmask;
				if (iTemp != iKnownShortMer) {
					trieVec.push_back(packedBigPairTrie(uint32_t(iKnownShortMer >> (_iMaxK - _iMaxLevel) * 5), iCount));
					iKnownShortMer = iTemp;
					iCount = 1;
				}
				else {
					++iCount;
				}
			}
			trieVec.push_back(packedBigPairTrie(uint32_t(iKnownShortMer >> (_iMaxK - _iMaxLevel) * 5), iCount));
			ofstream sizeFile(savePath + "_trie.txt");
			sizeFile << trieVec.size();
			trieVec.export_files("_");
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}


	///////////////////////////
	inline void LoadFromStxxlVec(const string& loadPath) {
		try {
			_root.reset(new Node(_iMaxLevel));
			ios::sync_with_stdio(false);
			ifstream sizeFile(loadPath + "_trie.txt");
			uint64_t iSizeOfTrieVec;
			sizeFile >> iSizeOfTrieVec;
			stxxlFile trieFile(loadPath + "_trie", stxxl::file::RDONLY);
			const trieVector trieVec(&trieFile, iSizeOfTrieVec);
			uint64_t iStart = 0, iEnd = 0, iCount = 0;
			uint32_t iReducedkMer;
			//ofstream smallAA("D:/tmp/paper/trie.txt");
			//uint64_t iMaxRange = 0, iSumOfRanges = 0, iNumOfRanges = 0;
			//vector<uint32_t> vRanges;
			for (const auto& entry : trieVec) {
				iReducedkMer = entry.first;
				iCount = entry.second;
				//vRanges.push_back(iCount);
				//if (iCount > iMaxRange) {
				//	iMaxRange = iCount;
				//}
				//iSumOfRanges += iCount;
				//++iNumOfRanges;
				//cout << kASA::kASA::kMerToAminoacid(iReducedkMer, 6) << " " << iCount << endl;
				const uint64_t& iSum = iCount + iEnd;
				const int32_t iDiff = static_cast<int32_t>(int64_t(iSum) - 1 - iStart);
				assert(iDiff >= 0);
				if (static_cast<uint32_t>(iDiff) > _iMaxRange) {
					_iMaxRange = iDiff;
				}
				if (_vAllPrefixes.size() > 0) {
					_vAllPrefixes[iReducedkMer] = make_tuple(iStart, iDiff);
				}
				else {
					_root->IncreaseIndex(iReducedkMer, iStart, static_cast<uint32_t>(iDiff), iSizeOfTrie);
				}
				iStart = iSum;
				iEnd += iCount;
			}
			//sort(vRanges.begin(), vRanges.end());
			//cout << vRanges[vRanges.size() / 2] << " " << vRanges.back() << endl;
			//cout << iMaxRange << " " << iSumOfRanges/double(iNumOfRanges) << " " << iNumOfRanges <<  endl;
			//cout << "Nodes per Level: 1:" << iNumOfNodes[0] << " 2:" << iNumOfNodes[1] << " 3:" << iNumOfNodes[2] << " 4:" << iNumOfNodes[3] << " 5:" << iNumOfNodes[4] << endl;
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Create without saving
	template<typename T>
	inline void Create(const T vKMerVec) { //TODO!!!

		_root.reset(new Node(_iMaxLevel));
		uint64_t iCount = 0, iStart = 0; //iEnd = 0;
		uint64_t iKnownShortMer = get<0>(vKMerVec->at(0)) & _Bitmask;
		uint64_t iTemp = 0;
		for (const auto& elem : *vKMerVec) {
			const auto& entry = get<0>(elem);
			//cout << kASA::kASA::kMerToAminoacid(entry, 12) << endl;
			iTemp = entry & _Bitmask;
			if (iTemp != iKnownShortMer) {
				const uint64_t& iSum = iCount + iStart;
				_root->IncreaseIndex(static_cast<uint32_t>(iKnownShortMer >> (_iMaxK - _iMaxLevel) * 5), iStart, static_cast<uint32_t>( iSum - 1), iSizeOfTrie);
				//cout << kASA::kASA::kMerToAminoacid(entry, 12) << " " << kASA::kASA::kMerToAminoacid((iKnownShortMer >> (_iMaxK - _iMaxLevel) * 5), 12) << endl;
				//iEnd += iCount;
				iKnownShortMer = iTemp;
				iStart = iSum;
				iCount = 1;
			}
			else {
				++iCount;
			}
		}
		_root->IncreaseIndex(static_cast<uint32_t>(iTemp >> (_iMaxK - _iMaxLevel) * 5), iStart, static_cast<uint32_t>(iCount + iStart - 1), iSizeOfTrie);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Access
	inline std::tuple<uint64_t, uint32_t> GetIndexRange(const uint64_t& kMer, const int8_t& iCurrentK, const uint64_t& iLatestStartIndex, const int32_t& iShift, const int32_t& iThreadID = 1) const {
		//cout << "tempmr: " << kASA::kASA::kMerToAminoacid(tempMer, 12) << endl;
		const uint64_t& iShiftedkMer = kMer << iShift;
		if (iShiftedkMer != _iTempKMer[iThreadID]) {
			_iTempKMer[iThreadID] = iShiftedkMer;
			_tempTuple[iThreadID] = _root->GetRange(kMer, iCurrentK, _iMinK, _iMaxLevel);
			if (get<0>(_tempTuple[iThreadID]) < iLatestStartIndex && iLatestStartIndex <= get<0>(_tempTuple[iThreadID]) + get<1>(_tempTuple[iThreadID])) {
				get<0>(_tempTuple[iThreadID]) = iLatestStartIndex;
			}
			else {
				if (iLatestStartIndex > get<0>(_tempTuple[iThreadID]) + get<1>(_tempTuple[iThreadID])) {
					_tempTuple[iThreadID] = make_tuple(numeric_limits<uint64_t>::max(), 0);
				}
			}
			return _tempTuple[iThreadID];
		}
		else {
			// This class can't know where the last search stopped so it must be given that information to save redundant calls if necessary
			if (iLatestStartIndex > get<0>(_tempTuple[iThreadID]) + get<1>(_tempTuple[iThreadID])) {
				_tempTuple[iThreadID] = make_tuple(numeric_limits<uint64_t>::max(), 0);
			}
			return std::tuple<uint64_t, uint32_t>(iLatestStartIndex, std::get<1>(_tempTuple[iThreadID]));
		}
	}

	inline std::tuple<uint64_t, uint32_t> GetIndexRange(const uint64_t& kMer, const int8_t& iCurrentK) const {
		if (_vAllPrefixes.size() > 0) {
			return make_tuple(_vAllPrefixes[kMer].first, _vAllPrefixes[kMer].second);
		}
		else {
			return _root->GetRange(kMer, iCurrentK, _iMinK, _iMaxLevel);
		}
	}
	
	///////////////////////////
	inline bool IsInTrie(const uint64_t kMer) const {
		return _root->isInTrie(static_cast<uint32_t>(kMer >> 30), _ikForIsInTrie, _iMaxLevel);
	}
};
