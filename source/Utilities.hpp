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

#include <initializer_list>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <utility>
#include <unordered_map>

namespace Utilities {

	///////////////////////////////////////////////////////
	// check if value is in list of values
	template <typename T>
	inline bool is_in(const T& val, const initializer_list<T>& list)
	{
		for (const auto& i : list) {
			if (val == i) {
				return true;
			}
		}
		return false;
	}

	///////////////////////////////////////////////////////
	struct rangeContainer {
		vector<pair<uint32_t, uint32_t>> kMers_GT6;
		vector<pair<uint64_t, uint32_t>> kMers_ST6;
		uint32_t range;
	};

	///////////////////////////////////////////////////////
	inline string removeSpaceAndEndline(const string& sIn) {
		string sOut = sIn;
		sOut.erase(find_if(sOut.rbegin(), sOut.rend(), [](char ch) {
			return !(ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t');
		}).base(), sOut.end());
		return sOut;
	}

	///////////////////////////////////////////////////////
	inline string removeCharFromString(const string& sIn, const char& c) {
		string sOut = sIn;
		sOut.erase(find_if(sOut.rbegin(), sOut.rend(), [&c](char ch) {
			return !(ch == c);
		}).base(), sOut.end());
		return sOut;
	}

	///////////////////////////////////////////////////////
	inline string lstrip(const string& sIn) {
		string sOut = sIn;
		auto it = sOut.begin();
		while (*it == ' ' || *it == '\n' || *it == '\r' || *it == '\t') {
			it = sOut.erase(it);
		}
		return sOut;
	}

	///////////////////////////////////////////////////////
	inline string rstrip(const string& sIn) {
		string sOut = sIn;
		auto it = sOut.end() - 1;
		while (*it == ' ' || *it == '\n' || *it == '\r' || *it == '\t') {
			it = sOut.erase(it);
			it--;
		}
		return sOut;
	}
	///////////////////////////////////////////////////////
	inline string lstrip(const string& sIn, const char& c) {
		string sOut = sIn;
		auto it = sOut.begin();
		while (*it == c) {
			it = sOut.erase(it);
		}
		return sOut;
	}

	///////////////////////////////////////////////////////
	inline string rstrip(const string& sIn, const char& c) {
		string sOut = sIn;
		auto it = sOut.end() - 1;
		while (*it == c) {
			it = sOut.erase(it);
			it--;
		}
		return sOut;
	}

	///////////////////////////////////////////////////////
	inline pair<vector<string>, size_t> gatherFilesFromPath(const string& sPath) {
		try {
			vector<string> files;
			size_t overallFileSize = 0;
#if _WIN32 || _WIN64 
			if (sPath.back() == '/') {
				for (auto& fsPath : filesystem::directory_iterator(sPath)) {
					const string& fileName = (fsPath.path()).string();
					files.push_back(fileName);

					bool isGzipped = (fileName[fileName.length() - 3] == '.' && fileName[fileName.length() - 2] == 'g' && fileName[fileName.length() - 1] == 'z');
					if (!isGzipped) {
						ifstream fast_q_a_File;
						//fast_q_a_File.exceptions(std::ifstream::failbit | std::ifstream::badbit);
						fast_q_a_File.open(fileName);
						fast_q_a_File.seekg(0, fast_q_a_File.end);
						overallFileSize += fast_q_a_File.tellg();
					}
				}
			}
#else
#if __GNUC__ || defined(__llvm__)
			if (sPath.back() == '/') {
				DIR           *dirp;
				struct dirent *directory;

				dirp = opendir(sPath.c_str());
				if (dirp) {
					while ((directory = readdir(dirp)) != NULL)
					{
						const string fName(directory->d_name);
						if (fName != "." && fName != "..") {
							const string& fileName = sPath + fName;
							files.push_back(fileName);

							bool isGzipped = (fileName[fileName.length() - 3] == '.' && fileName[fileName.length() - 2] == 'g' && fileName[fileName.length() - 1] == 'z');
							if (!isGzipped) {
								ifstream fast_q_a_File;
								//fast_q_a_File.exceptions(std::ifstream::failbit | std::ifstream::badbit);
								fast_q_a_File.open(fileName);
								fast_q_a_File.seekg(0, fast_q_a_File.end);
								overallFileSize += fast_q_a_File.tellg();
		}
						}
					}

					closedir(dirp);
				}
			}
#endif
#endif
			else {
				bool isGzipped = (sPath[sPath.length() - 3] == '.' && sPath[sPath.length() - 2] == 'g' && sPath[sPath.length() - 1] == 'z');
				if (!isGzipped) {
					ifstream fast_q_a_File;
					//fast_q_a_File.exceptions(std::ifstream::failbit | std::ifstream::badbit);
					fast_q_a_File.open(sPath);
					fast_q_a_File.seekg(0, fast_q_a_File.end);
					overallFileSize += fast_q_a_File.tellg();
				}
				files.push_back(sPath);
			}

			return make_pair(files,overallFileSize);
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}


	///////////////////////////////////////////////////////
	template<typename Out>
	inline void split(const string &s, char delim, Out result) {
		stringstream ss;
		ss.str(s);
		string item;
		while (getline(ss, item, delim)) {
			*(result++) = item;
		}
	}

	///////////////////////////////////////////////////////
	inline vector<string> split(const string &s, char delim) {
		vector<string> elems;
		split(s, delim, back_inserter(elems));
		return elems;
	}

	///////////////////////////////////////////////////////
	inline vector<unsigned int> splitToUInts(const string &s, char delim) {
		vector<unsigned int> elems;
		stringstream ss;
		ss.str(s);
		string item;
		while (getline(ss, item, delim)) {
			elems.push_back(stoul(item));
		}
		return elems;
	}

	///////////////////////////////////////////////////////
	template<typename T>
	inline pair<string,bool> getChunk(T&& input, uint64_t&& iNumOfChars = 0) {
		char sTempCArr[100000]; //100000
		input.get(sTempCArr, 100000);
		iNumOfChars = input.gcount();
		bool bLineFinished = false;
		if (iNumOfChars != 99999 && input) {
			bLineFinished = true;
			input.get();
		}
		string sTempString(sTempCArr);
		if (sTempString == "") {
			if (input.fail() && !input.eof()) {
				input.clear();
				getline(input, sTempString);
			}
		}
		return make_pair(sTempString, bLineFinished);
	}

	///////////////////////////////////////////////////////

	inline void createFile(const string& s) {
		ofstream derp(s);
		if (derp.fail()) {
			throw runtime_error("File couldn't be created, maybe a wrong path was used?");
		}
	}

	///////////////////////////////////////////////////////
	inline uint64_t countBits(uint64_t val) {
		// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
		if (val == 0) {
			return 0;
		}
		val = val - ((val >> 1) & (uint64_t)~(uint64_t)0 / 3);
		val = (val & (uint64_t)~(uint64_t)0 / 15 * 3) + ((val >> 2) & (uint64_t)~(uint64_t)0 / 15 * 3);
		val = (val + (val >> 4)) & (uint64_t)~(uint64_t)0 / 255 * 15;
		return ((uint64_t)(val * ((uint64_t)~(uint64_t)0 / 255)) >> (sizeof(uint64_t) - 1) * CHAR_BIT);
	}

	///////////////////////////////////////////////////////
	const int tab64[64] = {
		63,  0, 58,  1, 59, 47, 53,  2,
		60, 39, 48, 27, 54, 33, 42,  3,
		61, 51, 37, 40, 49, 18, 28, 20,
		55, 30, 34, 11, 43, 14, 22,  4,
		62, 57, 46, 52, 38, 26, 32, 41,
		50, 36, 17, 19, 29, 10, 13, 21,
		56, 45, 25, 31, 35, 16,  9, 12,
		44, 24, 15,  8, 23,  7,  6,  5 };

	inline uint64_t log2OfMSB(uint64_t val) {
		// https://www.geeksforgeeks.org/find-significant-set-bit-number/
		// https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
		if (val == 0) {
			return 0;
		}

		val |= val >> 1;
		val |= val >> 2;
		val |= val >> 4;
		val |= val >> 8;
		val |= val >> 16;
		val |= val >> 32;
		return tab64[((uint64_t)((val - (val >> 1)) * 0x07EDD5E59A4E28C2)) >> 58];
	}

	///////////////////////////////////////////////////////
	template<typename T1, typename T2, typename E>
	inline typename unordered_map<T1,T2>::const_iterator checkIfInMap(const unordered_map<T1,T2>& map, const E& element) {
		try {
			auto res = map.find(element);
			if (res != map.end()) {
				return res;
			}
			else {
				cerr << "ERROR: in " << element << " " << typeid(map).name() << endl;
				throw runtime_error("ID not found! Maybe you got the wrong database or contentfile?");
			}
			return res;
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}

	///////////////////////////////////////////////////////
	template<typename T1, typename T2>
	inline void mapSetValue(unordered_map<T1, T2>& map, const T1& search, const T2& element) {
		try {
			auto res = map.find(search);
			if (res != map.end()) {
				res->second += element;
			}
			else {
				res->second = element;
			}
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}

	///////////////////////////////////////////////////////
	template<typename T1, typename T2, typename E>
	inline T2 checkIfInMapReturnZero(const unordered_map<T1, T2>& map, const E& element) {
		auto res = map.find(element);
		if (res != map.end()) {
			return res->second;
		}
		else {
			return static_cast<T2>(0);
		}
		return static_cast<T2>(0);
	}

	///////////////////////////////////////////////////////
	// Field of ones to mask the bits up until (2^31 - 1)
	const uint64_t _bitMasks[32] = { 0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535, 131071, 262143, 524287, 1048575, 2097151, 4194303, 8388607, 16777215, 33554431, 67108863, 134217727, 268435455, 536870911, 1073741823, 2147483647 };
	const uint64_t _bitmasksSingle_64[64] = {
		0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000,
		0x10000, 0x20000, 0x40000, 0x80000, 0x100000, 0x200000, 0x400000, 0x800000, 0x1000000, 0x2000000, 0x4000000, 0x8000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000,
		0x100000000, 0x200000000, 0x400000000, 0x800000000, 0x1000000000, 0x2000000000, 0x4000000000, 0x8000000000, 0x10000000000, 0x20000000000, 0x40000000000, 0x80000000000, 0x100000000000, 0x200000000000, 0x400000000000, 0x800000000000,
		0x1000000000000, 0x2000000000000, 0x4000000000000, 0x8000000000000, 0x10000000000000, 0x20000000000000, 0x40000000000000, 0x80000000000000, 0x100000000000000, 0x200000000000000, 0x400000000000000, 0x800000000000000, 0x1000000000000000, 0x2000000000000000, 0x4000000000000000, 0x8000000000000000
	};
	const uint32_t _bitmasksSingle_32[32] = {
		0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000,
		0x10000, 0x20000, 0x40000, 0x80000, 0x100000, 0x200000, 0x400000, 0x800000, 0x1000000, 0x2000000, 0x4000000, 0x8000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000
	};

	///////////////////////////////////////////////////////
	class sBitArray {
	private:
		unique_ptr<uint64_t[]> _arr;
		uint64_t _iNumOfBlocks = 0;
		uint64_t _iNumOfElements = 0;
		pair<uint64_t, uint8_t> _pPosOfFirstElem = make_pair(numeric_limits<uint64_t>::max(), static_cast<uint8_t>(64));

	public:
		////////////////////////////////
		sBitArray(const uint64_t& iSize) {
			_iNumOfBlocks = (iSize + 63) / 64;
			_arr.reset(new uint64_t[_iNumOfBlocks]);
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				_arr[i] = 0ULL;
			}
		}


		////////////////////////////////
		inline sBitArray(const sBitArray& ba) {
			_iNumOfBlocks = ba._iNumOfBlocks;
			_iNumOfElements = ba._iNumOfElements;
			_pPosOfFirstElem = ba._pPosOfFirstElem;
			if (_arr == nullptr) {
				_arr.reset(new uint64_t[_iNumOfBlocks]);
			}
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				_arr[i] = ba._arr[i];
			}
		}

		////////////////////////////////
		inline sBitArray(sBitArray&& ba) {
			_iNumOfBlocks = ba._iNumOfBlocks;
			_iNumOfElements = ba._iNumOfElements;
			_pPosOfFirstElem = ba._pPosOfFirstElem;
			if (_arr == nullptr) {
				/*_arr = ba._arr;
				ba._arr = nullptr;*/
				_arr.reset(new uint64_t[_iNumOfBlocks]);
			}
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				_arr[i] = ba._arr[i];
			}
		}

		////////////////////////////////
		inline sBitArray& operator = (const sBitArray& ba) {
			if (this != &ba) {
				_iNumOfBlocks = ba._iNumOfBlocks;
				_iNumOfElements = ba._iNumOfElements;
				_pPosOfFirstElem = ba._pPosOfFirstElem;
				if (_arr == nullptr) {
					_arr.reset(new uint64_t[_iNumOfBlocks]);
				}
				for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
					_arr[i] = ba._arr[i];
				}
			}
			return *this;
		}

		////////////////////////////////
		inline sBitArray& operator = (sBitArray&& ba) {
			if (this != &ba) {
				_iNumOfBlocks = ba._iNumOfBlocks;
				_iNumOfElements = ba._iNumOfElements;
				_pPosOfFirstElem = ba._pPosOfFirstElem;
				_arr.swap(ba._arr);
				/*if (_arr == nullptr) {
					_arr.reset(new uint64_t[_iNumOfBlocks]);
				}
				for (uint64_t i = 0; i < _iNumOfBlocks && (_arr[i] != 0 || ba._arr[i] != 0); ++i) {
					_arr[i] = ba._arr[i];
				}*/
			}
			return *this;
		}

		////////////////////////////////
		inline bool at(const uint64_t& idx) {
			const uint64_t& blockID = idx / 64;
			const uint64_t& fieldID = idx % 64;
			return ((_arr[blockID] >> fieldID) & 1);
		}

		////////////////////////////////
		/*inline uint64_t getUnique() {
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				//cout << _arr[i] << endl;
				const auto& bitPos = Utilities::log2OfMSB(_arr[i]);
				if (bitPos || _arr[i] == 1) {
					return bitPos + i * 64;
				}
			}
			return 0;
		}*/

		////////////////////////////////
		/*inline void setArr(const sBitArray& ba) {
			if (ba._pPosOfFirstElem < _pPosOfFirstElem) {
				_pPosOfFirstElem.first = ba._pPosOfFirstElem.first;
				_pPosOfFirstElem.second = ba._pPosOfFirstElem.second;
			}
			_iNumOfElements = ba._iNumOfElements;
			assert(ba._iNumOfBlocks == _iNumOfBlocks);
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				if (_arr[i] != 0 || ba._arr[i]) {
					_arr[i] |= ba._arr[i];
				}
			}
		}*/

		////////////////////////////////
		inline void set(const uint64_t& idx) {
			const uint64_t& blockID = idx / 64;
			const uint64_t& fieldID = idx % 64;
			assert(blockID < _iNumOfBlocks);
			// if there already is a 1, the counter should not increase
			_iNumOfElements += !((_arr[blockID] >> fieldID) & 1);
			_arr[blockID] |= (1ULL << fieldID);

			if (blockID * 64 + fieldID < _pPosOfFirstElem.first * 64 + _pPosOfFirstElem.second || _pPosOfFirstElem.first == numeric_limits<uint64_t>::max()) {
				_pPosOfFirstElem.first = blockID;
				_pPosOfFirstElem.second = static_cast<uint8_t>(fieldID);
			}
		}

		////////////////////////////////
		inline uint64_t numOfEntries() {
			/*if (_iNumOfElements == 0) {
				uint64_t iSum = 0;
				for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
					iSum += countBits(_arr[i]);
				}
				return iSum;
			}*/
			return _iNumOfElements;
		}

		inline void clear() {
			memset(_arr.get(), 0, _iNumOfBlocks * sizeof(uint64_t));
			_iNumOfElements = 0;
			_pPosOfFirstElem = make_pair(numeric_limits<uint64_t>::max(), static_cast<uint8_t>(64));
		}

		////////////////////////////////
		inline void spillInsides() {
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				if (_arr[i]) {
					for (uint8_t j = 0; j < 64; ++j) {
						if ((_arr[i] >> j) & 1) {
							cout << j + i * 64 << " ";
						}
					}
				}
			}
			cout << endl;
		}

		////////////////////////////////
		inline uint32_t sizeInBytes() const {
			return static_cast<uint32_t>(_iNumOfBlocks * 8 + sizeof(sBitArray));
		}

		// ranged based for loops
		class Iterator {
			uint64_t _CurrentVal = 0;
			uint64_t* _arrForIt;
			uint64_t _iNumOfBlocksIt = 0;
			uint64_t _iNumOfElements = 0, _iCounterOfElements = 0;
			pair<uint64_t, uint8_t> _CurrentPos = make_pair(0, 0);


		public:
			inline Iterator(unique_ptr<uint64_t[]>& arr, const uint64_t& blocks, const uint64_t& iNumOfElements, const pair<uint64_t, uint8_t>& pFirstPosition, const bool& end = 0) : _iNumOfBlocksIt(blocks), _iNumOfElements(iNumOfElements), _CurrentPos(pFirstPosition) {
				_arrForIt = arr.get();

				if (pFirstPosition.first != numeric_limits<uint64_t>::max() && !end) {
					_iCounterOfElements = 1;
					_CurrentVal = pFirstPosition.first * 64 + pFirstPosition.second;
				}
				else {
					_CurrentPos = make_pair(blocks, static_cast<uint8_t>(64));
				}
			}

			////////////////////////////////
			inline void SetNumOfEntries(const uint64_t& iNumOfElements) {
				_iNumOfElements = iNumOfElements;
			}

			////////////////////////////////
			inline void next() {
				if (_iCounterOfElements < _iNumOfElements) {
					uint8_t j = _CurrentPos.second + 1;
					for (uint64_t i = _CurrentPos.first; i < _iNumOfBlocksIt; ++i) {
						if (_arrForIt[i]) {
							for (; j < 64; ++j) {
								if (_arrForIt[i] & Utilities::_bitmasksSingle_64[j]) {
									_CurrentPos = make_pair(i, j);
									_CurrentVal = i * 64 + j;
									++_iCounterOfElements;
									return;
								}
							}
						}
						j = 0;
					}
				}
				_CurrentPos = make_pair(_iNumOfBlocksIt, static_cast<uint8_t>(64));
			}

			////////////////////////////////
			inline uint64_t& operator*() { return _CurrentVal; }

			////////////////////////////////
			inline Iterator& operator++() {
				next();
				return *this;
			}

			////////////////////////////////
			inline bool operator!=(const Iterator& it) {
				return _CurrentPos != it._CurrentPos;
			}
		};

		////////////////////////////////
		inline Iterator begin() {
			return Iterator(_arr, _iNumOfBlocks, _iNumOfElements, _pPosOfFirstElem);
		}

		////////////////////////////////
		inline Iterator end() {
			return Iterator(_arr, _iNumOfBlocks, _iNumOfElements, _pPosOfFirstElem, true);
		}
	};

	///////////////////////////////////////////////////////
	template<typename T>
	struct vector_set
	{
		using vec_type = vector<T>;
		using const_iterator = typename vec_type::const_iterator;

		vector_set() {
			_v.reserve(3);
		}
		vector_set(const uint64_t&) {
			_v.reserve(3);
			iSeenElem = -1;
		}
		vector_set(const pair<uint64_t, float>&) {
			_v.reserve(3);
			iSeenElem = make_pair(static_cast<uint64_t>(-1), 0.f);
		}

		vector_set(size_t max_size) {
			_v.reserve(max_size);
		}
		vector_set(size_t max_size, const uint64_t&) {
			_v.reserve(max_size);
			iSeenElem = static_cast<T>(-1);
		}
		vector_set(size_t max_size, const pair<uint64_t, float>&) {
			_v.reserve(max_size);
			iSeenElem = make_pair(static_cast<uint64_t>(-1), 0.f);
		}

		vec_type getV() {
			return _v;
		}

		bool insert(const T& elem) {
			if (elem == iSeenElem) {
				return false;
			}
			else {
				iSeenElem = elem;
				const auto it = lower_bound(_v.begin(), _v.end(), elem);
				itSeenDistance = it - _v.begin();
				if (_v.end() == it || *it != elem) {
					_v.insert(it, elem);
					return true;
				}
			}
			return false;
		}

		void insert_pair(const pair<uint64_t, float>& elem) {
			if (elem.first == iSeenElem.first) {
				//(_v.begin() + itSeenDistance)->second += elem.second;
				_v.at(itSeenDistance).second += elem.second;
			}
			else {
				iSeenElem.first = elem.first;
				const auto& it = lower_bound(_v.begin(), _v.end(), elem, [](const pair<uint64_t, float>& a, const pair<uint64_t, float>& val) { return a.first < val.first; });
				itSeenDistance = it - _v.begin();
				if (_v.end() == it || it->first != elem.first) {
					_v.insert(it, elem);
				}
				else {
					if (it->first == elem.first) {
						it->second += elem.second;
					}
				}
			}
		}

		void insert_pair_first(const pair<uint64_t, float>& elem) {
			iSeenElem.first = elem.first;
			itSeenDistance = 0;
			_v.insert(_v.begin(), elem);
		}

		void insert(const_iterator it, const T& elem) {
			_v.insert(it, elem);
		}

		void insert_first(const T& elem) {
			_v.insert(_v.begin(), elem);
		}

		auto find(const T& elem) const -> const_iterator {
			auto vend = _v.end();
			auto it = lower_bound(_v.begin(), vend, elem);
			if (it != vend && *it != elem) {
				it = vend;
			}
			return it;
		}

		float findScore(const uint64_t& elem) const {
			auto vend = _v.end();
			auto it = lower_bound(_v.begin(), vend, elem, [](const pair<uint64_t, float>& a, const uint64_t& val) { return a.first < val; });
			if (it != vend && it->first != elem) {
				it = vend; // signal that this shouldn't be
			}
			return it->second;
		}

		//bool contains(const T& elem) const {
		//	return find(elem) != _v.end();
		//}

		const_iterator begin() const {
			return _v.begin();
		}

		const_iterator end() const {
			return _v.end();
		}

		size_t size() const {
			return _v.size();
		}

		void clear(const uint64_t&) {
			if (_v.size()) {
				_v.clear();
			}
			itSeenDistance = 0;
			iSeenElem = static_cast<T>(-1);
		}

		void clear(const pair<uint64_t, float>&) {
			if (_v.size()) {
				_v.clear();
			}
			itSeenDistance = 0;
			iSeenElem = make_pair(static_cast<uint64_t>(-1), 0.f);
		}

	private:
		vec_type _v;
		size_t itSeenDistance;
		T iSeenElem;
	};
} // namespace