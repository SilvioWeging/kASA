/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*
*  Copyright (C) 2020 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
*  
*  Some parts of this file are taken from http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
*  Parts from Stackoverflow are marked as such.
*
**************************************************************************/
#pragma once

#include <initializer_list>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <utility>
#include <unordered_map>

#include "WorkerThread.hpp"
#include "BitArray.hpp"
#include "../MetaHeader.h"

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
	inline string replaceCharacter(const string& sIn, const char& c, const char& rep) {
		string sOut = sIn;
		std::replace(sOut.begin(), sOut.end(), c, rep);
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
		if (sIn.size() == 0) {
			return sIn;
		}
		if (sIn.size() == 1 && sIn[0] == c) {
			return "";
		}
		string sOut = sIn;
		auto it = sOut.begin();
		while (*it == c) {
			it = sOut.erase(it);
		}
		return sOut;
	}

	///////////////////////////////////////////////////////
	inline string rstrip(const string& sIn, const char& c) {
		if (sIn.size() == 0) {
			return sIn;
		}
		if (sIn.size() == 1 && sIn[0] == c) {
			return "";
		}
		string sOut = sIn;
		auto it = sOut.end() - 1;
		while (*it == c) {
			it = sOut.erase(it);
			it--;
		}
		return sOut;
	}

	///////////////////////////////////////////////////////
	inline pair<vector<pair<string,int8_t>>, size_t> gatherFilesFromPath(const string& sPath) {
		try {
			vector<pair<string, int8_t>> files; // 0: not zipped; 1: gzipped; 2: bzip2'ed
			size_t overallFileSize = 0;
#if _WIN32 || _WIN64 
			if (sPath.back() == '/' || sPath.back() == '\\') {
				for (auto& fsPath : filesystem::directory_iterator(sPath)) {
					const string& fileName = (fsPath.path()).string();
					
					// read magic numbers at the beginning
					ifstream prelimFile(fileName);
					char firstByte;
					prelimFile.read(&firstByte, 1);
					prelimFile.close();
					bool isGzipped = false;
					if ((firstByte == 0x1f)) {
						isGzipped = true;
					}
					if (firstByte == 0x42) {
						throw runtime_error("ERROR: File was compressed with bzip2. kASA can't handle this at the moment, sorry.");
					}

					if (!isGzipped) {
						ifstream fast_q_a_File;
						//fast_q_a_File.exceptions(std::ifstream::failbit | std::ifstream::badbit);
						fast_q_a_File.open(fileName);
						fast_q_a_File.seekg(0, fast_q_a_File.end);
						overallFileSize += fast_q_a_File.tellg();
						files.push_back(make_pair(fileName,int8_t(0)));
					}
					else {
						files.push_back(make_pair(fileName, int8_t(1)));
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

							ifstream prelimFile(fileName);
							char firstByte;
							prelimFile.read(&firstByte, 1);
							prelimFile.close();
							bool isGzipped = false;
							if (firstByte == 0x1f) {
								isGzipped = true;
							}
							if (firstByte == 0x42) {
								throw runtime_error("ERROR: File was compressed with bzip2. kASA can't handle this at the moment, sorry.");
							}

							if (!isGzipped) {
								ifstream fast_q_a_File;
								//fast_q_a_File.exceptions(std::ifstream::failbit | std::ifstream::badbit);
								fast_q_a_File.open(fileName);
								fast_q_a_File.seekg(0, fast_q_a_File.end);
								overallFileSize += fast_q_a_File.tellg();
								files.push_back(make_pair(fileName, int8_t(0)));
							}
							else {
								files.push_back(make_pair(fileName, int8_t(1)));
							}
						}
					}

					closedir(dirp);
				}
			}
#endif
#endif
			else {
				ifstream prelimFile(sPath);
				char firstByte;
				prelimFile.read(&firstByte, 1);
				prelimFile.close();
				bool isGzipped = false;
				if (firstByte == 0x1f) {
					isGzipped = true;
				}
				if (firstByte == 0x42) {
					throw runtime_error("ERROR: File was compressed with bzip2. kASA can't handle this at the moment, sorry.");
				}
				if (!isGzipped) {
					ifstream fast_q_a_File;
					//fast_q_a_File.exceptions(std::ifstream::failbit | std::ifstream::badbit);
					fast_q_a_File.open(sPath);
					fast_q_a_File.seekg(0, fast_q_a_File.end);
					overallFileSize += fast_q_a_File.tellg();
					files.push_back(make_pair(sPath, int8_t(0)));
				}
				else {
					files.push_back(make_pair(sPath, int8_t(1)));
				}
			}

			return make_pair(files,overallFileSize);
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}

	///////////////////////////////////////////////////////
	inline vector<pair<string, size_t>> gatherFilesAndSizesFromPath(const string& sPath) {
		try {
			vector<pair<string, size_t>> files;
#if _WIN32 || _WIN64 
			if (sPath.back() == '/' || sPath.back() == '\\') {
				for (auto& fsPath : filesystem::directory_iterator(sPath)) {
					const string& fileName = (fsPath.path()).string();

					// read magic numbers at the beginning
					ifstream prelimFile(fileName);
					char firstByte;
					prelimFile.read(&firstByte, 1);
					prelimFile.close();
					bool isGzipped = false;
					if ((firstByte == 0x1f)) {
						isGzipped = true;
					}
					if (firstByte == 0x42) {
						throw runtime_error("ERROR: File was compressed with bzip2. kASA can't handle this at the moment, sorry.");
					}

					if (!isGzipped) {
						files.push_back(make_pair(fileName, filesystem::file_size(fsPath)));
					}
					else {
						files.push_back(make_pair(fileName, filesystem::file_size(fsPath)));
					}
				}
			}
#else
#if __GNUC__ || defined(__llvm__)
			if (sPath.back() == '/') {
				DIR* dirp;
				struct dirent* directory;

				dirp = opendir(sPath.c_str());
				if (dirp) {
					while ((directory = readdir(dirp)) != NULL)
					{
						const string fName(directory->d_name);
						if (fName != "." && fName != "..") {
							const string& fileName = sPath + fName;

							ifstream prelimFile(fileName);
							char firstByte;
							prelimFile.read(&firstByte, 1);
							prelimFile.close();
							bool isGzipped = false;
							if (firstByte == 0x1f) {
								isGzipped = true;
							}
							if (firstByte == 0x42) {
								throw runtime_error("ERROR: File was compressed with bzip2. kASA can't handle this at the moment, sorry.");
							}

							if (!isGzipped) {
								struct stat stat_buf;
								stat(sPath.c_str(), &stat_buf);
								files.push_back(make_pair(fileName, stat_buf.st_size));
							}
							else {
								struct stat stat_buf;
								stat(sPath.c_str(), &stat_buf);
								files.push_back(make_pair(fileName, stat_buf.st_size));
							}
						}
					}

					closedir(dirp);
				}
			}
#endif
#endif
			else {
				ifstream prelimFile(sPath);
				char firstByte;
				prelimFile.read(&firstByte, 1);
				prelimFile.close();
				bool isGzipped = false;
				if (firstByte == 0x1f) {
					isGzipped = true;
				}
				if (firstByte == 0x42) {
					throw runtime_error("ERROR: File was compressed with bzip2. kASA can't handle this at the moment, sorry.");
				}
				if (!isGzipped) {
					struct stat stat_buf;
					stat(sPath.c_str(), &stat_buf);
					files.push_back(make_pair(sPath, stat_buf.st_size));
				}
				else {
					struct stat stat_buf;
					stat(sPath.c_str(), &stat_buf);
					files.push_back(make_pair(sPath, stat_buf.st_size));
				}
			}

			return files;
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
	class FileReader {
	private:
		static const size_t _bufferSize = 2048;
		char _bufferArrayForInput[_bufferSize];
		char* _startOfBuffer = &_bufferArrayForInput[0];
		const string _newlineCharacters = "\n";
		size_t _currendPosition = 0, _maxCharsRead = 0;
		T* _file = nullptr;

		inline bool readFromFile() {
			if (_currendPosition >= _maxCharsRead) {
				_file->read(_bufferArrayForInput, _bufferSize);
				_maxCharsRead = _file->gcount();
				_currendPosition = 0;

				// Nothing to read, file is either corrupted or at the end
				if (_maxCharsRead == 0) {
					if (_file->fail() && !_file->eof()) {
						_file->clear();
						string dummyString = "";
						getline(*_file, dummyString);
					}
					else if (_file->eof()) {
						return false;
					}
				}
				//  Reset state to continue reading characters from buffer
				if (_maxCharsRead < _bufferSize) {
					_file->clear();
					_bufferArrayForInput[_maxCharsRead] = '\n';
				}
			}

			return true;
		}

	public:

		// Standard contructor is sufficient

		inline bool notNull() {
			return _file != nullptr;
		}

		inline void setFile(T* input) {
			_file = input;
			memset(_bufferArrayForInput, ' ', _bufferSize);
		}

		inline void ignore(uint64_t& iNumOfChars) {
			pair<std::string, bool> dummyObject("",false);
			uint64_t iNumOfDummyChars = 0; 
			while (!dummyObject.second) {
				getChunk(dummyObject, iNumOfDummyChars);
				iNumOfChars += iNumOfDummyChars;
			}
		}

		inline bool eof() {
			return _file->eof();
		}

		inline void getChunk(pair<std::string, bool>& outPair, uint64_t& iNumOfChars) {
			iNumOfChars = 0;
			if (readFromFile()) {
				auto newlineCharPos = find_first_of(_startOfBuffer + _currendPosition, _startOfBuffer + _bufferSize, _newlineCharacters.cbegin(), _newlineCharacters.cend());
				if (newlineCharPos != _startOfBuffer + _bufferSize) {
					outPair.first.append(_startOfBuffer + _currendPosition, newlineCharPos);
					iNumOfChars = newlineCharPos - &_bufferArrayForInput[_currendPosition] + 1;
					_currendPosition = newlineCharPos - _startOfBuffer + 1;
					outPair.second = true;
				}
				else {
					outPair.first.append(_startOfBuffer + _currendPosition, _startOfBuffer + _bufferSize);
					iNumOfChars = newlineCharPos - &_bufferArrayForInput[_currendPosition];
					_currendPosition = _bufferSize;
					outPair.second = false;
				}
			}
		}
	};

	///////////////////////////////////////////////////////
	// Class to simulate 2D Arrays but with contiguous memory and constant element access
	// DONT USE WITH ANYTHING ELSE THAN BASIC TYPES LIKE int, float, ...
	template<typename T>
	class Vector2D {
	private:
		vector<T> _vValues;
		vector<T*> _vPointers;

	public:

		inline void setSize(const size_t& rows, const size_t& cols) {
			setZero();
			_vValues.resize(rows * cols);
			_vPointers.clear();
			for (size_t i = 0; i < rows; ++i) {
				_vPointers.push_back(&_vValues[0] + i * cols);
			}
		}

		inline void setZero() {
			if (_vValues.size()) {
				memset(&_vValues[0], 0, sizeof(T)*_vValues.size());
			}
		}

		inline size_t size() {
			return _vValues.size();
		}

		inline T* operator[](const size_t& row) {
			return _vPointers[row];
		}

		inline const T* operator[](const size_t& row) const {
			return _vPointers[row];
		}

	};

	class Non_contiguousArray {
	private:
		vector<unique_ptr<float[]>> _vPointers;

	public:
		// we must assume, that cols stays constant
		inline void generate(const size_t& rows, const size_t& cols) {
			if (_vPointers.size()) {
				while (rows >= _vPointers.size()) {
					_vPointers.push_back(unique_ptr<float[]>(new float[cols]));
				}
				for (size_t i = 0; i < rows; ++i) {
					memset(_vPointers[i].get(), 0, sizeof(float) * cols);
				}
			}
			else {
				_vPointers.resize(rows);
				for (size_t i = 0; i < rows; ++i) {
					_vPointers[i].reset(new float[cols]);
					memset(_vPointers[i].get(), 0, sizeof(float) * cols);
				}
			}
		}

		inline size_t size() {
			return _vPointers.size();
		}

		inline float* operator[](const size_t& row) {
			return _vPointers[row].get();
		}

		inline const float* operator[](const size_t& row) const {
			return _vPointers[row].get();
		}

	};

	///////////////////////////////////////////////////////
	// Write strings to file, which are large enough
	class BufferedWriter {
	private:
		ofstream* _file = nullptr;
		size_t _iMaxBufferSize = 0;
		string _sStringToBeWritten = "";
	
	public:

		BufferedWriter(ofstream& ofstr, const size_t& bufferSize) {
			_file = &ofstr;
			_iMaxBufferSize = bufferSize;
			_sStringToBeWritten.reserve(bufferSize + bufferSize/1024);
		}

		~BufferedWriter() {
			writeToFile();
		}

		inline void writeToFile() {
			_file->write(&_sStringToBeWritten[0], _sStringToBeWritten.size());
			_sStringToBeWritten = "";
		}

		inline void operator+=(const string& str) {
			_sStringToBeWritten.append(str);
			if (_sStringToBeWritten.size() >= _iMaxBufferSize) {
				writeToFile();
			}
		}

		inline void operator+=(const char* const charArr) {
			_sStringToBeWritten.append(charArr);
			if (_sStringToBeWritten.size() >= _iMaxBufferSize) {
				writeToFile();
			}
		}

		inline void operator+=(const char& _Ch) {
			_sStringToBeWritten.push_back(_Ch);
			if (_sStringToBeWritten.size() >= _iMaxBufferSize) {
				writeToFile();
			}
		}

		inline string& getString() {
			return _sStringToBeWritten;
		}
	};

	///////////////////////////////////////////////////////
	template<typename T>
	inline void getChunk(T&& input, pair<std::string, bool>&& outPair, char* buffer, uint64_t&& iNumOfChars = 0) {
		input.get(buffer, 1000); // .get and .getline are way slower than .read but .read would ignore the newline character. only mmap would be faster albeit more memory intensive
		iNumOfChars = input.gcount();
		bool bLineFinished = false;
		if (iNumOfChars != 999 && input) {
			bLineFinished = true;
			input.get();
		}

		if (iNumOfChars == 0) {
			if (input.fail() && !input.eof()) {
				input.clear();
				getline(input, outPair.first);
			}
		}
		outPair.first.append(buffer, buffer + iNumOfChars);
		outPair.second = bLineFinished;
	}


	///////////////////////////////////////////////////////

	inline void createFile(const string& s) {
		if (!ofstream(s)) {
			throw runtime_error("File couldn't be created, maybe a wrong path was used?");
		}
	}

	///////////////////////////////////////////////////////

	inline void copyFile(const string& sIn, const string& sOut) {
		// rename is bugging under linux in that it creates a segmentation fault if the file is on another drive
		// so we just take this solution from https://stackoverflow.com/users/1054324/peter: https://stackoverflow.com/questions/10195343/copy-a-file-in-a-sane-safe-and-efficient-way
		createFile(sOut);
		ifstream source(sIn, ios::binary);
		ofstream dest(sOut, ios::binary);

		istreambuf_iterator<char> begin_source(source);
		istreambuf_iterator<char> end_source;
		ostreambuf_iterator<char> begin_dest(dest);
		copy(begin_source, end_source, begin_dest);

		source.close();
		dest.close();
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
	// https://www.geeksforgeeks.org/find-significant-set-bit-number/
	// https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers (https://stackoverflow.com/users/944687/desmond-hume)
	// License for this code: https://creativecommons.org/licenses/by-sa/4.0/
	
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


} // namespace