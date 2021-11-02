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
#include <chrono>
#include <array>

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
		sOut.erase(std::remove(sOut.begin(), sOut.end(), c), sOut.end());
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
	template<typename FileType>
	inline string getFirstSequenceOfFile(FileType& sequenceFile) {
		string sLine = "";
		getline(sequenceFile, sLine); // In fasta as well as in fastq, the first line is uninteresting
		getline(sequenceFile, sLine); // This should now contain enough information to figure it out
		sequenceFile.seekg(0);
		return sLine;
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
						files.push_back(make_pair(fileName,int8_t(0)));
					}
					else {
						files.push_back(make_pair(fileName, int8_t(1)));
					}
				}
			}
#else
#if (__GNUC__ || defined(__llvm__)) && !_MSC_VER
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
	inline size_t getSizeOfFile(const string& file) {
		try {
#if (_WIN32 || _WIN64) && !__APPLE__ 

			return filesystem::file_size(file);
#else
#if (__GNUC__ || defined(__llvm__)) && !_MSC_VER

			struct stat stat_buf;
			stat(file.c_str(), &stat_buf);
			return stat_buf.st_size;
#endif
#endif
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}

	///////////////////////////////////////////////////////
	inline vector<pair<string, size_t>> gatherFilesAndSizesFromPath(const string& sPath) {
		try {
			vector<pair<string, size_t>> files;
#if (_WIN32 || _WIN64) && !__APPLE__ 
			if (sPath.back() == '/' || sPath.back() == '\\') {
				for (auto& fsPath : filesystem::directory_iterator(sPath)) {
					const string& fileName = (fsPath.path()).string();

					// read magic numbers at the beginning
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
						files.push_back(make_pair(fileName, filesystem::file_size(fsPath)));
					}
					else {
						files.push_back(make_pair(fileName, filesystem::file_size(fsPath)));
					}
				}
			}
#else
#if (__GNUC__ || defined(__llvm__)) && !_MSC_VER
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
	inline vector<string> splitAtFirstOccurence(const string& s, const char& delim) {
		vector<string> elems;
		auto res = find(s.begin(), s.end(), delim);
		if (res != s.end()) {
			elems.push_back(string(s.begin(), res));
			elems.push_back(string(res + 1, s.end()));
		}
		else {
			elems.push_back(s);
		}
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
		const string _newlineCharacters = string("\n");
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

		inline size_t size() const {
			return _vPointers.size();
		}

		inline float* operator[](const size_t& row) {
			return _vPointers[row].get();
		}

		inline const float* operator[](const size_t& row) const {
			return _vPointers[row].get();
		}

	};

#include <deque>
	class SparseMatrix {
	private:
		uint64_t _iNumOfRows;

		deque<float> _valArr;
		deque<uint64_t> _colsArr;
		vector<uint64_t> _rowsArr;

		void insert(uint64_t index, uint64_t row, uint64_t col, float val) {
			if (_valArr.empty()) {
				_valArr.push_back(val);
				_colsArr.push_back(col);
			}
			else {
				_valArr.insert(_valArr.begin() + index, val);
				_colsArr.insert(_colsArr.begin() + index, col);
			}

			for_each(_rowsArr.begin() + row + 1, _rowsArr.end(), [](uint64_t& n) { ++n; });
		}

	public:
		inline void generate(const uint64_t& rows) {
			if (_iNumOfRows != 0) {
				_valArr.clear();
				_colsArr.clear();
				_rowsArr.clear();
			}
			_iNumOfRows = rows;
			_rowsArr.resize(rows + 1);
			
		}

		inline size_t size() {
			return _rowsArr.size() - 1;
		}

		inline size_t size() const {
			return _rowsArr.size() - 1;
		}

		inline void set(const uint64_t& row, const uint64_t& col, const float& val) {

			uint64_t pos = _rowsArr.at(row);
			uint64_t currCol = numeric_limits<uint64_t>::max();

			for (; pos < _rowsArr.at(row + 1); pos++) {
				currCol = _colsArr.at(pos);

				if (currCol >= col) {
					break;
				}
			}

			if (currCol != col) {
				insert(pos, row, col, val);
			}
			else {
				_valArr.at(pos) += val;
			}
		}

		inline float get(const uint64_t& row, const uint64_t& col) {
			uint64_t currCol = 0;

			for (uint64_t pos = _rowsArr.at(row); pos < _rowsArr.at(row + 1); pos++) {
				currCol = _colsArr.at(pos);

				if (currCol == col) {
					return _valArr.at(pos);
				}
				else {
					if (currCol > col) {
						break;
					}
				}
			}

			return 0;
		}

		inline float get(const uint64_t& row, const uint64_t& col) const {
			uint64_t currCol = 0;

			for (uint64_t pos = _rowsArr.at(row); pos < _rowsArr.at(row + 1); pos++) {
				currCol = _colsArr.at(pos);

				if (currCol == col) {
					return _valArr.at(pos);
				}
				else {
					if (currCol > col) {
						break;
					}
				}
			}

			return 0;
		}
	};

	///////////////////////////////////////////////////////
	// Write strings to file, which are large enough
	class BufferedWriter {
	private:
		ofstream* _file = nullptr;
		size_t _iMaxBufferSize = 0;
		string _sStringToBeWritten = string("");
	
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
			//auto start = std::chrono::high_resolution_clock::now();
			_file->write(&_sStringToBeWritten[0], _sStringToBeWritten.size());
			//auto end = std::chrono::high_resolution_clock::now();
			//cout << "Written with " << _sStringToBeWritten.size() / (1024.0*1024.0*(chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(1000000000))) << " MB/s" << endl;
			_sStringToBeWritten = "";
		}

		inline void operator+=(const string& str) {
			_sStringToBeWritten += str;
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

	inline void checkIfFileCanBeCreated(const string& s) {
		if (!ofstream(s)) {
			throw runtime_error("File couldn't be created, maybe a wrong path was used?");
		}
	}

	inline ofstream createFileAndGetIt(const string& s) {
		ofstream f(s);
		if (!f) {
			throw runtime_error("File couldn't be created, maybe a wrong path was used?");
		}
		return f;
	}

	///////////////////////////////////////////////////////

	inline void copyFile(const string& sIn, const string& sOut) {
		// if rename is not working 
		// I took this solution from https://stackoverflow.com/users/1054324/peter: https://stackoverflow.com/questions/10195343/copy-a-file-in-a-sane-safe-and-efficient-way
		checkIfFileCanBeCreated(sOut);
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
	inline void moveFile(const string& sIn, const string& sOut) {
		if (rename(sIn.c_str(), sOut.c_str())) { // in C++17 we can change that to filesystem::rename
			cerr << "WARNING: Moving the file did not work, attempting to copy..." << endl;
			remove(sOut.c_str());
			Utilities::copyFile(sIn, sOut);
			remove(sIn.c_str());
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
			typename unordered_map<T1, T2>::const_iterator res = map.find(element);
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
	// Check if duplicate entries exist in the content file. 
	// WARNING: This will read the entire content file into RAM!
	inline void checkIfContentFileIsCorrupted(const string& sContentFile, const string& sFixedContentFileOut) {
		try {
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

			ifstream fContent;
			bool bTaxIdsAsStrings = false;
			//fContent.exceptions(std::ifstream::failbit | std::ifstream::badbit); 
			fContent.open(sContentFile);
			if (!fContent) {
				throw runtime_error("Content file couldn't be read!");
			}
			string sTempLine = "";
			unordered_map<string, tuple<string, string, string, string>> mOrganisms;
			uint64_t iLargestLineIndex = 1;
			while (getline(fContent, sTempLine)) {
				if (sTempLine != "") {
					const auto& tempLineContent = Utilities::split(sTempLine, '\t');
					if (tempLineContent.size() >= 5 && !bTaxIdsAsStrings) {
						bTaxIdsAsStrings = true;
					}

					bool bDummy = false;
					// check for dummys to get the current counter
					if (tempLineContent[0].find("EWAN") != string::npos) {
						bDummy = true;
					}
					auto entry = mOrganisms.find(tempLineContent[1]);
					if (entry != mOrganisms.end()) {
						if (!bDummy) {
							// content file is corrupted
							cout << "OUT: Content file is corrupted, duplicate entries" << tempLineContent[0] << " and " << get<0>(entry->second) << " were found. Merging them now..." << endl;

							const auto& res = func(tempLineContent[2], get<1>(entry->second), tempLineContent[3], get<2>(entry->second));

							entry->second = make_tuple(get<0>(entry->second), res.first, res.second, tempLineContent[4]);
						}
					}
					else {
						if (bTaxIdsAsStrings) {
							const auto& currLineIdx = std::stoull(tempLineContent[4]);
							iLargestLineIndex = (iLargestLineIndex < currLineIdx) ? currLineIdx : iLargestLineIndex;
							mOrganisms[tempLineContent[1]] = make_tuple(tempLineContent[0], tempLineContent[2], tempLineContent[3], tempLineContent[4]); // [taxID] -> (Name, species ID, Acc Nrs., lineIdx) 
						}
						else {
							mOrganisms[tempLineContent[1]] = make_tuple(tempLineContent[0], tempLineContent[2], tempLineContent[3], ""); // [taxID] -> (Name, species ID, Acc Nrs.) 
						}
					}
				}
			}

			// write content file
			ofstream fContentOS(sFixedContentFileOut);
			if (!fContentOS) {
				throw runtime_error("Content file couldn't be opened for writing!");
			}
			for (const auto& entry : mOrganisms) {
				fContentOS << get<0>(entry.second) << "\t" << entry.first << "\t" << get<1>(entry.second) << "\t" << get<2>(entry.second) << ((bTaxIdsAsStrings) ? ("\t" + get<3>(entry.second)) : "") << endl;
			}
		}
		catch (...) {
			cerr << "ERROR: in: " << __PRETTY_FUNCTION__ << endl; throw;
		}
	}

	///////////////////////////////////////////////////////
	template<typename T1, typename T2>
	inline uint64_t calculateSizeInByteOfUnorderedMap(const unordered_map<T1,T2>& map) {
		uint64_t iSizeInBytes = 0;
		for (size_t i = 0; i < map.bucket_count(); ++i) {
			const size_t& bucket_size = map.bucket_size(i);
			if (bucket_size == 0) {
				iSizeInBytes += sizeof(pair<T1, T2>);
			}
			else {
				iSizeInBytes += sizeof(pair<T1, T2>) * bucket_size;
			}
		}
		return iSizeInBytes;
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
	// Shitty conversion function because using real json is too much of a hassle in C++
	// Maybe, someday, I'll use https://github.com/nlohmann/json
	inline void readParametersFromYaml(const string& sYamlFile) {
		if (!ifstream(sYamlFile)) {
			throw runtime_error("Config file not found!");
		}
		ifstream jsonFile(sYamlFile);
		string sLine = "";
		while (getline(jsonFile, sLine)) {
			if (sLine.empty()) {
				continue;
			}

			sLine = Utilities::lstrip(sLine);
			if (sLine.front() == '#') {
				continue;
			}
			
			auto parameterPair = Utilities::splitAtFirstOccurence(sLine, ':');
			if (parameterPair.size() > 1) {
				// long ugly list 
				parameterPair[0] = Utilities::removeCharFromString(parameterPair[0], '"');
				parameterPair[1] = Utilities::removeCharFromString(parameterPair[1], '"');
				parameterPair[1] = Utilities::lstrip(parameterPair[1]);

				if (parameterPair[0] == "Mode") {
					GlobalInputParameters.cMode = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "Index") {
					GlobalInputParameters.indexFile = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "ContentFile") {
					GlobalInputParameters.contentFileIn = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "kHigh") {
					GlobalInputParameters.iHigherK = stoi(parameterPair[1]);
					GlobalInputParameters.iHigherK = (GlobalInputParameters.iHigherK > HIGHESTPOSSIBLEK) ? HIGHESTPOSSIBLEK : GlobalInputParameters.iHigherK;
					GlobalInputParameters.bHighKSetByUser = true;
					continue;
				}
				if (parameterPair[0] == "kLow") {
					GlobalInputParameters.iLowerK = stoi(parameterPair[1]);
					GlobalInputParameters.iLowerK = (GlobalInputParameters.iLowerK < 1) ? 1 : GlobalInputParameters.iLowerK;
					continue;
				}
				if (parameterPair[0] == "NumberOfThreads") {
					GlobalInputParameters.iNumOfThreads = stoi(parameterPair[1]);
					continue;
				}
				if (parameterPair[0] == "AvailableRAMinGB") {
					GlobalInputParameters.iMemorySizeAvail = 1024ull * stoi(parameterPair[1]);
					GlobalInputParameters.bCustomMemorySet = true;
					continue;
				}
				if (parameterPair[0] == "FilePathForTemporaryFiles") {
					GlobalInputParameters.sTempPath = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "CallIndex") {
					GlobalInputParameters.iNumOfCall = stoi(parameterPair[1]);
					continue;
				}
				if (parameterPair[0] == "Verbose") {
					GlobalInputParameters.bVerbose = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "AlphabetFile") {
					GlobalInputParameters.codonTable = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "AlphabetIndex") {
					GlobalInputParameters.sCodonID = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "InputFileOrFolder") {
					GlobalInputParameters.sInput = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "PairedEnd-First") {
					GlobalInputParameters.sPairedEnd1 = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "PairedEnd-Second") {
					GlobalInputParameters.sPairedEnd2 = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "AlreadyTranslated") {
					GlobalInputParameters.bTranslated = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "TaxonomicLevel") {
					GlobalInputParameters.sTaxLevel = parameterPair[1];
					if (GlobalInputParameters.sTaxLevel == "sequence") {
						GlobalInputParameters.sTaxLevel = "lowest";
					}
					continue;
				}
				if (parameterPair[0] == "AccessionToTaxIDFileOrFolder") {
					GlobalInputParameters.sAccToTaxFiles = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "TaxonomyFolder") {
					GlobalInputParameters.sTaxonomyPath = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "TaxIDsAreStrings") {
					GlobalInputParameters.bTaxIdsAsStrings = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "One") {
					GlobalInputParameters.bOnlyOneFrame = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "Three") {
					GlobalInputParameters.bThreeFrames = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "Six") {
					GlobalInputParameters.bSixFrames = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "ProfileOutputfile") {
					GlobalInputParameters.tableFile = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "ReadIDtoTaxIDOutputfile") {
					GlobalInputParameters.readToTaxaFile = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "ReadIDtoTaxIDOutputFormat") {
					if (parameterPair[1] == "") {
						continue;
					}
					if (parameterPair[1] == "json") {
						GlobalInputParameters.outputFormat = 1;
						continue;
					}
					if (parameterPair[1] == "jsonl") {
						GlobalInputParameters.outputFormat = 2;
						continue;
					}
					if (parameterPair[1] == "kraken") {
						GlobalInputParameters.outputFormat = 0;
						continue;
					}
					if (parameterPair[1] == "tsv") {
						GlobalInputParameters.outputFormat = 3;
						continue;
					}
					continue;
				}
				if (parameterPair[0] == "UseRAMOnly") {
					GlobalInputParameters.bRAM = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "NumberOfTaxaPerRead") {
					GlobalInputParameters.iNumOfBeasts = stoi(parameterPair[1]);
					continue;
				}
				if (parameterPair[0] == "UniqueKmersOnly") {
					GlobalInputParameters.bUnique = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "ThresholdForScore") {
					GlobalInputParameters.threshold = stof(parameterPair[1]);
					continue;
				}
				if (parameterPair[0] == "PrintCoverage") {
					GlobalInputParameters.bCoverage = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "Filter") {
					auto filterFiles = Utilities::split(parameterPair[1], ' ');
					if (filterFiles[0] != "_" || filterFiles[1] != "_") {
						GlobalInputParameters.bFilter = true;
						GlobalInputParameters.sFilteredCleanOut = Utilities::removeSpaceAndEndline(filterFiles[0]);
						GlobalInputParameters.sFilteredContaminantsOut = Utilities::removeSpaceAndEndline(filterFiles[1]);
					}
					continue;
				}
				if (parameterPair[0] == "ErrorThreshold") {
					GlobalInputParameters.fErrorThreshold = stof(parameterPair[1]);
					continue;
				}
				if (parameterPair[0] == "Gzip") {
					GlobalInputParameters.bGzipOut = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "FileWithDeletedTaxa") {
					GlobalInputParameters.delnodesFile = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "ShrinkingStrategy") {
					GlobalInputParameters.shrinkStrategy = stoi(parameterPair[1]);
					continue;
				}
				if (parameterPair[0] == "ShrinkPercentage") {
					GlobalInputParameters.fPercentageOfThrowAway = stof(parameterPair[1]);
					continue;
				}
				if (parameterPair[0] == "ContentFile-First") {
					GlobalInputParameters.contentFile1 = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "ContentFile-Second") {
					GlobalInputParameters.contentFile2 = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "ContentFile-Out") {
					GlobalInputParameters.contentFileAfterUpdate = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "FirstOldIndex") {
					GlobalInputParameters.firstOldIndex = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "SecondOldIndex") {
					GlobalInputParameters.secondOldIndex = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "NewIndex") {
					GlobalInputParameters.sDBPathOut = parameterPair[1];
					continue;
				}
				if (parameterPair[0] == "Debug") {
					_bShowLineDebugMode = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "Visualize") {
					GlobalInputParameters.bVisualize = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "Spaced") {
					GlobalInputParameters.bSpaced = (parameterPair[1] == "true") ? true : false;
					continue;
				}
				if (parameterPair[0] == "SpacedMaskIdx") {
					GlobalInputParameters.iNumOfMask = stoi(parameterPair[1]);
					continue;
				}
				//if (parameterPair[0] == "Mode") {
				//	GlobalInputParameters.cMode = parameterPair[1];
				//	continue;
				//}

			}
		}
	}

} // namespace