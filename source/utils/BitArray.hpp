/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*
*  Copyright (C) 2020 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
*
**************************************************************************/
#pragma once

#include "../MetaHeader.h"

namespace Utilities {

	///////////////////////////////////////////////////////
	class sBitArray {
	private:
		unique_ptr<uint64_t[]> _arr;
		uint64_t _iNumOfBlocks = 0;
		vector<uint64_t> _vElementsArray;
		size_t _currentPositionInArray = 0;

	public:
		////////////////////////////////
		sBitArray(const uint64_t& iSize) {
			_iNumOfBlocks = (iSize + 63) >> 6; // idx / 64
			_arr.reset(new uint64_t[_iNumOfBlocks]);
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				_arr[i] = 0ULL;
			}
			_vElementsArray.resize(iSize);
		}


		////////////////////////////////
		inline sBitArray(const sBitArray& ba) {
			_iNumOfBlocks = ba._iNumOfBlocks;
			if (_arr == nullptr) {
				_arr.reset(new uint64_t[_iNumOfBlocks]);
			}
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				_arr[i] = ba._arr[i];
			}
			_vElementsArray = ba._vElementsArray;
			_currentPositionInArray = ba._currentPositionInArray;
		}

		////////////////////////////////
		inline sBitArray(sBitArray&& ba) {
			_iNumOfBlocks = ba._iNumOfBlocks;
			if (_arr == nullptr) {
				/*_arr = ba._arr;
				ba._arr = nullptr;*/
				_arr.reset(new uint64_t[_iNumOfBlocks]);
			}
			for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
				_arr[i] = ba._arr[i];
			}
			_vElementsArray = ba._vElementsArray;
			_currentPositionInArray = ba._currentPositionInArray;
		}

		////////////////////////////////
		inline sBitArray& operator = (const sBitArray& ba) {
			if (this != &ba) {
				_iNumOfBlocks = ba._iNumOfBlocks;
				if (_arr == nullptr) {
					_arr.reset(new uint64_t[_iNumOfBlocks]);
				}
				for (uint64_t i = 0; i < _iNumOfBlocks; ++i) {
					_arr[i] = ba._arr[i];
				}
				_vElementsArray = ba._vElementsArray;
				_currentPositionInArray = ba._currentPositionInArray;
			}
			return *this;
		}

		////////////////////////////////
		inline sBitArray& operator = (sBitArray&& ba) {
			if (this != &ba) {
				_iNumOfBlocks = ba._iNumOfBlocks;
				_arr.swap(ba._arr);
				_vElementsArray = ba._vElementsArray;
				_currentPositionInArray = ba._currentPositionInArray;
			}
			return *this;
		}

		////////////////////////////////
		inline bool at(const uint64_t& idx) {
			const uint64_t& blockID = idx >> 6; // idx / 64
			const uint64_t& fieldID = idx & 63; // idx % 64
			return ((_arr[blockID] >> fieldID) & 1);
		}

		////////////////////////////////
		inline void set(const uint64_t& idx) {
			const uint64_t& blockID = idx >> 6; // idx / 64
			const uint64_t& fieldID = idx & 63; // idx % 64

			// if there already is a 1, the counter should not increase
			if (!((_arr[blockID] >> fieldID) & 1)) {
				_vElementsArray[_currentPositionInArray++] = idx;
			}

			_arr[blockID] |= (1ULL << fieldID);
		}

		////////////////////////////////
		inline uint64_t numOfEntries() {
			return _currentPositionInArray;
		}

		inline void clear() {
			memset(_arr.get(), 0, _iNumOfBlocks * sizeof(uint64_t));
			_currentPositionInArray = 0;
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
			return static_cast<uint32_t>(_iNumOfBlocks * 8 + sizeof(sBitArray) + sizeof(uint64_t)*_vElementsArray.size());
		}

		////////////////////////////////
		inline vector<uint64_t>::const_iterator begin() const {
			return _vElementsArray.cbegin();
		}

		////////////////////////////////
		inline vector<uint64_t>::const_iterator end() const {
			return _vElementsArray.cbegin() + _currentPositionInArray;
		}
	};

}