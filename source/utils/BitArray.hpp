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
			_iNumOfBlocks = (iSize + 63) >> 6; // idx / 64
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
			const uint64_t& blockID = idx >> 6; // idx / 64
			const uint64_t& fieldID = idx & 63; // idx % 64
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
			const uint64_t& blockID = idx >> 6; // idx / 64
			const uint64_t& fieldID = idx & 63; // idx % 64
			assert(blockID < _iNumOfBlocks);
			// if there already is a 1, the counter should not increase
			_iNumOfElements += !((_arr[blockID] >> fieldID) & 1);
			_arr[blockID] |= (1ULL << fieldID);

			// blockID * 64 = blockID << 6
			if ((blockID << 6) + fieldID < (_pPosOfFirstElem.first << 6) + _pPosOfFirstElem.second || _pPosOfFirstElem.first == numeric_limits<uint64_t>::max()) {
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
					_CurrentVal = (pFirstPosition.first << 6) + pFirstPosition.second;
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
									_CurrentVal = (i << 6) + j;
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

}