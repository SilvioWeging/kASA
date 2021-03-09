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

#include <cstdint>
#include <tuple>
#include "uint128_t.hpp"

// Check windows
#if (_WIN32 || _WIN64) && !__APPLE__

#pragma pack(push, 1)
struct packedPair {
	uint32_t first = 0;
	uint16_t second = 0;
	packedPair() {}
	packedPair(const uint32_t& a, const uint16_t& b) : first(a), second(b) {}
};
#pragma pack(pop)

#pragma pack(push, 1)
struct packedBigPair {
	uint64_t first = 0; uint32_t second = 0;

	packedBigPair() {}
	packedBigPair(const uint64_t& a, const uint32_t& b) : first(a), second(b) {}

	packedBigPair& operator=(const std::tuple<uint64_t, uint32_t>& a) {
		first = std::get<0>(a);
		second = std::get<1>(a);
		return *this;
	}

	bool operator==(const packedBigPair& a) const {
		return a.first == this->first && a.second == this->second;
	}

	bool operator<(const packedBigPair& b) const {
		return (this->first < b.first || (!(b.first < this->first) && this->second < b.second));
	}
};
#pragma pack(pop)

#pragma pack(push, 1)
struct packedLargePair {
	uint128_t first = 0; uint32_t second = 0;

	packedLargePair() {}
	packedLargePair(const uint128_t& a, const uint32_t& b) : first(a), second(b) {}

	packedLargePair& operator=(const std::tuple<uint128_t, uint32_t>& a) {
		first = std::get<0>(a);
		second = std::get<1>(a);
		return *this;
	}

	bool operator==(const packedLargePair& a) const {
		return a.first == this->first && a.second == this->second;
	}

	bool operator<(const packedLargePair& b) const {
		return (this->first < b.first || (!(b.first < this->first) && this->second < b.second));
	}

};
#pragma pack(pop)

#pragma pack(push, 1)
struct packedBigPairTrie {
	uint64_t second = 0; uint32_t first = 0;

	packedBigPairTrie() {
		second = 0; first = 0;
	}

	packedBigPairTrie(const uint32_t& f, const uint64_t& s) : second(s), first(f) {}
};
#pragma pack(pop)

#endif


// Check GCC or Clang
#if (__GNUC__ || __clang__) && !_MSC_VER
struct __attribute__((packed)) packedPair {
	uint32_t first = 0;
	uint16_t second = 0;
	packedPair() {}
	packedPair(const uint32_t & a, const uint16_t & b) : first(a), second(b) {}
};

struct __attribute__((packed)) packedBigPair {
	uint64_t first = 0; uint32_t second = 0;

	packedBigPair() {}
	packedBigPair(const uint64_t & a, const uint32_t & b) : first(a), second(b) {}

	packedBigPair& operator=(const std::tuple<uint64_t, uint32_t> & a) {
		first = std::get<0>(a);
		second = std::get<1>(a);
		return *this;
	}

	bool operator==(const packedBigPair & a) const {
		return a.first == this->first && a.second == this->second;
	}

	bool operator<(const packedBigPair & b) const {
		return (this->first < b.first || (!(b.first < this->first) && this->second < b.second));
	}
};

struct __attribute__((packed)) packedLargePair {
	uint128_t first; uint32_t second = 0;

	packedLargePair() {}
	packedLargePair(const uint128_t& a, const uint32_t& b) : first(a), second(b) {}

	packedLargePair& operator=(const std::tuple<uint128_t, uint32_t>& a) {
		first = std::get<0>(a);
		second = std::get<1>(a);
		return *this;
	}

	bool operator==(const packedLargePair& a) const {
		return a.first == this->first && a.second == this->second;
	}

	bool operator<(const packedLargePair& b) const {
		return (this->first < b.first || (!(b.first < this->first) && this->second < b.second));
	}
};

struct __attribute__((packed)) packedBigPairTrie {
	uint64_t second = 0; uint32_t first = 0;

	packedBigPairTrie() {
		second = 0;
		first = 0;
	}

	packedBigPairTrie(const uint32_t & f, const uint64_t & s) : second(s), first(f) {}
};

#endif