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
///////////////////////////////////////////////////////
// shamelessly stolen from "Timo"(https://stackoverflow.com/users/458825/timo): https://stackoverflow.com/questions/4351371/c-performance-challenge-integer-to-stdstring-conversion
// License for this code: https://creativecommons.org/licenses/by-sa/4.0/
#pragma once

#include "../MetaHeader.h"

namespace Utilities {

	const char digit_pairs[201] = {
	"00010203040506070809"
	"10111213141516171819"
	"20212223242526272829"
	"30313233343536373839"
	"40414243444546474849"
	"50515253545556575859"
	"60616263646566676869"
	"70717273747576777879"
	"80818283848586878889"
	"90919293949596979899"
	};

	static const int BUFFER_SIZE = 11;

	inline void itostr(int32_t val, string& s)
	{
		char buf[BUFFER_SIZE];
		char *it = &buf[BUFFER_SIZE - 2];

		if (val >= 0) {
			int div = val / 100;
			while (div) {
				memcpy(it, &digit_pairs[2 * (val - div * 100)], 2);
				val = div;
				it -= 2;
				div = val / 100;
			}
			memcpy(it, &digit_pairs[2 * val], 2);
			if (val < 10)
				it++;
		}
		else {
			int div = val / 100;
			while (div) {
				memcpy(it, &digit_pairs[-2 * (val - div * 100)], 2);
				val = div;
				it -= 2;
				div = val / 100;
			}
			memcpy(it, &digit_pairs[-2 * val], 2);
			if (val <= -10)
				it--;
			*it = '-';
		}

		s.append(it, &buf[BUFFER_SIZE] - it);
	}

	///////////////////////////////////////////////////////
	inline void itostr(uint32_t val, string& s)
	{
		char buf[BUFFER_SIZE];
		char *it = (char*)&buf[BUFFER_SIZE - 2];

		int div = val / 100;
		while (div) {
			memcpy(it, &digit_pairs[2 * (val - div * 100)], 2);
			val = div;
			it -= 2;
			div = val / 100;
		}
		memcpy(it, &digit_pairs[2 * val], 2);
		if (val < 10)
			it++;

		s.append((char*)it, (char*)&buf[BUFFER_SIZE] - (char*)it);
	}

	///////////////////////////////////////////////////////
	// modified version of the naive algorithm from https://github.com/miloyip/itoa-benchmark (MIT License)
	inline void u64toa_naive(uint64_t value, string& s) {
		char temp[20];
		char *p = temp;
		do {
			*p++ = char(value % 10) + '0';
			value /= 10;
		} while (value > 0);

		for (auto i = p - temp; i > 0; --i) {
			s += *--p;
		}
	}

	///////////////////////////////////////////////////////
	inline void itostr(const uint64_t& val, string& s) {
		if (val <= 4294967296) {
			itostr(static_cast<uint32_t>(val), s);
		}
		else {
			u64toa_naive(val, s);
		}
	}

	inline string itostr(const uint64_t& val) {
		string s = "";
		itostr(val, s);
		return s;
	}
}