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
// Copied from https://github.com/calccrypto/uint128_t and modified for use in this software
// License for this code: https://opensource.org/licenses/MIT
// 
#pragma once

#ifndef __UINT128_T__
#define __UINT128_T__

#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/*#if defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN || \
    defined(__BIG_ENDIAN__) ||                               \
    defined(__ARMEB__) ||                                    \
    defined(__THUMBEB__) ||                                  \
    defined(__AARCH64EB__) ||                                \
    defined(_MIBSEB) || defined(__MIBSEB) || defined(__MIBSEB__)
#ifndef __BIG_ENDIAN__
#define __BIG_ENDIAN__
#endif
#elif defined(__BYTE_ORDER) && __BYTE_ORDER == __LITTLE_ENDIAN || \
    defined(__LITTLE_ENDIAN__) ||                                 \
    defined(__ARMEL__) ||                                         \
    defined(__THUMBEL__) ||                                       \
    defined(__AARCH64EL__) ||                                     \
    defined(_MIPSEL) || defined(__MIPSEL) || defined(__MIPSEL__)
#ifndef __LITTLE_ENDIAN__
#define __LITTLE_ENDIAN__
#endif
#endif*/ // this doesn't work

#ifndef __LITTLE_ENDIAN__
#define __LITTLE_ENDIAN__ // as if kASA would compile on a big-endian :DD
#endif

// !=, ==, >>, <<, < , >

class uint128_t;

// Give uint128_t type traits
namespace std {  // This is probably not a good idea
    template <> struct is_arithmetic <uint128_t> : std::true_type {};
    template <> struct is_integral   <uint128_t> : std::true_type {};
    template <> struct is_unsigned   <uint128_t> : std::true_type {};
}

#if __GNUC__ || __clang__
class __attribute__((packed)) uint128_t {
#else
class uint128_t {
#endif

private:
#ifdef __BIG_ENDIAN__
    uint64_t UPPER, LOWER;
#endif
#ifdef __LITTLE_ENDIAN__
    uint64_t LOWER = 0, UPPER = 0;
#endif

public:
    // Constructors
    uint128_t() = default;
    uint128_t(const uint128_t& rhs) = default;
    uint128_t(uint128_t&& rhs) = default;

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t(const T& rhs)
#ifdef __BIG_ENDIAN__
        : UPPER(0), LOWER(rhs)
#endif
#ifdef __LITTLE_ENDIAN__
        : LOWER(rhs), UPPER(0)
#endif
    {
        if (std::is_signed<T>::value) {
            if (rhs < 0) {
                UPPER = -1;
            }
        }
    }

    template <typename S, typename T, typename = typename std::enable_if <std::is_integral<S>::value&& std::is_integral<T>::value, void>::type>
    uint128_t(const S& upper_rhs, const T& lower_rhs)
#ifdef __BIG_ENDIAN__
        : UPPER(upper_rhs), LOWER(lower_rhs)
#endif
#ifdef __LITTLE_ENDIAN__
        : LOWER(lower_rhs), UPPER(upper_rhs)
#endif
    {}

    //  RHS input args only

    // Assignment Operator
    uint128_t& operator=(const uint128_t& rhs) = default;
    uint128_t& operator=(uint128_t&& rhs) = default;

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator=(const T& rhs) {
        UPPER = 0;

        if (std::is_signed<T>::value) {
            if (rhs < 0) {
                UPPER = -1;
            }
        }

        LOWER = rhs;
        return *this;
    }

    // Typecast Operators
    uint128_t(std::string& s) {
        init(s.c_str());
    }

    uint128_t(const char* s) {
        init(s);
    }

    void init(const char* s) {
        if (s == NULL || s[0] == 0) { uint128_t(); return; }
        if (s[1] == 'x')
            s += 2;
        else if (*s == 'x')
            s++;

        UPPER = ConvertToUint64(s);
        LOWER = ConvertToUint64(s + 16);
    }

    uint64_t ConvertToUint64(const char* s) const {
        int count = 0;
        uint64_t val = 0;
        uint8_t hv = HexToInt(s++);
        while (hv != 0xFF && count < 16) {
            val = (val << 4) | hv;
            hv = HexToInt(&s[count]);
            count++;
        }
        return val;
    }

    uint8_t HexToInt(const char* s) const {
        uint8_t ret = 0xFF;
        if (*s >= '0' && *s <= '9') {
            ret = uint8_t(*s - '0');
        }
        else if (*s >= 'a' && *s <= 'f') {
            ret = uint8_t(*s - 'a' + 10);
        }
        else if (*s >= 'A' && *s <= 'F') {
            ret = uint8_t(*s - 'A' + 10);
        }
        return ret;
    }

    operator bool() const {
        return static_cast<bool>(UPPER | LOWER);
    }

    operator uint8_t() const {
        return static_cast<uint8_t>(LOWER);
    }

    operator uint16_t() const {
        return static_cast<uint16_t>(LOWER);
    }

    operator uint32_t() const {
        return static_cast<uint32_t>(LOWER);
    }

    operator uint64_t() const {
        return LOWER;
    }

    // Bitwise Operators


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator&(const T& rhs) const {
        return uint128_t(0, LOWER & (uint64_t)rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator&=(const T& rhs) {
        UPPER = 0;
        LOWER &= rhs;
        return *this;
    }

  

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator|(const T& rhs) const {
        return uint128_t(UPPER, LOWER | (uint64_t)rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator|=(const T& rhs) {
        LOWER |= (uint64_t)rhs;
        return *this;
    }



    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator^(const T& rhs) const {
        return uint128_t(UPPER, LOWER ^ (uint64_t)rhs);
    }



    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator^=(const T& rhs) {
        LOWER ^= (uint64_t)rhs;
        return *this;
    }


    // Bit Shift Operators


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator<<(const T& rhs) const {
        return *this << uint128_t(rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator<<=(const T& rhs) {
        *this = *this << uint128_t(rhs);
        return *this;
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator>>(const T& rhs) const {
        return *this >> uint128_t(rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator>>=(const T& rhs) {
        *this = *this >> uint128_t(rhs);
        return *this;
    }

    // Logical Operators

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    bool operator&&(const T& rhs) {
        return static_cast <bool> (*this && rhs);
    }

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    bool operator||(const T& rhs) {
        return static_cast <bool> (*this || rhs);
    }

    // Comparison Operators

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    bool operator==(const T& rhs) const {
        return (!UPPER && (LOWER == (uint64_t)rhs));
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    bool operator!=(const T& rhs) const {
        return (UPPER | (LOWER != (uint64_t)rhs));
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    bool operator>(const T& rhs) const {
        return (UPPER || (LOWER > (uint64_t)rhs));
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    bool operator<(const T& rhs) const {
        return (!UPPER) ? (LOWER < (uint64_t)rhs) : false;
    }

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    friend bool operator<(const T& lhs, const uint128_t& rhs) {
        if (rhs.upper()) {
            return true;
        }
        return ((uint64_t)lhs < rhs.lower());
    }

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    bool operator>=(const T& rhs) const {
        return ((*this > rhs) | (*this == rhs));
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    bool operator<=(const T& rhs) const {
        return ((*this < rhs) | (*this == rhs));
    }

    // Arithmetic Operators

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator+(const T& rhs) const {
        return uint128_t(UPPER + ((LOWER + (uint64_t)rhs) < LOWER), LOWER + (uint64_t)rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator+=(const T& rhs) {
        return *this += uint128_t(rhs);
    }



    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator-(const T& rhs) const {
        return uint128_t((uint64_t)(UPPER - ((LOWER - rhs) > LOWER)), (uint64_t)(LOWER - rhs));
    }



    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator-=(const T& rhs) {
        return *this = *this - uint128_t(rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator*(const T& rhs) const {
        return *this * uint128_t(rhs);
    }

    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator*=(const T& rhs) {
        return *this = *this * uint128_t(rhs);
    }


public:


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator/(const T& rhs) const {
        return *this / uint128_t(rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator/=(const T& rhs) {
        return *this = *this / uint128_t(rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t operator%(const T& rhs) const {
        return *this % uint128_t(rhs);
    }


    template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
    uint128_t& operator%=(const T& rhs) {
        return *this = *this % uint128_t(rhs);
    }



    uint128_t operator&(const uint128_t& rhs) const {
        return uint128_t(UPPER & rhs.UPPER, LOWER & rhs.LOWER);
    }

    uint128_t& operator&=(const uint128_t& rhs) {
        UPPER &= rhs.UPPER;
        LOWER &= rhs.LOWER;
        return *this;
    }

    uint128_t operator|(const uint128_t& rhs) const {
        return uint128_t(UPPER | rhs.UPPER, LOWER | rhs.LOWER);
    }

    uint128_t& operator|=(const uint128_t& rhs) {
        UPPER |= rhs.UPPER;
        LOWER |= rhs.LOWER;
        return *this;
    }

    uint128_t operator^(const uint128_t& rhs) const {
        return uint128_t(UPPER ^ rhs.UPPER, LOWER ^ rhs.LOWER);
    }

    uint128_t& operator^=(const uint128_t& rhs) {
        UPPER ^= rhs.UPPER;
        LOWER ^= rhs.LOWER;
        return *this;
    }

    uint128_t operator~() const {
        return uint128_t(~UPPER, ~LOWER);
    }

    uint128_t operator<<(const uint128_t& rhs) const {
        const uint64_t shift = rhs.LOWER;
        if (((bool)rhs.UPPER) || (shift >= 128)) {
            return 0;
        }
        else if (shift == 64) {
            return uint128_t(LOWER, 0);
        }
        else if (shift == 0) {
            return *this;
        }
        else if (shift < 64) {
            return uint128_t((UPPER << shift) + (LOWER >> (64 - shift)), LOWER << shift);
        }
        else if ((128 > shift) && (shift > 64)) {
            return uint128_t(LOWER << (shift - 64), 0);
        }
        else {
            return 0;
        }
    }

    uint128_t& operator<<=(const uint128_t& rhs) {
        *this = *this << rhs;
        return *this;
    }

    uint128_t operator>>(const uint128_t& rhs) const {
        const uint64_t shift = rhs.LOWER;
        if (((bool)rhs.UPPER) || (shift >= 128)) {
            return 0;
        }
        else if (shift == 64) {
            return uint128_t(0, UPPER);
        }
        else if (shift == 0) {
            return *this;
        }
        else if (shift < 64) {
            return uint128_t(UPPER >> shift, (UPPER << (64 - shift)) + (LOWER >> shift));
        }
        else if ((128 > shift) && (shift > 64)) {
            return uint128_t(0, (UPPER >> (shift - 64)));
        }
        else {
            return 0;
        }
    }

    uint128_t& operator>>=(const uint128_t& rhs) {
        *this = *this >> rhs;
        return *this;
    }

    bool operator!() const {
        return !(bool)(UPPER | LOWER);
    }

    bool operator&&(const uint128_t& rhs) const {
        return ((bool)*this && rhs);
    }

    bool operator||(const uint128_t& rhs) const {
        return ((bool)*this || rhs);
    }

    bool operator==(const uint128_t& rhs) const {
        return ((UPPER == rhs.UPPER) && (LOWER == rhs.LOWER));
    }

    bool operator!=(const uint128_t& rhs) const {
        return ((UPPER != rhs.UPPER) | (LOWER != rhs.LOWER));
    }

    bool operator>(const uint128_t& rhs) const {
        if (UPPER == rhs.UPPER) {
            return (LOWER > rhs.LOWER);
        }
        return (UPPER > rhs.UPPER);
    }

    bool operator<(const uint128_t& rhs) const {
        if (UPPER == rhs.UPPER) {
            return (LOWER < rhs.LOWER);
        }
        return (UPPER < rhs.UPPER);
    }

    bool operator>=(const uint128_t& rhs) const {
        return ((*this > rhs) | (*this == rhs));
    }

    bool operator<=(const uint128_t& rhs) const {
        return ((*this < rhs) | (*this == rhs));
    }

    uint128_t operator+(const uint128_t& rhs) const {
        return uint128_t(UPPER + rhs.UPPER + ((LOWER + rhs.LOWER) < LOWER), LOWER + rhs.LOWER);
    }

    uint128_t& operator+=(const uint128_t& rhs) {
        UPPER += rhs.UPPER + ((LOWER + rhs.LOWER) < LOWER);
        LOWER += rhs.LOWER;
        return *this;
    }

    uint128_t operator-(const uint128_t& rhs) const {
        return uint128_t(UPPER - rhs.UPPER - ((LOWER - rhs.LOWER) > LOWER), LOWER - rhs.LOWER);
    }

    uint128_t& operator-=(const uint128_t& rhs) {
        *this = *this - rhs;
        return *this;
    }

    uint128_t operator*(const uint128_t& rhs) const {
        // split values into 4 32-bit parts
        uint64_t top[4] = { UPPER >> 32, UPPER & 0xffffffff, LOWER >> 32, LOWER & 0xffffffff };
        uint64_t bottom[4] = { rhs.UPPER >> 32, rhs.UPPER & 0xffffffff, rhs.LOWER >> 32, rhs.LOWER & 0xffffffff };
        uint64_t products[4][4];

        // multiply each component of the values
        for (int y = 3; y > -1; y--) {
            for (int x = 3; x > -1; x--) {
                products[3 - x][y] = top[x] * bottom[y];
            }
        }

        // first row
        uint64_t fourth32 = (products[0][3] & 0xffffffff);
        uint64_t third32 = (products[0][2] & 0xffffffff) + (products[0][3] >> 32);
        uint64_t second32 = (products[0][1] & 0xffffffff) + (products[0][2] >> 32);
        uint64_t first32 = (products[0][0] & 0xffffffff) + (products[0][1] >> 32);

        // second row
        third32 += (products[1][3] & 0xffffffff);
        second32 += (products[1][2] & 0xffffffff) + (products[1][3] >> 32);
        first32 += (products[1][1] & 0xffffffff) + (products[1][2] >> 32);

        // third row
        second32 += (products[2][3] & 0xffffffff);
        first32 += (products[2][2] & 0xffffffff) + (products[2][3] >> 32);

        // fourth row
        first32 += (products[3][3] & 0xffffffff);

        // move carry to next digit
        third32 += fourth32 >> 32;
        second32 += third32 >> 32;
        first32 += second32 >> 32;

        // remove carry from current digit
        fourth32 &= 0xffffffff;
        third32 &= 0xffffffff;
        second32 &= 0xffffffff;
        first32 &= 0xffffffff;

        // combine components
        return uint128_t((first32 << 32) | second32, (third32 << 32) | fourth32);
    }

    uint128_t& operator*=(const uint128_t& rhs) {
        *this = *this * rhs;
        return *this;
    }

    void ConvertToVector(std::vector<uint8_t>& ret, const uint64_t& val) const {
        ret.push_back(static_cast<uint8_t>(val >> 56));
        ret.push_back(static_cast<uint8_t>(val >> 48));
        ret.push_back(static_cast<uint8_t>(val >> 40));
        ret.push_back(static_cast<uint8_t>(val >> 32));
        ret.push_back(static_cast<uint8_t>(val >> 24));
        ret.push_back(static_cast<uint8_t>(val >> 16));
        ret.push_back(static_cast<uint8_t>(val >> 8));
        ret.push_back(static_cast<uint8_t>(val));
    }

    void export_bits(std::vector<uint8_t>& ret) const {
        ConvertToVector(ret, const_cast<const uint64_t&>(UPPER));
        ConvertToVector(ret, const_cast<const uint64_t&>(LOWER));
    }

    std::pair <uint128_t, uint128_t> divmod(const uint128_t& lhs, const uint128_t& rhs) const {
        // Save some calculations /////////////////////
        if (rhs == uint128_t(0)) {
            throw std::domain_error("Error: division or modulus by 0");
        }
        else if (rhs == uint128_t(1)) {
            return std::pair <uint128_t, uint128_t>(lhs, 0);
        }
        else if (lhs == rhs) {
            return std::pair <uint128_t, uint128_t>(1, 0);
        }
        else if ((lhs == uint128_t(0)) || (lhs < rhs)) {
            return std::pair <uint128_t, uint128_t>(uint128_t(0), lhs);
        }

        std::pair <uint128_t, uint128_t> qr(uint128_t(0), uint128_t(0));
        for (uint8_t x = lhs.bits(); x > 0; x--) {
            qr.first <<= 1;
            qr.second <<= uint128_t(1);

            if ((lhs >> (x - 1U)) & uint128_t(1)) {
                ++qr.second;
            }

            if (qr.second >= rhs) {
                qr.second -= rhs;
                ++qr.first;
            }
        }
        return qr;
    }

    uint128_t operator/(const uint128_t& rhs) const {
        return divmod(*this, rhs).first;
    }

    uint128_t& operator/=(const uint128_t& rhs) {
        *this = *this / rhs;
        return *this;
    }

    uint128_t operator%(const uint128_t& rhs) const {
        return divmod(*this, rhs).second;
    }

    uint128_t& operator%=(const uint128_t& rhs) {
        *this = *this % rhs;
        return *this;
    }

    uint128_t& operator++() {
        return *this += uint128_t(1);
    }

    uint128_t operator++(int) {
        uint128_t temp(*this);
        ++* this;
        return temp;
    }

    uint128_t& operator--() {
        return *this -= uint128_t(1);
    }

    uint128_t operator--(int) {
        uint128_t temp(*this);
        --* this;
        return temp;
    }

    uint128_t operator+() const {
        return *this;
    }

    uint128_t operator-() const {
        return ~*this + uint128_t(1);
    }

    uint64_t upper() const {
        return UPPER;
    }

    uint64_t lower() const {
        return LOWER;
    }

    uint8_t bits() const {
        uint8_t out = 0;
        if (UPPER) {
            out = 64;
            uint64_t up = UPPER;
            while (up) {
                up >>= 1;
                out++;
            }
        }
        else {
            uint64_t low = LOWER;
            while (low) {
                low >>= 1;
                out++;
            }
        }
        return out;
    }

    std::string str(uint8_t base, const unsigned int& len) const {
        if ((base < 2) || (base > 16)) {
            throw std::invalid_argument("Base must be in the range [2, 16]");
        }
        std::string out = "";
        if (!(*this)) {
            out = "0";
        }
        else {
            std::pair <uint128_t, uint128_t> qr(*this, uint128_t(0));
            do {
                qr = divmod(qr.first, base);
                out = "0123456789abcdef"[static_cast<uint8_t>(qr.second)] + out;
            } while (qr.first);
        }
        if (out.size() < len) {
            out = std::string(len - out.size(), '0') + out;
        }
        return out;
    }

    friend std::ostream& operator<<(std::ostream& stream, const uint128_t& rhs) {
        if (stream.flags() & stream.oct) {
            stream << rhs.str(8,40);
        }
        else if (stream.flags() & stream.dec) {
            stream << rhs.str(10,40);
        }
        else if (stream.flags() & stream.hex) {
            stream << rhs.str(16,40);
        }
        return stream;
    }
};

#endif

namespace std {
    template<> class numeric_limits<uint128_t> {
    public:
        static uint128_t max() { return uint128_t(numeric_limits<uint64_t>::max(), numeric_limits<uint64_t>::max()); };
    };

   /* uint64_t operator&(uint128_t& a, uint64_t& b) {
        return a.lower() & b;
    }*/
}