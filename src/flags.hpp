#ifndef __TS_FLAGS__
#define __TS_FLAGS__

#include <cstdint>
#include <climits>
#ifdef DEBUG
    #include <iostream>
    #include <bitset>
    #include <vector>
#endif

namespace ts {

/// Type to be used for record flags.
using flag_type = uint8_t;

/// Size (in bits) for a flag instance.
constexpr std::size_t flag_bit_sz { CHAR_BIT * sizeof(flag_type) };

#ifdef DEBUG
const char* abc = "abcdefghijklmnopqrstuvwxyz";
void bs_print(flag_type mt)
{
    std::bitset<flag_bit_sz> x(mt);

    std::cout << "\n|";
    for (int i=flag_bit_sz-1; i>=0; --i) { 
        std::cout << *(abc+i); 
    }
    std::cout << "|\n";

    std::cout << " " << x;
    std::cout << " (" << (unsigned)mt << ")";
}
#endif

/// Enumeration to hold all possible flags, a record can have.
enum class ts_flag : flag_type {
    outlier = 0b00000001,
    skip    = 0b00000010,
    ethq    = 0b00000100,
    velchg  = 0b00001000,
    jump    = 0b00010000
    /*
    f = 0b00100000,
    g = 0b01000000,
    h = 0b10000000,
    */
};

/// Underlying type of ts_flag enumeration.
using UT = typename std::underlying_type<ts_flag>::type;

/// Check if a flag is set.
inline
flag_type
check_flag(flag_type& lhs, ts_flag rhs)
noexcept
{
    return ( lhs >> static_cast<UT>(rhs) ) & 1;
}

/// Overload operator bitwise 'OR'.
inline 
flag_type& operator|=(flag_type& lhs, ts_flag rhs)
noexcept
{
    return lhs |= static_cast<UT>(rhs);
}

/// Overload operator bitwise 'OR'.
inline 
flag_type operator|(ts_flag lhs, ts_flag rhs)
noexcept
{
    return static_cast<UT>(lhs) | static_cast<UT>(rhs);
}

/// Overload operator bitwise 'XOR'.
/// Usage:
inline 
flag_type& operator^=(flag_type& lhs, ts_flag rhs)
noexcept
{
    return lhs ^= static_cast<UT>(rhs);
}

/// Overload operator bitwise 'AND'.
///// Usage:
inline
flag_type& operator&=(flag_type& lhs, ts_flag rhs)
noexcept
{
    return lhs &= static_cast<UT>(rhs);
}

/// Overload operator bitwise complement.
/// Usage:
inline 
flag_type operator~(ts_flag lhs)
noexcept
{
    return ~static_cast<UT>(lhs);
}

} // end namespace

#endif
