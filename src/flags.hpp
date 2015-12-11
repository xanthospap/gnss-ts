#ifndef __TS_FLAFS__
#define __TS_FLAGS__

#include <cstdint>
#include <climits>

#include <vector>
#ifdef DEBUG
    #include <iostream>
#endif

using bttp = uint8_t;

enum bp_flag {
    outlier  = /*0x01*/ /*1 << 0*/ 0b10000,
    skip     = /*0x02*/ /*1 << 1*/ 0b01000,
    erthq    = /*0x04*/ /*1 << 2*/ 0b00100,
    velchg   = /*0x08*/ /*1 << 3*/ 0b00010,
    jump     = /*0x10*/ /*1 << 4*/ 0b00001
};

static_assert(CHAR_BIT * sizeof(bttp) >= 8,
    "Can't store 8 bits in this shit!");

/*
using T = std::underlying_type_t<bp_flag>;

inline bp_flag& operator|=(bp_flag& lhs, bp_flag rhs)
{
    lhs = (bp_flag)(static_cast<T>(lhs) | static_cast<T>(rhs));
    return lhs;
}
*/

class flag {

private:
    inline bttp 
    _test_bit_(bp_flag n) 
    const noexcept
    {
        return bits_ & (1 << static_cast<bttp>(n));
    }

    inline bttp&
    _set_bit_(bp_flag n) 
    noexcept
    {
        return bits_ |= 1 << static_cast<bttp>(n);
    }
    
    inline bttp&
    _unset_bit_(bp_flag n) 
    noexcept
    {
        return bits_ &= ~(1 << static_cast<bttp>(n));
    }

public:
    flag() noexcept 
        : bits_{0x00} 
    {}
    
    bttp is_outlier() const noexcept
    {
        /*return _test_bit_(bp_flag::outlier);*/
        return bits_ & ( 1 << 0b1 );
    }
    bttp& set_outlier() noexcept
    {
        /*return _set_bit_(bp_flag::outlier);*/
        return bits_ |= 1 << 0b1;
    }
    bttp& unset_outlier() noexcept
    {
        /*return _unset_bit_(bp_flag::outlier);*/
        return bits_ &= ~(1 << 0b1);
    }

    bttp is_skip() const noexcept
    {
        return _test_bit_(bp_flag::skip);
    }
    bttp& set_skip() noexcept
    {
        return _set_bit_(bp_flag::skip);
    }
    bttp& unset_skip() noexcept
    {
        return _unset_bit_(bp_flag::skip);
    }
    
    bttp is_erthq() const noexcept
    {
        return _test_bit_(bp_flag::erthq);
    }
    bttp& set_erthq() noexcept
    {
        return _set_bit_(bp_flag::erthq);
    }
    bttp& unset_erthq() noexcept
    {
        return _unset_bit_(bp_flag::erthq);
    }
    
    bttp is_velchg() const noexcept
    {
        return _test_bit_(bp_flag::velchg);
    }
    bttp& set_velchg() noexcept
    {
        return _set_bit_(bp_flag::velchg);
    }
    bttp& unset_velchg() noexcept
    {
        return _unset_bit_(bp_flag::velchg);
    }
    
    bttp is_jump() const noexcept
    {
        return _test_bit_(bp_flag::jump);
    }
    bttp& set_jump() noexcept
    {
        return _set_bit_(bp_flag::jump);
    }
    bttp& unset_jump() noexcept
    {
        return _unset_bit_(bp_flag::jump);
    }

#ifdef DEBUG
    void debug_print() const
    {
        std::vector<bp_flag>  vf {  bp_flag::outlier, bp_flag::skip, 
            bp_flag::erthq, bp_flag::velchg, bp_flag::jump };

        for ( auto i : vf ) {
          if ( this->_test_bit_(i) ) {
            std::cout << "1";
          } else {
            std::cout << "0";
          }
        }
        std::cout << "(" << static_cast<unsigned>(this->bits_) << ")";
    }
#endif

private:
    bttp bits_;

};

#endif
