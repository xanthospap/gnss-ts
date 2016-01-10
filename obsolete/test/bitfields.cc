#include <cstdint>
#include <climits>
#include <bitset>
#include <vector>
#include <iostream>

/*
 * Μου 'σπασε τα αρχίδια. Να μην το διαγράψω, δεν ξέρεις ποτέ.
 */

using mytype = uint8_t;

constexpr std::size_t numb { CHAR_BIT * sizeof(mytype) };

const char* abc = "abcdefghijklmnopqrstuvwxyz";

void bs_print(mytype mt)
{
    std::bitset<numb> x(mt);

    std::cout << "\n|";
    for (int i=numb-1; i>=0; --i) { std::cout << *(abc+i); }
    std::cout << "|\n";

    std::cout << " " << x;
    std::cout << " (" << (unsigned)mt << ")\n";
}

enum class my_fags : mytype {
    a = 0b00000001,
    b = 0b00000010,
    c = 0b00000100,
    d = 0b00001000,
    e = 0b00010000,
    f = 0b00100000,
    g = 0b01000000,
    h = 0b10000000,
};

using UT = typename std::underlying_type<my_fags>::type;

inline mytype& operator|=(mytype& lhs, my_fags rhs) noexcept
{
    return lhs |= static_cast<UT>(rhs);
}
inline mytype& operator^=(mytype& lhs, my_fags rhs) noexcept
{
    return lhs ^= static_cast<UT>(rhs);
}
inline mytype& operator&=(mytype& lhs, my_fags rhs) noexcept
{
    return lhs &= static_cast<UT>(rhs);
}
inline mytype operator~(my_fags lhs) noexcept
{
    return ~static_cast<UT>(lhs);
}

int main ()
{

    std::cout << "\nSize of given type is : " << numb << " bits.";
    std::cout << "\nSize of flag       is : " 
              << CHAR_BIT * sizeof(my_fags) << " bits.";
    
    mytype bf;
    std::cout << "\nDefault initialize";
    bs_print( bf );

    /* set the bits ... */
    std::cout << "\nIncrementaly adding bits";
    bf |= my_fags::a;
    bf |= my_fags::b;
    bf |= my_fags::c;
    bf |= my_fags::d;
    bf |= my_fags::e;
    bf |= my_fags::f;
    bf |= my_fags::g;
    bf |= my_fags::h;
    bs_print( bf );

    /* clear some bits ... */
    std::cout << "\nClearing bits e, f, g";
    bf &= ~( my_fags::e );
    bf &= ~( my_fags::f );
    bf &= ~( my_fags::g );
    bs_print( bf );

    /* toggle bits ... */
    std::cout << "\nToggling bits d and e";
    bf ^= my_fags::d;
    bf ^= my_fags::e;
    bs_print( bf );

    /*
    std::vector<mytype> nums = { 0b0000001,
                                 0b0000011,
                                 0b0000111,
                                 0b0001111,
                                 0b0011111,
                                 0b0111111 };

    for ( auto i : nums ) { bs_print( i ); std::cout << "\n"; }
    */

    std::cout << "\n";
    return 0;
}
