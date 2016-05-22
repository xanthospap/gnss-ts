#ifndef __GENERIC_NGPT_FLAGS__
#define __GENERIC_NGPT_FLAGS__

#include <type_traits>

namespace ngpt {

/*
 * ///////////////////////////////////////////////////////////////////////////
 *  A generic way to implement the above for whatever enum class (i.e. strongly
 *  typed). This works for C++14.
 * ///////////////////////////////////////////////////////////////////////////
 */

/* wrapper class around an enum to make it act like a flag class */
template<typename FlagEnum>
    class flag
{
public:
    using ft = typename std::underlying_type<FlagEnum>::type;

private:
    ft _f;

public:
    
    flag() noexcept : _f{} {}
    
    explicit flag(FlagEnum f) noexcept : _f{static_cast<ft>(f)} {}
    
    explicit flag(std::initializer_list<FlagEnum>&& fs) : _f{}
    {
        for (auto& f : fs) { _f |= static_cast<ft>(f); }
    }
    
    /// Set a flag
    void set(FlagEnum f) noexcept { _f |= static_cast<ft>(f); }
    
    /// Clear a flag
    void clear(FlagEnum f) noexcept { _f &= ~(static_cast<ft>(f)); }
    
    /// Check if a flag is set
    bool check(FlagEnum f) const noexcept { return _f & static_cast<ft>(f); }

}; // class flag<FlagEnum>

} // ngpt

#endif
