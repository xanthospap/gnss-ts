///
/// @file genflags.hpp
///
/// @brief A generic wrapper class to provide flag-like functionality for a
///        (strongly-typed) enumeration type.
///
#ifndef __GENERIC_NGPT_FLAGS__
#define __GENERIC_NGPT_FLAGS__

#include <type_traits>
#include <vector>

namespace ngpt
{

/// Wrapper class around an enum to make it act like a flag class.
/// The class only provides basic functionality like, construction of flags,
/// set, check and clear flags.
///
/// @tparam FlagEnum A (strongly typed) enumeration class
///
template<typename FlagEnum>
    class flag
{
public:
    /// ft is the underlying enum class type (normaly int or char).
    using ft = typename std::underlying_type<FlagEnum>::type;

private:
    ft _f; ///< the flag as underlying type (ft).

public:
    
    /// Default constructor
    flag() noexcept : _f{} {}
    
    /// Constructor from a flag enumeration type (i.e. from a enumeration of the
    /// FlagEnum class).
    explicit flag(FlagEnum f) noexcept : _f{static_cast<ft>(f)} {}
    
    /// Constructor from multiple flags (via an std::initializer_list).
    explicit flag(std::initializer_list<FlagEnum>&& fs)
    : _f{}
    { for (auto& f : fs) { _f |= static_cast<ft>(f); } }
    
    /// Set a FlagEnum on.
    void set(FlagEnum f) noexcept { _f |= static_cast<ft>(f); }
    
    /// Clear a FlagEnum.
    void clear(FlagEnum f) noexcept { _f &= ~(static_cast<ft>(f)); }
    
    /// Clear the instance from all FlagEnums.
    void clear() noexcept { _f = static_cast<ft>(0); }
    
    /// Check if a FlagEnum is set (i.e. on).
    bool check(FlagEnum f) const noexcept { return _f & static_cast<ft>(f); }

    /// Check if a flag is clean (nothing is set).
    bool is_clean() const noexcept { return !static_cast<ft>(_f); }

    /// Equality operator
    bool operator==(flag f) const noexcept { return _f == f._f; }

    /// InEquality operator
    bool operator!=(flag f) const noexcept { return !(this->operator==(f)); }

}; // class flag<FlagEnum>

} // ngpt

#endif
