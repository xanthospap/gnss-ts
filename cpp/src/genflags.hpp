///
/// \file genflags.hpp
///
/// \brief A generic wrapper class to provide flag-like functionality for a
///        (strongly-typed) enumeration type.
///
#ifndef __GENERIC_NGPT_FLAGS__
#define __GENERIC_NGPT_FLAGS__

#include <type_traits>
#include <vector>

namespace ngpt
{

/// \brief Wrapper class around an enum to make it act like a flag class.
///
/// The class only provides basic functionality like, construction of flags,
/// set, check and clear flags.
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
    
    /// Constructor from multiple flags.
    explicit flag(std::initializer_list<FlagEnum>&& fs) : _f{}
    {
        for (auto& f : fs) { _f |= static_cast<ft>(f); }
    }
    
    /// Set a flag.
    void set(FlagEnum f) noexcept { _f |= static_cast<ft>(f); }
    
    /// Clear a flag.
    void clear(FlagEnum f) noexcept { _f &= ~(static_cast<ft>(f)); }
    
    /// Check if a flag is set.
    bool check(FlagEnum f) const noexcept { return _f & static_cast<ft>(f); }

}; // class flag<FlagEnum>

/// An enumeration type to hold possible flags for coordinate time-series.
enum class ts_event : char
{
    jump,
    earthquake,
    velocity_change,
    outlier,
    skip
};

/// Convert a ts_event to its identifing character.
char event2char(ts_event event) noexcept
{
    switch (event) {
        case ngpt::ts_event::jump: return 'j';
        case ngpt::ts_event::earthquake: return 'e';
        case ngpt::ts_event::velocity_change: return 'v';
        case ngpt::ts_event::outlier: return 'o';
        case ngpt::ts_event::skip: return 's';
        default: return 'x';
    }
}

/// For any enumerationtype that can be wrapped around the flag (template) class,
/// there should be a function called skip that determines if a data point with
/// a certain flag should be ignored.
bool
__skip__(flag<ts_event> p) noexcept
{
    return p.check(ts_event::outlier) || p.check(ts_event::skip);
}

/// For any enumerationtype that can be wrapped around the flag (template) class,
/// there should be an overload for the '<<' operator.
std::ostream&
operator<<(std::ostream& os, const flag<ts_event>& event)
{
    std::vector<ts_event> events = {
        ts_event::jump,
        ts_event::earthquake,
        ts_event::velocity_change,
        ts_event::outlier,
        ts_event::skip };

    for (auto i : events) if ( event.check(i) ) os << event2char(i);

    return os;
}

} // ngpt

#endif
