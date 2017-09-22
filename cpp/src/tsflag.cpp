#include "tsflag.hpp"

using ngpt::ts_event;
using ngpt::pt_marker;

/// Convert a ts_event to its identifing character.
/// \param[in] event A ts-event enumeration
/// \return          A char, identifying the ts-event. This can be:
///                  - 'j' for a jump
///                  - 'e' for an earthquake
///                  - 'v' for a velocity change
/// \warning Make sure all possible enum types (of ts_event class) are treated
///          here.
char
ngpt::event2char(ts_event event) noexcept
{
    switch (event)
    {
        case ts_event::jump:            return 'j';
        case ts_event::earthquake:      return 'e';
        case ts_event::velocity_change: return 'v';
        default: return 'x';
    }
}

/// For any enumeration type that can be wrapped around the flag (template)
/// class, there should be a function called 'skip' that determines if a data
/// point with a certain flag should be ignored.
bool
ngpt::__skip__(ngpt::flag<pt_marker> p) noexcept
{
    return /* p.check(pt_marker::outlier) || p.check(pt_marker::skip);
              or more simply ... */
           !p.is_clean();
}
