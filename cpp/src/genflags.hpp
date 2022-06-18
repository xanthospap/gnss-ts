#ifndef __GENERIC_NGPT_FLAGS_HPP__
#define __GENERIC_NGPT_FLAGS_HPP__

///
/// @file genflags.hpp
///
/// @brief A generic wrapper class to provide flag-like functionality for a
///        (strongly-typed) enumeration type.
///

#include <type_traits>
#include <vector>

namespace dso {

/// @class flag
///
/// @brief Wrapper class around an enum to make it act like a flag class.
///
/// This class is meant to act like a wrapper around a (strongly-typed) enum
/// type. That is, a user can declare an enum type and wrap it around this
/// class to make a flag-like class.
/// Within this project, this class is meant to act like a time-series data
/// point flag; aka, all data points in a time-series can have an instance
/// of this kind to mark their quality (outlier, ok, jump, etc, where outlier,
/// ok and jump are declared within the enumeration). Hence, a data point can
/// be flaged with none, one or more of these instances (it can be both an
/// outlier and a jump).
/// The class only provides basic functionality like, construction of flags,
/// set, check and clear flags.
/// Here is a simple example:
/// @code
/// enum class pt_marker : int { outlier = 1, skip = 2 };
/// flag<pt_marker> p{};
/// p.set(outlier);
/// @endcode
///
/// @see tsflag.hpp
///
/// @tparam FlagEnum A (strongly typed) enumeration class
template <typename FlagEnum> class flag {
public:
  /// ft is the underlying enum class type (normaly int or char).
  using ft = typename std::underlying_type<FlagEnum>::type;

private:
  ft _f; ///< the flag as underlying type (ft).

public:
  /// Default constructor (no FlagEnum is set).
  flag() noexcept : _f{} {}

  /// Constructor from a flag enumeration type (i.e. from a enumeration of
  /// the FlagEnum class).
  /// @param[in] f  A FlagEnum instance.
  explicit flag(FlagEnum f) noexcept : _f{static_cast<ft>(f)} {}

  /// Constructor from multiple flags (via an std::initializer_list).
  /// @param[in] fs  An std::initializer_list<FlagEnum> instance; all elements
  ///                of this list will be set on on the produced instance.
  explicit flag(std::initializer_list<FlagEnum> &&fs) : _f{} {
    for (auto &f : fs) {
      _f |= static_cast<ft>(f);
    }
  }

  /// Set a FlagEnum on.
  /// @param[in] f  A FlagEnum instance; this will be "switched on".
  void set(FlagEnum f) noexcept { _f |= static_cast<ft>(f); }

  /// Clear a FlagEnum.
  /// @param[in] f  Clear (i.e. "switch off") this FlagEnum.
  void clear(FlagEnum f) noexcept { _f &= ~(static_cast<ft>(f)); }

  /// Clear the instance from all FlagEnums. All FlagEnum are "switched off".
  void clear() noexcept { _f = static_cast<ft>(0); }

  /// Check if a FlagEnum is set (i.e. on).
  /// @param[in] f  Check if this FlagEnum is "switched on".
  /// @return true if FlagEnum f is set; false otherwise.
  bool is_set(FlagEnum f) const noexcept { return _f & static_cast<ft>(f); }
  /* TODO the following is obsolete */
  bool check(FlagEnum f) const noexcept { return this->is_set(f); }

  /// Check if a flag is clean (nothing is set).
  /// @return true if instance has no FlagEnum set; false otherwise.
  bool is_clean() const noexcept { return !static_cast<ft>(_f); }

  /// Equality operator (two instances have the same FlagEnum s on)
  /// @param[in] f Instance to check against
  /// @return true if the two instances have exactly the same FlagEnum s set;
  ///         false otherwise.
  bool operator==(flag f) const noexcept { return _f == f._f; }

  /// InEquality operator.
  /// @param[in] f Instance to check against
  /// @return false if the two instances have exactly the same FlagEnum s set;
  ///         else true.
  bool operator!=(flag f) const noexcept { return !(this->operator==(f)); }

}; // flag<FlagEnum>

} // namespace dso

#endif
