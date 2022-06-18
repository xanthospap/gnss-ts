#include "earthquake.hpp"
#include "datetime/datetime_read.hpp"
#include "datetime/datetime_write.hpp"
#include "geodesy/units.hpp"
#include <cstring>
#include <fstream>
#include <stdexcept>

/// @brief Format a datetime<T> instance based on the NOA catalogue files.
///
/// Given a datetime<T> instance, format it as a string of type:
/// YYYY OOO DD  HH MM SS.S
/// Where 'OOO' is the 3 first chars of the month, in uppercase.
/// @tparam    T  datetime instance resolution
/// @param[in] t  datetime<T> instance to be transformed to string
/// @return       string; the input datetime instance formated as:
///               'YYYY OOO DD HH MM SS.S'
std::string strfdt_as_noa(const dso::datetime<dso::milliseconds> &t) {
  using dso::_d2s_;
  using dso::_i2s_;

  auto ymd = t.as_ymd();
  auto hmsf = t.as_hmsf();

  double secs = std::get<2>(hmsf).as_underlying_type() +
                std::get<3>(hmsf) / dso::milliseconds::sec_factor<double>();

  const std::string wspace_str(1, ' ');

  return _i2s_((ymd.__year).as_underlying_type(), 4) + wspace_str +
         std::string((ymd.__month).short_name()) + wspace_str +
         _i2s_((ymd.__dom).as_underlying_type(), 2) + wspace_str +
         _i2s_(std::get<0>(hmsf).as_underlying_type(), 2) + wspace_str +
         _i2s_(std::get<1>(hmsf).as_underlying_type(), 2) + wspace_str +
         _d2s_(secs, 1);
}

std::string dso::earthquake::to_string() const noexcept {
  using dso::_d2s_;
  using dso::_i2s_;

  std::string evnt_str = strfdt_as_noa(m_epoch);
  evnt_str += "   " + _d2s_(dso::rad2deg<double>(m_lat), 2);
  evnt_str += "   " + _d2s_(dso::rad2deg<double>(m_lon), 2);
  evnt_str += "   " + _i2s_((static_cast<int>(m_depth / 1e3)), 3);
  evnt_str += "   " + _d2s_(m_magnitude, 1);
  return evnt_str;
}

/// @brief Format a datetime<T> instance based on the NOA catalogue files.
///
/// Given a datetime<T> instance, format it as a string of type:
/// YYYY OOO DD  HH MM SS.S
/// Where 'OOO' is the 3 first chars of the month, in uppercase.
/// @tparam    T  datetime instance resolution
/// @param[in] t  datetime<T> instance to be transformed to string
/// @return       string; the input datetime instance formated as:
///               'YYYY OOO DD HH MM SS.S'
std::string dso::strfdt_as_noa(const dso::datetime<dso::milliseconds> &t) {
  using dso::_d2s_;
  using dso::_i2s_;

  auto ymd = t.as_ymd();
  auto hmsf = t.as_hmsf();

  double secs = std::get<2>(hmsf).as_underlying_type() +
                std::get<3>(hmsf) / dso::milliseconds::sec_factor<double>();

  const std::string wspace_str(1, ' ');

  return _i2s_((ymd.__year).as_underlying_type(), 4) + wspace_str +
         std::string((ymd.__month).short_name()) + wspace_str +
         _i2s_((ymd.__dom).as_underlying_type(), 2) + wspace_str +
         _i2s_(std::get<0>(hmsf).as_underlying_type(), 2) + wspace_str +
         _i2s_(std::get<1>(hmsf).as_underlying_type(), 2) + wspace_str +
         _d2s_(secs, 1);
}
