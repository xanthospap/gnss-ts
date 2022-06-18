#ifndef __NGPT_PSD_HPP__
#define __NGPT_PSD_HPP__

#ifdef DEBUG
#include <cfenv>
#include <cstdlib>
#endif
#include "event_list.hpp"
#include "datetime/dtcalendar.hpp"

namespace dso {

enum class psd_model : int {
  pwl,    ///< Piece-Wise Linear
  log,    ///< Logarithmic
  exp,    ///< Exponential
  logexp, ///< Logarithmic and Exponential
  expexp  ///< Two Exponential functions
};

int _psd2int_(psd_model x) noexcept {
  switch (x) {
  case psd_model::pwl:
    return 0;
  case psd_model::log:
    return 1;
  case psd_model::exp:
    return 2;
  case psd_model::logexp:
    return 3;
  case psd_model::expexp:
    return 4;
  default:
    return -1;
  }
}

class md_earthquake {
public:
  explicit md_earthquake(dso::datetime<dso::milliseconds> t, psd_model md,
                         double a1 = 0e0, double t1 = 1.0e0, double a2 = 0e0,
                         double t2 = 1.0e0) noexcept
      : m_model{md}, m_start{t}, m_a1{a1}, m_t1{t1}, m_a2{a2}, m_t2{t2} {}

  void set_psd_type(psd_model md) noexcept { m_model = md; }

  psd_model psd_type() const noexcept { return m_model; }

  double &a1() noexcept { return m_a1; }

  double &a2() noexcept { return m_a2; }

  double &t1() noexcept { return m_t1; }

  double &t2() noexcept { return m_t2; }

  dso::datetime<dso::milliseconds> &start() noexcept { return m_start; }

  double a1() const noexcept { return m_a1; }

  double a2() const noexcept { return m_a2; }

  double t1() const noexcept { return m_t1; }

  double t2() const noexcept { return m_t2; }

  dso::datetime<dso::milliseconds> start() const noexcept { return m_start; }

  std::size_t parameters() const noexcept {
    switch (m_model) {
    case psd_model::pwl:
      return 1;
    case psd_model::log:
      return 2;
    case psd_model::exp:
      return 2;
    case psd_model::logexp:
      return 4;
    case psd_model::expexp:
      return 4;
    }
    return 0;
  }

  double offset() const noexcept {
    switch (m_model) {
    case psd_model::pwl:
      return m_a1;
    case psd_model::log:
      return m_a1;
    case psd_model::exp:
      return m_a1;
    case psd_model::logexp:
      return m_a1 + m_a2;
    case psd_model::expexp:
      return m_a1 + m_a2;
    }
    return 0;
  }

  double abs_offset() const noexcept {
    switch (m_model) {
    case psd_model::pwl:
      return std::abs(m_a1);
    case psd_model::log:
      return std::abs(m_a1);
    case psd_model::exp:
      return std::abs(m_a1);
    case psd_model::logexp:
      return std::sqrt(m_a1 * m_a1 + m_a2 * m_a2);
    case psd_model::expexp:
      return std::sqrt(m_a1 * m_a1 + m_a2 * m_a2);
    }
    return 0;
  }

  double value_at(const dso::datetime<dso::milliseconds> &t) const {
    if (t <= this->start())
      return 0e0;

    auto dt_intvr = dso::delta_date(t, start());
    double dtq = dt_intvr.as_mjd() / 365.25e0;
#ifdef DEBUG
    double dtq_conf = (t.as_mjd() - start().as_mjd()) / 365.25e0;
    assert(std::abs(dtq - dtq_conf) < 1e-10);
#endif
    double d{0e0}, te1, te2;

    switch (m_model) {
    case psd_model::pwl:
      d = m_a1;
      break;
    case psd_model::log:
      te1 = dtq / m_t1;
      d = m_a1 * std::log10(1e0 + te1);
      break;
    case psd_model::exp:
      te1 = dtq / m_t1;
      d = -m_a1 * std::expm1(-te1);
      break;
    case psd_model::logexp:
      te1 = dtq / m_t1;
      te2 = dtq / m_t2;
      d = m_a1 * std::log(1e0 + te1) + m_a2 * (1e0 - std::exp(-te2));
      break;
    case psd_model::expexp:
      te1 = dtq / m_t1;
      te2 = dtq / m_t2;
      d = m_a1 * (1e0 - std::exp(-te1)) + m_a2 * (1e0 - std::exp(-te2));
    }
    return d;
  }

  void diriv_at(const dso::datetime<dso::milliseconds> &t, double &da1,
                double &dt1, double &da2, double &dt2) const {
    auto dt_intvr = dso::delta_date(t, this->start());
    double dtq = dt_intvr.as_mjd() / 365.25e0;
#ifdef DEBUG
    assert(dtq >= 0e0);
    double dtq_conf = (t.as_mjd() - start().as_mjd()) / 365.25e0;
    assert(std::abs(dtq - dtq_conf) < 1e-10);
#endif
    double te1, te2;

    switch (m_model) {
    case psd_model::pwl:
      da1 = 1e0;
      break;
    case psd_model::log:
      te1 = dtq / m_t1;
      da1 = std::log10(1e0 + te1);
      dt1 = -(m_a1 * dtq) / (m_t1 * m_t1 * (1e0 + te1));
      break;
    case psd_model::exp:
#ifdef DEBUG
      assert(std::abs(te1) < 708.0e0);
#endif
      te1 = dtq / m_t1;
      da1 = -std::expm1(-te1);
      dt1 = -m_a1 * std::exp(-te1) * dtq / (m_t1 * m_t1);
      break;
    case psd_model::logexp:
      te1 = dtq / m_t1;
      te2 = dtq / m_t2;
      da1 = std::log(1e0 + te1);
      dt1 = m_a1 * (-te1 / m_t1) / (1e0 + te1);
      da2 = 1e0 - std::exp(-te2);
      dt2 = m_a2 * std::exp(-te2) * (-te2) / m_t2;
      break;
    case psd_model::expexp:
      te1 = dtq / m_t1;
      te2 = dtq / m_t2;
      da1 = 1e0 - std::exp(-te1);
      dt1 = m_a1 * std::exp(-te1) * (-te1) / m_t1;
      da2 = 1e0 - std::exp(-te2);
      dt2 = m_a2 * std::exp(-te2) * (-te2) / m_t2;
    }
    return;
  }

  auto to_earthquake(const dso::event_list<T> &lst) const {
    auto evn = lst.event_at(m_start);
    assert(evn.event_type() == dso::ts_event::earthquake);
    // return event_string2earthquake<T>(evn.info_str());
    return resolve_noa_earthquake_line<T>(evn.info_str().c_str());
  }

  double magnitude(const dso::event_list<T> &lst) const {
    auto ethq{to_earthquake(lst)};
    return ethq.magnitude();
  }

private:
  psd_model m_model;
  dso::datetime<dso::milliseconds> m_start;
  double m_a1, m_t1, m_a2, m_t2;
}; // md_earthquake

} // namespace dso

#endif
