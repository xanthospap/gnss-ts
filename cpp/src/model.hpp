#ifndef __NGPT_TS_MODEL_HPP__
#define __NGPT_TS_MODEL_HPP__

#include <iostream>
#ifdef DEBUG
#include <cfenv>
#include <cstdlib>
#endif
#include "eigen3/Eigen/Core"
#include "event_list.hpp"
#include "ggdatetime/dtcalendar.hpp"
#include "md_harmonics.hpp"
#include "md_jump.hpp"
#Include "md_velocity_change.hpp"
// #include "psd.hpp"

namespace ngpt {

/// This class represents a (time-series) model.
///
/// The purpose of this class is to describe a (mathematical) model of a
/// time-series. So it should serve several requirements, e.g.
/// - given an epoch (i.e. datetime), compute the value of the model,
/// - be used within Least-Squares modules to perform adjustments (i.e. fit)
/// The model can consist of any number of md_jump, md_harmonics and
/// md_velocity_change, plus a simple linear model (i.e. x0 and Vx0 parameters);
/// so if no md_jump, md_harmonics or md_velocity_change are added, then the
/// model defaults to a linear one.
/// Every instance also hold am epoch (i.e. datetime<ngpt::milliseconds>) member
/// variable, which is the central epoch of computation/estimation.
class ts_model {
public:
  /// Constructor; this default to a linear model.
  ts_model() noexcept : m_x0{0e0}, m_vx{0e0} {};

  /// Constructor using an event_list instance.
  /// @param[in] events An instance of type event_list; all events recorded
  ///                   in this instance, are going to be included in the
  ///                   constructed model. This includes:
  ///                   - jumps
  ///                   - velocity changes
  ///                   - earthquakes
  explicit ts_model(const event_list &events) noexcept;

  /// This function adds new md_harmonics instances to this (model)
  /// instance. The new md_harmonics added, are constructed using only
  /// their period (i.e. all other parametrs are default-constructed).
  /// @param[in] periods A vector<double> containing period values for the
  ///                    harmonics to be added.
  /// @note If you wish to add a harmonic signal in the model, using more info
  /// (e.g. include in- and out-of-phase amplitudes, see the function
  /// add_period).
  void add_periods(const std::vector<double> &periods) noexcept {
    for (auto i : periods)
      m_harmonics.emplace_back(i);
  }

  /// This function adds a new md_harmonics instance to this (model) instance.
  /// The new harmonic to be added, can be fully specified.
  /// @param[in] period    The period of the harmonic.
  /// @param[in] in_phase  Amplitude of the in-phase component.
  /// @param[in] out_phase Amplitude of the out-of-phase component.
  /// @param[in] start     Start of the validity interval (default value
  ///                      datetime<ngpt::milliseconds>::min()).
  /// @param[in] stop      End of the validity interval (default value
  ///                      datetime<ngpt::milliseconds>::max()).
  void add_period(
      double period, double in_phase = 0e0, double out_of_phase = 0e0,
      datetime<ngpt::milliseconds> start = datetime<ngpt::milliseconds>::min(),
      datetime<ngpt::milliseconds> stop =
          datetime<ngpt::milliseconds>::max()) noexcept {
    m_harmonics.emplace_back(period, start, stop, in_phase, out_of_phase);
  }

  /// Add a jump (i.e. offset) in this model instance. This function will
  /// create a new md_jump instance based on the input parametrs and add it
  /// to this md_model.
  /// @param[in] at  The epoch the jump happened.
  /// @param[in] val The value of the offset (default is 0).
  void add_jump(datetime<ngpt::milliseconds> at, double val = 0e0) noexcept {
    m_jumps.emplace_back(at, val);
  }

  /// Add a velocity change in this model instance. This function will
  /// create a new md_velocity_change instance based on the input parametrs
  /// and add it to this md_model.
  /// @param[in] start   The starting epoch of this velocity change.
  /// @param[in] val     The magnitude of the velocity change.
  /// @param[in] stop    The end epoch of this velocity change; default is
  ///                    ngpt::datetime<ngpt::milliseconds>::max(), i.e. this
  ///                    velocity change starts at epoch start and never stops.
  void
  add_velocity_change(datetime<ngpt::milliseconds> start, double val = 0e0,
                      ngpt::datetime<ngpt::milliseconds> stop =
                          ngpt::datetime<ngpt::milliseconds>::max()) noexcept {
    m_vel_changes.emplace_back(start, stop, val);
  }

  void add_earthquake(datetime<ngpt::milliseconds> start, psd_model m,
                      double a1 = 1e-3, double t1 = 1e-1, double a2 = 1e-3,
                      double t2 = 1e-1) noexcept {
    m_earthqs.emplace_back(start, m, a1, t1, a2, t2);
  }

  void add_earthquake(const md_earthquake<T> &e) noexcept {
    m_earthqs.emplace_back(e);
  }

  /// Get the mean epoch of the model (aka the central epoch of computation).
  /// Const version.
  datetime<ngpt::milliseconds> mean_epoch() const noexcept {
    return m_mean_epoch;
  }

  /// Get the mean epoch of the model (akak the central epoch of computation).
  /// Non-Const version.
  datetime<ngpt::milliseconds> &mean_epoch() noexcept { return m_mean_epoch; }

  /// Get the constant linear term (const version).
  double x0() const noexcept { return m_x0; }

  /// Get the constant linear term (non-const version).
  double &x0() noexcept { return m_x0; }

  /// Get the velocity (const version).
  double vx() const noexcept { return m_vx; }

  /// Get the velocity (non-const version).
  double &vx() noexcept { return m_vx; }

private:
  double m_x0, ///< const linear term.
      m_vx,    ///< linear velocity
      m_x0_stddev, m_vx_stddev;
  std::vector<md_jump> m_jumps;                  ///< vector of jumps
  std::vector<md_harmonics> m_harmonics;         ///< vector of harmonics
  std::vector<md_velocity_change> m_vel_changes; ///< vector of velocity changes
  std::vector<md_earthquake> m_earthqs;          ///< vector of earthquakes
  datetime<ngpt::milliseconds>
      m_mean_epoch; ///< mean epoch (aka central computation epoch)

}; // ts_model

} // namespace ngpt

#endif
