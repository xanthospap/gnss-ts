#ifndef __NGPT_VELCHNG_MODEL_HPP__
#define __NGPT_VELCHNG_MODEL_HPP__

namespace dso {

/// A class to represent a velocity change (in a time-series). The class is
/// dead simple, just holding a start and stop date (i.e. the validity interval
/// of this velocity change) and the magnitude of the velocity (change).
class md_velocity_change {
public:
  /// Constructor. The stop date defaults to
  /// dso::datetime<dso::milliseconds>::max().
  /// @param[in] start   The starting epoch of this velocity change.
  /// @param[in] stop    The end epoch of this velocity change; default is
  ///                    dso::datetime<dso::milliseconds>::max(), i.e. this
  ///                    velocity change starts at epoch start and never stops.
  /// @param[in] new_vel The magnitude of the velocity change.
  /// @param[in] new_vel_stddev The std. deviation of the magnitude of the
  ///                    velocity change.
  md_velocity_change(dso::datetime<dso::milliseconds> start,
                     dso::datetime<dso::milliseconds> stop =
                         dso::datetime<dso::milliseconds>::max(),
                     double new_vel = 0e0, double new_vel_stddev = 0e0) noexcept
      : m_start{start}, m_stop{stop}, m_newvel{new_vel}, m_newvel_stddev{
                                                             new_vel_stddev} {};

  /// Get the start date (start of validity interval), const version.
  dso::datetime<dso::milliseconds> start() const noexcept { return m_start; }

  /// Get the stop date (end of validity interval), const version.
  dso::datetime<dso::milliseconds> stop() const noexcept { return m_stop; }

  /// Get the velocity change magnitude (const version).
  double value() const noexcept { return m_newvel; }

  /// Get the velocity change magnitude (non-const version).
  double &value() noexcept { return m_newvel; }

  /// Get the velocity change std. deviation (const version).
  double stddev() const noexcept { return m_newvel_stddev; }

  /// Get the velocity change std. deviation (non-const version).
  double &stddev() noexcept { return m_newvel_stddev; }

private:
  dso::datetime<dso::milliseconds> m_start, ///< Start of validity niterval.
      m_stop;                                 ///< End of validity interval.
  double m_newvel,                            ///< Magnitude of velocity change.
      m_newvel_stddev; ///< Velocity change std. deviation.

}; // md_velocity_change

} // namespace dso
#endif
