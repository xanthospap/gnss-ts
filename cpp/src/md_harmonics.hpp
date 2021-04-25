#ifndef __NGPT_HARMONICS_MODEL_HPP__
#define __NGPT_HARMONICS_MODEL_HPP__

#include "ggeodesy/geoconst.hpp"

namespace ngpt {

/// A class to represent a harmonic signal (in a time-series). The harmonic
/// signal is described by in and out-of-phase amplitudes, a period/frequency,
/// and the starting and ending times (i.e. its validity interval).
///@note      Angular frequency (i.e. ω) is given by: ω = 2πf = ω = 2π/T
class md_harmonics {
public:
  /// Constructor.
  /// @param[in] period    The period of the harmonic.
  /// @param[in] start     Start of the validity interval (default value
  ///                      datetime<T>::min()).
  /// @param[in] stop      End of the validity interval (default value
  ///                      datetime<T>::max()).
  /// @param[in] in_phase  Amplitude of the in-phase component.
  /// @param[in] out_phase Amplitude of the out-of-phase component.
  /// @param[in] in_phase_stddev  Std. deviation of the amplitude value/estimate
  ///                      of the in-phase component.
  /// @param[in] out_phase_stddev Std. deviation of the amplitude value/estimate
  ///                      of the out-of-phase component.
  md_harmonics(double period,
               ngpt::datetime<ngpt::milliseconds> start =
                   ngpt::datetime<ngpt::milliseconds>::min(),
               ngpt::datetime<ngpt::milliseconds> stop =
                   ngpt::datetime<ngpt::milliseconds>::max(),
               double in_phase = 0e0, double out_phase = 0e0,
               double in_phase_stddev = 0e0,
               double out_phase_stddev = 0e0) noexcept
      : m_start{start}, m_stop{stop}, m_afreq{ngpt::D2PI / period},
        m_in_phase{in_phase}, m_in_stddev{in_phase_stddev},
        m_out_phase{out_phase}, m_out_stddev{out_phase_stddev} {};

  /// Get the starting epoch (of the validity interval).
  ngpt::datetime<ngpt::milliseconds> start() const noexcept { return m_start; }

  /// Get the ending epoch (of the validity interval).
  ngpt::datetime<ngpt::milliseconds> stop() const noexcept { return m_stop; }

  /// Get the angular frequency (ω = 2πf).
  double angular_frequency() const noexcept { return m_afreq; }

  /// Get the period (i.e. T = 2π/ω).
  double period() const noexcept { return ngpt::D2PI / m_afreq; }

  /// Get the in-phase component amplitude (const version).
  double in_phase() const noexcept { return m_in_phase; }

  /// Get the out-of-phase component amplitude (const version).
  double out_of_phase() const noexcept { return m_out_phase; }

  /// Get the in-phase component amplitude (non-const version).
  double &in_phase() noexcept { return m_in_phase; }

  /// Get the out-of-phase component amplitude (non-const version).
  double &out_of_phase() noexcept { return m_out_phase; }

  /// Get the in-phase component std. deviation (const version).
  double in_phase_stddev() const noexcept { return m_in_stddev; }

  /// Get the out-of-phase component std. deviation (const version).
  double out_phase_stddev() const noexcept { return m_out_stddev; }

  /// Get the in-phase component std. deviation (non-const version).
  double &in_phase_stddev() noexcept { return m_in_stddev; }

  /// Get the out-of-phase component std. deviation (non-const version).
  double &out_phase_stddev() noexcept { return m_out_stddev; }

  /// Get the amplitude of the harmonic signal, i.e. if the harmonic is:
  /// A*sin(ωt) + B*cos(ωt), then its aplitude is sqrt(A^2 + B^2).
  double amplitude() const noexcept {
    return std::sqrt(m_in_phase * m_in_phase + m_out_phase * m_out_phase);
  }

private:
  ngpt::datetime<ngpt::milliseconds>
      m_start,      ///< Starting epoch of the harmonic
      m_stop;       ///< Ending epoch of the harmonic
  double m_afreq,   ///< Angular frequency (i.e. ω)
      m_in_phase,   ///< In-Phase component (amplitude)
      m_in_stddev,  ///< Std. deviation of In-Phase component
      m_out_phase,  ///< Out-Of-Phase component (amplitude)
      m_out_stddev; ///< Std. deviation of Out-Of-Phase component

}; // md_harmonics

} // namespace ngpt
#endif
