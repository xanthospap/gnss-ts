#ifndef __NGPT_JUMP_MODEL_HPP__
#define __NGPT_JULP_MODEL_HPP__

namespace ngpt {
/// A class to represent a "jump" (i.e. offset) of a time-series. A 'jump' is
/// determined by the value of the offset, the epoch (time tag) at which it
/// occured and a possible value for its std. deviation.
class md_jump {
public:
  /// Constructor.
  /// @param[in] start  The epoch the jump happened.
  /// @param[in] offset The value of the offset (default is 0)
  /// @param[in] stdd   The std. deviation of the offset value/estimate
  md_jump(ngpt::datetime<ngpt::milliseconds> start, double offset = 0e0,
          double stdd = 0e0) noexcept
      : m_start{start}, m_offset{offset}, m_stddev{stdd} {};

  /// Get the epoch the jump happened at (const).
  datetime<ngpt::milliseconds> start() const noexcept { return m_start; }

  /// Get the value of the offset (const).
  double value() const noexcept { return m_offset; }

  /// Get the value of the offset (non-const).
  double &value() noexcept { return m_offset; }

  /// Get the std. deviation value of the offset (const).
  double stddev() const noexcept { return m_stddev; }

  /// Get the value of the offset (non-const).
  double &stddev() noexcept { return m_stddev; }

private:
  ngpt::datetime<ngpt::milliseconds> m_start; ///< When the "jump" occured.
  double m_offset; ///< Value (i.e. offset amplitude).
  double m_stddev; ///< Std. deviation of the estimated value.

}; // md_jump

} // namespace ngpt
#endif
