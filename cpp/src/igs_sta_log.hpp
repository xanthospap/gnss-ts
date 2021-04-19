#ifndef __NGPT_IGS_LOG_FILE_HPP__
#define __NGPT_IGS_LOG_FILE_HPP_

#include <cassert>
#include <fstream>
#include <ggdatetime/dtcalendar.hpp>
#include <string>
#include <vector>

namespace ngpt {

class igs_log {
public:
  explicit igs_log(const std::string &filename);

  /// No copy constructor.
  igs_log(const igs_log &) = delete;

  /// No assignment operator.
  igs_log &operator=(const igs_log &) = delete;

  /// Move constructor.
  igs_log(igs_log &&ec)
      : m_filename{std::move(ec.m_filename)}, m_ifs{std::move(ec.m_ifs)} {}

  /// Move assignment operator.
  igs_log &operator=(igs_log &&ec) {
    if (this != &ec) {
      m_filename = std::move(ec.m_filename);
      m_ifs = std::move(ec.m_ifs);
    }
    return *this;
  }

  std::vector<ngpt::datetime<ngpt::milliseconds>> receiver_changes();

  std::vector<ngpt::datetime<ngpt::milliseconds>> antenna_changes();

private:
  /// The name of the log file.
  std::string m_filename;
  /// The input file stream
  std::ifstream m_ifs;

  /// @brief Go to the begining of the stream (input file).
  void rewind() noexcept {
    m_ifs.seekg(0, std::ios::beg);
    return;
  }

  int read_receiver_block(std::string &receiver_type,
                          ngpt::datetime<ngpt::milliseconds> &start,
                          ngpt::datetime<ngpt::milliseconds> &stop);

  int read_antenna_block(std::string &antenna_type,
                         ngpt::datetime<ngpt::milliseconds> &start,
                         ngpt::datetime<ngpt::milliseconds> &stop);
}; // igs_log

} // namespace ngpt

#endif
