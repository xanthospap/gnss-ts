#include "timeseries.hpp"
#include "ggdatetime/dtfund.hpp"
#include <algorithm>

ngpt::timeseries::timeseries(const ngpt::timeseries &ts, std::size_t start,
                             std::size_t end)
    : m_epochs(ts.m_epochs) {
  if (start || end) {
    if (!end || end > ts.m_data.size())
      end = ts.m_data.size();
    if (end < start) {
      throw std::domain_error(
          "timeseries: Invalid start/stop indexes for copy c'tor");
    }
    m_data = std::vector<ngpt::data_point>{ts.m_data.cbegin() + start,
                                           ts.m_data.cbegin() + end};
  } else {
    m_data = ts.m_data;
  }
}

std::size_t
ngpt::timeseries::mark(ngpt::pt_marker type,
                       ngpt::datetime<ngpt::milliseconds> start,
                       ngpt::datetime<ngpt::milliseconds> end) noexcept {
  if (!m_epochs)
    return -1;

  std::size_t start_idx, stop_idx;
  if (start == ngpt::datetime<ngpt::milliseconds>::min()) {
    start_idx = 0;
  } else {
    auto it = std::lower_bound(m_epochs->cbegin(), m_epochs->cend(), start);
    start_idx = std::distance(m_epochs->cbegin(), it);
  }

  if (end == ngpt::datetime<ngpt::milliseconds>::max()) {
    stop_idx = m_epochs->size();
  } else {
    auto it = std::lower_bound(m_epochs->cbegin(), m_epochs->cend(), end);
    stop_idx = std::distance(m_epochs->cbegin(), it);
  }

  std::transform(m_data.begin() + start_idx, m_data.begin() + stop_idx,
                 m_data.begin() + start_idx, [=](const data_point &p) {
                   data_point p2{p};
                   p2.m_flag.set(type);
                   return p2;
                 });
  return stop_idx - start_idx;
}
