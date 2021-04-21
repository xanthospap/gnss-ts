#include "timeseries.hpp"
#include "ggdatetime/dtfund.hpp"
#include <algorithm>

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

long ngpt::timeseries::cut(
    ngpt::datetime<ngpt::milliseconds> start,
    ngpt::datetime<ngpt::milliseconds> end,
    std::vector<ngpt::datetime<ngpt::milliseconds>> *epoch_vec) noexcept {

  if (!m_epochs || m_epochs->size() != m_data.size()) {
    return -1;
  }
  if (epoch_vec && epoch_vec->size() != m_data.size())
    return -1;

  auto it_start = std::lower_bound(m_epochs->cbegin(), m_epochs->cend(), start);
  auto start_idx = std::distance(m_epochs->cbegin(), it_start);
  auto it_end = std::lower_bound(m_epochs->cbegin(), m_epochs->cend(), end);
  auto end_idx = std::distance(m_epochs->cbegin(), it_end);

  m_data = std::vector<data_point>{m_data.begin() + start_idx,
                                   m_data.begin() + end_idx};

  if (epoch_vec) {
    std::rotate(epoch_vec->begin(), epoch_vec->begin() + start_idx,
                epoch_vec->end());
    epoch_vec->erase(epoch_vec->begin() + end_idx - start_idx,
                     epoch_vec->end());
  }
  return (long)(end_idx - start_idx);
}
  
ngpt::datetime<ngpt::milliseconds> ngpt::timeseries::central_epoch() const noexcept {
  auto first_epoch = (*m_epochs)[0], last_epoch = (*m_epochs)[m_epochs->size()-1];
  auto delta_dt = ngpt::delta_date(last_epoch, first_epoch);
  auto central_epoch{first_epoch};
  central_epoch += (delta_dt / 2);
  return central_epoch;
}
