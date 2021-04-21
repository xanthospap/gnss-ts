#include "model.hpp"

std::size_t ngpt::model::num_parameters() const noexcept {
  std::size_t parameters = 2;
  parameters += md_jumps.size();
  parameters+=md_harmonics.size()*2;
  parameters+=m_vel_changes.size();
  return parameters;
}
