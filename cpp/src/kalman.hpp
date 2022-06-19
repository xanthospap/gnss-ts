#ifndef __NGPT_KALMAN_EST_HPP__
#define __NGPT_KALMAN_EST_HPP__
#include "eigen3/Eigen/Core"
#include "model.hpp"
#include "timeseries.hpp"
namespace dso {
Eigen::VectorXd kalman(const dso::timeseries &ts, const dso::ts_model &mdl,
                       double sigma0);
} // namespace dso
#endif