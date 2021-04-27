#ifndef __NGPT_KALMAN_EST_HPP__
#define __NGPT_KALMAN_EST_HPP__
#include "eigen3/Eigen/Core"
#include "model.hpp"
#include "timeseries.hpp"
namespace ngpt {
Eigen::VectorXd kalman(const ngpt::timeseries &ts, const ngpt::ts_model &mdl,
                       double sigma0);
} // namespace ngpt
#endif
