/*
 * Copyright (c) Helio Perroni Filho <xperroni@gmail.com>
 *
 * This file is part of KalmOn.
 *
 * KalmOn is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KalmOn is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KalmOn. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef KALMON_KALMAN_FILTER_H
#define KALMON_KALMAN_FILTER_H

#include "model.h"
#include "sensors.h"
#include "state.h"

namespace ukan {

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @brief Kalman filter implementation.
 */
struct UKF {
  /** @brief Current process state. */
  State x;

  /** @brief State covariance matrix. */
  MatrixXd P;

  /**
   * @brief Default constructor.
   */
  UKF(Model model);

  /**
   * @brief Update the filter and return the current state.
   */
  State operator () (const Measurement z);

private:
  /** @brief Process model specification. */
  Model model_;

  /** @brief Matrix of state samples. */
  MatrixXd X_;

  /** @brief Timestamp of last received measurement. */
  double t_;

  /**
   * @brief Compute mean and covariance of a distribution from a set of samples.
   */
  void estimate(const MatrixXd &X, VectorXd &m, MatrixXd &C);

  /**
   * @brief Predict the state (and its covariance) of the process at the next step.
   *
   * @param dt Time difference between current and next steps.
   */
  void predict(double dt);

  /**
   * @brief Update the system state based on latest measurement.
   *
   * @param z Latest measurement.
   */
  void update(const Measurement z);
};

} // namespace ukan

#endif
