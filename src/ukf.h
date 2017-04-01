/*
 * Copyright (c) Helio Perroni Filho <xperroni@gmail.com>
 *
 * This file is part of U-KAN.
 *
 * U-KAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * U-KAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with U-KAN. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UKAN_UKF_H
#define UKAN_UKF_H

#include "model.h"
#include "measurement.h"

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

  /** @brief Normalized Innovation Squared (NIS) */
  VectorXd e;

  /**
   * @brief Create a new filter for the given process model.
   */
  UKF(Model model);

  /**
   * @brief Update the filter and return the current state.
   */
  State operator () (const Measurement z);

private:
  /** @brief Process model specification. */
  Model model_;

  /** @brief Augmented state vector. */
  VectorXd x_aug_;

  /** @brief Augmented covariance matrix. */
  MatrixXd P_aug_;

  /** @brief Vector of prediction weights. */
  VectorXd w_;

  /** @brief Matrix of state samples. */
  MatrixXd X_;

  /** @brief Timestamp of last received measurement. */
  double t_;

  /**
   * @brief Estimate mean and covariance matrix of a sample set.
   */
  void estimate(const MatrixXd &X, Vector &m, MatrixXd &C);

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
