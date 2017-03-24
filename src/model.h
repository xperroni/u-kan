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

#ifndef UKAN_MODEL_H
#define UKAN_MODEL_H

#include "Eigen/Dense"

#include <memory>

namespace ukan {

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace process {

/**
 * @brief Process model specification.
 */
struct Model {
  /** @brief Dimension of the state vector. */
  int n_x;

  /** @brief Dimension of the augmented state vector. */
  int n_aug;

  /** @brief Sigma point factor. */
  double l;

  /** @brief Initial state parameter variance. */
  double s2_P0;

  /**
   * @brief Create a new model with given parameters.
   */
  Model(int n, int n_q, int p0):
    n_x(n),
    n_aug(n + n_q),
    l(3 - n_aug),
    s2_P0(p0)
  {
    // Nothing to do.
  }

  /**
   * @brief Augment the mean vector and covariance matrix with the process noise parameters.
   */
  virtual void augment(const VectorXd &x, const MatrixXd &P, VectorXd &x_aug, MatrixXd &P_aug) = 0;

  /**
   * @brief Normalize the state vector `x`.
   */
  virtual VectorXd normalize(VectorXd x) = 0;

  /**
   * @brief Compute the updated state using time difference `dt` and column `j` of matrix `S` as state.
   *
   * The updated state is written to column `j` of matrix `X`.
   */
  virtual void iterate(double dt, int j, const MatrixXd &S, MatrixXd &X) = 0;
};

} // namespace process

/**
 * @brief Reference-counting wrapper for process model objects.
 */
struct Model: std::shared_ptr<process::Model> {
  /**
   * @brief Default constructor.
   */
  Model() {
    // Nothing to do.
  }

  /**
   * @brief Create a new wrapper for a sensor measurement object.
   */
  Model(process::Model *m) {
    reset(m);
  }
};

} // namespace ukan

#endif
