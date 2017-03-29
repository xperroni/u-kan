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

#ifndef UKAN_MODEL_H
#define UKAN_MODEL_H

#include "state.h"

#include <memory>

namespace ukan {

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @brief Normalize an angle parameter to the range (-pi, pi).
 */
void normalize_angle(VectorXd &x, int j);

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

  /** @brief Process covariance matrix. */
  MatrixXd Q;

  /**
   * @brief Create a new model with given parameters.
   */
  Model(int n, int n_Q, double p0);

  /**
   * @brief Compute the updated state using time difference `dt` and column `j` of matrix `S` as state.
   *
   * The updated state is written to column `j` of matrix `X`.
   */
  virtual void iterate(double dt, int j, const MatrixXd &S, MatrixXd &X) = 0;

  /**
   * @brief Create a new, "empty" process state.
   */
  virtual State *newState() = 0;
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
