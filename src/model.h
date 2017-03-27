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

#include "state.h"
#include "sensors.h"

#include "Eigen/Dense"

#include <memory>

namespace ukan {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

/**
 * @brief Normalize angle parameters of a vector to the range `[0, pi)`.
 *
 * @param x Vector to be normalized.
 *
 * @param k Indexes of angle parameters in the vector.
 */
VectorXd normalize_angles(VectorXd x, const VectorXi &k);

/**
 * @brief Compute the weighted sum (and accompanying covariance matrix) of the column vectors in matrix `X`.
 *
 * @param X Matrix of column vectors to be summed.
 *
 * @param w Weights of the individual column vectors.
 *
 * @param k Indexes of angle parameters in the vectors (so they can be normalized).
 *
 * @param m Summed vector.
 *
 * @param C Covariance matrix of the weighted sum.
 */
void weighted_sum(const MatrixXd &X, const VectorXd &w, const VectorXi &k, VectorXd &m, MatrixXd &C);

namespace base {

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

  /** @brief Vector of prediction weights. */
  VectorXd w;

  /**
   * @brief Create a new model with given parameters.
   */
  Model(int n, int n_q, double p0);

  /**
   * @brief Create a new model with given parameters.
   */
  template<class ...Indexes> Model(int n, int n_q, double p0, int j, Indexes... indexes):
    Model(n, n_q, p0)
  {
    k_.resize(1 + sizeof...(indexes));
    initIndexes(0, j, indexes...);
  }

  /**
   * @brief Compute mean and covariance of a distribution from a set of samples.
   */
  void estimate(const MatrixXd &X, VectorXd &m, MatrixXd &C);

  /**
   * @brief Augment the mean vector and covariance matrix with the process noise parameters.
   */
  virtual void augment(const VectorXd &x, const MatrixXd &P, VectorXd &x_aug, MatrixXd &P_aug) = 0;

  MatrixXd correlate(const State &x, const MatrixXd &X, const Measurement y, const MatrixXd &Z);

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

private:
  /** @brief Indexes of state angle parameters. */
  VectorXi k_;

  /**
   * @brief Leaf state of the recursive initializer (see below).
   */
  void initIndexes(int i, int j) {
    k_(i) = j;
  }

  /**
   * @brief Recursively initialize the index vector using the given parameter pack.
   */
  template<class ...Indexes> void initIndexes(int i, int j, Indexes... indexes) {
    k_(i) = j;
    initIndexes(i + 1, indexes...);
  }
};

} // namespace base

/**
 * @brief Reference-counting wrapper for process model objects.
 */
struct Model: std::shared_ptr<base::Model> {
  /**
   * @brief Default constructor.
   */
  Model() {
    // Nothing to do.
  }

  /**
   * @brief Create a new wrapper for a sensor measurement object.
   */
  Model(base::Model *m) {
    reset(m);
  }
};

} // namespace ukan

#endif
