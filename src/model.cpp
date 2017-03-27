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

#include "model.h"

namespace ukan {

VectorXd normalize_angles(VectorXd x, const VectorXi &k) {
  static double pi2 = 2.0 * M_PI;

  for (int i = 0, n = k.rows(); i < n; ++i) {
    int j = k(i);

    // Normalize positive angle.
    while (x(j) > M_PI) {
      x(j) -= pi2;
    }

    // Normalize negative angle.
    while (x(j) < -M_PI) {
      x(j) += pi2;
    }
  }

  return x;
}

void weighted_sum(const MatrixXd &X, const VectorXd &w, const VectorXi &k, VectorXd &m, MatrixXd &C) {
  int n = X.cols();

  // Estimate state.
  m.fill(0.0);
  for (int j = 0; j < n; ++j) {
      m += X.col(j) * w(j);
  }

  // Estimate covariance matrix.
  C.fill(0.0);
  for (int j = 0; j < n; ++j) {
      VectorXd c = normalize_angles(X.col(j) - m, k);
      C += c * c.transpose() * w(j);
  }
}

namespace base {

Model::Model(int n, int n_q, double p0):
  n_x(n),
  n_aug(n + n_q),
  l(3 - n_aug),
  s2_P0(p0),
  w(2 * n_aug + 1)
{
  w(0) = l / (l + n_aug);
  w.tail(w.rows() - 1).fill(0.5 / (l + n_aug));
}


void Model::estimate(const MatrixXd &X, VectorXd &m, MatrixXd &C) {
  weighted_sum(X, w, k_, m, C);
}

MatrixXd Model::correlate(const State &x, const MatrixXd &X, const Measurement y, const MatrixXd &Z) {
  MatrixXd T = MatrixXd::Constant(n_x, y->rows(), 0.0);
  for (int j = 0, n = Z.cols(); j < n; ++j) {
    T += w(j) * normalize_angles(X.col(j) - x, k_) * (Z.col(j) - y)->transpose();
  }

  return T;
}

} // namespace base

} // namespace ukan
