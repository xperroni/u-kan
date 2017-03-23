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

#include "model_ctrv.h"

namespace ukan {

ModelCTRV::ModelCTRV(double s2_a, double s2_u, double s2_P0):
  Model(5, 2, s2_P0),
  s2_a_(s2_a),
  s2_u_(s2_u)
{
  // Nothing to do.
}

void ModelCTRV::augment(const VectorXd &x, const MatrixXd &P, VectorXd &x_aug, MatrixXd &P_aug) {
  // Create augmented mean vector.
  x_aug = VectorXd::Constant(n_aug, 0.0);
  x_aug.head(n_x) = x;

  // Create augmented state covariance matrix.
  P_aug = MatrixXd::Constant(n_aug, n_aug, 0.0);
  P_aug.block(0, 0, n_x, n_x) = P;
  P_aug.block(n_x, n_x, 2, 2) <<
    s2_a_, 0,
    0, s2_u_;
}

VectorXd ModelCTRV::difference(int j, const MatrixXd &X, const VectorXd &x) {
  static double pi2 = 2.0 * M_PI;

  VectorXd d = X.col(j) - x;

  // Perform angle normalization.
  while (d(3) > M_PI) {
    d(3) -= pi2;
  }

  // Perform angle normalization.
  while (d(3) < -M_PI) {
    d(3) += pi2;
  }

  return d;
}

void ModelCTRV::iterate(double dt, int j, const MatrixXd &S, MatrixXd &X) {
  double v = S(2, j);
  double o = S(3, j);
  double w = S(4, j);
  double a = S(5, j);
  double u = S(6, j);

  double dt2 = dt * dt;
  double cos_o = cos(o);
  double sin_o = sin(o);

  if (w != 0) {
    double w = S(4, j);
    X(0, j) = S(0, j) + (v / w) * (sin(o + w * dt) - sin_o) + 0.5 * dt2 * cos_o * a;
    X(1, j) = S(1, j) + (v / w) * (cos_o - cos(o + w * dt)) + 0.5 * dt2 * sin_o * a;
  }
  else {
    X(0, j) = S(0, j) + v * cos_o * dt + 0.5 * dt2 * cos_o * a;
    X(1, j) = S(1, j) + v * sin_o * dt + 0.5 * dt2 * sin_o * a;
  }

  X(2, j) = S(2, j) + dt * a;
  X(3, j) = S(3, j) + w * dt + 0.5 * dt2 * u;
  X(4, j) = S(4, j) + dt * u;
}

} // namespace ukan
