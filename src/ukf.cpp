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

#include "ukf.h"

#include "settings.h"

using namespace std;

namespace ukan {

// Size of the filter state.
static const int len_x = 5;

// Size of the augmented filter state.
static const int len_a = len_x + 2;

// Identity matrix.
static const MatrixXd I = MatrixXd::Identity(len_x, len_x);

UKF::UKF():
  P(0, 0)
{
    // Nothing to do.
}

State UKF::operator () (const Measurement z) {
  if (P.rows() == 0) {
    // Initialize state with first measurement.
    x = z->x();

    // Initialize covariance matrix.
    P.resize(len_x, len_x);
    P = I * getSettings().s2_P0;

    // Initialize timestamp.
    t_ = z->timestamp;

    return x;
  }

  // Compute elapsed time in seconds.
  double t = z->timestamp;
  double dt = (t - t_) / 1000000.0;

  // Predict new state.
  predict(dt);

  // Update predictions.
  update(z);

  // Update timestamp.
  t_ = t;

  return x;
}

inline void model(double dt, int j, const MatrixXd &S, MatrixXd &X) {
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

void UKF::predict(double dt) {
  // Create augmented mean vector.
  VectorXd x_aug = VectorXd::Constant(len_a, 0.0);
  x_aug.head(len_x) = x;

  // Create augmented state covariance matrix.
  double s2_a = getSettings().s2_a;
  double s2_u = getSettings().s2_u;
  MatrixXd P_aug = MatrixXd::Constant(len_a, len_a, 0.0);
  P_aug.block(0, 0, len_x, len_x) = P;
  P_aug.block(len_x, len_x, 2, 2) <<
    s2_a, 0,
    0, s2_u;

  // Create square root matrix.
  double l = getSettings().lambda;
  MatrixXd A = P_aug.llt().matrixL();
  A *= sqrt(l + len_a);

  // Create sigma point matrix.
  int n = 2 * len_a + 1;
  MatrixXd S = MatrixXd(len_a, n);
  S.col(0) = x_aug;
  S.block(0, 1, len_a, len_a) = A.colwise() + x_aug;
  S.block(0, len_a + 1, len_a, len_a) = (-A).colwise() + x_aug;

  // Create matrix with predicted sigma points as columns.
  MatrixXd X = MatrixXd(len_x, n);
  for (int j = 0; j < n; ++j) {
      model(dt, j, S, X);
  }

  // Compute prediction weights.
  double w0 = l / (l + len_a);
  double wj = 0.5 / (l + len_a);

  // Predict state.
  x = X.col(0) * w0;
  for (int j = 1; j < n; ++j) {
      x += X.col(j) * wj;
  }

  VectorXd c0 = X.col(0) - x;
  P = c0 * c0.transpose() * w0;
  for (int j = 1; j < n; ++j) {
      VectorXd cj = X.col(j) - x;
      P += cj * cj.transpose() * wj;
  }
}

void UKF::update(const Measurement z) {
  // Retrieve the measurement model matrix.
  MatrixXd H = z->H(x);
  MatrixXd Ht = H.transpose();

  VectorXd y = *z - H * x;
  MatrixXd S = H * P * Ht + z->R();
  MatrixXd K = P * Ht * S.inverse();

  // Update next state and covariance matrix.
  x += K * y;
  P = (I - K * H) * P;
}

} // namespace ukan
