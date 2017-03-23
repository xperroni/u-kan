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

using namespace std;

namespace ukan {

UKF::UKF(Model model):
  P(0, 0),
  model_(model),
  X_(model->n_x, 2 * model->n_aug + 1)
{
  // Nothing to do.
}

State UKF::operator () (const Measurement z) {
  if (P.rows() == 0) {
    // Initialize state with first measurement.
    x = z->x();

    // Initialize covariance matrix.
    int n = model_->n_x;
    P.resize(n, n);
    P = MatrixXd::Identity(n, n) * model_->s2_P0;

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

void UKF::estimate(const MatrixXd &X, VectorXd &m, MatrixXd &C) {
  // Compute prediction weights.
  double l = model_->l;
  double n_aug = model_->n_aug;
  double w0 = l / (l + n_aug);
  double wj = 0.5 / (l + n_aug);

  // Estimate state.
  m = X.col(0) * w0;
  int n = X.cols();
  for (int j = 1; j < n; ++j) {
      m += X.col(j) * wj;
  }

  // Estimate covariance matrix.
  VectorXd c0 = model_->difference(0, X, m);
  C = c0 * c0.transpose() * w0;
  for (int j = 1; j < n; ++j) {
      VectorXd cj = model_->difference(j, X, m);
      C += cj * cj.transpose() * wj;
  }
}

void UKF::predict(double dt) {
  // Create augmented mean vector and covariance matrix.
  VectorXd x_aug;
  MatrixXd P_aug;
  model_->augment(x, P, x_aug, P_aug);

  // Create square root matrix.
  double l = model_->l;
  int n_aug = model_->n_aug;
  MatrixXd A = P_aug.llt().matrixL();
  A *= sqrt(l + n_aug);

  // Create sigma point matrix.
  int n = X_.cols();
  MatrixXd S(n_aug, n);
  S.col(0) = x_aug;
  S.block(0, 1, n_aug, n_aug) = A.colwise() + x_aug;
  S.block(0, n_aug + 1, n_aug, n_aug) = (-A).colwise() + x_aug;

  // Create matrix with predicted sigma points as columns.
  for (int j = 0; j < n; ++j) {
      model_->iterate(dt, j, S, X_);
  }

  // Estimate process state and covariance matrix.
  estimate(X_, x, P);
}

void UKF::update(const Measurement z) {
  VectorXd s;
  MatrixXd S;
  estimate(z->transform(X_), s, S);
  S += z->R();

  // Retrieve the measurement model matrix.
  MatrixXd H = z->H(x);
  MatrixXd Ht = H.transpose();

  VectorXd y = *z - s;
  MatrixXd K = P * Ht * S.inverse();

  // Update next state and covariance matrix.
  x += K * y;
  //P = (I - K * H) * P;
}

} // namespace ukan
