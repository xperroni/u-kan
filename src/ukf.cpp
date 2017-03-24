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
  X_(model->n_x, 2 * model->n_aug + 1),
  w_(X_.cols())
{
  double l = model_->l;
  double n_aug = model_->n_aug;
  w_(0) = l / (l + n_aug);
  w_.tail(w_.rows() - 1).fill(0.5 / (l + n_aug));
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
  int n = X.cols();

  // Estimate state.
  m.fill(0.0);
  for (int j = 0; j < n; ++j) {
      m += X.col(j) * w_(j);
  }

  // Estimate covariance matrix.
  C.fill(0.0);
  for (int j = 0; j < n; ++j) {
      VectorXd c = model_->normalize(X.col(j) - m);
      C += c * c.transpose() * w_(j);
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
  MatrixXd Z = z->transform(X_);

  // Estimate measurement mean and covariance.
  VectorXd y;
  MatrixXd S;
  estimate(Z, y, S);
  S += z->R();

  // Compute cross-correlation matrix between state and measurement sigma points.
  MatrixXd T = MatrixXd::Constant(model_->n_x, z->rows(), 0.0);
  for (int j = 0, n = Z.cols(); j < n; ++j) {
    T += w_(j) * model_->normalize(X_.col(j) - x) * model_->normalize(Z.col(j) - y).transpose();
  }

  // Compute Kalman gain.
  MatrixXd K = T * S.inverse();

  //update state mean and covariance matrix
  x += K * model_->normalize(*z - y);
  P -= K * S * K.transpose();
}

} // namespace ukan
