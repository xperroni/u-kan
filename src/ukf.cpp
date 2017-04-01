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

#include "ukf.h"

using namespace std;

namespace ukan {

UKF::UKF(Model model):
  x(model->newState()),
  P(0, 0),
  e(1),
  model_(model),
  w_(2 * model->n_aug + 1),
  X_(model->n_x, 2 * model->n_aug + 1)
{
  // Model parameters.
  int l = model_->l;
  int n_x = model_->n_x;
  int n_aug = model_->n_aug;
  int n_Q = n_aug - n_x;

  // Initialize augmented state and covariance matrices.
  x_aug_ = VectorXd::Constant(n_aug, 0.0);
  P_aug_ = MatrixXd::Constant(n_aug, n_aug, 0.0);
  P_aug_.block(n_x, n_x, n_Q, n_Q) = model->Q;

  // Initialize prediction weight vector.
  w_.fill(0.5 / (l + n_aug));
  w_(0) = l / (l + n_aug);

}

State UKF::operator () (const Measurement z) {
  if (P.rows() == 0) {
    // Initialize state with first measurement.
    x = z->state();

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

  return x->copy();
}

void UKF::estimate(const MatrixXd &X, Vector &m, MatrixXd &C) {
  int n = X.cols();

  // Estimate state.
  m->fill(0.0);
  for (int j = 0; j < n; ++j) {
      m += w_(j) * X.col(j);
  }

  // Estimate covariance matrix.
  C.fill(0.0);
  for (int j = 0; j < n; ++j) {
      Vector c = X.col(j) - m;
      C += w_(j) * c * c->transpose();
  }
}

void UKF::predict(double dt) {
  int n_x = model_->n_x;
  int n_aug = model_->n_aug;
  double l = model_->l;

  // Update augmented state vector and covariance matrix.
  x_aug_.head(n_x) = *x;
  P_aug_.block(0, 0, n_x, n_x) = P;

  // Create weighted square root matrix.
  MatrixXd A = P_aug_.llt().matrixL();
  A *= sqrt(l + n_aug);

  // Create sigma point matrix.
  int n = X_.cols();
  MatrixXd S(n_aug, n);
  S.col(0) = x_aug_;
  S.block(0, 1, n_aug, n_aug) = A.colwise() + x_aug_;
  S.block(0, 1 + n_aug, n_aug, n_aug) = (-A).colwise() + x_aug_;

  // Create matrix with predicted sigma points as columns.
  for (int j = 0; j < n; ++j) {
      model_->iterate(dt, j, S, X_);
  }

  // Estimate process state and covariance matrix.
  estimate(X_, x, P);
}

void UKF::update(const Measurement z) {
  // Compute sigma points in measurement space.
  MatrixXd Z = z->transform(X_);

  // Compute measurement cross-correlation matrix.
  Measurement y = z->copy();
  MatrixXd S(z->rows(), z->rows());
  estimate(Z, y, S);
  S += z->R();

  // Compute cross-correlation matrix between state and measurement sigma points.
  MatrixXd T = MatrixXd::Constant(x->rows(), y->rows(), 0.0);
  for (int j = 0, n = Z.cols(); j < n; ++j) {
    T += w_(j) * (X_.col(j) - x) * (Z.col(j) - y)->transpose();
  }

  // Compute Kalman gain.
  MatrixXd Si = S.inverse();
  MatrixXd K = T * Si;

  // Update state mean and covariance matrix.
  Measurement d = z - y;
  x += K * d;
  P -= K * S * K.transpose();

  // Compute NIS.
  e = d->transpose() * Si * d;
}

} // namespace ukan
