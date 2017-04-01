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

#include "model_ctrv.h"

namespace ukan {

namespace ctrv {

/**
 * @brief An state of the CTRV model.
 */
struct State: process::State {
  /**
   * @brief Default constructor.
   */
  State():
    process::State(5)
  {
    // Nothing to do.
  }

  // See vector.h for documentation.
  base::Vector *copy() const {
    return new State(*this);
  }

  // See vector.h for documentation.
  virtual void iadd(const VectorXd &b) {
    *this += b;
    normalize_angle(*this, 3);
  }

  // See vector.h for documentation.
  virtual void imul(double b) {
    *this *= b;
    normalize_angle(*this, 3);
  }
};

Model::Model(double s2_a, double s2_u, double s2_P0):
  process::Model(5, 2, s2_P0)
{
  Q <<
    s2_a, 0,
    0, s2_u;
}

void Model::iterate(double dt, int j, const MatrixXd &S, MatrixXd &X) {
  double v = S(2, j);
  double o = S(3, j);
  double w = S(4, j);
  double a = S(5, j);
  double u = S(6, j);

  double dt2 = dt * dt;
  double cos_o = cos(o);
  double sin_o = sin(o);

  if (w != 0) {
    X(0, j) = S(0, j) + (v / w) * (sin(o + w * dt) - sin_o) + 0.5 * dt2 * cos_o * a;
    X(1, j) = S(1, j) + (v / w) * (cos_o - cos(o + w * dt)) + 0.5 * dt2 * sin_o * a;
  }
  else {
    X(0, j) = S(0, j) + v * cos_o * dt + 0.5 * dt2 * cos_o * a;
    X(1, j) = S(1, j) + v * sin_o * dt + 0.5 * dt2 * sin_o * a;
  }

  X(2, j) = v + dt * a;
  X(3, j) = o + w * dt + 0.5 * dt2 * u;
  X(4, j) = w + dt * u;
}

process::State *Model::newState() {
  return new ctrv::State();
}

} // namespace ctrv

} // namespace ukan
