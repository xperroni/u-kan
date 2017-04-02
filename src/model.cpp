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

#include "model.h"

namespace ukan {

void normalize_angle(VectorXd &x, int j) {
  static double pi2 = 2.0 * M_PI;

  // Normalize positive angle.
  while (x(j) > M_PI) {
    x(j) -= pi2;
  }

  // Normalize negative angle.
  while (x(j) < -M_PI) {
    x(j) += pi2;
  }
}

namespace process {

Model::Model(int n, int n_Q, double p0):
  n_x(n),
  n_aug(n + n_Q),
  l(1 - n_aug),
  s2_P0(p0),
  Q(n_Q, n_Q)
{
  // Nothing to do.
}

} // namespace process

} // namespace ukan
