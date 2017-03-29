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

#include "state.h"

namespace ukan {

namespace process {

State::State(int n):
  Vector(n)
{
  // Nothing to do.
}

base::Vector *State::copy() const {
  return new State(*this);
}

} // namespace process

VectorXd RMSE(const vector<State> &estimates, const vector<State> &ground_truth) {
  int n = estimates.size();
  VectorXd rmse = VectorXd::Constant(estimates[0]->rows(), 0.0);
  for(int i = 0; i < n; ++i) {
    const VectorXd &a = *estimates[i];
    const VectorXd &b = *ground_truth[i];

    VectorXd d = a - b;
    VectorXd d2 = d.array() * d.array();
    rmse += d2;
  }

  rmse /= n;

  rmse = rmse.array().sqrt();

  return rmse;
}

} // namespace ukan
