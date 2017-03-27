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

#include "state.h"

namespace ukan {

State::State():
  VectorXd()
{
  // Nothing to do.
}

State::State(int n):
  VectorXd(n)
{
  VectorXd &x = *this;
  for (int i = 0; i < n; ++i) {
    x(i) = 0;
  }
}

istream &operator >> (istream &in, State &x) {
  for (int i = 0, n = x.rows(); i < n; ++i) {
    in >> x(i);
  }

  return in;
}

ostream &operator << (ostream &out, const State &x) {
  for (int i = 0, n = x.rows(); i < n;) {
    out << x(i);
    if (++i >= n) {
      break;
    }

    out << '\t';
  }

  return out;
}

State RMSE(const vector<State> &estimates, const vector<State> &ground_truth) {
  int n = estimates.size();
  State rmse(estimates[0].rows());
  for(int i = 0; i < n; ++i) {
    const State &a = estimates[i];
    const State &b = ground_truth[i];

    State c = a - b;
    State d2 = c.array() * c.array();
    rmse += d2;
  }

  rmse /= n;

  rmse = rmse.array().sqrt();

  return rmse;
}

} // namespace ukan
