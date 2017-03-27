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

#ifndef KALMON_STATE_H
#define KALMON_STATE_H

#include "Eigen/Dense"

#include <iostream>
#include <vector>

namespace ukan {

using Eigen::EigenBase;
using Eigen::VectorXd;

using std::istream;
using std::ostream;
using std::vector;

/**
 * @brief An state of the process handled by the Kalman filter.
 */
struct State: VectorXd {
  /**
   * @brief Default constructor.
   */
  State();

  /**
   * @brief Create new state of given length.
   */
  State(int n);

  /**
   * @brief Create a new state from a linear algebra object or operation.
   */
  template<class T> State(const EigenBase<T> &values):
    VectorXd()
  {
    VectorXd &x = *this;
    x = values;
  }

  /**
   * @brief Create a new state with given parameters.
   */
  template<class ...Values> State(double v, Values... values):
    VectorXd(1 + sizeof...(values))
  {
    init(0, v, values...);
  }

private:
  /**
   * @brief Leaf state of the recursive initializer (see below).
   */
  void init(int i, double v) {
    (*this)(i) = v;
  }

  /**
   * @brief Recursively initialize this state using the given parameter pack.
   */
  template<class ...Values> void init(int i, double v, Values... values) {
    (*this)(i) = v;
    init(i + 1, values...);
  }
};

/**
 * @brief Read a state from the given input stream.
 */
istream &operator >> (istream &in, State &x);

/**
 * @brief Write a state to the given output stream.
 */
ostream &operator << (ostream &out, const State &x);

/**
 * @brief Compute the Root Mean Square Error (RMSE) between two state sequences.
 */
State RMSE(const vector<State> &estimates, const vector<State> &ground_truth);

} // namespace ukan

#endif
