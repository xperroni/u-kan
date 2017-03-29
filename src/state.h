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

#ifndef UKAN_STATE_H
#define UKAN_STATE_H

#include "vector.h"

#include <vector>

namespace ukan {

using Eigen::VectorXd;

using std::vector;

namespace process {

/**
 * @brief An state of the process handled by the Kalman filter.
 */
struct State: base::Vector {
  /**
   * @brief Create new state of given length.
   */
  State(int n);

  // See vector.h for documentation.
  virtual base::Vector *copy() const;
};

} // namespace process

/**
 * @brief Reference-counting wrapper for states.
 */
struct State: Vector {
  /**
   * @brief Default constructor.
   */
  State() {
    // Nothing to do.
  }

  /**
   * @brief Create a new wrapper for a sensor measurement object.
   */
  State(process::State *x) {
    reset(x);
  }
};

/**
 * @brief Compute the Root Mean Square Error (RMSE) between two state sequences.
 */
VectorXd RMSE(const vector<State> &estimates, const vector<State> &ground_truth);

} // namespace ukan

#endif
