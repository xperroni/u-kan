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

#ifndef UKAN_MODEL_CTRV_H
#define UKAN_MODEL_CTRV_H

#include "model.h"

namespace ukan {

namespace ctrv {

/**
 * @brief The CTRV model describes non-linear motion in terms of position, speed, direction and ratio of turning.
 */
struct Model: process::Model {
  /**
   * @brief Create a new CTRV model with given parameters.
   */
  Model(double s2_a, double s2_u, double s2_P0);

  // See model.h for documentation.
  virtual void iterate(double dt, int j, const MatrixXd &S, MatrixXd &X);

  // See model.h for documentation.
  virtual process::State *newState();
};

} // namespace ctrv

} // namespace ukan

#endif
