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

#ifndef UKAN_MODEL_CTRV_H
#define UKAN_MODEL_CTRV_H

#include "model.h"

namespace ukan {

struct ModelCTRV: base::Model {
  /**
   * @brief Create a new CTRV model with given parameters.
   */
  ModelCTRV(double s2_a, double s2_u, double s2_P0);

  // See model.h for documentation.
  virtual void augment(const VectorXd &x, const MatrixXd &P, VectorXd &x_aug, MatrixXd &P_aug);

  // See model.h for documentation.
  virtual VectorXd normalize(VectorXd x);

  // See model.h for documentation.
  virtual void iterate(double dt, int j, const MatrixXd &S, MatrixXd &X);

private:
  /** @brief Linear acceleration variance. */
  double s2_a_;

  /** @brief Radial acceleration variance. */
  double s2_u_;
};

} // namespace ukan

#endif
