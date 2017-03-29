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

#ifndef UKAN_SENSORS_H
#define UKAN_SENSORS_H

#include "measurement.h"
#include "model.h"

#include <iostream>

namespace ukan {

using Eigen::MatrixXd;

using std::istream;
/**
 * @brief The sensor array used to retrieve location measurements.
 */
struct Sensors {
  /**
   * @brief Create a new sensor array with given parameters.
   */
  Sensors(
    Model model,
    double s2_px,
    double s2_py,
    double s2_d,
    double s2_r,
    double s2_v
  );

  /**
   * @brief Return a new sensor measurement from the given input stream.
   */
  Measurement operator () (istream &data) const;

private:
  /** @brief Process model. */
  Model model_;

  /** @brief Laser measurement covariance matrix. */
  MatrixXd laserR_;

  /** @brief Radar measurement covariance matrix. */
  MatrixXd radarR_;
};

} // namespace ukan

#endif
