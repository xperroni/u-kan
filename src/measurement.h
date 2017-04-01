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

#ifndef UKAN_MEASUREMENT_H
#define UKAN_MEASUREMENT_H

#include "state.h"

#include <memory>
#include <stdexcept>

namespace ukan {

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace sensor {

/**
  * @brief A measurement retrieved from one of the sensors in the array.
  */
struct Measurement: base::Vector {
  /** @brief Time of measurement retrieval. */
  double timestamp;

  /**
   * @brief Create new measurement of given length.
   */
  Measurement(int n);

  /**
    * @brief Transform the given matrix of column state vectors into a matrix of column measurement vectors.
    */
  virtual MatrixXd transform(const MatrixXd &X) const = 0;

  /**
    * @brief Convert this measurement into an state estimate.
    */
  virtual State state() const = 0;

  /**
   * @brief Compute the measurement covariance matrix.
   */
  virtual MatrixXd R() const = 0;
};

} // namespace sensor

/**
 * @brief Reference-counting wrapper for sensor measurement objects.
 */
struct Measurement: Vector {
  /**
   * @brief Default constructor.
   */
  Measurement() {
    // Nothing to do.
  }

  /**
   * @brief Create a new wrapper for a sensor measurement object.
   */
  Measurement(sensor::Measurement *z):
    Vector(z)
  {
    // Nothing to do.
  }

  /**
   * @brief Create a new wrapper for a sensor measurement object.
   */
  Measurement(base::Vector *z):
    Vector(z)
  {
    // Nothing to do.
  }

  /**
   * @brief Create a new wrapper for a sensor measurement object.
   */
  Measurement(Vector x):
    Vector(x)
  {
    if (dynamic_cast<sensor::Measurement*>(get()) == nullptr) {
      throw std::runtime_error("Failed converting Vector input to Measurement");
    }
  }


  /**
   * @brief Return a pointer to the wrapped sensor measurement object.
   */
  sensor::Measurement *operator -> () {
    return (sensor::Measurement*) get();
  }

  /**
   * @brief Return a pointer to the wrapped sensor measurement object.
   */
  const sensor::Measurement *operator -> () const {
    return (const sensor::Measurement*) get();
  }
};

} // namespace ukan

#endif
