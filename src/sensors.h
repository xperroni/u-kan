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

#ifndef KALMON_SENSORS_H
#define KALMON_SENSORS_H

#include "state.h"

#include <iostream>
#include <memory>

namespace ukan {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

using std::istream;
using std::ostream;
using std::shared_ptr;

/**
 * @brief The sensor array used to retrieve location measurements.
 */
struct Sensors {
  /**
   * @brief A measurement retrieved from one of the sensors in the array.
   */
  struct Measurement: VectorXd {
    /** @brief Time of measurement retrieval. */
    double timestamp;

    virtual Measurement *subtractFrom(const VectorXd &y) const = 0;

    /**
     * @brief Write this measurement to the given output stream.
     */
    virtual ostream &write(ostream &out) const = 0;

    /**
    * @brief Compute mean and covariance of a distribution from a set of samples.
    */
    virtual Measurement *estimate(const MatrixXd &Z, MatrixXd &C) const = 0;

    /**
     * @brief Transform the given matrix of column state vectors into a matrix of column measurement vectors.
     */
    virtual MatrixXd transform(const MatrixXd &X) const = 0;

    /**
     * @brief Convert this measurement into an state estimate.
     */
    virtual State state() const = 0;
  };

  /**
   * @brief Create a new sensor array with given parameters.
   */
  Sensors(
    double s2_px,
    double s2_py,
    double s2_d,
    double s2_r,
    double s2_v,
    const VectorXd &w
  );

  /**
   * @brief Return a new sensor measurement from the given input stream.
   */
  Measurement *operator () (istream &data) const;

private:
  /** @brief Laser measurement covariance matrix. */
  MatrixXd laserR_;

  /** @brief Radar measurement covariance matrix. */
  MatrixXd radarR_;

  /** @brief Radar angle parameters. */
  VectorXi radar_k_;

  /** @brief Vector of prediction weights. */
  VectorXd w_;
};

/**
 * @brief Reference-counting wrapper for sensor measurement objects.
 */
struct Measurement: shared_ptr<Sensors::Measurement> {
  /**
   * @brief Default constructor.
   */
  Measurement() {
    // Nothing to do.
  }

  /**
   * @brief Create a new wrapper for a sensor measurement object.
   */
  Measurement(Sensors::Measurement *z) {
    reset(z);
  }
};

Measurement operator - (const VectorXd &y, const Measurement z);

Measurement operator - (const Measurement y, const Measurement z);

/**
 * @brief Write a sensor measurement to the given output stream.
 */
ostream &operator << (ostream &data, Measurement &z);


} // namespace ukan

#endif
