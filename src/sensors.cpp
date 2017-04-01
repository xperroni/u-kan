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

#include "sensors.h"

#include "model.h"

using std::string;

namespace ukan {

/**
 * @brief A measurement retrieved from the laser sensor.
 */
struct MeasurementLaser: sensor::Measurement {
  /**
   * @brief Retrieve a new laser measurement from the given input stream.
   */
  MeasurementLaser(istream &data, Model model, const MatrixXd &R):
    sensor::Measurement(2),
    model_(model),
    R_(R)
  {
    sensor::Measurement &z = *this;
    data >> z(0) >> z(1) >> z.timestamp;
    if (z(0) == 0 && z(1) == 0)
      z(0) = z(1) = 0.0001;
  }

  // See vector.h for documentation.
  base::Vector *copy() const {
    return new MeasurementLaser(*this);
  }

  // See measurement.h for documentation.
  virtual ostream &write(ostream &out) const {
      const sensor::Measurement &z = *this;
      return out << z(0) << "\t" << z(1);
  }

  // See measurement.h for documentation.
  virtual MatrixXd transform(const MatrixXd &X) const {
    int m = rows();
    int n = X.cols();
    MatrixXd Z(m, n);
    for (int j = 0; j < n; ++j) {
      Z(0, j) = X(0, j);
      Z(1, j) = X(1, j);
    }

    return Z;
  }

  // See measurement.h for documentation.
  virtual State state() const {
    State x = model_->newState();
    const sensor::Measurement &z = *this;
    *x << z(0), z(1), 0.0, 0.0, 0.0;
    return x;
  }

  // See measurement.h for documentation.
  virtual MatrixXd R() const {
    return R_;
  }

private:
  /** @brief Process model. */
  Model model_;

  /** @brief Covariance matrix for this measurement. */
  const MatrixXd &R_;
};

/**
 * @brief A measurement retrieved from the radar sensor.
 */
struct MeasurementRadar: sensor::Measurement {
  /**
   * @brief Retrieve a new radar measurement from the given input stream.
   */
  MeasurementRadar(istream &data, Model model, const MatrixXd &R):
    sensor::Measurement(3),
    model_(model),
    R_(R)
  {
    sensor::Measurement &z = *this;
    data >> z(0) >> z(1) >> z(2) >> z.timestamp;
    if (z(0) == 0 && z(1) == 0)
      z(0) = z(1) = 0.0001;
  }

  // See vector.h for documentation.
  base::Vector *copy() const {
    return new MeasurementRadar(*this);
  }

  // See vector.h for documentation.
  virtual void iadd(const VectorXd &b) {
    *this += b;
    normalize_angle(*this, 1);
  }

  // See vector.h for documentation.
  virtual void imul(const VectorXd &b) {
    *this *= b;
    normalize_angle(*this, 1);
  }

  // See measurement.h for documentation.
  virtual ostream &write(ostream &out) const {
    const sensor::Measurement &z = *this;
    double r = z(0);
    double p = z(1);
    return out << r * ::cos(p) << "\t" << r * ::sin(p);
  }

  // See measurement.h for documentation.
  virtual MatrixXd transform(const MatrixXd &X) const {
    int m = rows();
    int n = X.cols();
    MatrixXd Z(m, n);
    for (int j = 0; j < n; ++j) {
      double x = X(0, j);
      double y = X(1, j);
      double v = X(2, j);
      double o = X(3, j);

      double d = ::sqrt(x*x + y*y);

      Z(0, j) = d;
      Z(1, j) = atan2(y, x);
      Z(2, j) = (x * ::cos(o) + y * ::sin(o)) * v / d;
    }

    return Z;
  }

  // See measurement.h for documentation.
  virtual State state() const {
    const sensor::Measurement &z = *this;
    double d = z(0);
    double r = z(1);
    double v = z(2);

    double cos_r = ::cos(r);
    double sin_r = ::sin(r);

    // r is not really the same value as the psi state parameter,
    // but it's a reasonable approximation.
    State x = model_->newState();
    *x << d * cos_r, d * sin_r, v, r, 0.0;
    return x;
  }

  // See measurement.h for documentation.
  virtual MatrixXd R() const {
    return R_;
  }

private:
  /** @brief Process model. */
  Model model_;

  /** @brief Covariance matrix for this measurement. */
  const MatrixXd &R_;
};

Sensors::Sensors(
  Model model,
  double s2_px,
  double s2_py,
  double s2_d,
  double s2_r,
  double s2_v
):
  model_(model),
  laserR_(2, 2),
  radarR_(3, 3)
{
  laserR_ <<
    s2_px, 0,
    0, s2_py;

  radarR_ <<
    s2_d, 0, 0,
    0, s2_r, 0,
    0, 0, s2_v;
}

Measurement Sensors::operator () (istream &data) const {
    string sensor;
    data >> sensor;

    if (sensor == "L") {
      return new MeasurementLaser(data, model_, laserR_);
    }
    else /* if (sensor_type == "R") */ {
      return new MeasurementRadar(data, model_, radarR_);
    }
}

} // namespace ukan
