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

#include "sensors.h"

#include "model.h"

using std::string;

namespace ukan {

/**
 * @brief A measurement retrieved from the laser sensor.
 */
struct MeasurementLaser: Sensors::Measurement {
  /**
   * @brief Retrieve a new laser measurement from the given input stream.
   */
  MeasurementLaser(istream &data, const VectorXd &w, const MatrixXd &R):
    w_(w),
    R_(R)
  {
    resize(2);
    Sensors::Measurement &z = *this;
    data >> z(0) >> z(1) >> z.timestamp;
  }

  virtual Sensors::Measurement *subtractFrom(const VectorXd &y) const {
    Sensors::Measurement *d = new MeasurementLaser(*this);
    VectorXd &v = *d;
    v = y - *this;
    return d;
  }

  virtual ostream &write(ostream &out) const {
      const Sensors::Measurement &z = *this;
      return out << z(0) << "\t" << z(1);
  }

  virtual State state() const {
    const Sensors::Measurement &z = *this;
    return State(z(0), z(1), 0.0, 0.0, 0.0);
  }

  virtual Sensors::Measurement *estimate(const MatrixXd &Z, MatrixXd &S) const {
    Sensors::Measurement *z = new MeasurementLaser(*this);
    int m = rows();
    S.resize(m, m);
    weighted_sum(Z, w_, VectorXi(), *z, S);
    S += R_;
    return z;
  }

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

private:
  /** @brief Vector of prediction weights. */
  const VectorXd &w_;

  /** @brief Covariance matrix for this measurement. */
  const MatrixXd &R_;
};

/**
 * @brief A measurement retrieved from the radar sensor.
 */
struct MeasurementRadar: Sensors::Measurement {
  /**
   * @brief Retrieve a new radar measurement from the given input stream.
   */
  MeasurementRadar(istream &data, const VectorXd &w, const VectorXi &k, const MatrixXd &R):
    w_(w),
    k_(k),
    R_(R)
  {
    resize(3);
    Sensors::Measurement &z = *this;
    data >> z(0) >> z(1) >> z(2) >> z.timestamp;
  }

  virtual Sensors::Measurement *subtractFrom(const VectorXd &y) const {
    Sensors::Measurement *d = new MeasurementRadar(*this);
    VectorXd &v = *d;
    v = normalize_angles(y - *this, k_);
    return d;
  }

  virtual ostream &write(ostream &out) const {
    const Sensors::Measurement &z = *this;
    double r = z(0);
    double p = z(1);
    return out << r * ::cos(p) << "\t" << r * ::sin(p);
  }

  virtual State state() const {
    const Sensors::Measurement &z = *this;
    double d = z(0);
    double r = z(1);
    double v = z(2);

    double cos_r = ::cos(r);
    double sin_r = ::sin(r);

    // r is not really the same value as the psi state parameter,
    // but it's a reasonable approximation.
    return State(d * cos_r, d * sin_r, v, r, 0.0);
  }

  virtual Sensors::Measurement *estimate(const MatrixXd &Z, MatrixXd &S) const {
    Sensors::Measurement *z = new MeasurementRadar(*this);
    int m = rows();
    S.resize(m, m);
    weighted_sum(Z, w_, k_, *z, S);
    S += R_;
    return z;
  }

  virtual MatrixXd transform(const MatrixXd &X) const {
    int m = rows();
    int n = X.cols();
    MatrixXd Z(m, n);
    for (int j = 0; j < n; ++j) {
      double x = X(0, j);
      double y = X(1, j);
      double v = X(2, j);
      double o = X(3, j);

      Z(0, j) = ::sqrt(x*x + y*y);
      Z(1, j) = atan2(y, x);
      Z(2, j) = (x * ::cos(o) + y * ::sin(o)) * v / Z(0, j);
    }

    return Z;
  }

private:
  /** @brief Vector of prediction weights. */
  const VectorXd &w_;

  /** @brief Vector of angle indexes. */
  const VectorXi &k_;

  /** @brief Covariance matrix for this measurement. */
  const MatrixXd &R_;
};

Sensors::Sensors(
  double s2_px,
  double s2_py,
  double s2_d,
  double s2_r,
  double s2_v,
  const VectorXd &w
):
  laserR_(2, 2),
  radarR_(3, 3),
  radar_k_(1),
  w_(w)
{
  laserR_ <<
    s2_px, 0,
    0, s2_py;

  radarR_ <<
    s2_d, 0, 0,
    0, s2_r, 0,
    0, 0, s2_v;

  radar_k_(0) = 1;
}

Sensors::Measurement *Sensors::operator () (istream &data) const {
    string sensor;
    data >> sensor;

    if (sensor == "L") {
      return new MeasurementLaser(data, w_, laserR_);
    }
    else /* if (sensor_type == "R") */ {
      return new MeasurementRadar(data, w_, radar_k_, radarR_);
    }
}

Measurement operator - (const Measurement y, const Measurement z) {
  return z->subtractFrom(*y);
}

Measurement operator - (const VectorXd &y, const Measurement z) {
  return z->subtractFrom(y);
}

ostream &operator << (ostream &data, Measurement &z) {
  return z->write(data);
}

} // namespace ukan
