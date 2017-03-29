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

#ifndef UKAN_VECTOR_H
#define UKAN_VECTOR_H

#include "Eigen/Dense"

#include <iostream>
#include <memory>
#include <vector>

namespace ukan {

using Eigen::EigenBase;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using std::istream;
using std::ostream;
using std::vector;

namespace base {

/**
 * @brief A multidimensional value.
 */
struct Vector: VectorXd {
  /**
   * @brief Create new vector of given length.
   */
  Vector(int n);

  /**
   * @brief Create a copy of this vector.
   */
  virtual Vector *copy() const;

  /**
   * @brief Compute the addition between this vector and another one.
   */
  virtual Vector *add(const VectorXd &b) const;

  /**
   * @brief Compute the parameter-wise product between this vector and a scalar.
   */
  virtual Vector *mul(double b) const;

  /**
   * @brief Increment this vector by the given value.
   */
  virtual void iadd(const VectorXd &b);

  /**
   * @brief Multiply this vector by the given value.
   */
  virtual void imul(double b);

  /**
    * @brief Write this vector to the given output stream.
    */
  virtual ostream &write(ostream &out) const;
};

} // namespace base

/**
 * @brief Reference-counting wrapper for states.
 */
struct Vector: std::shared_ptr<base::Vector> {
  /**
   * @brief Default constructor.
   */
  Vector() {
    // Nothing to do.
  }

  /**
   * @brief Create a new wrapper for a sensor measurement object.
   */
  Vector(base::Vector *x) {
    reset(x);
  }

  /**
   * @brief Return the state parameter at index `k`.
   */
  double &operator () (int k) {
    return (*get())(k);
  }

  /**
   * @brief Return the state parameter at index `k`.
   */
  const double &operator () (int k) const {
    return (*get())(k);
  }
};

Vector operator += (Vector a, const VectorXd &b);

Vector operator - (const Vector a, const Vector b);

Vector operator - (const VectorXd &a, const Vector b);

Vector operator * (const Vector a, double b);

Vector operator * (double a, const Vector b);

MatrixXd operator * (const Vector a, MatrixXd b);

/**
 * @brief Read a state from the given input stream.
 */
istream &operator >> (istream &in, Vector &x);

/**
 * @brief Write a state to the given output stream.
 */
ostream &operator << (ostream &out, const Vector &x);

} // namespace ukan

#endif
