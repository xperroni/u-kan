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

#include "vector.h"

namespace ukan {

namespace base {

Vector::Vector(int n):
  VectorXd(n)
{
  VectorXd &x = *this;
  for (int i = 0; i < n; ++i) {
    x(i) = 0;
  }
}

Vector *Vector::copy() const {
  return new Vector(*this);
}

Vector *Vector::add(const VectorXd &b) const {
  Vector *a = copy();
  a->iadd(b);
  return a;
}

Vector *Vector::mul(double b) const {
  Vector *a = copy();
  a->imul(b);
  return a;
}

void Vector::iadd(const VectorXd &b) {
  *this += b;
}

void Vector::imul(double b) {
  *this *= b;
}

ostream &Vector::write(ostream &out) const {
  for (int i = 0, n = rows(); i < n;) {
    out << (*this)(i);
    if (++i >= n) {
      break;
    }

    out << '\t';
  }

  return out;
}

} // namespace base

Vector operator += (Vector a, const VectorXd &b) {
  a->iadd(b);
  return a;
}

Vector operator - (const Vector a, const Vector b) {
  return a->add(-(*b));
}

Vector operator - (const VectorXd &a, const Vector b) {
  Vector x = b->add(-a);
  *x *= -1;
  return x;
}

Vector operator * (const Vector a, double b) {
  return a->mul(b);
}

Vector operator * (double a, const Vector b) {
  return b->mul(a);
}

MatrixXd operator * (const Vector a, MatrixXd b) {
  return *a * b;
}

istream &operator >> (istream &in, Vector &x) {
  for (int i = 0, n = x->rows(); i < n; ++i) {
    in >> x(i);
  }

  return in;
}

ostream &operator << (ostream &out, const Vector &x) {
  return x->write(out);
}

} // namespace ukan
