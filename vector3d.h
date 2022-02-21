/*
 * vector3d.h
 *
 *  Created on: 6 апр. 2021 г.
 *      Author: User
 */
#pragma once
#ifndef NAVIGINE_VECTOR3D_H
#define NAVIGINE_VECTOR3D_H

#include <cmath>
static const double DOUBLE_EPSILON = 1e-8;

class Vector3d
{
public:
  double x;
  double y;
  double z;
  Vector3d():x(0),y(0),z(0) {}
  Vector3d(double _x, double _y, double _z);
  ~Vector3d() {}
  double    length()  const;
  double    norm() const;
  Vector3d  normalized()  const;
  Vector3d& normalize();
  bool      isNull() const;


  static Vector3d crossProduct(const Vector3d& v1, const Vector3d& v2);
  static double   dotProduct  (const Vector3d& v1, const Vector3d& v2);
  static Vector3d rotate      (const Vector3d& r, const Vector3d& axis, double fi);

  Vector3d& operator+=(const Vector3d& v2);
  Vector3d& operator-=(const Vector3d& v2);
  Vector3d& operator*=(double multiplier);
  Vector3d& operator/=(double divisor);
};


inline bool operator==(const Vector3d& v1, const Vector3d& v2)
{
  return std::fabs(v1.x - v2.x) < DOUBLE_EPSILON &&
         std::fabs(v1.y - v2.y) < DOUBLE_EPSILON &&
         std::fabs(v1.z - v2.z) < DOUBLE_EPSILON;
}

inline bool operator!=(const Vector3d& v1, const Vector3d& v2)
{
  return !(v1 == v2);
}

inline Vector3d operator*(double multiplier, const Vector3d& v)
{
  return Vector3d(v.x * multiplier, v.y * multiplier, v.z * multiplier);
}

inline Vector3d operator*(const Vector3d& v, double multiplier)
{
  return Vector3d(v.x * multiplier, v.y * multiplier, v.z * multiplier);
}

inline Vector3d operator+(const Vector3d& v1, const Vector3d& v2)
{
  return Vector3d(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector3d operator-(const Vector3d& v1, const Vector3d& v2)
{
  return Vector3d(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector3d operator-(const Vector3d& v)
{
  return Vector3d(-v.x, -v.y, -v.z);
}

inline Vector3d operator/(const Vector3d& v1, double divisor)
{
  return Vector3d(v1.x / divisor, v1.y / divisor, v1.z / divisor);
}

inline Vector3d operator*(const double Matr[3][3], const Vector3d& v)
{
  double x = Matr[0][0] * v.x + Matr[0][1] * v.y + Matr[0][2] * v.z;
  double y = Matr[1][0] * v.x + Matr[1][1] * v.y + Matr[1][2] * v.z;
  double z = Matr[2][0] * v.x + Matr[2][1] * v.y + Matr[2][2] * v.z;
  return Vector3d(x,y,z);
}

#endif // NAVIGINE_VECTOR3D_H
