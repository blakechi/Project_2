#include <iostream>
#include <cmath>

#include "Vector3.hpp"


Vector3::Vector3(): _data{0.0f, 0.0f, 0.0f}
{
}

Vector3::Vector3(const float &f): _data{f, f, f}
{
}


Vector3::Vector3(const float &x, const float &y, const float &z): _data{x, y, z}
{
}

Vector3::Vector3(const Vector3 &v): _data{v.x(), v.y(), v.z()}
{
}

Vector3::~Vector3() 
{
}


Vector3 Vector3::operator - () const
{
    return Vector3(-_data[0], -_data[1], -_data[2]);
}

float Vector3::operator [] (const int &idx) const 
{
    return _data[idx];
}

float& Vector3::operator [] (const int &idx)
{
    return _data[idx];
}

Vector3& Vector3::operator = (const Vector3 &v)
{
    _data[0] = v.x();
    _data[1] = v.y();
    _data[2] = v.z();

    return *this;
}

Vector3& Vector3::operator += (const Vector3 &v)
{
    _data[0] += v.x();
    _data[1] += v.y();
    _data[2] += v.z();

    return *this;
}

Vector3& Vector3::operator += (const float &f)
{
    _data[0] += f;
    _data[1] += f;
    _data[2] += f;

    return *this;
}

Vector3& Vector3::operator -= (const Vector3 &v)
{
    _data[0] -= v.x();
    _data[1] -= v.y();
    _data[2] -= v.z();

    return *this;
}

Vector3& Vector3::operator -= (const float &f)
{
    _data[0] -= f;
    _data[1] -= f;
    _data[2] -= f;

    return *this;
}

Vector3& Vector3::operator *= (const Vector3 &v)
{
    _data[0] *= v.x();
    _data[1] *= v.y();
    _data[2] *= v.z();

    return *this;
}

Vector3& Vector3::operator *= (const float &f)
{
    _data[0] *= f;
    _data[1] *= f;
    _data[2] *= f;

    return *this;
}

Vector3& Vector3::operator /= (const Vector3 &v)
{
    _data[0] /= v.x();
    _data[1] /= v.y();
    _data[2] /= v.z();

    return *this;
}

Vector3& Vector3::operator /= (const float &f)
{
    _data[0] /= f;
    _data[1] /= f;
    _data[2] /= f;

    return *this;
}


float Vector3::x() const
{
    return _data[0];
}

float Vector3::y() const
{
    return _data[1];
}

float Vector3::z() const
{
    return _data[2];
}

float& Vector3::x_()
{
    return _data[0];
}

float& Vector3::y_()
{
    return _data[1];
}

float& Vector3::z_()
{
    return _data[2];
}

void Vector3::setX(const float &f)
{
    _data[0] = f;
}

void Vector3::setY(const float &f)
{
    _data[1] = f;
}

void Vector3::setZ(const float &f)
{
    _data[2] = f;
}

void Vector3::set(const float &x, const float &y, const float &z)
{
    _data[0] = x;
    _data[1] = y;
    _data[2] = z;
}

void Vector3::set(const float &f)
{
    _data[0] = f;
    _data[1] = f;
    _data[2] = f;
}

float Vector3::squareMagnitude()
{
    return _data[0]*_data[0] + _data[1]*_data[1] + _data[2]*_data[2];
}

float Vector3::magnitude()
{
    return std::sqrt(squareMagnitude());
}

float Vector3::normalize()
{
    float m = magnitude();
    if(m != 0)
    {
        *this /= m;
    }

    return m;
}

Vector3 Vector3::normalized() const
{
    Vector3 r(*this);
    r.normalize();

    return r;
}

float Vector3::dot(const Vector3 &v) const
{
    return _data[0]*v.x() + _data[1]*v.y() + _data[2]*v.z();
}

Vector3 Vector3::cross(const Vector3 &v) const
{
    return Vector3(
        _data[1]*v.z() - _data[2]*v.y(),
        _data[2]*v.x() - _data[0]*v.z(),
        _data[0]*v.y() - _data[1]*v.x()
    );
}


Vector3 operator + (const Vector3 &v1, const Vector3 &v2)
{
    return Vector3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

Vector3 operator + (const Vector3 &v, const float &f)
{
    return Vector3(v.x() + f, v.y() + f, v.z() + f);
}

Vector3 operator + (const float &f, const Vector3 &v)
{
    return Vector3(v.x() + f, v.y() + f, v.z() + f);
}

Vector3 operator - (const Vector3 &v1, const Vector3 &v2)
{
    return Vector3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

Vector3 operator - (const Vector3 &v, const float &f)
{
    return Vector3(v.x() - f, v.y() - f, v.z() - f);
}

Vector3 operator - (const float &f, const Vector3 &v)
{
    return Vector3(f - v.x(), f - v.y(), f - v.z());
}

Vector3 operator * (const Vector3 &v1, const Vector3 &v2)
{
    return Vector3(v1.x()*v2.x(), v1.y()*v2.y(), v1.z()*v2.z());
}

Vector3 operator * (const Vector3 &v, const float &f)
{
    return Vector3(v.x()*f, v.y()*f, v.z()*f);
}

Vector3 operator * (const float &f, const Vector3 &v)
{
    return Vector3(v.x()*f, v.y()*f, v.z()*f);
}

Vector3 operator / (const Vector3 &v1, const Vector3 &v2) 
{
    return Vector3(v1.x()/v2.x(), v1.y()/v2.y(), v1.z()/v2.z());
}

Vector3 operator / (const Vector3 &v, const float &f) 
{
    return Vector3(v.x()/f, v.y()/f, v.z()/f);
}

Vector3 operator / (const float &f, const Vector3 &v)
{
    return Vector3(f/v.x(), f/v.y(), f/v.z());
}

std::istream &operator >> (std::istream &s, Vector3 &v)
{
    s >> v.x_() >> v.y_() >> v.z_();

    return s;
}

std::ostream &operator << (std::ostream &s, const Vector3 &v)
{
    s << v.x() << " " << v.y() << " " << v.z() << std::endl;

    return s;
}