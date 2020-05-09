#pragma once


#include <stdio.h>
#include <iostream>


class Color;


class Vector3
{
private:
    float _data[3];

public:
    Vector3();
    Vector3(const float &f);
    Vector3(const float &x, const float &y, const float &z);
    Vector3(const Vector3 &v);
    ~Vector3();

    Vector3 operator - () const;
    float operator [] (const int &idx) const;
    float& operator [] (const int &idx);
    Vector3& operator = (const Vector3 &v);
    Vector3& operator += (const Vector3 &v);
    Vector3& operator += (const float &f);
    Vector3& operator -= (const Vector3 &v);
    Vector3& operator -= (const float &f);
    Vector3& operator *= (const Vector3 &v);
    Vector3& operator *= (const float &f);
    Vector3& operator /= (const Vector3 &v);
    Vector3& operator /= (const float &f);

    float x() const;
    float y() const;
    float z() const;
    float& x_();
    float& y_();
    float& z_();
    void setX(const float &f);
    void setY(const float &f);
    void setZ(const float &f);
    void set(const float &x, const float &y, const float &z);
    void set(const float &f);
    float squareMagnitude();
    float magnitude();
    float normalize();
    Vector3 normalized() const;
    float dot(const Vector3 &v) const;
    Vector3 cross(const Vector3 &v) const;
};

Vector3 operator + (const Vector3 &v1, const Vector3 &v2);
Vector3 operator + (const Vector3 &v, const float &f);
Vector3 operator + (const float &f, const Vector3 &v);
Vector3 operator - (const Vector3 &v1, const Vector3 &v2);
Vector3 operator - (const Vector3 &v, const float &f);
Vector3 operator - (const float &f, const Vector3 &v);
Vector3 operator * (const Vector3 &v1, const Vector3 &v2);
Vector3 operator * (const Vector3 &v, const float &f);
Vector3 operator * (const float &f, const Vector3 &v);
Vector3 operator / (const Vector3 &v1, const Vector3 &v2);
Vector3 operator / (const Vector3 &v, const float &f);
Vector3 operator / (const float &f, const Vector3 &v);
std::istream &operator >> (std::istream &s, Vector3 &v);
std::ostream &operator << (std::ostream &s, const Vector3 &v);

typedef Vector3 Point;
