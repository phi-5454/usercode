#ifndef _Vectors_h_
#define _Vectors_h_

#include <cmath>
#include <iostream>

//
struct Vector1
{
  double val; // value
  double sig; // sigma
};

//
struct Vector2
{
  double x,y;         // value
  double cxx,cxy,cyy; // covariance

  Vector2 operator+=(const Vector2 & a) // increment
  { 
    this->x += a.x;
    this->y += a.y;

    this->cxx += a.cxx;
    this->cxy += a.cxy;
    this->cyy += a.cyy;

    return *this;
  }
};

//
struct Vector3
{
  double x,y,z;

  Vector3 operator+(const Vector3 & a) const // add
  { return {x + a.x, y + a.y, z + a.z}; }
  Vector3 operator-(const Vector3 & a) const // subst
  { return {x - a.x, y - a.y, z - a.z}; }

  Vector3 operator+=(const Vector3 & a) // increment
  { this->x += a.x; this->y += a.y; this->z += a.z; return *this; }
  Vector3 operator-=(const Vector3 & a) // decrement
  { this->x -= a.x; this->y -= a.y; this->z -= a.z; return *this; }

  double operator*(const Vector3 & a) const // scalar
  { return {x*a.x + y*a.y + z*a.z}; }
  Vector3 operator^(const Vector3 & a) const // cross
  { return {y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x}; }

  Vector3 operator*(double a) const
  { return {x*a, y*a, z*a}; }

  //
  double phi() const
  { return atan2(this->y, this->x); }

  double dphi(const Vector3 & a) const
  {
    // https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors
    double dot = this->x*a.x + this->y*a.y;
    double det = this->x*a.y - this->y*a.x;

    return fabs(atan2(det, dot));
/*
    double s = (this->x * a.x + this->y * a.y) /
               (this->trans() * a.trans());

    if(s >= 1.) return 0;
    if(s <=-1.) return M_PI;
    return acos(s);
*/
  }

  //
  double trans() const
  { return sqrt(x*x + y*y); }

  double trans2() const
  { return      x*x + y*y; }

  double length() const
  { return sqrt((*this)*(*this)); }

  double length2() const
  { return     ((*this)*(*this)); }

  Vector3 norm() const
  { return ((*this) * (1/this->length())); }

  //
  double eta() const
  { return 1./2 * log( (this->length() + this->z)
                     / (this->length() - this->z) ); }

//  void print() const
//  { std::cerr << " p = ( " << x << ", " << y << ", " << z << ")" << std::endl; }
};

//
typedef std::pair<double,Vector3> Boost;

//
struct Vector4
{
  double E;
  Vector3 p;

  Vector4 operator+(const Vector4 & a) const // add
  { return {E + a.E, p + a.p}; }
  Vector4 operator-(const Vector4 & a) const // subst
  { return {E - a.E, p - a.p}; }

  double operator*(const Vector4 & a) const // scalar
  { return {E*a.E - p*a.p}; }

  double mass() const
  { return sqrt((*this)*(*this)); }

  double mass2() const
  { return ((*this)*(*this)); }

  double rap() const
  { return 1./2 * log( (this->E + this->p.z)
                      /(this->E - this->p.z) ); }

  static Vector4 fourVector(double m, const Vector3 & p)
  {
    double E = sqrt(m*m + p.length2());

    Vector4 v4; v4.E = E; v4.p = p; return v4;
    // return {E, p.x, p.y, p.z};
  }

  static Boost getBoost(const Vector4 & a)
  {
    double beta = a.p.length() / a.E;
    Vector3 n   = a.p.norm();
    Boost b = std::pair<double,Vector3>(beta,n);

    return b;
  }

  Vector4 operator^(const Boost & b) // boost
  {
    const double & beta = b.first;
    const double   gamm = 1/sqrt(1 - beta*beta);
    const Vector3 & n   = b.second;

    double  & E = this->E;
    Vector3 & p = this->p;

    double pl = n*p;

    double E_  = gamm*( E - beta*pl);
    double pl_ = gamm*(-beta*E + pl);

    this->E = E_;
    this->p = (n*pl_) + (p - n*pl); // pl + pt

    return *this;
  }

  void print() const
  { std::cerr << " m = " << mass()
              << " p = ( " << p.x << ", " << p.y << ", " << p.z << ")" << std::endl; }
};

#endif
