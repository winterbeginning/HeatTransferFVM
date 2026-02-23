#ifndef _Point3D_
#define _Point3D_

#include <ostream>
#include <cmath>

template <typename ValueType>
class Tensor3D;

using namespace std;

#define sqr(num) pow(num, 2)

template <typename ValueType>
class Point3D
{
public:
    ValueType x, y, z;
    Point3D() : x(0.0), y(0.0), z(0.0){};
    Point3D(ValueType x, ValueType y, ValueType z) : x(x), y(y), z(z){};

    double getMag()
    {
        return sqrt(sqr(x) + sqr(y) + sqr(z));
    }

    double getDistance(const Point3D& v) const
    {
        return sqrt(sqr(x - v.x) + sqr(y - v.y) + sqr(z - v.z));
    };
    double dotWith(const Point3D& v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
    Point3D crossWith(const Point3D& v) const
    {
        double cx = y * v.z - z * v.y;
        double cy = z * v.x - x * v.z;
        double cz = x * v.y - y * v.x;

        return Point3D(cx, cy, cz);
    }

    Tensor3D<ValueType> outProductWith(const Point3D& v) const
    {
        return Tensor3D(x * v.x,
                        x * v.y,
                        x * v.z,
                        y * v.x,
                        y * v.y,
                        y * v.z,
                        z * v.x,
                        z * v.y,
                        z * v.z);
    }

    friend ostream& operator<<(ostream& out, const Point3D& rhs)
    {
        out << " " << rhs.x << "  " << rhs.y << "  " << rhs.z << "  ";
        return out;
    }

    friend Point3D operator*(const double a, const Point3D& rhs)
    {
        return Point3D(a * rhs.x, a * rhs.y, a * rhs.z);
    }

    Point3D operator*(const double a) const
    {
        return Point3D(a * x, a * y, a * z);
    }

    Point3D operator/(const double a) const
    {
        return Point3D(x / a, y / a, z / a);
    }

    Point3D operator-() const
    {
        return Point3D(-x, -y, -z);
    }

    Point3D operator+(const Point3D& rhs) const
    {
        return Point3D(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    Point3D operator-(const Point3D& rhs) const
    {
        return Point3D(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    Point3D operator*(const Point3D& rhs) const
    {
        return Point3D(x * rhs.x, y * rhs.y, z * rhs.z);
    }

    Point3D operator/(const Point3D& rhs) const
    {
        return Point3D(x / rhs.x, y / rhs.y, z / rhs.z);
    }

    Point3D& operator+=(const double& scalar)
    {
        x += scalar;
        y += scalar;
        z += scalar;
        return *this;
    }
    Point3D& operator+=(const Point3D& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
};

using Vector = Point3D<double>;

#endif