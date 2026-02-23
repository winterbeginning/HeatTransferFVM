#ifndef _Tensor_
#define _Tensor_

#include <ostream>
#include <cmath>
#include "Point3D.hpp"

using namespace std;

#ifndef sqr
#define sqr(num) pow(num, 2)
#endif

template <typename ValueType>
class Tensor3D
{
    // 3x3 张量的9个分量（按行存储：xx xy xz; yx yy yz; zx zy zz）
    ValueType xx, xy, xz;
    ValueType yx, yy, yz;
    ValueType zx, zy, zz;

    // 构造函数
    // 1. 默认构造：零张量
    Tensor3D() : xx(0), xy(0), xz(0), yx(0), yy(0), yz(0), zx(0), zy(0), zz(0)
    {
    }

    // 2. 全参数构造
    Tensor3D(ValueType xx,
             ValueType xy,
             ValueType xz,
             ValueType yx,
             ValueType yy,
             ValueType yz,
             ValueType zx,
             ValueType zy,
             ValueType zz)
        : xx(xx), xy(xy), xz(xz), yx(yx), yy(yy), yz(yz), zx(zx), zy(zy), zz(zz)
    {
    }

    // 3. 从对角分量构造（对角张量，如单位张量、标量张量）
    Tensor3D(ValueType diag)
        : xx(diag), xy(0), xz(0), yx(0), yy(diag), yz(0), zx(0), zy(0), zz(diag)
    {
    }

    // ========== 核心数学操作 ==========
    // 1. 转置
    Tensor3D transpose() const
    {
        return Tensor3D(xx, yx, zx, xy, yy, zy, xz, yz, zz);
    }

    // 2. 迹（对角元素和）
    ValueType trace() const
    {
        return xx + yy + zz;
    }

    // 3. 行列式（3x3张量行列式）
    double determinant() const
    {
        return xx * (yy * zz - yz * zy) - xy * (yx * zz - yz * zx) +
               xz * (yx * zy - yy * zx);
    }

    // 4. 张量-矢量乘法（核心：张量作用于矢量）
    Point3D<ValueType> multiply(const Point3D<ValueType>& v) const
    {
        ValueType x = xx * v.x + xy * v.y + xz * v.z;
        ValueType y = yx * v.x + yy * v.y + yz * v.z;
        ValueType z = zx * v.x + zy * v.y + zz * v.z;
        return Point3D<ValueType>(x, y, z);
    }

    // 5. 张量-张量乘法
    Tensor3D multiply(const Tensor3D& t) const
    {
        return Tensor3D(xx * t.xx + xy * t.yx + xz * t.zx,
                        xx * t.xy + xy * t.yy + xz * t.zy,
                        xx * t.xz + xy * t.yz + xz * t.zz,
                        yx * t.xx + yy * t.yx + yz * t.zx,
                        yx * t.xy + yy * t.yy + yz * t.zy,
                        yx * t.xz + yy * t.yz + yz * t.zz,
                        zx * t.xx + zy * t.yx + zz * t.zx,
                        zx * t.xy + zy * t.yy + zz * t.zy,
                        zx * t.xz + zy * t.yz + zz * t.zz);
    }

    // 6. 张量加法
    Tensor3D operator+(const Tensor3D& rhs) const
    {
        return Tensor3D(xx + rhs.xx,
                        xy + rhs.xy,
                        xz + rhs.xz,
                        yx + rhs.yx,
                        yy + rhs.yy,
                        yz + rhs.yz,
                        zx + rhs.zx,
                        zy + rhs.zy,
                        zz + rhs.zz);
    }

    // 7. 张量减法
    Tensor3D operator-(const Tensor3D& rhs) const
    {
        return Tensor3D(xx - rhs.xx,
                        xy - rhs.xy,
                        xz - rhs.xz,
                        yx - rhs.yx,
                        yy - rhs.yy,
                        yz - rhs.yz,
                        zx - rhs.zx,
                        zy - rhs.zy,
                        zz - rhs.zz);
    }

    // 8. 标量乘法（张量*标量）
    Tensor3D operator*(const double a) const
    {
        return Tensor3D(a * xx,
                        a * xy,
                        a * xz,
                        a * yx,
                        a * yy,
                        a * yz,
                        a * zx,
                        a * zy,
                        a * zz);
    }

    // 9. 标量除法
    Tensor3D operator/(const double a) const
    {
        return Tensor3D(xx / a,
                        xy / a,
                        xz / a,
                        yx / a,
                        yy / a,
                        yz / a,
                        zx / a,
                        zy / a,
                        zz / a);
    }

    // ========== 友元函数 ==========
    // 标量*张量（左乘）
    friend Tensor3D operator*(const double a, const Tensor3D& rhs)
    {
        return rhs * a; // 复用张量*标量
    }

    // 张量-矢量乘法的运算符重载（更直观：tensor * vector）
    friend Point3D<ValueType> operator*(const Tensor3D& t,
                                        const Point3D<ValueType>& v)
    {
        return t.multiply(v);
    }

    // 输出运算符重载
    friend std::ostream& operator<<(std::ostream& out, const Tensor3D& rhs)
    {
        out << "[" << rhs.xx << ", " << rhs.xy << ", " << rhs.xz << "]\n"
            << "[" << rhs.yx << ", " << rhs.yy << ", " << rhs.yz << "]\n"
            << "[" << rhs.zx << ", " << rhs.zy << ", " << rhs.zz << "]";
        return out;
    }

    // ========== 常用张量构造 ==========
    // 单位张量
    static Tensor3D identity()
    {
        return Tensor3D(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    }

    // 从矢量外积构造张量（v1 ⊗ v2）
    static Tensor3D outerProduct(const Point3D<ValueType>& v1,
                                 const Point3D<ValueType>& v2)
    {
        return Tensor3D(v1.x * v2.x,
                        v1.x * v2.y,
                        v1.x * v2.z,
                        v1.y * v2.x,
                        v1.y * v2.y,
                        v1.y * v2.z,
                        v1.z * v2.x,
                        v1.z * v2.y,
                        v1.z * v2.z);
    }
};

// 定义常用类型别名（和 Vector 对应）
typedef Tensor3D<double> Tensor;

#endif