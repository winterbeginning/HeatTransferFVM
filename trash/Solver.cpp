#include "Solver.hpp"
#include <iostream>

template <typename ValueType>
std::vector<ValueType>
Solver<ValueType>::solve(SpaceMatrix<ValueType>& Eqn,
                         const std::vector<ValueType>& x0)
{
    Eqn.compress();
    switch (solverType)
    {
    case SolverType::JACOBI:
        return jacobi(Eqn, x0);
    case SolverType::GAUSS_SEIDEL:
        return gaussSeidel(Eqn, x0);
    case SolverType::CG:
        return conjugateGradient(Eqn, x0);
    case SolverType::BICGSTAB:
        return bicgstab(Eqn, x0);
    default:
        return conjugateGradient(Eqn, x0);
    }
}

// dot 函数需要针对不同类型特化
template <typename ValueType>
inline ValueType Solver<ValueType>::dot(const std::vector<ValueType>& v1,
                                        const std::vector<ValueType>& v2)
{
    ValueType res = 0;
    for (size_t i = 0; i < v1.size(); ++i)
        res += v1[i] * v2[i];
    return res;
}

// template <>
// inline Vector Solver<Vector>::dot(const std::vector<Vector>& v1,
//                                   const std::vector<Vector>& v2)
// {
//     Vector res(0, 0, 0);
//     for (size_t i = 0; i < v1.size(); ++i)
//         res = res + v1[i] % v2[i]; // 分量乘法
//     return res;
// }

// norm 对所有类型通用
template <typename ValueType>
double Solver<ValueType>::norm(const std::vector<ValueType>& v)
{
    auto d = dot(v, v);
    // 对 double: 直接 sqrt
    // 对 Vector: 需要先取模
    if constexpr (std::is_same_v<ValueType, double>)
        return std::sqrt(d);
    else
        return std::sqrt(d.getMag()); // Vector 的模
}

template <typename ValueType>
std::vector<ValueType>
Solver<ValueType>::jacobi(const SpaceMatrix<ValueType>& Eqn,
                          const std::vector<ValueType>& x0)
{
    std::vector<double> x =
        (x0.size() == (size_t)Eqn.n) ? x0 : std::vector<double>(Eqn.n, 0.0);
    std::vector<double> x_new(Eqn.n, 0.0);
    for (int iter = 0; iter < maxIter; ++iter)
    {
        for (int i = 0; i < Eqn.n; ++i)
        {
            double diag = 0.0;
            double sum = 0.0;
            for (int k = Eqn.csr_ptr[i]; k < Eqn.csr_ptr[i + 1]; ++k)
            {
                if (Eqn.csr_col[k] == i)
                    diag = Eqn.csr_val[k];
                else
                    sum += Eqn.csr_val[k] * x[Eqn.csr_col[k]];
            }
            x_new[i] = (Eqn.b[i] - sum) / diag; // Vector 需要重载 operator/
        }
        x = x_new;
        std::vector<double> r = Eqn.multiply(x);
        for (int i = 0; i < Eqn.n; ++i)
            r[i] -= Eqn.b[i];
        double res = norm(r);
        if (verbose && iter % 100 == 0)
            std::cout << "[Jacobi] Iter " << iter << ", Residual: " << res
                      << std::endl;
        if (res < tol)
        {
            if (verbose)
                std::cout << "[Jacobi] converged in " << iter
                          << " iterations, Final Res: " << res << std::endl;
            return x;
        }
    }
    return x;
}

// Gauss-Seidel (通用模板)
template <typename ValueType>
std::vector<ValueType>
Solver<ValueType>::gaussSeidel(const SpaceMatrix<ValueType>& Eqn,
                               const std::vector<ValueType>& x0)
{
    std::vector<ValueType> x = (x0.size() == (size_t)Eqn.n)
                                   ? x0
                                   : std::vector<ValueType>(Eqn.n, ValueType{});

    for (int iter = 0; iter < maxIter; ++iter)
    {
        for (int i = 0; i < Eqn.n; ++i)
        {
            ValueType diag{};
            ValueType sum{};

            for (int k = Eqn.csr_ptr[i]; k < Eqn.csr_ptr[i + 1]; ++k)
            {
                if (Eqn.csr_col[k] == i)
                    diag = Eqn.csr_val[k];
                else
                    sum = sum + Eqn.csr_val[k] * x[Eqn.csr_col[k]];
            }
            x[i] = (Eqn.b[i] - sum) / diag; // Vector 需要重载 operator/
        }

        auto r = Eqn.multiply(x);
        for (int i = 0; i < Eqn.n; ++i)
            r[i] = r[i] - Eqn.b[i];

        double res = norm(r);
        if (verbose && iter % 100 == 0)
            std::cout << "[GS] Iter " << iter << ", Res: " << res << std::endl;

        if (res < tol)
        {
            if (verbose)
                std::cout << "[GS] converged in " << iter << " iterations\n";
            return x;
        }
    }
    return x;
}

template <typename ValueType>
std::vector<ValueType>
Solver<ValueType>::conjugateGradient(const SpaceMatrix<ValueType>& Eqn,
                                     const std::vector<ValueType>& x0)
{
    std::vector<double> x =
        (x0.size() == (size_t)Eqn.n) ? x0 : std::vector<double>(Eqn.n, 0.0);
    std::vector<double> r = Eqn.b;
    std::vector<double> Ax = Eqn.multiply(x);
    for (int i = 0; i < Eqn.n; ++i)
        r[i] -= Ax[i];
    std::vector<double> p = r;
    double rsold = dot(r, r);

    for (int iter = 0; iter < maxIter; ++iter)
    {
        std::vector<double> Ap = Eqn.multiply(p);
        double check_denom = dot(p, Ap);
        if (std::abs(check_denom) < 1e-20)
            break;

        double alpha = rsold / check_denom;
        for (int i = 0; i < Eqn.n; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        double rsnew = dot(r, r);
        double res = std::sqrt(rsnew);
        if (verbose && iter % 100 == 0)
            std::cout << "[CG] Iter " << iter << ", Residual: " << res
                      << std::endl;
        if (res < tol)
        {
            if (verbose)
                std::cout << "[CG] converged in " << iter
                          << " iterations, Final Res: " << res << std::endl;
            return x;
        }
        for (int i = 0; i < Eqn.n; ++i)
        {
            p[i] = r[i] + (rsnew / rsold) * p[i];
        }
        rsold = rsnew;
    }
    return x;
}

template <typename ValueType>
std::vector<ValueType>
Solver<ValueType>::bicgstab(const SpaceMatrix<ValueType>& Eqn,
                            const std::vector<ValueType>& x0)
{
    std::vector<double> x =
        (x0.size() == (size_t)Eqn.n) ? x0 : std::vector<double>(Eqn.n, 0.0);
    std::vector<double> r = Eqn.b;
    std::vector<double> Ax = Eqn.multiply(x);
    for (int i = 0; i < Eqn.n; ++i)
        r[i] -= Ax[i];
    std::vector<double> r_star = r;
    std::vector<double> p = r;

    double rho = 1.0, alpha = 1.0, omega = 1.0;
    std::vector<double> v(Eqn.n, 0.0), s(Eqn.n, 0.0);

    for (int iter = 0; iter < maxIter; ++iter)
    {
        double rho_prev = rho;
        rho = dot(r_star, r);
        if (std::abs(rho) < 1e-20)
            break;

        if (iter > 0)
        {
            double beta = (rho / rho_prev) * (alpha / omega);
            for (int i = 0; i < Eqn.n; ++i)
            {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }

        v = Eqn.multiply(p);
        double denom = dot(r_star, v);
        if (std::abs(denom) < 1e-20)
            break;

        alpha = rho / denom;
        for (int i = 0; i < Eqn.n; ++i)
            s[i] = r[i] - alpha * v[i];

        double res_s = norm(s);
        if (res_s < tol)
        {
            for (int i = 0; i < Eqn.n; ++i)
                x[i] += alpha * p[i];
            if (verbose)
                std::cout << "[BiCGSTAB] (early) converged in " << iter
                          << " iterations, Final Res: " << res_s << std::endl;
            return x;
        }

        std::vector<double> t = Eqn.multiply(s);
        double t_norm_sq = dot(t, t);
        if (std::abs(t_norm_sq) < 1e-20)
            break;

        omega = dot(t, s) / t_norm_sq;
        for (int i = 0; i < Eqn.n; ++i)
        {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }

        double res = norm(r);
        if (verbose && iter % 100 == 0)
            std::cout << "[BiCGSTAB] Iter " << iter << ", Residual: " << res
                      << std::endl;
        if (res < tol)
        {
            if (verbose)
                std::cout << "[BiCGSTAB] converged in " << iter
                          << " iterations, Final Res: " << res << std::endl;

            return x;
        }
    }
    return x;
}

template class Solver<double>;
