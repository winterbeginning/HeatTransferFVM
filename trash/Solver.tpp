#ifndef _Solver_TPP_
#define _Solver_TPP_

#include <iostream>
#include <type_traits>

// ============ 辅助函数：提取标量值 ============
template <typename T>
inline double extractScalar(const T& val)
{
    if constexpr (std::is_same_v<T, double>)
        return val;
    else
        return val.getMag(); // 对 Vector/Tensor 调用 getMag()
}

// ============ solve 主函数 ============
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
        return gaussSeidel(Eqn, x0);
    }
}

// ============ dot 函数 ============
template <typename ValueType>
inline ValueType Solver<ValueType>::dot(const std::vector<ValueType>& v1,
                                        const std::vector<ValueType>& v2)
{
    ValueType res{};
    for (size_t i = 0; i < v1.size(); ++i)
        res = res + v1[i] * v2[i];
    return res;
}

// ============ norm 函数 ============
template <typename ValueType>
double Solver<ValueType>::norm(const std::vector<ValueType>& v)
{
    auto d = dot(v, v);
    return std::sqrt(extractScalar(d)); // 使用辅助函数
}

// ============ Jacobi ============
template <typename ValueType>
std::vector<ValueType>
Solver<ValueType>::jacobi(const SpaceMatrix<ValueType>& Eqn,
                          const std::vector<ValueType>& x0)
{
    std::vector<ValueType> x = (x0.size() == (size_t)Eqn.n)
                                   ? x0
                                   : std::vector<ValueType>(Eqn.n, ValueType{});
    std::vector<ValueType> x_new(Eqn.n, ValueType{});

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
            x_new[i] = (Eqn.b[i] - sum) / diag;
        }
        x = x_new;

        auto r = Eqn.multiply(x);
        for (int i = 0; i < Eqn.n; ++i)
            r[i] = r[i] - Eqn.b[i];

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

// ============ Gauss-Seidel ============
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
            x[i] = (Eqn.b[i] - sum) / diag;
        }

        auto r = Eqn.multiply(x);
        for (int i = 0; i < Eqn.n; ++i)
            r[i] = r[i] - Eqn.b[i];

        double res = norm(r);
        if (verbose && iter % 100 == 0)
            std::cout << "[Gauss-Seidel] Iter " << iter << ", Res: " << res
                      << std::endl;

        if (res < tol)
        {
            if (verbose)
                std::cout << "[Gauss-Seidel] converged in " << iter
                          << " iterations\n";
            return x;
        }
    }
    return x;
}

// ============ Conjugate Gradient ============
template <typename ValueType>
std::vector<ValueType>
Solver<ValueType>::conjugateGradient(const SpaceMatrix<ValueType>& Eqn,
                                     const std::vector<ValueType>& x0)
{
    std::vector<ValueType> x = (x0.size() == (size_t)Eqn.n)
                                   ? x0
                                   : std::vector<ValueType>(Eqn.n, ValueType{});
    std::vector<ValueType> r = Eqn.b;
    std::vector<ValueType> Ax = Eqn.multiply(x);
    for (int i = 0; i < Eqn.n; ++i)
        r[i] = r[i] - Ax[i];
    std::vector<ValueType> p = r;

    double rsold = extractScalar(dot(r, r)); // 使用辅助函数

    for (int iter = 0; iter < maxIter; ++iter)
    {
        std::vector<ValueType> Ap = Eqn.multiply(p);
        double check_denom = extractScalar(dot(p, Ap)); // 使用辅助函数

        if (std::abs(check_denom) < 1e-20)
            break;

        double alpha = rsold / check_denom;
        for (int i = 0; i < Eqn.n; ++i)
        {
            x[i] = x[i] + p[i] * alpha;
            r[i] = r[i] - Ap[i] * alpha;
        }

        double rsnew = extractScalar(dot(r, r)); // 使用辅助函数
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
            p[i] = r[i] + p[i] * (rsnew / rsold);
        }
        rsold = rsnew;
    }
    return x;
}

// ============ BiCGSTAB ============
template <typename ValueType>
std::vector<ValueType>
Solver<ValueType>::bicgstab(const SpaceMatrix<ValueType>& Eqn,
                            const std::vector<ValueType>& x0)
{
    std::vector<ValueType> x = (x0.size() == (size_t)Eqn.n)
                                   ? x0
                                   : std::vector<ValueType>(Eqn.n, ValueType{});
    std::vector<ValueType> r = Eqn.b;
    std::vector<ValueType> Ax = Eqn.multiply(x);
    for (int i = 0; i < Eqn.n; ++i)
        r[i] = r[i] - Ax[i];
    std::vector<ValueType> r_star = r;
    std::vector<ValueType> p = r;

    double rho = 1.0, alpha = 1.0, omega = 1.0;
    std::vector<ValueType> v(Eqn.n, ValueType{}), s(Eqn.n, ValueType{});

    for (int iter = 0; iter < maxIter; ++iter)
    {
        double rho_prev = rho;
        rho = extractScalar(dot(r_star, r)); // 使用辅助函数

        if (std::abs(rho) < 1e-20)
            break;

        if (iter > 0)
        {
            double beta = (rho / rho_prev) * (alpha / omega);
            for (int i = 0; i < Eqn.n; ++i)
            {
                p[i] = r[i] + (p[i] - v[i] * omega) * beta;
            }
        }

        v = Eqn.multiply(p);
        double denom = extractScalar(dot(r_star, v)); // 使用辅助函数

        if (std::abs(denom) < 1e-20)
            break;

        alpha = rho / denom;
        for (int i = 0; i < Eqn.n; ++i)
            s[i] = r[i] - v[i] * alpha;

        double res_s = norm(s);
        if (res_s < tol)
        {
            for (int i = 0; i < Eqn.n; ++i)
                x[i] = x[i] + p[i] * alpha;
            if (verbose)
                std::cout << "[BiCGSTAB] (early) converged in " << iter
                          << " iterations, Final Res: " << res_s << std::endl;
            return x;
        }

        std::vector<ValueType> t = Eqn.multiply(s);
        double t_norm_sq = extractScalar(dot(t, t)); // 使用辅助函数

        if (std::abs(t_norm_sq) < 1e-20)
            break;

        omega = extractScalar(dot(t, s)) / t_norm_sq; // 使用辅助函数

        for (int i = 0; i < Eqn.n; ++i)
        {
            x[i] = x[i] + p[i] * alpha + s[i] * omega;
            r[i] = s[i] - t[i] * omega;
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

#endif