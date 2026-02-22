#include "Solver.hpp"
#include <iostream>

std::vector<double> Solver::solve(SpaceMatrix& Eqn,
                                  const std::vector<double>& x0)
{
    Eqn.compress(); // 求解前将组装好的 A 压缩为 CSR
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

double Solver::dot(const std::vector<double>& v1, const std::vector<double>& v2)
{
    double res = 0;
    for (size_t i = 0; i < v1.size(); ++i)
        res += v1[i] * v2[i];
    return res;
}

double Solver::norm(const std::vector<double>& v)
{
    return std::sqrt(dot(v, v));
}

std::vector<double> Solver::jacobi(const SpaceMatrix& Eqn,
                                   const std::vector<double>& x0)
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
            x_new[i] = (Eqn.b[i] - sum) / diag;
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

std::vector<double> Solver::gaussSeidel(const SpaceMatrix& Eqn,
                                        const std::vector<double>& x0)
{
    std::vector<double> x =
        (x0.size() == (size_t)Eqn.n) ? x0 : std::vector<double>(Eqn.n, 0.0);
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
            x[i] = (Eqn.b[i] - sum) / diag;
        }
        std::vector<double> r = Eqn.multiply(x);
        for (int i = 0; i < Eqn.n; ++i)
            r[i] -= Eqn.b[i];
        double res = norm(r);
        if (verbose && iter % 100 == 0)
            std::cout << "[Gauss-Seidel] Iter " << iter << ", Residual: " << res
                      << std::endl;
        if (res < tol)
        {
            if (verbose)
                std::cout << "[Gauss-Seidel] converged in " << iter
                          << " iterations, Final Res: " << res << std::endl;
            return x;
        }
    }
    return x;
}

std::vector<double> Solver::conjugateGradient(const SpaceMatrix& Eqn,
                                              const std::vector<double>& x0)
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

std::vector<double> Solver::bicgstab(const SpaceMatrix& Eqn,
                                     const std::vector<double>& x0)
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
