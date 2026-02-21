#include "Solver.hpp"
#include <iostream>

std::vector<double> Solver::solve(SolverType type,
                                  const SpaceMatrix& A,
                                  const std::vector<double>& b,
                                  const std::vector<double>& x0,
                                  int maxIter,
                                  double tol)
{
    switch (type)
    {
    case SolverType::JACOBI:
        return jacobi(A, b, x0, maxIter, tol);
    case SolverType::GAUSS_SEIDEL:
        return gaussSeidel(A, b, x0, maxIter, tol);
    case SolverType::CG:
        return conjugateGradient(A, b, x0, maxIter, tol);
    case SolverType::BICGSTAB:
        return bicgstab(A, b, x0, maxIter, tol);
    default:
        return conjugateGradient(A, b, x0, maxIter, tol);
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

std::vector<double> Solver::jacobi(const SpaceMatrix& A,
                                   const std::vector<double>& b,
                                   const std::vector<double>& x0,
                                   int maxIter,
                                   double tol)
{
    std::vector<double> x =
        (x0.size() == A.n) ? x0 : std::vector<double>(A.n, 0.0);
    std::vector<double> x_new(A.n, 0.0);
    for (int iter = 0; iter < maxIter; ++iter)
    {
        for (int i = 0; i < A.n; ++i)
        {
            double diag = A.diag(i);
            double sum = 0.0;
            for (const auto& [j, val] : A.data[i])
            {
                if (j != i)
                    sum += val * x[j];
            }
            x_new[i] = (b[i] - sum) / diag;
        }
        x = x_new;
        std::vector<double> r = A.multiply(x);
        for (int i = 0; i < A.n; ++i)
            r[i] -= b[i];
        if (norm(r) < tol)
        {
            std::cout << "Jacobi converged in " << iter << " iterations."
                      << std::endl;
            return x;
        }
    }
    return x;
}

std::vector<double> Solver::gaussSeidel(const SpaceMatrix& A,
                                        const std::vector<double>& b,
                                        const std::vector<double>& x0,
                                        int maxIter,
                                        double tol)
{
    std::vector<double> x =
        (x0.size() == A.n) ? x0 : std::vector<double>(A.n, 0.0);
    for (int iter = 0; iter < maxIter; ++iter)
    {
        for (int i = 0; i < A.n; ++i)
        {
            double diag = A.diag(i);
            double sum = 0.0;
            for (const auto& [j, val] : A.data[i])
            {
                if (j != i)
                    sum += val * x[j];
            }
            x[i] = (b[i] - sum) / diag;
        }
        std::vector<double> r = A.multiply(x);
        for (int i = 0; i < A.n; ++i)
            r[i] -= b[i];
        if (norm(r) < tol)
        {
            std::cout << "Gauss-Seidel converged in " << iter << " iterations."
                      << std::endl;
            return x;
        }
    }
    return x;
}

std::vector<double> Solver::conjugateGradient(const SpaceMatrix& A,
                                              const std::vector<double>& b,
                                              const std::vector<double>& x0,
                                              int maxIter,
                                              double tol)
{
    std::vector<double> x =
        (x0.size() == A.n) ? x0 : std::vector<double>(A.n, 0.0);
    std::vector<double> r = b;
    std::vector<double> Ax = A.multiply(x);
    for (int i = 0; i < A.n; ++i)
        r[i] -= Ax[i];
    std::vector<double> p = r;
    double rsold = dot(r, r);

    for (int iter = 0; iter < maxIter; ++iter)
    {
        std::vector<double> Ap = A.multiply(p);
        double check_denom = dot(p, Ap);
        if (std::abs(check_denom) < 1e-20)
            break;

        double alpha = rsold / check_denom;
        for (int i = 0; i < A.n; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        double rsnew = dot(r, r);
        if (std::sqrt(rsnew) < tol)
        {
            std::cout << "CG converged in " << iter << " iterations."
                      << std::endl;
            return x;
        }
        for (int i = 0; i < A.n; ++i)
        {
            p[i] = r[i] + (rsnew / rsold) * p[i];
        }
        rsold = rsnew;
    }
    return x;
}

std::vector<double> Solver::bicgstab(const SpaceMatrix& A,
                                     const std::vector<double>& b,
                                     const std::vector<double>& x0,
                                     int maxIter,
                                     double tol)
{
    std::vector<double> x =
        (x0.size() == A.n) ? x0 : std::vector<double>(A.n, 0.0);
    std::vector<double> r = b;
    std::vector<double> Ax = A.multiply(x);
    for (int i = 0; i < A.n; ++i)
        r[i] -= Ax[i];
    std::vector<double> r_star = r;
    std::vector<double> p = r;

    double rho = 1.0, alpha = 1.0, omega = 1.0;
    std::vector<double> v(A.n, 0.0), s(A.n, 0.0);

    for (int iter = 0; iter < maxIter; ++iter)
    {
        double rho_prev = rho;
        rho = dot(r_star, r);
        if (std::abs(rho) < 1e-20)
            break;

        if (iter > 0)
        {
            double beta = (rho / rho_prev) * (alpha / omega);
            for (int i = 0; i < A.n; ++i)
            {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }

        v = A.multiply(p);
        double denom = dot(r_star, v);
        if (std::abs(denom) < 1e-20)
            break;

        alpha = rho / denom;
        for (int i = 0; i < A.n; ++i)
            s[i] = r[i] - alpha * v[i];

        if (norm(s) < tol)
        {
            for (int i = 0; i < A.n; ++i)
                x[i] += alpha * p[i];
            std::cout << "BiCGSTAB (early) converged in " << iter
                      << " iterations." << std::endl;
            return x;
        }

        std::vector<double> t = A.multiply(s);
        double t_norm_sq = dot(t, t);
        if (std::abs(t_norm_sq) < 1e-20)
            break;

        omega = dot(t, s) / t_norm_sq;
        for (int i = 0; i < A.n; ++i)
        {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }

        if (norm(r) < tol)
        {
            std::cout << "BiCGSTAB converged in " << iter << " iterations."
                      << std::endl;
            return x;
        }
    }
    return x;
}
