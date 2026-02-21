#pragma once

#include <vector>
#include <cmath>
#include "SpaceMatrix.hpp"

enum class SolverType
{
    JACOBI,
    GAUSS_SEIDEL,
    CG,
    BICGSTAB
};

class Solver
{
public:
    static std::vector<double> solve(SolverType type,
                                     const SpaceMatrix& A,
                                     const std::vector<double>& b,
                                     const std::vector<double>& x0 = {},
                                     int maxIter = 1000,
                                     double tol = 1e-6);
    static std::vector<double> jacobi(const SpaceMatrix& A,
                                      const std::vector<double>& b,
                                      const std::vector<double>& x0 = {},
                                      int maxIter = 1000,
                                      double tol = 1e-6);
    static std::vector<double> gaussSeidel(const SpaceMatrix& A,
                                           const std::vector<double>& b,
                                           const std::vector<double>& x0 = {},
                                           int maxIter = 1000,
                                           double tol = 1e-6);
    static std::vector<double>
    conjugateGradient(const SpaceMatrix& A,
                      const std::vector<double>& b,
                      const std::vector<double>& x0 = {},
                      int maxIter = 1000,
                      double tol = 1e-6);
    static std::vector<double> bicgstab(const SpaceMatrix& A,
                                        const std::vector<double>& b,
                                        const std::vector<double>& x0 = {},
                                        int maxIter = 1000,
                                        double tol = 1e-6);

private:
    static double dot(const std::vector<double>& v1,
                      const std::vector<double>& v2);
    static double norm(const std::vector<double>& v);
};
