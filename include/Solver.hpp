#ifndef _Solver_
#define _Solver_

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
private:
    int maxIter;
    double tol;
    bool verbose;
    SolverType solverType;

public:
    Solver(int maxIter = 1000,
           double tol = 1e-6,
           bool verbose = true,
           SolverType solverType = SolverType::GAUSS_SEIDEL)
        : maxIter(maxIter),
          tol(tol),
          verbose(verbose),
          solverType(solverType){};

    std::vector<double> solve(SpaceMatrix& Eqn,
                              const std::vector<double>& x0 = {});

    std::vector<double> jacobi(const SpaceMatrix& Eqn,
                               const std::vector<double>& x0 = {});

    std::vector<double> gaussSeidel(const SpaceMatrix& Eqn,
                                    const std::vector<double>& x0 = {});

    std::vector<double> conjugateGradient(const SpaceMatrix& Eqn,
                                          const std::vector<double>& x0 = {});

    std::vector<double> bicgstab(const SpaceMatrix& Eqn,
                                 const std::vector<double>& x0 = {});

private:
    double dot(const std::vector<double>& v1, const std::vector<double>& v2);
    double norm(const std::vector<double>& v);
};

#endif