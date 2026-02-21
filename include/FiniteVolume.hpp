#pragma once

#include "Mesh.hpp"
#include "Solver.hpp"
#include "Field.hpp"

enum class TimeScheme
{
    STEADY,
    EXPLICIT,
    IMPLICIT
};

class FiniteVolume
{
public:
    const Mesh& mesh;
    double k;     // Thermal conductivity
    double rhoCp; // Density * heat capacity

    Field<double> T;
    Field<double> T_old;

    FiniteVolume(const Mesh& mesh, double k, double rhoCp)
        : mesh(mesh), k(k), rhoCp(rhoCp), T("T", mesh), T_old("Told", mesh)
    {
    }

    void solve(TimeScheme ts,
               SolverType st = SolverType::CG,
               double dt = 0.0,
               int maxSteps = 1);

    void writeToTec(const string& filePath = "out.dat", double time = 0.0);

private:
    std::vector<std::vector<int>> cachedConnectivity;
    void prepareConnectivity();

    void assembleSteady(SpaceMatrix& A, std::vector<double>& b);
    void assembleImplicit(SpaceMatrix& A, std::vector<double>& b, double dt);
    void stepExplicit(double dt);
};
