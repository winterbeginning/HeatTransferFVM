#include <iostream>
#include "Mesh.hpp"
#include "FiniteVolume.hpp"
#include "Solver.hpp"

int main()
{
    std::cout << "--- Steady State Heat Conduction ---" << std::endl;
    int nx = 50, ny = 50;
    Mesh mesh;
    mesh.createSquareMesh(nx, ny, 1.0, 1.0);
    std::cout << "--- Mesh Cells " << mesh.numCells << " ---" << std::endl;

    FiniteVolume fvm(mesh);

    fvm.setSolveOption(true, false, false);
    fvm.setProperties(1.0, 1.0, 1.0, 0.001);

    Field<double>& T = fvm.T;
    Field<double>& SourceT = fvm.SourceT;

    // 设置初始值
    T.fill(0.0);
    // 设置边界条件：现在与 Mesh 的 Patch 名字强关联
    T.setBoundary("left", 0.0, 0.0, 1.0);
    T.setBoundary("down", 100.0, 0.0, 1.0);
    T.setBoundary("right", 0.0, 0.0, 0.0);
    T.setBoundary("top", 0.0, 0.0, 0.0);

    Field<Vector>& U = fvm.U;

    U.fill(Vector(1.0, 1.0, 0));
    U.setBoundary("left", Vector(1.0, 1.0, 0), Vector(), 1.0);
    U.setBoundary("down", Vector(1.0, 1.0, 0), Vector(), 1.0);
    U.setBoundary("right", Vector(1.0, 1.0, 0), Vector(), 0.0);
    U.setBoundary("top", Vector(1.0, 1.0, 0), Vector(), 0.0);

    Solver solver(1000, 1e-6, true, SolverType::GAUSS_SEIDEL);

    // fvm.solve(TimeScheme::STEADY, solver);

    fvm.solve(TimeScheme::IMPLICIT, solver, 10, 0.2);

    return 0;
}
