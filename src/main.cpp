#include <iostream>
#include "Mesh.hpp"
#include "FiniteVolume.hpp"
#include "Solver.hpp"

int main()
{
    int nx = 10, ny = 10;
    Mesh mesh;
    mesh.createSquareMesh(nx, ny, 1.0, 1.0);

    FiniteVolume fvm(mesh, 1.0, 1.0);

    Field<double>& T = fvm.T;

    // 设置边界条件：现在与 Mesh 的 Patch 名字强关联
    T.setBoundary("left", 100.0, 0.0, 1.0);
    T.setBoundary("right", 50.0, 0.0, 1.0);
    T.setBoundary("top", 20.0, 0.0, 1.0);
    T.setBoundary("down", 0.0, 0.0, 0.0);

    // 设置初始值
    T.fill(20.0);

    std::cout << "--- Steady State Heat Conduction (CG) ---" << std::endl;
    // fvm.solve(TimeScheme::STEADY, SolverType::CG);

    fvm.solve(TimeScheme::IMPLICIT, SolverType::BICGSTAB, 0.02, 20);

    return 0;
}
