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

    FiniteVolume fvm(mesh, 1.0, 1.0);

    Field<double>& T = fvm.T;
    Field<double>& SourceT = fvm.SourceT;

    // 设置边界条件：现在与 Mesh 的 Patch 名字强关联
    T.setBoundary("left", 100.0, 0.0, 1.0);
    T.setBoundary("right", 50.0, 0.0, 1.0);
    T.setBoundary("top", 20.0, 0.0, 1.0);
    T.setBoundary("down", 0.0, 0.0, 0.0);

    // 设置初始值
    T.fill(20.0);

    for (size_t i = 0; i < mesh.numCells; i++)
    {
        if (mesh.cellCentres[i].x > 0.4 && mesh.cellCentres[i].x < 0.6 &&
            mesh.cellCentres[i].y > 0.4 && mesh.cellCentres[i].y < 0.6)
        {
            SourceT[i] = 1500.0;
        }
    }

    fvm.writeToTec("out_init.plt");

    fvm.solve(TimeScheme::STEADY, SolverType::BICGSTAB);

    // fvm.solve(TimeScheme::IMPLICIT, SolverType::BICGSTAB, 0.02, 20);

    return 0;
}
