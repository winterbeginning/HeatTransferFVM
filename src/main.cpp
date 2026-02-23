#include <iostream>
#include "Mesh.hpp"
#include "FiniteVolume.hpp"
#include "Solver.hpp"

int main()
{
    std::cout << "--- Steady State Heat Conduction ---" << std::endl;
    int nx = 50, ny = 50;
    Mesh mesh;

    // 选择网格类型（取消注释其中一个）
    // mesh.ReadFromFile("/home/winter/HeatTransferFVM/polyMesh");
    mesh.createSquareMesh(nx, ny, 1.0, 1.0);

    std::cout << "--- Mesh Cells " << mesh.numCells << " ---" << std::endl;

    FiniteVolume fvm(mesh);

    // 求解选项：对流、扩散、源项
    fvm.setSolveOption(true, true, false); // 开启扩散求解

    // 非正交修正开关（仅对非结构网格有影响）
    fvm.NonOrthogonalCorrection = true; // true=开启, false=关闭

    fvm.properties.setProperties(0.001, 1.0, 1.0, 0.001);

    Field<double>& T = fvm.T;

    // 设置初始值
    T.fill(0.0);

    // 设置边界条件
    T.setBoundary("left", 0.0, 0.0, 1.0);  // Dirichlet: 0°C
    T.setBoundary("down", 1.0, 0.0, 1.0);  // Dirichlet: 100°C
    T.setBoundary("right", 0.0, 0.0, 0.0); // Neumann: 绝热
    T.setBoundary("top", 0.0, 0.0, 0.0);   // Neumann: 绝热

    Field<Vector>& U = fvm.U;

    Vector velocity(1.0, 1.0, 0);
    Vector zero(0.0, 0.0, 0.0);

    U.fill(velocity);
    U.setBoundary("left", velocity, zero, 1.0);
    U.setBoundary("down", velocity, zero, 1.0);
    U.setBoundary("right", velocity, zero, 0.0);
    U.setBoundary("top", velocity, zero, 0.0);

    Solver solver(5000, 1e-6, true, SolverType::GAUSS_SEIDEL);

    fvm.solve(TimeScheme::STEADY, solver);

    // fvm.solve(TimeScheme::IMPLICIT, solver, 10, 0.2);

    return 0;
}
