#include <iostream>
#include "Mesh.hpp"
#include "PostProcessor.hpp"
#include "Solver.hpp"
#include "Properties.hpp"
#include "fvMatrix.hpp"
#include "FiniteVolume.hpp"

int main()
{
    std::cout << "--- Steady State Heat Conduction (OpenFOAM Style) ---"
              << std::endl;

    // 1. 网格初始化
    Mesh mesh;
    mesh.ReadFromFile("/home/winter/HeatTransferFVM/polyMesh");
    // mesh.createSquareMesh(50, 50, 1.0, 1.0);
    std::cout << "--- Mesh Cells " << mesh.numCells << " ---" << std::endl;

    // 2. 场与物性初始化
    Field<double> T(mesh, 0.0);
    Field<Vector> U(mesh, Vector(1.0, 1.0, 0.0));
    Properties properties(mesh, 0.001, 1.0, 1.0, 0.001); // kappa, rho, Cp, mu

    // 设置边界条件
    T.setBoundary("left", 0.0, 0.0, 1.0);  // Dirichlet: 0°C
    T.setBoundary("down", 1.0, 0.0, 1.0);  // Dirichlet: 100°C
    T.setBoundary("right", 0.0, 0.0, 0.0); // Neumann: 绝热
    T.setBoundary("top", 0.0, 0.0, 0.0);   // Neumann: 绝热
    T.correctBoundaryField();

    Vector Veocity(1.0, 1.0, 0.0);
    Vector Zero(0.0, 0.0, 0.0);

    U.setBoundary("left", Veocity, Zero, 1.0);
    U.setBoundary("down", Veocity, Zero, 1.0);
    U.setBoundary("right", Veocity, Zero, 0.0);
    U.setBoundary("top", Veocity, Zero, 0.0);
    U.correctBoundaryField();

    // 3. 求解器与工具类
    PostProcessor post(mesh);
    Solver<double> solver(5000, 1e-6, true, SolverType::GAUSS_SEIDEL);

    // 将矩阵、求解器、场绑定
    fvMatrix<double> TEqn(T);
    TEqn.setSolver(solver);

    // 4. 稳态求解流程
    int nOuterIter = 30;
    bool enableNonOrthogonalCorrection = true;

    post.writeToTecplot("out_init.plt", 0.0, {{"T", &T}}, {{"U", &U}});

    for (int outerIter = 0; outerIter < nOuterIter; ++outerIter)
    {
        TEqn.clear();

        // --- 组装方程式 (fvm 风格) ---

        // 计算通量
        FaceField<double> faceFlux = fvc::flux(U);

        // 对流项
        fvm::Div(TEqn, properties, faceFlux, T);

        // 扩散项
        fvm::Laplacian(
            TEqn, properties.kappa(), T, enableNonOrthogonalCorrection);

        // 求解 (直接在矩阵对象上调用 solve)
        TEqn.solve();

        T.correctBoundaryField();
        T.storeOldField();

        if (outerIter % 5 == 0)
        {
            std::cout << "Iteration " << outerIter << " done." << std::endl;
        }
    }

    post.writeToTecplot("out_steady.plt", 1.0, {{"T", &T}}, {{"U", &U}});

    return 0;
}
