#include <iostream>
#include "Mesh.hpp"
#include "PostProcessor.hpp"
#include "Solver.hpp"
#include "Properties.hpp"
#include "fvMatrix.hpp"
#include "FiniteVolume.hpp"

int main()
{
    std::cout << "--- Transient Heat Conduction (OpenFOAM Style) ---"
              << std::endl;

    // 1. 网格初始化
    Mesh mesh;
    // 使用非结构网格/OpenFOAM网格
    // mesh.ReadFromFile("/home/winter/HeatTransferFVM/polyMesh");
    mesh.createSquareMesh(50, 50, 1.0, 1.0);
    std::cout << "--- Mesh Cells " << mesh.numCells << " ---" << std::endl;

    // 2. 场与物性初始化
    Field<double> T(mesh, 0.0); // 初始温度 20°C
    Field<Vector> U(mesh, Vector(1.0, 1.0, 0.0));
    Properties properties(mesh, 0.01, 1.0, 1.0, 0.001); // kappa, rho, Cp, mu

    // 设置边界条件
    T.setBoundary("left", 0.0, 0.0, 1.0);   // Dirichlet: 0°C
    T.setBoundary("down", 100.0, 0.0, 1.0); // Dirichlet: 100°C
    T.setBoundary("right", 0.0, 0.0, 0.0);  // Neumann: 绝热
    T.setBoundary("top", 0.0, 0.0, 0.0);    // Neumann: 绝热
    T.correctBoundaryField();
    T.storeOldField(); // 存储初始场作为 Time 0

    Vector Veocity(1.0, 1.0, 0.0);
    Vector Zero(0.0, 0.0, 0.0);

    U.setBoundary("left", Veocity, Zero, 1.0);
    U.setBoundary("down", Veocity, Zero, 1.0);
    U.setBoundary("right", Veocity, Zero, 0.0);
    U.setBoundary("top", Veocity, Zero, 0.0);
    U.correctBoundaryField();

    // 3. 求解器与工具
    PostProcessor post(mesh);
    Solver<double> solver(
        1000,
        1e-8,
        false,
        SolverType::GAUSS_SEIDEL); // 瞬态通常用更快的线性求解器

    fvMatrix<double> TEqn(T);
    TEqn.setSolver(solver);

    // 4. 时间步参数
    double dt = 0.01;
    double endTime = 2.0;
    int maxSteps = static_cast<int>(endTime / dt);
    bool enableNonOrthogonalCorrection = true;

    post.writeToTecplot("out_transient_0.plt", 0.0, {{"T", &T}}, {{"U", &U}});

    // 5. 时间步循环
    for (int step = 1; step <= maxSteps; ++step)
    {
        double currentTime = step * dt;

        // 每个时间步清理并更新矩阵
        TEqn.clear();

        // 计算通量
        FaceField<double> faceFlux = fvc::flux(U);

        // --- 组装瞬态方程式 (fvm 风格) ---

        // 时间项: d(rho*Cp*T)/dt
        fvm::Ddt(TEqn, properties, T, dt);

        // 对流项
        fvm::Div(TEqn, properties, faceFlux, T);

        // 扩散项
        // fvm::Laplacian(
        //     TEqn, properties.kappa(), T, enableNonOrthogonalCorrection);

        // 求解
        TEqn.solve();

        // 边界修正与旧场更新 (非常重要)
        T.correctBoundaryField();
        T.storeOldField();

        if (step % 20 == 0)
        {
            std::cout << "Time = " << currentTime << "s (Step " << step
                      << ") done." << std::endl;
            post.writeToTecplot("out_transient_" + std::to_string(step) +
                                    ".plt",
                                currentTime,
                                {{"T", &T}},
                                {{"U", &U}});
        }
    }

    std::cout << "Transient simulation completed." << std::endl;

    return 0;
}
