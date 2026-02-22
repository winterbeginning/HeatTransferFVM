#include "FiniteVolume.hpp"
#include "FVMDiscrete.hpp"
#include <iostream>
#include <cmath>

void FiniteVolume::assembleSource(SpaceMatrix& Eqn)
{
    for (size_t i = 0; i < mesh.numCells; ++i)
    {
        Eqn.addTob(i, SourceT[i] * mesh.cellVolumes[i]);
    }
}

void FiniteVolume::assembleMatrix(SpaceMatrix& Eqn, TimeScheme st, double dt)
{
    if (Convective)
    {
        SurfaceField<double> faceFlux = fvc::flux(U);
        fvm::Div(Eqn, properties, faceFlux, T);
    }
    if (Diffusive)
        fvm::Laplacian(Eqn, properties, T);
    if (Source)
        assembleSource(Eqn);
    if (st == TimeScheme::IMPLICIT)
        fvm::Ddt(Eqn, properties, T, dt);
}

void FiniteVolume::solve(TimeScheme ts, Solver& solver, int maxSteps, double dt)
{
    prepareConnectivity();

    // 初始修正边界
    T.correctBoundaryField();
    U.correctBoundaryField();

    if (ts == TimeScheme::STEADY)
    {
        writeToTecplot("out_init.plt", 0.0, {{"T", &T}}, {{"U", &U}});

        // 外层迭代：让非正交修正逐步收敛
        int nOuterIter = 3;
        for (int outerIter = 0; outerIter < nOuterIter; ++outerIter)
        {
            TEqn.clear();
            assembleMatrix(TEqn, ts);
            T.internalField = solver.solve(TEqn, T.internalField);
            T.correctBoundaryField();
            T.storeOldField();
        }

        writeToTecplot("out_steady.plt", 50.0, {{"T", &T}}, {{"U", &U}});
    }
    else if (ts == TimeScheme::IMPLICIT)
    {
        for (int step = 0; step <= maxSteps; ++step)
        {
            T.correctBoundaryField();
            TEqn.clear();
            assembleMatrix(TEqn, ts, dt);
            T.internalField = solver.solve(TEqn, T.internalField);
            T.correctBoundaryField();
            T.storeOldField();
            if (step % 5 == 0)
            {
                std::cout << "Time step " << step << " done." << std::endl;
                writeToTecplot("out_" + std::to_string(step) + ".plt",
                               (double)step * dt,
                               {{"T", &T}},
                               {{"U", &U}});
            }
        }
    }
}
