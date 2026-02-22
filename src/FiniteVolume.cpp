#include "FiniteVolume.hpp"
#include "FVMDiscrete.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
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
        fvm::Div(Eqn, properties, U, T);
    if (Diffusive)
        fvm::Laplacian(Eqn, properties, T);
    if (Source)
        assembleSource(Eqn);
    if (st == TimeScheme::IMPLICIT)
        fvm::Ddt(Eqn, properties, T, dt);
}

void FiniteVolume::prepareConnectivity()
{
    if (!cachedConnectivity.empty())
        return;

    cachedConnectivity.assign(mesh.numCells, std::vector<int>());
    for (int f = 0; f < (int)mesh.facePoints.size(); ++f)
    {
        int o = mesh.owner[f];
        for (int p : mesh.facePoints[f])
        {
            if (std::find(cachedConnectivity[o].begin(),
                          cachedConnectivity[o].end(),
                          p) == cachedConnectivity[o].end())
                cachedConnectivity[o].push_back(p);
        }

        // 仅处理内部面对应的邻居单元
        if (f < mesh.nInternalFace)
        {
            int n = mesh.neighbour[f];
            for (int p : mesh.facePoints[f])
            {
                if (std::find(cachedConnectivity[n].begin(),
                              cachedConnectivity[n].end(),
                              p) == cachedConnectivity[n].end())
                    cachedConnectivity[n].push_back(p);
            }
        }
    }

    // 排序
    for (int i = 0; i < mesh.numCells; ++i)
    {
        auto& nodes = cachedConnectivity[i];
        if (nodes.size() >= 3)
        {
            Vector center = mesh.cellCentres[i];
            std::sort(nodes.begin(),
                      nodes.end(),
                      [&](int a, int b)
                      {
                          return atan2(mesh.points[a].y - center.y,
                                       mesh.points[a].x - center.x) <
                                 atan2(mesh.points[b].y - center.y,
                                       mesh.points[b].x - center.x);
                      });
        }
    }
}

void FiniteVolume::solve(TimeScheme ts, Solver& solver, int maxSteps, double dt)
{
    prepareConnectivity();

    if (ts == TimeScheme::STEADY)
    {
        writeToTec("out_init.plt");
        TEqn.clear();
        assembleMatrix(TEqn, ts);
        T.internalField = solver.solve(TEqn, T.internalField);
        T.correctBoundaryField();
        T.storeOldField();
        writeToTec("out_steady.plt", 50.0);
    }
    else if (ts == TimeScheme::IMPLICIT)
    {
        for (int step = 0; step <= maxSteps; ++step)
        {
            TEqn.clear();
            assembleMatrix(TEqn, ts, dt);
            T.internalField = solver.solve(TEqn, T.internalField);
            T.correctBoundaryField();
            T.storeOldField();
            if (step % 5 == 0)
            {
                std::cout << "Time step " << step << " done." << std::endl;
                writeToTec("out_" + std::to_string(step) + ".plt",
                           (double)step * dt);
            }
        }
    }
}

void FiniteVolume::writeToTec(const string& filePath, double time)
{
    prepareConnectivity();

    std::ofstream outfile(filePath);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open " << filePath << " for writing!"
                  << std::endl;
        return;
    }

    int nNodes = mesh.points.size();
    int nElems = mesh.numCells;

    // 1. 写入文件头
    outfile << "TITLE = \"Heat Transfer FVM Results\"" << std::endl;
    outfile << "VARIABLES = \"X\", \"Y\", \"T\", \"Vx\", \"Vy\", \"Vz\""
            << std::endl;

    // 2. 写入区域头 (ZONE)
    outfile
        << "ZONE T=\"InternalField\", NODES=" << nNodes
        << ", ELEMENTS=" << nElems << ", DATAPACKING=BLOCK, "
        << "ZONETYPE=FEQUADRILATERAL, VARLOCATION=([3,4,5,6]=CELLCENTERED), "
        << "STRANDID=1, SOLUTIONTIME=" << time << std::endl;

    // 3. 写入节点 X 坐标
    for (int i = 0; i < nNodes; ++i)
        outfile << mesh.points[i].x << ((i + 1) % 10 == 0 ? "\n" : " ");
    outfile << "\n";

    // 4. 写入节点 Y 坐标
    for (int i = 0; i < nNodes; ++i)
        outfile << mesh.points[i].y << ((i + 1) % 10 == 0 ? "\n" : " ");
    outfile << "\n";

    // 5. 写入变量 T (CELLCENTERED)
    writeField(outfile, nElems, T);

    writeField(outfile, nElems, U);

    // 6. 使用缓存的 Connectivity
    for (int i = 0; i < nElems; ++i)
    {
        const auto& nodes = cachedConnectivity[i];
        for (int k = 0; k < 4; ++k)
        {
            int idx = (k < nodes.size()) ? nodes[k] : nodes.back();
            outfile << idx + 1 << (k == 3 ? "" : " ");
        }
        outfile << "\n";
    }

    outfile.close();
}
