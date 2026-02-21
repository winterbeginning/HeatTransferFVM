#include "FiniteVolume.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

void FiniteVolume::assembleSource(SpaceMatrix& A, std::vector<double>& b)
{
    for (size_t i = 0; i < mesh.numCells; i++)
    {
        b[i] += SourceT[i] * mesh.cellVolumes[i];
    }
}

void FiniteVolume::assembleSteady(SpaceMatrix& A, std::vector<double>& b)
{
    // 1. 处理内部面
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];
        double area = mesh.faceAreas[i];
        double dist = (mesh.cellCentres[n] - mesh.cellCentres[o]).getMag();
        double coeff = k * area / dist;

        A.addValue(o, o, coeff);
        A.addValue(o, n, -coeff);
        A.addValue(n, n, coeff);
        A.addValue(n, o, -coeff);
    }

    // 2. 处理边界
    for (int pIdx = 0; pIdx < mesh.boundary.size(); ++pIdx)
    {
        const auto& patch = mesh.boundary[pIdx];
        const auto& bc = T.boundaryField[pIdx];
        for (int i = patch.firstFaceIdx; i < patch.firstFaceIdx + patch.nFaces;
             ++i)
        {
            int o = mesh.owner[i];
            double area = mesh.faceAreas[i];
            double d = (mesh.faceCentres[i] - mesh.cellCentres[o]).getMag();

            double coeff_p = k * area * bc.fraction / d;
            double source = k * area * bc.fraction * bc.refValue / d +
                            k * area * (1.0 - bc.fraction) * bc.refGrad;

            A.addValue(o, o, coeff_p);
            b[o] += source;
        }
    }
}

void FiniteVolume::assembleImplicit(SpaceMatrix& A,
                                    std::vector<double>& b,
                                    double dt)
{
    assembleSteady(A, b);
    for (int i = 0; i < mesh.numCells; ++i)
    {
        double trans = mesh.cellVolumes[i] * rhoCp / dt;
        A.addValue(i, i, trans);
        b[i] += trans * T_old[i];
    }
}

void FiniteVolume::stepExplicit(double dt)
{
    std::vector<double> T_new = T.internalField;
    std::vector<double> cellFlux(mesh.numCells, 0.0);

    // 1. 内部面通量
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];
        double area = mesh.faceAreas[i];
        double dist = (mesh.cellCentres[n] - mesh.cellCentres[o]).getMag();

        double flux = k * area * (T[n] - T[o]) / dist;
        cellFlux[o] += flux;
        cellFlux[n] -= flux;
    }

    // 2. 边界面通量
    for (int pIdx = 0; pIdx < mesh.boundary.size(); ++pIdx)
    {
        const auto& patch = mesh.boundary[pIdx];
        const auto& bc = T.boundaryField[pIdx];
        for (int i = patch.firstFaceIdx; i < patch.firstFaceIdx + patch.nFaces;
             ++i)
        {
            int o = mesh.owner[i];
            double area = mesh.faceAreas[i];
            double d = (mesh.faceCentres[i] - mesh.cellCentres[o]).getMag();

            double Tf = bc.fraction * bc.refValue +
                        (1.0 - bc.fraction) * (T[o] + bc.refGrad * d);
            double flux = k * area * (Tf - T[o]) / d;
            cellFlux[o] += flux;
        }
    }

    for (int i = 0; i < mesh.numCells; ++i)
    {
        if (mesh.cellVolumes[i] > 1e-20)
            T_new[i] =
                T[i] + (dt / (rhoCp * mesh.cellVolumes[i])) * cellFlux[i];
    }
    T.internalField = T_new;
    T_old.internalField = T_new;
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

void FiniteVolume::solve(TimeScheme ts, SolverType st, double dt, int maxSteps)
{
    prepareConnectivity();

    if (ts == TimeScheme::STEADY)
    {
        SpaceMatrix A(mesh.numCells);
        std::vector<double> b(mesh.numCells, 0.0);
        assembleSteady(A, b);
        assembleSource(A, b);
        T.internalField = Solver::solve(st, A, b, T.internalField);
        writeToTec("out_steady.plt", 50.0);
    }
    else if (ts == TimeScheme::IMPLICIT)
    {
        for (int step = 0; step <= maxSteps; ++step)
        {
            SpaceMatrix A(mesh.numCells);
            std::vector<double> b(mesh.numCells, 0.0);
            assembleImplicit(A, b, dt);
            assembleSource(A, b);
            T.internalField = Solver::solve(st, A, b, T.internalField);
            T_old.internalField = T.internalField;
            if (step % 5 == 0)
            {
                std::cout << "Time step " << step << " done." << std::endl;
                writeToTec("out_" + std::to_string(step) + ".plt",
                           (double)step * dt);
            }
        }
    }
    else if (ts == TimeScheme::EXPLICIT)
    {
        for (int step = 0; step <= maxSteps; ++step)
        {
            stepExplicit(dt);
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
    outfile << "VARIABLES = \"X\", \"Y\", \"T\"" << std::endl;

    // 2. 写入区域头 (ZONE)
    outfile << "ZONE T=\"InternalField\", NODES=" << nNodes
            << ", ELEMENTS=" << nElems << ", DATAPACKING=BLOCK, "
            << "ZONETYPE=FEQUADRILATERAL, VARLOCATION=([3]=CELLCENTERED), "
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
    for (int i = 0; i < nElems; ++i)
        outfile << T[i] << ((i + 1) % 10 == 0 ? "\n" : " ");
    outfile << "\n";

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
