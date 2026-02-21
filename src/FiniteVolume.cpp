#include "FiniteVolume.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

void FiniteVolume::assembleDiv(SpaceMatrix& A, std::vector<double>& b)
{
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];
        double area = mesh.faceAreas[i];
        Vector areaDirect = mesh.faceNormals[i];
        Vector U_face = (U[o] + U[n]) / 2.0;

        double flux = rho * U_face.dotWith(areaDirect) * area;

        if (flux > 0) // 流动方向 o --> n   面变量等于o单元
        {
            A.addValue(o, o, flux); // o单元方程，面变量来自o，流出为正
            A.addValue(n, o, -flux); // n单元方程，面变量来自o，流入为负
        }
        else // 流动方向 n --> o   面变量等于n单元
        {
            A.addValue(o, n, flux); // o单元方程，面变量来自n，流入为正
            A.addValue(n, n, -flux); // n单元方程，面变量来自n，流出为负
        };
    }

    for (int pIdx = 0; pIdx < mesh.boundary.size(); ++pIdx)
    {
        const auto& patch = mesh.boundary[pIdx];
        const auto& Ubc = U.boundaryField[pIdx];
        const auto& phi_bc = T.boundaryField[pIdx];
        for (int i = patch.firstFaceIdx; i < patch.firstFaceIdx + patch.nFaces;
             ++i)
        {
            int o = mesh.owner[i];
            double area = mesh.faceAreas[i];
            Vector areaDirect = mesh.faceNormals[i];
            double d = (mesh.faceCentres[i] - mesh.cellCentres[o]).getMag();
            Vector U_face = Ubc.fraction * Ubc.refValue +
                            (1.0 - Ubc.fraction) * (U[o] + Ubc.refGrad * d);

            double phi_face =
                phi_bc.fraction * phi_bc.refValue +
                (1.0 - phi_bc.fraction) * (T[o] + phi_bc.refGrad * d);

            double flux = rho * U_face.dotWith(areaDirect) * area;

            if (flux > 0) // 流出
            {
                A.addValue(o, o, flux);
            }
            else // 流入
            {
                b[o] -= flux * phi_face;
            };
        }
    }
}

void FiniteVolume::assembleLaplacian(SpaceMatrix& A, std::vector<double>& b)
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

void FiniteVolume::assembleSource(SpaceMatrix& A, std::vector<double>& b)
{
    for (size_t i = 0; i < mesh.numCells; i++)
    {
        b[i] += SourceT[i] * mesh.cellVolumes[i];
    }
}

void FiniteVolume::assembleSpace(SpaceMatrix& A, std::vector<double>& b)
{
    if (Convective)
        assembleDiv(A, b);
    if (Diffusive)
        assembleLaplacian(A, b);
    if (Source)
        assembleSource(A, b);
}

void FiniteVolume::assembleTime(SpaceMatrix& A,
                                std::vector<double>& b,
                                double dt)
{
    assembleSpace(A, b);
    for (int i = 0; i < mesh.numCells; ++i)
    {
        double trans = mesh.cellVolumes[i] * rho * Cp / dt;
        A.addValue(i, i, trans);
        b[i] += trans * T_old[i];
    }
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
        SpaceMatrix A(mesh.numCells);
        std::vector<double> b(mesh.numCells, 0.0);
        assembleSpace(A, b);
        T.internalField = solver.solve(A, b, T.internalField);
        writeToTec("out_steady.plt", 50.0);
    }
    else if (ts == TimeScheme::IMPLICIT)
    {
        for (int step = 0; step <= maxSteps; ++step)
        {
            SpaceMatrix A(mesh.numCells);
            std::vector<double> b(mesh.numCells, 0.0);
            assembleTime(A, b, dt);
            T.internalField = solver.solve(A, b, T.internalField);
            T_old.internalField = T.internalField;
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
