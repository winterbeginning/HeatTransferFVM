#ifndef _FiniteVolume_
#define _FiniteVolume_

#include "Mesh.hpp"
#include "Solver.hpp"
#include "Field.hpp"
#include "Point3D.hpp"
#include "Properties.hpp"
#include <fstream>
#include <algorithm>

#include "FVMDiscrete.hpp"
#include <iostream>
#include <cmath>

enum class TimeScheme
{
    STEADY,
    IMPLICIT
};

class FiniteVolume
{
public:
    const Mesh& mesh;
    Properties properties;
    bool Convective;
    bool Diffusive;
    bool Source;
    bool NonOrthogonalCorrection; // 非正交修正开关

    Field<double> T;
    Field<double> SourceT;
    Field<Vector> U;

    // 对矩阵实例进行复用，避免每一时间步或迭代都重新申请内存
    SpaceMatrix<double> TEqn;

    FiniteVolume(const Mesh& mesh)
        : mesh(mesh),
          properties(mesh, 1.0, 1.0, 1.0, 1.0),
          Convective(false),
          Diffusive(true),
          Source(true),
          NonOrthogonalCorrection(true), // 默认开启非正交修正
          T(mesh),
          SourceT(mesh, 0.0),
          U(mesh, Vector(1.0, 1.0, 0)),
          TEqn(mesh.numCells)
    {
    }

    template <typename ValueType>
    void solve(TimeScheme ts,
               Solver<ValueType>& solver,
               int maxSteps = 20,
               double dt = 0.1)
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

    // 通用 Tecplot 输出函数
    void
    writeToTecplot(const string& filePath,
                   double time = 0.0,
                   const std::vector<std::pair<std::string, Field<double>*>>&
                       scalarFields = {},
                   const std::vector<std::pair<std::string, Field<Vector>*>>&
                       vectorFields = {});

    void setSolveOption(bool Convective, bool Diffusive, bool Source)
    {
        this->Convective = Convective;
        this->Diffusive = Diffusive;
        this->Source = Source;
    }

private:
    std::vector<std::vector<int>> cachedConnectivity;
    void prepareConnectivity();

    template <typename ValueType>
    void assembleSource(SpaceMatrix<ValueType>& Eqn)
    {
        for (size_t i = 0; i < mesh.numCells; ++i)
        {
            Eqn.addTob(i, SourceT[i] * mesh.cellVolumes[i]);
        }
    }

    template <typename ValueType>
    void
    assembleMatrix(SpaceMatrix<ValueType>& Eqn, TimeScheme st, double dt = 1.0)
    {
        if (Convective)
        {
            FaceField<double> faceFlux = fvc::flux(U);
            fvm::Div(Eqn, properties, faceFlux, T);
        }
        if (Diffusive)
            fvm::Laplacian(Eqn, properties, T, NonOrthogonalCorrection);
        if (Source)
            assembleSource(Eqn);
        if (st == TimeScheme::IMPLICIT)
            fvm::Ddt(Eqn, properties, T, dt);
    }

    void writeField(ofstream& outfile, int nElements, Field<double> field)
    {
        for (int i = 0; i < nElements; ++i)
            outfile << field[i] << ((i + 1) % 10 == 0 ? "\n" : " ");
        outfile << "\n";
    };
};

inline void FiniteVolume::writeToTecplot(
    const string& filePath,
    double time,
    const std::vector<std::pair<std::string, Field<double>*>>& scalarFields,
    const std::vector<std::pair<std::string, Field<Vector>*>>& vectorFields)
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

    // 1. 检测网格维度
    bool is3D = false;
    for (const auto& pt : mesh.points)
    {
        if (std::abs(pt.z) > 1e-10)
        {
            is3D = true;
            break;
        }
    }

    // 2. 检测单元类型（通过第一个单元的节点数判断）
    int nodesPerElement =
        cachedConnectivity.empty() ? 0 : cachedConnectivity[0].size();
    std::string zoneType;
    bool needsDegenerateConversion = false; // 标记是否需要退化转换

    if (is3D)
    {
        if (nodesPerElement == 8)
            zoneType = "FEBRICK"; // 六面体
        else if (nodesPerElement == 6)
        {
            // 楔形（棱柱）：Tecplot 不直接支持，转换为退化的六面体
            zoneType = "FEBRICK";
            needsDegenerateConversion = true;
            nodesPerElement = 8; // 输出时需要 8 个节点
        }
        else if (nodesPerElement == 5)
        {
            // 金字塔：也转换为退化的六面体
            zoneType = "FEBRICK";
            needsDegenerateConversion = true;
            nodesPerElement = 8;
        }
        else if (nodesPerElement == 4)
            zoneType = "FETETRAHEDRON"; // 四面体
        else
        {
            // 其他多面体：尝试使用 FEPOLYHEDRON（需要 Tecplot 360）
            zoneType = "FEBRICK";
            needsDegenerateConversion = true;
            nodesPerElement = 8;
            std::cerr << "Warning: Cell has " << cachedConnectivity[0].size()
                      << " nodes. Converting to degenerate hexahedron."
                      << std::endl;
        }
    }
    else
    {
        if (nodesPerElement == 4)
            zoneType = "FEQUADRILATERAL"; // 四边形
        else if (nodesPerElement == 3)
            zoneType = "FETRIANGLE"; // 三角形
        else
            zoneType = "FEPOLYGON"; // 多边形
    }

    // 3. 构建变量列表
    std::vector<std::string> varNames;
    varNames.push_back("X");
    varNames.push_back("Y");
    if (is3D)
        varNames.push_back("Z");

    // 添加标量场名称
    for (const auto& sf : scalarFields)
        varNames.push_back(sf.first);

    // 添加矢量场分量名称
    for (const auto& vf : vectorFields)
    {
        varNames.push_back(vf.first + "_X");
        varNames.push_back(vf.first + "_Y");
        if (is3D)
            varNames.push_back(vf.first + "_Z");
    }

    // 4. 写入文件头
    outfile << "TITLE = \"FVM Results\"" << std::endl;
    outfile << "VARIABLES = ";
    for (size_t i = 0; i < varNames.size(); ++i)
    {
        outfile << "\"" << varNames[i] << "\"";
        if (i < varNames.size() - 1)
            outfile << ", ";
    }
    outfile << std::endl;

    // 5. 构建 VARLOCATION 字符串（场变量在单元中心）
    std::string varLocation;
    int coordCount = is3D ? 3 : 2;
    int totalVars = varNames.size();
    if (totalVars > coordCount)
    {
        varLocation = "VARLOCATION=([" + std::to_string(coordCount + 1) + "-" +
                      std::to_string(totalVars) + "]=CELLCENTERED)";
    }

    // 6. 写入区域头
    outfile << "ZONE T=\"InternalField\", " << "NODES=" << nNodes << ", "
            << "ELEMENTS=" << nElems << ", " << "DATAPACKING=BLOCK, "
            << "ZONETYPE=" << zoneType;
    if (!varLocation.empty())
        outfile << ", " << varLocation;
    outfile << ", STRANDID=1, SOLUTIONTIME=" << time << std::endl;

    // 7. 写入节点坐标（BLOCK 格式）
    // X 坐标
    for (int i = 0; i < nNodes; ++i)
        outfile << mesh.points[i].x << ((i + 1) % 10 == 0 ? "\n" : " ");
    outfile << "\n";

    // Y 坐标
    for (int i = 0; i < nNodes; ++i)
        outfile << mesh.points[i].y << ((i + 1) % 10 == 0 ? "\n" : " ");
    outfile << "\n";

    // Z 坐标（3D）
    if (is3D)
    {
        for (int i = 0; i < nNodes; ++i)
            outfile << mesh.points[i].z << ((i + 1) % 10 == 0 ? "\n" : " ");
        outfile << "\n";
    }

    // 8. 写入标量场（CELLCENTERED）
    for (const auto& sf : scalarFields)
    {
        writeField(outfile, nElems, *sf.second);
    }

    // 9. 写入矢量场（CELLCENTERED）
    for (const auto& vf : vectorFields)
    {
        // X 分量
        for (int i = 0; i < nElems; ++i)
            outfile << (*vf.second)[i].x << ((i + 1) % 10 == 0 ? "\n" : " ");
        outfile << "\n";

        // Y 分量
        for (int i = 0; i < nElems; ++i)
            outfile << (*vf.second)[i].y << ((i + 1) % 10 == 0 ? "\n" : " ");
        outfile << "\n";

        // Z 分量（3D）
        if (is3D)
        {
            for (int i = 0; i < nElems; ++i)
                outfile << (*vf.second)[i].z
                        << ((i + 1) % 10 == 0 ? "\n" : " ");
            outfile << "\n";
        }
    }

    // 10. 写入连接性（Connectivity）
    for (int i = 0; i < nElems; ++i)
    {
        const auto& nodes = cachedConnectivity[i];

        if (needsDegenerateConversion && nodes.size() == 6)
        {
            // 楔形（6节点）-> 退化六面体（8节点）
            // 底面三角形: 0,1,2 -> 四边形: 0,1,2,2
            // 顶面三角形: 3,4,5 -> 四边形: 3,4,5,5
            outfile << nodes[0] + 1 << " " << nodes[1] + 1 << " "
                    << nodes[2] + 1 << " " << nodes[2] + 1 << " "
                    << nodes[3] + 1 << " " << nodes[4] + 1 << " "
                    << nodes[5] + 1 << " " << nodes[5] + 1 << "\n";
        }
        else if (needsDegenerateConversion && nodes.size() == 5)
        {
            // 金字塔（5节点）-> 退化六面体（8节点）
            // 底面四边形: 0,1,2,3
            // 顶点: 4 (重复4次)
            outfile << nodes[0] + 1 << " " << nodes[1] + 1 << " "
                    << nodes[2] + 1 << " " << nodes[3] + 1 << " "
                    << nodes[4] + 1 << " " << nodes[4] + 1 << " "
                    << nodes[4] + 1 << " " << nodes[4] + 1 << "\n";
        }
        else
        {
            // 标准单元类型（直接输出）
            int requiredNodes = nodesPerElement;
            for (int k = 0; k < requiredNodes; ++k)
            {
                int idx = (k < (int)nodes.size()) ? nodes[k] : nodes.back();
                outfile << idx + 1; // Tecplot 使用 1-based 索引
                if (k < requiredNodes - 1)
                    outfile << " ";
            }
            outfile << "\n";
        }
    }

    outfile.close();
    std::cout << "Tecplot file written: " << filePath
              << " (ZoneType: " << zoneType;
    if (needsDegenerateConversion)
        std::cout << ", Degenerate conversion applied";
    std::cout << ")" << std::endl;
}

inline void FiniteVolume::prepareConnectivity()
{
    if (!cachedConnectivity.empty())
        return;

    cachedConnectivity.assign(mesh.numCells, std::vector<int>());

    // 为每个单元收集所有唯一的节点
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

    // 检测是否为 3D 网格
    bool is3D = false;
    for (const auto& pt : mesh.points)
    {
        if (std::abs(pt.z) > 1e-10)
        {
            is3D = true;
            break;
        }
    }

    // 对节点进行排序
    for (int i = 0; i < mesh.numCells; ++i)
    {
        auto& nodes = cachedConnectivity[i];

        if (is3D)
        {
            // 3D：按 Z 坐标分层，然后在每层内按极角排序
            if (nodes.size() >= 4)
            {
                Vector center = mesh.cellCentres[i];

                // 按 Z 坐标排序（分为底面和顶面）
                std::sort(nodes.begin(),
                          nodes.end(),
                          [&](int a, int b)
                          { return mesh.points[a].z < mesh.points[b].z; });

                // 六面体（8节点）：分别对底面和顶面排序
                if (nodes.size() == 8)
                {
                    // 底面（前4个节点）按极角排序
                    std::sort(nodes.begin(),
                              nodes.begin() + 4,
                              [&](int a, int b)
                              {
                                  return atan2(mesh.points[a].y - center.y,
                                               mesh.points[a].x - center.x) <
                                         atan2(mesh.points[b].y - center.y,
                                               mesh.points[b].x - center.x);
                              });

                    // 顶面（后4个节点）按极角排序
                    std::sort(nodes.begin() + 4,
                              nodes.end(),
                              [&](int a, int b)
                              {
                                  return atan2(mesh.points[a].y - center.y,
                                               mesh.points[a].x - center.x) <
                                         atan2(mesh.points[b].y - center.y,
                                               mesh.points[b].x - center.x);
                              });
                }
                // 楔形（6节点）：底面3个节点，顶面3个节点
                else if (nodes.size() == 6)
                {
                    // 底面（前3个节点）按极角排序
                    std::sort(nodes.begin(),
                              nodes.begin() + 3,
                              [&](int a, int b)
                              {
                                  return atan2(mesh.points[a].y - center.y,
                                               mesh.points[a].x - center.x) <
                                         atan2(mesh.points[b].y - center.y,
                                               mesh.points[b].x - center.x);
                              });

                    // 顶面（后3个节点）按极角排序
                    std::sort(nodes.begin() + 3,
                              nodes.end(),
                              [&](int a, int b)
                              {
                                  return atan2(mesh.points[a].y - center.y,
                                               mesh.points[a].x - center.x) <
                                         atan2(mesh.points[b].y - center.y,
                                               mesh.points[b].x - center.x);
                              });
                }
                // 金字塔（5节点）：底面4个节点 + 顶点1个
                else if (nodes.size() == 5)
                {
                    // 底面（前4个节点）按极角排序
                    std::sort(nodes.begin(),
                              nodes.begin() + 4,
                              [&](int a, int b)
                              {
                                  return atan2(mesh.points[a].y - center.y,
                                               mesh.points[a].x - center.x) <
                                         atan2(mesh.points[b].y - center.y,
                                               mesh.points[b].x - center.x);
                              });
                    // 顶点（第5个节点）不需要排序
                }
            }
        }
        else
        {
            // 2D：按极角排序（原逻辑）
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
}

#endif