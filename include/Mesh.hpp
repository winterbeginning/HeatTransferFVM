#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include "Point3D.hpp"

// 定义边界组（Boundary Patches）
struct BoundaryPatch
{
    std::string name;
    std::vector<int> faceIndices;
    int firstFaceIdx;
    int nFaces;
};

class Mesh
{
public:
    // --- 拓扑描述 (Topological description) ---
    std::vector<Vector> points;               // 所有顶点坐标
    std::vector<std::vector<int>> facePoints; // 每个面包含的顶点ID
    std::vector<int> owner;                   // 每个面所属的主单元ID
    std::vector<int> neighbour; // 每个面的相邻单元ID (边界为-1)
    int nInternalFace;

    // --- 几何属性 (Pre-calculated Geometry) ---
    std::vector<Vector> cellCentres; // 单元中心
    std::vector<double> cellVolumes; // 单元体积
    std::vector<Vector> faceCentres; // 面中心
    std::vector<double> faceAreas;   // 面面积
    std::vector<Vector> faceNormals; // 面法向量 (指向外侧/neighbour)

    std::vector<BoundaryPatch> boundary;
    int numCells;

    Mesh() : numCells(0)
    {
    }

    void ReadFromFile(const string& meshDic)
    {
        readOwnerFile(meshDic + "/owner");
        readNeighbourFile(meshDic + "/neighbour");
        readPointsFile(meshDic + "/points");
        readFacesFile(meshDic + "/faces");
        readBoundaryFile(meshDic + "/boundary");
        nInternalFace = neighbour.size();
        calculateGeometry();
    };

    void readOwnerFile(const string& fileName)
    {
        ifstream file;
        string line;
        string old_line;
        file.open(fileName, ios_base::in);
        bool inDataBlock = false;

        if (!file.is_open())
        {
            std::cerr << "Error: 无法打开owner文件: " << fileName << std::endl;
            return;
        }

        while (std::getline(file, line))
        {
            if (line.find("(") != std::string::npos && !inDataBlock)
            {
                inDataBlock = true;
                owner.reserve(std::stoi(old_line));
                continue;
            }
            if (line.find(")") != std::string::npos && inDataBlock)
            {
                break;
            }
            if (inDataBlock)
            {
                if (line.empty())
                    continue;
                owner.push_back(std::stoi(line));
            }
            else
            {
                old_line = line;
            }
        }
        file.close();
    };

    void readNeighbourFile(const string& fileName)
    {
        ifstream file;
        string line;
        string old_line;
        file.open(fileName, ios_base::in);
        bool inDataBlock = false;

        if (!file.is_open())
        {
            std::cerr << "Error: 无法打开neighbour文件: " << fileName
                      << std::endl;
            return;
        }

        while (std::getline(file, line))
        {
            if (line.find("(") != std::string::npos && !inDataBlock)
            {
                inDataBlock = true;
                neighbour.reserve(std::stoi(old_line));
                continue;
            }
            if (line.find(")") != std::string::npos && inDataBlock)
            {
                break;
            }
            if (inDataBlock)
            {
                if (line.empty())
                    continue;
                neighbour.push_back(std::stoi(line));
            }
            else
            {
                old_line = line;
            }
        }
        file.close();
    };

    void readPointsFile(const string& fileName)
    {
        ifstream file;
        string line;
        string old_line;
        file.open(fileName, ios_base::in);
        bool inDataBlock = false;

        if (!file.is_open())
        {
            std::cerr << "Error: 无法打开points文件: " << fileName << std::endl;
            return;
        }

        while (std::getline(file, line))
        {
            // 1. 找到外层左括号，初始化points数组大小
            if (line.find("(") != std::string::npos && !inDataBlock)
            {
                inDataBlock = true;
                points.reserve(std::stoi(old_line));
                continue;
            }

            // 2. 找到外层右括号，退出读取
            if (line.find(")") != std::string::npos && inDataBlock &&
                line.find("(") == std::string::npos)
            {
                break;
            }

            if (inDataBlock)
            {
                if (line.empty())
                    continue;

                size_t left = line.find("(");
                size_t right = line.find(")");
                if (left == std::string::npos || right == std::string::npos)
                {
                    continue;
                }
                std::string coordStr = line.substr(left + 1, right - left - 1);

                // 解析x/y/z三个浮点数值
                std::istringstream iss(coordStr);
                double x, y, z;
                if (iss >> x >> y >> z)
                {
                    points.push_back(Vector(x, y, z));
                }
            }
            else
            {
                old_line = line;
            }
        }

        file.close();
    };

    void readFacesFile(const std::string& fileName)
    {
        std::ifstream file;
        std::string line;
        std::string old_line;
        file.open(fileName, std::ios_base::in);
        bool inDataBlock = false;
        int faceIdx = 0;

        if (!file.is_open())
        {
            std::cerr << "Error: 无法打开faces文件: " << fileName << std::endl;
            return;
        }

        while (std::getline(file, line))
        {
            if (line.find("(") != std::string::npos && !inDataBlock)
            {
                inDataBlock = true;
                int faceNum = std::stoi(old_line);
                facePoints.resize(faceNum);
                faceIdx = 0;
                continue;
            }

            if (line.find(")") != std::string::npos && inDataBlock &&
                line.find("(") == std::string::npos)
            {
                break;
            }

            if (inDataBlock && faceIdx < facePoints.size())
            {
                if (line.empty())
                    continue;

                size_t left = line.find("(");
                size_t right = line.find(")");
                if (left == std::string::npos || right == std::string::npos)
                    continue;

                std::string vertexNumStr = line.substr(0, left);
                int vertexNum = std::stoi(vertexNumStr);

                facePoints[faceIdx].reserve(vertexNum);
                std::string vertexIdStr =
                    line.substr(left + 1, right - left - 1);
                std::istringstream iss(vertexIdStr);

                int vertexIdx = 0;
                int vertexId;
                while (iss >> vertexId && vertexIdx < vertexNum)
                {
                    facePoints[faceIdx].push_back(vertexId);
                    vertexIdx++;
                }
                faceIdx++;
            }
            else
            {
                old_line = line;
            }
        }

        file.close();
    }

    void readBoundaryFile(const string& fileName)
    {
        std::ifstream file(fileName);
        if (!file.is_open())
        {
            std::cerr << "Error: 无法打开boundary文件: " << fileName
                      << std::endl;
            return;
        }

        std::string line;
        bool inDataBlock = false;
        while (std::getline(file, line))
        {
            if (line.find("(") != std::string::npos && !inDataBlock)
            {
                inDataBlock = true;
                continue;
            }

            if (inDataBlock)
            {
                if (line.empty() || line.find(")") != std::string::npos)
                    continue;

                // 假设行首是边界名称
                std::stringstream ss(line);
                std::string patchName;
                ss >> patchName;
                if (patchName.empty())
                    continue;

                BoundaryPatch patch;
                patch.name = patchName;

                // 读取直到找到 {
                while (std::getline(file, line) &&
                       line.find("{") == std::string::npos)
                    ;

                // 读取内部内容
                while (std::getline(file, line) &&
                       line.find("}") == std::string::npos)
                {
                    if (line.find("nFaces") != std::string::npos)
                    {
                        std::string sub = line.substr(line.find("nFaces") + 6);
                        patch.nFaces = std::stoi(sub.substr(0, sub.find(";")));
                    }
                    else if (line.find("startFace") != std::string::npos)
                    {
                        std::string sub =
                            line.substr(line.find("startFace") + 9);
                        patch.firstFaceIdx =
                            std::stoi(sub.substr(0, sub.find(";")));
                    }
                }
                boundary.push_back(patch);
            }
        }
        file.close();
    };

    // 通用的几何属性计算函数
    void calculateGeometry()
    {
        numCells = 0;
        for (int o : owner)
            numCells = std::max(numCells, o + 1);
        for (int n : neighbour)
            numCells = std::max(numCells, n + 1);

        int nFaces = facePoints.size();
        faceCentres.assign(nFaces, Vector(0, 0, 0));
        faceAreas.assign(nFaces, 0.0);
        faceNormals.assign(nFaces, Vector(0, 0, 0));

        for (int i = 0; i < nFaces; ++i)
        {
            Vector center(0, 0, 0);
            for (int vid : facePoints[i])
                center = center + points[vid];
            faceCentres[i] = center / (double)facePoints[i].size();

            // 2D 简化：计算面长度和法向量 (假设面由两个点组成)
            if (facePoints[i].size() >= 2)
            {
                Vector p1 = points[facePoints[i][0]];
                Vector p2 = points[facePoints[i][1]];
                Vector dP = p2 - p1;
                faceAreas[i] = dP.getMag();
                // 默认法向量 (dy, -dx, 0)
                faceNormals[i] = Vector(dP.y, -dP.x, 0.0);
                if (faceNormals[i].getMag() > 1e-15)
                    faceNormals[i] = faceNormals[i] / faceNormals[i].getMag();
            }
        }

        cellCentres.assign(numCells, Vector(0, 0, 0));
        cellVolumes.assign(numCells, 0.0);

        // 计算单元中心（平均面中心，简单处理）
        std::vector<int> faceCount(numCells, 0);
        for (int i = 0; i < nFaces; ++i)
        {
            cellCentres[owner[i]] = cellCentres[owner[i]] + faceCentres[i];
            faceCount[owner[i]]++;
            if (neighbour[i] != -1)
            {
                cellCentres[neighbour[i]] =
                    cellCentres[neighbour[i]] + faceCentres[i];
                faceCount[neighbour[i]]++;
            }
        }
        for (int i = 0; i < numCells; ++i)
        {
            if (faceCount[i] > 0)
                cellCentres[i] = cellCentres[i] / (double)faceCount[i];
        }

        // 修正法向量方向并计算单元体积 (2D 为面积)
        for (int i = 0; i < nFaces; ++i)
        {
            Vector vec = faceCentres[i] - cellCentres[owner[i]];
            if (vec.dotWith(faceNormals[i]) < 0)
                faceNormals[i] = faceNormals[i] * -1.0;

            // 累积单元体积贡献 (散度定理: V = 1/d * sum(x * n * A))
            double dot = (faceCentres[i] - cellCentres[owner[i]])
                             .dotWith(faceNormals[i]);
            cellVolumes[owner[i]] += 0.5 * dot * faceAreas[i];

            if (neighbour[i] != -1)
            {
                double dot_n = (faceCentres[i] - cellCentres[neighbour[i]])
                                   .dotWith(faceNormals[i] * -1.0);
                cellVolumes[neighbour[i]] += 0.5 * dot_n * faceAreas[i];
            }
        }
    }

    // 创建一个简单的 2D 正方形网格
    void createSquareMesh(int nx, int ny, double Lx, double Ly)
    {
        numCells = nx * ny;
        double dx = Lx / nx;
        double dy = Ly / ny;

        points.clear();
        facePoints.clear();
        owner.clear();
        neighbour.clear();
        boundary.clear();

        for (int j = 0; j <= ny; ++j)
        {
            for (int i = 0; i <= nx; ++i)
            {
                points.push_back(Vector(i * dx, j * dy, 0.0));
            }
        }

        auto getPointIdx = [&](int i, int j) { return j * (nx + 1) + i; };
        auto getCellIdx = [&](int i, int j) { return j * nx + i; };

        // 1. 内部面
        // 垂直内部面
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 1; i < nx; ++i)
            {
                facePoints.push_back(
                    {getPointIdx(i, j), getPointIdx(i, j + 1)});
                owner.push_back(getCellIdx(i - 1, j));
                neighbour.push_back(getCellIdx(i, j));
            }
        }
        // 水平内部面
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                facePoints.push_back(
                    {getPointIdx(i, j), getPointIdx(i + 1, j)});
                owner.push_back(getCellIdx(i, j - 1));
                neighbour.push_back(getCellIdx(i, j));
            }
        }
        nInternalFace = facePoints.size();

        // 2. 边界
        auto addPatch = [&](std::string name,
                            std::vector<std::vector<int>> fP,
                            std::vector<int> own)
        {
            BoundaryPatch patch;
            patch.name = name;
            patch.firstFaceIdx = facePoints.size();
            patch.nFaces = fP.size();
            for (int i = 0; i < fP.size(); ++i)
            {
                facePoints.push_back(fP[i]);
                owner.push_back(own[i]);
                // neighbour.push_back(-1);
            }
            boundary.push_back(patch);
        };

        // Left
        std::vector<std::vector<int>> leftFP;
        std::vector<int> leftOwn;
        for (int j = 0; j < ny; ++j)
        {
            leftFP.push_back({getPointIdx(0, j + 1), getPointIdx(0, j)});
            leftOwn.push_back(getCellIdx(0, j));
        }
        addPatch("left", leftFP, leftOwn);

        // Right
        std::vector<std::vector<int>> rightFP;
        std::vector<int> rightOwn;
        for (int j = 0; j < ny; ++j)
        {
            rightFP.push_back({getPointIdx(nx, j), getPointIdx(nx, j + 1)});
            rightOwn.push_back(getCellIdx(nx - 1, j));
        }
        addPatch("right", rightFP, rightOwn);

        // Down
        std::vector<std::vector<int>> downFP;
        std::vector<int> downOwn;
        for (int i = 0; i < nx; ++i)
        {
            downFP.push_back({getPointIdx(i + 1, 0), getPointIdx(i, 0)});
            downOwn.push_back(getCellIdx(i, 0));
        }
        addPatch("down", downFP, downOwn);

        // Top
        std::vector<std::vector<int>> topFP;
        std::vector<int> topOwn;
        for (int i = 0; i < nx; ++i)
        {
            topFP.push_back({getPointIdx(i, ny), getPointIdx(i + 1, ny)});
            topOwn.push_back(getCellIdx(i, ny - 1));
        }
        addPatch("top", topFP, topOwn);

        calculateGeometry();
    }
};
