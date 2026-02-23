#ifndef _Mesh_
#define _Mesh_

#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include "Point3D.hpp"

enum BoundaryType
{
    PATCH,
    WALL,
    EMPTY
};

// 定义边界组（Boundary Patches）
struct BoundaryPatch
{
    std::string name;
    int firstFaceIdx;
    int nFaces;
    BoundaryType boundaryType;

    BoundaryPatch(const std::string& n, int first, int num)
        : name(n), firstFaceIdx(first), nFaces(num)
    {
    }
    BoundaryPatch() : firstFaceIdx(0), nFaces(0)
    {
    }
};

inline std::string trim(const std::string& s)
{
    auto start = s.find_first_not_of(" \t\n\r");
    auto end = s.find_last_not_of(" \t\n\r");
    if (start == std::string::npos)
        return "";
    return s.substr(start, end - start + 1);
}

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
            // 去除行首尾的空白字符（空格、制表符、换行符等），避免空行误判
            line.erase(0, line.find_first_not_of(" \t\n\r"));
            line.erase(line.find_last_not_of(" \t\n\r") + 1);

            // 检测数据块开始（左括号）
            if (line.find("(") != std::string::npos && !inDataBlock)
            {
                inDataBlock = true;
                // 预分配空间（使用括号前一行的数字）
                if (!old_line.empty())
                {
                    owner.reserve(std::stoi(old_line));
                }
                continue;
            }

            // 检测数据块结束（右括号）
            if (line.find(")") != std::string::npos && inDataBlock)
            {
                break;
            }

            // 处理数据块内的内容
            if (inDataBlock)
            {
                if (line.empty())
                    continue;

                std::istringstream iss(line);
                std::string numStr;
                while (iss >> numStr)
                {
                    try
                    {
                        owner.push_back(std::stoi(numStr));
                    }
                    catch (const std::invalid_argument& e)
                    {
                        std::cerr << "警告: 无效的数字格式: " << numStr
                                  << "，已跳过" << std::endl;
                    }
                    catch (const std::out_of_range& e)
                    {
                        std::cerr << "警告: 数字超出范围: " << numStr
                                  << "，已跳过" << std::endl;
                    }
                }
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
            line = trim(line);

            if (line.find("(") != std::string::npos && !inDataBlock)
            {
                inDataBlock = true;
                // 预分配空间（使用括号前一行的数字）
                if (!old_line.empty())
                {
                    owner.reserve(std::stoi(old_line));
                }
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

                std::istringstream iss(line);
                std::string numStr;
                while (iss >> numStr)
                {
                    try
                    {
                        neighbour.push_back(std::stoi(numStr));
                    }
                    catch (const std::invalid_argument& e)
                    {
                        std::cerr << "警告: 无效的数字格式: " << numStr
                                  << "，已跳过" << std::endl;
                    }
                    catch (const std::out_of_range& e)
                    {
                        std::cerr << "警告: 数字超出范围: " << numStr
                                  << "，已跳过" << std::endl;
                    }
                }
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
                    if (line.find("type") != std::string::npos)
                    {
                        std::string sub = line.substr(line.find("type") + 4);
                        std::string ValueTypeue =
                            trim(sub.substr(0, sub.find(";")));
                        if (ValueTypeue == "patch")
                            patch.boundaryType = BoundaryType::PATCH;
                        else if (ValueTypeue == "wall")
                            patch.boundaryType = BoundaryType::WALL;
                        else
                            patch.boundaryType = BoundaryType::EMPTY;
                    }
                    else if (line.find("nFaces") != std::string::npos)
                    {
                        std::string sub = line.substr(line.find("nFaces") + 6);
                        patch.nFaces =
                            std::stoi(trim(sub.substr(0, sub.find(";"))));
                    }
                    else if (line.find("startFace") != std::string::npos)
                    {
                        std::string sub =
                            line.substr(line.find("startFace") + 9);
                        patch.firstFaceIdx =
                            std::stoi(trim(sub.substr(0, sub.find(";"))));
                    }
                }
                boundary.push_back(patch);
            }
        }
        file.close();
    };

    // 通用的几何属性计算函数，支持二维(Z=0)和三维网格
    // 改进版：使用迭代优化方法计算单元中心，提高对复杂网格的鲁棒性
    void calculateGeometry()
    {
        numCells = 0;
        for (int o : owner)
            numCells = std::max(numCells, o + 1);
        for (int n : neighbour)
            if (n != -1)
                numCells = std::max(numCells, n + 1);

        // 1. 自动判断维度
        bool is2D = true;
        for (const auto& fp : facePoints)
        {
            if (fp.size() > 2)
            {
                is2D = false;
                break;
            }
        }
        if (is2D)
        {
            for (const auto& p : points)
            {
                if (std::abs(p.z) > 1e-12)
                {
                    is2D = false;
                    break;
                }
            }
        }

        int nFaces = facePoints.size();
        faceCentres.assign(nFaces, Vector(0, 0, 0));
        faceAreas.assign(nFaces, 0.0);
        faceNormals.assign(nFaces, Vector(0, 0, 0));

        // 2. 计算面属性（面中心、面积、法向量）
        for (int i = 0; i < nFaces; ++i)
        {
            const auto& verts = facePoints[i];

            // 面中心：顶点平均
            Vector center(0, 0, 0);
            for (int vid : verts)
                center = center + points[vid];
            faceCentres[i] = center / (double)verts.size();

            if (is2D && verts.size() == 2)
            {
                // 二维：面是线段
                Vector p1 = points[verts[0]];
                Vector p2 = points[verts[1]];
                Vector dP = p2 - p1;
                faceAreas[i] = dP.getMag();
                faceNormals[i] = Vector(dP.y, -dP.x, 0.0);
            }
            else
            {
                // 三维：多边形面，使用扇形三角化计算矢量面积
                Vector vectorArea(0, 0, 0);
                for (size_t v = 0; v < verts.size(); ++v)
                {
                    Vector v1 = points[verts[v]] - faceCentres[i];
                    Vector v2 =
                        points[verts[(v + 1) % verts.size()]] - faceCentres[i];
                    vectorArea = vectorArea + v1.crossWith(v2);
                }
                vectorArea = vectorArea * 0.5;
                faceAreas[i] = vectorArea.getMag();
                faceNormals[i] = vectorArea / (faceAreas[i] + 1e-20);
            }

            // 归一化法向量
            if (faceNormals[i].getMag() > 1e-15)
                faceNormals[i] = faceNormals[i] / faceNormals[i].getMag();
        }

        // 3. 建立单元-面关系表（用于迭代优化单元中心）
        cellCentres.assign(numCells, Vector(0, 0, 0));
        cellVolumes.assign(numCells, 0.0);
        std::vector<std::vector<int>> cellFaces(numCells);

        for (int i = 0; i < nFaces; ++i)
        {
            cellFaces[owner[i]].push_back(i);
        }
        for (size_t i = 0; i < neighbour.size(); ++i)
        {
            cellFaces[neighbour[i]].push_back(i);
        }

        // 4. 初始估计：面中心平均
        for (int cellId = 0; cellId < numCells; ++cellId)
        {
            Vector sum(0, 0, 0);
            for (int faceId : cellFaces[cellId])
            {
                sum = sum + faceCentres[faceId];
            }
            if (cellFaces[cellId].size() > 0)
                cellCentres[cellId] = sum / (double)cellFaces[cellId].size();
        }

        // 5. 迭代优化单元中心：使用体积加权的金字塔中心（3次迭代）
        for (int iter = 0; iter < 3; ++iter)
        {
            for (int cellId = 0; cellId < numCells; ++cellId)
            {
                Vector volWeightedCenter(0, 0, 0);
                double totalVol = 0.0;

                for (int faceId : cellFaces[cellId])
                {
                    const auto& verts = facePoints[faceId];

                    if (is2D && verts.size() == 2)
                    {
                        // 2D：计算三角形（单元中心+两个顶点）
                        Vector p1 = points[verts[0]];
                        Vector p2 = points[verts[1]];
                        Vector triCenter =
                            (cellCentres[cellId] + p1 + p2) / 3.0;

                        // 三角形面积 = 0.5 * |叉积的Z分量|
                        double triArea =
                            0.5 *
                            std::abs((p2 - p1)
                                         .crossWith(cellCentres[cellId] - p1)
                                         .z);

                        volWeightedCenter =
                            volWeightedCenter + triCenter * triArea;
                        totalVol += triArea;
                    }
                    else
                    {
                        // 3D：将面分解为三角形，每个三角形与单元中心组成四面体
                        for (size_t v = 0; v < verts.size(); ++v)
                        {
                            Vector p1 = points[verts[v]];
                            Vector p2 = points[verts[(v + 1) % verts.size()]];
                            Vector p3 = faceCentres[faceId];

                            // 四面体中心 = (顶点0 + 顶点1 + 顶点2 + 顶点3) / 4
                            Vector tetCenter =
                                (cellCentres[cellId] + p1 + p2 + p3) * 0.25;

                            // 四面体体积 = |det(v1, v2, v3)| / 6
                            Vector v1 = p1 - cellCentres[cellId];
                            Vector v2 = p2 - cellCentres[cellId];
                            Vector v3 = p3 - cellCentres[cellId];
                            double tetVol =
                                std::abs(v1 * v2.crossWith(v3)) / 6.0;

                            volWeightedCenter =
                                volWeightedCenter + tetCenter * tetVol;
                            totalVol += tetVol;
                        }
                    }
                }

                if (totalVol > 1e-20)
                    cellCentres[cellId] = volWeightedCenter / totalVol;
            }
        }

        // 6. 修正面法向量方向（确保从owner指向neighbour）
        for (int i = 0; i < nFaces; ++i)
        {
            Vector dC = faceCentres[i] - cellCentres[owner[i]];
            if (dC * faceNormals[i] < 0)
                faceNormals[i] = faceNormals[i] * -1.0;
        }

        // 7. 计算单元体积：使用散度定理
        double dimInv = is2D ? 0.5 : (1.0 / 3.0);
        cellVolumes.assign(numCells, 0.0);

        for (int i = 0; i < nFaces; ++i)
        {
            // Owner单元的体积贡献
            double dot_o =
                (faceCentres[i] - cellCentres[owner[i]]) * faceNormals[i];
            cellVolumes[owner[i]] += dimInv * dot_o * faceAreas[i];

            // Neighbour单元的体积贡献（法向量反向）
            if (i < neighbour.size())
            {
                double dot_n = (faceCentres[i] - cellCentres[neighbour[i]]) *
                               faceNormals[i] * -1.0;
                cellVolumes[neighbour[i]] += dimInv * dot_n * faceAreas[i];
            }
        }

        // 8. 几何验证：检查异常情况
        int negVolCount = 0;
        for (int i = 0; i < numCells; ++i)
        {
            if (cellVolumes[i] < 0)
            {
                cellVolumes[i] = std::abs(cellVolumes[i]);
                negVolCount++;
            }
            if (cellVolumes[i] < 1e-20)
            {
                std::cerr << "警告: 单元 " << i
                          << " 体积过小: " << cellVolumes[i] << std::endl;
            }
        }

        if (negVolCount > 0)
        {
            std::cerr << "警告: 检测到 " << negVolCount
                      << " 个负体积单元（已修正为绝对值）" << std::endl;
            std::cerr << "       可能原因：网格质量差、法向量方向错误或强烈扭曲"
                      << std::endl;
        }

        // 验证面积
        for (int i = 0; i < nFaces; ++i)
        {
            if (faceAreas[i] < 1e-20)
            {
                std::cerr << "警告: 面 " << i << " 面积过小: " << faceAreas[i]
                          << std::endl;
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
            BoundaryPatch patch(name, facePoints.size(), fP.size());
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

#endif