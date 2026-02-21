#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include "Mesh.hpp"

// 边界条件结构

// 场类：模拟 OpenFOAM 的 volField
template <typename valType>
class Field
{
    struct BoundaryCondition
    {
        valType refValue; // 参考值 (如 T_left)
        valType refGrad;  // 参考梯度 (如 q/k)
        double fraction;  // 1.0 = Dirichlet, 0.0 = Neumann
    };

public:
    std::string name;
    const Mesh& mesh;

    // 内部场 (Internal field): cell-centered values
    std::vector<valType> internalField;

    // 边界场 (Boundary field): 与 Mesh 中的 boundary 一一对应
    std::vector<BoundaryCondition> boundaryField;

    Field(std::string name, const Mesh& m, valType initVal = 0.0)
        : name(name), mesh(m)
    {
        internalField.assign(mesh.numCells, initVal);
        boundaryField.resize(mesh.boundary.size());
    }

    void fill(const valType& val)
    {
        std::fill(internalField.begin(), internalField.end(), val);
    };

    // 设置指定边界的条件
    void setBoundary(const std::string& patchName,
                     valType refV,
                     valType refG,
                     double frac)
    {
        for (size_t i = 0; i < mesh.boundary.size(); ++i)
        {
            if (mesh.boundary[i].name == patchName)
            {
                boundaryField[i] = {refV, refG, frac};
                return;
            }
        }
        throw std::runtime_error("Patch " + patchName + " not found in mesh!");
    }

    // 获取特定单元的值 (辅助函数)
    valType& operator[](int idx)
    {
        return internalField[idx];
    }
    const valType& operator[](int idx) const
    {
        return internalField[idx];
    }

    // 获取场的大小
    size_t size() const
    {
        return internalField.size();
    }
};
