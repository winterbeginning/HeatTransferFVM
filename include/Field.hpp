#ifndef _Field_
#define _Field_

#include <vector>
#include <string>
#include <algorithm>
#include "Mesh.hpp"

// 边界条件结构

// 场类：模拟 OpenFOAM 的 volField
template <typename valType>
class Field
{
    struct BoundaryCondition : public BoundaryPatch
    {
        valType refValue; // 参考值 (如 T_left)
        valType refGrad;  // 参考梯度 (如 q/k)
        double fraction;  // 1.0 = Dirichlet, 0.0 = Neumann
        BoundaryCondition(const BoundaryPatch& patch) : BoundaryPatch(patch)
        {
        }
    };

public:
    const Mesh& mesh;

    // 内部场 (Internal field): cell-centered values
    std::vector<valType> internalField;

    std::vector<valType> oldInternalField;

    // 边界场 (Boundary field): 与 Mesh 中的 boundary 一一对应
    std::vector<BoundaryCondition> boundaryList;

    std::vector<vector<valType>> boundaryField;

    Field(const Mesh& m, valType initVal = valType{}) : mesh(m)
    {
        internalField.assign(mesh.numCells, initVal);
        oldInternalField = internalField;

        boundaryList.reserve(mesh.boundary.size());
        boundaryField.reserve(mesh.boundary.size());

        // 遍历网格边界，逐个初始化
        for (const auto& patch : mesh.boundary)
        {
            boundaryList.emplace_back(patch);
            boundaryField.emplace_back(patch.nFaces, initVal);
        }
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
            if (boundaryList[i].name == patchName)
            {
                boundaryList[i].refValue = refV;
                boundaryList[i].refGrad = refG;
                boundaryList[i].fraction = frac;
                boundaryField[i] = getBoundaryValue(i);
                return;
            }
        }
        throw std::runtime_error("Patch " + patchName + " not found in mesh!");
    }

    void correctBoundaryField()
    {
        for (size_t i = 0; i < boundaryList.size(); ++i)
        {
            boundaryField[i] = getBoundaryValue(i);
        }
    }

    vector<valType> getBoundaryValue(int i) const
    {
        int nFaces = boundaryList[i].nFaces;
        vector<valType> tmpValue(nFaces);

        for (int j = 0; j < nFaces; ++j)
        {
            int globalFaceId = j + boundaryList[i].firstFaceIdx;
            int o = mesh.owner[globalFaceId];
            double distance =
                (mesh.faceCentres[globalFaceId] - mesh.cellCentres[o]).getMag();
            valType phi_face =
                boundaryList[i].fraction * boundaryList[i].refValue +
                (1.0 - boundaryList[i].fraction) *
                    (this->internalField[o] +
                     boundaryList[i].refGrad * distance);
            tmpValue[j] = phi_face;
        }
        return tmpValue;
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

    void storeOldField()
    {
        oldInternalField = internalField;
    }

    const std::vector<valType>& oldField() const
    {
        return oldInternalField;
    }
};
#endif