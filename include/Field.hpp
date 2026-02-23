#ifndef _Field_
#define _Field_

#include <vector>
#include <string>
#include <algorithm>
#include "Mesh.hpp"

// 边界条件结构

// 场类：模拟 OpenFOAM 的 volField
template <typename ValueType>
class Field
{
public:
    // 边界条件定义 (Boundary Condition Definition)
    struct BoundaryPatchDefinition : public BoundaryPatch
    {
        ValueType refValue; // 给定的参考值 (Dirichlet)
        ValueType refGrad;  // 给定的参考梯度 (Neumann)
        double fraction;    // 1: Dirichlet, 0: Neumann
        BoundaryPatchDefinition(const BoundaryPatch& patch)
            : BoundaryPatch(patch),
              refValue(ValueType{}),
              refGrad(ValueType{}),
              fraction(0.0)
        {
        }

        void setBoundary(ValueType refV, ValueType refG, double frac)
        {
            this->refValue = refV;
            this->refGrad = refG;
            this->fraction = frac;
        }
    };

    const Mesh& mesh;

    // 内部场 (Internal Field): 单元值
    std::vector<ValueType> internalField;
    std::vector<ValueType> oldInternalField;

    // 边界条件参数
    std::vector<BoundaryPatchDefinition> boundaryList;

    // 缓存的边界面值 (Boundary Face Values): size = nBoundaryFaces
    // 对应索引从 mesh.nInternalFace 到 mesh.facePoints.size() - 1
    std::vector<ValueType> boundaryField;

    Field(const Mesh& m, ValueType initVal = ValueType{}) : mesh(m)
    {
        internalField.assign(mesh.numCells, initVal);
        oldInternalField = internalField;

        int nBoundaryFaces = (int)mesh.facePoints.size() - mesh.nInternalFace;
        boundaryField.assign(nBoundaryFaces, initVal);

        boundaryList.reserve(mesh.boundary.size());
        for (auto& patch : mesh.boundary)
        {
            BoundaryPatchDefinition fvPatch = BoundaryPatchDefinition(patch);
            if (fvPatch.boundaryType == BoundaryType::EMPTY)
                setBoundary(fvPatch, ValueType{}, ValueType{}, 0.0);
            boundaryList.emplace_back(fvPatch);
        }
    }

    void fill(const ValueType& value)
    {
        std::fill(internalField.begin(), internalField.end(), value);
    };

    // 设置边界条件的参数
    void setBoundary(const std::string& patchName,
                     ValueType refV,
                     ValueType refG,
                     double frac)
    {
        for (auto& bc : boundaryList)
        {
            if (bc.name == patchName)
            {
                bc.refValue = refV;
                bc.refGrad = refG;
                bc.fraction = frac;
                correctBoundaryField(); // 更新该边界的值
                return;
            }
        }
        throw std::runtime_error("Patch " + patchName + " not found!");
    }

    void setBoundary(BoundaryPatchDefinition& patch,
                     ValueType refV,
                     ValueType refG,
                     double frac)
    {
        patch.refValue = refV;
        patch.refGrad = refG;
        patch.fraction = frac;
        correctBoundaryField();
    }

    // 更新整个内部边界场 (基于当前单元中心值和 BC 参数)
    void correctBoundaryField()
    {
        for (const auto& bc : boundaryList)
        {
            for (int i = 0; i < bc.nFaces; ++i)
            {
                int globalFaceId = bc.firstFaceIdx + i;
                int o = mesh.owner[globalFaceId];
                double dist = mesh.faceCentres[globalFaceId].getDistance(
                    mesh.cellCentres[o]);

                // 线性外推基础公式: phi_f = f * phi_ref + (1-f) * (phi_cell +
                // grad_ref * d)
                ValueType phi_face = bc.fraction * bc.refValue +
                                     (1.0 - bc.fraction) *
                                         (internalField[o] + bc.refGrad * dist);

                boundaryField[globalFaceId - mesh.nInternalFace] = phi_face;
            }
        }
    }

    // 获取特定边界 Patch 的 face values (用于插值/输出)
    std::vector<ValueType> getPatchValues(int patchIdx) const
    {
        const auto& bc = boundaryList[patchIdx];
        std::vector<ValueType> vals(bc.nFaces);
        for (int i = 0; i < bc.nFaces; ++i)
        {
            vals[i] = boundaryField[bc.firstFaceIdx + i - mesh.nInternalFace];
        }
        return vals;
    }

    // 获取特定全局面索引的边界值
    const ValueType& getBoundaryFaceValue(int globalFaceId) const
    {
        return boundaryField[globalFaceId - mesh.nInternalFace];
    }

    // 获取特定单元的值 (辅助函数)
    ValueType& operator[](int idx)
    {
        return internalField[idx];
    }
    const ValueType& operator[](int idx) const
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

    const std::vector<ValueType>& oldField() const
    {
        return oldInternalField;
    }
};
#endif