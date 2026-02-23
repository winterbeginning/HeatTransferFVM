#ifndef _FVMDiscrete_
#define _FVMDiscrete_

#include "SpaceMatrix.hpp"
#include "Field.hpp"
#include "Properties.hpp"
#include "FaceField.hpp"
#include "Tensor3D.hpp"

namespace fvc
{
// 基础线性插值（用于梯度计算和非正交网格）
template <typename ValueType>
inline void interpolate(const Field<ValueType>& cellField,
                        FaceField<ValueType>& faceField)
{
    const Mesh& mesh = cellField.mesh;
    // 内部面：距离加权线性插值
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];

        // 计算插值因子：基于面中心到单元中心的距离
        Vector d = mesh.cellCentres[n] - mesh.cellCentres[o];
        Vector f = mesh.faceCentres[i] - mesh.cellCentres[o];
        double gf = (f * d) / (d * d);         // 插值权重
        gf = std::max(0.0, std::min(1.0, gf)); // 限制在[0,1]

        faceField[i] = cellField[o] * (1.0 - gf) + cellField[n] * gf;
    }

    // 边界：读取已维护的边界场值
    for (int i = mesh.nInternalFace; i < mesh.numFaces; ++i)
    {
        faceField[i] = cellField.getBoundaryFaceValue(i);
    }
}

// 计算面通量 phi = U.dot(S)
inline FaceField<double> flux(const Field<Vector>& U)
{
    const Mesh& mesh = U.mesh;
    FaceField<double> faceFlow(mesh);
    FaceField<Vector> Uf(mesh);

    interpolate(U, Uf);

    for (size_t i = 0; i < (size_t)mesh.facePoints.size(); ++i)
    {
        faceFlow[i] = Uf[i] * mesh.faceNormals[i] * mesh.faceAreas[i];
    }
    return faceFlow;
}

inline Field<Vector> Grad(const Field<double>& phi)
{
    const Mesh& mesh = phi.mesh;
    Field<Vector> g(mesh, Vector(0, 0, 0));

    FaceField<double> phi_f(mesh);
    interpolate(phi, phi_f);

    for (size_t i = 0; i < (size_t)mesh.facePoints.size(); ++i)
    {
        int o = mesh.owner[i];
        Vector flux = mesh.faceNormals[i] * mesh.faceAreas[i] * phi_f[i];

        g.internalField[o] = g.internalField[o] + flux;
        if (i < mesh.nInternalFace)
        {
            int n = mesh.neighbour[i];
            g.internalField[n] = g.internalField[n] - flux;
        }
    }

    for (int i = 0; i < mesh.numCells; ++i)
    {
        g.internalField[i] = g.internalField[i] / (mesh.cellVolumes[i] + 1e-20);
    }
    return g;
}

inline Field<Tensor> Grad(const Field<Vector>& phi)
{
    const Mesh& mesh = phi.mesh;
    Field<Tensor> g(mesh, Tensor());

    FaceField<Vector> phi_f(mesh);
    interpolate(phi, phi_f);

    for (int i = 0; i < mesh.numFaces; ++i)
    {
        int o = mesh.owner[i];
        Vector faceVector =
            mesh.faceNormals[i] * mesh.faceAreas[i]; // S_f = n_f * A_f
        Tensor flux = faceVector.outProductWith(phi_f[i]);

        g.internalField[o] = g.internalField[o] + flux;
        if (i < mesh.nInternalFace)
        {
            int n = mesh.neighbour[i];
            g.internalField[n] = g.internalField[n] - flux;
        }
    }

    for (int i = 0; i < mesh.numCells; ++i)
    {
        g.internalField[i] = g.internalField[i] / (mesh.cellVolumes[i] + 1e-20);
    }
    return g;
}
} // namespace fvc

namespace fvm
{
inline void Ddt(SpaceMatrix& phiEqn,
                const Properties& properties,
                const Field<double>& phi,
                double dt)
{
    const Mesh& mesh = phi.mesh;
    for (int i = 0; i < mesh.numCells; ++i)
    {
        double trans =
            mesh.cellVolumes[i] * properties.rho()[i] * properties.Cp()[i] / dt;
        phiEqn.addToA(i, i, trans);
        phiEqn.addTob(i, trans * phi.oldField()[i]);
    }
}

// 修改 Div 以支持单元场和面质量通量
inline void Div(SpaceMatrix& phiEqn,
                const Properties& properties,
                const FaceField<double>& faceFlux,
                const Field<double>& phi)
{
    const Mesh& mesh = phi.mesh;

    FaceField<double> faceRho(mesh);
    fvc::interpolate(properties.rho(), faceRho);
    // 1. 内部面 (Upwind)
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];
        double flux = faceFlux[i] * faceRho[i];

        if (flux > 0) // 流动方向 o --> n   面变量等于o单元
        {
            phiEqn.addToA(o, o, flux); // o单元方程，面变量来自o，流出为正
            phiEqn.addToA(n, o, -flux); // n单元方程，面变量来自o，流入为负
        }
        else // 流动方向 n --> o   面变量等于n单元
        {
            phiEqn.addToA(o, n, flux); // o单元方程，面变量来自n，流入为负
            phiEqn.addToA(n, n, -flux); // n单元方程，面变量来自n，流出为正
        };
    }

    // 2. 边界面 (根据边界流向决定使用单元值还是边界值)
    for (int pIdx = 0; pIdx < (int)phi.boundaryList.size(); ++pIdx)
    {
        const auto& fvPatch = phi.boundaryList[pIdx];

        if (fvPatch.boundaryType == BoundaryType::EMPTY)
            continue;

        for (int i = 0; i < fvPatch.nFaces; ++i)
        {
            int globalFaceId = fvPatch.firstFaceIdx + i;
            int o = mesh.owner[globalFaceId];
            double flux = faceFlux[globalFaceId] * faceRho[i];

            if (flux > 0) // 流出
            {
                phiEqn.addToA(o, o, flux);
            }
            else // 流入 (通过 b 累加背景源项)
            {
                double phi_face = phi.getBoundaryFaceValue(globalFaceId);
                phiEqn.addTob(o, -flux * phi_face);
            }
        }
    }
}

inline void Laplacian(SpaceMatrix& phiEqn,
                      const Properties& properties,
                      const Field<double>& phi,
                      bool enableNonOrthCorrection = true)
{
    const Mesh& mesh = phi.mesh;

    FaceField<double> faceKappa(mesh);
    fvc::interpolate(properties.kappa(), faceKappa);

    // 计算梯度场（用于非正交修正）
    FaceField<Vector> phiGrad(mesh);
    if (enableNonOrthCorrection)
        fvc::interpolate(fvc::Grad(phi), phiGrad);

    // 1. 内部面扩散
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];

        // 几何向量
        Vector d = mesh.cellCentres[n] - mesh.cellCentres[o]; // 单元中心连线
        Vector Sf = mesh.faceNormals[i] * mesh.faceAreas[i]; // 面向量
        double magSf = mesh.faceAreas[i];                    // 面积

        // === 正交部分（隐式处理，进入系数矩阵）===
        double SfdotD = Sf * d;
        double deltaCoeff;
        if (enableNonOrthCorrection)
        {
            // 非正交修正模式：deltaCoeff = κ|Sf|² / (Sf·d)
            deltaCoeff = faceKappa[i] * magSf * magSf / (SfdotD + 1e-30);
        }
        else
        {
            // 简化正交模式：deltaCoeff = κ|Sf| / |d|
            deltaCoeff = faceKappa[i] * magSf / d.getMag();
        }

        // 添加到系数矩阵（隐式）
        phiEqn.addToA(o, o, deltaCoeff);
        phiEqn.addToA(o, n, -deltaCoeff);
        phiEqn.addToA(n, n, deltaCoeff);
        phiEqn.addToA(n, o, -deltaCoeff);

        // === 非正交修正（显式处理，作为源项）===
        if (enableNonOrthCorrection)
        {
            // 计算非正交向量：T = Sf - Δ，其中 Δ = (|Sf|² / (Sf · d)) * d
            Vector Delta = d * (magSf * magSf / (SfdotD + 1e-30));
            Vector T = Sf - Delta;

            // 非正交修正项：κ * (∇φ)_f · T （显式，使用插值梯度）
            double nonOrthCorr = faceKappa[i] * (phiGrad[i] * T);

            // 添加到源项（扩散项在方程左边，源项需要取负）
            phiEqn.addTob(o, nonOrthCorr);
            phiEqn.addTob(n, -nonOrthCorr);
        }
    }

    // 2. 边界扩散
    for (int pIdx = 0; pIdx < (int)mesh.boundary.size(); ++pIdx)
    {
        const auto& fvPatch = phi.boundaryList[pIdx];

        if (fvPatch.boundaryType == BoundaryType::EMPTY)
            continue;

        for (int i = fvPatch.firstFaceIdx;
             i < fvPatch.firstFaceIdx + fvPatch.nFaces;
             ++i)
        {
            int o = mesh.owner[i];

            // 几何向量
            Vector d = mesh.faceCentres[i] - mesh.cellCentres[o];
            Vector Sf = mesh.faceNormals[i] * mesh.faceAreas[i];
            double magSf = mesh.faceAreas[i];

            // 边界的正交系数：deltaCoeff = |Sf|² / (Sf · d)
            double SfdotD = Sf * d;
            double deltaCoeff = faceKappa[i] * magSf * magSf / (SfdotD + 1e-30);

            // 应用混合边界条件：fraction * fixedValue + (1-fraction) *
            // fixedGradient
            double coeff_implicit = deltaCoeff * fvPatch.fraction;
            double source_value =
                deltaCoeff * fvPatch.fraction * fvPatch.refValue;
            double source_grad = faceKappa[i] * magSf *
                                 (1.0 - fvPatch.fraction) * fvPatch.refGrad;

            phiEqn.addToA(o, o, coeff_implicit);
            phiEqn.addTob(o, source_value + source_grad);
        }
    }
}

} // namespace fvm

#endif