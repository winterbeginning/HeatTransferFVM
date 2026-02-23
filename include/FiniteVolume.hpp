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
        double gf = (f.dotWith(d)) / (d.dotWith(d)); // 插值权重
        gf = std::max(0.0, std::min(1.0, gf));       // 限制在[0,1]

        faceField[i] = cellField[o] * (1.0 - gf) + cellField[n] * gf;
    }

    auto field = const_cast<Field<ValueType>&>(cellField);
    field.correctBoundaryField();

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
        faceFlow[i] = Uf[i].dotWith(mesh.faceNormals[i]) * mesh.faceAreas[i];
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
template <typename ValueType>
inline void Ddt(SpaceMatrix<ValueType>& phiEqn,
                const Properties& properties,
                const Field<ValueType>& phi,
                double dt)
{
    const Mesh& mesh = phi.mesh;
    for (int i = 0; i < mesh.numCells; ++i)
    {
        ValueType trans =
            mesh.cellVolumes[i] * properties.rho()[i] * properties.Cp()[i] / dt;
        phiEqn.addToA(i, i, trans);
        phiEqn.addTob(i, trans * phi.oldField()[i]);
    }
}

template <typename ValueType>
inline void Source(SpaceMatrix<ValueType>& Eqn, const Field<ValueType>& source)
{
    const Mesh& mesh = source.mesh;
    for (size_t i = 0; i < mesh.numCells; ++i)
    {
        Eqn.addTob(i, source[i] * mesh.cellVolumes[i]);
    }
}

inline double limiter_minmod(double r)
{
    return std::max(0.0, std::min(1.0, r));
}

inline double limiter_vanleer(double r)
{
    double r_abs = std::abs(r);
    return (r + r_abs) / (1 + r_abs);
}

inline double limiter_superbee(double r)
{
    return std::max({0.0, std::min(1.0, 2.0 * r), std::min(2.0, r)});
}

inline double
get_correction(Vector& gradUp, const Vector& rCf, double phiUp, double phiDown)
{
    double dPhi_linear = gradUp.dotWith(rCf);
    double dPhi_jump = phiDown - phiUp;

    // 避免除以极小值
    if (std::abs(dPhi_jump) < 1e-15)
        return 0.0;

    // r 的定义：上游预测斜率与实际台阶斜率的比值
    double r = dPhi_jump / (2.0 * dPhi_linear + 1e-15);

    // 使用限制器，注意：在梯度重构中，我们需要限制的是外推量
    // 这里采用一种更稳健的写法：直接限制面中心值在 [Up, Down] 之间
    double psi = limiter_minmod(r);

    // 如果预测方向与跳跃方向相反，psi 必为 0
    // 如果预测过大，psi 也会通过 r 减小
    return psi * dPhi_linear;
}

// template <typename ValueType>
//  修改 Div 以支持单元场和面质量通量
inline void Div(SpaceMatrix<double>& phiEqn,
                const Properties& properties,
                const FaceField<double>& faceFlux,
                const Field<double>& phi)
{
    const Mesh& mesh = phi.mesh;

    FaceField<double> faceRho(mesh);
    fvc::interpolate(properties.rho(), faceRho);

    Field<Vector> phiGrad = fvc::Grad(phi);

    // 1. 内部面 (Upwind)
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];
        double weightFlux = faceFlux[i] * faceRho[i];

        int up = (weightFlux > 0) ? o : n;
        int down = (weightFlux > 0) ? n : o;

        // 核心：基于梯度的线性重构
        // phi_f = phi_up + psi * (grad_up . r_cf)
        Vector rCf = mesh.faceCentres[i] - mesh.cellCentres[up];
        double correction =
            get_correction(phiGrad[up], rCf, phi[up], phi[down]);

        // A 矩阵保持一阶上风，保证对角优势
        phiEqn.addToA(o, up, weightFlux);  // 对o单元，外法相与Sf一致
        phiEqn.addToA(n, up, -weightFlux); // 对o单元，外法相与Sf相反

        // 高阶修正项放进 b 向量（延迟修正），放进 b 移项变为负
        phiEqn.addTob(o, weightFlux * (-correction));
        // 对o单元，外法相与Sf一致
        phiEqn.addTob(n, -weightFlux * (-correction));
        // 对o单元，外法相与Sf相反
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
            double flux = faceFlux[globalFaceId] * faceRho[globalFaceId];

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

inline void Laplacian(SpaceMatrix<double>& phiEqn,
                      const Field<double>& diffusionCoeff,
                      const Field<double>& phi,
                      bool enableNonOrthCorrection = true)
{
    const Mesh& mesh = phi.mesh;

    FaceField<double> faceDiffusionCoeff(mesh);
    fvc::interpolate(diffusionCoeff, faceDiffusionCoeff);

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
        double SfdotD = Sf.dotWith(d);
        double deltaCoeff;
        if (enableNonOrthCorrection)
        {
            // 非正交修正模式：deltaCoeff = κ|Sf|² / (Sf·d)
            deltaCoeff =
                faceDiffusionCoeff[i] * magSf * magSf / (SfdotD + 1e-30);
        }
        else
        {
            // 简化正交模式：deltaCoeff = κ|Sf| / |d|
            deltaCoeff = faceDiffusionCoeff[i] * magSf / d.getMag();
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
            Vector Tf = Sf - Delta;

            // 非正交修正项：κ * (∇φ)_f · T （显式，使用插值梯度）
            double nonOrthCorr =
                faceDiffusionCoeff[i] * (phiGrad[i].dotWith(Tf));

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
            double SfdotD = Sf.dotWith(d);
            double deltaCoeff =
                faceDiffusionCoeff[i] * magSf * magSf / (SfdotD + 1e-30);

            // 应用混合边界条件：fraction * fixedValue + (1-fraction) *
            // fixedGradient
            double coeff_implicit = deltaCoeff * fvPatch.fraction;
            double source_value =
                deltaCoeff * fvPatch.fraction * fvPatch.refValue;
            double source_grad = faceDiffusionCoeff[i] * magSf *
                                 (1.0 - fvPatch.fraction) * fvPatch.refGrad;

            phiEqn.addToA(o, o, coeff_implicit);
            phiEqn.addTob(o, source_value + source_grad);
        }
    }
}

inline void Laplacian(SpaceMatrix<Vector>& phiEqn,
                      const Field<double>& diffusionCoeff,
                      const Field<Vector>& phi,
                      bool enableNonOrthCorrection = true)
{
    const Mesh& mesh = phi.mesh;

    FaceField<double> faceDiffusionCoeff(mesh);
    fvc::interpolate(diffusionCoeff, faceDiffusionCoeff);

    // 计算梯度场（用于非正交修正）
    FaceField<Tensor> phiGrad(mesh);
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
        double SfdotD = Sf.dotWith(d);
        Vector deltaCoeff;
        if (enableNonOrthCorrection)
        {
            // 非正交修正模式：deltaCoeff = κ|Sf|² / (Sf·d)
            deltaCoeff = Vector(1, 1, 1) * faceDiffusionCoeff[i] * magSf *
                         magSf / (SfdotD + 1e-30);
        }
        else
        {
            // 简化正交模式：deltaCoeff = κ|Sf| / |d|
            deltaCoeff =
                Vector(1, 1, 1) * faceDiffusionCoeff[i] * magSf / d.getMag();
        }

        // 添加到系数矩阵（隐式）
        phiEqn.addToA(o, o, deltaCoeff);
        phiEqn.addToA(o, n, -1.0 * deltaCoeff);
        phiEqn.addToA(n, n, deltaCoeff);
        phiEqn.addToA(n, o, -1.0 * deltaCoeff);

        // === 非正交修正（显式处理，作为源项）===
        if (enableNonOrthCorrection)
        {
            // 计算非正交向量：T = Sf - Δ，其中 Δ = (|Sf|² / (Sf · d)) * d
            Vector Delta = d * (magSf * magSf / (SfdotD + 1e-30));
            Vector Tf = Sf - Delta;

            // 非正交修正项：κ * (∇φ)_f · T （显式，使用插值梯度）
            Vector nonOrthCorr =
                faceDiffusionCoeff[i] * (phiGrad[i].multiply(Tf));

            // 添加到源项（扩散项在方程左边，源项需要取负）
            phiEqn.addTob(o, nonOrthCorr);
            phiEqn.addTob(n, -1.0 * nonOrthCorr);
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
            double SfdotD = Sf.dotWith(d);
            Vector deltaCoeff = Vector(1, 1, 1) * faceDiffusionCoeff[i] *
                                magSf * magSf / (SfdotD + 1e-30);

            // 应用混合边界条件：fraction * fixedValue + (1-fraction) *
            // fixedGradient
            Vector coeff_implicit = deltaCoeff * fvPatch.fraction;
            Vector source_value =
                deltaCoeff * fvPatch.fraction * fvPatch.refValue;
            Vector source_grad = faceDiffusionCoeff[i] * magSf *
                                 (1.0 - fvPatch.fraction) * fvPatch.refGrad;

            phiEqn.addToA(o, o, coeff_implicit);
            phiEqn.addTob(o, source_value + source_grad);
        }
    }
}

} // namespace fvm

#endif