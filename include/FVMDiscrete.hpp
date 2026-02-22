#ifndef _FVMDiscrete_
#define _FVMDiscrete_

#include "SpaceMatrix.hpp"
#include "Field.hpp"
#include "Properties.hpp"
#include "SurfaceField.hpp"

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
            mesh.cellVolumes[i] * properties.rho() * properties.Cp() / dt;
        phiEqn.addToA(i, i, trans);
        phiEqn.addTob(i, trans * phi.oldField()[i]);
    }
}
inline void Div(SpaceMatrix& phiEqn,
                const Properties& properties,
                const Field<Vector>& U,
                const Field<double>& phi)
{
    const Mesh& mesh = phi.mesh;
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];
        double area = mesh.faceAreas[i];
        Vector areaDirect = mesh.faceNormals[i];
        Vector U_face = (U[o] + U[n]) / 2.0;

        double flux = properties.rho() * U_face.dotWith(areaDirect) * area;

        if (flux > 0) // 流动方向 o --> n   面变量等于o单元
        {
            phiEqn.addToA(o, o, flux); // o单元方程，面变量来自o，流出为正
            phiEqn.addToA(n, o, -flux); // n单元方程，面变量来自o，流入为负
        }
        else // 流动方向 n --> o   面变量等于n单元
        {
            phiEqn.addToA(o, n, flux); // o单元方程，面变量来自n，流入为正
            phiEqn.addToA(n, n, -flux); // n单元方程，面变量来自n，流出为负
        };
    }

    for (int pIdx = 0; pIdx < mesh.boundary.size(); ++pIdx)
    {
        const auto& patch = mesh.boundary[pIdx];
        const auto& Ubc = U.boundaryField[pIdx];
        const auto& phi_bc = phi.boundaryField[pIdx];
        for (int i = patch.firstFaceIdx; i < patch.firstFaceIdx + patch.nFaces;
             ++i)
        {
            int o = mesh.owner[i];
            double area = mesh.faceAreas[i];
            Vector areaDirect = mesh.faceNormals[i];

            Vector U_face = Ubc[i - patch.firstFaceIdx];
            double phi_face = phi_bc[i - patch.firstFaceIdx];

            double flux = properties.rho() * U_face.dotWith(areaDirect) * area;

            if (flux > 0) // 流出
            {
                phiEqn.addToA(o, o, flux);
            }
            else // 流入
            {
                phiEqn.addTob(o, -flux * phi_face);
            };
        }
    }
}

inline void Laplacian(SpaceMatrix& phiEqn,
                      const Properties& properties,
                      const Field<double> phi)
{ // 1. 处理内部面
    const Mesh& mesh = phi.mesh;
    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
        int o = mesh.owner[i];
        int n = mesh.neighbour[i];
        double area = mesh.faceAreas[i];
        double dist = (mesh.cellCentres[n] - mesh.cellCentres[o]).getMag();
        double coeff = properties.kappa() * area / dist;

        phiEqn.addToA(o, o, coeff);
        phiEqn.addToA(o, n, -coeff);
        phiEqn.addToA(n, n, coeff);
        phiEqn.addToA(n, o, -coeff);
    }

    // 2. 处理边界
    for (int pIdx = 0; pIdx < mesh.boundary.size(); ++pIdx)
    {
        const auto& patch = mesh.boundary[pIdx];
        const auto& bc = phi.boundaryList[pIdx];
        for (int i = patch.firstFaceIdx; i < patch.firstFaceIdx + patch.nFaces;
             ++i)
        {
            int o = mesh.owner[i];
            double area = mesh.faceAreas[i];
            double d = (mesh.faceCentres[i] - mesh.cellCentres[o]).getMag();

            double coeff_p = properties.kappa() * area * bc.fraction / d;
            double source =
                properties.kappa() * area * bc.fraction * bc.refValue / d +
                properties.kappa() * area * (1.0 - bc.fraction) * bc.refGrad;

            phiEqn.addToA(o, o, coeff_p);
            phiEqn.addTob(o, source);
        }
    }
}

}; // namespace fvm

namespace fvc
{
template <typename valType>
inline SurfaceField<valType> CellToSurface(const Field<valType>& cellField)
{
    const Mesh& mesh = cellField.mesh;
    SurfaceField<valType> tempField(mesh);

    const std::vector<Vector>& cellCentres = mesh.cellCentres;
    const std::vector<Vector>& faceCentres = mesh.faceCentres;
    const std::vector<Vector>& faceNormals = mesh.faceNormals;

    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
    }

    return tempField;
}

inline Field<Vector> Grad(Field<double> phi)
{
    const Mesh& mesh = phi.mesh;
    Field<Vector> tmpGrad(mesh, Vector());

    for (int i = 0; i < mesh.nInternalFace; ++i)
    {
    }
    return tmpGrad;
}
} // namespace fvc

#endif