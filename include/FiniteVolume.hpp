#ifndef _FiniteVolume_
#define _FiniteVolume_

#include "Mesh.hpp"
#include "Solver.hpp"
#include "Field.hpp"
#include "Point3D.hpp"
#include "Properties.hpp"

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

    Field<double> T;
    Field<double> T_old;
    Field<double> SourceT;
    Field<Vector> U;

    // 对矩阵实例进行复用，避免每一时间步或迭代都重新申请内存
    SpaceMatrix TEqn;

    FiniteVolume(const Mesh& mesh)
        : mesh(mesh),
          properties(1.0, 1.0, 1.0, 1.0),
          Convective(false),
          Diffusive(true),
          Source(true),
          T(mesh),
          T_old(mesh),
          SourceT(mesh, 0.0),
          U(mesh, Vector(1.0, 1.0, 0)),
          TEqn(mesh.numCells)
    {
    }

    void
    solve(TimeScheme ts, Solver& solver, int maxSteps = 20, double dt = 0.1);

    void writeToTec(const string& filePath = "out.dat", double time = 0.0);

    void setSolveOption(bool Convective, bool Diffusive, bool Source)
    {
        this->Convective = Convective;
        this->Diffusive = Diffusive;
        this->Source = Source;
    }

private:
    std::vector<std::vector<int>> cachedConnectivity;
    void prepareConnectivity();

    void assembleSource(SpaceMatrix& A);
    void assembleMatrix(SpaceMatrix& A, TimeScheme ts, double dt = 1.0);

    void writeField(ofstream& outfile, int nElements, Field<double> field)
    {
        for (int i = 0; i < nElements; ++i)
            outfile << field[i] << ((i + 1) % 10 == 0 ? "\n" : " ");
        outfile << "\n";
    };

    void writeField(ofstream& outfile, int nElements, Field<Vector> field)
    {
        for (int i = 0; i < nElements; ++i)
            outfile << field[i].x << ((i + 1) % 10 == 0 ? "\n" : " ");
        outfile << "\n";
        for (int i = 0; i < nElements; ++i)
            outfile << field[i].y << ((i + 1) % 10 == 0 ? "\n" : " ");
        outfile << "\n";
        for (int i = 0; i < nElements; ++i)
            outfile << field[i].z << ((i + 1) % 10 == 0 ? "\n" : " ");
        outfile << "\n";
    };
};

#endif