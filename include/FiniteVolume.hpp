#pragma once

#include "Mesh.hpp"
#include "Solver.hpp"
#include "Field.hpp"
#include "Point3D.hpp"

enum class TimeScheme
{
    STEADY,
    IMPLICIT
};

class FiniteVolume
{
public:
    const Mesh& mesh;
    double k; // Thermal conductivity
    double rho;
    double Cp;
    double mu;
    bool Convective;
    bool Diffusive;
    bool Source;

    Field<double> T;
    Field<double> T_old;
    Field<double> SourceT;
    Field<Vector> U;

    FiniteVolume(const Mesh& mesh)
        : mesh(mesh),
          k(1.0),
          rho(1.0),
          Cp(1.0),
          mu(1.0),
          Convective(false),
          Diffusive(true),
          Source(true),
          T("T", mesh),
          T_old("Told", mesh),
          SourceT("sourceT", mesh, 0.0),
          U("U", mesh, Vector(1.0, 1.0, 0))
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

    void setProperties(double k, double rho, double Cp, double mu)
    {
        this->k = k;
        this->rho = rho;
        this->Cp = Cp;
        this->mu = mu;
    }

private:
    std::vector<std::vector<int>> cachedConnectivity;
    void prepareConnectivity();

    void assembleDiv(SpaceMatrix& A);
    void assembleLaplacian(SpaceMatrix& A);
    void assembleSource(SpaceMatrix& A);
    void assembleSpace(SpaceMatrix& A);
    void assembleTime(SpaceMatrix& A, double dt);

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
