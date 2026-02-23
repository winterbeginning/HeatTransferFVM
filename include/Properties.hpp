#ifndef _Properties_
#define _Properties_

#include "Field.hpp"
#include "Mesh.hpp"

class Properties
{
private:
    Field<double> kappa_; // Thermal conductivity
    Field<double> rho_;
    Field<double> Cp_;
    Field<double> mu_;

public:
    Properties(const Mesh& mesh)
        : kappa_(mesh, 1.0), rho_(mesh, 1.0), Cp_(mesh, 1.0), mu_(mesh, 1.0)
    {
        initBCs();
    }
    Properties(const Mesh& mesh, double k, double rho, double Cp, double mu)
        : kappa_(mesh, k), rho_(mesh, rho), Cp_(mesh, Cp), mu_(mesh, mu)
    {
        initBCs();
    }

    void initBCs()
    {
        setFieldBoundary(kappa_);
        setFieldBoundary(rho_);
        setFieldBoundary(Cp_);
        setFieldBoundary(mu_);
        correctAll();
    }
    void setKappa(double kappa)
    {
        this->kappa_.fill(kappa);
        this->kappa_.correctBoundaryField();
    }
    void setRho(double rho)
    {
        this->rho_.fill(rho);
        this->rho_.correctBoundaryField();
    }
    void setCp(double Cp)
    {
        this->Cp_.fill(Cp);
        this->Cp_.correctBoundaryField();
    }
    void setMu(double mu)
    {
        this->mu_.fill(mu);
        this->mu_.correctBoundaryField();
    }

    const Field<double>& kappa() const
    {
        return kappa_;
    }
    const Field<double>& rho() const
    {
        return rho_;
    }
    const Field<double>& Cp() const
    {
        return Cp_;
    }
    const Field<double>& mu() const
    {
        return mu_;
    }
    void setProperties(double kappa, double rho, double Cp, double mu)
    {
        this->kappa_.fill(kappa);
        this->rho_.fill(rho);
        this->Cp_.fill(Cp);
        this->mu_.fill(mu);
        correctAll();
    }
    void setFieldBoundary(Field<double>& prop)
    {
        for (auto& patch : prop.boundaryList)
        {
            patch.setBoundary(0.0, 0.0, 0.0);
        }
    }
    void correctAll()
    {
        this->kappa_.correctBoundaryField();
        this->rho_.correctBoundaryField();
        this->Cp_.correctBoundaryField();
        this->mu_.correctBoundaryField();
    }
};

#endif