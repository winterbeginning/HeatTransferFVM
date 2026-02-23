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
    }
    Properties(const Mesh& mesh, double k, double rho, double Cp, double mu)
        : kappa_(mesh, k), rho_(mesh, rho), Cp_(mesh, Cp), mu_(mesh, mu)
    {
    }
    void setKappa(double kappa)
    {
        this->kappa_.fill(kappa);
    }
    void setRho(double rho)
    {
        this->rho_.fill(rho);
    }
    void setCp(double Cp)
    {
        this->Cp_.fill(Cp);
    }
    void setMu(double mu)
    {
        this->mu_.fill(mu);
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
    }
    void correct()
    {
    }
};

#endif