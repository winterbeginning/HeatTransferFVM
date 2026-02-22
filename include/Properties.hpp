#ifndef _Properties_
#define _Properties_

class Properties
{
private:
    double kappa_; // Thermal conductivity
    double rho_;
    double Cp_;
    double mu_;

public:
    Properties() : kappa_(1.0), rho_(1.0), Cp_(1.0), mu_(1.0)
    {
    }
    Properties(double k, double rho, double Cp, double mu)
        : kappa_(k), rho_(rho), Cp_(Cp), mu_(mu)
    {
    }
    void setKappa(double kappa)
    {
        this->kappa_ = kappa;
    }
    void setRho(double rho)
    {
        this->rho_ = rho;
    }
    void setCp(double Cp)
    {
        this->Cp_ = Cp;
    }
    void setMu(double mu)
    {
        this->mu_ = mu;
    }

    const double& kappa() const
    {
        return kappa_;
    }
    const double& rho() const
    {
        return rho_;
    }
    const double& Cp() const
    {
        return Cp_;
    }
    const double& mu() const
    {
        return mu_;
    }
    void setProperties(double kappa, double rho, double Cp, double mu)
    {
        this->kappa_ = kappa;
        this->rho_ = rho;
        this->Cp_ = Cp;
        this->mu_ = mu;
    }
};

#endif