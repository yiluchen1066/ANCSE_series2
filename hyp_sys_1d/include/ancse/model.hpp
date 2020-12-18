#ifndef HYPSYS1D_MODEL_HPP
#define HYPSYS1D_MODEL_HPP

#include <cmath>
#include <memory>

#include <Eigen/Dense>
#include <ancse/config.hpp>

/// Interface for implementing different models,
/// eg. Euler equations, Shallow-water equations
///
/// Add more functions to this interface if needed.
class Model {
  public:
    virtual ~Model() = default;

    virtual Eigen::VectorXd flux(const Eigen::VectorXd &u) const = 0;
    virtual Eigen::VectorXd eigenvalues(const Eigen::VectorXd &u) const = 0;
    virtual Eigen::MatrixXd eigenvectors(const Eigen::VectorXd &u) const = 0;
    virtual double max_eigenvalue(const Eigen::VectorXd &u) const = 0;

    virtual Eigen::VectorXd cons_to_prim(const Eigen::VectorXd &u) const = 0;
    virtual Eigen::VectorXd prim_to_cons(const Eigen::VectorXd &u) const = 0;

    virtual Eigen::VectorXd roe_avg(const Eigen::VectorXd &uL,
                                    const Eigen::VectorXd &uR) const = 0;

    virtual int get_nvars() const = 0;
    virtual std::string get_name() const = 0;
};

class Burgers : public Model {
  public:

    Eigen::VectorXd flux(const Eigen::VectorXd &u) const override
    {
        Eigen::VectorXd f(n_vars);
        f(0) = 0.5*u(0)*u(0);

        return f;
    }

    Eigen::VectorXd eigenvalues(const Eigen::VectorXd &u) const override
    {
        Eigen::VectorXd eigvals(n_vars);
        eigvals(0) = u(0);

        return eigvals;
    }

    Eigen::MatrixXd eigenvectors(const Eigen::VectorXd &) const override
    {
        Eigen::MatrixXd eigvecs(n_vars, n_vars);
        eigvecs (0,0) = 1;

        return eigvecs;
    }

    double max_eigenvalue(const Eigen::VectorXd &u) const override {
        return (eigenvalues(u).cwiseAbs()).maxCoeff();
    }

    Eigen::VectorXd cons_to_prim(const Eigen::VectorXd &u) const override {
        return u;
    }

    Eigen::VectorXd prim_to_cons(const Eigen::VectorXd &u) const override {
        return u;
    }

    Eigen::VectorXd roe_avg(const Eigen::VectorXd &,
                            const Eigen::VectorXd &) const override
    {return Eigen::VectorXd::Zero(n_vars);}

    int get_nvars() const override
    {
        return n_vars;
    }

    std::string get_name() const override
    {
        return "burgers";
    }

  private:
    int n_vars = 1;
};

/// Euler equations
class Euler : public Model {
  public:

    Eigen::VectorXd flux(const Eigen::VectorXd &u) const override {
        Eigen::VectorXd f(n_vars);
        f(0) = u(1);
        f(1)=u(1)*u(1)/u(0)+u(0);
        f(2)=u(2)*u(1)/u(0)+u(1);

        return f;

    }
    
    Eigen::VectorXd eigenvalues(const Eigen::VectorXd &u) const override {
        double v = u(1)/u(0);
        double p = (u(2)-0.5*u(0) * u(1)/u(0)*u(1)/u(0))*(gamma-1);

        double c = sqrt(gamma*p/u(0));

        Eigen::VectorXd eigenvals=eigenvalues(v,c);
        return eigenvals;

    }

    Eigen::MatrixXd eigenvectors(const Eigen::VectorXd &u) const override{
        double v= u(1)/u(0);
        double p = u(2)-0.5*u(0)*u(1)/u(0)*u(1)/u(0)*(gamma-1);
        double H = (u(2)+p)/u(0);

        Eigen::MatrixXd eigvecs = eigenvectors(v,H);
        return eigvecs;
    };
    
    double max_eigenvalue(const Eigen::VectorXd &u) const override {
        return (eigenvalues(u).cwiseAbs()).maxCoeff();

    }

    Eigen::VectorXd cons_to_prim(const Eigen::VectorXd &u_cons) const override;

    Eigen::VectorXd prim_to_cons(const Eigen::VectorXd &u_prim) const override;

    Eigen::VectorXd roe_avg(const Eigen::VectorXd &,
                            const Eigen::VectorXd &) const override;

    ///  ANCSE_COMMENT Add more functions if needed.
    //----------------ModelEulerBegin----------------
    inline std::tuple<double, double, double>
    primitive(const Eigen::VectorXd &u_cons) const
    {
        double rho, v, p;

        rho = u_cons(0);
        v = u_cons(1)/u_cons(0);
        p = (u_cons(2)-0.5*u_cons(0) * u_cons(1)/u_cons(0)*u_cons(1)/u_cons(0) )*(gamma-1);

        return std::make_tuple(rho, v, p);

    }

    inline std::tuple<double, double, double>
    consertive(const Eigen::VectorXd &u_prim) const{
        double rho,rhov, E;
        rho = u_prim(0);
        rhov = u_prim(0)*u_prim(1);
        E = u_prim(2)/(gamma-1)+0.5*u_prim(0)*u_prim(1)*u_prim(1);

        return std::make_tuple(rho,rhov,E);
    }

    inline double pressure(double rho, double v, double E) const
    {
        return (gamma-1)*(E - 0.5*rho*v*v);
    }

    inline double sound_speed(double rho, double p) const
    {
        return sqrt(gamma*p/rho);
    }

    inline double enthalpy(double rho, double E, double p) const
    {
        return (E+p)/rho;
    }

    inline double energy(double rho, double v, double p) const
    {
        return (p/(gamma-1) + 0.5*rho*v*v);
    }


    inline Eigen::VectorXd eigenvalues(double v, double c) const
    {

        Eigen::VectorXd eigvals(n_vars);
        eigvals(0) = v-c;
        eigvals(1) = v;
        eigvals(2) = v+c;
        ///  ANCSE_COMMENT Compute eigenvalues
        return eigvals;
    }

    inline Eigen::MatrixXd eigenvectors(double v, double H) const
    {
        Eigen::MatrixXd eigvecs(n_vars, n_vars);
        double c = sqrt(gamma * (H - 0.5*v*v - 1/(gamma-1)));

        eigvecs (0,0) = 1;
        eigvecs(0,1) = 1;
        eigvecs(0,2)= 1;
        eigvecs(1,0) = v-c;
        eigvecs(1,1) = v;
        eigvecs(1,2)= v+c;
        eigvecs(2,0)= H-v*c;
        eigvecs(2,1) = v*v;
        eigvecs(2,2) = H + v*c;
        ///  ANCSE_COMMENT Compute eigenvectors

        return eigvecs;
    }

    //----------------ModelEulerEnd----------------


    void set_gamma(double gamma_)
    {
        gamma = gamma_;
    }

    double get_gamma() const
    {
        return gamma;
    }

    int get_nvars() const override
    {
        return n_vars;
    }

    std::string get_name() const override
    {
        return "euler";
    }

  private:
    int n_vars = 3;
    double gamma = 5./3.;
};


std::shared_ptr<Model> make_model (const nlohmann::json &config);

#endif // HYPSYS1D_MODEL_HPP
