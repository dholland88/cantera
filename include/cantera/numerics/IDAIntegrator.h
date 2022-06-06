/**
 *  @file IDAIntegrator.h
 *  Header file for class IDAIntegrator
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_IDAIntegrator_H
#define CT_IDAIntegrator_H

#include "cantera/numerics/Integrator.h"
#include "cantera/base/ctexceptions.h"
#include "sundials/sundials_nvector.h"

namespace Cantera
{


/**
 * Wrapper for Sundials IDA solver
 * @see FuncEval.h. Classes that use IDAIntegrator:
 * FlowReactor
 */
class IDAIntegrator : public Integrator
{
public:
    /**
     *  Constructor. Default settings: dense Jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     */
    IDAIntegrator();
    virtual ~IDAIntegrator();
    virtual void setTolerances(double reltol, size_t n, double* abstol);
    virtual void setTolerances(double reltol, double abstol);
    virtual void setSensitivityTolerances(double reltol, double abstol);
    virtual void setProblemType(int probtype);
    virtual void initialize(double t0, FuncEval& func);
    virtual void reinitialize(double t0, FuncEval& func);
    virtual void integrate(double tout);
    virtual double step(double tout);
    virtual double& solution(size_t k);
    virtual double* solution();
    virtual int nEquations() const {
        return static_cast<int>(m_neq);
    }
    virtual int nEvals() const;
    virtual void setMaxOrder(int n) {
        m_maxord = n;
    }
    virtual void setMaxStepSize(double hmax);
    virtual void setMaxSteps(int nmax);
    virtual int maxSteps();
    virtual void setMaxErrTestFails(int n);
    virtual void setBandwidth(int N_Upper, int N_Lower) {
        m_mupper = N_Upper;
        m_mlower = N_Lower;
    }
    virtual int nSensParams() {
        return static_cast<int>(m_np);
    }
    virtual double sensitivity(size_t k, size_t p);

    //! Returns a string listing the weighted error estimates associated
    //! with each solution component.
    //! This information can be used to identify which variables are
    //! responsible for integrator failures or unexpected small timesteps.
    virtual std::string getErrorInfo(int N);

    //! Error message information provide by IDAS
    std::string m_error_message;

    virtual void setMaxNonlinIterations(int n);
    virtual void setMaxNonlinConvFailures(int n);
    virtual void inclAlgebraicInErrorTest(bool yesno);
    virtual void setMethod(MethodType t);
    virtual void setIterator(IterType t);

protected:
    protected:
    //! Applies user-specified options to the underlying IDAS solver. Called
    //! during integrator initialization or reinitialization.
    void applyOptions();

private:
    void sensInit(double t0, FuncEval& func);
    size_t m_neq;
    void* m_ida_mem; //!< Pointer to the IDA memory for the problem
    void* m_linsol; //!< Sundials linear solver object
    void* m_linsol_matrix; //!< matrix used by Sundials
    void * m_ctx;           //!< contex object used by Sundials
    FuncEval* m_func;
    double m_t0;
    double m_time; //!< The current integrator time
    N_Vector m_y, m_ydot, m_abstol;
    int m_type;
    int m_itol;
    int m_maxord;
    double m_reltol;
    double m_abstols;
    double m_reltolsens, m_abstolsens;
    size_t m_nabs;
    double m_hmax, m_hmin;
    int m_maxsteps;
    int m_maxErrTestFails;
    N_Vector* m_yS;
    N_Vector* m_ySdot;
    size_t m_np;
    int m_mupper, m_mlower;
    N_Vector m_constraints;

    //! Indicates whether the sensitivities stored in m_yS have been updated
    //! for at the current integrator time.
    bool m_sens_ok;

    //! Maximum number of nonlinear solver iterations at one solution
    /*!
     *  If zero, this is the default of 4.
     */
    int m_maxNonlinIters;

    //! Maximum number of nonlinear convergence failures
    int m_maxNonlinConvFails;

    //! If true, the algebraic variables don't contribute to error tolerances
    bool m_setSuppressAlg;

    //! Initial IDA stepsize
    double m_init_step;
};

}

#endif