#ifndef LD_INFGENS_HH
#define LD_INFGENS_HH

#include "LD_framework.hh"

struct rinfgen_workspace
{
  rinfgen_workspace(int perm_len_,int nvar_,int ord_len_): lamvec(new double[perm_len_]), vxu_wkspc(vxu_workspace(nvar_,ord_len_)) {}
  ~rinfgen_workspace() {delete [] lamvec;}
  double * const lamvec;
  vxu_workspace vxu_wkspc;
};

class r_xu_infgen: public LD_vector_field, private rinfgen_workspace
{

  double  * const xu,
          * const vx_vec;

  public:

    r_xu_infgen(function_space &fspc_,LD_vector_space **Vspaces_):
      LD_vector_field(fspc_,Vspaces_,fspc_.ndep,1), rinfgen_workspace(fspc_.perm_len,fspc_.nvar,fspc_.comp_ord_len()),
      xu(new double[fspc_.nvar]), vx_vec(new double[fspc_.perm_len]) {}
    r_xu_infgen(r_xu_infgen &rxu_): r_xu_infgen(rxu_.fspc,rxu_.Vspaces) {}
    ~r_xu_infgen() {delete [] xu; delete [] vx_vec;}

    virtual void dudx_eval(double x_,double *u_,double *dudx_)
    {
      xu[0] = x_;
      for (size_t i = 0; i < ndep; i++) xu[i+1] = u_[i] + (dudx_[i] = 0.0);
      fspc.lamvec_eval(xu,lamvec,vxu_wkspc);
      const int len_lam = fspc.perm_len,
                nvec_use = Vspaces[ispc]->nvec_use;
      double ** const Theta_mat = Vspaces[ispc]->Vmat,
                      Vx2 = 0.0;
      for (size_t itheta = 0; itheta < nvec_use; itheta++)
      {
        double &vx_itheta = vx_vec[itheta] = 0.0;
        for (size_t i = 0; i < len_lam; i++) vx_itheta += Theta_mat[itheta][i]*lamvec[i];
        for (size_t idep = 0,i_theta = len_lam; idep < ndep; idep++)
          for (size_t ilam = 0; ilam < len_lam; ilam++,i_theta++)
            dudx_[idep] += vx_itheta*Theta_mat[itheta][i_theta]*lamvec[ilam];
        Vx2 += vx_itheta*vx_itheta;
      }
      for (size_t i = 0; i < ndep; i++) dudx_[i] /= Vx2;
    }
};

class rn_infgen: public LD_prolonged_vfield, private rinfgen_workspace
{

  double  * const s_local,
          * const v_local,
          * const theta_local;

  public:

    rn_infgen(function_space_basis &fbse_,LD_vector_space **Vspaces_):
      LD_prolonged_vfield(fbse_,Vspaces_,fbse_.ndep,fbse_.eor), rinfgen_workspace(fspc.perm_len,fspc.nvar,fspc.comp_ord_len()),
      s_local(new double[fspc.ndim]), v_local(new double[fspc.ndim]), theta_local(new double[fspc.ndof_full]) {}
    rn_infgen(rn_infgen &rn_, function_space_basis &fbse_): rn_infgen(fbse_,rn_.Vspaces) {}
    ~rn_infgen() {delete [] s_local; delete [] v_local; delete [] theta_local;}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      s_local[0] = x_;
      for (size_t i = 0; i < ndep; i++) s_local[i+1] = u_[i];
      fspc.lamvec_eval(s_local,lamvec,vxu_wkspc);
      const int len_lam = fspc.perm_len,
                len_theta = fspc.ndof_full,
                nvec_use = Vspaces[ispc]->nvec_use,
                ndof_ODE = eor*ndep;
      double ** const Theta_mat = Vspaces[ispc]->Vmat,
                      Vx2 = 0.0;
      LD_linalg::fill_vec<double>(theta_local,len_theta,0.0);
      for (size_t itheta = 0; itheta < nvec_use; itheta++) // compute local parameter values
      {
        double  vx_itheta = 0.0;
        for (size_t i = 0; i < len_lam; i++) vx_itheta += Theta_mat[itheta][i]*lamvec[i];
        for (size_t i = 0; i < len_theta; i++) theta_local[i] += vx_itheta*Theta_mat[itheta][i];
        Vx2 += vx_itheta*vx_itheta;
      }
      for (size_t i = 0; i < len_theta; i++) theta_local[i] /= Vx2; // normalize local parameters
      for (size_t i = ndep; i < ndof_ODE; i++) s_local[i+1] = dudx_[i-ndep] = u_[i]; // load 1 : n-1 derivative terms
      for (size_t i = ndof_ODE+1; i < ndim; i++) s_local[i] = 0.0; // pad n'th order terms with zeroes
      fbse.v_eval(s_local,v_local,theta_local,eor-1);
      for (size_t i = ndof_ODE-ndep; i < ndof_ODE; i++) dudx_[i] = v_local[i+1];
    }
};

#endif
