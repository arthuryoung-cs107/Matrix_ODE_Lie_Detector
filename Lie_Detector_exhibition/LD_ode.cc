#include "LD_ode.hh"

ode_solspc_subset::ode_solspc_subset(ode_solspc_meta &meta_,int nobs_,bool palloc_,bool Jalloc_) :
  ode_solspc_element(meta_),
  local_data(true), nobs(nobs_),
  pts_mat(new double*[nobs]),
  sols(new ode_solution*[nobs]),
  dnp1xu_mat(palloc_?(new double*[nobs]):NULL),
  JFs_tns(Jalloc_?(new double**[nobs]):NULL)
  {for (int i = 0; i < nobs; i++) sols[i] = NULL;}
ode_solspc_subset::ode_solspc_subset(ode_solspc_meta &meta_,int nobs_,double **pts_mat_,ode_solution **sols_,double **dnp1xu_mat_,double ***JFs_tns_):
  ode_solspc_element(meta_),
  local_data(false), nobs(nobs_),
  pts_mat(pts_mat_), sols(sols_),
  dnp1xu_mat(dnp1xu_mat_), JFs_tns(JFs_tns_)
  {}
ode_solspc_subset::~ode_solspc_subset()
{
  if (local_data)
  {
    delete [] pts_mat;
    for (size_t i = 0; i < nobs; i++) if (sols[i] != NULL) delete sols[i];
    delete [] sols;
    if (dnp1xu_mat != NULL) delete [] dnp1xu_mat;
    if (JFs_tns != NULL) delete [] JFs_tns;
  }
}

solspc_data_chunk::solspc_data_chunk(ode_solspc_meta &meta_, int nobs_,bool palloc_,bool Jalloc_):
ode_solspc_subset(meta_,nobs_,
Tmatrix<double>(nobs_,meta_.ndim),
new ode_solution*[nobs_],
(palloc_)?(Tmatrix<double>(nobs_,meta_.ndep)):(NULL),
(Jalloc_)?(T3tensor<double>(nobs_,meta_.ndep,meta_.ndim)):(NULL)),
data_owner(true)
{initialize_solutions(true); initialize_additional_data();}
solspc_data_chunk::solspc_data_chunk(ode_solspc_meta &meta_,int nobs_,double **pts_mat_,ode_solution **sols_,double **dnp1xu_mat_,double ***JFs_tns_):
ode_solspc_subset(meta_,nobs_,pts_mat_,sols_,dnp1xu_mat_,JFs_tns_), data_owner(false)
{initialize_additional_data();}
solspc_data_chunk::~solspc_data_chunk()
{
  if ((data_owner)&&(!local_data))
  {
    free_Tmatrix<double>(pts_mat);
    for (size_t i = 0; i < nobs; i++) if (sols[i] != NULL) delete sols[i];
    delete [] sols;
    if (dnp1xu_mat != NULL) free_Tmatrix<double>(dnp1xu_mat);
    if (JFs_tns != NULL) free_T3tensor<double>(JFs_tns);
  }
}

ode_solcurve::ode_solcurve(int icrv_,ode_solspc_meta &meta_,int nobs_,bool palloc_,bool Jalloc_):
solspc_data_chunk(meta_,nobs_,palloc_,Jalloc_),
icrv(icrv),eps_vec(new double[nobs]) {}
ode_solcurve::ode_solcurve(int icrv_,ode_solspc_meta &meta_,int nobs_,double **pts_mat_,ode_solution **sols_,double **dnp1xu_mat_,double ***JFs_tns_):
solspc_data_chunk(meta_,nobs_,pts_mat_,sols_,dnp1xu_mat_,JFs_tns_),
icrv(icrv),eps_vec(new double[nobs]) {}
ode_solcurve::~ode_solcurve() {delete [] eps_vec;}
