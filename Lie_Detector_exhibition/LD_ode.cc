#include "LD_ode.hh"
#include <cstdio>
#include <cstring>

ode_system::ode_system(int eor_, int ndep_, const char name_[]): ode_solspc_meta(eor_,ndep_), name(new char[strlen(name_)+1])
{strcpy(name,name_);}
ode_system::~ode_system() {delete [] name;}

ode_solspc_subset::ode_solspc_subset(ode_solspc_meta &meta_, int nobs_):
  ode_solspc_element(meta_), nobs(nobs_),
  pts_mat(new double*[nobs]), sols(new ode_solution*[nobs]) {}
ode_solspc_subset::~ode_solspc_subset() {delete [] pts_mat; delete [] sols;}

ode_solcurve::ode_solcurve(ode_solspc_meta &meta_, int nobs_, double *pts_chunk_, int icrv_):
ode_solspc_subset(meta_,nobs_), icrv(icrv_), pts_chunk(pts_chunk_)
{
  for (size_t i = 0, i_dim = 0; i < nobs; i++, i_dim+=ndim)
    sols[i] = new ode_solution(meta_,pts_mat[i] = pts_chunk + i_dim);
}
ode_solcurve::~ode_solcurve()
  {for (size_t i = 0; i < nobs; i++) delete sols[i];}

void ode_solution::print_sol()
{
  for (size_t idim = 0; idim < ndim; idim++) printf("%.2e ", pts[idim]);
  printf("\n");
}
