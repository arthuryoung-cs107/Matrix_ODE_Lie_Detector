#include "LD_framework.hh"
#include "LD_aux.hh"
#include "LD_io.hh"

ode_curve_observations::ode_curve_observations(int eor_, int ndep_, int ncrv_, int nobs_):
eor(eor_), ndep(ndep_), ncrv(ncrv_), nobs(nobs_),
npts_per_crv(new int[ncrv_]), pts_in(new double[(1+ndep_*(eor_+1))*nobs_])
{for (int i = 0, np = nobs_/ncrv_; i < ncrv_; i++) npts_per_crv[i] = np;}
ode_curve_observations::~ode_curve_observations()
{
  if (npts_per_crv != NULL) delete [] npts_per_crv;
  if (pts_in != NULL) delete [] pts_in;
  if (JFs_in != NULL) delete [] JFs_in;
  if (dnp1xu_in != NULL) delete [] dnp1xu_in;
}

generated_ode_observations::generated_ode_observations(ode_system &ode_, int nc_, int np_):
ode_curve_observations(ode_,nc_,nc_*np_), ode(ode_), npts(np_),
pts_IC(new double*[nc_])
{
  for (size_t i = 0, iic = 0, idelic = npts*ndim; i < nc_; i++, iic+=idelic)
    pts_IC[i] = pts_in + iic;
}

generated_ode_observations::~generated_ode_observations()
{
  delete [] pts_IC;
}

void generated_ode_observations::set_random_ICs(LD_rng rng_, const double *ICR_)
{
  for (size_t i = 0; i < ncrv; i++)
    for (size_t i_dim = 1, iicr = 0; i_dim <= ndof_ODE; i_dim++, iicr+=2)
      pts_IC[i][i_dim] = rng_.rand_uni(ICR_[iicr],ICR_[iicr+1]);
}

void generated_ode_observations::generate_solution_curves(ode_integrator &integrator_, const double * indep_range_)
{
  integrator_.del_t = (indep_range_[1]-indep_range_[0])/((double)(npts-1));
  double  * const wvec = integrator_.get_wvec(),
          ** const integr_wkspc = Tmatrix<double>(npts,ndof_ODE);
  for (size_t icrv = 0; icrv < ncrv; icrv++)
  {
    for (size_t i_dof = 0; i_dof < ndof_ODE; i_dof++)
      wvec[i_dof] = pts_IC[icrv][i_dof+1];
    integrator_.init_curve_integration(icrv);
    integrator_.set_and_solve_time(indep_range_[0],indep_range_[1],npts,integr_wkspc);
    integrator_.unpack_time_sol(indep_range_[0],npts,integr_wkspc,pts_IC[icrv]);
  }
  free_Tmatrix<double>(integr_wkspc);
}

void generated_ode_observations::write_solution_curves(const char name_[])
{
  FILE * file_out = LD_io::fopen_SAFE(name_,"wb");
  fwrite(ode_meta,sizeof(int),2,file_out);
  fwrite(obs_meta,sizeof(int),2,file_out);
  fwrite(npts_per_crv,sizeof(int),ncrv,file_out);
  fwrite(pts_in,sizeof(double),ndim*nobs,file_out);
  LD_io::fclose_SAFE(file_out);
  printf("(generated_ode_observations::write_solution_curves) wrote %s\n",name_);
}

input_ode_observations::input_ode_observations(const char name_[]): ode_curve_observations(), name(LD_io::duplicate_string(name_))
{
  FILE * file_in = LD_io::fopen_SAFE(name,"r");
  LD_io::fread_SAFE(ode_meta,sizeof(int),2,file_in);
  LD_io::fread_SAFE(obs_meta,sizeof(int),2,file_in);
  npts_per_crv = new int[ncrv];
  int ndim = 1+(ndep*(eor+1));
  pts_in = new double[ndim*nobs];
  LD_io::fread_SAFE(npts_per_crv,sizeof(int),ncrv,file_in);
  LD_io::fread_SAFE(pts_in,sizeof(double),ndim*nobs,file_in);
  LD_io::fclose_SAFE(file_in);
}
input_ode_observations::~input_ode_observations() {delete [] name;}

void input_ode_observations::print_details()
{
  printf("(input_ode_observations::print_details) name: %s\n  eor = %d, ndep = %d, ncrv = %d, nobs = %d\n", name,eor,ndep,ncrv,nobs);
}

// LD_observations_set::LD_observations_set(ode_solspc_meta &meta_, input_ode_observations &input_):
LD_observations_set::LD_observations_set(ode_solspc_meta &meta_, input_ode_observations input_):
ode_solspc_element(meta_),
ncrvs_tot(input_.ncrv), nobs_full(input_.nobs),
// npts_per_crv(input_.npts_per_crv), pts_chunk_full(input_.pts_in),
npts_per_crv(new int[ncrvs_tot]), pts_chunk_full(new double[nobs_full*ndim]),
pts_mat_full(new double*[nobs_full]), pts_tns_full(new double**[ncrvs_tot]),
curves(new ode_solcurve*[ncrvs_tot]), sols_full(new ode_solution*[nobs_full]),
indep_range(new double[2])
{
  LD_linalg::copy_x(input_.npts_per_crv,npts_per_crv,ncrvs_tot);
  LD_linalg::copy_x(input_.pts_in,pts_chunk_full,nobs_full*ndim);
  for (size_t icrv = 0, ipts=0, idim=0; icrv < ncrvs_tot; icrv++)
  {
    pts_tns_full[icrv] = pts_mat_full+ipts;
    curves[icrv] = new ode_solcurve(meta_,npts_per_crv[icrv],pts_chunk_full+idim,icrv);
    for (size_t ipts_crv = 0; ipts_crv < npts_per_crv[icrv]; ipts_crv++, ipts++, idim+=ndim)
    {
      pts_mat_full[ipts] = pts_chunk_full + idim;
      sols_full[ipts] = curves[icrv]->sols[ipts_crv];
    }
  }
}

LD_observations_set::~LD_observations_set()
{
  for (size_t i = 0; i < ncrvs_tot; i++) delete curves[i];
  delete [] curves; delete [] sols_full;

  delete [] pts_tns_full;
  delete [] pts_mat_full;

  delete [] npts_per_crv;
  delete [] pts_chunk_full;

  delete [] indep_range;
}

void LD_observations_set::get_solspace_val_extrema(double **sve_g_)
{
  bool first2finish=true;
  #pragma omp parallel
  {
    double sve_t[2][ndim];
    for (int l = 0; l < ndim; l++)
    {
      sve_t[0][l] = DBL_MAX;
      sve_t[1][l] = -1.0*DBL_MAX;
    }
    #pragma omp for nowait
    for (int j = 0; j < nobs_full; j++)
    {
      double * const pts_j = sols_full[j]->pts;
      for (int l = 0; l < ndim; l++)
      {
        if (sve_t[0][l]>pts_j[l]) sve_t[0][l]=pts_j[l];
        if (sve_t[1][l]<pts_j[l]) sve_t[1][l]=pts_j[l];
      }
    }
    #pragma omp critical
    {
      if (first2finish)
      {
        for (size_t i = 0; i < ndim; i++)
        {
          sve_g_[0][i] = sve_t[0][i];
          sve_g_[1][i] = sve_t[1][i];
        }
        first2finish = false;
      }
      else
        for (size_t i = 0; i < ndim; i++)
        {
          if (sve_g_[0][i]>sve_t[0][i]) sve_g_[0][i]=sve_t[0][i];
          if (sve_g_[1][i]<sve_t[1][i]) sve_g_[1][i]=sve_t[1][i];
        }
    }
  }
}

void LD_observations_set::get_solspace_mag_extrema(double **sme_g_)
{
  bool first2finish=true;
  #pragma omp parallel
  {
    double sme_t[2][ndim];
    for (int l = 0; l < ndim; l++)
    {
      sme_t[0][l] = DBL_MAX;
      sme_t[1][l] = 0.0;
    }
    #pragma omp for nowait
    for (int j = 0; j < nobs_full; j++)
    {
      double * const pts_j = sols_full[j]->pts;
      for (int l = 0; l < ndim; l++)
      {
        double mpjl = fabs(pts_j[l]);
        if (sme_t[0][l]>mpjl) sme_t[0][l]=mpjl;
        if (sme_t[1][l]<mpjl) sme_t[1][l]=mpjl;
      }
    }
    #pragma omp critical
    {
      if (first2finish)
      {
        for (size_t i = 0; i < ndim; i++)
        {
          sme_g_[0][i] = sme_t[0][i];
          sme_g_[1][i] = sme_t[1][i];
        }
        first2finish = false;
      }
      else
        for (size_t i = 0; i < ndim; i++)
        {
          if (sme_g_[0][i]>sme_t[0][i]) sme_g_[0][i]=sme_t[0][i];
          if (sme_g_[1][i]<sme_t[1][i]) sme_g_[1][i]=sme_t[1][i];
        }
    }
  }
}

LD_matrix::LD_matrix(function_space &fspc_, LD_observations_set &Sset_, int dim_cnstr_):
function_space_element(fspc_), LD_experiment(Sset_), dim_cnstr(dim_cnstr_),
Attns(new double***[ncrvs_tot]), Atns(T3tensor<double>(nobs_full,dim_cnstr,ndof_full))
{
  for (size_t icrv = 0, idim_crv = 0; icrv < ncrvs_tot; idim_crv+=npts_per_crv[icrv++])
    Attns[icrv] = Atns + idim_crv;
}
LD_matrix::~LD_matrix()
{
  free_T3tensor<double>(Atns);
  delete [] Attns;
}
void LD_matrix::write_matrix(const char name_[])
{
  FILE * file = LD_io::fopen_SAFE(name_,"wb");
  int hlen = 2,
      header[] = {hlen,net_rows,ndof_full};
  fwrite(header, sizeof(int), hlen+1, file);
  fwrite(Avec, sizeof(double), net_rows*ndof_full, file);
  fclose(file);
  printf("(LD_matrix::write_matrix) wrote %s\n",name_);
}
void LD_matrix::read_matrix(const char name_[])
{
  const int hlen_check = 2,
            Alen_check = net_rows*ndof_full;
  int hlen_in,
      Alen_in,
      header_in[hlen_check],
      &net_rows_in = header_in[0],
      &ndof_full_in = header_in[1];
  FILE * file_in = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(&hlen_in,sizeof(int),1,file_in);
  if (hlen_in==hlen_check)
  {
    LD_io::fread_SAFE(header_in,sizeof(int),hlen_in,file_in);
    Alen_in=net_rows_in*ndof_full_in;
  }
  else
  {
    printf("(LD_matrix::read_matrix) ERROR: hlen_in != hlen_check (%d vs. %d)\n", hlen_in, hlen_check);
    LD_io::fclose_SAFE(file_in);
    return;
  }
  if (Alen_in<=Alen_check)
  {
    if (Alen_in<Alen_check)
      printf("(LD_matrix::read_matrix) WARNING: Alen_in < Alen_check (%d vs %d). ", Alen_in, Alen_check);
    LD_io::fread_SAFE(Avec,sizeof(double),Alen_in,file_in);
    LD_io::fclose_SAFE(file_in);
    printf("(LD_matrix::read_matrix) read %s\n",name_);
  }
  else
  {
    printf("(LD_matrix::read_matrix) ERROR: Alen_in > Alen_check (%d vs %d) ", Alen_in, Alen_check);
    LD_io::fclose_SAFE(file_in);
  }
}
