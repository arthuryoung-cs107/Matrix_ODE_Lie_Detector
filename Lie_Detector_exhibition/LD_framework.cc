#include "LD_framework.hh"

ode_curve_observations::ode_curve_observations(int eor_, int ndep_, int nobs_):
  eor(eor_), ndep(ndep_), nobs(nobs_),
  pts_in(new double[(1+ndep_*(eor_+1))*nobs_]) {}
ode_curve_observations::ode_curve_observations(int eor_, int ndep_, int ncrv_, int nobs_):
  eor(eor_), ndep(ndep_), ncrv(ncrv_), nobs(nobs_),
  npts_per_crv(new int[ncrv_]), pts_in(new double[(1+ndep_*(eor_+1))*nobs_])
  {for (int i = 0, np = nobs_/ncrv_; i < ncrv_; i++) npts_per_crv[i] = np;}
ode_curve_observations::ode_curve_observations(int eor_, int ndep_, int ncrv_, int *npts_per_crv_):
  eor(eor_), ndep(ndep_), ncrv(ncrv_), nobs(LD_linalg::sum_vec<int>(npts_per_crv_,ncrv_)),
  npts_per_crv(new int[ncrv_]), pts_in(new double[(1+ndep_*(eor_+1))*nobs])
  {LD_linalg::copy_vec<int>(npts_per_crv,npts_per_crv_,ncrv_);}
ode_curve_observations::ode_curve_observations(const char name_[])
{
  FILE * file_in = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(ode_meta,sizeof(int),2,file_in);
  LD_io::fread_SAFE(obs_meta,sizeof(int),2,file_in);
  npts_per_crv = new int[ncrv];
  int ndim = 1+(ndep*(eor+1));
  pts_in = new double[ndim*nobs];
  LD_io::fread_SAFE(npts_per_crv,sizeof(int),ncrv,file_in);
  LD_io::fread_SAFE(pts_in,sizeof(double),ndim*nobs,file_in);
  LD_io::fclose_SAFE(file_in);
  printf("(ode_curve_observations::ode_curve_observations) read %s\n",name_);
}
ode_curve_observations::~ode_curve_observations()
{
  if (npts_per_crv != NULL) delete [] npts_per_crv;
  if (pts_in != NULL) delete [] pts_in;
  if (JFs_in != NULL) delete [] JFs_in;
  if (dnp1xu_in != NULL) delete [] dnp1xu_in;
}

void ode_curve_observations::read_basic_observations(const char name_[], bool force_overwrite_)
{
  int ode_meta_in[2],
      &eor_in = ode_meta_in[0],
      &ndep_in = ode_meta_in[1],
      obs_meta_in[2],
      &ncrv_in = obs_meta_in[0],
      &nobs_in = obs_meta_in[1];

  FILE * file_in = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(ode_meta_in,sizeof(int),2,file_in);
  LD_io::fread_SAFE(obs_meta_in,sizeof(int),2,file_in);
  if ((eor_in==eor)&&(ndep_in==ndep)&&(ncrv_in==ncrv)&&(nobs_in==nobs)) // everything matches, can read directly
  {
    int ndim = 1+(ndep*(eor+1));
    if (npts_per_crv==NULL) npts_per_crv = new int[ncrv];
    LD_io::fread_SAFE(npts_per_crv,sizeof(int),ncrv,file_in);
    if (pts_in==NULL) pts_in = new double[ndim*nobs];
    LD_io::fread_SAFE(pts_in,sizeof(double),ndim*nobs,file_in);
  }
  else if ((pts_in!=NULL)&&( (comp_ndim(eor,ndep)*nobs) == (comp_ndim(eor_in,ndep_in)*nobs_in) )) // if the pts buffers are same size
  {
    if (npts_per_crv==NULL) npts_per_crv = new int[ncrv=ncrv_in];
    if (ncrv_in<=ncrv) LD_io::fread_SAFE(npts_per_crv,sizeof(int),ncrv=ncrv_in,file_in);
    LD_io::fread_SAFE(pts_in,sizeof(double),comp_ndim(eor=eor_in,ndep=ndep_in)*(nobs=nobs_in),file_in); // reconfigure system dims
  }
  else if (force_overwrite_)
  {
    eor=eor_in; ndep=ndep_in; ncrv=ncrv_in; nobs=nobs_in;
    int ndim_in = 1 + ndep*(eor+1);
    if ((npts_per_crv==NULL)&&(pts_in==NULL)) // can do a clean initialization
    {
      npts_per_crv = new int[ncrv];
      pts_in = new double[ndim_in*nobs];
    }
    LD_io::fread_SAFE(npts_per_crv,sizeof(int),ncrv,file_in);
    LD_io::fread_SAFE(pts_in,sizeof(double),ndim_in*nobs,file_in);
  }
  LD_io::fclose_SAFE(file_in);
  printf("(ode_curve_observations::read_basic_observations) read %s\n",name_);
}

void ode_curve_observations::read_additional_observations(const char name_addtl_[])
{
  const char fcn_name[] = "ode_curve_observations::read_additional_observations";
  const int hlen_max = 4,
            ndim = (1+(ndep*(eor+1))),
            chunk_len_max = nobs*ndep*ndim;
  int hlen,
      header[hlen_max],
      chunk_len_read;
  FILE * file_in = LD_io::fopen_SAFE(name_addtl_,"r");
  LD_io::fread_SAFE(&hlen,sizeof(int),1,file_in);
  LD_io::fread_SAFE(header,sizeof(int),(hlen<=hlen_max)?(hlen):(hlen_max),file_in);
  double * const chunk_in = new double[chunk_len_read = (header[0]<=chunk_len_max)?(header[0]):(chunk_len_max)];
  LD_io::fread_SAFE(chunk_in,sizeof(double),chunk_len_read,file_in);
  switch (hlen)
  {
    case 3:
      if (dnp1xu_in==NULL)
      {
        dnp1xu_in = chunk_in;
        if (!( (header[0]==(ndep*nobs)) && (header[1]==nobs) && (header[2]==ndep) ))
          printf("(%s) WARNING - inconsistent dnp1xu data dimensions in %s \n",fcn_name,name_addtl_);
      }
      else
      {
        printf("(%s) ERROR - attemping to overwrite dnp1xu_in with %s \n",fcn_name,name_addtl_);
        delete [] chunk_in;
      }
      break;
    case 4:
      if (JFs_in==NULL)
      {
        JFs_in = chunk_in;
        if (!( (header[0]==(chunk_len_max)) && (header[1]==nobs) && (header[2]==ndep) && (header[3]==ndim)))
          printf("(%s) WARNING - inconsistent JFs data dimensions in %s \n",fcn_name,name_addtl_);
      }
      else
      {
        printf("(%s) ERROR - attemping to overwrite JFs_in with %s \n",fcn_name,name_addtl_);
        delete [] chunk_in;
      }
      break;
    default:
      printf("(%s) FAILED to read %s (hlen=%d).\n",fcn_name,name_addtl_,hlen);
      delete [] chunk_in;
  }
  LD_io::fclose_SAFE(file_in);
  printf("(%s) read %s\n",fcn_name,name_addtl_);
}

void ode_curve_observations::write_observed_solutions(const char name_[])
{
  FILE * file_out = LD_io::fopen_SAFE(name_,"wb");
  fwrite(ode_meta,sizeof(int),2,file_out);
  fwrite(obs_meta,sizeof(int),2,file_out);
  fwrite(npts_per_crv,sizeof(int),ncrv,file_out);
  fwrite(pts_in,sizeof(double),(1 + ndep*(eor+1))*nobs,file_out);
  LD_io::fclose_SAFE(file_out);
  printf("(ode_curve_observations::write_observed_solutions) wrote %s\n",name_);
}

generated_ode_observations::generated_ode_observations(ode_system &ode_, int nc_, int np_):
  ode_curve_observations(ode_,nc_,nc_*np_), ode(ode_), pts_IC(new double*[nc_]), xrange_mat(Tmatrix<double>(nc_,2))
  {for (size_t i = 0, iic = 0, idelic = np_*ndim; i < nc_; i++, iic+=idelic) pts_IC[i] = pts_in + iic;}
generated_ode_observations::generated_ode_observations(ode_system &ode_, int nc_, int *npts_per_crv_):
  ode_curve_observations(ode_,nc_,npts_per_crv_), ode(ode_), pts_IC(new double*[nc_]), xrange_mat(Tmatrix<double>(nc_,2))
  {for (size_t ic = 0, iic = 0; ic < nc_; iic+=ndim*(npts_per_crv_[ic++])) pts_IC[ic] = pts_in + iic;}

void generated_ode_observations::set_Gaussian_random_ICs(LD_rng rng_, const double *ICr_, const double *xrange_)
{
  for (size_t i = 0; i < ncrv; i++)
  {
    xrange_mat[i][0] = xrange_[0]; xrange_mat[i][1] = xrange_[1];
    for (size_t i_dim = 1, iicr = 0; i_dim <= ndof_ODE; i_dim++, iicr+=2)
      pts_IC[i][i_dim] = rng_.rand_gau( 0.5*(ICr_[iicr+1]+ICr_[iicr]),
                                        (ICr_[iicr+1]-ICr_[iicr])/10 );
  }
}
void generated_ode_observations::set_random_ICs(LD_rng rng_, const double *ICr_, const double *xrange_)
{
  for (size_t i = 0; i < ncrv; i++)
  {
    xrange_mat[i][0] = xrange_[0]; xrange_mat[i][1] = xrange_[1];
    for (size_t i_dim = 1, iicr = 0; i_dim <= ndof_ODE; i_dim++, iicr+=2)
      pts_IC[i][i_dim] = rng_.rand_uni(ICr_[iicr],ICr_[iicr+1]);
  }
}
void generated_ode_observations::generate_solution_curves(ode_integrator &integrator_)
{
  double  * const u_state = integrator_.get_u_state(),
          ** const integr_wkspc = Tmatrix<double>(LD_linalg::max_val<int>(npts_per_crv,ncrv),ndof_ODE),
          t0 = LD_threads::tic();
  for (size_t icrv = 0; icrv < ncrv; icrv++)
  {
    for (size_t i_dof = 0; i_dof < ndof_ODE; i_dof++) u_state[i_dof] = pts_IC[icrv][i_dof+1];
    integrator_.init_curve_integration((xrange_mat[icrv][1]-xrange_mat[icrv][0])/((double) (npts_per_crv[icrv]-1)),icrv);
    integrator_.set_and_solve_time(xrange_mat[icrv][0],xrange_mat[icrv][1],npts_per_crv[icrv],integr_wkspc);
    integrator_.unpack_time_sol(xrange_mat[icrv][0],npts_per_crv[icrv],integr_wkspc,pts_IC[icrv]);
  }
  free_Tmatrix<double>(integr_wkspc);
  printf("(generated_ode_observations::generate_solution_curves) exponentiated %d integral curves (%d net snaps, %d degrees of freedom) in %.4f seconds (SINGLE thread)\n",
    ncrv, nobs, integrator_.ndof_ODE, LD_threads::toc(t0));
}

void generated_ode_observations::generate_JFs()
{
  const int len_JFs_i = ndep*ndim;
  double  * const Jac_chunk = JFs_in = new double[len_JFs_i*nobs],
          ** const JFs_mat_i = new double*[ndep];
  for (size_t iobs = 0; iobs < nobs; iobs++)
  {
    double  * const Jac_mat_i_start = Jac_chunk + (iobs*len_JFs_i),
            * const pts_i = pts_in + (iobs*ndim);
    for (size_t idep = 0; idep < ndep; idep++) JFs_mat_i[idep] = Jac_mat_i_start + (idep*ndim);
    ode.JacF_eval(pts_i[0],pts_i+1,JFs_mat_i);
  }
  delete [] JFs_mat_i;
}
void generated_ode_observations::write_JFs(const char name_[])
{
  int hlen = 4,
      len_JFs = ndep*ndim*nobs,
      header[] = {hlen,len_JFs,nobs,ndep,ndim};
  FILE * file_out = LD_io::fopen_SAFE(name_,"wb");
  fwrite(header,sizeof(int),hlen+1,file_out);
  fwrite(JFs_in,sizeof(double),len_JFs,file_out);
  LD_io::fclose_SAFE(file_out);
  printf("(generated_ode_observations::write_JFs) wrote %s\n",name_);
}

void generated_ode_observations::generate_dnp1xu()
{
  dnp1xu_in = new double[ndep*nobs];
  for (size_t iobs = 0; iobs < nobs; iobs++)
  {
    double  * const pts_i = pts_in + (iobs*ndim);
    ode.dnp1xu_eval(pts_i[0],pts_i+1,dnp1xu_in + (iobs*ndep));
  }
}
void generated_ode_observations::write_dnp1xu(const char name_[])
{
  int hlen = 3,
      len_dnp1xu = ndep*nobs,
      header[] = {hlen,len_dnp1xu,nobs,ndep};
  FILE * file_out = LD_io::fopen_SAFE(name_,"wb");
  fwrite(header,sizeof(int),hlen+1,file_out);
  fwrite(dnp1xu_in,sizeof(double),len_dnp1xu,file_out);
  LD_io::fclose_SAFE(file_out);
  printf("(generated_ode_observations::write_dnp1xu) wrote %s\n",name_);
}

LD_observations_set::LD_observations_set(ode_solspc_meta &meta_, ode_curve_observations input_):
solspc_data_chunk(meta_,input_.nobs,input_.palloc(),input_.Jalloc()),
ncrvs_tot(input_.ncrv), npts_per_crv(new int[ncrvs_tot]), pts_tns(new double**[ncrvs_tot]), curves(new ode_solcurve*[ncrvs_tot]),
dnp1xu_tns((input_.palloc())?( new double**[ncrvs_tot] ):NULL),
JFs_crv((input_.Jalloc())?( new double***[ncrvs_tot] ):NULL)
{
  if (input_.npts_per_crv != NULL) LD_linalg::copy_vec<int>(npts_per_crv,input_.npts_per_crv,ncrvs_tot);
  if (input_.pts_in != NULL) LD_linalg::copy_vec<double>(pts_chunk,input_.pts_in,nobs*ndim);
  if (input_.dnp1xu_in != NULL) LD_linalg::copy_vec<double>(dnp1xu_chunk,input_.dnp1xu_in,nobs*ndep);
  if (input_.JFs_in != NULL) LD_linalg::copy_vec<double>(JFs_chunk,input_.JFs_in,nobs*ndep*ndim);
  for (size_t icrv = 0, ipts=0; icrv < ncrvs_tot; ipts+=npts_per_crv[icrv++])
    curves[icrv] = new ode_solcurve(icrv,meta,npts_per_crv[icrv], pts_tns[icrv] = pts_mat+ipts,
                                                                  sols+ipts,
                                                                  (dnp1xu_tns!=NULL)?(dnp1xu_tns[icrv] = dnp1xu_mat+ipts):NULL,
                                                                  (JFs_crv!=NULL)?(JFs_crv[icrv] = JFs_tns+ipts):NULL);
}
LD_observations_set::LD_observations_set(ode_solspc_meta &meta_, int ncrv_, int npts_, bool palloc_, bool Jalloc_):
solspc_data_chunk(meta_,ncrv_*npts_,palloc_,Jalloc_),
ncrvs_tot(ncrv_), npts_per_crv(new int[ncrvs_tot]), pts_tns(new double**[ncrvs_tot]), curves(new ode_solcurve*[ncrvs_tot]),
dnp1xu_tns((palloc_)?( new double**[ncrvs_tot] ):NULL),
JFs_crv((Jalloc_)?( new double***[ncrvs_tot] ):NULL)
{
  for (size_t icrv = 0, ipts=0; icrv < ncrvs_tot; ipts+=npts_per_crv[icrv++])
    curves[icrv] = new ode_solcurve(icrv,meta,  npts_per_crv[icrv] = npts_,
                                                pts_tns[icrv] = pts_mat+ipts,
                                                sols+ipts,
                                                (dnp1xu_tns!=NULL)?(dnp1xu_tns[icrv] = dnp1xu_mat+ipts):NULL,
                                                (JFs_crv!=NULL)?(JFs_crv[icrv] = JFs_tns+ipts):NULL);
}
LD_observations_set::~LD_observations_set()
{
  for (size_t i = 0; i < ncrvs_tot; i++) delete curves[i];
  delete [] curves;
  delete [] npts_per_crv;
  delete [] pts_tns;
  if (dnp1xu_tns != NULL) delete [] dnp1xu_tns;
  if (JFs_crv != NULL) delete [] JFs_crv;
}
void LD_observations_set::load_additional_inputs(ode_curve_observations input_, bool overwrite_basics_)
{
  if (overwrite_basics_)
  {
    if (input_.npts_per_crv != NULL) LD_linalg::copy_vec<int>(npts_per_crv,input_.npts_per_crv,ncrvs_tot);
    if (input_.pts_in != NULL) LD_linalg::copy_vec<double>(pts_chunk,input_.pts_in,nobs*ndim);
  }
  if (input_.dnp1xu_in != NULL)
  {
    if (dnp1xu_tns==NULL)
    {
      alloc_dnp1xu_safe(input_.dnp1xu_in); // allocate contiguous dnp1xu data, set pointers for solutions
      dnp1xu_tns = new double**[ncrvs_tot]; // allocate dnp1xu pointers for each curve
      for (size_t icrv = 0, ipts=0; icrv < ncrvs_tot; ipts+=npts_per_crv[icrv++])
      {
        curves[icrv]->dnp1xu_mat = dnp1xu_tns[icrv] = dnp1xu_mat+ipts; // set each curve's pointer to dnp1xu data
        curves[icrv]->initialize_additional_data(); // set each curve's pointer to contiguous dnp1xu data
      }
    }
    else LD_linalg::copy_vec<double>(dnp1xu_chunk,input_.dnp1xu_in,nobs*ndep); // just copy data, pointers already initialized
  }
  if (input_.JFs_in != NULL)
  {
    if (JFs_crv==NULL)
    {
      alloc_JFs_safe(input_.JFs_in); // allocate contiguous JFs data, set pointers for solutions
      JFs_crv = new double***[ncrvs_tot]; // allocate JFs pointers for each curve
      for (size_t icrv = 0, ipts=0; icrv < ncrvs_tot; ipts+=npts_per_crv[icrv++])
      {
        curves[icrv]->JFs_tns = JFs_crv[icrv] = JFs_tns+ipts; // set each curve's pointer to JFs data
        curves[icrv]->initialize_additional_data(); // set each curve's pointer to contiguous JFs data
      }
    }
    else LD_linalg::copy_vec<double>(JFs_chunk,input_.JFs_in,nobs*ndep*ndim); // just copy data, pointers already initialized
  }
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
    for (int j = 0; j < nobs; j++)
    {
      double * const pts_j = sols[j]->pts;
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
    for (int j = 0; j < nobs; j++)
    {
      double * const pts_j = sols[j]->pts;
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

void LD_observations_set::get_solspace_center_mass(double *scm_g_)
{
  const int ind_1 = (nobs%2)?(nobs/2):((nobs/2) - 1),
            ind_2 = (nobs%2)?(nobs/2):(nobs/2);
  #pragma omp parallel
  {
    int ind_work;
    double scm_t[nobs];
    #pragma omp for
    for (size_t idim = 0; idim < ndim; idim++)
    {
      for (size_t iobs = 0, iind = idim; iobs < nobs; iobs++, iind+=ndim) scm_t[iobs] = pts_chunk[iind];
      LD_linalg::sort_vec_inc<double>(scm_t,nobs,ind_work);
      scm_g_[idim] = 0.5*(scm_t[ind_1] + scm_t[ind_1]);
    }
  }
}

LD_matrix_file::LD_matrix_file(const char name_[])
{
  FILE * file_in = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(&hlen_in,sizeof(int),1,file_in);
  if (hlen_in==hlen_check)
  {
    LD_io::fread_SAFE(header_in,sizeof(int),hlen_check,file_in);
    Amat_in = Tmatrix<double>(nrows_in,ncols_in);
    LD_io::fread_SAFE(Amat_in[0],sizeof(double),nrows_in*ncols_in,file_in);
    printf("(LD_matrix_file::LD_matrix_file) read %s \n", name_);
  }
  else printf("(LD_matrix_file::LD_matrix_file) ERROR: hlen_in != hlen_check (%d vs. %d)\n", hlen_in, hlen_check);
  LD_io::fclose_SAFE(file_in);
}

LD_matrix::LD_matrix(function_space &fspc_, LD_observations_set &Sset_, int dim_cnstr_, int net_cols_):
function_space_element(fspc_), LD_experiment(Sset_), dim_cnstr(dim_cnstr_), net_cols(net_cols_),
Attns(new double***[ncrvs_tot]), Atns(T3tensor<double>(nobs_full,dim_cnstr,net_cols))
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
      header[] = {hlen,net_rows,net_cols};
  fwrite(header, sizeof(int), hlen+1, file);
  fwrite(Avec, sizeof(double), net_rows*net_cols, file);
  fclose(file);
  printf("(LD_matrix::write_matrix) wrote %s\n",name_);
}
void LD_matrix::read_matrix(const char name_[])
{
  const int hlen_check = 2,
            Alen_check = net_rows*net_cols;
  int hlen_in,
      Alen_in,
      header_in[hlen_check],
      &net_rows_in = header_in[0],
      &net_cols_in = header_in[1];
  FILE * file_in = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(&hlen_in,sizeof(int),1,file_in);
  if (hlen_in==hlen_check)
  {
    LD_io::fread_SAFE(header_in,sizeof(int),hlen_in,file_in);
    Alen_in=net_rows_in*net_cols_in;
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
