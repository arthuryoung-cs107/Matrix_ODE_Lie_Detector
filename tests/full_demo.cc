#include "LD_matrices.hh"
#include "matrix_Lie_detector.hh"
#include "LD_ode_system.hh"
#include "LD_integrators.hh"
#include <sys/stat.h>

// specify data directory for writing binary files
const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";
const char addtl_prefix[] = "";

// specify subject ordinary differential equation for tests
Duffing_ode ode;
// VanDerPol_ode ode;
// Pendulum_ode ode;
// Bessel_ode ode;
// Riccati_ode ode;
// Brusselator_ode ode;

ode_solspc_meta meta0(ode.eor,ode.ndep);

// number of curves and uniform number of points per curve for dataset
// const int nc = 30, // number of curves
const int nc = 50, // number of curves
// const int nc = 75, // number of curves
// const int nc = 100, // number of curves
// const int nc = 150, // number of curves
// const int nc = 200, // number of curves
          np_min = 300, // min number of points for extrapolation experiment
          np = np_min; // points per curve
          // np = 350; // points per curve

// identifier for independent variable range of data
const int xrange = 0;
// const int xrange = 1;
// const int xrange = 2;

const bool  write_dnp1xu = true,
            write_JFs = true;

// level of noise applied to observational data. If <0, then unnoised
const int noise_level = -1;
// const int noise_level = 0;

LD_observations_set Sobs(meta0,nc,np,write_dnp1xu,write_JFs);

// specify runge kutta integrator for the generation of synthetic data
DoP853_settings integrator_settings; DoP853 ode_integrator(ode,integrator_settings);
// DoPri5_settings integrator_settings; DoPri5 ode_integrator(ode,integrator_settings);

// specify order of embedding function space
const int bor = 10;
// const int bor = 9;
// const int bor = 8;
// const int bor = 7;
// const int bor = 6;
// const int bor = 5;
// const int bor = 4;
// const int bor = 3;
// const int bor = 2;


orthopolynomial_space fspace0(meta0,bor); // class of embedding function space
orthopolynomial_basis fbasis0(fspace0); // prolongation workspace based on function space

// orthogonal polynomial family for function space configuration
// const char fam_name[] = "Legendre";
const char fam_name[] = "Chebyshev1";
// const char fam_name[] = "Chebyshev2";

// names of encoded data matrices
const char Lmat_name[] = "Lmat";
const char Omat_name[] = "Omat";
const char Rmat_name[] = "Rmat";
const char Pmat_name[] = "Pmat";
const char Qmat_name[] = "Qmat";
const char Gmat_name[] = "Gmat";

LD_L_matrix Lmat(fspace0,Sobs);

LD_R_matrix Amat(fspace0,Sobs); const char * const Amat_name = Rmat_name;
// LD_Q_matrix Amat(fspace0,Sobs); const char * const Amat_name = Qmat_name;

LD_O_matrix Bmat(fspace0,Sobs); const char * const Bmat_name = Omat_name;
// LD_R_matrix Bmat(fspace0,Sobs); const char * const Bmat_name = Rmat_name;
// LD_Q_matrix Bmat(fspace0,Sobs); const char * const Bmat_name = Qmat_name;

LD_matrix_svd_result Lmat_svd(nc,fspace0.perm_len);
LD_matrix_svd_result Amat_svd(nc,fspace0.ndof_full);
LD_matrix_svd_result Bmat_svd(nc,fspace0.ndof_full);
LD_alternate_svd_result AYmat_svd(nc,fspace0.ndof_full);
LD_alternate_svd_result BYmat_svd(nc,fspace0.ndof_full);

const bool  write_gen_obs_data = true,
            write_fspace_config = true,
            write_encoded_mats = true,
            write_decoded_mats = true,
            write_all_svds = true,
            write_recon_data = true,
            write_extnd_data = true;

r1space_infinitesimal_generator Ar1infgen0(fspace0,Amat_svd.VTtns),
                                AYr1infgen0(fspace0,AYmat_svd.VTtns_alt);
rnspace_infinitesimal_generator Brninfgen0(fbasis0,Bmat_svd.VTtns),
                                BYrninfgen0(fbasis0,BYmat_svd.VTtns_alt);

// runge kutta integrator for the reconstruction of observational data
DoP853_settings intgr_rec_settings; DoP853 intgr_rec(Ar1infgen0,intgr_rec_settings); // 8th order accurate, 7th order interpolation
// DoPri5_settings intgr_rec_settings; DoPri5 intgr_rec(Ar1infgen0,intgr_rec_settings); // 5th order accurate, 4th order interpolation

// buffers for writing file names via meta strings
char  eqn_name[50],
      intgen_name[50],
      noise_name[20],
      data_name[100],
      bse_name[75],
      intrec_name[50],
      recon_mat_name[50];

// character buffer for writing / reading data
char  name_buffer[250],
      name_dnp1xu_buffer[250],
      name_JFs_buffer[250];

bool  dnp1xu_empty = true,
      JFs_empty = true;

int generate_observational_data()
{
  mkdir(dir_name, S_IRWXU);

  if (np>np_min) sprintf(eqn_name,"%s_xrange%d_extrap",ode.name,xrange);
  else sprintf(eqn_name,"%s_xrange%d",ode.name,xrange);

  strcpy(intgen_name,ode_integrator.name);

  generated_ode_observations inputs_gen(ode,nc,np);
  inputs_gen.set_random_ICs(LD_rng(9365),ode.get_default_IC_range(),ode.get_default_IC_indep_range(xrange));  // shay's number
  // inputs_gen.set_random_ICs(LD_rng(8426),ode.get_default_IC_range(),ode.get_default_IC_indep_range(xrange));  // my number
  // inputs_gen.set_random_ICs(LD_rng(86427),ode.get_default_IC_range(),ode.get_default_IC_indep_range(xrange));  // another number
  inputs_gen.generate_solution_curves(ode_integrator);

  if (noise_level>=0)
  {
    // inputs_gen.apply_noise_pts();
    sprintf(noise_name,"noise%d",noise_level);
  }
  else sprintf(noise_name,"true");

  sprintf(data_name,"%s%s_%s_%sgen",addtl_prefix,eqn_name,noise_name,intgen_name);

  if (write_gen_obs_data)
  {
    sprintf(name_buffer,"%s/%s.%s",dir_name,data_name,dat_suff);
    inputs_gen.write_observed_solutions(name_buffer);
    if (write_dnp1xu)
    {
      inputs_gen.generate_dnp1xu();
      sprintf(name_buffer,"%s/%s_dnp1xu.%s",dir_name,data_name,dat_suff);
      inputs_gen.write_dnp1xu(name_buffer);
    }
    if (write_JFs)
    {
      inputs_gen.generate_JFs();
      sprintf(name_buffer,"%s/%s_JFs.%s",dir_name,data_name,dat_suff);
      inputs_gen.write_JFs(name_buffer);
    }
  }
  return 0;
}

int configure_function_space(bool check_fspaces_=false)
{
  sprintf(name_buffer,"%s/%s.%s",dir_name,data_name,dat_suff);
  Sobs.load_additional_inputs(ode_curve_observations(name_buffer),true);

  if (write_fspace_config)
  {
    fspace0.set_Legendre_coeffs(); Sobs.configure_centered_domain(fspace0);
    if (check_fspaces_) fspace0.debugging_description();
    sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s",dir_name,data_name,"Legendre",bor,dat_suff);
    fspace0.write_configuration_file(name_buffer);

    fspace0.set_Chebyshev1_coeffs(); Sobs.configure_0maxmag_0pi05_domain(fspace0);
    if (check_fspaces_) fspace0.debugging_description();
    sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s",dir_name,data_name,"Chebyshev1",bor,dat_suff);
    fspace0.write_configuration_file(name_buffer);

    // fspace0.set_Chebyshev2_coeffs(); Sobs.configure_0maxmag_0pi05_domain(fspace0);
    // fspace0.set_Chebyshev2_coeffs(); Sobs.configure_centered_domain(fspace0);
    fspace0.set_Chebyshev2_coeffs(); Sobs.configure_center_mass_domain(fspace0);
    if (check_fspaces_) fspace0.debugging_description();
    sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s",dir_name,data_name,"Chebyshev2",bor,dat_suff);
    fspace0.write_configuration_file(name_buffer);
  }
  sprintf(bse_name, "%s.%d",fam_name,bor);
  sprintf(name_buffer, "%s/%s_%s.domain_config.%s",dir_name,data_name,bse_name,dat_suff);
  fspace0.configure_self(name_buffer);
  fspace0.debugging_description();

  return 0;
}

template <class BSIS> int encode_data_matrices(BSIS **bases_)
{
  bases_[LD_threads::thread_id()]->debugging_description();

  sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Lmat_name,dat_suff);

  if (write_encoded_mats)
  {
    Lmat.populate_L_matrix<BSIS>(bases_); Lmat.write_matrix(name_buffer);

    LD_O_matrix Omat(fspace0,Sobs); Omat.populate_O_matrix<BSIS>(bases_);
    sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Omat_name,dat_suff);
    Omat.write_matrix(name_buffer);

    LD_R_matrix Rmat(fspace0,Sobs); Rmat.populate_R_matrix<BSIS>(bases_);
    sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Rmat_name,dat_suff);
    Rmat.write_matrix(name_buffer);

    if (write_dnp1xu)
    {
      if (dnp1xu_empty)
      {
        sprintf(name_buffer, "%s/%s.%s",dir_name,data_name,dat_suff);
        sprintf(name_dnp1xu_buffer, "%s/%s_dnp1xu.%s",dir_name,data_name,dat_suff);
        Sobs.load_additional_inputs(ode_curve_observations(name_buffer,name_dnp1xu_buffer),false);
        dnp1xu_empty = false;
      }

      LD_P_matrix Pmat(fspace0,Sobs); Pmat.populate_P_matrix<BSIS>(bases_);
      sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Pmat_name,dat_suff);
      Pmat.write_matrix(name_buffer);

      LD_Q_matrix Qmat(fspace0,Sobs); Qmat.populate_Q_matrix<BSIS>(bases_);
      sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Qmat_name,dat_suff);
      Qmat.write_matrix(name_buffer);

      // additional matrix / svd result
      LD_matrix OPmat(fspace0,Sobs,Omat.dim_cnstr + Pmat.dim_cnstr,fspace0.ndof_full);
      LD_matrix::concatenate_matrices(OPmat,Omat,Pmat);
      sprintf(name_buffer, "%s/%s_%s.OPmat.%s", dir_name,data_name,bse_name,dat_suff);
      OPmat.write_matrix(name_buffer);
      if (write_decoded_mats)
      {
        matrix_Lie_detector::compute_curve_svds(OPmat,Amat_svd,OPmat.nrow_curve(np_min));
        sprintf(name_buffer, "%s/%s_%s.OPmat_svd.%s", dir_name,data_name,bse_name,dat_suff);
        Amat_svd.write_svd_results(name_buffer);
      }
    }
    if (write_JFs)
    {
      if (JFs_empty)
      {
        sprintf(name_buffer, "%s/%s.%s",dir_name,data_name,dat_suff);
        sprintf(name_JFs_buffer, "%s/%s_JFs.%s",dir_name,data_name,dat_suff);
        Sobs.load_additional_inputs(ode_curve_observations(name_buffer,name_JFs_buffer),false);
        JFs_empty = false;
      }

      LD_G_matrix Gmat(fspace0,Sobs); Gmat.populate_G_matrix<BSIS>(bases_);
      sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Gmat_name,dat_suff);
      Gmat.write_matrix(name_buffer);

      // additional matrix / svd result
      LD_matrix OGmat(fspace0,Sobs,Omat.dim_cnstr + Gmat.dim_cnstr, fspace0.ndof_full);
      LD_matrix::concatenate_matrices(OGmat,Omat,Gmat);
      sprintf(name_buffer, "%s/%s_%s.OGmat.%s", dir_name,data_name,bse_name,dat_suff);
      OGmat.write_matrix(name_buffer);
      if (write_decoded_mats)
      {
        matrix_Lie_detector::compute_curve_svds(OGmat,Amat_svd,OGmat.nrow_curve(np_min));
        sprintf(name_buffer, "%s/%s_%s.OGmat_svd.%s", dir_name,data_name,bse_name,dat_suff);
        Amat_svd.write_svd_results(name_buffer);
      }
    }
  }
  else Lmat.read_matrix(name_buffer);

  sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Amat_name,dat_suff);
  Amat.read_matrix(name_buffer);
  sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Bmat_name,dat_suff);
  Bmat.read_matrix(name_buffer);

  return 0;
}

template <class MAT> void comp_crv_svd(MAT mat_, const char * const mat_name_)
{
  sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,mat_name_,dat_suff);
  mat_.read_matrix(name_buffer);
  matrix_Lie_detector::compute_curve_svds(mat_,Amat_svd,mat_.nrow_curve(np_min));
  sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,mat_name_,dat_suff);
  Amat_svd.write_svd_results(name_buffer);
}

int decode_data_matrices()
{
  sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,Lmat_name,dat_suff);
  if (write_decoded_mats)
  {
    matrix_Lie_detector::compute_curve_svds(Lmat,Lmat_svd,Lmat.nrow_curve(np_min));
    Lmat_svd.write_svd_results(name_buffer);

    AYmat_svd.compute_restricted_curve_svds(Amat,Lmat_svd,Amat.nrow_curve(np_min));
    sprintf(name_buffer, "%s/%s_%s.%sYL_svd.%s", dir_name,data_name,bse_name,Amat_name,dat_suff);
    AYmat_svd.write_svd_results(name_buffer);

    BYmat_svd.compute_restricted_curve_svds(Bmat,Lmat_svd,Bmat.nrow_curve(np_min));
    sprintf(name_buffer, "%s/%s_%s.%sYL_svd.%s", dir_name,data_name,bse_name,Bmat_name,dat_suff);
    BYmat_svd.write_svd_results(name_buffer);

    if (write_all_svds)
    {
      comp_crv_svd<LD_O_matrix>(LD_O_matrix(fspace0,Sobs),Omat_name);
      comp_crv_svd<LD_R_matrix>(LD_R_matrix(fspace0,Sobs),Rmat_name);
      if (write_dnp1xu)
      {
        if (dnp1xu_empty)
        {
          sprintf(name_buffer, "%s/%s.%s",dir_name,data_name,dat_suff);
          sprintf(name_dnp1xu_buffer, "%s/%s_dnp1xu.%s",dir_name,data_name,dat_suff);
          Sobs.load_additional_inputs(ode_curve_observations(name_buffer,name_dnp1xu_buffer),false);
          dnp1xu_empty = false;
        }

        comp_crv_svd<LD_P_matrix>(LD_P_matrix(fspace0,Sobs),Pmat_name);
        comp_crv_svd<LD_Q_matrix>(LD_Q_matrix(fspace0,Sobs),Qmat_name);
      }
      if (write_JFs)
      {
        if (JFs_empty)
        {
          sprintf(name_buffer, "%s/%s.%s",dir_name,data_name,dat_suff);
          sprintf(name_JFs_buffer, "%s/%s_JFs.%s",dir_name,data_name,dat_suff);
          Sobs.load_additional_inputs(ode_curve_observations(name_buffer,name_JFs_buffer),false);
          JFs_empty = false;
        }
        comp_crv_svd<LD_G_matrix>(LD_G_matrix(fspace0,Sobs),Gmat_name);
      }

      sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,Amat_name,dat_suff);
      Amat_svd.read_svd_results(name_buffer);
      sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,Bmat_name,dat_suff);
      Bmat_svd.read_svd_results(name_buffer);
    }
    else
    {
      matrix_Lie_detector::compute_curve_svds(Amat,Amat_svd,Amat.nrow_curve(np_min));
      sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,Amat_name,dat_suff);
      Amat_svd.write_svd_results(name_buffer);

      matrix_Lie_detector::compute_curve_svds(Bmat,Bmat_svd,Bmat.nrow_curve(np_min));
      sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,Bmat_name,dat_suff);
      Bmat_svd.write_svd_results(name_buffer);
    }
  }
  else
  {
    Lmat_svd.read_svd_results(name_buffer);

    sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,Amat_name,dat_suff);
    Amat_svd.read_svd_results(name_buffer);
    sprintf(name_buffer, "%s/%s_%s.%sYL_svd.%s", dir_name,data_name,bse_name,Amat_name,dat_suff);
    AYmat_svd.read_svd_results(name_buffer); AYmat_svd.set_VTtns_alt_Y_VAY(Lmat_svd);

    sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,Bmat_name,dat_suff);
    Bmat_svd.read_svd_results(name_buffer);
    sprintf(name_buffer, "%s/%s_%s.%sYL_svd.%s", dir_name,data_name,bse_name,Bmat_name,dat_suff);
    BYmat_svd.read_svd_results(name_buffer); BYmat_svd.set_VTtns_alt_Y_VAY(Lmat_svd);
  }
  Lmat_svd.print_details("Lmat_svd");
  Amat_svd.print_details("Amat_svd");
  AYmat_svd.print_details("AYmat_svd");
  Bmat_svd.print_details("Bmat_svd");
  BYmat_svd.print_details("BYmat_svd");

  return 0;
}

template <class BSIS, class INFGN, class INTGR> void infgen_reconstruct(BSIS **bases_,INFGN &infgn_,INTGR &intgr_)
{
  generated_ode_observations inputs_recon(infgn_,Sobs.ncrvs_tot,Sobs.npts_per_crv);
  sprintf(name_buffer, "%s/%s_%s.%s.%srec.%s",dir_name,data_name,bse_name,recon_mat_name,intrec_name,dat_suff);
  inputs_recon.set_solcurve_ICs(Sobs.curves);
  inputs_recon.parallel_generate_solution_curves<BSIS,INFGN,INTGR>(bases_,infgn_,intgr_);
  inputs_recon.write_observed_solutions(name_buffer);
}

template <class INFGN, class INTGR> void infgen_reconstruct(INFGN &infgn_,INTGR &intgr_)
{
  generated_ode_observations inputs_recon(infgn_,Sobs.ncrvs_tot,Sobs.npts_per_crv);
  sprintf(name_buffer, "%s/%s_%s.%s.%srec.%s",dir_name,data_name,bse_name,recon_mat_name,intrec_name,dat_suff);
  inputs_recon.set_solcurve_ICs(Sobs.curves);
  inputs_recon.parallel_generate_solution_curves<INFGN,INTGR>(infgn_,intgr_);
  inputs_recon.write_observed_solutions(name_buffer);
}

template <class BSIS, class INFGN> void infgen_extend(BSIS **bases_,INFGN &infgn_)
{
  generated_ode_observations inputs_recon(infgn_,Sobs.ncrvs_tot,Sobs.npts_per_crv);
  sprintf(name_buffer, "%s/%s_%s.%s.%srec.%s", dir_name,data_name,bse_name,recon_mat_name,intrec_name,dat_suff);
  inputs_recon.read_basic_observations(name_buffer);

  ode_curve_observations inputs_extnd(meta0.eor,meta0.ndep,inputs_recon.nobs);
  INFGN::init_extended_observations(inputs_extnd,inputs_recon);
  matrix_Lie_detector::extend_ode_observations<INFGN,BSIS>(inputs_extnd,infgn_,bases_);
  sprintf(name_buffer, "%s/%s_%s.%s.%sext.%s", dir_name,data_name,bse_name,recon_mat_name,intrec_name,dat_suff);
  inputs_extnd.write_observed_solutions(name_buffer);
}

template <class BSIS, class INFGN, class INTGR> void infgen_reconstruct_extend(BSIS **bases_,INFGN &infgn_,INTGR &intgr_)
{
  generated_ode_observations inputs_recon(infgn_,Sobs.ncrvs_tot,Sobs.npts_per_crv);
  if (write_recon_data)
  {
    sprintf(name_buffer, "%s/%s_%s.%s.%srec.%s",dir_name,data_name,bse_name,recon_mat_name,intrec_name,dat_suff);
    inputs_recon.set_solcurve_ICs(Sobs.curves);
    inputs_recon.parallel_generate_solution_curves<INFGN,INTGR>(infgn_,intgr_);
    inputs_recon.write_observed_solutions(name_buffer);
  }
  else
  {
    sprintf(name_buffer, "%s/%s_%s.%s.%srec.%s", dir_name,data_name,bse_name,recon_mat_name,intrec_name,dat_suff);
    inputs_recon.read_basic_observations(name_buffer);
  }
  if (write_extnd_data)
  {
    ode_curve_observations inputs_extnd(meta0.eor,meta0.ndep,inputs_recon.nobs);
    INFGN::init_extended_observations(inputs_extnd,inputs_recon);
    matrix_Lie_detector::extend_ode_observations<INFGN,BSIS>(inputs_extnd,infgn_,bases_);
    sprintf(name_buffer, "%s/%s_%s.%s.%sext.%s", dir_name,data_name,bse_name,recon_mat_name,intrec_name,dat_suff);
    inputs_extnd.write_observed_solutions(name_buffer);
  }
}

template <class BSIS, class INTGR> int reconstruct_data(BSIS ** bases_, INTGR &intgr_)
{
  sprintf(intrec_name, "%sr1",intgr_.name);
  strcpy(recon_mat_name,Amat_name); Ar1infgen0.init_svd_default(Amat_svd);
  infgen_reconstruct_extend<BSIS,r1space_infinitesimal_generator,INTGR>(bases_,Ar1infgen0,intgr_);
  sprintf(recon_mat_name, "%sYL",Amat_name); AYr1infgen0.init_svd_default(AYmat_svd);
  infgen_reconstruct_extend<BSIS,r1space_infinitesimal_generator,INTGR>(bases_,AYr1infgen0,intgr_);

  sprintf(intrec_name, "%srn",intgr_.name);
  strcpy(recon_mat_name,Bmat_name); Brninfgen0.init_svd_default(Bmat_svd);
  infgen_reconstruct<BSIS,rnspace_infinitesimal_generator,INTGR>(bases_,Brninfgen0,intgr_);
  sprintf(recon_mat_name, "%sYL",Bmat_name); BYrninfgen0.init_svd_default(BYmat_svd);
  infgen_reconstruct<BSIS,rnspace_infinitesimal_generator,INTGR>(bases_,BYrninfgen0,intgr_);

  return 0;
}

int main()
{
  orthopolynomial_basis **bases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0);

  int gen_check = generate_observational_data(),
      cnf_check = configure_function_space(),
      enc_check = encode_data_matrices<orthopolynomial_basis>(bases0),
      dec_check = decode_data_matrices(),
      rec_check = reconstruct_data(bases0,intgr_rec);

  free_evaluation_bases<orthopolynomial_basis>(bases0);
  return 0;
}
