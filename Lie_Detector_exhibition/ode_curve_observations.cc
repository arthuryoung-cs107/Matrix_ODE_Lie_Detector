#include "ode_curve_observations.hh"

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
