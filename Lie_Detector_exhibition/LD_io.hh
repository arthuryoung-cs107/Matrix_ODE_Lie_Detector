#ifndef LD_IO_HH
#define LD_IO_HH

#include <cstring>
#include <cstdio>
#include <cstdlib>

struct LD_io
{
  LD_io() {}
  ~LD_io() {}

  static char * duplicate_string(const char in_[])
    {char * out = new char[strlen(in_)+1]; return strcpy(out,in_);}

  static FILE * fopen_SAFE(const char name_[], const char tags_[])
  {
    FILE * file_out = fopen(name_,tags_);
    if (file_out==NULL)
    {
      printf("(LD_io::fopen_SAFE) fopen error with %s (%s attempted). Exiting\n", name_, tags_);
      exit(1);
    }
    else return file_out;
  }
  static void fclose_SAFE(FILE *file_) {fclose(file_);}
  static void fread_SAFE(void *ptr_, size_t size_, size_t m_, FILE *file_)
  {
    if (fread(ptr_,size_,m_,file_)!=m_)
    {
      printf("(LD_io::fread_SAFE): failed to read %zu bytes (size = %zu, count = %zu). Exiting\n", size_*m_, size_, m_);
      exit(1);
    }
  }
  static void fseek_SAFE(FILE * file_, long int offset_, int start_)
  {
    if (fseek(file_,offset_,start_)!=0)
    {
      printf("(LD_io::fseek_SAFE): failed to shift file position by %ld bytes from %d. Exiting\n", offset_,start_);
      exit(1);
    }
  }

  template <typename T> static void write_Tmat(const char name_[], T **mat_, int nrows_, int ncols_)
  {
    const int hlen_c = 3;
    int header[hlen_c+1],
        &hlen = header[0] = hlen_c,
        &nrows = header[1] = nrows_,
        &ncols = header[2] = ncols_,
        &sze = header[3] = (int) sizeof(T);

    FILE * file_out = LD_io::fopen_SAFE(name_,"wb");
    fwrite(header,sizeof(int),hlen_c+1,file_out);
    for (size_t i = 0; i < nrows_; i++)
      fwrite(mat_[i],sze,ncols_,file_out);
    LD_io::fclose_SAFE(file_out);
    printf("(LD_io::write_Tmat): wrote %s\n", name_);
  }

};

class LD_name_buffer
{
  int check,
      len_buf_cmp;

  inline size_t max_size_t(size_t sz1_,size_t sz2_) {return (sz1_>sz2_)?(sz1_):(sz2_);}

  inline char * check_bad_encode(const char fcn_name_[])
  {
    if (check<0) printf("(LD_name_buffer::%s) ERROR - bad encoding for %s \n",fcn_name_,name);
    return name;
  }

  inline bool not_enough_space(int check_, bool alloc_more_=true, bool verbose_=false)
  {
    if (verbose_) printf("(LD_name_buffer::not_enough_space) WARNING - buffer overflow for %s \n", name);
    if ( (len_buf_cmp = check_ + 1) > len_buf )
    {
      if (alloc_more_)
      {
        delete [] name;
        name = new char[len_buf = (size_t)(len_buf_cmp)];
      }
      return true;
    }
    else return false;
  }

  inline void check_output_names() {check_dir_name(); check_dat_suf();}
  inline void check_dir_name() {if (dir_name==NULL) load_dir_name();} // initialize to default
  inline void check_dat_suf() {if (dat_suf==NULL) load_dat_suf();} // initialize to default

  public:

    LD_name_buffer(size_t len_buf_=10): len_buf(len_buf_), name(new char[len_buf_]) {}

    LD_name_buffer(const char name_[], size_t len_buf_=10):
      LD_name_buffer(max_size_t(strlen(name_)+1,len_buf_)) {load_name(name_);}

    LD_name_buffer(const char dir_name_[], const char dat_suf_[], size_t len_buf_=10):
      LD_name_buffer(len_buf_) {load_dir_name(dir_name_); load_dat_suf(dat_suf_);}

    LD_name_buffer(const char name_[], const char dir_name_[], const char dat_suf_[], size_t len_buf_=10):
     LD_name_buffer(name_,len_buf_) {load_dir_name(dir_name_); load_dat_suf(dat_suf_);}

    ~LD_name_buffer()
    {
      delete [] name;
      if (dir_name != NULL) delete [] dir_name;
      if (dat_suf != NULL) delete [] dat_suf;
    }

    size_t  len_buf;
    char  * name,
          * dir_name = NULL,
          * dat_suf = NULL;

    inline char * name_file(const char name_[])
    {
      check_output_names();
      return (not_enough_space(check = snprintf(name,len_buf,"%s/%s.%s",dir_name,name_,dat_suf)))
        ?(name_file(name_)) // try again
        :(check_bad_encode("name_file"));
    }
    inline char * name_file(LD_name_buffer &in_) {return name_file(in_.name);}
    inline char * name_file(const char name_[], const char suffix_[])
    {
      check_output_names();
      return (not_enough_space(check = snprintf(name,len_buf,"%s/%s%s.%s",dir_name,name_,suffix_,dat_suf)))
        ?(name_file(name_,suffix_)) // try again
        :(check_bad_encode("name_file"));
    }
    inline char * name_file(LD_name_buffer &in0_,const char in1_[]) {return name_file(in0_.name,in1_);}
    inline char * name_file(const char name1_[], const char name2_[], const char name3_[])
    {
      check_output_names();
      return (not_enough_space(check = snprintf(name,len_buf,"%s/%s%s%s.%s",dir_name,name1_,name2_,name3_,dat_suf)))
        ?(name_file(name1_,name2_,name3_)) // try again
        :(check_bad_encode("name_file"));
    }

    inline char * name_svd_file(const char dat_[], const char fam_[], const char a_[], const char suf_[] = "mat")
    {
      check_output_names();
      return (not_enough_space(check = snprintf(name,len_buf,"%s/%s_%s.%s%s_svd.%s",
                                                              dir_name,dat_,fam_,a_,suf_,dat_suf)))
        ?(name_svd_file(dat_,fam_,a_,suf_)) // try again
        :(check_bad_encode("name_svd_file"));
    }
    inline char * name_svd_file(LD_name_buffer &dat_,LD_name_buffer &fam_, const char a_[], const char suf_[] = "mat")
      {return name_svd_file(dat_.name,fam_.name,a_,suf_);}

    inline char * name_matrix_file(const char dat_[], const char fam_[], const char a_[], const char suf_[] = "mat")
    {
      check_output_names();
      return (not_enough_space(check = snprintf(name,len_buf,"%s/%s_%s.%s%s.%s",dir_name,dat_,fam_,a_,suf_,dat_suf)))
        ?(name_matrix_file(dat_,fam_,a_,suf_)) // try again
        :(check_bad_encode("name_matrix_file"));
    }
    inline char * name_matrix_file(LD_name_buffer &dat_,LD_name_buffer &fam_, const char a_[], const char suf_[] = "mat")
      {return name_matrix_file(dat_.name,fam_.name,a_,suf_);}

    inline char * name_domain_config_file(const char dat_[], const char fam_[])
    {
      check_output_names();
      return (not_enough_space( check = snprintf(name,len_buf,"%s/%s_%s.domain_config.%s",dir_name,dat_,fam_,dat_suf)))
        ?(name_domain_config_file(dat_,fam_)) // try again
        :(check_bad_encode("name_domain_config_file"));
    }
    inline char * name_domain_config_file(LD_name_buffer &dat_,LD_name_buffer &fam_)
      {return name_domain_config_file(dat_.name,fam_.name);}


    inline char * name_function_space(const char fam_[],int bor_)
    {
      return (not_enough_space( check = snprintf(name,len_buf,"%s.%d",fam_,bor_)))
        ?(name_function_space(fam_,bor_)) // try again
        :(check_bad_encode("name_function_space"));
    }
    inline char * name_observational_data(const char eqn_[],int xr_,int nse_,const char ign_[], const char adtl_[] = "")
    {
      return (not_enough_space( check = (nse_<0)
                                ?(snprintf(name,len_buf,"%s%s_xrange%d_true_%sgen",adtl_,eqn_,xr_,ign_))
                                :(snprintf(name,len_buf,"%s%s_xrange%d_nse%d_%sgen",adtl_,eqn_,xr_,nse_,ign_))))
        ?(name_observational_data(eqn_,xr_,nse_,ign_,adtl_)) // try again
        :(check_bad_encode("name_observational_data"));
    }


    inline char * load_name(const char name_[])
    {
      return (not_enough_space( check = snprintf(name,len_buf,"%s",name_)))
                                ?(load_name(name_))
                                :(check_bad_encode("load_name"));
    }
    inline char * load_name(const char name1_[],const char name2_[])
    {
      return (not_enough_space( check = snprintf(name,len_buf,"%s%s",name1_,name2_)))
                                ?(load_name(name1_,name2_))
                                :(check_bad_encode("load_name"));
    }
    inline char * load_name(const char name1_[],const char name2_[],const char name3_[])
    {
      return (not_enough_space( check = snprintf(name,len_buf,"%s%s%s",name1_,name2_,name3_)))
                                ?(load_name(name1_,name2_,name3_))
                                :(check_bad_encode("load_name"));
    }
    inline char * load_dir_name(const char dir_name_[] = "./data_directory")
    {
      if (dir_name != NULL)
      {
        if (strlen(dir_name) >= strlen(dir_name_)) return strcpy(dir_name,dir_name_);
        else
        {
          delete [] dir_name;
          return dir_name = LD_io::duplicate_string(dir_name_);
        }
      }
      else return dir_name = LD_io::duplicate_string(dir_name_);
    }
    inline char * load_dat_suf(const char dat_suf_[] = "lddat")
    {
      if (dat_suf != NULL)
      {
        if ( strlen(dat_suf) >= strlen(dat_suf_) ) return strcpy(dat_suf,dat_suf_);
        else
        {
          delete [] dat_suf;
          return dat_suf = LD_io::duplicate_string(dat_suf_);
        }
      }
      else return dat_suf = LD_io::duplicate_string(dat_suf_);
    }
};

#endif
