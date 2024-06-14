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
    {char * out = new char[strlen(in_)]; strcpy(out,in_); return out;}

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

#endif
