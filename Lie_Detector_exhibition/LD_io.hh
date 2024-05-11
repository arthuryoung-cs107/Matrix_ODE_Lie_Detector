#ifndef LD_IO_HH
#define LD_IO_HH

#include <cstring>
#include <cstdio>
#include <cstdlib>

struct LD_io
{
  LD_io() {}
  ~LD_io() {}
  // LD_io();
  // ~LD_io();

  static char * duplicate_string(const char in_[])
    {char * out = new char[strlen(in_)]; strcpy(out,in_); return out;}

  // static FILE * fopen_SAFE(char * const name_, const char tags_[])
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
};

#endif
