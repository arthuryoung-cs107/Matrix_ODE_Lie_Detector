#ifndef LD_COKERN_HH
#define LD_COKERN_HH

#include "LD_cokernals.hh"

struct cokernal_vfield_workspace
{
  cokernal_vfield_workspace(cokernal_refinement &rfne_) :
    nset0(rfne_.nset0), nsetf(rfne_.nset)

    {}
  ~cokernal_vfield_workspace() {}

  const int nset0,
            nsetf;

};

// class cokernal_vfield : public LD_vector_field, public cokernal_vfield_workspace
// {
//
// };

#endif
