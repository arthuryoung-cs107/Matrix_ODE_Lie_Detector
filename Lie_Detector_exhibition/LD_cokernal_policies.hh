#ifndef LD_CKRNPOL_HH
#define LD_CKRNPOL_HH

#include "LD_cokernals.hh"

class nullspace_ckrn_policy: public LD_cokernal_policy
{
  public:

    nullspace_ckrn_policy(LD_vspace_measure &msr_,bool wdistance_): LD_cokernal_policy(msr_,wdistance_) {}
    ~nullspace_ckrn_policy() {}

    virtual double ** init_cokernal_collapse(Jet_function_vector_space &jfvs_,bool verbose_=true)
    {
      jfvs_.svd0.compute_Acode_curve_svds(jfvs_.Acode,verbose_);
      jfvs_.svd0.evaluate_iscl_ranks(2*(jfvs_.vlen_full)); // using 2*N, debiasing tall matrices
      // jfvs_.svd0.evaluate_iscl_ranks(jfvs_.fspc.perm_len); // consider standardizing w perm_len
      if (verbose_) jfvs_.svd0.print_details("svd0");

      jfvs_.svd0.set_Vspaces_nullspc();
      if (verbose_) jfvs_.rec0.print_subspace_details("rec0",false,false);

      if (wdistance) msr.init_distances(jfvs_.Vbndle0.Vspaces,jfvs_.svd0.Smat,jfvs_.nset0,verbose_);
      else msr.init_distances(jfvs_.Vbndle0.Vspaces,jfvs_.nset0,verbose_);

      return jfvs_.svd0.Smat;
    }

  protected:

    int iterate_nullspace_cokernals(cokernal_bundle &cokern_,k_medoids_package &kmed_,bool verbose_);

};

class nullspace_clst_policy: public nullspace_ckrn_policy
{
  public:

    nullspace_clst_policy(LD_vspace_measure &msr_,bool wdistance_): nullspace_ckrn_policy(msr_,wdistance_) {}
    ~nullspace_clst_policy() {}

    virtual int reduce_cokernal_bundle(cokernal_bundle &cokern_, bool verbose_=true)
    {
      const int nset0 = cokern_.nset,
                nvK0 = cokern_.nvK();
      int &nset = cokern_.nset,
          &kSC = cokern_.kSC,
          &nred_succ = cokern_.nred_succ = 0,
          &generation = cokern_.generation = 0;

      double ** const dsym = msr.dsym;

      double t0 = LD_threads::tic();
      do
      {
        k_medoids_package kmed(dsym,nset);
        printf("\n");
        int kSC0_gen = kSC = kmed.comp_kSC_medoids(2,verbose_);
        if (iterate_nullspace_cokernals(cokern_,kmed,verbose_)) generation++;
        else
        {
          if (kSC0_gen > 2)
          {
            printf("\n(nullspace_clst_policy::iterate_cokernal_reduction) gen %d - IRREDUCABLE kSC0=%d clusters (nset=%d). Attempting to reduce 2 <= kSC < kSC0 = %d clusters.\n",generation,kSC0_gen,nset,kSC0_gen);
            int kSCm1 = kSC0_gen-1;
            kSC = kmed.comp_kSC_krange_medoids(2,kSCm1,verbose_);
            while ( (kSCm1 >= 2) && !(iterate_nullspace_cokernals(cokern_,kmed,verbose_)) )
              if (kSCm1 == 2) break;
              else kSC = kmed.comp_kSC_krange_medoids(2,--kSCm1,verbose_);
          }

          if (nred_succ) generation++;
          else
          {
            printf("\n(nullspace_clst_policy::iterate_cokernal_reduction) gen %d - IRREDUCABLE CLUSTERS (kSC0=%d, nset=%d). Attempting to consolidate remaining cokernal spaces\n",generation,kSC0_gen,nset);
            // if all else fails, just try consolidating one pair
            kmed.set_one_medoid();
            kSC = 1;
            if (iterate_nullspace_cokernals(cokern_,kmed,verbose_)) generation++;
            else kSC = kSC0_gen; // record the last authentic cluster count.
          }
        }
      } while (nred_succ);
      double work_time = LD_threads::toc(t0);

      const int nvKf = cokern_.nvK();

      if (verbose_)
        printf("\n\n(nullspace_clst_policy::reduce_cokernal_bundle) collapsed %d (nvK0 = %d) kernal spaces into %d (nvKf = %d, kSC = %d) in %.1f seconds (%d threads)\n",
          nset0, nvK0, nset, nvKf, kSC,
          work_time, LD_threads::numthreads() );

      return nvKf;
    }
};

class nullspace_near_policy: public nullspace_ckrn_policy
{
  public:

    nullspace_near_policy(LD_vspace_measure &msr_,bool wdistance_): nullspace_ckrn_policy(msr_,wdistance_) {}
    ~nullspace_near_policy() {}

    virtual int reduce_cokernal_bundle(cokernal_bundle &cokern_,bool verbose_=true)
    {
      const int nset0 = cokern_.nset,
                nvK0 = cokern_.nvK();
      int &nset = cokern_.nset,
          &kSC = cokern_.kSC = 1,
          &nred_succ = cokern_.nred_succ = 0,
          &generation = cokern_.generation = 0;

      double ** const dsym = msr.dsym;

      double t0 = LD_threads::tic();
      do
      {
        k_medoids_package kmed(dsym,nset);
        printf("\n");
        kmed.set_one_medoid();
        if (iterate_nullspace_cokernals(cokern_,kmed,verbose_)) generation++;
      } while (nred_succ);
      double work_time = LD_threads::toc(t0);

      k_medoids_package kmed_f(dsym,nset);
      kSC = kmed_f.comp_kSC_medoids(2,verbose_);

      const int nvKf = cokern_.nvK();

      if (verbose_)
        printf("\n\n(nullspace_near_policy::reduce_cokernal_bundle) collapsed %d (nvK0 = %d) kernal spaces into %d (nvKf = %d, kSC = %d) in %.1f seconds (%d threads)\n",
          nset0, nvK0, nset, nvKf, kSC,
          work_time, LD_threads::numthreads() );
      return nvKf;
    }

};


#endif
