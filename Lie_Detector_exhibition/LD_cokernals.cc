#include "LD_cokernals.hh"

cokernal_sub_bundle::cokernal_sub_bundle(Jet_function_vector_space &jfvs_,cokernal_policy &pol_,bool verbose_):
  vlen_full(jfvs_.vlen_full),
  nset(jfvs_.nset0),
  cokern(jfvs_.Vbndle0,cokern_dimvec,pol_.init_cokernal_collapse(jfvs_,verbose_))
{
  const int nset0 = nset;

  int generation = 0;
  double t0 = LD_threads::tic();
  while (pol_.iterate_cokernal_reduction(cokern,generation,verbose_)) generation++;
  double work_time = LD_threads::toc(t0);

  const int nset_f = nset;

  LD_encoded_matrix ** const Aenc_mats = jfvs_.Acode.Acodes;
  LD_vector_space ** const krn0_spcs = jfvs_.Vbndle0.Vspaces;

  int nKf = 0,
      ckrn_V0_assign[nset0];
  for (size_t ickrn = 0; ickrn < nset_f; nKf+=ckrn_spcs[ickrn++]->nvec_use)
    for (size_t iset = 0; iset < ckrn_spcs[ickrn]->n_0Vkrn; iset++)
      ckrn_V0_assign[ckrn_spcs[ickrn]->ivec_0Vkrn[iset]] = ickrn;

  if (verbose_)
    printf("\n(cokernal_sub_bundle::cokernal_sub_bundle) collapsed %d kernal spaces (%d x %d) into %d (nKf = %d) in %.1f seconds (%d threads)\n",
      nset0, vlen_full, vlen_full, nset_f, nKf,
      work_time, LD_threads::numthreads() );

  LD_svd  Kf_svd(vlen_full,nKf);
  for (size_t iset = 0, jK = 0; iset < nset_f; iset++)
    for (size_t lcol = 0; lcol < ckrn_spcs[iset]->nvec_use; lcol++, jK++)
      for (size_t i = 0; i < vlen_full; i++) Kf_svd.Umat[i][jK] = ckrn_spcs[iset]->Vmat[lcol][i];

  Kf_svd.decompose_U(vlen_full,nKf);

  Kf_svd.print_result("Kf");
  int rank_Kf = Kf_svd.rank();
  double  nuke_norm_Kf = Kf_svd.norm_nuke(),
          effrank_Kf = nuke_norm_Kf/Kf_svd.s0;
  printf("  nKf = %d, rank_Kf = %d, nuke_norm_Kf = %.2f (%.2f), effective rank_Kf = %.2f (%.2f)\n",
    nKf,rank_Kf,
    nuke_norm_Kf, nuke_norm_Kf/((double)nKf),
    effrank_Kf,effrank_Kf/((double)nKf));

  k_medoids_package kmed_f(pol_.msr.dsym,nset_f);
  const int kSC_f = kSC = kmed_f.comp_kSC_medoids(2,verbose_);
  int nmem_clst_vec_f[kSC_f],
      imem_clst_vec_f[nset_f],
      *imem_clst_mat_f[kSC_f];
  kmed_f.assign_clusters(nmem_clst_vec_f,imem_clst_vec_f,imem_clst_mat_f,kSC_f,false);
  int n_0V_clst_f[kSC_f],
      nvec_clst_f[kSC_f],
      i_0V_list_f[nset0],
      *i_0V_clst_f[kSC_f],
      *i_0V_clst_med_f[kSC_f];
  printf("k = %d cluster assignments:\n", kSC_f);
  for (size_t k = 0, ii_0V = 0; k < kSC_f; ii_0V += n_0V_clst_f[k++])
  {
    const int n_clst_f_k = nmem_clst_vec_f[k],
              imed_k = kmed_f.i_meds[k];
    printf("\n(k=%d) imed = %d, nmem = %d.\n  i_ckrns : ", k, imed_k, n_clst_f_k);

    i_0V_clst_med_f[k] = ckrn_spcs[imed_k]->ivec_0Vkrn;

    int * const i_0V_clst_k = i_0V_clst_f[k] = i_0V_list_f + ii_0V,
        * const imem_k = imem_clst_mat_f[k];
    int i_ckrn_ki,
        &nvec_clst_k = nvec_clst_f[k] = 0,
        &n_0V_clst_k = n_0V_clst_f[k] = 0;
    for (size_t i = 0; i < n_clst_f_k; i++)
    {
      printf("%s%d%s ",(imem_k[i]==imed_k)?("["):(""),i_ckrn_ki = imem_k[i],(imem_k[i]==imed_k)?("]"):(""));
      for (size_t j = 0; j < ckrn_spcs[i_ckrn_ki]->n_0Vkrn; j++)
        i_0V_clst_k[n_0V_clst_k++] = ckrn_spcs[i_ckrn_ki]->ivec_0Vkrn[j];
      nvec_clst_k += ckrn_spcs[i_ckrn_ki]->nvec_use;
    }
    LD_linalg::sort_vec_inc<int>(i_0V_clst_k,n_0V_clst_k);
    printf("--> iV0 : ");
    for (size_t i = 0; i < n_0V_clst_k; i++) printf("%d ", i_0V_clst_k[i]);
    printf("\n  (tot. nV0 = %d, nvec = %d)\n  [ ickrn : iV0 | (nV0, nvec0->nvecf) ] - \n", n_0V_clst_k, nvec_clst_k);
    for (size_t i = 0; i < n_clst_f_k; i++)
    {
      printf("  [ %d : ", imem_k[i]);
      int nvec_0_net = 0;
      for (size_t j = 0; j < ckrn_spcs[imem_k[i]]->n_0Vkrn; j++)
      {
        printf("%d ", ckrn_spcs[imem_k[i]]->ivec_0Vkrn[j]);
        nvec_0_net += jfvs_.svd0.nV_spcvec[ckrn_spcs[imem_k[i]]->ivec_0Vkrn[j]];
      }
      printf(" | (%d, %d->%d) ] %s\n", ckrn_spcs[imem_k[i]]->n_0Vkrn, nvec_0_net,ckrn_spcs[imem_k[i]]->nvec_use,
                                      (imem_k[i]==imed_k)?(" (<-MED)"):(""));
    }
  }

  int mrow_Amat[nset0],
      nvec_krn0[nset0],
      nvec_ckrn[nset0];
  double  ** A_mats[nset0],
          ** V_krn0[nset0],
          ** V_ckrn[nset0],
          krn0_stats[nset0][4],
          ckrn_stats[nset0][4];

  #pragma omp parallel for
  for (size_t iset = 0; iset < nset0; iset++)
  {
    const int iickrn = ckrn_V0_assign[iset],
              mr_Amat_i = mrow_Amat[iset] = Aenc_mats[iset]->nrow,
              nV_krn0_i = nvec_krn0[iset] = krn0_spcs[iset]->nvec_use,
              nV_ckrn_i = nvec_ckrn[iset] = ckrn_spcs[iickrn]->nvec_use,
              len_Akrn0_i = mr_Amat_i*nV_krn0_i,
              len_Ackrn_i = mr_Amat_i*nV_ckrn_i;
    double  Ai_wkspc[len_Akrn0_i],
           *A_krn0_i[mr_Amat_i],
           *A_ckrn_i[mr_Amat_i];
    for (size_t i = 0, ii_krn0 = 0, ii_ckrn = 0; i < mr_Amat_i; i++, ii_krn0+=nV_krn0_i, ii_ckrn+=nV_ckrn_i)
    {
      A_krn0_i[i] = Ai_wkspc + ii_krn0;
      A_ckrn_i[i] = Ai_wkspc + ii_ckrn;
    }

    LD_linalg::A__B_CT( A_krn0_i,
                        A_mats[iset]=Aenc_mats[iset]->Amat,
                        V_krn0[iset]=krn0_spcs[iset]->Vmat,
                        mr_Amat_i,vlen_full,nV_krn0_i,true);
    LD_linalg::abs_vec<double>(A_krn0_i[0],len_Akrn0_i);
    LD_linalg::comp_Tvec_stats<double>(krn0_stats[iset],A_krn0_i[0],len_Akrn0_i);

    LD_linalg::A__B_CT( A_ckrn_i,
                        A_mats[iset],
                        V_ckrn[iset]=ckrn_spcs[iickrn]->Vmat,
                        mr_Amat_i,vlen_full,nV_ckrn_i,true);
    LD_linalg::abs_vec<double>(A_ckrn_i[0],len_Ackrn_i);
    LD_linalg::comp_Tvec_stats<double>(ckrn_stats[iset],A_ckrn_i[0],len_Ackrn_i);
  }

  printf("\n  original kernal vs. cokernal inner product performance\n");
  for (size_t iset = 0; iset < nset0; iset++)
    printf("    (crv %d, iickrn %d, krn0 vs. ckrn):  nV (%d -> %d) ; min (%.2e -> %.2e) ; med (%.2e -> %.2e) ; avg (%.2e -> %.2e) ; max (%.2e -> %.2e) ; s0 (%.2e -> %.2e) ; |.|_* (%.2e -> %.2e) \n",
                iset, ckrn_V0_assign[iset], krn0_spcs[iset]->nvec_use, ckrn_spcs[ckrn_V0_assign[iset]]->nvec_use,
                krn0_stats[iset][0], ckrn_stats[iset][0],
                krn0_stats[iset][1], ckrn_stats[iset][1],
                krn0_stats[iset][2], ckrn_stats[iset][2],
                krn0_stats[iset][3], ckrn_stats[iset][3],
                jfvs_.svd0.Smat[iset][0], ckrn_spcs[ckrn_V0_assign[iset]]->wvec[0],
                LD_linalg::sum_vec<double>(jfvs_.svd0.Smat[iset],vlen_full),
                LD_linalg::sum_vec<double>(ckrn_spcs[ckrn_V0_assign[iset]]->wvec,vlen_full)
              );

}

int nullspace_ckrn_policy::iterate_nullspace_cokernals(cokernal_bundle &cokern_,int gen_, k_medoids_package &kmed_, bool verbose_)
{
  const int vlen_full = msr.vlen;
  int &nset = cokern_.nset,
      &kSC = cokern_.kSC,
      &nred_succ = cokern_.nred_succ,
      nmem_clst_vec[kSC],
      imem_clst_vec[nset],
      *imem_clst_mat[kSC];

  kmed_.assign_clusters(nmem_clst_vec,imem_clst_vec,imem_clst_mat,kSC,verbose_);

  cokernal_space ** const ckrn_spcs = cokern_.ckrn_spcs;
  bool f_reduce_success[kSC];
  int nsvd_cmps = nred_succ = 0,
      nuld_ckrn[kSC],
      isat_ckrn[kSC][vlen_full],
      i_n_pairs_ckrn[kSC][2],
      ij_pairs_ckrn[kSC][2],
      ij_nn_pairs[kSC][2];
  double  dpairs_ckrn[kSC][2],
          svec_ckrns[kSC][vlen_full],
          Vmat_ckrns[kSC][vlen_full][vlen_full],
          t0 = LD_threads::tic();
  #pragma omp parallel reduction(+:nsvd_cmps,nred_succ)
  {
    LD_svd  Vsvd_t(2*vlen_full,vlen_full);
    int i_setk,j_setk;
    double  * const svec_V_t = Vsvd_t.svec,
            ** const Vmat_V_t = Vsvd_t.Vmat,
            ** const Umat_V_t = Vsvd_t.Umat;
    #pragma omp for
    for (size_t k = 0; k < kSC; k++)
    {
      if (f_reduce_success[k] = (nmem_clst_vec[k] > 1))
      {
        const int npairs_k = i_n_pairs_ckrn[k][1] = (nmem_clst_vec[k])*(nmem_clst_vec[k]-1)/2;
        int ipair = i_n_pairs_ckrn[k][0] = 0,
            ipairs_k[npairs_k],
            iipairs_k[npairs_k][2];
        double dpairs_k[npairs_k];
        kmed_.comp_sort_cluster_distances(dpairs_k,ipairs_k,iipairs_k[0],imem_clst_mat[k],nmem_clst_vec[k]);
        ij_nn_pairs[k][0] = iipairs_k[0][0]; ij_nn_pairs[k][1] = iipairs_k[0][1];
        dpairs_ckrn[k][1] = dpairs_k[0];
        do
        {
          i_setk = iipairs_k[ipair][0]; j_setk = iipairs_k[ipair][1];

          for (size_t i = 0; i < vlen_full; i++)
            for (size_t j = 0; j < vlen_full; j++)
              Umat_V_t[i][j] = (ckrn_spcs[i_setk]->wvec[i])*(ckrn_spcs[i_setk]->Vmat_data[i][j]);

          for (size_t i = 0; i < vlen_full; i++)
            for (size_t j = 0; j < vlen_full; j++)
              Umat_V_t[i+vlen_full][j] = (ckrn_spcs[j_setk]->wvec[i])*(ckrn_spcs[j_setk]->Vmat_data[i][j]);

          Vsvd_t.decompose_U(2*vlen_full,vlen_full); nsvd_cmps++;

          if (f_reduce_success[k] = (bool)( nuld_ckrn[k] = (vlen_full-Vsvd_t.rank()) ))
          {
            nred_succ++;
            for (size_t i = 0, isat = vlen_full-nuld_ckrn[k]; i < nuld_ckrn[k]; i++, isat++) isat_ckrn[k][i] = isat;
            for (size_t i = nuld_ckrn[k], insat = 0; i < vlen_full; i++, insat++) isat_ckrn[k][i] = insat;
            i_n_pairs_ckrn[k][0] = ipair;
            ij_pairs_ckrn[k][0] = i_setk; ij_pairs_ckrn[k][1] = j_setk;
            dpairs_ckrn[k][0] = dpairs_k[ipair];
            for (size_t i = 0; i < vlen_full; i++) svec_ckrns[k][i] = svec_V_t[i];
            for (size_t i = 0; i < vlen_full; i++)
              for (size_t j = 0; j < vlen_full; j++)
                Vmat_ckrns[k][i][j] = Vmat_V_t[j][i];
            break;
          }
          else if (ipair == (npairs_k-1)) break;
          else ipair++;
        } while (true);
      }
      else
      {
        ij_nn_pairs[k][0] = ij_nn_pairs[k][1] = kmed_.i_meds[k];
        dpairs_ckrn[k][0] = dpairs_ckrn[k][1] = 0.0;
        i_n_pairs_ckrn[k][0] = -1;
        i_n_pairs_ckrn[k][1] = 0;
      }
    }
  }
  double work_time = LD_threads::toc(t0);

  if (verbose_)
  {
    printf("  (nullspace_ckrn_policy::iterate_nullspace_cokernals) gen %d: reduced %d (of kSC = %d) cokernals via %d svds in %.2f seconds (%d threads).\n  Reduction details: \n",
      gen_,nred_succ, kSC,
      nsvd_cmps,work_time,LD_threads::numthreads() );
    for (size_t k = 0; k < kSC; k++)
    {
      int imd=kmed_.i_meds[k], inn=ij_nn_pairs[k][0], jnn=ij_nn_pairs[k][1];
      printf("    k=%d: imed=%d (nv=%d), nclst=%d. inn=%d, jnn=%d (nv = %d, %d ; dij=%.2e) of %d pairs. Intersection %s",
        k,imd,ckrn_spcs[imd]->nvec_use,nmem_clst_vec[k],
        inn,jnn,ckrn_spcs[inn]->nvec_use,ckrn_spcs[jnn]->nvec_use,dpairs_ckrn[k][1],
        i_n_pairs_ckrn[k][1],
        (f_reduce_success[k])?( (i_n_pairs_ckrn[k][0])?("FOUND: "):("FOUND (default).\n") ):("NOT found.\n"));
      if (f_reduce_success[k] && i_n_pairs_ckrn[k][0])
        printf("i_cok=%d, j_cok=%d (d=%.2e, %d).\n",
          ij_pairs_ckrn[k][0],ij_pairs_ckrn[k][1],dpairs_ckrn[k][0],i_n_pairs_ckrn[k][0]);
    }
  }

  if (nred_succ)
    nset = cokern_.collapse_cokernals(msr,gen_+1,
      f_reduce_success,nuld_ckrn,ij_pairs_ckrn[0],svec_ckrns[0],Vmat_ckrns[0][0],isat_ckrn[0],wdistance,verbose_);

  return nred_succ;
}

int cokernal_bundle::collapse_cokernals(LD_vspace_measure &msr_, int gen_spawn_, bool *frdc_succ_, int *nsat_ckrn_, int *ij_prs_, double *w_ckrn_, double *vvec_ckrn_, int *isat_ckrn_,bool wdistance_,bool verbose_)
{
  const int vl2 = vlen*vlen;
  int i_i, i_j,
      ij_prnts[nred_succ][2],
      gen_iprnt[nred_succ],
      gen_jprnt[nred_succ],
      nvec_use_iprnt[nred_succ],
      nvec_use_jprnt[nred_succ],
      nV0krn_iprnt[nred_succ],
      nV0krn_jprnt[nred_succ];
  cokernal_space * ckrn_rdc[nred_succ];
  for (size_t i = 0, irdc = 0, i_ij = 0, iw = 0, iv = 0; irdc < nred_succ; irdc += (int)(frdc_succ_[i++]), i_ij+=2, iw+=vlen, iv+=vl2)
    if (frdc_succ_[i])
    {
      cokernal_space  * const ckrn_i_old = ckrn_spcs[i_i = ij_prnts[irdc][0] = ij_prs_[i_ij]],
                      * &ckrn_j_old = ckrn_spcs[i_j = ij_prnts[irdc][1] = ij_prs_[i_ij+1]];
      gen_iprnt[irdc] = ckrn_i_old->gen_spawn; gen_jprnt[irdc] = ckrn_j_old->gen_spawn;
      nvec_use_iprnt[irdc] = ckrn_i_old->nvec_use; nvec_use_jprnt[irdc] = ckrn_j_old->nvec_use;
      nV0krn_iprnt[irdc] = ckrn_i_old->n_0Vkrn; nV0krn_jprnt[irdc] = ckrn_j_old->n_0Vkrn;

      if ( ckrn_i_old->gen_spawn ) // if NOT an original kernal, can be overwritten
        ckrn_spcs[i_i] = new cokernal_space(gen_spawn_,ckrn_i_old,ckrn_j_old); // spawn new cokernal space by overwriting ckrn_i and its workspace
      else // if an original kernal, extra workspace must be overwritten or allocated
      {
        int i_check = 0;
        while ( (i_check<nxtra_wkspc) && (factive_ckrn_wkspc_xtra[i_check]) ) i_check++;
        if (i_check==nxtra_wkspc) // if all allocated workspaces are active, allocate a new one.
          ckrn_wkspcs_xtra[nxtra_wkspc++] = new cokernal_workspace(nset0+nxtra_wkspc,factive_ckrn_wkspc_xtra[nxtra_wkspc]=true,vlen);

        // spawn new cokernal space at ckrn_i using an allocated workspace
        ckrn_spcs[i_i] = new cokernal_space(gen_spawn_,ckrn_wkspcs_xtra[i_check],ckrn_i_old,ckrn_j_old);
      }
      // set the underlying vector space
      (ckrn_rdc[irdc] = ckrn_spcs[i_i])->init_cokernal_space(vlen,nsat_ckrn_[i],w_ckrn_+iw,vvec_ckrn_+iv,isat_ckrn_+iw);
    }

  const int nset_new = nset-nred_succ;
  int nV_setvec_new[nset_new];
  double * Wmat_new[nset_new];
  LD_vector_space * vspcs_new[nset_new];
  for (size_t i = 0, ickrn_new = 0; ickrn_new < nset_new; i++)
    if (ckrn_spcs[i] != NULL)
    {
      ckrn_spcs[i]->ickrn = ickrn_new;
      nV_setvec_new[ickrn_new] = ckrn_spcs[i]->nvec_use;
      Wmat_new[ickrn_new] = ckrn_spcs[i]->wvec;
      vspcs_new[ickrn_new] = ckrn_wkspcs[ckrn_spcs[i]->iwkspc]->vspc;
      ckrn_spcs[ickrn_new++] = ckrn_spcs[i];
    }
  for (size_t i = nset_new; i < nset0; i++) ckrn_spcs[i] = NULL;

  if (verbose_)
  {
    printf("  (cokernal_bundle::collapse_cokernals) reassigned %d cokernals:\n",nset_new);
    for (size_t i = 0; i < nset_new; i++)
    {
      cokernal_space &cki = *(ckrn_spcs[i]);
      printf("    (ickrn=%d, iwkspc=%d): nvec=%d, n_0Vkrn=%d, (gen,iprnt,jprnt)=(%d,%d,%d)",
        cki.ickrn, cki.iwkspc,
        cki.nvec_use, cki.n_0Vkrn,
        cki.gen_spawn,cki.ickrn_parent,cki.jckrn_parent);
      if (cki.n_0Vkrn > 1)
      {
        printf(", i_0Vkrn: ( ");
        for (size_t iV0 = 0; iV0 < cki.n_0Vkrn; iV0++) printf("%d ", cki.ivec_0Vkrn[iV0]);
        printf(")\n");
      }
      else printf(".\n");
    }

    printf("    New cokernals: \n");
    for (size_t irdc = 0; irdc < nred_succ; irdc++)
    {
      cokernal_space &ckrni = *(ckrn_rdc[irdc]);
      printf("      (irdc=%d) [iVprnt=%d, jVprnt=%d -> ickrn=%d] : gen [%d,%d -> %d] ; nV [%d,%d -> %d] ; nV0krn [%d,%d -> %d]. iV0krn = ( ",
        irdc,
        ckrni.ickrn_parent,ckrni.jckrn_parent,ckrni.ickrn,
        gen_iprnt[irdc],gen_jprnt[irdc],ckrni.gen_spawn,
        nvec_use_iprnt[irdc],nvec_use_jprnt[irdc],ckrni.nvec_use,
        nV0krn_iprnt[irdc],nV0krn_jprnt[irdc],ckrni.n_0Vkrn);
      for (size_t i = 0; i < ckrni.n_0Vkrn; i++) printf("%d ", ckrni.ivec_0Vkrn[i]);
      printf(").\n");
    }
  }

  // if (wdistance_) msr_.init_distances(vspcs_new,Wmat_new,nset_new,verbose_); // recompute (expensive), but definitely correct weighted distances
  // else msr_.init_distances(vspcs_new,nset_new,verbose_); // recompute (expensive), but definitely correct distances

  if (wdistance_) msr_.update_distances(nred_succ,ij_prnts[0],vspcs_new,Wmat_new,verbose_); // update weighted distances
  else msr_.update_distances(nred_succ,ij_prnts[0],vspcs_new,verbose_); // update distances

  return nset_new;
}
