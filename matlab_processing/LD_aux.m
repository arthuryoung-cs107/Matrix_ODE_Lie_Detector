classdef LD_aux
    properties

    end
    methods
        function obj = LD_aux()

        end
    end
    methods (Static)
        function mat_out = read_Tmat(name_)
            file = fopen(name_);
            hlen = fread(file,1,'int=>int');
            header = fread(file,hlen,'int=>int');
            switch header(3)
                case 8
                    mat_flat = fread(file,header(1)*header(2),'double=>double');
                case 4
                    mat_flat = fread(file,header(1)*header(2),'int=>int');
                otherwise
                    fprintf('(LD_aux::read_Tmat) ERROR - failed to read %s, non-standard byte size(%d) \n', name_,header(3));
            end
            fclose(file);

            mat_out = reshape(mat_flat,header(2),header(1))';
            fprintf('(LD_aux::read_Tmat) read %s \n', name_);
        end
        function AYLglbsvd_glb_out = make_global_restricted_svd_package(Amat_,Lglbsvd_,kappa_L_)
            [ncol,ncrv] = deal(Amat_.ncol,length(Lglbsvd_.rvec));
            AYLglb_mat = LD_aux.Atns_Ytns_mult(Amat_.matT',Lglbsvd_.Vmat_glb(:,1:(end-kappa_L_)));
            ncol_AYLglb = size(AYLglb_mat,2);
            AYLglb_tns = permute(reshape(AYLglb_mat',ncol_AYLglb,[],ncrv),[2 1 3]);

            [~,svec_glb_AYLglb,Vmat_glb_AYLglb] = svd(AYLglb_mat,'econ','vector');

            rvec_AYLglb = nan(ncrv,1);
            smat_AYLglb = nan(ncol_AYLglb,ncrv);
            Vtns_AYLglb = nan(ncol_AYLglb,ncol_AYLglb,ncrv);
            for i = 1:ncrv
                rvec_AYLglb(i) = rank(AYLglb_tns(:,:,i));
                [~,smat_AYLglb(:,i),Vtns_AYLglb(:,:,i)] = svd(AYLglb_tns(:,:,i),'econ','vector');
            end

            AYLglbsvd_out = Lglbsvd_;
            AYLglbsvd_glb_out.nrow = Amat_.nrow;
            AYLglbsvd_glb_out.ncol = ncol_AYLglb;
            AYLglbsvd_glb_out.matT = AYLglb_mat';
            AYLglbsvd_glb_out.rvec = rvec_AYLglb;
            AYLglbsvd_glb_out.smat = smat_AYLglb;
            AYLglbsvd_glb_out.Vtns = Vtns_AYLglb;
            AYLglbsvd_glb_out.svec_glb = svec_glb_AYLglb;
            AYLglbsvd_glb_out.Vmat_glb = Vmat_glb_AYLglb;
        end
        function pckgs_out = make_global_svd_package(pckgs_in_)
            len_pckgs = length(pckgs_in_);
            for i = 1:len_pckgs
                pckg_i_new = pckgs_in_(i);
                [~,pckg_i_new.svec_glb,pckg_i_new.Vmat_glb] = svd(pckg_i_new.matT','econ','vector');
                pckgs_out(i) = pckg_i_new;
            end
        end
        function Asvd_rstr_out = comp_restricted_space(meta_,Lsvd_,Asvd_)
            Asvd_rstr_out = Asvd_;

            nvar = meta_.ndep + 1;
            ncrv = length(Lsvd_.rvec);
            ndof = Asvd_.ncol;
            rho_L = min(Lsvd_.rvec);
            kappa_A = ndof - max(Asvd_.rvec);
            ncols_PL = nvar*rho_L;

            Atns = permute(reshape(Asvd_.matT,ndof,[],ncrv),[2 1 3]);

            AYtns = LD_aux.Atns_Ytns_mult(Atns,Lsvd_.Vtns(:,1:rho_L,:));
            VYtns = LD_aux.Atns_Ytns_mult(Asvd_.Vtns,Lsvd_.Vtns(:,1:rho_L,:));

            rvec_AY = nan(ncrv,1);
            smat_AY = nan(ncols_PL,ncrv);
            Vtns_AY = nan(ncols_PL,ncols_PL,ncrv);
            for i = 1:ncrv
                rvec_AY(i) = rank(AYtns(:,:,i));
                [~,smat_AY(:,i),Vtns_AY(:,:,i)] = svd(AYtns(:,:,i),'econ','vector');
            end

            Asvd_rstr_out.AYtns = AYtns;
            Asvd_rstr_out.rvec_AY = rvec_AY;
            Asvd_rstr_out.smat_AY = smat_AY;
            Asvd_rstr_out.Vtns_AY = Vtns_AY;
            Asvd_rstr_out.VYtns = VYtns;
        end
        function AYtns_out = Atns_Ytns_mult(Atns_,Ytns_)
            [nrow,ndof,ncrv] = size(Atns_);
            [perm_len,rho,sz3Y] = size(Ytns_);
            nvar = ndof/perm_len;
            N_PL = nvar*rho;

            switch sz3Y
                case ncrv
                    AYtns_out = nan(nrow,N_PL,ncrv);
                    for i = 1:ncrv
                        AYtns_out(:,:,i) = reshape(pagemtimes(reshape(Atns_(:,:,i), ...
                                        nrow,perm_len,nvar),Ytns_(:,:,i)),nrow,[]);
                    end
                case 1
                    AYtns_out = nan(nrow,N_PL,ncrv);
                    for i = 1:ncrv
                        AYtns_out(:,:,i) = reshape(pagemtimes(reshape(Atns_(:,:,i),nrow,perm_len,nvar),Ytns_),nrow,[]);
                    end
            end
        end
        function YVtns_out = Ytns_Vtns_mult(Ytns_,Vtns_)
            [perm_len,rho,sz3Y] = size(Ytns_);
            [N_PL,kappa,ncrv] = size(Vtns_);
            nvar = N_PL/rho;
            ndof = nvar*perm_len;

            switch sz3Y
                case ncrv
                    YVtns_out = nan(ndof,kappa,ncrv);
                    for i = 1:ncrv
                        YVtns_out(:,:,i) = reshape(pagemtimes(Ytns_(:,:,i), ...
                                        reshape(Vtns_(:,:,i),rho,nvar,kappa)),ndof,kappa);
                    end
                case 1
                    YVtns_out = nan(ndof,kappa,ncrv);
                    for i = 1:ncrv
                        YVtns_out(:,:,i) = reshape(pagemtimes(Ytns_,reshape(Vtns_(:,:,i),rho,nvar,kappa)),ndof,kappa);
                    end
            end
        end
        function str1_out = overwrite_struct(str1_,str2_)
            str1_out = str1_;
            str1_names = fieldnames(str1_);
            for i = 1:length(str1_names)
                name_i = str1_names{i};
                str1_out = setfield(str1_out,name_i,getfield(str2_,name_i));
            end
        end
        function struct_out = combine_structs(str1_,str2_)
            struct_out = cell2struct([struct2cell(str1_); struct2cell(str2_)], [fieldnames(str1_); fieldnames(str2_)]);
        end
        function [min_out,med_out,avg_out,max_out] = mat_colstats(mat_)
            [min_out,med_out,avg_out,max_out] = deal(min(mat_,[],1),median(mat_,1),mean(mat_,1),max(mat_,[],1));
        end
        function [min_out,med_out,avg_out,max_out] = mat_rowstats(mat_)
            [min_out,med_out,avg_out,max_out] = LD_aux.mat_colstats(mat_');
        end
        function [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d] = naive_k_medoid(dmat_,kfull_)
            [i_meds_0,n_meds_0,net_d_0,clst_mem_0,loc_d_0,med2med_d_0] = LD_aux.greedy_naive_k_medoids(dmat_,kfull_);
            [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d] = deal(i_meds_0,n_meds_0,net_d_0,clst_mem_0,loc_d_0,med2med_d_0);

            npts = size(dmat_,1);
            i_full = 1:npts;

            get_not_medoids = @(i_medoids_) i_full(logical(prod( i_full ~= i_medoids_', 1 )));
            augment_medoids = @(i_medoids_,i_) deal([i_medoids_,i_],get_not_medoids([i_medoids_,i_]));

            while (true)
                i_nmedoids = get_not_medoids(i_meds);
                [mswap_best,oswap_best] = deal(0);
                dswap_best = net_d;
                for k = 1:kfull_
                    i_nkmedoids = [i_meds(1:(k-1)) , i_meds((k+1):end)];
                    for inm = i_nmedoids
                        [i_medoids_try,i_nmedoids_try] = augment_medoids(i_nkmedoids,inm);
                        [d_mins_try,imed_assignments_try] = min(dmat_(i_medoids_try,i_nmedoids_try) ,[],1);
                        dswap_try = sum(d_mins_try);
                        if (dswap_best>dswap_try)
                            mswap_best = k;
                            oswap_best = inm;
                            dswap_best = dswap_try;
                        end
                    end
                end
                if (mswap_best == 0)
                    break;
                else
                    i_meds(mswap_best) = oswap_best;
                    net_d = dswap_best;
                end
            end
            if (nargout > 1)
                [d_mins,imed_assignments] = min( dmat_(i_meds,i_nmedoids) ,[],1);
                n_meds(1:kfull_) = sum(imed_assignments == (1:kfull_)',2);
                if (nargout > 2)
                    net_d = sum(d_mins);
                    [clst_mem,loc_d] = deal(zeros(1,npts));
                    clst_mem(i_nmedoids) = imed_assignments;
                    loc_d(i_nmedoids) = d_mins;
                    med2med_d = dmat_(i_meds,i_meds);
                end
            end
        end
        function [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d] = greedy_naive_k_medoids(dmat_,kfull_)
            npts = size(dmat_,1);
            i_full = 1:npts;
            i_meds = 1:kfull_;
            [~,i_meds(1)] = min(sum(dmat_,1),[],2); % first medoid

            get_not_medoids = @(i_medoids_) i_full(logical(prod( i_full ~= i_medoids_', 1 )));
            augment_medoids = @(i_medoids_,i_) deal([i_medoids_,i_],get_not_medoids([i_medoids_,i_]));

            for k = 2:kfull_
                [i_medoids_k,i_nmedoids_k] = augment_medoids(i_meds(1:(k-2)),i_meds(k-1));

                ik = i_nmedoids_k(1);
                imed_assignments_ik = ones(1,length(i_nmedoids_k));
                d_net_ik = realmax;

                for kk = 1:length(i_nmedoids_k)
                    ikk = i_nmedoids_k(kk);
                    [i_medoids_ikk,i_nmedoids_ikk] = augment_medoids(i_medoids_k,ikk);
                    [d_mins_ikk,imed_assignments_ikk] = min( dmat_(i_medoids_ikk,i_nmedoids_ikk) ,[],1);
                    d_net_ikk = sum(d_mins_ikk);
                    if (d_net_ik > d_net_ikk)
                        ik = ikk;
                        imed_assignments_ik = imed_assignments_ikk;
                        d_net_ik = d_net_ikk;
                    end
                end
                i_meds(k) = ik;
            end
            if (nargout > 1)
                i_nmedoids_greedy = get_not_medoids(i_meds);
                [d_mins_greedy,imed_assignments_greedy] = min( dmat_(i_meds,i_nmedoids_greedy) ,[],1);
                n_meds(1:kfull_) = sum(imed_assignments_greedy == (1:kfull_)',2);
                if (nargout > 2)
                    net_d = sum(d_mins_greedy);
                    [clst_mem,loc_d] = deal(zeros(1,npts));
                    clst_mem(i_nmedoids_greedy) = imed_assignments_greedy;
                    loc_d(i_nmedoids_greedy) = d_mins_greedy;
                    med2med_d = dmat_(i_meds,i_meds);
                end
            end
        end
    end
end
