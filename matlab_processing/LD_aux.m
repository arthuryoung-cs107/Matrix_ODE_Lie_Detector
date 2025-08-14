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
            if (file==-1)
                fprintf('(LD_aux::read_Tmat) : WARNING - failed to read %s. Defaulting to empty output \n',name_);
                mat_out = [];
            else
                hlen = fread(file,1,'int=>int');
                header = fread(file,hlen,'int=>int');
                switch header(3)
                    case 8
                        mat_flat = fread(file,header(1)*header(2),'double=>double');
                    case 4
                        mat_flat = fread(file,header(1)*header(2),'int=>int');
                    otherwise
                        fprintf('(LD_aux::read_Tmat) WARNING - non-standard byte size(%d) in %s \n',header(3), name_);
                        mat_flat = fread(file,header(1)*header(2)*header(3),'char*1=>char*1');
                end
                fclose(file);
                mat_out = reshape(mat_flat,header(2),header(1))';
                fprintf('(LD_aux::read_Tmat) read %s \n', name_);
            end
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
        function Ktns_A_YLrst_out = compute_Ktns_A_YLrst(Ytns_L_,Ktns_AYL_,rho_AYL_)
            if (nargin == 1)
                svd = Ytns_L_;
                [Ytns_L,Ktns_AYL] = deal(svd.Ytns_L,svd.Vtns(:,(max(svd.rvec)+1):end,:));
            elseif (nargin == 2)
                [Ytns_L,Ktns_AYL] = deal(Ytns_L_,Ktns_AYL_);
            elseif (nargin == 3)
                [Ytns_L,Ktns_AYL] = deal(Ytns_L_,Ktns_AYL_(:,(rho_AYL_+1):end,:));
            end
            Ktns_A_YLrst_out = LD_aux.Ytns_Vtns_mult(Ytns_L,Ktns_AYL);
        end
        % function pckg_out = compute_frobenius_closeness_matrix(Ktns_)
        function [KTK_mat,minds_lw] = compute_frobenius_closeness_matrix(Ktns_)
            [ncol,kappa,ncrv] = size(Ktns_);
            Kmat = reshape(Ktns_,ncol,[]);

            mflag_lw = logical(tril(ones(ncrv),-1)); % flags of lower triangular elements
            vflag_lw = reshape(mflag_lw,[],1); % column stack of lower triangular flags
            vinds_full = (1:(ncrv*ncrv))'; % 1 thru ncrv^2
            vinds_full_mat = reshape(vinds_full,ncrv,ncrv); % matrix of column stacked indices
            vinds_lw = vinds_full(vflag_lw); % column stacked indices of strictly lower triangular elements

            n_lw = length(vinds_lw);

            sqrt_kappa = sqrt(kappa);

            minds_lw = nan(2,n_lw);
            cvec_fro = nan(n_lw,1);
            for iv = 1:n_lw
                vind = vinds_lw(iv);
                [irow,icol] = find(vinds_full_mat==vind);
                minds_lw(:,iv) = [irow;icol];
                cvec_fro(iv) = norm((Ktns_(:,:,irow)')*Ktns_(:,:,icol), 'fro')/sqrt_kappa;
            end

            KTK_lwmat = zeros(ncrv);
            KTK_lwmat(vinds_lw) = cvec_fro;
            KTK_mat = KTK_lwmat + eye(ncrv) + KTK_lwmat';
        end
        % function [D_X_out,C_X_out] = compute_Zassenhaus_bases(A1_,A2_)
        function [pck_out1,pck_out2] = compute_Zassenhaus_bases(A1_,A2_)
            [D_X_out,C_X_out] = deal(0);
            [n_A,k_A1,k_A2] = deal(size(A1_,1),size(A1_,2),size(A2_,2));
            [rho_A1,rho_A2,rho_A1A2] = deal(rank(A1_),rank(A2_),rank([A1_,A2_]));
            min_rho_A = min([rho_A1,rho_A2]);
            X = [A1_' A1_'; A2_' (zeros(size(A2_)))'];
            % X_padded = [ X ; zeros( (rho_A1A2+min_rho_A) - (k_A1+k_A2), size(X,2) )];
            X_padded = [ X; A2_' (zeros(size(A2_)))' ];

            % [~,U_X] = lu(X);
            U_X = rref(X);
            % [D_X_out,C_X_out] = rref(X);

            % [~,U_X_padded] = lu(X_padded);
            U_X_padded = rref(X_padded);

            pck_out1 = struct(  'U_X', U_X, ...
                                'C_X', U_X(:, 1:n_A), ...
                                'D_X', U_X(:, (n_A+1):end));
            pck_out2 = struct(  'U_X', U_X_padded, ...
                                'C_X', U_X_padded(:, 1:n_A), ...
                                'D_X', U_X_padded(:, (n_A+1):end));

            % [~,U_X] = lu([A1_' A1_'; A2_' (zeros(size(A2_)))']);
            % U_X(:,1:n_A)
            % size(U_X)
            % C_X_out = U_X(1:rho_A1A2,1:n_A)';
            % D_X_raw = U_X( (rho_A1A2+1):(rho_A1A2+min([rho_A1,rho_A2])) , (n_A+1):end )';
            % D_X_out = D_X_raw(:,1:rank(D_X_raw));
        end
        function [med_pckg_out,det_pckg_out] = post_process_medoid_package(med_pckg_)
            det_pckg_out = LD_aux.compute_medoid_silhouette(med_pckg_);
            [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d,dmat_] = LD_aux.unpack_medoid_package(med_pckg_);

            [kfull,npts] = deal(length(i_meds), length(clst_mem));
            i_full = 1:npts;
            ii_nmeds = logical(prod( i_full ~= i_meds', 1 ));
            i_nmeds = i_full(ii_nmeds);
            iclst_mems = (clst_mem == (1:kfull)');

            get_member_data = @(data_mat_) data_mat_(iclst_mems);

            d_loc_mat = nan(kfull,npts);
            d_loc_mat(iclst_mems) = get_member_data(ones(kfull,npts).*loc_d);
            avg_d_loc = mean(d_loc_mat,2,'omitnan');

            det_pckg_out.avg_cluster_distances = avg_d_loc';
            det_pckg_out.i_full = i_full;
            det_pckg_out.i_nmeds = i_nmeds;

            if (nargout == 1)
                med_pckg_out = combine_structs(med_pckg_,det_pckg_out);
            else
                med_pckg_out = med_pckg_;
            end
        end
        function evl_pckg_out = compute_medoid_silhouette(med_pckg_)
            [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d,dmat_] = LD_aux.unpack_medoid_package(med_pckg_);
            [kfull,npts] = deal(length(i_meds), length(clst_mem));

            i_full = 1:npts;
            ii_nmeds = logical(prod( i_full ~= i_meds', 1 ));

            a_loc_d = loc_d;
            [b_loc_d,s_loc_d] = deal(nan(1,npts));
            for jmed = 1:kfull
                i_meds_nj = [i_meds(1:(jmed-1)),i_meds((jmed+1):end)];
                imems_jclst = i_full(clst_mem == jmed);
                b_loc_d(imems_jclst) = min(dmat_(i_meds_nj,imems_jclst),[],1);
                s_loc_d(imems_jclst) = 1 - a_loc_d(imems_jclst)./b_loc_d(imems_jclst);
            end

            evl_pckg_out = struct(  'silhouette_coeff', mean(s_loc_d,'omitnan'), ...
                                    'silhouette_values', s_loc_d, ...
                                    'b_cluster_distances', b_loc_d);
        end
        % function [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d] = naive_k_medoid(dmat_,kfull_)
        function [med_pckg_out,greedy_med_pckg_out,improved_on_greedy] = naive_k_medoid(dmat_,kfull_)
            greedy_med_pckg_out = LD_aux.greedy_naive_k_medoids(dmat_,kfull_);
            [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d] = LD_aux.unpack_medoid_package(greedy_med_pckg_out);
            net_d_0 = net_d;

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
            improved_on_greedy = net_d < net_d_0;
            [d_mins,imed_assignments] = min( dmat_(i_meds,i_nmedoids) ,[],1);
            n_meds(1:kfull_) = sum(imed_assignments == (1:kfull_)',2);
            net_d = sum(d_mins);
            [clst_mem,loc_d] = deal(zeros(1,npts));
            clst_mem(i_nmedoids) = imed_assignments;
            loc_d(i_nmedoids) = d_mins;
            med2med_d = dmat_(i_meds,i_meds);

            med_pckg_out = LD_aux.pack_medoid_package(i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d,dmat_);
        end
        % function [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d] = greedy_naive_k_medoids(dmat_,kfull_)
        function med_pckg_out = greedy_naive_k_medoids(dmat_,kfull_)
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
            i_nmedoids_greedy = get_not_medoids(i_meds);
            [d_mins_greedy,imed_assignments_greedy] = min( dmat_(i_meds,i_nmedoids_greedy) ,[],1);
            n_meds(1:kfull_) = sum(imed_assignments_greedy == (1:kfull_)',2);
            net_d = sum(d_mins_greedy);
            [clst_mem,loc_d] = deal(zeros(1,npts));
            clst_mem(i_nmedoids_greedy) = imed_assignments_greedy;
            loc_d(i_nmedoids_greedy) = d_mins_greedy;
            med2med_d = dmat_(i_meds,i_meds);

            med_pckg_out = LD_aux.pack_medoid_package(i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d,dmat_);
        end
        function clst_info_out = get_cluster_info(pckg_in_)
            clst_mem = pckg_in_.cluster_membership;
            clst_dst = pckg_in_.cluster_distances;
            kvec = 1:(length(pckg_in_.i_medoids));
            i_full = 1:length(pckg_in_.cluster_membership);

            iclst_mems = clst_mem == kvec';

            clst_info_out = cell(length(kvec),1);
            for ik = kvec
                iclst_mems_k = iclst_mems(ik,:);
                clst_info_out{ik} = [i_full(iclst_mems_k) ; clst_dst(iclst_mems_k)];
            end
        end
        function med_pckg_out = pack_medoid_package(i_meds_,n_meds_,net_d_,clst_mem_,loc_d_,med2med_d_,dmat_)
            med_pckg_out = struct(  'i_medoids', i_meds_, ...
                                    'n_medoids', n_meds_, ...
                                    'net_dist', net_d_, ...
                                    'cluster_membership', clst_mem_, ...
                                    'cluster_distances', loc_d_, ...
                                    'medoid_distances', med2med_d_);
            if (nargin == 7)
                med_pckg_out.distance_matrix = dmat_;
            end
        end
        function [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d,dmat] = unpack_medoid_package(pckg_)
            [i_meds,n_meds,net_d,clst_mem,loc_d,med2med_d] = deal(  pckg_.i_medoids, ...
                                                                    pckg_.n_medoids, ...
                                                                    pckg_.net_dist, ...
                                                                    pckg_.cluster_membership, ...
                                                                    pckg_.cluster_distances, ...
                                                                    pckg_.medoid_distances);
            if (nargout==7)
                dmat = pckg_.distance_matrix;
            end
        end
    end
end
