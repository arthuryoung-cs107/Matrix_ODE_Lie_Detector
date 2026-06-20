classdef mvp_1D_jspc_model
    properties (Constant)

    end

    properties

        Smat;
        npts_per_crv;
        fspace;

        lambdas;

        Rksvds_glb;
        vth_k_glb;
        tau_uk_glb;
        g_tau_dNm1xu_glb;

        Rk_glb_sub_svds;
        vth_k_glb_sub;
        tau_uk_glb_sub;

        Rksvds_crv_cell;
        RNsvds_crv;
        vth_N_crv;
        tau_uN_crv;
        g_tau_uN_crv;

        RN_crv_sub_svds;
        vth_N_crv_sub;
        tau_uN_crv_sub;
    end

    methods (Static)

        function [mods,fspaces] = compute_models(Sobs_,fspace_)
            nvar = length(fspace_.Omap_b);
            ndep = nvar - 1;
            bor_Qjspc = fspace_.bor;

            [Smat,nobs,ncrv,kor_obs,ndim_obs,npts_per_crv,ipts_crv] = ldaux.unpack_Scell(Sobs_,ndep);
            ncrv = length(npts_per_crv);

            xumat = Smat(1:nvar,:);
            xvec = xumat(1,:);
            umat = xumat(2:end,:);
            untns = reshape( Smat(2:end,:), ndep,kor_obs+1,nobs );
            dxuntns = reshape(untns(:,2:end,:),ndep,kor_obs,nobs);

            npts_min = min(npts_per_crv);
            bor_max = sqrt(npts_min/2) - 1; % max C = (1+Q)*P for tall N=Q=1 R matrix
            bor = floor( bor_max ) - 1; % maximally expressive overdetermined N=Q=1 Lambda space
            % bor = min([bor,8]); % 8-10 should be enough
            % bor = 3; % cubic model

            [Plen_1D, Pmat_1D, ~] = ldaux.count_set_P_len(bor,2);
            Omap_b_i = @(i_) [ fspace_.Omap_b(1) ; fspace_.Omap_b(1+i_) ];
            Omap_A_i = @(i_) [ fspace_.Omap_A(1,1) 0 ; 0 fspace_.Omap_A(1+i_,1+i_) ];
            fspace_1D = struct( ...
                'bor', bor, ...
                'Plen', Plen_1D, ...
                'Pmat', Pmat_1D ...
            );
            if (isfield(fspace_,'fam'))
                fspace_1D.fam = fspace_.fam;
            end

            mods = cell([ndep,1]);
            fspaces = cell([ndep,1]);
            for idep = 1:ndep

                s_i_mat = [ ...
                    xvec ; ...
                    reshape(untns(idep,:,:),kor_obs+1,nobs) ...
                ];

                fspace_i = fspace_1D;
                fspace_i.A = Omap_A_i(idep);
                fspace_i.b = Omap_b_i(idep);
                fspace_i = adlam.init_fspace_family(fspace_i);

                mods{idep} = mvp_1D_jspc_model( ...
                    s_i_mat, ...
                    npts_per_crv, ...
                    fspace_i ...
                );

            end

        end

    end

    methods

        function obj = mvp_1D_jspc_model(Smat_,npts_per_crv_,fspace_)
            ndep = 1;
            ncrv = length(npts_per_crv_(:));
            nobs = sum(npts_per_crv_);
            bor = fspace_.bor;
            Pmat = fspace_.Pmat;
            Plen = fspace_.Plen;
            ntheta = 2*Plen;
            csum_npts_per_crv = cumsum(npts_per_crv_);
            ipts_crv = [ [ 1 , (csum_npts_per_crv(1:(end-1))'+1) ] ; csum_npts_per_crv' ];

            xumat = Smat_(1:2,:);
            unmat = Smat_(2:end,:);
            dxunmat = Smat_(3:end,:);
            kor_obs = size(dxunmat,1);

            %% overarching package for all SVD usage
            function Asvd_out = Asvd_package(A_);
                dim_scl = max(size(A_));
                [~,S_,V_] = svd(A_,'econ');
                s_ = reshape(diag(S_),[],1); % output s_ is col vector
                r_ = sum(double(s_ > dim_scl*eps(S_(1)))); % default matlab tol

                Asvd_out = struct( ...
                    'dim', length(s_), ...
                    'r', r_, ...
                    's', s_, ...
                    'V', V_, ...
                    'D', ( s_(:)  )' .* V_, ... % rowspace (domain) Y approx. basis
                    'W', ( s_(end) ./ s_(:)  )' .* V_ ... % nullspace (kernal) K approx. basis
                );
            end
            normalize_Renc = @(r_) r_/norm(r_); % unit length Rmat rows
            % normalize_Renc = @(r_) r_; % unnormalized Rmat rows

            lambdas = cell([nobs,1]);
            lvs = nan(Plen,nobs);
            dkxls = nan(kor_obs,Plen,nobs);
            lkxs = nan(Plen,kor_obs,nobs);
            Renc_cell = cell(kor_obs,nobs);
            LamN_T_tns = nan(ntheta,1 + (kor_obs+1), nobs);
            for iobs = 1:nobs

                lambdas{iobs} = adlam(fspace_, xumat(1:2,iobs), Smat_(3:end,iobs) );
                l_i = lambdas{iobs};
                lv_i = l_i.lrow_vals;

                lvs(:,iobs) = lv_i;
                dkxls(:,:,iobs) = l_i.dkxl;
                lkxs(:,:,iobs) = reshape(l_i.lkx(1,:,:),Plen,kor_obs);

                LamN_T_tns(:,1:2,iobs) = [ [lv_i(:) ; zeros(Plen,1)] , [zeros(Plen,1) ; lv_i(:)] ];
                Renc_cell{1,iobs} = normalize_Renc([ -l_i.dxu(1,1)*lv_i , lv_i ]);
                iLam_km1 = 3;
                LamN_T_tns(:,iLam_km1,iobs) = [ -l_i.lkx(1,:,1) , l_i.dkxl(1,:) ]';
                for k = 2:kor_obs
                    Renc_cell{k,iobs} = normalize_Renc( ...
                    [ -LamN_T_tns(1:Plen,iLam_km1,iobs)'-l_i.dxu(1,k)*lv_i , LamN_T_tns((Plen+1):end,iLam_km1,iobs)' ] ...
                    );
                    iLam_km1 = iLam_km1 + 1;
                    LamN_T_tns(:,iLam_km1,iobs) = [ -l_i.lkx(1,:,k) , l_i.dkxl(k,:) ]';
                end
            end
            % debug_script
            Lam_dkxu_tns = reshape(LamN_T_tns(:,3:end,:),ntheta,kor_obs,nobs);
            LamN_tns = permute(LamN_T_tns,[2 1 3]);

            %% prepare local utility functions
            inds_icrv = @(icrv_) ipts_crv(1,icrv_):ipts_crv(2,icrv_);
            vertcat_cell = @(c_) vertcat(c_{:});
            horzcat_cell = @(c_) horzcat(c_{:});
            shape_Rkm_glb =  @(RkcT_) reshape( horzcat_cell(RkcT_),ntheta,nobs)';
            Renc = struct( ...
            'Rkc_2_Rkm_glb', @(Rkc_) shape_Rkm_glb( cellfun( @(R_) R_',Rkc_,'UniformOutput', false) ), ...
            'Rk_icrv', @(k_,icrv_) vertcat_cell(reshape(Renc_cell(k_,inds_icrv(icrv_)),[],1))  ...
            );
            function vth_out = comp_vartheta(W_,lvs_)
                vx_W = (W_(1:(size(W_,1)/2),:))' * lvs_;
                Deltas_W = sum(vx_W.*vx_W,1);
                vth_out = W_ * (vx_W./Deltas_W);
            end
            function tuk_out = comp_tau_uk_glb(th_,k_,iP_)
                Plen_b = size(th_,1)/2;
                tuk_out = nan(nobs,1);
                if (k_==1)
                    tuk_out = sum(th_((Plen_b+1):end,:) .* lvs(iP_,:), 1)';
                else
                    for i = 1:nobs
                        tuk_out(i) =  ...
                        dkxls(k_-1,iP_,i) * th_((Plen_b+1):end,i)  ...
                         - ( lkxs(iP_,k_-1,i)' * th_(1:Plen_b,i) );
                    end
                end
            end
            function tuN_out = comp_tau_uN_crv(th_,inds_,iP_)
                Plen_b = size(th_,1)/2;
                nobs_crv = size(th_,2);
                tuN_out = nan(kor_obs+1,nobs_crv);
                for i = 1:nobs_crv
                    ii = inds_(i);
                    tuN_out(1,i) = th_((Plen_b+1):end,i)' * lvs(iP_,ii);
                    for k = 2:(kor_obs+1)
                        tuN_out(k,i) = ...
                        dkxls(k-1,iP_,ii) * th_((Plen_b+1):end,i)  ...
                         - ( lkxs(iP_,k-1,ii)' * th_(1:Plen_b,i) );
                    end
                end
            end
            shape_pPbs = @(pPbs_) reshape(double(pPbs_(:)).*ones( length(pPbs_),ndep+1 ), [],1)';
            is_prjtheta_b_sub = @(b_) logical(shape_pPbs(mvp_jspc_model.is_prjP_b_sub(Pmat,b_)));
            % prepare b'th suborder model indexing
            is_pPbs_mat = zeros(bor,Plen);
            is_pTbs_mat = zeros(bor,ntheta);
            for b = 1:bor
                is_pPbs_mat(b,:) = mvp_jspc_model.is_prjP_b_sub(Pmat,b);
                is_pTbs_mat(b,:) = is_prjtheta_b_sub(b);
            end
            is_pPbs_mat = logical(is_pPbs_mat);
            is_pTbs_mat = logical(is_pTbs_mat);
            bor_to_one = flip(1:bor);

            %% compute global R models of each jetspace order
            Rksvds_glb = cell([1,kor_obs]);
            vth_k_glb = nan(ntheta,nobs,kor_obs);
            tau_uk_glb = nan(nobs,kor_obs+1);
            Rk_glb_sub_svds = cell([bor,kor_obs]);
            vth_k_glb_sub = cell([bor,kor_obs]);
            tau_uk_glb_sub = nan(nobs,bor,kor_obs+1);
            for k = 1:kor_obs
                Rkm_glb = Renc.Rkc_2_Rkm_glb(Renc_cell(k,:));
                Rksvds_glb{k} = Asvd_package(Rkm_glb);
                vth_k_glb(:,:,k) = comp_vartheta( ...
                    Rksvds_glb{k}.W, lvs ...
                );
                % estimate of k'th derivative
                tau_uk_glb(:,k) = comp_tau_uk_glb(vth_k_glb(:,:,k),k,1:Plen);
                if (k == kor_obs)
                    % finish by estimating N+1 derivative and restricted N'th derivative
                    tau_uk_glb(:,k+1) = comp_tau_uk_glb(vth_k_glb(:,:,k),k+1,1:Plen);
                end
                % execute the same procedure on ordered subspaces of P permutation model
                for b = bor_to_one
                    Rk_glb_sub_svds{b,k} = Asvd_package( ...
                        Rkm_glb( :, is_pTbs_mat(b,:) ) ...
                    );
                    vth_k_glb_sub{b,k} = comp_vartheta(Rk_glb_sub_svds{b,k}.W, ...
                        lvs(is_pPbs_mat(b,:),:) ...
                    );
                    tau_uk_glb_sub(:,b,k) = comp_tau_uk_glb(vth_k_glb_sub{b,k},k,is_pPbs_mat(b,:));
                    if (k == kor_obs)
                        tau_uk_glb_sub(:,b,k+1) = comp_tau_uk_glb(vth_k_glb_sub{b,k},k+1,is_pPbs_mat(b,:));
                    end
                end
            end
            tau_uk_glb = tau_uk_glb'; % --> kor_obs+1,nobs
            tau_uk_glb_sub = permute(tau_uk_glb_sub, [3 1 2] ); % --> kor_obs+1,nobs,bor

            g_tau_dNm1xu_glb = nan(1+(kor_obs+1),nobs);
            for iobs = 1:nobs
                g_tau_dNm1xu_glb(:,iobs) = lambdas{iobs}.J_tau_dNm1xu( ...
                    vth_k_glb(:,iobs,kor_obs), ...
                    tau_uk_glb(kor_obs,iobs) ...
                );
            end

            %% compute local (curve) models
            Rksvds_crv_cell = cell([kor_obs,ncrv]);
            RNsvds_crv = cell([1,ncrv]);
            vth_N_crv = nan(ntheta,nobs);
            tau_uN_crv = nan(kor_obs+1,nobs);
            g_tau_uN = nan(1+(kor_obs+1),kor_obs+1,nobs);
            RN_crv_sub_svds = cell([bor,ncrv]);
            vth_N_crv_sub = cell([bor,ncrv]);
            tau_uN_crv_sub = nan(kor_obs+1,nobs,bor);
            for icrv = 1:ncrv
                inds_i = ipts_crv(1,icrv):ipts_crv(2,icrv);
                Rkenc_cell_i = Renc_cell(:, inds_i);
                RNmat_agg_i = vertcat( Rkenc_cell_i{:} ); % true aggregate RN matrix
                RNtns_agg_i = reshape(RNmat_agg_i',ntheta,kor_obs,[]);
                for k = 1:kor_obs
                    Rksvds_crv_cell{k,icrv} = Asvd_package( reshape(RNtns_agg_i(:,k,:),ntheta,[])' );
                end
                RN_svd_i_net = reshape(Rksvds_crv_cell(:,icrv),[],1);
                RN_mat_i_net = horzcat( vertcat(RN_svd_i_net{:}).D )';
                RN_mat_i = RN_mat_i_net;

                % compute concatenated curve model, all information together
                RNsvds_crv{icrv} = Asvd_package(RN_mat_i); % compute P permutation (full) local curve model
                % evaluate local parameters predicted by local curve model
                vth_N_crv(:,inds_i) = comp_vartheta( ...
                    RNsvds_crv{icrv}.W,lvs(:,inds_i) ...
                );
                % estimate all curve derivatives
                % tau_uN_crv(:,inds_i) = comp_tau_uN_crv( ...
                %     vth_N_crv(:,inds_i), inds_i, 1:Plen ...
                % );
                for iiobs = inds_i
                    [g_tau_uN_crv(:,:,iiobs),tau_uN_crv(:,iiobs)] = ...
                        lambdas{iiobs}.J_tau_uN(vth_N_crv(:,iiobs));
                end

                % execute the same procedure on ordered subspaces of P permutation model
                for b = bor_to_one
                    RN_crv_sub_svds{b,icrv} = Asvd_package( ...
                        RN_mat_i( :, is_pTbs_mat(b,:) ) ...
                    );
                    vth_N_crv_sub{b,icrv} = comp_vartheta( ...
                        RN_crv_sub_svds{b,icrv}.W, lvs(is_pPbs_mat(b,:),inds_i) ...
                    );
                    tau_uN_crv_sub(:,inds_i,b) = comp_tau_uN_crv( ...
                        vth_N_crv_sub{b,icrv},inds_i,is_pPbs_mat(b,:) ...
                    );
                end
            end

            %% core assignments
            obj.Smat = Smat_;
            obj.npts_per_crv = npts_per_crv_;
            obj.fspace = fspace_;

            obj.lambdas = lambdas;

            obj.Rksvds_glb = Rksvds_glb;
            obj.vth_k_glb = vth_k_glb;
            obj.tau_uk_glb = tau_uk_glb;
            obj.g_tau_dNm1xu_glb = g_tau_dNm1xu_glb;

            obj.Rk_glb_sub_svds = Rk_glb_sub_svds;
            obj.vth_k_glb_sub = vth_k_glb_sub;
            obj.tau_uk_glb_sub = tau_uk_glb_sub;

            obj.Rksvds_crv_cell = Rksvds_crv_cell;
            obj.RNsvds_crv = RNsvds_crv;
            obj.vth_N_crv = vth_N_crv;
            obj.tau_uN_crv = tau_uN_crv;
            obj.g_tau_uN_crv = g_tau_uN_crv;

            obj.RN_crv_sub_svds = RN_crv_sub_svds;
            obj.vth_N_crv_sub = vth_N_crv_sub;
            obj.tau_uN_crv_sub = tau_uN_crv_sub;

        end

    end

end
