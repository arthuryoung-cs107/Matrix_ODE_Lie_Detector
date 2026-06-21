classdef mvp_jspc_model
    properties (Constant)

        is_prjP_b_sub = @(P_,b_) sum(P_,1) <= b_;
        pick_prjP_b_sub = @(is_,ifull_) ifull_(is_);
        inds_prjP_b_sub = @(P_,b_) mvp_jspc_model.pick_prjP_b_sub(sum(P_,1) <= b_, 1:size(P_,2));
        prjP_b_sub = @(P_,b_) P_( : , sum(P_,1) <= b_ );

        vertcat_glb_SVDs = @(o_) vertcat(o_.Rqksvd_glb_cell(:),o_.Rksvds_glb(:),o_.Rk_glb_sub_svds(:));

    end

    properties

        Sobs;
        Sdat;
        fdat;

        fspace;

        jspc_1Dmods;
        jspc_N1mod;

        Smat;
        npts_per_crv;
        ipts_crv;

        lambdas;
        lvs;
        LamN_tns;
        Renc_cell;

        Gsvd;

        lvs_svd;
        LamN_svd;

        Lam_dkxuq_svd_glb_cell;
        Lam_dkxu_svd_glb;

        Rqksvd_glb_cell;
        vth_qk_glb;
        tau_uqk_glb;

        Rksvds_glb;
        vth_k_glb;
        tau_uk_glb;

        Rk_glb_sub_svds;
        vth_k_glb_sub;
        tau_uk_glb_sub;

        RN_YLdNm1xu_glb_svd;
        vth_RN_YLdNm1xu_glb;
        tau_dNm1xu_RN_YLdNm1xu_glb;
        vthx_N_glb_svd;

        Rqksvds_crv_cell;
        RNsvds_crv;
        vth_N_crv;
        tau_uN_crv;

        RNsvds_crv_Yvthxglb;
        vth_N_crv_Yvthxglb;
        tau_uN_crv_Yvthxglb;

        RN_crv_sub_svds;
        vth_N_crv_sub;
        tau_uN_crv_sub;

    end

    methods
        function obj = mvp_jspc_model(Sobs_,Sdat_,fdat_)
            [ndep,nvar] = deal(Sdat_.ndep, Sdat_.ndep+1);
            bor = fdat_.bor;
            if ( ~isfield(fdat_,'Omap_A') )
                if (isfield(fdat_,'Omap_a'))
                    fdat_.Omap_A = fdat_.Omap_a(:) .* eye(ndep+1);
                else
                    fdat_.Omap_A = eye(ndep+1);
                end
            end
            if ( ~isfield(fdat_,'Omap_b') )
                Omap_b = zeros(1+ndep,1);
            end
            if ( ~isfield(fdat_,'Pmat') )
                [Plen, Pmat, ~] = ldaux.count_set_P_len(bor,ndep+1);
                fdat_.Pmat = Pmat;
            else
                Pmat = fdat_.Pmat;
                Plen = size(fdat_.Pmat,2);
            end

            %% overarching package for all SVD usage
            function Asvd_out = Asvd_package(A_);
                % [r_,s_,V_,U_] = mvp_jspc_model.rsV_unpack(Ri_);
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
                % 'A', @(o_) o_.U * ( o_.s(:) .* (o_.V') ), ...
                % 'U', U_, ...
            end
            ntheta = nvar*Plen;
            %% initialize injection into Lambda column space
            i_imm = zeros(ndep*Plen,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            function l_imm = immerse_lambda(l_)
                l_imm = zeros((ndep*Plen)*ndep,1);
                l_imm(logical( i_imm(:) )) = reshape( l_(:) * ones(1,ndep), [], 1);
                l_imm = (reshape(l_imm,ndep*Plen,ndep))';
                % l_ -> [ l_ 0 ... 0 ; 0 l_ ... 0 ; ... ; 0 0 ... l_ ]
            end
            fspace = fdat_;
            fspace.Plen = Plen;
            fspace.ntheta = ntheta;
            fspace.imm_l = @(lambda_) immerse_lambda(lambda_);
            % fspace.fam = 'Hermite1';
            % fspace.fam = 'Hermite2';
            % fspace.fam = 'Legendre';
            % fspace.fam = 'Chebyshev1';
            % fspace.fam = 'Chebyshev2';
            fspace = adlam.init_fspace_family(fspace);
            % fspace

            jspc_1Dmods = mvp_1D_jspc_model.compute_models(Sobs_,fspace);
            jspc_N1mod = mvp_N1_jspc_model(Sobs_,fspace);
            JF_N1mod = jspc_N1mod.JF_glb;

            [Smat,nobs,nset,kor_obs,ndim_obs,npts_per_crv,ipts_crv] = ldaux.unpack_Scell(Sobs_,ndep);
            ncrv = length(npts_per_crv);

            xumat = Smat(1:nvar,:);
            untns = reshape( Smat(2:end,:), ndep,kor_obs+1,nobs );
            dxuntns = reshape(untns(:,2:end,:),ndep,kor_obs,nobs);

            %% Process observations point by point. Compute prolongation quantities of function space
            lambdas = adlam.empty(nobs,0);
            lvs = nan(Plen,nobs);
            LamN_T_tns = nan(ntheta,1 + ndep*(kor_obs+1), nobs);
            Renc_cell = cell(kor_obs,nobs);
            Genc_tns = nan(ndep,ntheta,nobs);
            tic_prN = tic;
            for iobs = 1:nobs

                lambdas(iobs) = adlam(fdat_, xumat(:,iobs), Smat((nvar+1):end,iobs) );

                l_i = lambdas(iobs);
                % lv_i = adlam.lrowP(l_i);
                lv_i = l_i.lrow_vals;

                lvs(:,iobs) = lv_i;

                l_imm_i = fspace.imm_l(lv_i);

                LamN_T_tns(:,1:nvar,iobs) = [ [lv_i(:) ; zeros(ntheta-Plen,1) ] , [ zeros(Plen,ndep) ; l_imm_i' ] ];
                Renc_cell{1,iobs} = [ -l_i.dxu(:,1)*lv_i , l_imm_i ];
                iLam_km1 = (nvar+1):(nvar+ndep);
                LamN_T_tns(:,iLam_km1,iobs) = [ -l_i.lkx(:,:,1) , fspace.imm_l(l_i.dkxl(1,:)) ]';
                for k = 2:kor_obs
                    % Renc_cell{k,iobs} = [ -l_i.lkx(:,:,k-1) - l_i.dxu(:,k)*lv_i , fspace.imm_l(l_i.dkxl(k-1,:)) ];
                    Renc_cell{k,iobs} = [ ...
                    -LamN_T_tns(1:Plen,iLam_km1,iobs)'-l_i.dxu(:,k)*lv_i , LamN_T_tns((Plen+1):end,iLam_km1,iobs)' ...
                    ];
                    iLam_km1 = iLam_km1 + ndep;
                    LamN_T_tns(:,iLam_km1,iobs) = [ -l_i.lkx(:,:,k) , fspace.imm_l(l_i.dkxl(k,:)) ]';
                end
                Genc_tns(:,:,iobs) = JF_N1mod(:,:,iobs)*LamN_T_tns(:,:,iobs)';
            end
            toc_prN = toc(tic_prN);
            fprintf('(mvp_jspc_model) Prolonged %d observations over Q=%d, N=%d (B=%d) jet space with order O=%d mvpolynomials (C=%d) in %.3f seconds \n', ...
            nobs, ndep, kor_obs, ndim_obs, fspace.bor, ntheta, ...
            toc_prN);

            Gmat = (reshape(permute(Genc_tns,[2 1 3]),ntheta,ndep*nobs))';
            Gsvd = Asvd_package(Gmat);

            % debug_script
            Lam_dkxu_ttns = reshape(LamN_T_tns(:,(nvar+1):end,:),ntheta,ndep,kor_obs,nobs);
            LamN_tns = permute(LamN_T_tns,[2 1 3]); % --> ndim, ntheta, nobs

            %% prepare local utility functions
            inds_icrv = @(icrv_) ipts_crv(1,icrv_):ipts_crv(2,icrv_);
            vertcat_cell = @(c_) vertcat(c_{:});
            horzcat_cell = @(c_) horzcat(c_{:});
            shape_Rkt_glb =  @(RkcT_) permute(reshape( horzcat_cell(RkcT_),ntheta,ndep,nobs),[3 1 2]);
            Renc = struct( ...
            'Rkc_2_Rkt_glb', @(Rkc_) shape_Rkt_glb( cellfun( @(R_) R_',Rkc_,'UniformOutput', false) ), ...
            'Rk_icrv', @(k_,icrv_) vertcat_cell(reshape(Renc_cell(k_,inds_icrv(icrv_)),[],1))  ...
            );
            function vth_out = comp_vartheta(W_,lvs_)
                vx_W = (W_(1:(size(W_,1)/nvar),:))' * lvs_;
                Deltas_W = sum(vx_W.*vx_W,1);
                vth_out = W_ * (vx_W./Deltas_W);
            end
            function tuk_out = comp_tau_uk_glb(th_,k_,iP_)
                Plen_b = size(th_,1)/(ndep+1);
                tuk_out = nan(ndep,nobs);
                if (k_==1)
                    for i = 1:nobs
                        tuk_out(:,i) =  reshape(th_((Plen_b+1):end,i),Plen_b,ndep)' * lvs(iP_,i);
                    end
                else
                    for i = 1:nobs
                        tuk_out(:,i) =  ...
                        ( lambdas(i).dkxl(k_-1,iP_) * reshape(th_((Plen_b+1):end,i),Plen_b,ndep) )' ...
                         - ( lambdas(i).lkx(:,iP_,k_-1)*th_(1:Plen_b,i) );
                    end
                end
            end
            function tuqk_out = comp_tau_uqk_glb(th_,q_,k_,iP_)
                Plen_b = size(th_,1)/(ndep+1);
                th_u = reshape(th_((Plen_b+1):end,:), Plen_b,ndep,nobs );
                th_uq = reshape(th_u(:,q_,:),Plen_b,nobs);
                if (k_==1)
                    tuqk_out = sum(th_uq.*lvs(iP_,:),1);
                else
                    tuqk_out = nan(nobs,1);
                    for i = 1:nobs
                        tuqk_out(i) =  ...
                            lambdas(i).dkxl(k_-1,iP_) * th_uq(:,i) ...
                             - ( lambdas(i).lkx(q_,iP_,k_-1)*th_(1:Plen_b,i) );
                    end
                end
            end
            function tuN_out = comp_tau_uN_crv(th_,inds_,iP_)
                Plen_b = size(th_,1)/(ndep+1);
                nobs_crv = size(th_,2);
                tuN_out = nan(ndep,kor_obs+1,nobs_crv);
                for i = 1:nobs_crv
                    ii = inds_(i);
                    th_u_i_mat = reshape(th_((Plen_b+1):end,i),Plen_b,ndep);
                    tuN_out(:,1,i) = th_u_i_mat' * lvs(iP_,ii);
                    for k = 2:(kor_obs+1)
                        tuN_out(:,k,i) = ...
                        ( lambdas(ii).dkxl(k-1,iP_) * th_u_i_mat )' ...
                         - ( lambdas(ii).lkx(:,iP_,k-1) * th_(1:Plen_b,i) );
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

            %% use utility functions and decomposed R and Lambda matrices to model data

            lvs_svd = Asvd_package(lvs');
            LamN_svd = Asvd_package( reshape(LamN_tns,ntheta,[])' );

            %% compute global models
            Lam_dkxuq_svd_glb_cell = cell([ndep,kor_obs]);
            Lam_dkxu_svd_glb = cell([1,kor_obs]);

            Rqksvd_glb_cell = cell([ndep,kor_obs]);
            vth_qk_glb = nan(ntheta,nobs,ndep,kor_obs);
            tau_uqk_glb = nan(nobs,ndep,kor_obs);

            Rksvds_glb = cell([1,kor_obs]);
            vth_k_glb = nan(ntheta,nobs,kor_obs);
            tau_uk_glb = nan(ndep,nobs,kor_obs+1);

            Rk_glb_sub_svds = cell([bor,kor_obs]);
            vth_k_glb_sub = cell([bor,kor_obs]);
            tau_uk_glb_sub = nan(ndep,nobs,bor,kor_obs+1);
            tic_glb = tic;
            for k = 1:kor_obs
                Rkt_glb = Renc.Rkc_2_Rkt_glb(Renc_cell(k,:));
                [Rk_glb_net,Lam_dkxu_net] = deal(zeros(ntheta,ntheta,ndep));
                for idep = 1:ndep
                    Lam_dkxuq_svd_glb_cell{idep,k} = Asvd_package( reshape(Lam_dkxu_ttns(:,idep,k,:),ntheta,[])' );
                    Lam_dkxu_net(:,:,idep) = Lam_dkxuq_svd_glb_cell{idep,k}.D;

                    Rqksvd_glb_cell{idep,k} = Asvd_package( Rkt_glb(:,:,idep) );
                    Rk_glb_net(:,:,idep) = Rqksvd_glb_cell{idep,k}.D;
                    vth_qk_glb(:,:,idep,k) = comp_vartheta( ...
                        Rqksvd_glb_cell{idep,k}.W, lvs ...
                    );
                    tau_uqk_glb(:,idep,k) = comp_tau_uqk_glb(vth_qk_glb(:,:,idep,k),idep,k,1:Plen);
                end
                Lam_dkxu_svd_glb{k} = Asvd_package( reshape(Rk_glb_net,ntheta,ntheta*ndep)' );
                Rk_glb_net = reshape(Rk_glb_net,ntheta,[])';
                Rksvds_glb{k} = Asvd_package(Rk_glb_net);
                vth_k_glb(:,:,k) = comp_vartheta( ...
                    Rksvds_glb{k}.W, lvs ...
                );
                % estimate of k'th derivative, global model
                tau_uk_glb(:,:,k) = comp_tau_uk_glb(vth_k_glb(:,:,k),k,1:Plen);
                if (k == kor_obs)
                    % finish by estimating N+1 derivative and restricted N'th derivative
                    tau_uk_glb(:,:,k+1) = comp_tau_uk_glb(vth_k_glb(:,:,k),k+1,1:Plen);
                end
                % execute the same procedure on ordered subspaces of P permutation model
                for b = bor_to_one
                    Rk_glb_sub_svds{b,k} = Asvd_package( ...
                        Rk_glb_net( :, is_pTbs_mat(b,:) ) ...
                    );
                    vth_k_glb_sub{b,k} = comp_vartheta(Rk_glb_sub_svds{b,k}.W, ...
                        lvs(is_pPbs_mat(b,:),:) ...
                    );
                    tau_uk_glb_sub(:,:,b,k) = comp_tau_uk_glb(vth_k_glb_sub{b,k},k,is_pPbs_mat(b,:));
                    if (k == kor_obs)
                        tau_uk_glb_sub(:,:,b,k+1) = comp_tau_uk_glb(vth_k_glb_sub{b,k},k+1,is_pPbs_mat(b,:));
                    end
                end
            end
            toc_glb = toc(tic_glb);
            fprintf('(mvp_jspc_model) computed global svds in %.3f seconds \n', toc_glb);

            tau_uqk_glb = permute(tau_uqk_glb,[2 3 1]); % --> ndep,kor_obs,nobs
            tau_uk_glb = permute(tau_uk_glb, [1 3 2] ); % --> ndep,kor_obs+1,nobs
            tau_uk_glb_sub = permute(tau_uk_glb_sub, [1 4 2 3] ); % --> ndep,kor_obs+1,nobs,bor

            vthx_N_glb = reshape(vth_k_glb(1:Plen,:,end),Plen,nobs)';
            vthx_N_glb_svd = Asvd_package(vthx_N_glb);
            dim_Yvthx = min([ vthx_N_glb_svd.r max([ (min(npts_per_crv)*ndep*kor_obs + Plen - ntheta - 1) 2 ]) ]);
            iYvthx = 1:dim_Yvthx;
            % iYvthx = 1:(vthx_N_glb_svd.r);
            % iYvthx = 1:(vthx_N_glb_svd.r - 1);
            % iYvthx = 1:(Plen-1);
            % iYvthx = 1:16;
            C_Yvthx = length(iYvthx);
            Yvthx_N_glb = vthx_N_glb_svd.V(:,iYvthx); % restricted orthonormal basis for curve R^(N) kernels, induced by global data

            % compute local (curve) models
            Rqksvds_crv_cell = cell([ndep,kor_obs,ncrv]);
            RNsvds_crv = cell([1,ncrv]);
            vth_N_crv = nan(ntheta,nobs);
            tau_uN_crv = nan(ndep,kor_obs+1,nobs);
            RNsvds_crv_Yvthxglb = cell([1,ncrv]);
            vth_N_crv_Yvthxglb = nan(ntheta,nobs);
            tau_uN_crv_Yvthxglb = nan(ndep,kor_obs+1,nobs);
            RN_crv_sub_svds = cell([bor,ncrv]);
            vth_N_crv_sub = cell([bor,ncrv]);
            tau_uN_crv_sub = nan(ndep,kor_obs+1,nobs,bor);
            tic_crv = tic;
            for icrv = 1:ncrv
                inds_i = ipts_crv(1,icrv):ipts_crv(2,icrv);
                Rkenc_cell_i = Renc_cell(:, inds_i);
                RNmat_agg_i = vertcat( Rkenc_cell_i{:} ); % true aggregate RN matrix
                RNttns_agg_i = reshape(RNmat_agg_i',ntheta,ndep,kor_obs,[]);
                for k = 1:kor_obs
                    for idep = 1:ndep
                        Rqksvds_crv_cell{idep,k,icrv} = Asvd_package( reshape(RNttns_agg_i(:,idep,k,:),ntheta,[])' );
                    end
                end

                % RN_svd_i_net = Rqksvds_crv_cell(:,:,icrv)
                RN_svd_i_net = reshape(Rqksvds_crv_cell(:,:,icrv),[],1);
                RN_mat_i_net = horzcat( vertcat(RN_svd_i_net{:}).D )';
                % RN_net_svd_i = Asvd_package(RN_mat_i_net)
                % RN_agg_svd_i = Asvd_package( RNmat_agg_i ) % SVD of true aggregate RN matrix

                % RN_svdarray_i_net(1).D
                % RN_mat_i = RNmat_agg_i;
                RN_mat_i = RN_mat_i_net;
                % size(RN_mat_i)
                % pause

                % compute concatenated curve model, all information together
                RNsvds_crv{icrv} = Asvd_package(RN_mat_i); % compute P permutation (full) local curve model
                % evaluate local parameters predicted by local curve model
                vth_N_crv(:,inds_i) = comp_vartheta( ...
                    RNsvds_crv{icrv}.W,lvs(:,inds_i) ...
                );
                % estimate all curve derivatives
                tau_uN_crv(:,:,inds_i) = comp_tau_uN_crv( ...
                    vth_N_crv(:,inds_i), inds_i, 1:Plen ...
                );

                % compute restricted SVD to span of global vartheta set (by R^N kernel)
                RNsvds_crv_Yvthxglb{icrv} = Asvd_package([ RN_mat_i(:,1:Plen)*Yvthx_N_glb , RN_mat_i(:,(Plen+1):end) ]);
                % evaluate local parameters predicted by local restricted curve model, embedded back in full P permutation space
                vth_N_crv_Yvthxglb(:,inds_i) = comp_vartheta( ...
                    [ zeros(ntheta,Plen-C_Yvthx) , ...
                        [Yvthx_N_glb*RNsvds_crv_Yvthxglb{icrv}.W(1:C_Yvthx,:) ; RNsvds_crv_Yvthxglb{icrv}.W((C_Yvthx+1):end,:)] ] , ...
                    lvs(:,inds_i) ...
                );
                % estimate all curve derivatives predicted by restricted curve model
                tau_uN_crv_Yvthxglb(:,:,inds_i) = comp_tau_uN_crv( ...
                    vth_N_crv_Yvthxglb(:,inds_i), inds_i, 1:Plen ...
                );

                % execute the same procedure on ordered subspaces of P permutation model
                for b = bor_to_one
                    RN_crv_sub_svds{b,icrv} = Asvd_package( ...
                        RN_mat_i( :, is_pTbs_mat(b,:) ) ...
                    );
                    vth_N_crv_sub{b,icrv} = comp_vartheta( ...
                        RN_crv_sub_svds{b,icrv}.W, lvs(is_pPbs_mat(b,:),inds_i) ...
                    );
                    tau_uN_crv_sub(:,:,inds_i,b) = comp_tau_uN_crv( ...
                        vth_N_crv_sub{b,icrv},inds_i,is_pPbs_mat(b,:) ...
                    );
                end
            end
            toc_crv = toc(tic_crv);
            fprintf('(mvp_jspc_model) computed local (curve) svds in %.3f seconds \n', toc_crv);

            % fspc.print_vshort_polynomial_theta(Rksvds_glb{1}.W(:,end),Pmat)
            % fspc.print_vshort_polynomial_theta(Rksvds_glb{1}.W(:,end-1),Pmat)

            %% assignments

            obj.Sobs = Sobs_;
            obj.Sdat = Sdat_;
            obj.fdat =  fdat_;

            obj.fspace = fspace;

            obj.Smat = Smat;

            obj.npts_per_crv = npts_per_crv;
            obj.ipts_crv = ipts_crv;

            obj.jspc_1Dmods = jspc_1Dmods;
            obj.jspc_N1mod = jspc_N1mod;

            obj.lambdas = lambdas;
            obj.lvs = lvs;
            obj.LamN_tns = LamN_tns;
            obj.Renc_cell = Renc_cell;

            obj.Gsvd = Gsvd;

            obj.lvs_svd = lvs_svd;
            obj.LamN_svd = LamN_svd;
            obj.Lam_dkxuq_svd_glb_cell = Lam_dkxuq_svd_glb_cell;
            obj.Lam_dkxu_svd_glb = Lam_dkxu_svd_glb;

            obj.Rqksvd_glb_cell = Rqksvd_glb_cell;
            obj.vth_qk_glb = vth_qk_glb;
            obj.tau_uqk_glb = tau_uqk_glb;

            obj.Rksvds_glb = Rksvds_glb;
            obj.vth_k_glb = vth_k_glb;
            obj.tau_uk_glb = tau_uk_glb;

            obj.Rk_glb_sub_svds = Rk_glb_sub_svds;
            obj.vth_k_glb_sub = vth_k_glb_sub;
            obj.tau_uk_glb_sub = tau_uk_glb_sub;

            % obj.RN_YLdNm1xu_glb_svd = RN_YLdNm1xu_glb_svd;
            % obj.vth_RN_YLdNm1xu_glb = vth_RN_YLdNm1xu_glb;
            % obj.tau_dNm1xu_RN_YLdNm1xu_glb = tau_dNm1xu_RN_YLdNm1xu_glb;

            obj.vthx_N_glb_svd = vthx_N_glb_svd;

            obj.Rqksvds_crv_cell = Rqksvds_crv_cell;
            obj.RNsvds_crv = RNsvds_crv;
            obj.RNsvds_crv_Yvthxglb = RNsvds_crv_Yvthxglb;
            obj.RN_crv_sub_svds = RN_crv_sub_svds;
            obj.vth_N_crv = vth_N_crv;
            obj.tau_uN_crv = tau_uN_crv;
            obj.vth_N_crv_Yvthxglb = vth_N_crv_Yvthxglb;
            obj.tau_uN_crv_Yvthxglb = tau_uN_crv_Yvthxglb;
            obj.vth_N_crv_sub = vth_N_crv_sub;
            obj.tau_uN_crv_sub = tau_uN_crv_sub;

        end
    end

    methods (Static)

        function obj_out = verify(obj)

            ndep = obj.Sdat.ndep;
            nvar = ndep + 1;
            [ndim_obs,nobs] = size(obj.Smat);
            kor_obs = (ndim_obs-1)/ndep - 1;

            mods_1D = obj.jspc_1Dmods;

            Smat = obj.Smat;
            fspace = obj.fspace;
            lvs = obj.lvs;

            xumat = Smat(1:nvar,:);
            untns = reshape( Smat(2:end,:), ndep,kor_obs+1,nobs );
            dxuntns = reshape(untns(:,2:end,:),ndep,kor_obs,nobs);

            function err_stats(prefix_,err_)
                err_tol = 1e-3;

                err_inc = sort(err_(:));
                err_min = err_inc(1);
                err_med = median(err_inc);
                err_avg = mean(err_inc);
                err_max = err_inc(end);
fprintf( '(%s err) [min,med,avg,max]=[%.1e,%.1e,%.1e,%.1e]. Success: [med,max] = [%d %d] \n', ...
            prefix_, ...
            err_min,err_med,err_avg,err_max, ...
            err_med < err_tol, err_max < err_tol ...
            );
                err_mat_min = min(err_,[],3);
                err_mat_med = median(err_,3);
                err_mat_avg = mean(err_,3);
                err_mat_max = max(err_,[],3);
                for q = 1:ndep
                    fprintf('   (q=%d [min,med,avg,max]) ',q);
                    for k = 1:kor_obs
                        fprintf('   [%.2e %.2e %.2e %.2e]', ...
                            err_mat_min(q,k), err_mat_med(q,k), err_mat_avg(q,k), err_mat_max(q,k)  );
                    end
                    fprintf('\n');
                end
            end
            comp_err_msr = @(u_) abs((u_-dxuntns)./dxuntns);

            tau_uqk_glb = obj.tau_uqk_glb;
            tau_uk_glb = obj.tau_uk_glb(:,1:(end-1),:);
            tau_uN_crv = obj.tau_uN_crv(:,1:(end-1),:);
            tau_uN_crv_Yvthxglb = obj.tau_uN_crv_Yvthxglb(:,1:(end-1),:);
            tau_uN_crv_1Dmod = nan(size(tau_uN_crv));
            for i = 1:ndep
                tau_uN_crv_1Dmod(i,:,:) = mods_1D{i}.tau_uN_crv(1:(end-1),:);
            end
            tau_uN_glb_N1mod = reshape(obj.jspc_N1mod.tau_uN_net(:,1,:),ndep,kor_obs,[]);

            size(dxuntns);
            size(tau_uqk_glb);
            size(tau_uk_glb);
            size(tau_uN_crv);
            size(tau_uN_crv_Yvthxglb);

            err_stats( 'tau_uqk_glb', comp_err_msr(tau_uqk_glb) );
            err_stats( 'tau_uk_glb', comp_err_msr(tau_uk_glb) );
            err_stats( 'tau_uN_crv', comp_err_msr(tau_uN_crv) );
            err_stats( 'tau_uN_crv_Yvthxglb', comp_err_msr(tau_uN_crv_Yvthxglb) );
            err_stats( 'tau_uN_crv_1D', comp_err_msr(tau_uN_crv_1Dmod) );
            err_stats( 'tau_uN_crv_N1', comp_err_msr(tau_uN_glb_N1mod) );

            Gsvd = obj.Gsvd;
            if ( Gsvd.r < Gsvd.dim )
                r_start = Gsvd.dim - min([5,Gsvd.dim-Gsvd.r]) + 1;
                for r = r_start:(Gsvd.dim)
                    fspc.print_vshort_polynomial_theta( Gsvd.V(:,r), obj.fspace.Pmat, ['\n K_G ' num2str(r) '/' num2str(Gsvd.dim) ' : '] )
                end
            end

            % err_dkxuq_tns_glb = nan(size(dxuntns));

            % tvf_glb = struct( ...
            %     'dkxuq_tns', dkxuq_tns_glb ...
            % );

            fprintf('\n\n (ending verification) \n\n');
            % debug_script
            obj_out = obj;

        end

        function fspace_o = Omap_solutions(fspace_i_,Sobs_,O_)
            nvar = length(fspace_i_.Omap_b);
            ndep = nvar-1;
            [Smat,nobs,nset,kor_obs,ndim_obs] = ldaux.unpack_Scell(Sobs_,ndep);
            xumat = Smat(1:nvar,:);

            if (O_(1) == 0)
                y_max = max(abs(xumat),[],2);
                y_min = zeros(size(y_max));
            else
                y_max = max(xumat,[],2);
                y_min = min(xumat,[],2);
            end

            a_out = reshape( (O_(2) - O_(1)) ./ (y_max-y_min), nvar,1);
            b_out = 0.5*((O_(2)-(a_out.*y_max)) + (O_(1)-(a_out.*y_min)));
            A_out = a_out.*eye(nvar);

            fspace_o = fspace_i_;
            fspace_o.Omap_b = b_out;
            fspace_o.Omap_A = A_out;
            if (isfield(fspace_o,'Omap_a'))
                fspace_o.Omap_a = a_out;
            end
        end

    end
end
