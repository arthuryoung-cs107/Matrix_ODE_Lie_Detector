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

        Smat;

        lambdas;
        lvs;

        Renc_cell;

        Rqksvd_glb_cell;
        Rksvds_glb;
        Rk_glb_sub_svds;
        vth_k_glb;
        vth_k_glb_sub;
        tau_uk_glb;
        tau_uk_glb_sub;

        RNsvds_crv;
        RN_crv_sub_svds;
        vth_N_crv;
        vth_N_crv_sub;
        tau_uN_crv;
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

            [Smat,nobs,nset,kor_obs,ndim_obs,npts_per_crv,ipts_crv] = ldaux.unpack_Scell(Sobs_,ndep);
            ncrv = length(npts_per_crv);

            xumat = Smat(1:nvar,:);
            untns = reshape( Smat(2:end,:), ndep,kor_obs+1,nobs );
            dxuntns = reshape(untns(:,2:end,:),ndep,kor_obs,nobs);

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
                % l_ -> [l_ 0 ... 0 ; 0 l_ ... 0 ; ...]
            end

            lambdas = adlam.empty(nobs,0);
            lvs = nan(Plen,nobs);
            Renc_cell = cell(kor_obs,nobs);
            for iobs = 1:nobs

                lambdas(iobs) = adlam(fdat_, xumat(:,iobs), Smat((nvar+1):end,iobs) );

                l_i = lambdas(iobs);
                lv_i = adlam.lrowP(l_i);

                lvs(:,iobs) = lv_i;

                Renc_cell{1,iobs} = [ -l_i.dxu(:,1)*lv_i , immerse_lambda(lv_i) ];
                for k = 2:kor_obs
                    Renc_cell{k,iobs} = [ ...
                    -l_i.lkx(:,:,k-1) - l_i.dxu(:,k)*lv_i , immerse_lambda(l_i.dkxl(k-1,:)) ...
                    ];
                end
                % Renc_cell{1,iobs}
                % Renc_cell{2,iobs}
                % pause
            end

            shape_Rkt_glb = @(RkcT_) permute(reshape( horzcat(RkcT_{:}),ntheta,ndep,nobs),[3 1 2]);
            Rkc_2_Rkt_glb = @(Rkc_) shape_Rkt_glb( ...
            cellfun( @(R_) R_',Rkc_,'UniformOutput', false) ...
            );
            function Asvd_out = Asvd_package(Ri_);
                [r_,s_,V_,U_] = fspc.rsV_unpack(Ri_);
                Asvd_out = struct( ...
                    'A', @(o_) o_.U * ( o_.s(:) .* (o_.V') ), ...
                    'r', r_, ...
                    's', s_, ...
                    'V', V_, ...
                    'U', U_, ...
                    'D', ( s_(:)  )' .* V_, ... % rowspace (domain) Y approx. fiber
                    'W', ( s_(end) ./ s_(:)  )' .* V_ ... % nullspace (kernal) approx. fiber
                );
            end
            function vth_out = comp_vartheta(W_,lvs_)
                vx_W = (W_(1:( length(W_)/(ndep+1) ),:))' * lvs_;
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

            is_pPbs_mat = zeros(bor,Plen);
            is_pTbs_mat = zeros(bor,ntheta);
            for b = 1:bor
                is_pPbs_mat(b,:) = mvp_jspc_model.is_prjP_b_sub(Pmat,b);
                is_pTbs_mat(b,:) = is_prjtheta_b_sub(b);
            end
            is_pPbs_mat = logical(is_pPbs_mat);
            is_pTbs_mat = logical(is_pTbs_mat);
            bor_to_one = flip(1:bor);

            Rqksvd_glb_cell = cell([ndep,kor_obs]);
            Rksvds_glb = cell([1,kor_obs]);
            Rk_glb_sub_svds = cell([bor,kor_obs]);
            vth_k_glb = nan(ntheta,nobs,kor_obs);
            vth_k_glb_sub = cell([bor,kor_obs]);
            tau_uk_glb = nan(ndep,nobs,kor_obs+1);
            tau_uk_glb_sub = nan(ndep,nobs,bor,kor_obs+1);
            % compute global models
            for k = 1:kor_obs
                Rkt_glb = Rkc_2_Rkt_glb(Renc_cell(k,:));
                [Rk_glb_net,Rk_crv_net] = deal(zeros(ntheta,ntheta,ndep));
                for idep = 1:ndep
                    Rqksvd_glb_cell{idep,k} = Asvd_package( Rkt_glb(:,:,idep) );
                    Rk_glb_net(:,:,idep) = Rqksvd_glb_cell{idep,k}.D;
                end
                Rk_glb_net = reshape(Rk_glb_net,ntheta,[])';
                Rksvds_glb{k} = Asvd_package(Rk_glb_net);
                vth_k_glb(:,:,k) = comp_vartheta( ...
                    Rksvds_glb{k}.W, lvs ...
                );
                % estimate of k'th derivative, global model
                tau_uk_glb(:,:,k) = comp_tau_uk_glb(vth_k_glb(:,:,k),k,1:Plen);
                if (k == kor_obs)
                    tau_uk_glb(:,:,k+1) = comp_tau_uk_glb(vth_k_glb(:,:,k),k+1,1:Plen);
                end
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
            tau_uk_glb = permute(tau_uk_glb, [1 3 2] ); % ndep,kor_obs+1,nobs
            tau_uk_glb_sub = permute(tau_uk_glb_sub, [1 4 2 3] ); % ndep,kor_obs+1,nobs,bor

            RNsvds_crv = cell([1,ncrv]);
            RN_crv_sub_svds = cell([bor,ncrv]);
            vth_N_crv = nan(ntheta,nobs);
            vth_N_crv_sub = cell([bor,ncrv]);
            tau_uN_crv = nan(ndep,kor_obs+1,nobs);
            tau_uN_crv_sub = nan(ndep,kor_obs+1,nobs,bor);
            % compute local curve models
            for icrv = 1:ncrv
                inds_i = ipts_crv(1,icrv):ipts_crv(2,icrv);
                RN_enc_i = reshape(Renc_cell(:, inds_i), [], 1);
                RN_mat_i = vertcat( RN_enc_i{:} );

                % compute P permutation (full) local curve model
                RNsvds_crv{icrv} = Asvd_package(RN_mat_i);
                % evaluate local parameters predicted by local curve model
                vth_N_crv(:,inds_i) = comp_vartheta( ...
                    RNsvds_crv{icrv}.W,lvs(:,inds_i) ...
                );
                % estimate all curve derivatives
                tau_uN_crv(:,:,inds_i) = comp_tau_uN_crv( ...
                    vth_N_crv(:,inds_i), inds_i, 1:Plen ...
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

            % fspc.print_vshort_polynomial_theta(Rksvds_glb{1}.W(:,end),Pmat)
            % fspc.print_vshort_polynomial_theta(Rksvds_glb{1}.W(:,end-1),Pmat)

            %% assignments

            obj.Sobs = Sobs_;
            obj.Sdat = Sdat_;
            obj.fdat =  fdat_;

            obj.Smat = Smat;

            obj.lambdas = lambdas;
            obj.lvs = lvs;
            obj.Renc_cell = Renc_cell;

            obj.Rqksvd_glb_cell = Rqksvd_glb_cell;
            obj.Rksvds_glb = Rksvds_glb;
            obj.Rk_glb_sub_svds = Rk_glb_sub_svds;
            obj.vth_k_glb = vth_k_glb;
            obj.vth_k_glb_sub = vth_k_glb_sub;
            obj.tau_uk_glb = tau_uk_glb;
            obj.tau_uk_glb_sub = tau_uk_glb_sub;

            obj.RNsvds_crv = RNsvds_crv;
            obj.RN_crv_sub_svds = RN_crv_sub_svds;
            obj.vth_N_crv = vth_N_crv;
            obj.vth_N_crv_sub = vth_N_crv_sub;
            obj.tau_uN_crv = tau_uN_crv;
            obj.tau_uN_crv_sub = tau_uN_crv_sub;
        end
    end

    methods (Static)

        function obj_out = verify(obj)

            ndep = obj.Sdat.ndep;
            [ndim_obs,nobs] = size(obj.Smat);
            kor_obs = (ndim_obs-1)/ndep - 1;

            Smat = obj.Smat;
            tau_uk_glb = obj.tau_uk_glb;
            tau_uk_glb_sub = obj.tau_uk_glb_sub;

            xvec = Smat(1,:);
            unmat = Smat(2:end,:);
            untns = reshape( unmat, ndep, kor_obs+1, nobs );
            umat = reshape( untns(:,1,:), ndep, nobs );
            dxutns = reshape( untns(:,2:end,:), ndep, kor_obs, nobs );



            obj_out = obj;

        end
    end
end
