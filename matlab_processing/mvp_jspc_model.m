classdef mvp_jspc_model
    properties (Constant)

        is_prjP_b_sub = @(P_,b_) sum(P_,1) <= b_;
        pick_prjP_b_sub = @(is_,ifull_) ifull_(is_);
        inds_prjP_b_sub = @(P_,b_) mvp_jspc_model.pick_prjP_b_sub(sum(P_,1) <= b_, 1:size(P_,2));
        prjP_b_sub = @(P_,b_) P_( : , sum(P_,1) <= b_ );

    end

    properties

        Sobs;
        Sdat;
        fdat;

        lambdas;
        lvs;

        Renc_cell;
        Rsvd_cell;
        Rksvds;
        Rk_sub_svds;

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

            [Smat,nobs,nset,kor_obs,ndim_obs] = ldaux.unpack_Scell(Sobs_,ndep);
            xumat = Smat(1:nvar,:);
            untns = reshape( Smat(2:end,:), ndep,kor_obs+1,nobs );
            dxuntns = reshape(untns(:,2:end,:),ndep,kor_obs,nobs);

            ntheta = nvar*Plen;
            %% initialize injection into Lambda column space
            ncol_Lambda_u = ntheta-Plen;
            i_imm = zeros(ncol_Lambda_u,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            i_imm_vec = logical( i_imm(:) );
            len_Lambda_u = prod(size(i_imm));
            l_imm_init = zeros(len_Lambda_u,1);
            ones_imm = ones(ndep,1);
            function l_imm = immerse_lambda(l_)
                l_imm = l_imm_init;
                l_imm(i_imm_vec) = reshape( (ones_imm*l_)', [], 1);
                l_imm = reshape(l_imm,ndep,ncol_Lambda_u);
                % l_ -> [l_ 0 ... 0 ; 0 l_ ... 0 ; ...]
            end

            lambdas = adlam.empty(nobs,0);
            lvs = nan(Plen,nobs);
            Renc_cell = cell(kor_obs,nobs);
            for i = 1:nobs

                lambdas(i) = adlam(fdat_, xumat(:,i), Smat((nvar+1):end,i) );

                l_i = lambdas(i);
                lv_i = adlam.lrowP(l_i);

                lvs(:,i) = lv_i;

                Renc_cell{1,i} = [(-l_i.dxu(:,1).*lv_i),immerse_lambda(lv_i)];
                for k = 2:kor_obs
                    Renc_cell{k,i} = [ ...
                    -(l_i.lkx(:,:,k-1)+(l_i.dxu(:,k).*lv_i)),immerse_lambda(l_i.dkxl(k-1,:)) ...
                    ];
                end
            end

            shape_Rkt = @(RkcT_) permute(reshape( horzcat(RkcT_{:}),ntheta,ndep,nobs),[3 1 2]);
            Rkc_2_Rkt = @(Rkc_) shape_Rkt( ...
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
            function tuk_out = comp_tau_uk(th_,k_,iP_)
                Plen_b = size(th_,1)/(ndep+1);
                tuk_out = nan(ndep,nobs);
                for i = 1:nobs
                    tuk_out(:,i) =  ...
                    ( lambdas(i).dkxl(k_,iP_) * reshape(th_((Plen_b+1):end,i),Plen_b,ndep) )' ...
                     - ( lambdas(i).lkx(:,iP_,k_)*th_(1:Plen_b,i) );
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

            Rsvd_cell = cell([ndep,kor_obs]);
            Rksvds = cell([1,kor_obs]);
            Rk_sub_svds = cell([bor,kor_obs]);

            varthetas_k = nan(ntheta,nobs,kor_obs);
            varthetas_k_sub = cell([bor,kor_obs]);
            tau_uk = nan(ndep,nobs,kor_obs);
            tau_uk_sub = nan(ndep,nobs,bor,kor_obs);
            for k = 1:kor_obs
                Rkt = Rkc_2_Rkt(Renc_cell(k,:));

                Rk_net = zeros(ntheta,ntheta,ndep);
                for i = 1:ndep
                    Rsvd_cell{i,k} = Asvd_package( Rkt(:,:,i) );
                    Rk_net(:,:,i) = Rsvd_cell{i,k}.D;
                end
                Rk_net = reshape(Rk_net,ntheta,[])';
                Rksvds{k} = Asvd_package(Rk_net);
                varthetas_k(:,:,k) = comp_vartheta(Rksvds{k}.W, ...
                    lvs ...
                );
                tau_uk(:,:,k) = comp_tau_uk(varthetas_k(:,:,k),k,1:Plen);

                for b = bor_to_one
                    Rk_sub_svds{b,k} = Asvd_package( ...
                        Rk_net( :, is_pTbs_mat(b,:) ) ...
                    );
                    varthetas_k_sub{b,k} = comp_vartheta(Rk_sub_svds{b,k}.W, ...
                        lvs(is_pPbs_mat(b,:),:) ...
                    );
                    tau_uk_sub(:,:,b,k) = comp_tau_uk(varthetas_k_sub{b,k},k,is_pPbs_mat(b,:));
                end
            end
            R1svd = Rksvds{end,1} % computed by intersecting svds
            fspc.print_vshort_polynomial_theta(R1svd.W(:,end),Pmat)
            fspc.print_vshort_polynomial_theta(R1svd.W(:,end-1),Pmat)

            %% assignments

            obj.Sobs = Sobs_;
            obj.Sdat = Sdat_;
            obj.fdat =  fdat_;

            obj.lambdas = lambdas;
            obj.lvs = lvs;
            obj.Renc_cell = Renc_cell;
            obj.Rsvd_cell = Rsvd_cell;
            obj.Rksvds = Rksvds;
            obj.Rk_sub_svds = Rk_sub_svds;

        end
    end

    methods (Static)

    end

end
