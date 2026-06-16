classdef mvp_N1_jspc_model
    properties (Constant)

    end
    properties

        Sobs;
        fspace;

        lambdas;
        LamN_tns;
        Renc_cell;

        Rmat_glb;
        Rsvd_glb;
        vth_glb;
        tau_uN_glb;
    end
    methods
        function obj = mvp_N1_jspc_model(Sobs_,fspace_)

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
            % normalize_Renc = @(r_) r_./sqrt(sum(r_.*r_,2));
            normalize_Renc = @(r_) r_;

            % get original jet space specs
            ncrv = length(Sobs_(:));
            ndep0 = size(fspace_.Pmat,1)-1;
            ndim0 = size(Sobs_{1},1);
            kor0 = ((ndim0-1)/ndep0) - 1;
            nvar0 = 1+ndep0;

            % convert to corresponding first order system
            kor = 1;
            ndep = (ndim0-1)-ndep0;
            ndim = 1 + ndep*(kor+1);
            nvar = 1+ndep;

            Sobs = cell([ncrv,1]);
            for icrv = 1:ncrv
                Sobs{icrv} = [ ...
                    Sobs_{icrv}(1:(end-ndep0),:) ; ...
                    Sobs_{icrv}((nvar0+1):end,:) ; ...
                ];
            end
            [Smat,nobs,~,kor_obs,ndim_obs,npts_per_crv,ipts_crv] = ldaux.unpack_Scell(Sobs,ndep);
            xumat = Smat(1:nvar,:);

            bor = fspace_.bor;
            [~, Pmat, ~] = ldaux.count_set_P_len(bor,ndep+1);

            %% restrict new function space to first order in all derivative coordinates
            Pmat = Pmat(:, sum(Pmat((nvar0+1):end,:),1)<=1 );
            Plen = size(Pmat,2);
            ntheta = nvar*Plen;
            %% initialize injection into Lambda column space
            i_imm = zeros(ndep*Plen,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            function l_imm = immerse_lambda_N1(l_)
                l_imm = zeros((ndep*Plen)*ndep,1);
                l_imm(logical( i_imm(:) )) = reshape( l_(:) * ones(1,ndep), [], 1);
                l_imm = (reshape(l_imm,ndep*Plen,ndep))';
                % l_ -> [ l_ 0 ... 0 ; 0 l_ ... 0 ; ... ; 0 0 ... l_ ]
            end

            fspace = fspace_;
            fspace.Plen = Plen;
            fspace.Pmat = Pmat;
            fspace.ntheta = ntheta;
            fspace.Omap_b = [ fspace_.Omap_b ; zeros(ndep-ndep0,1) ];
            fspace.Omap_A = [ ...
                fspace_.Omap_A , zeros(nvar0, ndep-ndep0) ;  ...
                zeros(ndep-ndep0, nvar0) , eye(ndep-ndep0) ...
            ];
            if (isfield(fspace_,'Omap_a'))
                fspace.Omap_a = [ fspace_.Omap_a(:) ; ones(ndep-ndep0,1) ];
            end
            fspace.imm_l = @(lambda_) immerse_lambda_N1(lambda_);
            fspace = adlam.init_fspace_family(fspace);

            lambdas = adlam.empty(nobs,0);
            lvs = nan(Plen,nobs);
            LamN_T_tns = nan(ntheta,1 + ndep*(kor_obs+1), nobs);
            Renc_cell = cell(1,nobs);
            Renc_tns = nan(ntheta,ndep,nobs);
            tic_prN = tic;
            for iobs = 1:nobs

                lambdas(iobs) = adlam(fspace, xumat(:,iobs), Smat((nvar+1):end,iobs) );

                l_i = lambdas(iobs);
                lv_i = l_i.lrow_vals;
                lvs(:,iobs) = lv_i;
                l_imm_i = fspace.imm_l(lv_i);

                LamN_T_tns(:,1:nvar,iobs) = [ [lv_i(:) ; zeros(ntheta-Plen,1) ] , [ zeros(Plen,ndep) ; l_imm_i' ] ];
                % Renc_cell{1,iobs} = [ -l_i.dxu(:,1)*lv_i , l_imm_i ];
                Renc_cell{1,iobs} = normalize_Renc([ -l_i.dxu(:,1)*lv_i , l_imm_i ]);
                Renc_tns(:,:,iobs) = Renc_cell{1,iobs}';
                LamN_T_tns(:,(nvar+1):(nvar+ndep),iobs) = [ -l_i.lkx(:,:,1) , fspace.imm_l(l_i.dkxl(1,:)) ]';

            end
            toc_prN = toc(tic_prN);
            fprintf('(mvp_N1_jspc_model) Prolonged %d observations over Q=%d, N=%d (B=%d) jet space with order O=%d mvpolynomials (C=%d) in %.3f seconds \n', ...
            nobs, ndep, kor_obs, ndim_obs, fspace.bor, ntheta, ...
            toc_prN);

            LamN_tns = permute(LamN_T_tns,[2 1 3]);

            function vth_out = comp_vartheta(W_,lvs_)
                vx_W = (W_(1:(size(W_,1)/nvar),:))' * lvs_;
                Deltas_W = sum(vx_W.*vx_W,1);
                vth_out = W_ * (vx_W./Deltas_W);
            end
            function tuN_out = comp_tau_uN_crv(th_,inds_,iP_)
                Plen_b = size(th_,1)/(ndep+1);
                nobs_crv = size(th_,2);
                tuN_out = nan(ndep,2,nobs_crv);
                for i = 1:nobs_crv
                    ii = inds_(i);
                    th_u_i_mat = reshape(th_((Plen_b+1):end,i),Plen_b,ndep);
                    tuN_out(:,1,i) = th_u_i_mat' * lvs(iP_,ii);
                    tuN_out(:,2,i) = (lambdas(ii).dkxl(1,iP_) * th_u_i_mat )' ...
                                        - (lambdas(ii).lkx(:,iP_,1) * th_(1:Plen_b,i));
                end
            end

            Rmat_glb = (reshape(Renc_tns,ntheta,ndep*nobs))';
            Rsvd_glb = Asvd_package(Rmat_glb);
            vth_glb = comp_vartheta( Rsvd_glb.W, lvs );
            tau_uN_glb = comp_tau_uN_crv( vth_glb, 1:nobs, 1:Plen );

            %% core initializations

            obj.Sobs = Sobs;
            obj.fspace = fspace;

            obj.lambdas = lambdas;
            obj.LamN_tns = LamN_tns;
            obj.Renc_cell = Renc_cell;

            obj.Rmat_glb = Rmat_glb;
            obj.Rsvd_glb = Rsvd_glb;
            obj.vth_glb = vth_glb;
            obj.tau_uN_glb = tau_uN_glb;
        end
    end
end
