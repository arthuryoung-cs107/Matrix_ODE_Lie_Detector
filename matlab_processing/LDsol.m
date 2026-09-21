classdef LDsol
    properties (Constant)

    end
    properties
        xu;

        lamN1;
        lam0;
        lamRN1;

    end
    methods (Static)
        function mod_out = model_solspace(Sobs_,dat_,fmap_)
            if (isfield(dat_,'bor_max'))
                bor_max = dat_.bor_max;
            else
                % bor_max = 8;
                bor_max = 6;
            end

            %% overarching package for all SVD usage
            function [Asvd_out,U_] = Asvd_package(A_);
                dim_scl = max(size(A_));
                if (nargout==1)
                    [~,S_,V_] = svd(A_,'econ');
                else
                    [U_,S_,V_] = svd(A_,'econ');
                end
                s_ = reshape(diag(S_),[],1); % output s_ is col vector
                r_ = sum(double(s_ > dim_scl*eps(S_(1)))); % default matlab tol

                Asvd_out = struct( ...
                    'mrow', size(A_,1), ...
                    'ncol', size(A_,2), ...
                    'dim', length(s_), ...
                    'r', r_, ...
                    's', s_, ...
                    'V', V_, ...
                    'D', ( s_(:)  )' .* V_, ... % rowspace (domain) Y approx. basis
                    'W', ( s_(end) ./ s_(:)  )' .* V_ ... % nullspace (kernal) K approx. basis
                );
            end
            normalize_Renc = @(r_) r_./sqrt(sum(r_.*r_,2)); % unit length Rmat rows
            normalize_DprN = @(d_) d_./sqrt( sum(d_.*d_,2) ); % unit length DprN cols (rows after transposition)
            normalize_Genc = @(g_) g_./sqrt(sum(g_.*g_,2)); % unit length Gmat rows
            normalize_Benc = @(b_) b_./sqrt(sum(b_.*b_,2)); % unit length Bmat rows
            normalize_Henc = @(h_) h_./sqrt(sum(h_.*h_,2)); % unit length Hmat rows
            normalize_Tenc = @(t_) t_./sqrt(sum(t_.*t_,2)); % unit length Tmat rows
            % normalize_LamWtns = @(LW_) reshape();

            [ndep,nvar] = deal(dat_.ndep, dat_.ndep+1);
            [Smat,nobs,ncrv,kor,ndim,npts_per_crv,ipts_crv] = ldaux.unpack_Scell(Sobs_,ndep);
            xumat = Smat(1:nvar,:);
            xvec = xumat(1,:);
            umat = xumat(2:nvar,:);
            uNtns = reshape( Smat(2:end,:) , ndep , kor+1 , nobs );
            uNm1mat = reshape( uNtns(:,1:(end-1),:), ndep*kor, nobs );
            dxuNtns = uNtns(:,2:end,:);
            dxumat = reshape( dxuNtns, ndep*kor, nobs );

            %% X | U0 | ... | UN --> X | U0 ... UNm1 | U1 ... UN
            xumat_N1 = [ xvec ; uNm1mat ];
            Smat_N1 = [ xumat_N1 ;  dxumat ];

            fmap_.Omap_a = fmap_.Omap_a(:);
            fmap_.Omap_A = fmap_.Omap_a .* eye(length(fmap_.Omap_a));

            ndep_N1 = ndep*kor;
            nvar_N1 = 1+ndep_N1;
            ndim_N1 = 1 + 2*ndep_N1;
            iN1_2_N = [1:nvar_N1 (ndim_N1-ndep+1):ndim_N1];
            bor_N1 = 1;
            [Plen_N1, Pmat_N1, ~] = ldaux.count_set_P_len(bor_N1,ndep_N1+1);
            while ((Plen_N1 < nobs) && (bor_N1 < bor_max))
                bor_N1 = bor_N1 + 1;
                [Plen_N1, Pmat_N1, ~] = ldaux.count_set_P_len(bor_N1,ndep_N1+1);
            end
            if (Plen_N1>=nobs)
                %% pare down jet space model if fspace is too big for overdetermined H_N1 matrix
                if (kor==1)
                    %% pare down trivially
                    bor_N1 = bor_N1 - 1;
                    [Plen_N1, Pmat_N1, ~] = ldaux.count_set_P_len(bor_N1,ndep_N1+1);
                else
                    %% pare down by restricting polynomial order of derivative terms
                    ord_i = bor_N1;
                    Pmat_N1 = Pmat_N1(:, sum(Pmat_N1((nvar+1):end,:),1) <= ord_i ); % restrict to b'th order polys in dxu
                    while ((size(Pmat_N1,2) >= nobs)&&(ord_i > 1))
                        ord_i = ord_i-1;
                        Pmat_N1 = Pmat_N1(:, sum(Pmat_N1((nvar+1):end,:),1) <= ord_i );
                    end
                    Plen_N1 = size(Pmat_N1,2);
                end
            end

            fspace_N1 = fmap_;
            fspace_N1.Omap_a = [
                fmap_.Omap_a(:) ;  ...
                reshape( (fmap_.Omap_a(1)).^(-(1:(kor-1))) .* ( fmap_.Omap_a(2:end) .* ones(1,kor-1) ) ,[],1 ) ...
            ];
            fspace_N1.Omap_b = [ fspace_N1.Omap_b(:) ; zeros((kor-1).*ndep,1) ];
            fspace_N1.Omap_A = fspace_N1.Omap_a .* eye(length(fspace_N1.Omap_a));
            fspace_N1.bor = bor_N1;
            fspace_N1.Plen = Plen_N1;
            fspace_N1.Pmat = Pmat_N1;
            fspace_N1.ntheta = (1+ndep_N1)*Plen_N1;
            fspace_N1 = adlam.init_fspace_family(fspace_N1);

            if (kor>1)
                bor_0 = min([bor_max , max([ 1 , floor( ((nobs.*ndep)./nvar).^(1./nvar) - 1  ) ])]);
                [Plen_0, Pmat_0, ~] = ldaux.count_set_P_len(bor_0,ndep+1);

                fspace_0 = fmap_;
                fspace_0.Omap_A = fspace_0.Omap_a .* eye(length(fspace_0.Omap_a));
                fspace_0.bor = bor_0;
                fspace_0.Plen = Plen_0;
                fspace_0.Pmat = Pmat_0;
                fspace_0.ntheta = (1+ndep)*Plen_0;
                fspace_0 = adlam.init_fspace_family(fspace_0);

                fspace_RN1 = fspace_N1;
                % ord_i = max( reshape(Pmat_N1((nvar+1):end,:),[],1) );
                ord_i = max( sum(Pmat_N1((nvar+1):end,:),1) );
                mrow_RN1_net = nobs*(ndep*( 2*kor - 1 )); % = nobs*( ndep*kor + ndep*(kor-1) )
                mrow_GN1_net_full = nobs*ndep*kor;
                % while ( (mrow_RN1_net <= fspace_RN1.ntheta)&&(ord_i>1) )
                while ( (mrow_GN1_net_full <= fspace_RN1.ntheta)&&(ord_i>1) ) % R and G same Lambda space
                    Pmat_i = fspace_RN1.Pmat;
                    ord_i = ord_i - 1;
                    fspace_RN1.Pmat = Pmat_i( :, sum(Pmat_i((nvar+1):end,:),1) <= ord_i );
                    fspace_RN1.Plen = size(fspace_RN1.Pmat,2);
                    fspace_RN1.ntheta = nvar_N1*(fspace_RN1.Plen);
                end
                if (mrow_RN1_net <= fspace_RN1.ntheta)
                    fprintf('(LDsol::model_solspace) honestly you should get some more points (mrow_RN1_net = %d, ntheta=%d) \n', ...
                    mrow_RN1_net, fspace_RN1.ntheta );
                end
            else
                [bor_0,Plen_0,Pmat_0] = deal(bor_N1,Plen_N1,Pmat_N1);
                [fspace_0,fspace_RN1] = deal(fspace_N1); % all identical in the case of N = 1

                ord_i = max(sum(Pmat_N1(2:end,:),1));
                mrow_RN1_net = nobs*ndep; % = nobs*( ndep*kor + ndep*(kor-1) )
                while ( (mrow_RN1_net <= fspace_RN1.ntheta)&&(ord_i>1) )
                    Pmat_i = fspace_RN1.Pmat;
                    ord_i = ord_i - 1;
                    fspace_RN1.Pmat = Pmat_i( :, sum(Pmat_i(2:end,:),1) <= ord_i );
                    fspace_RN1.Plen = size(fspace_RN1.Pmat,2);
                    fspace_RN1.ntheta = nvar_N1*(fspace_RN1.Plen);
                end
            end
            Pmat_RN1_full = fspace_RN1.Pmat;
            Plen_RN1 = fspace_RN1.Plen;
            ntheta_RN1 = fspace_RN1.ntheta;
            %% initialize injection into RN1 Lambda column space
            i_imm_RN1 = zeros(ndep_N1*Plen_RN1,ndep_N1);
            for i = 1:ndep_N1
                idel = (i-1)*Plen_RN1;
                i_imm_RN1( (1+idel):(Plen_RN1+idel), i ) = 1;
            end
            function l_imm = immerse_lambda_RN1(l_)
                l_imm = zeros((ndep_N1*Plen_RN1)*ndep_N1,1);
                l_imm(logical( i_imm_RN1(:) )) = reshape( l_(:) * ones(1,ndep_N1), [], 1);
                l_imm = (reshape(l_imm,ndep_N1*Plen_RN1,ndep_N1))';
                % l_ -> [ l_ 0 ... 0 ; 0 l_ ... 0 ; ... ; 0 0 ... l_ ]
            end
            function Lam_out = immerse_l_Lamspc(l_,nvar_)
                Plen_ = length(l_(:));
                Lam_out = zeros(Plen_,nvar_,nvar_);
                for ivar = 1:nvar_
                    Lam_out(:,ivar,ivar) = l_(:);
                end
                Lam_out = reshape(Lam_out,Plen_*nvar_,nvar_)'; % result is nvar x ntheta Lambda matrix
            end
            fspace_RN1.imm_l = @(lambda_) immerse_lambda_RN1(lambda_);

            Hmat_N1 = nan(nobs,Plen_N1);
            Hmat_0 = nan(nobs,Plen_0);
            Hmat_RN1 = nan(nobs,Plen_RN1);
            lvs_N1 = nan(Plen_N1,nobs);
            Jltns_N1 = nan(nvar_N1,Plen_N1,nobs);
            lvs_RN1 = nan(Plen_RN1,nobs);
            Jltns_RN1 = nan(nvar_N1,Plen_RN1,nobs);
            dxl_RN1 = nan(Plen_RN1,nobs);
            lNx_RN1 = nan(ndep_N1,Plen_RN1,nobs);
            LamN_T_tns_RN1 = nan(ntheta_RN1,ndim_N1,nobs);
            Renc_tns_RN1 = nan(ntheta_RN1,ndep_N1,nobs);
            DprN_T_ttns = zeros(ntheta_RN1,ndep,kor-1,nobs);
            tic0 = tic;
            for iobs = 1:nobs
                sols(iobs) = LDsol(Smat(:,iobs));

                sols(iobs).lamN1 = adlam( fspace_N1, xumat_N1(:,iobs) );
                lvs_N1(:,iobs) = (sols(iobs).lamN1.lrow_vals)';
                Jltns_N1(:,:,iobs) = sols(iobs).lamN1.Jl;
                Hmat_N1(iobs,:) = [ 1.0 , dxumat(:,iobs)' ]*Jltns_N1(:,:,iobs);

                sols(iobs).lam0 = adlam( fspace_0, xumat(:,iobs), Smat((nvar+1):end,iobs) );
                Hmat_0(iobs,:) = [ 1.0 , dxumat(1:ndep,iobs)' ]*(sols(iobs).lam0.Jl);

                sols(iobs).lamRN1 = adlam( fspace_RN1, xumat_N1(:,iobs), Smat((nvar+1):end,iobs) );
                Hmat_RN1(iobs,:) = [ 1.0 , dxumat(:,iobs)' ]*(sols(iobs).lamRN1.Jl);

                lRN1_i = sols(iobs).lamRN1;
                lv_RN1_i = lRN1_i.lrow_vals;
                dxl_RN1_i = lRN1_i.dkxl(1,:);
                l_imm_i = fspace_RN1.imm_l(lv_RN1_i);

                lvs_RN1(:,iobs) = lv_RN1_i;
                Jltns_RN1(:,:,iobs) = lRN1_i.Jl;
                dxl_RN1(:,iobs) = dxl_RN1_i;
                lNx_RN1(:,:,iobs) = lRN1_i.lkx(:,:,1);
                LamN_T_tns_RN1(:,1:nvar_N1,iobs) = ...
                    [ [ lv_RN1_i(:) ; zeros(ntheta_RN1-Plen_RN1,1) ] , [ zeros(Plen_RN1,ndep_N1) ; l_imm_i' ] ];
                Renc_tns_RN1(:,:,iobs) = ([ -lRN1_i.dxu(:,1)*lv_RN1_i , l_imm_i ])';

                LamN_T_tns_RN1(:,(nvar_N1+1):(nvar_N1+ndep_N1),iobs) = ...
                    [ -lRN1_i.lkx(:,:,1) , fspace_RN1.imm_l(dxl_RN1_i) ]';
                iiLam0 = nvar + (1:ndep);
                iiLam1 = iiLam0-nvar+nvar_N1;
                for k = 2:kor
                    DprN_T_ttns(:,:,k-1,iobs) = LamN_T_tns_RN1(:,iiLam1,iobs) - LamN_T_tns_RN1(:,iiLam0,iobs);
                    iiLam0 = iiLam0 + ndep;
                    iiLam1 = iiLam1 + ndep;
                end
            end
            toc1 = toc(tic0);
            fprintf('(LDsol::model_solspace) Prolonged %d observations over Q=%d, N=%d (B=%d) jet space with order O=%d mvpolynomials (C_N1=%d, C_0=%d, C_RN1=%d), encoded %dx%d R matrix in %.3f seconds \n', ...
            nobs, ndep, kor, ndim, fspace_N1.bor, ...
            fspace_N1.ntheta, fspace_0.ntheta, fspace_RN1.ntheta, ...
            ndep_N1*nobs, ntheta_RN1, ...
            toc1);

            LamN_T_ttns_RN1 = reshape(LamN_T_tns_RN1,[Plen_RN1 nvar_N1 ndim_N1 nobs]);
            % Plen x nvar_N1 x ndim x nobs, Lambda tensors over first jet space
            LamN1_T_ttns = LamN_T_ttns_RN1(:,:,[ 1:nvar_N1 (ndim_N1-ndep+1):ndim_N1 ],:);
            LamN_tns_RN1 = permute(LamN_T_tns_RN1,[2 1 3]); % --> ndim (N1) x ntheta_N1 x nobs
            Lam_dNxu_RN1_i = @(i_) LamN_tns_RN1((end-ndep+1):end,:,i_); % ndep x ntheta_N1, acts on theta vectors

            Jl_N1_svd = Asvd_package( reshape(permute(Jltns_N1,[2 1 3]),Plen_N1,nvar_N1*nobs)' );
            Jl_N1_svd.YJl = Jl_N1_svd.V(:,1:(Jl_N1_svd.r));
            Jl_N1_svd.H_YJl_svd = Asvd_package( normalize_Henc(Hmat_N1*Jl_N1_svd.YJl) );

            H_N1_svd = Asvd_package(normalize_Henc(Hmat_N1));
            H_0_svd = Asvd_package(normalize_Henc(Hmat_0));
            H_RN1_svd = Asvd_package(normalize_Henc(Hmat_RN1));

            function vth_out = comp_vartheta_RN1(W_,lvs_)
                vx_W = (W_(1:(size(W_,1)/nvar_N1),:))' * lvs_;
                vth_out = W_ * (vx_W./sum(vx_W.*vx_W,1));
            end
            function tuN_out = comp_tau_uN_RN1(th_,iP_) % global N1 same as crv Nfull
                Plen_b = size(th_,1)/nvar_N1;
                tuN_out = nan(ndep_N1,2,nobs);
                for i = 1:nobs
                    th_u_i_mat = reshape(th_((Plen_b+1):end,i),Plen_b,ndep_N1);
                    tuN_out(:,1,i) = th_u_i_mat' * lvs_RN1(iP_,i);
                    tuN_out(:,2,i) = ( sols(i).lamRN1.dkxl(1,iP_) * th_u_i_mat )' ...
                                        - (sols(i).lamRN1.lkx(:,iP_,1) * th_(1:Plen_b,i) );
                end
            end
            function [f_s0_out,dxf_s0_out,vth_out,lam_out,Jtu_s0_out,JtuN_s0_out] = comp_f_s0(s0_,fspc_,W_,iP_)
                nvar_ = length(s0_(:));
                Plen_b = size(fspc_.Pmat(:,iP_),2);
                lam_out = adlam( fspc_, s0_ );
                vx_W = ( W_(1:Plen_b,:) )' * reshape( lam_out.lrow_vals(iP_) , [], 1 );
                vth_out = W_ * (vx_W./sum(vx_W.*vx_W,1));
                [vTh_x,vTh_u] = deal( vth_out(1:Plen_b) , reshape(vth_out((Plen_b+1):end),Plen_b,[]) );
                f_s0_out = vTh_u' * (lam_out.lrow_vals(iP_))';
                lam_out = lam_out.prolong_jet_space( reshape(f_s0_out,[],1) );
                dxf_s0_out = vTh_u'*lam_out.dkxl(1,iP_)' - lam_out.lkx(:,iP_,1)*vTh_x;
                if (nargout>=5)
                    JtuN_s0_out = lam_out.J_tau_uN(vth_out);
                    Jtu_s0_out = JtuN_s0_out(:,1:nvar_,1);
                end
            end
            function [inds_i,ncol_i,ord_i] = pare_colspc(Pmat_,mrow_)
                inds_ = 1:(size(Pmat_,2));
                inds_i = inds_; ncol_i = nvar_N1*length(inds_); ord_i = max(sum(Pmat_((nvar+1):end,:),1))-1;
                while ( ( mrow_<=ncol_i )&&(ord_i>1) )
                    inds_i = inds_(sum(Pmat_((nvar+1):end,:),1) <= ord_i);
                    ncol_i = nvar_N1*length(inds_i);
                    if ( ncol_i < mrow_ )
                        break;
                    else
                        ord_i = ord_i-1;
                    end
                end
            end
            function sO_basis_out = Tspc_package(VW_sO_,tO_)
            % function [sO_basis_out,sO_PCA_out] = Tspc_package(VW_sO_,tO_)
                [sO_basis_out,U_V_sO] = Asvd_package(VW_sO_); % V is ntheta x nvar, linearly combine W columns, yield phi
                sO_basis_out.U = U_V_sO; % nvar x nvar, orthogonal basis for tangent space of S at sO
                sO_basis_out.Tspc_basis = sO_basis_out.U; % U is a column orthonormal basis for T_sO S

                B_in = size(U_V_sO,1);
                tO_B = reshape(tO_(1:B_in),B_in,1);
                mag_tO_B = norm(tO_B);

                sO_basis_out.Uscl = mag_tO_B; % choose to scale w.r.t. norm of trivial tvector at sO
                % image of the trivial tvector over chosen basis for T_sO S
                sO_basis_out.x_Tspc_t_sO = (sO_basis_out.Tspc_basis)' * tO_B;

                % decompose section of tangent space orthogonal to tvf tangent space
                tO_B_unit = tO_B/mag_tO_B;
                VW_YVW_sO = VW_sO_*sO_basis_out.V(:,1:nvar_N1);
                [sO_nTVF,U_V_nTVF_sO] = Asvd_package((eye(B_in)-tO_B_unit*(tO_B_unit'))*VW_YVW_sO);
                sO_nTVF.U = U_V_nTVF_sO;
                sO_nTVF.Tspc_image = VW_YVW_sO * ( ( sO_nTVF.s(1:(nvar_N1-1))/sO_nTVF.s(1) )' .* sO_nTVF.V(:,1:(nvar_N1-1)) );
                sO_nTVF.nTVF_Tspc_image = (sO_nTVF.s(1:(nvar_N1-1)) / sO_nTVF.s(1))' .* U_V_nTVF_sO(:,1:(nvar_N1-1));

                sO_basis_out.sO_nTVF = sO_nTVF;
            end
            function [s0O_basis_out,sNO_basis_out] = compute_sO_W_tspc(W_,iP_,lam_sO_,tO_)
                lv_b_sO = lam_sO_.lrow_vals(iP_);
                dkxl_b_sO = lam_sO_.dkxl(:,iP_);
                lkx_b_sO = lam_sO_.lkx(:,iP_,:);
                nvar_ = 1+size(lkx_b_sO,1);
                Plen_b = length(lv_b_sO);
                ntheta_b = Plen_b*nvar_;
                % nvar x nvar matrix, columns span tangent space of S0 at s0O
                V0spc_image = @(thtns_) reshape( ...
                    pagemtimes( thtns_, lv_b_sO' ) , nvar_,nvar_ );
                % VdNxuspc_sO is ndep x nvar matrix, columns are vfield coeffs of dNxu in jet space
                VdNxuspc_image = @(thmat_,thtns_) reshape( ...
                    pagemtimes( thtns_((end-ndep+1):end,:,:) , dkxl_b_sO(1,:)' ) , [], nvar_ ...
                    ) - lkx_b_sO((end-ndep+1):end,:,1) * thmat_(1:Plen_b,:) ;
                % B x nvar matrix, columns are vfield coeffs, span tangent space of SN at sNO, lie algebra at origin
                VNspc_image = @(thmat_,thtns_) deal([V0spc_image(thtns_) ; VdNxuspc_image(thmat_,thtns_)],thtns_);
                % thtns_ is nvar x Plen x nvar, pages are coordinate vfield parameter matrices, act on lambda vectors, yield coeffs
                complete_sO_Tspc_image = @(thmat_) VNspc_image(thmat_,permute(reshape(thmat_,Plen_b,nvar_,nvar_),[2 1 3]));

                %% evaluate image of parameter space at the origin
                % nvar x Plen x ntheta, pages multiply column lambda vecs
                Wtns = permute(reshape(W_,Plen_b,nvar_,ntheta_b) , [2 1 3]);
                % V0W_sO is nvar x ntheta, cols are base space vector coefficients which span T_sO S0
                V0W_sO = reshape(pagemtimes( Wtns , lv_b_sO(:) ), nvar_, ntheta_b);

                % V of V0W_sO is an orthonormal basis for parameters of all vfields over base space, U an orthogonal basis for T_sO S0
                s0O_basis_out = Tspc_package( V0W_sO,tO_ );
                % [s0O_basis_out,s0O_PCA_out] = Tspc_package( V0W_sO,tO_ );
                s0O_basis_out.WYmu = W_ * s0O_basis_out.V(:,1:nvar_);
                % theta_WV_sO is ntheta x nvar parameter matrix, cols lincom lambda fcns, yield vfield coeffs
                s0O_basis_out.theta_WV_sO = ( s0O_basis_out.Uscl(:)' / s0O_basis_out.s(1)) .* s0O_basis_out.WYmu;
                % VNspc_sO is B x nvar matrix, columns are vfield coeffs, span tangent space of SN at sNO, lie algebra at origin
                [s0O_basis_out.Vspc_sO,s0O_basis_out.theta_tns_WV_sO] = complete_sO_Tspc_image( s0O_basis_out.theta_WV_sO );

                if (nargout == 2)
                    % VdNxuW_sNO is ndep x ntheta, cols are N'th derivative vector field coefficients in N'th jet space
                    VdNxuW_sNO = reshape( pagemtimes( Wtns((end-ndep+1):end,:,:), dkxl_b_sO(1,:)' ), ndep, ntheta_b )  ...
                        - ( lkx_b_sO((end-ndep+1):end,:,1) * W_(1:Plen_b,:) );
                    % V of VN=[V0W_sO ; VdNxuW_sNO] is an orthonormal basis for parameters of all vfields over SN, U an orthogonal basis for T_O SN
                    sNO_basis_out = Tspc_package( [V0W_sO ; VdNxuW_sNO],tO_ );
                    % [sNO_basis_out,sNO_PCA_out] = Tspc_package( [V0W_sO ; VdNxuW_sNO],tO_ );
                    sNO_basis_out.WYmu = W_ * sNO_basis_out.V(:,1:nvar_);
                    % theta_WV_sO is ntheta x nvar parameter matrix, cols lincom lambda fcns, yield vfield coeffs
                    sNO_basis_out.theta_WV_sO = ( sNO_basis_out.Uscl(:)' / sNO_basis_out.s(1)) .* sNO_basis_out.WYmu;
                    % VNspc_sO is B x nvar matrix, columns are vfield coeffs, span tangent space of SN at sNO, lie algebra at origin
                    [sNO_basis_out.Vspc_sO,sNO_basis_out.theta_tns_WV_sO] = complete_sO_Tspc_image( sNO_basis_out.theta_WV_sO );
                end
            end
            function [Nsvd_out,sO_NT_basis_out,s0O_NT_basis_out,sNO_NT_basis_out] = compute_sO_nontrivial_basis(W_,iPv_,tS_,lvs_,Lam_v_T_ttns_,tO_,lam_sO_,Jltns_xi_)
                [Plen_v,nvar_,ndim_,nobs_] = size(Lam_v_T_ttns_);
                ntheta_v = Plen_v*nvar_;

                lv_b_sO = lam_sO_.lrow_vals(iPv_);
                dkxl_b_sO = lam_sO_.dkxl(1,iPv_);
                lkx_b_sO = lam_sO_.lkx((end-ndep+1):end,iPv_,1);
                % nvar x C matrix, columns span tangent space of S0 at s0O
                V0spc_image = @(thtns_) reshape( ...
                    pagemtimes( thtns_, lv_b_sO(:) ) , [nvar_ size(thtns_,3)] );
                % VdNxuspc_sO is ndep x nvar matrix, columns are vfield coeffs of dNxu in jet space
                VdNxuspc_image = @(thmat_,thtns_) reshape( ...
                    pagemtimes( thtns_((end-ndep+1):end,:,:) , dkxl_b_sO(:) ) , [ndep size(thtns_,3)] ...
                    ) - lkx_b_sO((end-ndep+1):end,:,1) * thmat_(1:Plen_v,:) ;
                % B x nvar matrix, columns are vfield coeffs, span tangent space of SN at sNO, lie algebra at origin
                VNspc_image = @(thmat_,thtns_) deal([V0spc_image(thtns_) ; VdNxuspc_image(thmat_,thtns_)],thtns_);
                % thtns_ is nvar x Plen x nvar, pages are coordinate vfield parameter matrices, act on lambda vectors, yield coeffs
                complete_sO_Tspc_image = @(thmat_) VNspc_image( ...
                    thmat_ , permute(reshape(thmat_,Plen_v,nvar_,size(thmat_,2)),[2 1 3]) ...
                );

                tS_unit = tS_./sqrt( sum(tS_.*tS_,1) );
                % ndim x C x nobs, Lambda matrices projected over candidate vfield space (sample of tangent bundle)
                Mu_S = permute( ...
                    pagemtimes(W_',reshape(Lam_v_T_ttns_,[ntheta_v ndim_ nobs])), ...
                [2 1 3]);
                % ndim x ntheta x nobs, Mu tangent bundle stripped of component in the direction of tvf
                Nu_S = Mu_S-pagemtimes( reshape(tS_unit,[ndim_ 1 nobs_]) , ...
                                        pagemtimes(reshape(tS_unit,[1 ndim_ nobs_]) , Mu_S) );
                % the principle components of this matrix correspond to vector fields not parallel to the tvf everywhere.
                Nsvd_out = Asvd_package( reshape(permute( Nu_S,[2 1 3]), size(W_,2), ndim_*nobs_ )' );

                Theta_N = W_*(Nsvd_out.D/Nsvd_out.s(1));

                [s0O_NT_basis_out,sNO_NT_basis_out] = compute_sO_W_tspc( ...
                    Theta_N, iPv_, lam_sO_, tO_ ...
                );

                Theta_nu = Theta_N*sNO_NT_basis_out.V(:,1:nvar_);
                [Nu_O,Theta_nu_tns] = complete_sO_Tspc_image(Theta_nu);
                a_nu_tO = lsqminnorm(Nu_O,tO_);
                theta_nu_tO = Theta_nu*a_nu_tO;

                sNO_NT_basis_out.Theta_nu = Theta_nu;
                sNO_NT_basis_out.a_nu_tO = a_nu_tO;
                sNO_NT_basis_out.theta_nu_tO = theta_nu_tO;

                sO_NT_basis_out = sNO_NT_basis_out;

            end

            tic0 = tic;
            inds_P_RN1 = 1:Plen_RN1;
            Rmat_N1 = (reshape(Renc_tns_RN1,ntheta_RN1,ndep_N1*nobs))'; % nobs*Qhat by ntheta, where Qhat = kor*Q
            Rsvd_N1 = Asvd_package(normalize_Renc(Rmat_N1)); % global R1 svd
            mrow_R_net = size(Rmat_N1,1);
            %{
                Computed the svd of an "overdetermined" R matrix, which we actually expect to be rank deficient.
                Letting KR denote an orthonormal basis for the kernal of R, which has dimension kappa,
                for each k = 1, ... , N, there exists a subspace of ker(R) which has \theta \in \R^C also satisfying
                    D^k th = vdkxu - pr1 (vdkm1xu) = l . th_dkxu - dx (l . th_dkm1xu) + dkxu dx (l.th_x) = 0,
                that is, first order prolongation adherence for each k = 1, ..., N, automatically satisfied when N=1.
                Intersect these nullspaces for unique trivial vector field model.
            %}
            nsvd = 1;
            if (kor>1)
                DprN_mat = reshape(DprN_T_ttns,ntheta_RN1, ndep*(kor-1)*nobs )';
                DprN_svd = Asvd_package(normalize_DprN(DprN_mat));
                Rmat_N1_net = [ Rmat_N1 ; DprN_mat ];
                Rsvd_N1_net = Asvd_package(normalize_Renc(Rmat_N1_net));
                mrow_R_net = size(Rmat_N1_net,1); % redefine column space dimension to that of concatenated, transposed row space
                Rtns_T_net = reshape(Rmat_N1_net',Plen_RN1,nvar_N1,mrow_R_net);
                nsvd = 3;
            else
                Rmat_N1_net = Rmat_N1; % automatically equivalent
                % [Rsvd_N1_net,DprN_svd] = deal(Rsvd_N1);
                Rsvd_N1_net = Rsvd_N1;
                Rtns_T_net = reshape(Rmat_N1_net',Plen_RN1,nvar_N1,mrow_R_net);
            end
            toc1 = toc(tic0);
fprintf('(LDsol::model_solspace) Decomposed %d R + DprN matrices in %.2f seconds: %dx%d (r=%d,k=%d) -> %dx%d (r=%d,k=%d,o0=%d,oN=%d,o=%d) \n', ...
            nsvd, toc1, ...
            size(Rmat_N1,1), size(Rmat_N1,2), Rsvd_N1.r, Rsvd_N1.dim-Rsvd_N1.r, ...
            mrow_R_net, Rsvd_N1_net.dim, Rsvd_N1_net.r, Rsvd_N1_net.dim-Rsvd_N1_net.r, ...
            max(sum(fspace_RN1.Pmat(1:nvar,inds_P_RN1),1)), max(sum(fspace_RN1.Pmat((nvar+1):end,inds_P_RN1),1)), ...
            max(sum(fspace_RN1.Pmat(:,inds_P_RN1),1)) );

            vth_RN1_net = comp_vartheta_RN1( Rsvd_N1_net.W , lvs_RN1(inds_P_RN1,:) ); % ntheta by nobs
            tau_uN_RN1_net = comp_tau_uN_RN1( vth_RN1_net , inds_P_RN1 ); % ndep by 2 by nobs
            uNp1_tvf_mat = reshape(cat(1, tau_uN_RN1_net(:,1,:), tau_uN_RN1_net((end-ndep+1):end,2,:)),ndim-1,nobs);
            tauN_S_mat = [ ones(1, nobs) ; uNp1_tvf_mat ];

            %% prepare function space for representation of G kernal
            Pmat_GN1_full = fspace_RN1.Pmat;
            Plen_GN1 = fspace_RN1.Plen;
            ntheta_GN1 = fspace_RN1.ntheta; % ntheta x 1
            LamN_tns_GN1 = LamN_tns_RN1; % ndim x ntheta x nobs

            %% generate G matrices, kernel vector fields satisfy infinitesimal criterion
            J_tau_u_RN1 = nan(ndep_N1,ndim_N1,2,nobs);
            JF_N1 = nan(ndep,ndim,nobs);
            Jdxl_N1_tns = zeros(ndim_N1,Plen_GN1,nobs);
            Jl1x_N1_ttns = zeros(ndep_N1,Plen_GN1,ndim_N1,nobs);
            Gtns_N1 = nan(ndep,ntheta_GN1,nobs);
            Btns_t0_N1 = nan(nvar_N1,ntheta_GN1,nobs);
            Btns_tdxu_N1 = nan(ndep_N1,ntheta_GN1,nobs);
            tic0 = tic;
            for iobs = 1:nobs
                J_tau_u_RN1(:,:,:,iobs) = sols(iobs).lamRN1.J_tau_uN( vth_RN1_net(:,iobs) );
                % holds due to first order ratio condition and DprN enforcement
                JF_N1(:,:,iobs) = ...
                    [ J_tau_u_RN1((end-ndep+1):end,1:nvar_N1,1,iobs) , -eye(ndep) ];
                % induced inf criterion
                Gtns_N1(:,:,iobs) =  JF_N1(:,:,iobs) ...
                                    * [LamN_tns_GN1(1:nvar_N1,:,iobs) ; LamN_tns_GN1((end-ndep+1):end,:,iobs)];
                Btns_t0_N1(:,:,iobs) = immerse_l_Lamspc(dxl_RN1(:,iobs)',nvar_N1) ...
                    - [ zeros(1,nvar_N1) ;
                        J_tau_u_RN1(:,1:nvar_N1,1,iobs) ] * LamN_tns_GN1(1:nvar_N1,:,iobs);

                Jdxl_N1_tns(:,:,iobs) = sols(iobs).lamRN1.Jdkxl(:,:,1);
                Jl1x_N1_ttns(:,:,:,iobs) = reshape(sols(iobs).lamRN1.Jlkx(:,:,1,:),[ndep_N1 Plen_GN1 ndim_N1]);

                tN_N1_i = [1 ; reshape(tau_uN_RN1_net(:,:,iobs),[],1)];
                Btns_tdxu_N1(:,:,iobs) = ...
                    [   permute(pagemtimes(tN_N1_i',permute(-Jl1x_N1_ttns(:,:,:,iobs), [3 2 1])),[3 2 1]) , ...
                            immerse_l_Lamspc(tN_N1_i' * Jdxl_N1_tns(:,:,iobs),ndep_N1) ] ...
                    - (J_tau_u_RN1(:,:,2,iobs) * LamN_tns_GN1(:,:,iobs)) ;
            end
            toc1 = toc(tic0);
fprintf('(LDsol::model_solspace) encoded G, %dx%dx%d, in %.2f seconds.\n', ...
            size(Gtns_N1,1), size(Gtns_N1,2), size(Gtns_N1,3),  ...
            toc1);

            Bmat_t0_N1 = (reshape(permute(Btns_t0_N1,[2 1 3]),ntheta_GN1,nvar_N1*nobs))';
            Bsvd_t0_N1 = Asvd_package(normalize_Benc(Bmat_t0_N1));
            Bmat_tdxu_N1 = (reshape(permute(Btns_tdxu_N1,[2 1 3]),ntheta_GN1,ndep_N1*nobs))';
            Bsvd_tdxu_N1 = Asvd_package(normalize_Benc(Bmat_tdxu_N1));

            Bsvd_tN_N1 = Asvd_package(normalize_Benc([ Bmat_t0_N1 ; Bmat_tdxu_N1 ]));

            %% assemble net G matrix, obeying prolongation, as well as optional tvf commutativity constraints
            tic0 = tic;
            [inds_P_GN1_full,inds_P_GN1_net,inds_P_GN1_com] = deal(1:Plen_GN1);
            Gmat_N1 = (reshape(permute(Gtns_N1,[2 1 3]),ntheta_GN1,ndep*nobs))';
            Gsvd_N1 = Asvd_package(normalize_Genc(Gmat_N1));
            nsvd = 1;
            if (kor>1)

                % most general class of vector fields over solution space, likely underdetermined
                Gmat_N1_net_full = [ Gmat_N1 ; DprN_mat ];
                Gsvd_N1_net_full = Asvd_package(normalize_Genc(Gmat_N1_net_full));
                nsvd = nsvd+1;
                [mrow_GN1_net,ntheta_GN1_net] = size(Gmat_N1_net_full);
                % if Gmat_N1_net_full is underdetermined, pare down column space until overdetermined
                if ( mrow_GN1_net<=ntheta_GN1_net )
                    Gtns_T_net = reshape(Gmat_N1_net_full',Plen_GN1,nvar_N1,mrow_GN1_net);
                    [inds_P_GN1_net,ntheta_GN1_net] = pare_colspc(Pmat_GN1_full,mrow_GN1_net);
                    Gmat_N1_net = reshape(Gtns_T_net(inds_P_GN1_net,:,:),ntheta_GN1_net,mrow_GN1_net)';
                    Gsvd_N1_net = Asvd_package( ...
                        normalize_Genc(Gmat_N1_net) ...
                    );
                    nsvd = nsvd+1;
                else
                    Gmat_N1_net = Gmat_N1_net_full;
                    Gsvd_N1_net = Gsvd_N1_net_full;
                end

            else % N=1 => G matrix rows >= R matrix rows, nothing to do
                DprN_mat = Gmat_N1; % reinforce infinitesimal criterion, if no consistent prolongation constraint is in play
                DprN_svd = Gsvd_N1;
                [Gmat_N1_net,Gmat_N1_net_full] = deal(Gmat_N1);
                [Gsvd_N1_net,Gsvd_N1_net_full] = deal(Gsvd_N1);
                [mrow_GN1_net,ntheta_GN1_net] = size(Gmat_N1_net);
                nsvd = nsvd+1;
            end
            toc1 = toc(tic0);
fprintf('(LDsol::model_solspace) Decomposed %d G+DprN matrices in %.2f seconds: %dx%d (r=%d,k=%d) -> %dx%d (r=%d,k=%d,o0=%d,oN=%d,o=%d) \n', ...
            nsvd, toc1, ...
            size(Gmat_N1,1), size(Gmat_N1,2), Gsvd_N1.r, Gsvd_N1.dim - Gsvd_N1.r, ...
            mrow_GN1_net, Gsvd_N1_net.dim, Gsvd_N1_net.r, Gsvd_N1_net.dim - Gsvd_N1_net.r, ...
            max(sum(fspace_RN1.Pmat(1:nvar,inds_P_GN1_net),1)), max(sum(fspace_RN1.Pmat((nvar+1):end,inds_P_GN1_net),1)) , ...
            max(sum(fspace_RN1.Pmat(:,inds_P_GN1_net),1)) );

            %{
                SVDs of Gnet and Gcom reveal kernal vfields of S. The latter are guaranteed to commute with the TVF.

                Using the Gnet kernal vfields, we generate global canonical coordinates.
                Using Gcom kernal vfields, we generate local coordinates suitable for flow transformations, which perturb observed
                integral curves into those passing through arbitrary initial conditions in a neighborhood of the observations
            %}
            % refine the trivial vector field model by intersecting R matrix kernal with Gnet matrix kernal
            RGsvd_N1_net = Asvd_package([ normalize_Renc(Rmat_N1_net) ; normalize_Genc(Gmat_N1) ]);

            %% choose an arbitrary origin for the generation of an intrinsic coordinate system
            [icrv_sO,i_sO_0,i_sO_1] = deal(1,1,2); % mid point between first and second observed solutions on curve 1, w.l.o.g.
            jt_O = LDsol.compute_trivial_Hermite_jet( ...
                [Smat(:,i_sO_0) ; tau_uN_RN1_net((end-ndep+1):end,2,i_sO_0)], ...
                [Smat(:,i_sO_1) ; tau_uN_RN1_net((end-ndep+1):end,2,i_sO_1)], ...
                ndep ...
            );
            s_O0 = [ jt_O.xh ; reshape( jt_O.Amat(1:kor,:)', ndep*(kor), 1 ) ]; % extract fitted base space origin, s_O0
            lamN1_sO = adlam( fspace_N1, s_O0 );
            % pass s_O0 to tvf model, get s_N, N+1'th derivative, and Jacobian at the origin
            [f_O0,dxf_O0,vth_sO,lamRN1_sO,Jtu_sO,JtuN_sO] = comp_f_s0(s_O0,fspace_RN1,RGsvd_N1_net.W,inds_P_RN1);
            s_O = [ s_O0 ; f_O0((end-ndep+1):end) ]; % set the jet space origin as the graph of tvf on s_O0
            sNp1_O = [ s_O ; dxf_O0((end-ndep+1):end) ];
            t_O = [ 1 ; sNp1_O((nvar+1):end) ]; % tvf tangent vector in the N'th jet space at the origin

            % P x kappa+1, non globally constant function parameters w gradient orthogonal to tvf
            Jl_N1_svd.YJlW_Eta_tvf = Jl_N1_svd.YJl * Jl_N1_svd.H_YJl_svd.W;
            % svd of the candidate canonical independent coordinate variable gradients at origin
            Jl_N1_svd.gEta_tvf_YJl_sO_svd = Asvd_package(lamN1_sO.Jl * Jl_N1_svd.YJlW_Eta_tvf );
            % P x nvar, parameters of canonical independent coordinates with orthogonal gradients at the origin
            Jl_N1_svd.Theta_Eta_tvf = Jl_N1_svd.YJlW_Eta_tvf * Jl_N1_svd.gEta_tvf_YJl_sO_svd.V( :, 1:(nvar_N1-1) );
            % 1 x (nvar-1), image of tvf independent coordinates over sO (includes near null space solution)
            Jl_N1_svd.Eta_tvf_sO = lamN1_sO.lrow_vals * Jl_N1_svd.Theta_Eta_tvf;
            % nobs x (nvar-1), image of tvf independent coordinates over S (includes near null space solution)
            Jl_N1_svd.Eta_tvf_S = lvs_N1' * Jl_N1_svd.Theta_Eta_tvf;

            %% validate flow transformation technique
            inds_P_GN1_com = 1:Plen_GN1;
            [flow_pckg,Gsvd_N1_com] = verify_flow_transformation(inds_P_GN1_com,t_O,JtuN_sO,lamRN1_sO);

            function [flow_out,Gc_svd] = verify_flow_transformation(iPv_,tO_,JtuN_sO_,lam_sO_)
                lv_b_sO = lam_sO_.lrow_vals(iPv_);
                dkxl_b_sO = lam_sO_.dkxl(1,iPv_);
                lkx_b_sO = lam_sO_.lkx(:,iPv_,1);
                Jlv_b_sO = lam_sO_.Jl(:,iPv_);
                Jdxlv_b_sO = lam_sO_.Jdkxl(:,iPv_,1);
                Plen_v = length(lv_b_sO(:));
                ntheta_v = nvar_N1*Plen_v;
                Jl1xv_b_sO = reshape(lam_sO_.Jlkx(:,iPv_,1,:),ndep_N1,Plen_v,ndim_N1);

                % nvar x C matrix, columns span tangent space of S0 at s0O
                V0spc_image = @(thtns_) reshape( ...
                    pagemtimes( thtns_, lv_b_sO(:) ) , [nvar_N1 size(thtns_,3)] );
                % VdNxuspc_sO is ndep x nvar matrix, columns are vfield coeffs of dNxu in jet space
                VdNxuspc_image = @(thmat_,thtns_) reshape( ...
                    pagemtimes( thtns_((end-ndep+1):end,:,:) , dkxl_b_sO(:) ) , [ndep size(thtns_,3)] ...
                    ) - lkx_b_sO((end-ndep+1):end,:,1) * thmat_(1:Plen_v,:) ;
                % B x nvar matrix, columns are vfield coeffs, span tangent space of SN at sNO, lie algebra at origin
                VNspc_image = @(thmat_,thtns_) [ V0spc_image(thtns_) ; VdNxuspc_image(thmat_,thtns_) ];
                % thtns_ is nvar x Plen x nvar, pages are coordinate vfield parameter matrices, act on lambda vectors
                compute_sO_Tspc_image = @(thmat_) VNspc_image( ...
                    thmat_ , permute(reshape(thmat_,Plen_v,nvar_N1,size(thmat_,2)),[2 1 3]) ...
                );

                LamN_v_tns_sO = zeros(Plen_v,nvar_N1,ndim_N1);
                LamN_v_tns_sO(:,1,1) = lv_b_sO';
                for idep = 1:ndep_N1
                    LamN_v_tns_sO(:,idep+1,idep+1) = lv_b_sO';
                    LamN_v_tns_sO(:,idep+1,idep+nvar_N1) = dkxl_b_sO';
                    LamN_v_tns_sO(:,1,idep+nvar_N1) = -lkx_b_sO(idep,:)';
                end
                LamN_v_mat_sO = reshape(LamN_v_tns_sO,[ntheta_v ndim_N1])';
                DprN_mat_sO = zeros(ndep*(kor-1),ntheta_v);
                if (kor>1)
                    DprN_tns_sO = zeros(ntheta_v,ndep,kor-1);
                    iiLam0 = nvar + (1:ndep);
                    iiLam1 = iiLam0-nvar+nvar_N1;
                    for k = 2:kor
                        DprN_tns_sO(:,:,k-1) = LamN_v_mat_sO(iiLam1,:)' - LamN_v_mat_sO(iiLam0,:)';
                        iiLam0 = iiLam0 + ndep;
                        iiLam1 = iiLam1 + ndep;
                    end
                    DprN_mat_sO = reshape(DprN_tns_sO,ntheta_v,[])';
                end
                tN_O_unit = tO_ / norm(tO_);
                Jtu0_sO = [ zeros(1,ndim) ; [ JtuN_sO_(:,1:nvar_N1,1) , zeros((nvar_N1-1),ndep) ] ]; % ndep_N1 x nvar_N1
                Jf_sO = Jtu0_sO( (end-ndep+1):end,1:nvar_N1 ); % ndep x nvar_N1, grad of N'th derivatives at sO
                JF_sO = [ Jf_sO , -eye(ndep) ]; % ndep x ndim, JF = ( Jf , -I )
                [JF_sO_svd,U_JF_sO] =  Asvd_package(JF_sO');
                JF_sO_svd.U = U_JF_sO;
                P_JF_sO = eye(ndim) - U_JF_sO*U_JF_sO';

                % ndim x nobs x ndep, pages are unit gradient vectors (normalized row of Jacobian of F)
                JF_S_unit_tns = permute(JF_N1 ./ sqrt(sum(JF_N1.^2,2)),[2 3 1]);
                tauN_S_unit_mat = tauN_S_mat ./ sqrt( sum(tauN_S_mat.*tauN_S_mat,1) );

                lvs_v0_mat = lvs_RN1(iPv_,:); % Plen x nobs
                Jl_v0_tns = permute(Jltns_RN1(:,iPv_,:),[2 1 3]); % Plen x nvar_N1 x nobs
                dxl_vdxu_mat = dxl_RN1(iPv_,:); % Plen x nobs
                Jdxl_vN_tns = permute( Jdxl_N1_tns(:,iPv_,:) , [2 1 3] ); % Plen x ndim_N1 x nobs
                l1x_vdxu_tns = lNx_RN1(:,iPv_,:); % ndep x Plen x nobs
                Jl1x_vN_ttns = Jl1x_N1_ttns(:,iPv_,:,:); % ndep x Plen x ndim_N1 x nobs
                compute_vTh_S = @(thvmat_) [ ...
                    thvmat_' * lvs_v0_mat ; ...
                    thvmat_(:,2:end)' * dxl_vdxu_mat ...
                    - reshape(pagemtimes(l1x_vdxu_tns, thvmat_(:,1)), [ndep_N1 nobs] ) ...
                ];
                LamN_T_vN_ttns = LamN_T_ttns_RN1(iPv_,:,:,:); % Plen x nvar_N1 x ndim_N1 x nobs
                LamN_T_v1_ttns = LamN_T_vN_ttns(:,:,iN1_2_N,:); % Plen x nvar_N1 x ndim_N1 x nobs
                LamN_T_v0_ttns = LamN_T_vN_ttns(:,:,1:nvar_N1,:); % Plen x nvar_N1 x nvar_N1 x nobs
                function [vN_sO_vec vN_S_mat Btns_v0 Btns_vdxu Bmat_v0_sO Bmat_vdxu_sO] = compute_vth_sO_S_data(th_)
                    % vN_sO_vec = compute_sO_Tspc_image( th_(:,1) );

                    Th_v_mat = reshape(th_,[Plen_v nvar_N1]);
                    vN1_sO_vec = [  Th_v_mat' * lv_b_sO(:) ;
                                    Th_v_mat(:,2:end)' * dkxl_b_sO(:) - lkx_b_sO*Th_v_mat(:,1) ];
                    vN_sO_vec = [vN1_sO_vec(1:nvar_N1) ; vN1_sO_vec((end-ndep+1):end) ];

                    v_S_mat = compute_vTh_S(Th_v_mat);
                    vN_S_mat = v_S_mat([1:nvar_N1 (ndim_N1-ndep+1):ndim_N1],:);
                    v_S_tns = reshape(v_S_mat,[ndim_N1 1 nobs]);
                    % Plen x 1, v^(0) ( lambda )
                    de_l_sO_vec = vN1_sO_vec(1:nvar_N1)' * Jlv_b_sO;
                    % Plen x nobs, v^(0) ( lambda )
                    de_l_mat = reshape(pagemtimes(Jl_v0_tns,v_S_tns(1:nvar_N1,1,:)),[Plen_v nobs]);
                    de_Lam_0_sO_tns = zeros(Plen_v,nvar_N1,nvar_N1);
                    de_Lam_0_ttns = zeros(Plen_v,nvar_N1,nvar_N1,nobs);
                    for ivar = 1:nvar_N1
                        de_Lam_0_sO_tns(:,ivar,ivar) = de_l_sO_vec';
                        de_Lam_0_ttns(:,ivar,ivar,:) = de_l_mat;
                    end
                    % keyboard
                    % Plen x 1, v^(1) ( dx lambda )
                    de_dxl_sO_vec = vN1_sO_vec' * Jdxlv_b_sO;
                    % Plen x nobs, v^(1) ( dx lambda )
                    de_dxl_mat = reshape(pagemtimes(Jdxl_vN_tns,v_S_tns),[Plen_v nobs]);
                    de_Lam_dxu_sO_tns = zeros(Plen_v,nvar_N1,ndep_N1);
                    de_Lam_dxu_ttns = zeros(Plen_v,nvar_N1,ndep_N1,nobs);
                    for idep = 1:ndep_N1
                        de_Lam_dxu_sO_tns(:,1,idep) = -reshape(Jl1xv_b_sO(idep,:,:),[Plen_v ndim_N1])*vN1_sO_vec;
                        de_Lam_dxu_sO_tns(:,idep+1,idep) = de_dxl_sO_vec';
                        de_Lam_dxu_ttns(:,1,idep,:) = -pagemtimes( ...
                            reshape(Jl1x_vN_ttns(idep,:,:,:), [Plen_v ndim_N1 nobs]), v_S_tns );
                        de_Lam_dxu_ttns(:,idep+1,idep,:) = de_dxl_mat;
                    end
                    
                    % nvar_N1 x nvar_N1 , J ( v^(0) ), cols are partials
                    Jv_0_sO_mat = Th_v_mat' * Jlv_b_sO' ;
                    % ndep x ndim , J ( v_dxu ), cols are partials
                    Jv_dxu_sO_mat = Th_v_mat(:,2:end)' * Jdxlv_b_sO' ...
                                    - reshape(pagemtimes(Jl1xv_b_sO , Th_v_mat(:,1)),[ndep_N1,ndim_N1]);
                    % nvar_N1 x nvar_N1 x nobs, J ( v^(0) )
                    Jv_0_S_tns = pagemtimes( Th_v_mat' , Jl_v0_tns );
                    % ndep x ndim x nobs, J ( v_dxu )
                    Jv_dxu_S_tns = pagemtimes( Th_v_mat(:,2:end)', Jdxl_vN_tns ) ...
                                    - reshape(sum(Th_v_mat(:,1).*permute(Jl1x_vN_ttns,[2 1 3 4]), 1),[ndep_N1 ndim_N1 nobs]);

                    % nvar_N1 x ntheta, base space v Lie bracket commutativity condition encoded as matrix
                    Bmat_v0_sO = reshape(de_Lam_0_sO_tns,[ntheta_v nvar_N1])' - Jv_0_sO_mat * LamN_v_mat_sO(1:nvar_N1,:);
                    % ndep_N1 x ntheta, jet space v Lie bracket commutativity condition encoded as matrix
                    Bmat_vdxu_sO = reshape(de_Lam_dxu_sO_tns,[ntheta_v ndep_N1])' - Jv_dxu_sO_mat * LamN_v_mat_sO;
                    % nvar_N1 x ntheta x nobs, base space v Lie bracket commutativity condition encoded as matrices
                    Btns_v0 =  ...
                        permute(reshape(de_Lam_0_ttns,[ntheta_v nvar_N1 nobs]),[2 1 3]) ...
                        -pagemtimes(Jv_0_S_tns,permute(reshape(LamN_T_v0_ttns,[ntheta_v nvar_N1 nobs]),[2 1 3]));
                    % ndep_N1 x ntheta x nobs, jet space v Lie bracket commutativity condition encoded as matrices
                    Btns_vdxu = ...
                        permute(reshape(de_Lam_dxu_ttns,[ntheta_v ndep_N1 nobs]),[2 1 3]) ...
                        -pagemtimes(Jv_dxu_S_tns,permute(reshape(LamN_T_vN_ttns,[ntheta_v ndim_N1 nobs]),[2 1 3]));
                end

                function [VNspc_Th_i,VNspc_Th_i0] = compute_transversal_Tspace(Th_,VNi_unit_sO_)
                    VNspc_Th_i0 = compute_sO_Tspc_image(Th_);
                    VNspc_Th_i = P_JF_sO*VNspc_Th_i0; % project away Jacobian components
                    % Gramm-Schmidtt away current basis unit tangent vectors
                    for ibse = 1:size(VNi_unit_sO_,2)
                        VNspc_Th_i = VNspc_Th_i - VNi_unit_sO_(:,ibse) * ( VNi_unit_sO_(:,ibse)' * VNspc_Th_i );
                    end
                end
                function [VNspc_Th_i,VNspc_Th_i0] = compute_orthogonal_Tspace(Th_,VNi_unit_sO_,BNi_sO_)
                    VNspc_Th_i0 = compute_sO_Tspc_image(Th_);

                    VNspc_Th_i = P_JF_sO*VNspc_Th_i0; % project away Jacobian components
                    % Gramm-Schmidtt away current basis unit tangent vectors
                    for ibse = 1:size(VNi_unit_sO_,2)
                        VNspc_Th_i = VNspc_Th_i - VNi_unit_sO_(:,ibse) * ( VNi_unit_sO_(:,ibse)' * VNspc_Th_i );
                    end

                    PJF_VNspc_Th_svd = Asvd_package(VNspc_Th_i);

                end

                theta_v_coords = zeros(ntheta_v,ndep_N1);
                BD_vN_tns = zeros(ntheta_v,ntheta_v,ndep_N1);
                Tspc_Nv_sO_tns = zeros(ndim,ntheta_v,ndep_N1);
                VN_spc_sO = zeros(ndim,nvar_N1);
                VN_spc_S = zeros(ndim,nobs,ndep_N1);

                VN_spc_unit_sO = zeros(ndim,nvar_N1);
                VN_spc_unit_S = zeros(ndim,nobs,nvar_N1);

                % Gn_svd = Gsvd_N1_net;
                Gn_svd = Asvd_package([ ...
                    Gmat_N1_net_full
                ]);
                Mu_S_0 = permute(pagemtimes((Gn_svd.W)',reshape(LamN_T_v1_ttns,[ntheta_v ndim nobs])), [2 1 3]);
                % ndim x ntheta x nobs, Mu tangent bundle stripped of component in the direction of tvf
                Nu_S_0 = Mu_S_0-pagemtimes( reshape(tauN_S_unit_mat,[ndim 1 nobs]) , ...
                                    pagemtimes(reshape(tauN_S_unit_mat,[1 ndim nobs]) , Mu_S_0) );
                % Nu_S_0 = Mu_S_0;
                % the principle components of this matrix correspond to vfields not parallel tvf everywhere
                N_v_svd0 = Asvd_package( reshape(permute( Nu_S_0,[2 1 3]), ntheta_v, ndim*nobs )' );
                % parameters of vector fields not parallel to current basis everywhere
                Theta_N_0 = Gn_svd.W * (N_v_svd0.D / N_v_svd0.s(1));
                %% identify subspace of candidate vector fields transversal to current basis at origin
                [TTspc_0,Tspc_Nv_sO_mat0] = compute_transversal_Tspace(Theta_N_0, []);
                [TTspc_svd_0, U_TTspc_0] = Asvd_package(TTspc_0);
                TTspc_svd_0.U = U_TTspc_0;
                dimY0 = ndep;
                % dimY0 = nvar_N1;
                gspc_sO_svd_0 = Asvd_package( U_JF_sO'*Tspc_Nv_sO_mat0*TTspc_svd_0.V(:,1:dimY0) );
                %% extract null vector field, identify as 1st coordinate vector field
                theta_v_coord0 = Theta_N_0*TTspc_svd_0.V(:,1:dimY0)*gspc_sO_svd_0.V(:,end);
                [VN_spc_sO_0 VN_spc_S_0 Btns_v0_0 Btns_v0_dxu] = compute_vth_sO_S_data(theta_v_coord0);
                % B_v_svd_0 = Asvd_package(normalize_Benc([ ...
                B_v_svd_0 = Asvd_package([ ...
                    reshape(permute(Btns_v0_0,[2 1 3]),ntheta_v,nobs*nvar_N1)' ;
                    reshape(permute(Btns_v0_dxu,[2 1 3]),ntheta_v,nobs*ndep_N1)'
                ]);
                % reshape(permute(Btns_v0_0,[2 1 3]),ntheta_v,nobs*nvar_N1)'  ...
                % reshape(permute(Btns_v0_dxu,[2 1 3]),ntheta_v,nobs*ndep_N1)' ...
                Gnc_svd = Asvd_package([ ...
                    Gn_svd.D/Gn_svd.s(1), ...
                    B_v_svd_0.D/B_v_svd_0.s(1) ...
                ]');

                Gn_bse = N_v_svd0;
                Gn_bse.Theta_N_0 = Theta_N_0;
                Gn_bse.Tspc_Nv_sO_mat0 = Tspc_Nv_sO_mat0;
                Gn_bse.TTspc_svd_0 = TTspc_svd_0;
                Gn_bse.gspc_sO_svd_0 = gspc_sO_svd_0;
                Gn_bse.theta_v_coord0 = theta_v_coord0;

                %% compute SVD of Gc : nullspace consists of vector fields that commute with tvf
                Gc_svd = Asvd_package([ ...
                    Gmat_N1_net_full ;
                    Bmat_t0_N1
                ]);
                % Gmat_N1
                % DprN_mat
                % Gmat_N1_net_full
                % Bmat_t0_N1
                % Bmat_tdxu_N1
                % Gc_svd = Asvd_package([ ...
                %     Gsvd_N1.D/Gsvd_N1.s(1), ...
                %     DprN_svd.D/DprN_svd.s(1), ...
                %     Bsvd_t0_N1.D/Bsvd_t0_N1.s(1) ...
                % ]');
                % Gsvd_N1.D/Gsvd_N1.s(1), ...
                % DprN_svd.D/DprN_svd.s(1), ...
                % Bsvd_tN_N1.D/Bsvd_tN_N1.s(1) ...
                % Bsvd_t0_N1.D/Bsvd_t0_N1.s(1) ...
                % Bsvd_tdxu_N1.D/Bsvd_tdxu_N1.s(1) ...

                G0_svd = Gnc_svd;
                VN_spc_sO(:,1) = VN_spc_sO_0;
                VN_spc_unit_sO(:,1) = VN_spc_sO_0 / norm(VN_spc_sO_0);
                VN_spc_unit_S(:,:,1) = VN_spc_S_0 ./ sqrt(sum(VN_spc_S_0.^2,1));

                % G0_svd = Gc_svd;
                % VN_spc_sO(:,1) = tO_;
                % VN_spc_unit_sO(:,1) = tN_O_unit;
                % VN_spc_unit_S(:,:,1) = tauN_S_unit_mat;

                % Gsvd_N1_net_full.D/Gsvd_N1_net_full.s(1), ...
                %% parameters of candidate vector fields which commute with tvf
                WGc = G0_svd.W;
                WGc_i = WGc;

                function out = build_Tspc_sO_basis(W_,v1_,B1_)

                    [TTspc_i,Tspc_Nv_sO_tns(:,:,ivec)] = compute_orthogonal_Tspace(W_,VN_spc_unit_sO(:,1:ivec));


                    out = 0
                end

                for ivec = 1:ndep_N1
                    % ndim x ntheta x nobs, Lambda matrices projected over candidate vfield space (sample of tangent bundle)
                    Mu_S_i = permute(pagemtimes(WGc_i',reshape(LamN_T_v1_ttns,[ntheta_v ndim nobs])), [2 1 3]);
                    % % ndim x ntheta x nobs, Mu tangent bundle stripped of component in the direction of tvf
                    % Nu_S_i = Mu_S_i-pagemtimes( reshape(tauN_S_unit_mat,[ndim 1 nobs]) , ...
                    %                     pagemtimes(reshape(tauN_S_unit_mat,[1 ndim nobs]) , Mu_S_i) );
                    Nu_S_i = Mu_S_i;
                    % the principle components of this matrix correspond to vfield basis candidates not parallel to tvf everywhere
                    N_v_svds(ivec) = Asvd_package( reshape(permute( Nu_S_i,[2 1 3]), ntheta_v, ndim*nobs )' );
                    % parameters of vector fields not parallel to current basis everywhere
                    Theta_N_i = WGc_i * (N_v_svds(ivec).D / N_v_svds(ivec).s(1));
                    %% identify subspace of candidate vector fields transversal to current basis at origin

                    % [TTspc_i,Tspc_Nv_sO_tns(:,:,ivec)] = compute_transversal_Tspace(Theta_N_i,[]);
                    [TTspc_i,Tspc_Nv_sO_tns(:,:,ivec)] = compute_transversal_Tspace(Theta_N_i,VN_spc_unit_sO(:,1:ivec));
                    [TTspc_svd_i, U_TTspc_i] = Asvd_package(TTspc_i);
                    TTspc_svd_i.U = U_TTspc_i;
                    TTspc_svds(ivec) = TTspc_svd_i;
                    dimYi = min([nvar_N1, ndep+ivec]);
                    % dimYi = nvar_N1;
                    gspc_sO_svds(ivec) = Asvd_package( [ U_JF_sO' ;
                                VN_spc_unit_sO(:,1:ivec)']*Tspc_Nv_sO_tns(:,:,ivec)*TTspc_svd_i.V(:,1:dimYi) );
                    %% extract null vector field, identify as next coordinate vector field
                    theta_v_coords(:,ivec) = Theta_N_i*TTspc_svd_i.V(:,1:dimYi)*gspc_sO_svds(ivec).V(:,end);

                    [VN_spc_sO(:,ivec+1) VN_spc_S(:,:,ivec) Btns_vi_0 Btns_vi_dxu] = ...
                        compute_vth_sO_S_data(theta_v_coords(:,ivec));
                    VN_spc_unit_sO(:,ivec+1) = VN_spc_sO(:,ivec+1) / norm(VN_spc_sO(:,ivec+1));
                    VN_spc_unit_S(:,:,ivec+1) = VN_spc_S(:,:,ivec) ./ sqrt(sum(VN_spc_S(:,:,ivec).^2,1));
                    % B_v_svds(ivec) = Asvd_package(normalize_Benc([ ...
                    B_v_svds(ivec) = Asvd_package([ ...
                        reshape(permute(Btns_vi_dxu,[2 1 3]),ntheta_v,nobs*ndep_N1)' ;
                        reshape(permute(Btns_vi_0,[2 1 3]),ntheta_v,nobs*nvar_N1)'
                    ]);
                    BD_vN_tns(:,:,ivec) = B_v_svds(ivec).D / B_v_svds(ivec).s(1);
                    if (ivec>1)
                        BD_vnet_svd = Asvd_package(reshape(BD_vN_tns(:,:,1:ivec),ntheta_v,[])');
                    else
                        BD_vnet_svd = B_v_svds(ivec);
                    end
                    Gc_v_svds(ivec) = Asvd_package([ ...
                        G0_svd.D/G0_svd.s(1) , ...
                        BD_vnet_svd.D/BD_vnet_svd.s(1) ...
                    ]');
                    WGc_i = Gc_v_svds(ivec).W;
                end

                flow_out = struct( ...
                    'theta_v_coords', theta_v_coords, ...
                    'BD_vN_tns', BD_vN_tns, ...
                    'VN_spc_sO', VN_spc_sO, ...
                    'VN_spc_unit_sO', VN_spc_unit_sO, ...
                    'VN_spc_S', VN_spc_S, ...
                    'Tspc_Nv_sO_tns', Tspc_Nv_sO_tns ...
                );
                flow_out.Gn_bse = Gn_bse;
                flow_out.B_v_svds = B_v_svds;
                flow_out.N_v_svds = N_v_svds;
                flow_out.TTspc_svds = TTspc_svds;
                flow_out.gspc_sO_svds = gspc_sO_svds;
                flow_out.Gc_v_svds = Gc_v_svds;
                flow_out.JF_sO_svd = JF_sO_svd;
            end

            %% compute non-trivial vector field bases
            [~,Mcom_sO_basis] = compute_sO_W_tspc(Gsvd_N1_com.W,inds_P_GN1_com,lamRN1_sO,t_O);
            [~,Mnet_sO_basis] = compute_sO_W_tspc(Gsvd_N1_net.W,inds_P_GN1_net,lamRN1_sO,t_O);

            [Nsvd_N1_com,Ncom_sO_basis,~,~] = compute_sO_nontrivial_basis( ...
                Gsvd_N1_com.W,inds_P_GN1_com,tauN_S_mat,lvs_RN1(inds_P_GN1_com,:),LamN1_T_ttns(inds_P_GN1_com,:,:,:),t_O,lamRN1_sO,Jltns_N1 ...
            );
            [Nsvd_N1_net,Nnet_sO_basis,~,~] = compute_sO_nontrivial_basis( ...
                Gsvd_N1_net.W,inds_P_GN1_net,tauN_S_mat,lvs_RN1(inds_P_GN1_net,:),LamN1_T_ttns(inds_P_GN1_net,:,:,:),t_O,lamRN1_sO,Jltns_N1 ...
            );

            %% bonus computations

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
                    for k = 1:kor
                        fprintf('   [%.2e %.2e %.2e %.2e]', ...
                            err_mat_min(q,k), err_mat_med(q,k), err_mat_avg(q,k), err_mat_max(q,k)  );
                    end
                    fprintf('\n');
                end
            end
            comp_err_msr = @(u_) abs((u_-dxuNtns)./dxuNtns);

            tau_u_RN1_net = reshape(tau_uN_RN1_net(:,1,:),ndep,kor,[]);
            err_tau_u_RN1_net = tau_u_RN1_net-dxuNtns;
            % err_stats( 'tau_u_RN1_net', comp_err_msr(tau_u_RN1_net) );
            err_stats( 'tau_u_RN1_net', abs(err_tau_u_RN1_net./dxuNtns) );

            % keyboard
            mod_out = struct( ...
                'Smat', Smat, ...
                'ipts_crv', ipts_crv, ...
                'Smat_N1', Smat_N1, ...
                'Hmat_N1', Hmat_N1, ...
                'Hmat_0', Hmat_0, ...
                'Hmat_RN1', Hmat_RN1, ...
                'LamN_tns_RN1', LamN_tns_RN1, ...
                'inds_P_RN1', inds_P_RN1, ...
                'Rmat_N1', Rmat_N1, ...
                'inds_P_GN1_net', inds_P_GN1_net, ...
                'inds_P_GN1_com', inds_P_GN1_com, ...
                'Gmat_N1', Gmat_N1, ...
                'vth_RN1_net', vth_RN1_net, ...
                'tau_uN_RN1_net', tau_uN_RN1_net, ...
                'J_tau_u_RN1', J_tau_u_RN1, ...
                'JF_N1', JF_N1 ...
            );
            mod_out.Sobs = Sobs_;
            mod_out.fspace_N1 = fspace_N1;
            mod_out.fspace_0 = fspace_0;
            mod_out.fspace_RN1 = fspace_RN1;
            mod_out.sols = sols;

            mod_out.Jl_N1_svd = Jl_N1_svd;

            mod_out.H_N1_svd = H_N1_svd;
            mod_out.H_0_svd = H_0_svd;
            mod_out.H_RN1_svd = H_RN1_svd;
            mod_out.DprN_svd = DprN_svd;

            mod_out.Rsvd_N1 = Rsvd_N1;
            mod_out.Rsvd_N1_net = Rsvd_N1_net;

            mod_out.err_tau_u_tvf = err_tau_u_RN1_net;

            mod_out.jt_O = jt_O;

            mod_out.f_O0 = f_O0;
            mod_out.vth_sO = vth_sO;
            mod_out.lamRN1_sO = lamRN1_sO;
            mod_out.Jtu_sO = Jtu_sO;
            mod_out.JtuN_sO = JtuN_sO;
            mod_out.lamN1_sO = lamN1_sO;

            mod_out.icrv_sO = icrv_sO;
            mod_out.s_O = s_O;
            mod_out.sNp1_O = sNp1_O;
            mod_out.t_O = t_O;

            mod_out.Bsvd_t0_N1 = Bsvd_t0_N1;
            mod_out.Bsvd_tdxu_N1 = Bsvd_tdxu_N1;
            mod_out.Bsvd_tN_N1 = Bsvd_tN_N1;

            mod_out.Gsvd_N1 = Gsvd_N1;

            mod_out.Gsvd_N1_net = Gsvd_N1_net;
            mod_out.Gsvd_N1_com = Gsvd_N1_com;

            mod_out.Nsvd_N1_net = Nsvd_N1_net;
            mod_out.Nsvd_N1_com = Nsvd_N1_com;

            mod_out.Mnet_sO_basis = Mnet_sO_basis;
            mod_out.Mcom_sO_basis = Mcom_sO_basis;
            mod_out.Nnet_sO_basis = Nnet_sO_basis;
            mod_out.Ncom_sO_basis = Ncom_sO_basis;

            mod_out.flow_pckg = flow_pckg;
        end
        function jt_out = compute_trivial_Hermite_jet(s0_,s1_,ndep_)
            if (s0_(1) > s1_(1)) % enforce trivial flow in the positive direction
                s0 = s1_;
                s1 = s0_;
            else
                s0 = s0_;
                s1 = s1_;
            end
            x0 = s0(1);
            x1 = s1(1);
            u0mat = reshape(s0(2:end),ndep_,[]);
            u1mat = reshape(s1(2:end),ndep_,[]);
            kor = size(u0mat,2)-1;
            xh = 0.5*(x0+x1);

            Jor = 2*(kor+1) - 1;
            Jorp1 = Jor+1;
            Jorp1_h = Jorp1/2;

            p_z2J = 0:Jor;
            pm_Jor = (-1).^p_z2J;
            dh1_Jor = (x1-xh).^p_z2J; % strictly positive
            dh0_Jor = dh1_Jor.*pm_Jor; % alternating sign

            [V0,V1] = deal( zeros(Jorp1_h , Jorp1) );
            fJ = [1, cumprod( p_z2J(2:end) )];
            ch = fJ.^(-1);
            V1(1,:) = ch.*dh1_Jor(1:(end-1+1)); % = ch.*(d).^pdi
            V0(1,:) = ch.*dh0_Jor(1:(end-1+1)); % = ch.*(-d).^pdi
            for i = 2:Jorp1_h % first derivative onwards
                ch(i:end) = ch(i:end) .* p_z2J( 2:(end-(i-1)+1) );
                V1(i,i:end) = ch(i:end) .* dh1_Jor( 1:(end-i+1) );
                V0(i,i:end) = ch(i:end) .* dh0_Jor( 1:(end-i+1) );
            end
            Vmat = [ V0 ; V1 ];
            Amat = nan(Jorp1,ndep_);
            Umat = nan(Jorp1,ndep_);
            for i = 1:ndep_
                Umat(:,i) = [ u0mat(i,:)' ; u1mat(i,:)' ];
                % Amat(:,i) = linsolve(Vmat, Umat(:,i) ); % full rank, possibly ill conditioned
                Amat(:,i) = lsqminnorm(Vmat, Umat(:,i)); % Tikhinov regularized when rank deficient
            end

            jt_out = struct( ...
                'xh', xh, ...
                'fJ', fJ, ...
                'Vmat', Vmat, ...
                'Umat', Umat, ...
                'Amat', Amat ...
            );
        end
    end
    methods
        function obj = LDsol(xu_)
            obj.xu = xu_(:);
        end
    end
end
