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
            bor_N1 = 1;
            [Plen_N1, Pmat_N1, ~] = ldaux.count_set_P_len(bor_N1,ndep_N1+1);
            while ((Plen_N1 < nobs) && (bor_N1 < 10)) % order 10 should be more than enough
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
                bor_0 = min([10 , max([ 1 , floor( ((nobs.*ndep)./nvar).^(1./nvar) - 1  ) ])]);
                [Plen_0, Pmat_0, ~] = ldaux.count_set_P_len(bor_0,ndep+1);

                fspace_0 = fmap_;
                fspace_0.Omap_A = fspace_0.Omap_a .* eye(length(fspace_0.Omap_a));
                fspace_0.bor = bor_0;
                fspace_0.Plen = Plen_0;
                fspace_0.Pmat = Pmat_0;
                fspace_0.ntheta = (1+ndep)*Plen_0;
                fspace_0 = adlam.init_fspace_family(fspace_0);

                fspace_RN1 = fspace_N1;
                ord_i = max(reshape(Pmat_N1((nvar+1):end,:),[],1));
                mrow_RN1_DprN = nobs*(ndep*( 2*kor - 1 )); % = nobs*( ndep*kor + ndep*(kor-1) )
                while ( (mrow_RN1_DprN <= fspace_RN1.ntheta)&&(ord_i>1) )
                    Pmat_i = fspace_RN1.Pmat;
                    ord_i = ord_i - 1;
                    fspace_RN1.Pmat = Pmat_i( :, sum(Pmat_i((nvar+1):end,:),1) <= ord_i );
                    fspace_RN1.Plen = size(fspace_RN1.Pmat,2);
                    fspace_RN1.ntheta = nvar_N1*(fspace_RN1.Plen);
                end
                if (mrow_RN1_DprN <= fspace_RN1.ntheta)
                    fprintf('(LDsol::model_solspace) honestly you should get some more points (mrow_RN1_DprN = %d, ntheta=%d) \n', ...
                    mrow_RN1_DprN, fspace_RN1.ntheta );
                end
            else
                [bor_0,Plen_0,Pmat_0] = deal(bor_N1,Plen_N1,Pmat_N1);
                [fspace_0,fspace_RN1] = deal(fspace_N1); % all identical in the case of N = 1
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
            fspace_RN1.imm_l = @(lambda_) immerse_lambda_RN1(lambda_);

            Hmat_N1 = nan(nobs,Plen_N1);
            Hmat_0 = nan(nobs,Plen_0);
            Hmat_RN1 = nan(nobs,Plen_N1);
            lvs_N1 = nan(Plen_N1,nobs);
            Jltns_N1 = nan(nvar_N1,Plen_N1,nobs);
            lvs_RN1 = nan(Plen_RN1,nobs);
            dxl_RN1 = nan(Plen_RN1,nobs);
            lNx_RN1 = nan(ndep,Plen_RN1,nobs);
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
                dxl_RN1(:,iobs) = dxl_RN1_i;
                lNx_RN1(:,:,iobs) = lRN1_i.lkx((end-ndep+1):end,:,1);
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
            function [f_s0_out, dxf_s0_out, vth_out, lam_out] = comp_f_s0(s0_,fspc_,W_,iP_)
                nvar_ = length(s0_(:));
                Plen_b = size(fspc_.Pmat(:,iP_),2);
                lam_out = adlam( fspc_, s0_ );
                vx_W = ( W_(1:Plen_b,:) )' * reshape( lam_out.lrow_vals(iP_) , [], 1 );
                vth_out = W_ * (vx_W./sum(vx_W.*vx_W,1));
                [vTh_x,vTh_u] = deal( vth_out(1:Plen_b) , reshape(vth_out((Plen_b+1):end),Plen_b,[]) );
                f_s0_out = vTh_u' * (lam_out.lrow_vals(iP_))';
                lam_out = lam_out.prolong_jet_space( reshape(f_s0_out,[],1) );
                dxf_s0_out = vTh_u'*lam_out.dkxl(1,iP_)' - lam_out.lkx(:,iP_,1)*vTh_x;
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
                % if (nargout==2)
                %     [sO_PCA_out,U_PCA_sO] = Asvd_package( VW_sO_ - mean(VW_sO_,2) );
                %     sO_PCA_out.U = U_PCA_sO;
                % end

                % decompose section of tangent space orthogonal to tvf tangent space
                tO_B_unit = tO_B/mag_tO_B;
                % [sO_nTVF,U_V_nTVF_sO] = Asvd_package(( eye(B_in) - tO_B_unit*(tO_B_unit') )*VW_sO_);
                % sO_nTVF.U = U_V_nTVF_sO;
                % sO_nTVF.Tspc_image = VW_sO_ * (( sO_nTVF.s(:)/sO_nTVF.s(1) )' .* sO_nTVF.V);
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
                V0W_sO = reshape(pagemtimes( Wtns , lv_b_sO' ), nvar_, ntheta_b);
                % VdNxuW_sNO is ndep x ntheta, cols are N'th derivative vector field coefficients in N'th jet space
                VdNxuW_sNO = reshape( pagemtimes( Wtns((end-ndep+1):end,:,:), dkxl_b_sO(1,:)' ), ndep, ntheta_b )  ...
                    - ( lkx_b_sO((end-ndep+1):end,:,1) * W_(1:Plen_b,:) );

                % V of V0W_sO is an orthonormal basis for parameters of all vfields over base space, U an orthogonal basis for T_sO S0
                s0O_basis_out = Tspc_package( V0W_sO,tO_ );
                % [s0O_basis_out,s0O_PCA_out] = Tspc_package( V0W_sO,tO_ );
                s0O_basis_out.WYmu = W_ * s0O_basis_out.V(:,1:nvar_);
                % theta_WV_sO is ntheta x nvar parameter matrix, cols lincom lambda fcns, yield vfield coeffs
                s0O_basis_out.theta_WV_sO = ( s0O_basis_out.Uscl(:)' / s0O_basis_out.s(1)) .* s0O_basis_out.WYmu;
                % VNspc_sO is B x nvar matrix, columns are vfield coeffs, span tangent space of SN at sNO, lie algebra at origin
                [s0O_basis_out.Vspc_sO,s0O_basis_out.theta_tns_WV_sO] = complete_sO_Tspc_image( s0O_basis_out.theta_WV_sO );

                if (nargout == 2)
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
            function [sO_coords_out,s0O_basis_out,sNO_basis_out] = compute_sO_coords(Gsvd_,iPv_,lvs_,lam_sO_,tO_,Jltns_xi_,Lam_v_T_ttns_,Gmat_)
                [nvar_,Plen_xi,nobs_] = size(Jltns_xi_);
                Plen_v = size(lvs_,1);
                ntheta_v = Plen_v*nvar_;

                % [s0O_basis_out,sNO_basis_out] = compute_sO_W_tspc(Gsvd_.W,iP_,lam_sO_,tO_);
                % sO_coords_out = s0O_basis_out;
                % sO_coords_out = sNO_basis_out;
                [sO_coords_out,s0O_basis_out,sNO_basis_out] = compute_sO_basis(Gsvd_.W,iPv_,lvs_,lam_sO_,tO_);


                % nvar x nvar x nobs, page rows are transversal tangent vectors evaluated at observed solutions
                Vspc_S = permute( ...
                    reshape(pagemtimes( sO_coords_out.theta_tns_WV_sO , lvs_), nvar_, nobs_, nvar_), ...
                [3 1 2]);
                % nvar x Plen x nobs, pages are directional derivatives wrt coordinate vfields at s^(N-1) |_j
                H_LamTheta_S = pagemtimes( Vspc_S , Jltns_xi_ );
                % svd of nvar*nobs x Plen net H matrix. Kernel consists of globally constant functions (e.g. f(x,u) = 1)
                Hsvd_LamTheta_S = Asvd_package( ...
                    reshape(permute(H_LamTheta_S,[2 1 3]),Plen_xi,nvar_*nobs_)' ...
                );
                % normalize_Henc(reshape(permute(H_LamTheta_S,[2 1 3]),Plen_xi,nvar_*nobs_)') ...
                % Plen x r_H, row space of global H svd over orthogonal vfield basis (at s_O), non globally constant fcns
                YH = Hsvd_LamTheta_S.V(:,1:Hsvd_LamTheta_S.r);
                % nobs x Plen x nvar, pages are Hmats (directional derivatives) over S = { s|_j } of each coordinate vfield
                H_LamTheta_S = permute(H_LamTheta_S, [3 2 1]);

                theta_xicoords_TspcO = zeros(Plen_xi,nvar_);
                err_xicoords_TspcO = nan(nobs_,nvar_);
                %% compute dependent variables of orthogonal tangent space basis at the origin
                for ixi = 1:nvar_
                    theta_xicoords_TspcO(:,ixi) = YH * lsqminnorm(H_LamTheta_S(:,:,ixi)*YH, ones(nobs_,1));
                    err_xicoords_TspcO(:,ixi) = H_LamTheta_S(:,:,ixi)*theta_xicoords_TspcO(:,ixi) - ones(nobs_,1);
                end
                JXi_sO_TspcO = lamN1_sO.Jl * theta_xicoords_TspcO;
                JXi_svd_sO_TspcO = Asvd_package(JXi_sO_TspcO);
                theta_svd_xicoords_TspcO = Asvd_package(theta_xicoords_TspcO');

                %% generate functionally independent canonical coordinate system
                Theta_vcoords = zeros(Plen_v,nvar_,nvar_);
                theta_xicoords = zeros(Plen_xi,nvar_);
                err_xicoords = nan(nobs_,nvar_);

                % parameters associated with principle coordinate vfield
                Theta_vcoords(:,:,1) = reshape(sO_coords_out.theta_WV_sO(:,1),Plen_v,nvar_);
                % directional derivatives wrt principle coordinate vfield
                Hxi_i = H_LamTheta_S(:,:,1);
                % restricted to non globally constant functions
                Hxi_i_YH = Hxi_i*YH;
                % svd of Hxi_i (nobs x r_H) : directional derivatives of non globally constant fcns, wrt i'th coord vfield
                [Hxi_i_svd, U_Hxi_i] = Asvd_package( normalize_Henc(Hxi_i_YH) );
                Hxi_i_svd.U = U_Hxi_i;

                % % % r_H x 1, regularized least squares solution to Hxi_i * YH * a_xi = 1
                % % Hxi_i_svd.a_xi = lsqminnorm(Hxi_i_YH,ones(nobs_,1));
                % % r_H x 1, restricted, regularized least squares solution to Hxi_i * YH * YHYH aa_xi = 1
                % Hxi_i_svd.a_xi = Hxi_i_svd.V(:,1:(Hxi_i_svd.r)) * lsqminnorm(Hxi_i_YH*Hxi_i_svd.V(:,1:(Hxi_i_svd.r)),ones(nobs_,1));
                % Plen x 1, linear combination of YH almost surely satisfying Hxi_i * theta_xi = 1 everywhere
                % theta_xicoords(:,1) = YH * Hxi_i_svd.a_xi;
                theta_xicoords(:,1) = YH * Hxi_i_svd.V(:,1:(Hxi_i_svd.r))  ...
                                            * lsqminnorm(Hxi_i_YH*Hxi_i_svd.V(:,1:(Hxi_i_svd.r)),ones(nobs_,1));
                % theta_xicoords(:,1) = YH * lsqminnorm(Hxi_i_YH,ones(nobs_,1));

                err_xicoords(:,1) = Hxi_i*theta_xicoords(:,1) - ones(nobs_,1);

                H_Xi_svds(1) = Hxi_i_svd;
                G_Xi_svds(1) = Gsvd_; % initial G matrix is just input infinitesimal criterion
                GnetXi_i = [];

                % Plen x 1, linear combination of YH almost surely satisfying Hxi_i * theta_xi = 1 everywhere
                Lam_v_tns_ = permute(reshape(Lam_v_T_ttns_(iPv_,:,1:nvar_,:),ntheta_v,nvar_,nobs_),[2 1 3]);
                % diagonalize canonical dependent variables
                for ixi = 2:nvar_
                    % GXi is an accumulating matrix with rows corresponding to gradients of identified dependent coords over observations
                    GnetXi_i = [ GnetXi_i ; ...
                                 reshape( pagemtimes( ...
                                    reshape( pagemtimes(Jltns_xi_,theta_xicoords(:,ixi-1)) ,1,nvar_,nobs_), ...
                                    Lam_v_tns_ ), ntheta_v,nobs_ )' ];
                    % net infinitesimal criterion wrt xi gradients (and standard infinitesimal criterion)
                    % G_Xi_svds(ixi) = Asvd_package([ Gsvd_.D' ; normalize_Genc(GnetXi_i) ]);
                    G_Xi_svds(ixi) = Asvd_package(normalize_Genc([Gmat_ ; GnetXi_i]));

                    GnetXi_i_sO_basis = compute_sO_basis(G_Xi_svds(ixi).W,iPv_,lvs_,lam_sO_,tO_);
                    Theta_vcoords(:,:,ixi) = reshape(GnetXi_i_sO_basis.theta_WV_sO(:,1),Plen_v,nvar_);
                    % nvar x nvar x nobs, page rows are transversal tangent vectors evaluated at observed solutions
                    V1_GnetXi_i_S = reshape(Theta_vcoords(:,:,ixi)' * lvs_,[1 nvar_ nobs_]);

                    % nobs x Plen, rows are directional derivatives wrt i'th coordinate vfield
                    Hxi_i = reshape(pagemtimes( V1_GnetXi_i_S(1,:,:) , Jltns_xi_ ), Plen_xi, nobs_)';
                    % restricted to non globally constant functions
                    Hxi_i_YH = Hxi_i*YH;
                    % svd of Hxi_i (nobs x r_H) : directional derivatives of non globally constant fcns, wrt i'th coord vfield
                    [Hxi_i_svd, U_Hxi_i] = Asvd_package( normalize_Henc(Hxi_i_YH) );
                    Hxi_i_svd.U = U_Hxi_i;

                    % % % r_H x 1, regularized least squares solution to Hxi_i * YH * a_xi = 1
                    % % Hxi_i_svd.a_xi = lsqminnorm(Hxi_i_YH,ones(nobs_,1));
                    % % r_H x 1, restricted, regularized least squares solution to Hxi_i * YH * YHYH aa_xi = 1
                    % Hxi_i_svd.a_xi = Hxi_i_svd.V(:,1:(Hxi_i_svd.r)) * lsqminnorm(Hxi_i_YH*Hxi_i_svd.V(:,1:(Hxi_i_svd.r)),ones(nobs_,1));
                    % Plen x 1, linear combination of YH almost surely satisfying Hxi_i * theta_xi = 1 everywhere
                    % theta_xicoords(:,ixi) = YH * Hxi_i_svd.a_xi;
                    theta_xicoords(:,ixi) = YH * Hxi_i_svd.V(:,1:(Hxi_i_svd.r)) ...
                                                * lsqminnorm(Hxi_i_YH*Hxi_i_svd.V(:,1:(Hxi_i_svd.r)),ones(nobs_,1));
                    % theta_xicoords(:,ixi) = YH * lsqminnorm(Hxi_i_YH,ones(nobs_,1));

                    err_xicoords(:,ixi) = Hxi_i*theta_xicoords(:,ixi) - ones(nobs_,1);

                    H_Xi_svds(ixi) = Hxi_i_svd;
                end
                sO_coords_out.Hsvd_LamTheta_S = Hsvd_LamTheta_S;
                sO_coords_out.H_LamTheta_S = H_LamTheta_S;
                sO_coords_out.YH = YH;
                sO_coords_out.theta_xicoords_TspcO = theta_xicoords_TspcO;
                sO_coords_out.theta_svd_xicoords_TspcO = theta_svd_xicoords_TspcO
                sO_coords_out.err_xicoords_TspcO = err_xicoords_TspcO;
                sO_coords_out.JXi_sO_TspcO = JXi_sO_TspcO;
                sO_coords_out.JXi_svd_sO_TspcO = JXi_svd_sO_TspcO;
                sO_coords_out.Theta_vcoords = Theta_vcoords;
                sO_coords_out.theta_xicoords = theta_xicoords;
                sO_coords_out.err_xicoords = err_xicoords;
                sO_coords_out.H_Xi_svds = H_Xi_svds;
                sO_coords_out.G_Xi_svds = G_Xi_svds;
            end
            function [sO_basis_out,s0O_basis_out,sNO_basis_out] = compute_sO_basis(W_,iP_,lvs_,lam_sO_,tO_)
                [s0O_basis_out,sNO_basis_out] = compute_sO_W_tspc(W_,iP_,lam_sO_,tO_);

                % % nvar x nobs x nvar, evaluated vector field basis (pages) over base space at each observed solution
                % s0O_basis_out.LamTheta_S = reshape(pagemtimes( s0O_basis_out.theta_tns_WV_sO , lvs_ ), nvar_N1, nobs, nvar_N1);
                % % nvar x nobs x nvar, evaluated vector field basis (pages) over base space at each observed solution
                % sNO_basis_out.LamTheta_S = reshape(pagemtimes( sNO_basis_out.theta_tns_WV_sO , lvs_ ), nvar_N1, nobs, nvar_N1);

                %% choose either base space or jet space orthonormal parameter basis
                % sO_basis_out = s0O_basis_out;
                sO_basis_out = sNO_basis_out;
            end
            function sO_coords_out = complete_sO_tvf_invariants(sO_coords_in_,lam_xi_sO_,Hmat_tvf_,lvs_xi_,lvs_v_,Jltns_xi_,Jl_xi_svd_)
                %% identify differential invariants of the TVF by finding functionally independent conserved quantities at the origin
                Hsvd_tvf_YH_S = Asvd_package( normalize_Henc(Hmat_tvf_*sO_coords_in_.YH) );
                % P x kappa+1, non globally constant function parameters w gradient orthogonal to tvf
                Theta_Eta_tvf_full = sO_coords_in_.YH*Hsvd_tvf_YH_S.W;
                gEta_tvf_svd = Asvd_package(lam_xi_sO_.Jl*Theta_Eta_tvf_full);
                % P x nvar, parameters of canonical independent coordinates with orthogonal gradients at the origin
                gEta_tvf_svd.Theta_Eta_tvf = Theta_Eta_tvf_full*gEta_tvf_svd.V( :, 1:(size(lam_xi_sO_.Jl,1)-1) );

                % 1 x nvar, image of tvf independent coordinates over sO (includes near null space solution)
                Eta_tvf_sO = lam_xi_sO_.lrow_vals*(gEta_tvf_svd.Theta_Eta_tvf);
                % nobs x nvar, image of tvf independent coordinates over S (includes near null space solution)
                Eta_tvf_S = lvs_xi_'*(gEta_tvf_svd.Theta_Eta_tvf);

                YJl = Jl_xi_svd_.V(:,1:(Jl_xi_svd_.r));
                Hsvd_tvf_YJl_S = Asvd_package( normalize_Henc(Hmat_tvf_*YJl) );
                % P x kappa+1, non globally constant function parameters w gradient orthogonal to tvf
                Theta_Eta_tvf_YJl_full = YJl*Hsvd_tvf_YJl_S.W;
                gEta_tvf_YJl_svd = Asvd_package(lam_xi_sO_.Jl*Theta_Eta_tvf_YJl_full);
                % P x nvar, parameters of canonical independent coordinates with orthogonal gradients at the origin
                gEta_tvf_YJl_svd.Theta_Eta_tvf = Theta_Eta_tvf_YJl_full * gEta_tvf_YJl_svd.V( :, 1:(nvar_N1-1) );

                % 1 x nvar, image of tvf independent coordinates over sO (includes near null space solution)
                Eta_tvf_YJl_sO = lam_xi_sO_.lrow_vals*(gEta_tvf_YJl_svd.Theta_Eta_tvf);
                % nobs x nvar, image of tvf independent coordinates over S (includes near null space solution)
                Eta_tvf_YJl_S = lvs_xi_'*(gEta_tvf_YJl_svd.Theta_Eta_tvf);

                %% identify differential invariants of non trivial vfields
                % ntheta x QN, parameters of non trivial vfields over S (transversal to TVF at origin)
                Y_nTVF = sO_coords_in_.sO_nTVF.V(:,1:(nvar_N1-1)) / sO_coords_in_.sO_nTVF.s(1);
                % C x QN, parameters of non trivial vfields over S in full coordinates, scaled,
                theta_nTVF = sO_coords_in_.WYmu * Y_nTVF ...
                                                .* (sO_coords_in_.sO_nTVF.s(1:(nvar_N1-1))' / sO_coords_in_.sO_nTVF.s(1));
                % (QN+1) x Plen x QN, parameters of non-trivial vfields organized as matrices
                theta_nTVF_tns = permute(reshape(theta_nTVF,[],nvar_N1,nvar_N1-1),[2 1 3]);
                % QN x (QN+1) x nobs, non trivial vfield tangent vectors evaluated at every point
                Vspc_nTVF_S = permute( ...
                    reshape(pagemtimes( theta_nTVF_tns , lvs_v_), nvar_N1, nobs, nvar_N1-1 ), ...
                [3 1 2]);
                % QN x Plen x nobs, pages are directional derivatives wrt coordinate vfields at s^(N-1) |_j
                Vspc_nTVF_Jl_S = pagemtimes( Vspc_nTVF_S,Jltns_xi_ );
                Vspc_Jl_svd_S = Asvd_package( ...
                    [ reshape(permute(Vspc_nTVF_Jl_S,[2 1 3]),[],(nvar_N1-1)*nobs)' ; Hmat_tvf_ ]*YJl ...
                );
                % [ reshape(permute(Vspc_nTVF_Jl_S,[2 1 3]),[],(nvar_N1-1)*nobs)' ; Hmat_tvf_ ] ...
                % [ reshape(permute(Vspc_nTVF_Jl_S,[2 1 3]),[],(nvar_N1-1)*nobs)' ; Hmat_tvf_ ]*YJl ...
                % normalize_Henc([ reshape(permute(Vspc_nTVF_Jl_S,[2 1 3]),[],(nvar_N1-1)*nobs)' ; Hmat_tvf_ ]*YJl) ...

                % spans set of non-globally constant functions along transversal vector fields of S
                Y_Vspc_Jl = YJl*Vspc_Jl_svd_S.V(:,1:Vspc_Jl_svd_S.r);
                Hsvd_tvf_Y_VJl_S = Asvd_package( Hmat_tvf_*Y_Vspc_Jl );
                % P x kappa+1, non globally constant function parameters w gradient orthogonal to tvf
                Theta_Eta_tvf_Y_VJl_full = Y_Vspc_Jl*Hsvd_tvf_Y_VJl_S.W;
                gEta_tvf_Y_VJl_svd = Asvd_package(lam_xi_sO_.Jl*Theta_Eta_tvf_Y_VJl_full);
                % P x nvar, parameters of canonical independent coordinates with orthogonal gradients at the origin
                gEta_tvf_Y_VJl_svd.Theta_Eta_tvf = Theta_Eta_tvf_Y_VJl_full * gEta_tvf_Y_VJl_svd.V( :, 1:(nvar_N1-1) );

                % 1 x nvar, image of tvf independent coordinates over sO (includes near null space solution)
                Eta_tvf_Y_VJl_sO = lam_xi_sO_.lrow_vals*(gEta_tvf_Y_VJl_svd.Theta_Eta_tvf);
                % nobs x nvar, image of tvf independent coordinates over S (includes near null space solution)
                Eta_tvf_Y_VJl_S = lvs_xi_'*(gEta_tvf_Y_VJl_svd.Theta_Eta_tvf);

                %% assignments
                sO_coords_out = sO_coords_in_;
                sO_coords_out.Hsvd_tvf_YH_S = Hsvd_tvf_YH_S;

                sO_coords_out.gEta_tvf_svd = gEta_tvf_svd;
                sO_coords_out.Eta_tvf_sO = Eta_tvf_sO;
                sO_coords_out.Eta_tvf_S = Eta_tvf_S;

                sO_coords_out.gEta_tvf_YJl_svd = gEta_tvf_YJl_svd;
                sO_coords_out.Eta_tvf_YJl_sO = Eta_tvf_YJl_sO;
                sO_coords_out.Eta_tvf_YJl_S = Eta_tvf_YJl_S;

                sO_coords_out.Vspc_Jl_svd_S = Vspc_Jl_svd_S;

                sO_coords_out.gEta_tvf_Y_VJl_svd = gEta_tvf_Y_VJl_svd;
                sO_coords_out.Eta_tvf_Y_VJl_sO = Eta_tvf_Y_VJl_sO;
                sO_coords_out.Eta_tvf_Y_VJl_S = Eta_tvf_Y_VJl_S;
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
                [Rsvd_N1_net,DprN_svd] = deal(Rsvd_N1);
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

            %% prepare function space for representation of G kernal
            Pmat_GN1_full = fspace_RN1.Pmat;
            Plen_GN1 = fspace_RN1.Plen;
            ntheta_GN1 = fspace_RN1.ntheta; % ntheta x 1
            LamN_tns_GN1 = LamN_tns_RN1; % ndim x ntheta x nobs

            %% generate G matrices, kernel vector fields satisfy infinitesimal criterion
            J_tau_u_RN1 = nan(ndep_N1,ndim_N1,2,nobs);
            JF_N1 = nan(ndep,ndim,nobs);
            Gtns_N1 = nan(ndep,ntheta_GN1,nobs);
            tic0 = tic;
            for iobs = 1:nobs
                J_tau_u_RN1(:,:,:,iobs) = sols(iobs).lamRN1.J_tau_uN( vth_RN1_net(:,iobs) );
                % holds due to first order ratio condition and DprN enforcement
                JF_N1(:,:,iobs) = ...
                    [ J_tau_u_RN1((end-ndep+1):end,1:nvar_N1,1,iobs) , -eye(ndep) ];
                LamN_i = [ LamN_tns_GN1(1:nvar_N1,:,iobs) ; LamN_tns_GN1((end-ndep+1):end,:,iobs) ];
                % induced inf criterion
                Gtns_N1(:,:,iobs) =  JF_N1(:,:,iobs) * LamN_i;
            end
            toc1 = toc(tic0);
fprintf('(LDsol::model_solspace) encoded G, %dx%dx%d, in %.2f seconds.\n', ...
            size(Gtns_N1,1), size(Gtns_N1,2), size(Gtns_N1,3),  ...
            toc1);

            %% assemble net G matrix, obeying prolongation, as well as optional tvf commutativity constraints
            tic0 = tic;
            [inds_P_GN1_full,inds_P_GN1_net,inds_P_GN1_com] = deal(1:Plen_GN1);
            Gmat_N1 = (reshape(permute(Gtns_N1,[2 1 3]),ntheta_GN1,ndep*nobs))';
            Gsvd_N1 = Asvd_package(normalize_Genc(Gmat_N1));
            Gmat_N1_com = [ Gmat_N1 ; [ Hmat_RN1 , zeros(nobs,ntheta_GN1-Plen_RN1) ] ];
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
                % most restrictive class of vector fields, commute w tvf, possibly underdetermined
                Gmat_N1_com_full = [ Gmat_N1_com ; DprN_mat ];
                Gsvd_N1_com_full = Asvd_package(normalize_Genc(Gmat_N1_com_full));
                nsvd = nsvd + 1;
                [mrow_GN1_com,ntheta_GN1_com] = size(Gmat_N1_com_full);
                % if Gmat_N1_com_full is underdetermined, pare down column space until overdetermined
                if ( mrow_GN1_com<=ntheta_GN1_com )
                    Gtns_T_com = reshape(Gmat_N1_com_full',Plen_GN1,nvar_N1,mrow_GN1_com);
                    [inds_P_GN1_com,ntheta_GN1_com] = pare_colspc(Pmat_GN1_full,mrow_GN1_com);
                    Gmat_N1_com = reshape(Gtns_T_com(inds_P_GN1_com,:,:),ntheta_GN1_com,mrow_GN1_com)';
                    Gsvd_N1_com = Asvd_package( ...
                        normalize_Genc(Gmat_N1_com) ...
                    );
                    nsvd = nsvd+1;
                else
                    Gmat_N1_com = Gmat_N1_com_full;
                    Gsvd_N1_com = Gsvd_N1_com_full;
                end

            else % N=1 => G matrix rows >= R matrix rows, nothing to do
                [Gmat_N1_net,Gmat_N1_net_full] = deal(Gmat_N1);
                [Gsvd_N1_net,Gsvd_N1_net_full] = deal(Gsvd_N1);
                [mrow_GN1_net,ntheta_GN1_net] = size(Gmat_N1_net);

                Gmat_N1_com_full = Gmat_N1_com;
                [Gsvd_N1_com,Gsvd_N1_com_full] = deal( Asvd_package(normalize_Genc(Gmat_N1_com)) );
                nsvd = nsvd+1;
                [mrow_GN1_com,ntheta_GN1_com] = size(Gmat_N1_com);
            end
            toc1 = toc(tic0);
fprintf('(LDsol::model_solspace) Decomposed %d G+DprN matrices in %.2f seconds: %dx%d (r=%d,k=%d) -> %dx%d (r=%d,k=%d,o0=%d,oN=%d,o=%d) \n', ...
            nsvd, toc1, ...
            size(Gmat_N1,1), size(Gmat_N1,2), Gsvd_N1.r, Gsvd_N1.dim - Gsvd_N1.r, ...
            mrow_GN1_net, Gsvd_N1_net.dim, Gsvd_N1_net.r, Gsvd_N1_net.dim - Gsvd_N1_net.r, ...
            max(sum(fspace_RN1.Pmat(1:nvar,inds_P_GN1_net),1)), max(sum(fspace_RN1.Pmat((nvar+1):end,inds_P_GN1_net),1)) , ...
            max(sum(fspace_RN1.Pmat(:,inds_P_GN1_net),1)) );

            Plen_GN1_net = ntheta_GN1_net / nvar_N1;
            % LamN0_tns_GN1_net = permute(reshape(LamN_T_ttns_RN1(inds_P_GN1_net,:,1:nvar_N1,:),ntheta_GN1_net,nvar_N1,nobs),[2 1 3]);

            Plen_GN1_com = ntheta_GN1_com / nvar_N1;
            % LamN0_tns_GN1_com = permute(reshape(LamN_T_ttns_RN1(inds_P_GN1_com,:,1:nvar_N1,:),ntheta_GN1_com,nvar_N1,nobs),[2 1 3]);

            % permute(pagemtimes(Gsvd_N1_com.W',reshape(LamN1_T_ttns(inds_P_GN1_net,:,:,:),ntheta_GN1_net,ndim,nobs)),[2 1 3]);

            %{
                SVDs of Gnet and Gcom reveal kernal vfields of S. The latter are guaranteed to commute with the TVF.

                Using the Gnet kernal vfields, we generate global canonical coordinates.
                Using Gcom kernal vfields, we generate local coordinates suitable for flow transformations, which perturb observed
                integral curves into those passing through arbitrary initial conditions in a neighborhood of the observations
            %}
            % refine the trivial vector field model by intersecting R matrix kernal with Gnet matrix kernal
            % RGsvd_N1_net = Asvd_package([ Rsvd_N1_net.D , Gsvd_N1_net_full.D ]');
            % RGsvd_N1_net = Asvd_package([ normalize_Renc(Rmat_N1_net) ; normalize_Genc(Gmat_N1) ]);

            %% choose an arbitrary origin for the generation of an intrinsic coordinate system
            [icrv_sO,i_sO_0,i_sO_1] = deal(1,1,2); % mid point between first and second observed solutions on curve 1, w.l.o.g.
            jt_O = LDsol.compute_trivial_Hermite_jet( ...
                [Smat(:,i_sO_0) ; tau_uN_RN1_net((end-ndep+1):end,2,i_sO_0)], ...
                [Smat(:,i_sO_1) ; tau_uN_RN1_net((end-ndep+1):end,2,i_sO_1)], ...
                ndep ...
            );
            s_O0 = [ jt_O.xh ; reshape( jt_O.Amat(1:kor,:)', ndep*(kor), 1 ) ]; % extract fitted base space origin, s_O0
            lamN1_sO = adlam( fspace_N1, s_O0 );
            [f_O0,dxf_O0,vth_sO,lamRN1_sO] = comp_f_s0(s_O0,fspace_RN1,Rsvd_N1_net.W,inds_P_RN1); % pass s_O0 to tvf model
            s_O = [ s_O0 ; f_O0((end-ndep+1):end) ]; % set the jet space origin as the graph of tvf on s_O0
            sNp1_O = [ s_O ; dxf_O0((end-ndep+1):end) ];
            t_O = [ 1 ; sNp1_O((nvar+1):end) ]; % tvf tangent vector in the N'th jet space at the origin

            %% use WGcom to validate flow transformation technique

            tic0 = tic;
            [Gcom_sO_coords,Gcom_s0O_basis,Gcom_sNO_basis] = compute_sO_coords( ...
                Gsvd_N1_com,inds_P_GN1_com,lvs_RN1(inds_P_GN1_com,:),lamRN1_sO,t_O,Jltns_N1,LamN_T_ttns_RN1,Gmat_N1_com ...
            );
            toc1 = toc(tic0);
            fprintf('(LDsol::model_solspace) Generated %d canonical coordinates in  %.2f seconds\n', ...
            nvar_N1, toc1 ...
            );
            Gcom_sO_coords = complete_sO_tvf_invariants(Gcom_sO_coords,lamN1_sO,Hmat_N1,lvs_N1,lvs_RN1(inds_P_GN1_com,:),Jltns_N1,Jl_N1_svd);

            tic0 = tic;
            [Gnet_sO_coords,Gnet_s0O_basis,Gnet_sNO_basis] = compute_sO_coords( ...
                Gsvd_N1_net,inds_P_GN1_net,lvs_RN1(inds_P_GN1_net,:),lamRN1_sO,t_O,Jltns_N1,LamN_T_ttns_RN1,Gmat_N1_net ...
            );
            toc1 = toc(tic0);
            fprintf('(LDsol::model_solspace) Generated %d canonical coordinates in  %.2f seconds\n', ...
            nvar_N1, toc1 ...
            );
            Gnet_sO_coords = complete_sO_tvf_invariants(Gnet_sO_coords,lamN1_sO,Hmat_N1,lvs_N1,lvs_RN1(inds_P_GN1_net,:),Jltns_N1,Jl_N1_svd);
            GN1_sO_basis = Gnet_sO_coords;
            GN1_s0O_basis = Gnet_s0O_basis;
            GN1_sNO_basis = Gnet_sNO_basis;

            function flow_out = verify_flow_transformation(Gsvd_,coords_,iPv_)
                % theta_v1 = Gsvd_.W * coords_.V(:,1:nvar_N1) * lsqminnorm( coords_.U(:,1:nvar_N1) * diag(Gcom_sO_coords.s(1:nvar_N1)) , t_O);
                % theta_v1 = theta_v1/norm(t_O);
                % theta_v2 = Gsvd_N1_com.W * Gcom_sO_coords.sO_nTVF.V(:,1);

                flow_out = coords_;
                flow_out.theta_v1 = theta_v1;
            end


            GN1_sO_basis.Xi_sO = Gnet_sO_coords.theta_xicoords' * lamN1_sO.lrow_vals(:) ;
            GN1_sO_basis.Xi_S =  Gnet_sO_coords.theta_xicoords' * lvs_N1 ;
            % GN1_sO_basis.Xi_sO = Gcom_sO_coords.theta_xicoords' * lamN1_sO.lrow_vals(:) ;
            % GN1_sO_basis.Xi_S =  Gcom_sO_coords.theta_xicoords' * lvs_N1 ;

            dXi_S_sO = GN1_sO_basis.Xi_S - GN1_sO_basis.Xi_sO;
            [sortmags_dXi_S_sO,isrtmags_dXi_S_sO] = sort(sqrt(sum(dXi_S_sO.*dXi_S_sO,1)));
            SXi_sO_wgtsrt = sortmags_dXi_S_sO(1)./(sortmags_dXi_S_sO) .* dXi_S_sO;

            SXi_sO_svd = Asvd_package( SXi_sO_wgtsrt' );

            dXi_S_sO_cell = cell([ncrv,1]);
            minmags_dXi_S_sO_crv = nan(ncrv,1);
            for icrv = 1:ncrv
                dXi_S_sO_i = dXi_S_sO(:,ipts_crv(1,icrv):ipts_crv(2,icrv));
                dXi_S_sO_cell{icrv} = dXi_S_sO_i;
                minmags_dXi_S_sO_crv(icrv) = min(sum(dXi_S_sO_i.*dXi_S_sO_i,1));
            end
            [sortmags_dXi_S_sO_crv,isrtmags_dXi_S_sO_crv] = sort(sqrt(minmags_dXi_S_sO_crv));

            tau_S_mat = [ ones(1, nobs) ; uNp1_tvf_mat ];
            i_s_1 = isrtmags_dXi_S_sO(1);
            s_1 = Smat(:,i_s_1);
            t_1 = tau_S_mat(:,i_s_1);
            lam_s_1 = sols(i_s_1).lamRN1;
            GN1_s1_basis = compute_sO_basis(Gsvd_N1_net.W, inds_P_GN1_net, lvs_RN1(inds_P_GN1_net,:), lam_s_1, t_1);

            %% bonus computations

            %% refine G matrix kernel basis
            % nvar x Plen x ntheta, pages are theta row vectors in base space coords, act on lambda column vectors
            WG_tns = permute(reshape(Gsvd_N1_net.W, [Plen_GN1_net nvar_N1 Gsvd_N1_net.dim]),[2 1 3]);
            % nvar x ntheta x nobs, pages are column vectors spanning base space tangent space at each point
            Lam0_WG_tns = permute( pagemtimes(WG_tns,lvs_RN1(inds_P_GN1_net,:)),  [1 3 2] );
            % ndep x ntheta x nobs, pages are column vectors spanning N'th jet space tangent space SECTION at each point
            LamdNxu_WG_tns = permute(pagemtimes( WG_tns((end-ndep+1):end,:,:),dxl_RN1(inds_P_GN1_net,:) ), [1 3 2]) ...
                - pagemtimes( lNx_RN1(:,inds_P_GN1_net,:), reshape(WG_tns(1,:,:),Plen_GN1_net,Gsvd_N1_net.dim) );
            % ndim x ntheta x nobs, pages are column vectors spanning jet space tangent space at each point
            Lam_WG_tns = cat(1,Lam0_WG_tns,LamdNxu_WG_tns);
            % can be viewed as an SVD of a finite subset of the tangent bundle section
            LamWG_svd = Asvd_package(reshape(permute(Lam_WG_tns,[2 1 3]),[ntheta_GN1_net,ndim*nobs])'); % no renormalization
            %% find dominant WG columns parameterizing vfields which correlate with tvf
            % nobs x ntheta, matrix of tangent vector inner product values wrt tvf over observations
            tauT_LamWG_mat = reshape( ...
                pagemtimes( permute(reshape([ones(1,nobs) ; uNp1_tvf_mat],ndim,1,nobs),[2 1 3]), Lam_WG_tns  ), ...
                Gsvd_N1_net.dim, nobs ...
            )';
            % principle components have large component in the direction of tvf
            tauT_LamWG_svd = Asvd_package(tauT_LamWG_mat);

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

            % for i = 1:ndep_N1
            %     fspc.print_vshort_polynomial_theta_z(gEta_svd_tvf_N1.Theta_Eta_tvf(:,i),Pmat_N1,['e' num2str(i)],['\n'])
            % end

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
            mod_out.lamN1_sO = lamN1_sO;

            mod_out.icrv_sO = icrv_sO;
            mod_out.s_O = s_O;
            mod_out.sNp1_O = sNp1_O;
            mod_out.t_O = t_O;

            mod_out.Gsvd_N1 = Gsvd_N1;
            mod_out.Gsvd_N1_com = Gsvd_N1_com;
            mod_out.Gsvd_N1_net = Gsvd_N1_net;

            mod_out.LamWG_svd = LamWG_svd;
            mod_out.tauT_LamWG_svd = tauT_LamWG_svd;

            mod_out.GN1_sO_basis = GN1_sO_basis;
            mod_out.GN1_s0O_basis = GN1_s0O_basis;
            mod_out.GN1_sNO_basis = GN1_sNO_basis;
            mod_out.Gcom_sNO_basis = Gcom_sNO_basis;
            mod_out.Gcom_s0O_basis = Gcom_s0O_basis;
            mod_out.Gnet_sNO_basis = Gnet_sNO_basis;
            mod_out.Gnet_s0O_basis = Gnet_s0O_basis;

            mod_out.Gnet_sO_coords = Gnet_sO_coords;
            mod_out.Gcom_sO_coords = Gcom_sO_coords;
            % mod_out.GnetXi_svds = GnetXi_svds;
            % mod_out.Hxi_svds = Hxi_svds;
            % mod_out.Theta_v_coords = Theta_v_coords;
            % mod_out.theta_xi_coords = theta_xi_coords;
            % mod_out.err_xi_coords = err_xi_coords;

            mod_out.dXi_S_sO = dXi_S_sO;
            mod_out.sortmags_dXi_S_sO = sortmags_dXi_S_sO;
            mod_out.isrtmags_dXi_S_sO = isrtmags_dXi_S_sO;
            mod_out.SXi_sO_svd = SXi_sO_svd;

            mod_out.dXi_S_sO_cell = dXi_S_sO_cell;
            mod_out.sortmags_dXi_S_sO_crv = sortmags_dXi_S_sO_crv;
            mod_out.isrtmags_dXi_S_sO_crv = isrtmags_dXi_S_sO_crv;

            mod_out.s_1 = s_1;
            mod_out.t_1 = t_1;
            mod_out.GN1_s1_basis = GN1_s1_basis;
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
