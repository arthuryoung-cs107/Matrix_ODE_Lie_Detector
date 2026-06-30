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
            % normalize_Renc = @(r_) r_./sqrt(sum(r_.*r_,2)); % unit length Rmat rows
            normalize_Renc = @(r_) r_; % unnormalized Rmat rows
            % normalize_DprN = @(d_) d_./sqrt( sum(d_.*d_,1) ); % unit length DprN cols (rows after transposition)
            normalize_DprN = @(d_) d_; % unnormalized DprN cols (rows after transposition)

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
                % bor_0 = 3;
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
                % ord_i = max(reshape(sum(Pmat_N1((nvar+1):end,:),1),[],1));
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
                [fspace_0,fspace_RN1] = deal(fspace_N1);
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
            lvs_RN1 = nan(Plen_RN1,nobs);
            LamN_T_tns_RN1 = nan(ntheta_RN1,ndim_N1,nobs);
            Renc_tns_RN1 = nan(ntheta_RN1,ndep_N1,nobs);
            DprN_T_ttns = zeros(ntheta_RN1,ndep,kor-1,nobs);
            tic0 = tic;
            for iobs = 1:nobs
                sols(iobs) = LDsol(Smat(:,iobs));

                sols(iobs).lamN1 = adlam( fspace_N1, xumat_N1(:,iobs) );
                Hmat_N1(iobs,:) = [ 1.0 , dxumat(:,iobs)' ]*(sols(iobs).lamN1.Jl);

                sols(iobs).lam0 = adlam( fspace_0, xumat(:,iobs), Smat((nvar+1):end,iobs) );
                Hmat_0(iobs,:) = [ 1.0 , dxumat(1:ndep,iobs)' ]*(sols(iobs).lam0.Jl);

                sols(iobs).lamRN1 = adlam( fspace_RN1, xumat_N1(:,iobs), Smat((nvar+1):end,iobs) );
                Hmat_RN1(iobs,:) = [ 1.0 , dxumat(:,iobs)' ]*(sols(iobs).lamRN1.Jl);

                lRN1_i = sols(iobs).lamRN1;
                lv_RN1_i = lRN1_i.lrow_vals;
                l_imm_i = fspace_RN1.imm_l(lv_RN1_i);

                lvs_RN1(:,iobs) = lv_RN1_i;
                LamN_T_tns_RN1(:,1:nvar_N1,iobs) = ...
                    [ [ lv_RN1_i(:) ; zeros(ntheta_RN1-Plen_RN1,1) ] , [ zeros(Plen_RN1,ndep_N1) ; l_imm_i' ] ];
                Renc_tns_RN1(:,:,iobs) = (normalize_Renc([ -lRN1_i.dxu(:,1)*lv_RN1_i , l_imm_i ]))';

                LamN_T_tns_RN1(:,(nvar_N1+1):(nvar_N1+ndep_N1),iobs) = ...
                    [ -lRN1_i.lkx(:,:,1) , fspace_RN1.imm_l(lRN1_i.dkxl(1,:)) ]';
                iiLam0 = nvar + (1:ndep);
                iiLam1 = iiLam0-nvar+nvar_N1;
                for k = 2:kor
                    DprN_T_ttns(:,:,k-1,iobs) = normalize_DprN(LamN_T_tns_RN1(:,iiLam1,iobs) - LamN_T_tns_RN1(:,iiLam0,iobs));
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

            LamN_tns_RN1 = permute(LamN_T_tns_RN1,[2 1 3]); % --> ndim (N1), ntheta_N1, nobs

            H_N1_svd = Asvd_package(Hmat_N1);
            H_0_svd = Asvd_package(Hmat_0);
            H_RN1_svd = Asvd_package(Hmat_RN1);

            function vth_out = comp_vartheta_RN1(W_,lvs_)
                vx_W = (W_(1:(size(W_,1)/nvar_N1),:))' * lvs_;
                % Deltas_W = sum(vx_W.*vx_W,1);
                % vth_out = W_ * (vx_W./Deltas_W);
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
            function [f_s0_out, vth_out, lam_out] = comp_f_s0(s0_,fspc_,W_,iP_)
                nvar_ = length(s0_(:));
                Plen_b = size(fspc_.Pmat(:,iP_),2);
                lam_out = adlam( fspc_, s0_ );
                vx_W = ( W_(1:Plen_b,:) )' * reshape( lam_out.lrow_vals(iP_) , [], 1 );
                vth_out = W_ * (vx_W./sum(vx_W.*vx_W,1));
                f_s0_out = reshape(vth_out((Plen_b+1):end),Plen_b,[])' * (lam_out.lrow_vals(iP_))';

                lam_out = lam_out.prolong_jet_space( reshape(f_s0_out,[],1) );
            end

            tic0 = tic;
            inds_P_RN1 = 1:Plen_RN1;
            Rmat_N1 = (reshape(Renc_tns_RN1,ntheta_RN1,ndep_N1*nobs))'; % nobs*Qhat by ntheta, where Qhat = kor*Q
            Rsvd_N1 = Asvd_package(Rmat_N1); % global R1 svd
            mrow_R_net = size(Rmat_N1,1);
            %{
                Computed the svd of an "overdetermined" R matrix, which we actually expect to be rank deficient.
                Letting KR denote an orthonormal basis for the kernal of R, which has dimension kappa,
                for each k = 1, ... , N, there exists a subspace of ker(R) which has th \in \R^C also satisfying
                    D^k th = vdkxu - pr1 (vdkm1xu) = l . th_dkxu - dx (l . th_dkm1xu) + dkxu dx (l.th_x) = 0,
                that is, first order prolongation adherence for each k = 1, ..., N, automatically satisfied when N=1.
                Intersect these nullspaces, then find lowest order rank deficient R matrix for just sufficient
                complexity trivial vector field model.
            %}
            nsvd = 1;
            if (kor>1)
                DprN_mat = reshape(DprN_T_ttns,ntheta_RN1,ndep*(kor-1)*nobs)';
                DprN_svd = Asvd_package( DprN_mat );
                Rmat_N1_DprN = [ Rsvd_N1.D , DprN_svd.D ]';
                Rsvd_N1_net = Asvd_package( Rmat_N1_DprN );
                mrow_R_net = size(Rmat_N1_DprN,1); % redefine column space dimension to that of concatenated, transposed row space
                Rtns_T_net = reshape(Rmat_N1_DprN',Plen_RN1,nvar_N1,mrow_R_net);
                nsvd = 3;
                %% further pare down jet space model to smallest rank deficient RN1 space
                if (Rsvd_N1_net.r < Rsvd_N1_net.dim)
                    inds_P_RN1_full = inds_P_RN1;

                    %% pare down order of coordinate functions for derivatives
                    ord_i = max(reshape(Pmat_RN1_full((nvar+1):end,:),[],1));
                    while ( (Rsvd_N1_net.r < Rsvd_N1_net.dim)&&(ord_i>1) )
                        ord_i = ord_i - 1;
                        inds_P_i = inds_P_RN1_full(sum(Pmat_RN1_full((nvar+1):end,:),1) <= ord_i);
                        Plen_i = length(inds_P_i(:));
                        nsvd = nsvd + 1;
                        Rsvd_i = Asvd_package(reshape(Rtns_T_net(inds_P_i,:,:),nvar_N1*Plen_i,mrow_R_net)');
                        if (Rsvd_i.r < Rsvd_i.dim) % accept sub svd if rank deficient
                            inds_P_RN1 = inds_P_i;
                            Rsvd_N1_net = Rsvd_i;
                        else % reject, maintain previous svd, break
                            break;
                        end
                    end

                    %% pare down order of coordinate functions over base space
                    ord_i = max(reshape(Pmat_RN1_full(1:nvar,:),[],1));
                    while ( (Rsvd_N1_net.r < Rsvd_N1_net.dim)&&(ord_i>1) )
                        inds_P_i = inds_P_RN1_full(sum(Pmat_RN1_full(1:nvar,:),1) <= ord_i);
                        Plen_i = length(inds_P_i(:));
                        nsvd = nsvd + 1;
                        Rsvd_i = Asvd_package(reshape(Rtns_T_net(inds_P_i,:,:),nvar_N1*Plen_i,mrow_R_net)');
                        if (Rsvd_i.r < Rsvd_i.dim) % accept sub svd if rank deficient
                            inds_P_RN1 = inds_P_i;
                            Rsvd_N1_net = Rsvd_i;
                            ord_i = ord_i - 1;
                        else % reject, maintain previous svd, break
                            break;
                        end
                    end
                end
            else
                Rmat_N1_DprN = Rmat_N1; % automatically equivalent
                [Rsvd_N1_net,DprN_svd] = deal(Rsvd_N1);
                Rtns_T_net = reshape(Rmat_N1_DprN',Plen_RN1,nvar_N1,mrow_R_net);

                %% further pare down jet space model to smallest rank deficient R1 space
                if (Rsvd_N1_net.r < Rsvd_N1_net.dim)
                    inds_P_RN1_full = inds_P_RN1;
                    %% pare down order of coordinate functions over base space
                    ord_i = max(Pmat_RN1_full(:));
                    while ( (Rsvd_N1_net.r < Rsvd_N1_net.dim)&&(ord_i>1) )
                        inds_P_i = inds_P_RN1_full(sum(Pmat_RN1_full,1) <= ord_i);
                        Plen_i = length(inds_P_i(:));
                        nsvd = nsvd + 1;
                        Rsvd_i = Asvd_package(reshape(Rtns_T_net(inds_P_i,:,:),nvar_N1*Plen_i,mrow_R_net)');
                        if (Rsvd_i.r < Rsvd_i.dim) % Accept sub svd if rank deficient. Decrement ord_i
                            inds_P_RN1 = inds_P_i;
                            Rsvd_N1_net = Rsvd_i;
                            ord_i = ord_i - 1;
                        else % reject, maintain previous svd, break
                            break;
                        end
                    end
                end
            end
            toc1 = toc(tic0);
fprintf('(LDsol::model_solspace) Decomposed %d R + DprN matrices in %.2f seconds: %dx%d (r=%d,k=%d) -> %dx%d (r=%d,k=%d,o0=%d,oN=%d) \n', ...
            nsvd, toc1, ...
            size(Rmat_N1,1), size(Rmat_N1,2), Rsvd_N1.r, Rsvd_N1.dim-Rsvd_N1.r, ...
            mrow_R_net, Rsvd_N1_net.dim, Rsvd_N1_net.r, Rsvd_N1_net.dim-Rsvd_N1_net.r, ...
            max(sum(fspace_RN1.Pmat(1:nvar,inds_P_RN1),1)), max(sum(fspace_RN1.Pmat((nvar+1):end,inds_P_RN1),1)) );

            vth_RN1_net = comp_vartheta_RN1( Rsvd_N1_net.W , lvs_RN1(inds_P_RN1,:) ); % ntheta by nobs
            tau_uN_RN1_net = comp_tau_uN_RN1( vth_RN1_net , inds_P_RN1 ); % ndep by 2 by nobs

            %% choose an arbitrary origin for the generation of an intrinsic coordinate system
            [i_sO_0,i_sO_1] = deal(1,2); % mid point between first and second observed solutions on curve 1, w.l.o.g.
            % Smat(:,i_sO_0), ...
            % Smat(:,i_sO_1), ...
            jt_O = LDsol.compute_trivial_Hermite_jet( ...
                [Smat(:,i_sO_0) ; tau_uN_RN1_net((end-ndep+1):end,2,i_sO_0)], ...
                [Smat(:,i_sO_1) ; tau_uN_RN1_net((end-ndep+1):end,2,i_sO_1)], ...
                ndep ...
            );
            s_O0 = [ jt_O.xh ; reshape( jt_O.Amat(1:kor,:)', ndep*(kor), 1 ) ]; % extract fitted base space origin, s_O0
            [f_O0,vth_sO,lamRN1_sO] = comp_f_s0( s_O0 , fspace_RN1, Rsvd_N1_net.W, inds_P_RN1 ); % pass s_O0 to tvf model
            s_O = [ s_O0 ; f_O0((end-ndep+1):end) ]; % set the jet space origin as the graph of tvf on s_O0

            %% if R was found to be rank deficient over a subspace, lift vth back into full space
            if ( size(vth_RN1_net,1) ~= ntheta_RN1 )
                vth_RN1_tns = reshape(vth_RN1_net,[],nvar_N1,nobs);
                vth_RN1_full_tns = zeros(Plen_RN1,nvar_N1,nobs);
                vth_RN1_full_tns(inds_P_RN1,:,:) = vth_RN1_tns;
                vth_RN1_net = reshape(vth_RN1_full_tns,ntheta_RN1,nobs);

                vth_sO_mat = zeros(Plen_RN1,nvar_N1);
                vth_sO_mat(inds_P_RN1,:) = reshape(vth_sO,[],nvar_N1);
                vth_sO = reshape(vth_sO_mat,ntheta_RN1,1);
            end

            %% prepare function space for representation of G kernal
            Pmat_GN1_full = fspace_RN1.Pmat;
            Plen_GN1 = fspace_RN1.Plen;
            ntheta_GN1 = fspace_RN1.ntheta;

            J_tau_u_RN1 = nan(ndep_N1,ndim_N1,2,nobs);
            JF_N1 = nan(ndep,ndim,nobs);
            Gtns_N1 = nan(ndep,ntheta_GN1,nobs);
            Tmat_N1 = nan(nobs,ntheta_GN1);
            tic0 = tic;
            for iobs = 1:nobs
                J_tau_u_RN1(:,:,:,iobs) = sols(iobs).lamRN1.J_tau_uN( vth_RN1_net(:,iobs) );
                % holds due to first order ratio condition and DprN enforcement
                JF_N1(:,:,iobs) = ...
                    [ J_tau_u_RN1((end-ndep+1):end,1:nvar_N1,1,iobs) , -eye(ndep) ];
                LamN_i = [ LamN_tns_RN1(1:nvar_N1,:,iobs) ; LamN_tns_RN1((end-ndep+1):end,:,iobs) ];
                % induced inf criterion
                Gtns_N1(:,:,iobs) =  JF_N1(:,:,iobs) * LamN_i;
                % trivial vector field orthogonality condition
                Tmat_N1(iobs,:) = ...
                    [ 1 , tau_uN_RN1_net(:,1,iobs)', tau_uN_RN1_net((end-ndep+1):end,2,iobs)' ]*LamN_i;
            end
            toc1 = toc(tic0);
            Tsvd_N1 = Asvd_package(Tmat_N1);
            toc2 = toc(tic0);
fprintf('(LDsol::model_solspace) encoded G+T, %dx%dx%d + %dx%d, in %.2f seconds, decomposed T in %.2f (r=%d,k=%d).\n', ...
            size(Gtns_N1,1), size(Gtns_N1,2), size(Gtns_N1,3),  ...
            size(Tmat_N1,1), size(Tmat_N1,2),  ...
            toc1, ...
            toc2-toc1, ...
            Tsvd_N1.r, Tsvd_N1.dim - Tsvd_N1.r );

            tic0 = tic;
            inds_P_GN1 = 1:Plen_GN1;
            Gmat_N1 = (reshape(permute(Gtns_N1,[2 1 3]),ntheta_GN1,ndep*nobs))';
            Gsvd_N1 = Asvd_package(Gmat_N1);
            mrow_GN1_net = size(Gmat_N1,1);
            nsvd = 1;
            if (kor>1)
                Gmat_N1_net = [ Gsvd_N1.D , DprN_svd.D ]';
                Gsvd_N1_net = Asvd_package( Gmat_N1_net );
                mrow_GN1_net = size(Gmat_N1_net,1); % redefine column space dimension to that of concatenated, transposed row space

                Gtns_T_net = reshape(Gmat_N1_net',Plen_GN1,nvar_N1,mrow_GN1_net);
                nsvd = 2;
                if (Gsvd_N1_net.r < Gsvd_N1_net.dim)
                    Pmat_GN1_full = Pmat_RN1_full;
                    inds_P_GN1_full = inds_P_GN1;

                    %% pare down order of coordinate functions for derivatives
                    ord_i = max(reshape(Pmat_GN1_full((nvar+1):end,:),[],1));
                    if ( ord_i == max(reshape( sum(Pmat_GN1_full((nvar+1):end,:), 1) ,[],1)) )
                        ord_i = ord_i -1; % decrement if already at restricted order
                    end
                    while ( (Gsvd_N1_net.r < Gsvd_N1_net.dim)&&(ord_i>1) )
                        inds_P_i = inds_P_GN1_full(sum(Pmat_GN1_full((nvar+1):end,:),1) <= ord_i);
                        Plen_i = length(inds_P_i(:));
                        nsvd = nsvd + 1;
                        Gsvd_i = Asvd_package(reshape(Gtns_T_net(inds_P_i,:,:),nvar_N1*Plen_i,mrow_GN1_net)');
                        if (Gsvd_i.r < Gsvd_i.dim) % accept sub svd if rank deficient
                            inds_P_GN1 = inds_P_i;
                            Gsvd_N1_net = Gsvd_i;
                            ord_i = ord_i - 1;
                        else % reject, maintain previous svd, break
                            break;
                        end
                    end
                    %% pare down order of coordinate functions over base space
                    ord_i = max(reshape(Pmat_GN1_full(1:nvar,:),[],1));
                    while ( (Gsvd_N1_net.r < Gsvd_N1_net.dim)&&(ord_i>1) )
                        inds_P_i = inds_P_GN1_full(sum(Pmat_GN1_full(1:nvar,:),1) <= ord_i);
                        Plen_i = length(inds_P_i(:));
                        nsvd = nsvd + 1;
                        Gsvd_i = Asvd_package(reshape(Gtns_T_net(inds_P_i,:,:),nvar_N1*Plen_i,mrow_GN1_net)');
                        if (Rsvd_i.r < Rsvd_i.dim) % accept sub svd if rank deficient
                            inds_P_GN1 = inds_P_i;
                            Gsvd_N1_net = Gsvd_i;
                            ord_i = ord_i - 1;
                        else % reject, maintain previous svd, break
                            break;
                        end
                    end
                end
            else
                Gsvd_N1_net = Gsvd_N1;
            end
            toc1 = toc(tic0);
fprintf('(LDsol::model_solspace) Decomposed %d G+DprN matrices in %.2f seconds: %dx%d (r=%d,k=%d) -> %dx%d (r=%d,k=%d,o0=%d,oN=%d) \n', ...
            nsvd, toc1, ...
            size(Gmat_N1,1), size(Gmat_N1,2), Gsvd_N1.r, Gsvd_N1.dim - Gsvd_N1.r, ...
            mrow_GN1_net, Gsvd_N1_net.dim, Gsvd_N1_net.r, Gsvd_N1_net.dim - Gsvd_N1_net.r, ...
            max(sum(fspace_RN1.Pmat(1:nvar,inds_P_GN1),1)), max(sum(fspace_RN1.Pmat((nvar+1):end,inds_P_GN1),1)) );

            tic0 = tic;
            inds_P_G_N1_TV = 1:Plen_GN1;
            Gmat_N1_TV = [ Gsvd_N1.D , Tsvd_N1.D ]'; % kernal contains non-trivial vector fields, unrestricted by prolongation
            Gsvd_N1_TV = Asvd_package( Gmat_N1_TV );
            mrow_G_TV_net = size(Gmat_N1_TV,1);
            nsvd = 1;
            if (kor>1)
                Gmat_N1_TV_net = [ Gsvd_N1_TV.D , DprN_svd.D ]';  % kernal contains non-trivial vector fields obeying prolongation
                Gsvd_N1_TV_net = Asvd_package( Gmat_N1_TV_net );
                mrow_G_TV_net = size(Gmat_N1_TV_net,1); % redefine column space dimension to that of concatenated, transposed row space

                Gtns_T_net = reshape(Gmat_N1_TV_net',Plen_GN1,nvar_N1,mrow_G_TV_net);
                nsvd = 2;
                if (Gsvd_N1_TV_net.r < Gsvd_N1_TV_net.dim)
                    Pmat_GN1_full = Pmat_RN1_full;
                    inds_P_GN1_full = inds_P_G_N1_TV;

                    %% pare down order of coordinate functions for derivatives
                    ord_i = max(reshape(Pmat_GN1_full((nvar+1):end,:),[],1));
                    if ( ord_i == max(reshape( sum(Pmat_GN1_full((nvar+1):end,:), 1) ,[],1)) )
                        ord_i = ord_i -1; % decrement if already at restricted order
                    end
                    while ( (Gsvd_N1_TV_net.r < Gsvd_N1_TV_net.dim)&&(ord_i>1) )
                        inds_P_i = inds_P_GN1_full(sum(Pmat_GN1_full((nvar+1):end,:),1) <= ord_i);
                        Plen_i = length(inds_P_i(:));
                        nsvd = nsvd + 1;
                        Gsvd_i = Asvd_package(reshape(Gtns_T_net(inds_P_i,:,:),nvar_N1*Plen_i,mrow_G_TV_net)');
                        if (Gsvd_i.r < Gsvd_i.dim) % accept sub svd if rank deficient
                            inds_P_G_N1_TV = inds_P_i;
                            Gsvd_N1_TV_net = Gsvd_i;
                            ord_i = ord_i - 1;
                        else % reject, maintain previous svd, break
                            break;
                        end
                    end
                    %% pare down order of coordinate functions over base space
                    ord_i = max(reshape(Pmat_GN1_full(1:nvar,:),[],1));
                    while ( (Gsvd_N1_TV_net.r < Gsvd_N1_TV_net.dim)&&(ord_i>1) )
                        inds_P_i = inds_P_GN1_full(sum(Pmat_GN1_full(1:nvar,:),1) <= ord_i);
                        Plen_i = length(inds_P_i(:));
                        nsvd = nsvd + 1;
                        Gsvd_i = Asvd_package(reshape(Gtns_T_net(inds_P_i,:,:),nvar_N1*Plen_i,mrow_G_TV_net)');
                        if (Gsvd_i.r < Gsvd_i.dim) % accept sub svd if rank deficient
                            inds_P_G_N1_TV = inds_P_i;
                            Gsvd_N1_TV_net = Gsvd_i;
                            ord_i = ord_i - 1;
                        else % reject, maintain previous svd, break
                            break;
                        end
                    end
                end
            else
                Gsvd_N1_TV_net = Gsvd_N1_TV;
            end
            toc1 = toc(tic0);
fprintf('(LDsol::model_solspace) Decomposed %d G+DprN+T matrices in %.2f seconds: %dx%d (r=%d,k=%d) -> %dx%d (r=%d,k=%d,o0=%d,oN=%d) \n', ...
            nsvd, toc1, ...
            size(Gmat_N1_TV,1), size(Gmat_N1_TV,2), Gsvd_N1_TV.r, Gsvd_N1_TV.dim-Gsvd_N1_TV.r, ...
            mrow_G_TV_net, Gsvd_N1_TV_net.dim, Gsvd_N1_TV_net.r, Gsvd_N1_TV_net.dim-Gsvd_N1_TV_net.r, ...
            max(sum(fspace_RN1.Pmat(1:nvar,inds_P_G_N1_TV),1)), max(sum(fspace_RN1.Pmat((nvar+1):end,inds_P_G_N1_TV),1)) );
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
            err_stats( 'tau_u_RN1_net', comp_err_msr(tau_u_RN1_net) );

            % keyboard
            mod_out = struct( ...
                'Smat', Smat, ...
                'Smat_N1', Smat_N1, ...
                'Hmat_N1', Hmat_N1, ...
                'Hmat_0', Hmat_0, ...
                'Hmat_RN1', Hmat_RN1, ...
                'LamN_tns_RN1', LamN_tns_RN1, ...
                'inds_P_RN1', inds_P_RN1, ...
                'Rmat_N1', Rmat_N1, ...
                'inds_P_GN1', inds_P_GN1, ...
                'Gmat_N1', Gmat_N1, ...
                'inds_P_G_N1_TV', inds_P_G_N1_TV, ...
                'Gmat_N1_TV', Gmat_N1_TV, ...
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

            mod_out.H_N1_svd = H_N1_svd;
            mod_out.H_0_svd = H_0_svd;
            mod_out.H_RN1_svd = H_RN1_svd;
            mod_out.DprN_svd = DprN_svd;

            mod_out.Rsvd_N1 = Rsvd_N1;
            mod_out.Rsvd_N1_net = Rsvd_N1_net;

            mod_out.jt_O = jt_O;

            mod_out.f_O0 = f_O0;
            mod_out.vth_sO = vth_sO;
            mod_out.lamRN1_sO = lamRN1_sO;

            mod_out.s_O = s_O;

            mod_out.Tsvd_N1 = Tsvd_N1;

            mod_out.Gsvd_N1 = Gsvd_N1;
            mod_out.Gsvd_N1_net = Gsvd_N1_net;

            mod_out.Gsvd_N1_TV = Gsvd_N1_TV;
            mod_out.Gsvd_N1_TV_net = Gsvd_N1_TV_net;
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
                Amat(:,i) = linsolve(Vmat, Umat(:,i) );
            end

            jt_out = struct( ...
                'xh', xh, ...
                'fJ', fJ, ...
                'Vmat', Vmat, ...
                'Umat', Umat, ...
                'Amat', Amat ...
            );
            % uhmat = ( Amat(1:kor,:) )';
            % s_out = [ xh ; reshape( uhmat ,ndep_*kor,1) ];
        end
    end
    methods
        function obj = LDsol(xu_)
            obj.xu = xu_(:);
        end
    end
end
