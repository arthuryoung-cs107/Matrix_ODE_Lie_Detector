classdef apv < fspc

    properties

        obs;

        Omap_a;
        Omap_b;

        % S;
        % s;
        P_mat;

        L_mat;
        l_vec;

        pz_coeff;
        pz_Lz;

        Ck_pL_nz;
        L;
        l;
        gxu_l;
        veps_gradxu;
        dx_l;

        vxu;
        txu;

        Pbor;
        Pbor_sub;

    end

    methods

        % automatically prolonging vector
        function obj = apv(obs_)

            obj@fspc(obs_);

            obj.obs = obs_;

            obj.Omap_a = ones(1,obj.ndep+1);
            obj.Omap_b = zeros(1,obj.ndep+1);

        end
        function Plen_out = Plen(obj)
            Plen_out = size(obj.P_mat,2);
        end
        function [Pmat_out,Plen_out] = Pmat_Plen(obj)
            [Pmat_out,Plen_out] = deal(obj.P_mat,size(obj.P_mat,2));
        end
        function [Pbor_thin, Pbor_square] = get_polynomial_bor_projection(obj,bor_,sub_flag_)
            if (nargin==3)
                sub_flag = sub_flag_;
            else
                sub_flag = false;
            end
            if (bor_ < obj.bor)
                if (sub_flag)
                    iibor_vec = double( sum(obj.P_mat,1) <= bor_ );
                else
                    iibor_vec = prod( double(obj.P_mat <= bor_) ,1 );
                end
                iibor_long = (iibor_vec(:)).*( ones(1,size(obj.P_mat,1)) );
                ntheta_hot = sum( iibor_long );
                Pbor_square = eye( prod(size(obj.P_mat)) ) .* iibor_long(:);
                Pbor_thin = Pbor_square(:,logical(iibor_long));
            else
                [Pbor_thin,Pbor_square] = deal(eye( prod(size(obj.P_mat)) ));
            end
        end
        function [txu_out,vartheta_out] = compute_txu(obj,lmat_,V_,s_)
            if (nargin==4)
                W = V_.*fspc.s_spi(s_);
            else
                W = V_;
            end
            [eor,ndep,ndim,nvar] = jspc.unpack_jetspace_dims(obj);
            Plen = size(W,1)/nvar;
            [dlm1,dlm2] = size(lmat_);
            if (dlm1 == Plen)
                nobs = dlm2;
                lmat = lmat_';
            else
                nobs = dlm1;
                lmat = lmat_;
            end
            vartheta_evl = @(vx_)  ( W*vx_ )/sum(vx_.*vx_);
            Wx = W(1:Plen,:);
            if (nobs==1)
                vartheta_out = vartheta_evl( (lmat*Wx)' );
                txu_out = (lmat*reshape(vartheta_out,Plen,nvar))';
            else
                txu_out = nan(nvar,nobs);
                vartheta_out = nan(Plen*nvar,nobs);
                for i = 1:nobs
                    vartheta_out(:,i) = vartheta_evl( (lmat(i,:)*Wx)' );
                    txu_out(:,i) = (lmat(i,:)*reshape(vartheta_out(:,i),Plen,nvar))';
                end
            end
        end
        function obj_out = init_polynomial_fspace(obj, bor_)

            bor = bor_;
            [P_len order_mat dim0_lens] = count_set_P_len(bor,obj.ndep);
            P_mat = order_mat';

            obj_out = obj;

            obj_out.Omap_a = ones(1,obj.ndep+1);
            obj_out.Omap_b = zeros(1,obj.ndep+1);

            obj_out.bor = bor;

            obj_out.P_mat = P_mat;
            obj_out.L_mat = @(s_,P_) ( s_ ).^( P_ );
            obj_out.l_vec = @(s_,P_) prod(( s_ ).^( P_ ),1);

            obj_out.pz_coeff = @(iz_,pr_) double( pr_.Ck_pL(iz_,:)-pr_.pz_count(iz_) );
            obj_out.pz_Lz = @(iz_,pr_) pr_.Lvals( iz_, pr_.Ck_pL(iz_,:)-pr_.pz_count(iz_)+1 ) ;

            obj_out.Ck_pL_nz = @(obj_) obj_.P_mat;
            obj_out.L = @(obj_,s_) obj_.L_polynomial(s_);
            obj_out.l = @(obj_,s_) obj_.l_polynomial(s_);
            obj_out.gxu_l = @(obj_,xu_) obj_.gxu_l_polynomial(xu_);
            obj_out.dx_l = @(obj_,xu_,dxu_) obj_.dx_l_polynomial(xu_,dxu_);

            obj_out.veps_gradxu = @(obj_,xu_,th_) obj_.veps_gradxu_polynomial(xu_,th_);

            obj_out.vxu = @(obj_,xu_,theta_) obj_.vxu_polynomial(xu_,theta_);
            obj_out.txu = @(obj_,xu_,W_) obj_.compute_txu(obj_.l_polynomial(xu_),W_);

            obj_out.Pbor = @(obj_,bor_) obj_.get_polynomial_bor_projection(bor_);
            obj_out.Pbor_sub = @(obj_,bor_) obj_.get_polynomial_bor_projection(bor_,true);

        end
        function vxu_out = vxu_polynomial(obj,xu_,theta_)
            [xu,nobs,Pmat,Plen] = obj.unpack_xu_fspace(xu_);
            [eor,ndep,ndim,nvar] = jspc.unpack_jetspace_dims(obj);

            vxu_evl = @(s_,t_) ( prod(s_ .^ Pmat,1)*reshape(t_,Plen,[]) )';

            if (nobs==1)
                vxu_out = vxu_evl(xu,theta_);
            else
                ntheta = Plen*nvar;
                theta_cols = reshape(theta_,ntheta,[]);
                ncol_theta = size(theta_cols,2);
                vxu_out = nan(nvar,nobs);
                if (ncol_theta==1)
                    for i = 1:nobs
                        vxu_out(:,i) = vxu_evl(xu(:,i),theta_cols);
                    end
                else
                    for i = 1:nobs
                        vxu_out(:,i) = vxu_evl(xu(:,i),theta_cols(:,i));
                    end
                end
            end
        end
        function L_out = L_polynomial(obj,xu_)
            [xu,nobs,Pmat,Plen] = obj.unpack_xu_fspace(xu_);
            L_mat_evl = @(s_) s_ .^ Pmat;
            if (nobs==1)
                L_out = L_mat_evl(xu);
            else
                L_out = nan(obj.ndep+1,Plen,nobs);
                for i = 1:nobs
                    L_out(:,:,i) = L_mat_evl(xu(:,i));
                end
            end
        end
        function l_out = l_polynomial(obj,xu_)
            [xu,nobs,Pmat,Plen] = obj.unpack_xu_fspace(xu_);
            l_vec_evl = @(s_) prod(s_ .^ Pmat,1);
            if (nobs==1)
                l_out = l_vec_evl(xu);
            else
                l_out = nan(nobs,Plen);
                for i = 1:nobs
                    l_out(i,:) = l_vec_evl(xu(:,i));
                end
            end
        end
        function [gxu_out,L_out] = gxu_l_polynomial(obj,xu_)
            [xu,nobs,Pmat,Plen] = obj.unpack_xu_fspace(xu_);
            [eor,ndep,ndim,nvar] = jspc.unpack_jetspace_dims(obj);
            ip_tns = zeros(nvar,Plen,nvar);
            for i = 1:nvar
                ip_tns(i,:,i) = 1;
            end
            ip_vec = logical(ip_tns(:));
            ones_tns = ones(nvar,Plen,nvar);
            len_ones_tns = prod(size(ones_tns));
            len_hot_ip_vec = sum(ip_tns(:));
            function [ gxu_l_out,L_mat_out ] = comp_gxu_l(s_)

                L_mat_out = s_ .^ Pmat;
                L_vec_i = reshape(ones_tns.*L_mat_out, len_ones_tns,1);
                % substitute partial derivatives
                L_vec_i(ip_vec) = reshape((Pmat.*( s_ .^ (Pmat-1) ))', len_hot_ip_vec,1);
                gxu_l_out = reshape(prod(reshape(L_vec_i,nvar,Plen,nvar),1),Plen,nvar)';

            end
            if (nobs==1)
                [gxu_out,L_out] = comp_gxu_l(xu);
            else
                [gxu_out,L_out] = deal(nan( nvar,Plen,nobs ));
                for i = 1:nobs
                    [gxu_out(:,:,i),L_out(:,:,i)] = comp_gxu_l( xu(:,i) );
                end
            end
        end

function [veps_out,gxu_out,L_out] = veps_gradxu_polynomial(obj,xu_,theta_eps_)
    [xu,nobs,Pmat,Plen] = obj.unpack_xu_fspace(xu_);
    [eor,ndep,ndim,nvar] = jspc.unpack_jetspace_dims(obj);
    theta_eps = reshape(theta_eps_,Plen);
    ip_tns = zeros(nvar,Plen,nvar);
    for i = 1:nvar
        ip_tns(i,:,i) = 1;
    end
    ip_vec = logical(ip_tns(:));
    ones_tns = ones(nvar,Plen,nvar);
    len_ones_tns = prod(size(ones_tns));
    len_hot_ip_vec = sum(ip_tns(:));
    function [ veps_l_out,gxu_l_out,L_mat_out ] = comp_veps_l(s_)

        L_mat_out = s_ .^ Pmat;
        L_vec_i = reshape(ones_tns.*L_mat_out, len_ones_tns,1);
        % substitute partial derivatives
        L_vec_i(ip_vec) = reshape((Pmat.*( s_ .^ (Pmat-1) ))', len_hot_ip_vec,1);
        gxu_l_out = reshape(prod(reshape(L_vec_i,nvar,Plen,nvar),1),Plen,nvar)';
        veps_l_out = gxu_l_out*theta_eps;


    end
    if (nobs==1)
        [veps_out,gxu_out,L_out] = comp_veps_l(xu);
    else
        [gxu_out,L_out] = deal(nan( nvar,Plen,nobs ));
        veps_l_out = nan(nvar,nobs);
        for i = 1:nobs
            [veps_l_out(:,:,i),gxu_out(:,:,i),L_out(:,:,i)] = comp_veps_l( xu(:,i) );
        end
    end
end
        function [dx_out,gxu_out,L_out,l_out] = dx_l_polynomial(obj,xu_,dxu_)
            [xu,nobs,Pmat,Plen] = obj.unpack_xu_fspace(xu_);
            [eor,ndep,ndim,nvar] = jspc.unpack_jetspace_dims(obj);
            dxu = reshape(dxu_,ndep,nobs);
            ip_tns = zeros(nvar,Plen,nvar);
            for i = 1:nvar
                ip_tns(i,:,i) = 1;
            end
            ip_vec = logical(ip_tns(:));
            ones_tns = ones(nvar,Plen,nvar);
            len_ones_tns = prod(size(ones_tns));
            len_hot_ip_vec = sum(ip_tns(:));
            function [ dx_l_out, gxu_l_out, L_mat_out ] = comp_dx_l(s_,dxu_)
                L_mat_out = s_ .^ Pmat;
                L_vec_i = reshape(ones_tns.*L_mat_out, len_ones_tns,1);
                % substitute partial derivatives
                L_vec_i(ip_vec) = reshape((Pmat.*( s_ .^ (Pmat-1) ))', len_hot_ip_vec,1);
                gxu_l_out = reshape(prod(reshape(L_vec_i,nvar,Plen,nvar),1),Plen,nvar)';
                % finish by computing action of trivial vfield on gradient
                dx_l_out = sum( [ 1.0 ; dxu_ ].*gxu_l_out , 1 );
            end
            if (nobs==1)
                [dx_out,gxu_out,L_out] = comp_dx_l( xu,dxu );
            else
                [gxu_out,L_out] = deal(nan( nvar,Plen,nobs ));
                dx_out = nan(nobs,Plen);
                for i = 1:nobs
                    [dx_out(i,:),gxu_out(:,:,i),L_out(:,:,i)] = ...
                     comp_dx_l( xu(:,i),dxu(:,i) );
                end
            end
            if (nargout == 4)
                l_out = (reshape(prod(L_out,1),Plen,nobs))';
            else
                l_out = [];
            end
        end
        function [Vx_,Vu_] = split_Vxu_mat(obj,V_)
            [Vx_,Vu_] = deal( V_(1:obj.Plen,:) , V_((obj.Plen+1):end,:) );
        end
        function [Pmat_out,Plen_out] = unpack_P_fspace(obj)
            [Pmat_out,Plen_out] = deal(obj.P_mat,size(obj.P_mat,2));
        end
        function [xu_out,nobs_out,Pmat_out,Plen_out] = unpack_xu_fspace(obj,xu_)
            [Pmat_out,Plen_out] = deal(obj.P_mat,size(obj.P_mat,2));
            xu_out = reshape(xu_,obj.ndep+1,[]);
            nobs_out = size(xu_out,2);
        end
    end
    methods (Static)

        function [tvf_out,Renc_out,Renc_specs] = model_trivial_vfield(obj_,obs_,S_)

            [Renc_out,Renc_specs] = apv.encode_R( obj_,obs_,S_ );

            [eor,ndep] = deal(obj_.eor,obj_.ndep);
            ndim = 1 + ndep*(eor+1);
            nvar = ndep+1;
            [Pmat,Plen] = obj_.Pmat_Plen;

            ntheta = size(Renc_out{1},2);

            nobs = Renc_specs.nobs;
            nset = Renc_specs.nset;
            kor = Renc_specs.kor_obs;
            ndim_obs = Renc_specs.ndim_obs;
            Smat = jspc.Scell_2_Smat(S_,ndim_obs);

            L_S = Renc_specs.L_S;
            gxu_l_S = Renc_specs.gxu_l_S;
            dx_l_S = Renc_specs.dx_l_S;
            xu = Renc_specs.xu;
            un = Renc_specs.un;
            dxun = Renc_specs.dxun;
            d1xu = Renc_specs.d1xu;
            Lam_x_dxu = Renc_specs.Lam_x_dxu;
            Lam_u_dxu = Renc_specs.Lam_u_dxu;

            % l_S = (reshape(prod(L_S,1),Plen,[]))';
            l_S = Renc_specs.lambda_evl(L_S);
            Lambda_k_S = Renc_specs.prk_Lambda(l_S,Lam_x_dxu,Lam_u_dxu);
            Lambda_1_S = Lambda_k_S(1:(nvar+ndep),:,:);
            Lambda_dxu_S = Lambda_1_S( (end-(ndep-1)):end,:,: );

            %{
                RRk :
            %}
            rsV_Rk_cell = cell(3,kor,nset);
            RRk_full_ttns = nan(ntheta,ntheta,kor,nset);
            idel = 0;
            for i = 1:nset
                npts_i = size(S_{i},2);
                inds_i = (1+idel):(npts_i+idel);
                for k = 1:kor

                    % concatenated matrix of Rk constraints, all dependent variables together
                    Rkmat_i = reshape( permute(Renc_out{k,i},[2 1 3]) ,ntheta,[])';

                    [rsV_Rk_cell{1,k,i},rsV_Rk_cell{2,k,i},rsV_Rk_cell{3,k,i}] = ...
                        fspc.rsV_unpack(Rkmat_i);

                    % save row space representation of this Rk matrix for full null space problem
                    RRk_full_ttns(:,:,k,i) = rsV_Rk_cell{3,k,i}.*( rsV_Rk_cell{2,k,i}' );

                end
                idel = idel + npts_i;
            end

            % Evaluate aggregate R1 svd: always exists, always relevant.
            R1_svd = fspc.compute_svd_package(RRk_full_ttns(:,:,1,:));
            tvf_out = R1_svd;
            %{
                Evaluate spectral tvf model on observational data
                by computing image of local parameter map and passing
                to base space Lambda matrix.
            %}
            [tvf_out.txu_S,tvf_out.vartheta_S] = obj_.txu(obj_,xu,R1_svd.W);
            tvf_out.vartheta_svd = fspc.compute_svd_package((tvf_out.vartheta_S)');

            tvf_out.l_S = l_S;
            tvf_out.Lambda_k_S = Lambda_k_S;
            tvf_out.Lambda_1_S = Lambda_1_S;
            tvf_out.Lambda_dxu_S = Lambda_dxu_S;

            %% extra credit

            [tx_S,tu_S] = deal(tvf_out.txu_S(1,:), tvf_out.txu_S(2:end,:));
            tu_mat_S = reshape(tu_S,ndep,nobs);

            % [ txu_K_S , vartheta_K_S ] = obj_.txu(obj_,xu,K_R1); % pass back to txu
            vartheta_tns_S = reshape(tvf_out.vartheta_S,Plen,nvar,nobs);

            % T0_full_mat = nan(ntheta,nobs);
            % [T1_full_mat,H1_full_mat] = deal(nan(ntheta,nobs));
            % T1_full_mat = nan(ntheta,nobs);
            Jtu_S = nan(ndep,nvar,nobs);
            tdxu_S = nan(ndep,nobs);
            for i = 1:nobs

                vtxu_i = vartheta_tns_S(:,:,i);
                vtx_i = vtxu_i(:,1);
                vtu_i = vtxu_i(:,2:end);
                g_i = gxu_l_S(:,:,i)';

                % Jacobian of tvf over base space
                Jtu_S(:,:,i) = (vtu_i' - ( tu_mat_S(:,i) .* vtx_i' ))*g_i;

                % tvf estimate of the second derivative of u wrt x
                tdxu_S(:,i) = (tvf_out.Lambda_dxu_S(:,:,i))*(vtxu_i(:));

            end

            tvf_out.Jtu_S = Jtu_S;
            tvf_out.tdxu_S = tdxu_S;

            % tdxu_true_S = obj_.obs.dxf( xu(1,:) , [xu(2:end,:) ; d1xu ] );
            % norm( tdxu_true_S(:) - tdxu_S(:) )
        end
        function [mod_out,tvf_out,Renc_out,Renc_specs] = model_Fode_observations(obj_,obs_,S_)

            [tvf_out,Renc_out,Renc_specs] = apv.model_trivial_vfield( obj_,obs_,S_ )

            [eor,ndep] = deal(obj_.eor,obj_.ndep);
            ndim = 1 + ndep*(eor+1);
            nvar = ndep+1;
            [Pmat,Plen] = obj_.Pmat_Plen;

            ntheta = size(Renc_out{1},2);
            nobs = Renc_specs.nobs;
            nset = Renc_specs.nset;
            kor = Renc_specs.kor_obs;
            ndim_obs = Renc_specs.ndim_obs;
            Smat = jspc.Scell_2_Smat(S_,ndim_obs);

            L_S = Renc_specs.L_S;
            gxu_l_S = Renc_specs.gxu_l_S;
            dx_l_S = Renc_specs.dx_l_S;
            xu = Renc_specs.xu;
            un = Renc_specs.un;
            dxun = Renc_specs.dxun;
            d1xu = Renc_specs.d1xu;
            Lam_x_dxu = Renc_specs.Lam_x_dxu;
            Lam_u_dxu = Renc_specs.Lam_u_dxu;

            txu_S = tvf_out.txu_S;
            [tx_S,tu_S] = deal(txu_S(1,:), txu_S(2:end,:));
            tu_mat_S = reshape(tu_S,ndep,nobs);
            vartheta_tns_S = reshape(tvf_out.vartheta_S,Plen,nvar,nobs);

            l_S = tvf_out.l_S;
            Lambda_1_S = tvf_out.Lambda_1_S;
            Lambda_0_S = tvf_out.Lambda_1_S(1:nvar,:,:);

            Jtu_S = tvf_out.Jtu_S;
            tdxu_S = tvf_out.tdxu_S;

            [T1_full_mat,H1_full_mat] = deal(nan(ntheta,nobs));
            G1_full_ttns = zeros(ntheta,ndep,nobs);
            % G1_true_ttns = zeros(ntheta,ndep,nobs);
            for i = 1:nobs
                vtxu_i = vartheta_tns_S(:,:,i);
                vtx_i = vtxu_i(:,1);
                vtu_i = vtxu_i(:,2:end);
                g_i = gxu_l_S(:,:,i)';
                l_i = l_S(i,:)';

                Lam1_x_i = reshape(Lam_x_dxu(:,:,1,i),ndep,Plen)';
                Lam1_u_i = Lam_u_dxu(:,1,i);

                Lam_dxu_T_i = G1_full_ttns(:,:,i);
                Lam0_xu_T_i = zeros(ntheta,nvar);

                % T0_full_mat(1:Plen,i) = l_i; % tx*vx = 1*vx
                Lam0_xu_T_i(1:Plen,1) = l_i;
                idel_P = Plen;
                for idep = 1:ndep
                    Lam_dxu_T_i(1:Plen,idep) = Lam1_x_i(:,idep);

                    inds_ii = (1+idel_P):(Plen+idel_P);

                    Lam_dxu_T_i(inds_ii,idep) = Lam1_u_i;
                    Lam0_xu_T_i(inds_ii,idep+1) = l_i;

                    % T0_full_mat(inds_ii,i) = d1xu(:,idep)*l_i; % tuq*vuq = dxuq*li

                    idel_P = idel_P + Plen;
                end

                G1_full_ttns(:,:,i) = ...
                ( Jtu_S(:,:,i) * Lam0_xu_T_i' - Lam_dxu_T_i')';
                % ( (vtu_i' - (tu_mat_S(:,i).*vtx_i'))*g_i *Lam0_xu_T_i' - Lam_dxu_T_i')';

                T1_full_mat(:,i) = ( [ txu_S(:,i) ; tdxu_S(:,i) ]' * Lambda_1_S(:,:,i) )';

                % oracle verification
                % G1_true_ttns(:,:,i) = ...
                % (obj_.obs.gradf(Smat(1,i),Smat(2:nvar,i))'*Lam0_xu_T_i' - Lam_dxu_T_i')';
            end

            G1_svd = fspc.compute_svd_package(G1_full_ttns)
                % G1_true_svd = fspc.compute_svd_package(G1_true_ttns)

            bor_cap =  obj_.bor_cap(obj_);
            bor_range = 1:(bor_cap-1);
            Pbor_cell = cell(bor_cap-1,1);

            [Pbor_cell{1},~] = obj_.Pbor(obj_,1);
            G1Pbor_svd(1) = fspc.compute_svd_package(G1_svd.mat*Pbor_cell{1});
            KGnet = Pbor_cell{1} * G1Pbor_svd(1).K(G1Pbor_svd(1));
            dim_KGnet = size(KGnet,2);
            for ibor = 2:(bor_cap-1)
                % [Pbor_cell{ibor},~] = obj_.Pbor(obj_,ibor);
                [Pbor_cell{ibor},~] = obj_.Pbor_sub(obj_,ibor);
                G1Pbor_svd(ibor) = fspc.compute_svd_package(G1_svd.mat*Pbor_cell{ibor});
                if ( G1Pbor_svd(ibor).k(G1Pbor_svd(ibor)) > dim_KGnet )
                    [KGnet,dim_KGnet,s_KGnet] = ...
                    fspc.append_to_basis( KGnet, G1Pbor_svd(ibor).K(G1Pbor_svd(ibor)) );
                end
            end
            if ( G1_svd.k(G1_svd) > dim_KGnet )
                [KGnet,dim_KGnet,s_KGnet] = ...
                fspc.append_to_basis( KGnet, G1_svd.K(G1_svd), size(G1_svd.mat,1) );
            end

            % compute differential invariants of extracted vector field basis
            for i = 1:dim_KGnet
                fspace_i = obj_;
                ki = KGnet(:,i);

                %{
                    seek eta = lambda . theta^eta over base space such that
                    v (eta) = 0. These are the conserved quantities along this flow
                %}
                H_svd(i) = fspc.compute_svd_package( ...
                reshape( ...
                pagemtimes(permute(pagemtimes( Lambda_0_S , ki(:) ),[2 1 3]), Lambda_0_S), ...
                ntheta,nobs ...
                ));
            end

            % R1 G1 kernel intersection is in the kernel of their scaled concatenation
            R1G1_svd = fspc.compute_svd_package( [ (tvf_out.Y)' ; (G1_svd.Y)' ] );

            WGmWR = G1_svd.W - (R1G1_svd.W)*(R1G1_svd.W')*(G1_svd.W);
            KGmKR = G1_svd.K(G1_svd) - ...
             (R1G1_svd.K(R1G1_svd))*(R1G1_svd.K(R1G1_svd)')*(G1_svd.K(G1_svd));

            WGmWR_svd = fspc.compute_svd_package( WGmWR,true );
            KGmKR_svd = fspc.compute_svd_package( KGmKR,true );

            T1_svd = fspc.compute_svd_package(T1_full_mat);

            mod_out = G1_svd;
            mod_out.bor_range = bor_range;
            mod_out.Pbor_cell = Pbor_cell;
            mod_out.G1Pbor_svd = G1Pbor_svd;
            mod_out.KGnet = KGnet;
            mod_out.H_svd = H_svd;
            mod_out.s_KGnet = s_KGnet;
            mod_out.R1G1_svd = R1G1_svd;
            mod_out.T1_svd = T1_svd;
            mod_out.KGmKR_svd = KGmKR_svd;
            mod_out.WGmWR_svd = WGmWR_svd;

            % mod_out = 0
            fprintf('\n(model_Fode_observations end)\n');
            % pause
        end

        function [Renc_out,enc_specs] = encode_R(obj_,obs_,S_)

            [eor,ndep,ndim,nvar] = jspc.unpack_jetspace_dims(obj_);
            [Smat,nobs,nset,kor_obs,ndim_obs] = jspc.unpack_Scell(S_,ndep);

            Pmat = obj_.P_mat;
            Plen = size(Pmat,2);

            xumat = Smat(1:nvar,:);
            untns = reshape( Smat(2:end,:), ndep,eor+1,nobs );
            dxuntns = reshape(untns(:,2:end,:),ndep,eor,nobs);
            d1xumat = reshape(dxuntns(:,1,:),ndep,nobs);
            %{
                First ratio condition needs no prolongation
                First prolongation is free after computing first total derivative
            %}
            [dx_l_S,gxu_l_S,L_S,l_S] = obj_.dx_l(obj_,xumat,d1xumat);
            ntheta = nvar*Plen;

            %% stage gradient evaluation by initializing injection into Lambda column space
            ncol_Lambda_u = ntheta-Plen;
            i_imm = zeros(ncol_Lambda_u,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            i_imm_vec = logical( i_imm(:) );
            len_Lambda_u = prod(size(i_imm));
            l_imm_init = zeros(len_Lambda_u,1);
            % ones_imm = ones(ndep,Plen);
            ones_imm = ones(ndep,1);
            function l_imm = immerse_lambda(l_)
                l_imm = l_imm_init;
                l_imm(i_imm_vec) = reshape( (ones_imm.*l_)', len_Lambda_u ,1);
                l_imm = reshape(l_imm,ndep,ncol_Lambda_u);
                % l_ -> [l_ 0 ... 0 ; 0 l_ ... 0 ; ...]
            end

            %{
                The first prolongation block of the Lambda matrix is free as
                soon as one has computed the first total derivative
            %}
            R1_tns = nan(ndep,ntheta,nobs);
            Lam_x_dxu_ttns = nan(ndep,Plen,kor_obs,nobs);
            Lam_u_dxu_tns = nan(Plen,kor_obs,nobs);
            for i = 1:nobs
                % check1 = (-d1xumat(:,i)*l_S(i,:)) % also works
                R1_tns(:,:,i) = [ (-d1xumat(:,i).*l_S(i,:)) , immerse_lambda(l_S(i,:)) ];
                % free evaluations of Lambda 1
                Lam_x_dxu_ttns(:,:,1,i) = -d1xumat(:,i).*dx_l_S(i,:);
                Lam_u_dxu_tns(:,1,i) = dx_l_S(i,:);
            end

            Renc_out = cell(kor_obs,nset);
            R_ttns = nan(ndep,ntheta,kor_obs,nobs);
            idel = 0;
            for i = 1:nset
                nobs_set_i = size(S_{i},2);
                iinds_i = (1+idel):(nobs_set_i+idel);

                [Renc_out{1,i},R_ttns(:,:,1,iinds_i)] = deal(R1_tns( :,:,iinds_i ));

                idel = idel+nobs_set_i;
            end

            Lam_x_ttns = zeros(size(Lam_x_dxu_ttns));
            Lam_u_tns = zeros(size(Lam_u_dxu_tns));
            g_vdxu_ttns = zeros(ndep,ndim,kor_obs,nobs);
            % if (kor_obs>1)

                iLv = ones(1,Plen);
                for i = 2:nvar
                    iLv = iLv.*double( Pmat(i-1,:)==Pmat(i,:) );
                end
                iLv = logical(iLv);

                prn_spc = struct( ...
                    'kor', kor_obs, ...
                    'ndep', ndep, ...
                    'ndim', ndim, ...
                    'nvar',  nvar, ...
                    'iobs', 1, ...
                    'iLv', iLv, ...
                    'P', Pmat, ...
                    'S', Smat, ...
                    'untns_S', untns, ... % ndep,eor+1,nobs
                    'dxuntns_S', dxuntns, ... % ndep,eor,nobs
                    'l_S', l_S, ...
                    'L_S', L_S, ...
                    'gxu_l_S', gxu_l_S, ...
                    'dx_l_S', dx_l_S, ...
                    'Ck_pL_0', obj_.Ck_pL_nz(obj_), ...
                    'inds_dkxu', @(spc_,k_) (1+(spc_.ndep*k_))+(1:spc_.ndep), ...
                    'xun', @(spc_) spc_.S(:,spc_.iobs), ...
                    'dxu', @(spc_) spc_.S( spc_.nvar+(1:spc_.ndep) , spc_.iobs), ...
                    'dxuq', @(spc_,q_) spc_.S( spc_.nvar+q_ , spc_.iobs ), ...
                    'txuk', @(spc_,k_) [ 1 ; spc_.untns_S(:,2:(k_+1),spc_.iobs) ], ...
                    'l', @(spc_) spc_.l_S(spc_.iobs,:), ...
                    'L', @(spc_) spc_.L_S(:,:,spc_.iobs), ...
                    'gxu_l', @(spc_) spc_.gxu_l_S(:,:,spc_.iobs), ...
                    'dx_l', @(spc_) spc_.dx_l_S(spc_.iobs,:), ...
                    'i_pz', @(pr_,iz_) (pr_.Ck_pL(iz_,:)-pr_.pz_count(iz_)) > 0, ...
                    'pz_ncoeff', @(pr_,iz_) obj_.pz_coeff(iz_,pr_), ...
                    'evl_pz_Lz', @(pr_,iz_) obj_.pz_Lz(iz_,pr_), ...
                    'evl_vz_hat', @(pr_) prod(pr_.Ldxu_hat(:))*pr_.cl_hat ...
                );
                % prn_str = cell(kor_obs+1,Plen,nobs);

                function [pro_, vz_hat] = init_pr_pz_ipz(spc_,pri_,iz_,ipz_)
                    pro_ = pri_;
                    % restrict to non-zero indices
                    pro_.iPl = pri_.iPl(:,ipz_);
                    pro_.Ck_pL = pri_.Ck_pL(:,ipz_);
                    pro_.coeffs = pri_.coeffs(ipz_);
                    pro_.L_hat = pri_.L_hat(:,ipz_);

                    pro_.coeffs = pro_.coeffs .* spc_.pz_ncoeff(pro_,iz_); % compute new coefficients

                    pro_.pz_count(iz_) = pro_.pz_count(iz_) + 1; % increment partial derivative count
                    pro_.L_hat(iz_,:) = spc_.evl_pz_Lz(pro_,iz_); % update L_z function value history
                    if (iz_>1)
                        % increase the power of the L polynomial by 1
                        idep = iz_-1;
                        pro_.P_dxu(idep,1) = pro_.P_dxu(idep,1)+1;
                        pro_.Ldxu_hat(idep,1) = pro_.Ldxu_hat(idep,1)*pro_.dxu(idep,1);
                    end
                    % recompute lambda maps, product of L functions
                    pro_.l_hat = prod(pro_.L_hat , 1 );
                    pro_.cl_hat = pro_.coeffs(:) .* pro_.l_hat(:);
                    if (nargout == 2)
                        vz_hat = spc_.evl_vz_hat(pro_);
                    end
                end
                function Lamu = prk_vu(k_,spc_,prkm1_vu_)
                    Lamu = zeros(length(prkm1_vu_.iPl),spc_.kor);
                    % del_Lamu = zeros(size(Lamu_));
                    if (k_<=spc_.kor)
                        for ivar = 1:nvar
                            % indices of lambda maps which have non-zero partial der. in z
                            i_pz = spc_.i_pz(prkm1_vu_,ivar);
                            if (sum(double(i_pz)))
                                prk_dxz = init_pr_pz_ipz(spc_,prkm1_vu_,ivar,i_pz);
                                Lamu(prk_dxz.iPl,k_) = ...
                                    Lamu(prk_dxz.iPl,k_) + prod(prk_dxz.Ldxu_hat(:))*prk_dxz.vz_hat;
                                if (k_ < spc_.kor)
                                    % pass into next prolongation
                                    Lamu(prk_dxz.iPl,:) = prk_vu(k_+1,spc_,prk_dxz);
                                end
                            end
                        end
                        ivar = nvar+1;
                        for kord = 1:(k_-1)
                            for idep = 1:ndep
                                if ( prkm1_vu_.P_dxu(idep,kord) )
                                    prk_dx_dxuq = prkm1_vu_;

                                    prk_dx_dxuq.coeffs = prk_dx_dxuq.coeffs * pro_.P_dxu(idep,kord);

                                    prk_dx_dxuq.P_dxu(idep,kord) = pro_.P_dxu(idep,kord)-1;
                                    prk_dx_dxuq.Ldxu_hat(idep,kord) = ...
                                        prk_dx_dxuq.dxu(idep,kord)^prk_dx_dxuq.P_dxu(idep,kord);

                                    prk_dx_dxuq.P_dxu(idep,kord+1) = pro_.P_dxu(idep,kord+1)+1;
                                    prk_dx_dxuq.Ldxu_hat(idep,kord+1) = ...
                                        prk_dx_dxuq.Ldxu_hat(idep,kord+1)*prk_dx_dxuq.dxu(idep,kord+1);

                                    Lamu(prk_dxz.iPl,k_) = ...
                                        Lamu(prk_dxz.iPl,k_) + spc_.evl_vz_hat(prk_dxz);
                                end
                                ivar = ivar + 1;
                            end
                        end
                    end
                end
                function Lamu = prk_vx(k_,k_src_,spc_,prkm1_vx_)

                    %% To do. Pending results from operator overloaded autodiff

                end
                function [Lamx,Lamu,prn_vxu,g_vdxu] = prn_obs_i(spc_,iobs_)
                    spc_.iobs = iobs_;

                    l_i = spc_.l(spc_);
                    dxl_i = spc_.dx_l(spc_) ;
                    s_i = spc_.xun(spc_);
                    L0 = spc_.L(spc_);
                    %% initialize base space prolongation state.
                    pr0_ = struct( ...
                    'Lam_x_1_val', -spc_.dxu(spc_) .* dxl_i, ...
                    'Lam_u_1_val', dxl_i(:), ...
                    's', s_i, ...
                    'x', s_i(1), ...
                    'u', reshape( s_i(2:end) , ndep, kor_obs+1 ), ...
                    'dxu', reshape( s_i((nvar+1):ndim),ndep,kor_obs ), ...
                    'Lvals', L0(:, spc_.iLv ), ...
                    'iPl', 1:Plen, ...
                    'Ck_pL', spc_.Ck_pL_0, ...
                    'pz_count', zeros(nvar,1), ...
                    'P_dxu', zeros(ndep,spc_.kor), ...
                    'coeffs', ones(1,Plen), ...
                    'L_hat', L0, ...
                    'l_hat', l_i, ...
                    'cl_hat', l_i, ...
                    'Ldxu_hat', ones(ndep,spc_.kor) ...
                    );

                    % allocate full prolongation block for a point
                    Lamu = zeros(Plen,spc_.kor);
                    Lamx = zeros(ndep,Plen,spc_.kor);
                    g_vdxu = zeros(ndep,spc_.ndim,spc_.kor);
                    for ivar = 1:nvar
                        % indices of lambda maps which have non-zero partial der. in z
                        i_pz = spc_.i_pz(pr0_,ivar);
                        if (sum(double(i_pz)))
                            % del_Lamu = compute_pr_pz_ipz(spc_,pr0_,ivar,i_pz);
                            [pr1_gxu(ivar),vz_hat] = init_pr_pz_ipz(spc_,pr0_,ivar,i_pz);

                            Lamu(pr1_gxu(ivar).iPl,1) = ...
                                Lamu(pr1_gxu(ivar).iPl,1) + vz_hat;

                            % pass into second prolongation, if only to compute Jacobian
                            Lamu(pr1_gxu(ivar).iPl,:) = Lamu(pr1_gxu(ivar).iPl,1) ...
                                + prk_vu(2,spc_,pr1_gxu(ivar));

                            %% accumulate vx contributions to k = 1, ..., kor prolongation coords
                            for kord = 1:spc_.kor
                                Lamx(:,pr1_gxu(ivar).iPl,kord) = ...
                                    Lamx(:,pr1_gxu(ivar).iPl,kord) ...
                                    + ( -pr0_.dxu(:,kord) .* (vz_hat(:))' );
                                if (kord < spc_.kor)
                                    Lamx(:,pr1_gxu(ivar).iPl,:) = ...
                                        Lamx(:,pr1_gxu(ivar).iPl,kord) ...
                                        + prk_vx(kord+1,kord,spc_,pr1_gxu(ivar));
                                end
                            end
                        end
                    end
                    prn_vxu = pr1_gxu;

                end

                %% populate output Lambda matrix data
                for iobs = 1:nobs
                    % [Lam_x_ttns(:,:,:,iobs), Lam_u_tns(:,:,iobs), prn_vxu_obs_i ] ...
                    [Lam_x_ttns(:,:,:,iobs), Lam_u_tns(:,:,iobs), ~ ] ...
                        = prn_obs_i(prn_spc,iobs);

                    % prn_vx_obs_i = prn_vxu_obs_i(1);
                    % prn_vu_obs_i = prn_vxu_obs_i(2:end);
                    % Lam_x_ttns(:,:,:,iobs), Lam_u_tns(:,:,iobs);
                    % prn_vx_obs_i.Lam_x_1_val;
                    % prn_vu_obs_i.Lam_u_1_val;
                end

                norm(Lam_x_dxu_ttns(:)-Lam_x_ttns(:))
                norm(Lam_u_dxu_tns(:)-Lam_u_tns(:))
                % pause

            % end

            enc_specs = struct( ...
            'R_ttns', R_ttns, ...
            'R1_tns', R1_tns, ...
            'nobs',nobs, ...
            'nset',nset, ...
            'kor_obs',kor_obs, ...
            'ndim_obs',ndim_obs, ...
            'dx_l_S',dx_l_S, ...
            'gxu_l_S',gxu_l_S, ...
            'L_S', L_S, ...
            'xu', xumat, ...
            'un', untns, ...
            'dxun', dxuntns, ...
            'd1xu', d1xumat, ...
            'lambda_evl', @(L_) (reshape(prod(L_,1),size(L_,2),[]))', ...
            'prk_Lambda', @(l_,Lx_,Lu_) fspc.comp_prk_Lambda(l_,Lx_,Lu_), ...
            'Lam_x_dxu', Lam_x_dxu_ttns, ...
            'Lam_u_dxu', Lam_u_dxu_tns ...
            );

        end

    end

end

function [P_len_out pow_mat_out dim0_lens_out] = count_set_P_len(bor_,ndep_)
    d_ = ndep_ + 1;
    P_len_out = (double(bor_+1))^(double(d_));

    dim0_lens_out = nan(bor_+1,1);
    pow_mat_out = nan(P_len_out,d_);

    k_perm = 0;
    for i = 0:bor_
        [delk pow_mat_out] = set_powmat_full_recursive(pow_mat_out,bor_,d_,1,k_perm);
        dim0_lens_out(i+1) = delk;
        for k = k_perm:(k_perm+delk-1)
            pow_mat_out(k+1,1) = i;
        end
        k_perm = k_perm + delk;
    end
end

function [delk_out pow_mat_out] = set_powmat_full_recursive(mat_, bor_, d_, ilevel_, k_perm_)
    pow_mat_out = mat_;
    if (ilevel_>=d_)
        delk_out = 0;
    else
        delk_level = 0;
        if (ilevel_ == (d_-1))
            for i = 0:bor_
                pow_mat_out(k_perm_+delk_level+1,ilevel_+1) = i;
                delk_level = delk_level + 1;
            end
        else
            for i = 0:bor_
                [delk pow_mat_out] = set_powmat_full_recursive(pow_mat_out,bor_,d_,ilevel_+1,delk_level + k_perm_);
                for k = 0:(delk-1)
                    pow_mat_out(k_perm_+delk_level+1,ilevel_+1) = i;
                    delk_level = delk_level + 1;
                end
            end
        end
        delk_out = delk_level;
    end
end
