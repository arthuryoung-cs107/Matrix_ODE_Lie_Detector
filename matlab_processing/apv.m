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

        L;
        l;
        gxu_l;
        dx_l;

        vxu;
        txu;

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

            obj_out.L = @(obj_,s_) obj_.L_polynomial(s_);
            obj_out.l = @(obj_,s_) obj_.l_polynomial(s_);
            obj_out.gxu_l = @(obj_,xu_) obj_.gxu_l_polynomial(xu_);
            obj_out.dx_l = @(obj_,xu_,dxu_) obj_.dx_l_polynomial(xu_,dxu_);

            obj_out.vxu = @(obj_,xu_,theta_) obj_.vxu_polynomial(xu_,theta_);
            obj_out.txu = @(obj_,xu_,W_) obj_.compute_txu(obj_.l_polynomial(xu_),W_);
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
                [gxu_l_out,L_mat_out] = comp_gxu_l(xu);
            else
                [L_out,gxu_out] = deal(nan( nvar,Plen,nobs ));
                for i = 1:nobs
                    [gxu_out(:,:,i),L_out(:,:,i)] = comp_gxu_l( xu(:,i) );
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
                [L_out,gxu_out] = deal(nan( nvar,Plen,nobs ));
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

        function [Renc_out,enc_specs,coord_specs] = model_trivial_vfield(obj_,obs_,S_)

            [Renc_out,Renc_specs] = apv.encode_R( obj_,obs_,S_ )

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
            l_S = (reshape(prod(L_S,1),Plen,[]))';
            gxu_l_S = Renc_specs.gxu_l_S;
            dx_l_S = Renc_specs.dx_l_S;

            xu = Renc_specs.xu;
            un = Renc_specs.un;
            dxun = Renc_specs.dxun;
            d1xu = Renc_specs.d1xu;

            Lam_x_dxu = Renc_specs.Lam_x_dxu;
            Lam_u_dxu = Renc_specs.Lam_u_dxu;

            %{
                H1 :
                RRk :
            %}

            % T1_full_tns = nan(ntheta,ntheta,nset);
            H1_full_tns = nan(Plen,Plen,nset);
            rsV_H1_cell = cell(3,nset);
            rsV_Rk_cell = cell(3,kor,nset);
            RRk_full_ttns = nan(ntheta,ntheta,kor,nset);
            idel = 0;
            for i = 1:nset

                npts_i = size(S_{i},2);
                inds_i = (1+idel):(npts_i+idel);

                H1_i = dx_l_S(inds_i,:);
                % T1_i = fspc.Lambda_xu_immersion( l_S(i        ) );

                [rsV_H1_cell{1,i},rsV_H1_cell{2,i},rsV_H1_cell{3,i}] = ...
                    fspc.rsV_unpack(H1_i);
                H1_full_tns(:,:,i) = rsV_H1_cell{3,i}.*( rsV_H1_cell{2,i}' );

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

            H1_svd = fspc.compute_svd_package(H1_full_tns)

            % Evaluate aggregate R1 svd: always exists, always relevant.
            R1_svd = fspc.compute_svd_package(RRk_full_ttns(:,:,1,:))
            % [W_R1_mat,rank_R1_full,s_R1_full,V_R1_full,R1_full_mat] = ...
            %     fspc.safely_process_net_svd(RRk_full_ttns(:,:,1,:));

            fspc.print_vshort_polynomial_theta(R1_svd.W(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(R1_svd.W(:,end-1),obj_.P_mat);

            [R1_svd.Wx,R1_svd.Wu] = obj_.split_Vxu_mat(R1_svd.W);
            WR1x_svd = fspc.compute_svd_package(R1_svd.Wx)

            % these are the scaled singular vectors of the Vx submatrix.
            fspc.print_vshort_polynomial_theta_z(WR1x_svd.Y(:,1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta_z(WR1x_svd.Y(:,2),obj_.P_mat);

            %% model base space trivial vector field
            [txu_S,vartheta_S] = obj_.txu(obj_,xu,R1_svd.W);

            [rank_vtheta_S,s_vtheta_S,V_vtheta_S] = fspc.rsV_unpack(vartheta_S');
            rank_vtheta_S
            s_vtheta_S_row = s_vtheta_S'
            Y_vtheta_S = V_vtheta_S.*(s_vtheta_S'/s_vtheta_S(1));

            rstats = @(r_) deal( sqrt(sum(r_)), max(r_(:)), median(r_(:)) );

            [tx_S,tu_S] = deal(txu_S(1,:), txu_S(2:end,:));
            res_tu_S = sum((tu_S-d1xu).^2,1);
            [err_tu_S_net,maxr_tu_S,medr_tu_S] = rstats(res_tu_S)

            % these are the scaled singular vectors of the Vx submatrix.
            fspc.print_vshort_polynomial_theta(Y_vtheta_S(:,1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(Y_vtheta_S(:,2),obj_.P_mat);

            K_R1 = V_vtheta_S(:,1:rank_vtheta_S); % an orthonormal basis for ker(R1)
            [txu_K_S,vartheta_K_S] = obj_.txu(obj_,xu,K_R1); % pass back to txu

            [tx_K_S,tu_K_S] = deal(txu_K_S(1,:), txu_K_S(2:end,:));
            res_tu_K_S = sum((tu_K_S-d1xu).^2,1);
            [err_tu_K_S_net,maxr_tu_K_S,medr_tu_K_S] = rstats(res_tu_K_S)

            tu_mat_S = reshape(tu_S,ndep,nobs);
            vartheta_tns_S = reshape(vartheta_K_S,Plen,nvar,nobs);

            T1_full_mat = nan(ntheta,nobs);
            G1_full_ttns = zeros(ntheta,ndep,nobs);
            G1_true_ttns = zeros(ntheta,ndep,nobs);
            for i = 1:nobs
                vtx_i = vartheta_tns_S(:,1,i);
                vtu_i = vartheta_tns_S(:,2:end,i);

                l_i = l_S(i,:)';
                g_i = gxu_l_S(:,:,i)';
                Lam1_x_i = reshape(Lam_x_dxu(:,:,1,i),ndep,Plen)';
                Lam1_u_i = Lam_u_dxu(:,1,i);

                Lam_dxu_T_i = G1_full_ttns(:,:,i);
                Lam0_xu_T_i = zeros(ntheta,nvar);

                T1_full_mat(1:Plen,i) = l_i; % tx*vx = 1*vx
                Lam0_xu_T_i(1:Plen,1) = l_i;
                idel_P = Plen;
                for idep = 1:ndep
                    Lam_dxu_T_i(1:Plen,idep) = Lam1_x_i(:,idep);

                    inds_ii = (1+idel_P):(Plen+idel_P);
                    Lam_dxu_T_i(inds_ii,idep) = Lam1_u_i;
                    Lam0_xu_T_i(inds_ii,idep+1) = l_i;
                    T1_full_mat(inds_ii,i) = d1xu(:,idep)*l_i; % tuq*vuq = dxuq*li
                    idel_P = idel_P + Plen;
                end
                G1_full_ttns(:,:,i) = ...
                ((vtu_i' - (tu_mat_S(:,i).*vtx_i'))*g_i*Lam0_xu_T_i' - Lam_dxu_T_i')';

                % oracle verification
                G1_true_ttns(:,:,i) = ...
                (obj_.obs.gradf(Smat(1,i),Smat(2:nvar,i))'*Lam0_xu_T_i' - Lam_dxu_T_i')';
            end
            T1_svd = fspc.compute_svd_package(T1_full_mat)

            fspc.print_vshort_polynomial_theta(T1_svd.V(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(T1_svd.V(:,end-1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(T1_svd.V(:,1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(T1_svd.V(:,2),obj_.P_mat);

            G1_svd = fspc.compute_svd_package(G1_full_ttns)

            fspc.print_vshort_polynomial_theta(G1_svd.W(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(G1_svd.W(:,end-1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(G1_svd.W(:,end-2),obj_.P_mat);

            [W_G1_true,rank_G1_true,s_G1_true,V_G1_true,G1_true_mat] = ...
                fspc.safely_process_net_svd( G1_true_ttns );
            rank_G1_true
            s_G1_true_row = s_G1_true'

            check = sqrt(sum( (G1_true_mat*V_G1_full(:,(end-3):end)).^2 ,1))

            fspc.print_vshort_polynomial_theta(W_G1_true(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(W_G1_true(:,end-1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(W_G1_true(:,end-2),obj_.P_mat);

            size(s_T1_full)
            YT1_full = V_T1_full.*(s_T1_full');
            YG1_full = V_G1_full.*(s_G1_full');

            size(YT1_full)
            size(YG1_full)

            [W_T1G1,rank_T1G1,s_T1G1,V_T1G1,T1G1_full_mat] = ...
                fspc.safely_process_net_svd( [ YT1_full, YG1_full ] );
            size(T1G1_full_mat)
            rank_T1G1
            s_T1G1_row = s_T1G1'

            fspc.print_vshort_polynomial_theta(W_T1G1(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(W_T1G1(:,end-1),obj_.P_mat);


            fprintf('\n(model_trivial_vfield) end\n');
            pause

            coord_specs = 0;
        end
        function [Rn_encoding,mod_specs] = model_Fode_observations(obj_,obs_,Sobs_)

            % [Rn_encoding,enc_specs] = apv.encode_R( obj_,obs_,Sobs_ )
            [Rn_encoding,enc_specs,coord_specs] = apv.model_trivial_vfield( obj_,obs_,Sobs_ )

            pause

            [eor,ndep] = deal(obj_.eor,obj_.ndep);
            [Pmat,Plen] = obj_.Pmat_Plen;

            ntheta = size(Rn_encoding{1},2);

            nobs = enc_specs.nobs;
            nset = enc_specs.nset;
            kor = enc_specs.kor_obs;
            ndim = enc_specs.ndim_obs;

            s_inv_scl = @(s_) reshape(1.0./(s_/s_(end)),1,[]);

            SVD_Rk_cell = deal(cell(3,kor,nset));
            RRk_full_ttns = nan(ntheta,ntheta,kor,nset);
            for i = 1:nset
                for k = 1:kor

                    % concatenated matrix of Rk constraints, all dependent variables together
                    Rkmat_i = reshape( permute(Rn_encoding{k,i},[2 1 3]) ,ntheta,[])';

                    % [~,SVD_Rk_cell{2,k,i},SVD_Rk_cell{3,k,i}] = svd(Rkmat_i,'econ','vector');
                    % SVD_Rk_cell{1,k,i} = rank(Rkmat_i);
                    [SVD_Rk_cell{1,k,i},SVD_Rk_cell{2,k,i},SVD_Rk_cell{3,k,i}] = ...
                        fspc.rsV_unpack(Rkmat_i);

                    % save row space representation of this Rk matrix for full null space problem
                    RRk_full_ttns(:,:,k,i) = SVD_Rk_cell{3,k,i}.*( SVD_Rk_cell{2,k,i}' );

                end
            end

            function [WA_,rA_,sA_,VA_,Afull_] = process_net_svd(Atns_)
                Afull_ = ( reshape( Atns_, size(Atns_,1),[] ) )';
                [rA_,sA_,VA_] = fspc.rsV_unpack(Afull_);
                WA_ = VA_.*s_inv_scl(sA_);
            end

            % Evaluate aggregate R1 svd: always exists, always relevant.
            [W_R1_mat,rank_R1_full,s_R1_full,V_R1_full,R1_full_mat] = ...
                process_net_svd(RRk_full_ttns(:,:,1,:));
            rank_R1_full
            s_R1_full'

            fspc.print_vshort_polynomial_theta(W_R1_mat(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(W_R1_mat(:,end-1),obj_.P_mat);

            [WR1x,WR1u] = obj_.split_Vxu_mat(W_R1_mat);
            [rank_WR1x,s_WR1x,V_WR1x] = fspc.rsV_unpack(WR1x');
            rank_WR1x
            s_WR1x'

            YW1_x = V_WR1x.*(s_WR1x'); % these are the scaled singular vectors of the Vx submatrix.
            rank_YW1_x = rank(YW1_x)

            % these are the scaled singular vectors of the Vx submatrix.
            fspc.print_vshort_polynomial_theta_z(YW1_x(:,1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta_z(YW1_x(:,2),obj_.P_mat);



            %% set outputs

            mod_specs = enc_specs;
            mod_specs.R1_full_mat = R1_full_mat;
            mod_specs.s_R1_full = s_R1_full;
            mod_specs.V_R1_full = V_R1_full;
            mod_specs.rank_R1_full = rank(R1_full_mat);

            fprintf('(model_Fode_observations end)\n');
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
                The
            %}
            [dx_l_S,gxu_l_S,L_S,l_S] = obj_.dx_l(obj_,xumat,d1xumat);

            ntheta = nvar*Plen;

            ncol_Lambda_u = ntheta-Plen;
            i_imm = zeros(ncol_Lambda_u,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            i_imm_vec = logical( i_imm(:) );

            len_Lambda_u = prod(size(i_imm));
            l_imm_init = zeros(len_Lambda_u,1);
            ones_imm = ones(ndep,Plen);

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

            if (kor_obs>1)

                % To do. Second order ratio condition is almost free, only needs first total deriv. x

            end

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
