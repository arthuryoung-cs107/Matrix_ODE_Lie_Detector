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
        function [dx_out,gxu_out,L_out] = dx_l_polynomial(obj,xu_,dxu_)
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
        end
    end
    methods (Static)

        function [Rn_encoding,mod_specs] = model_Fode_observations(obj_,obs_,Sobs_)

            [Rn_encoding,enc_specs] = apv.encode_R( obj_,obs_,Sobs_ )

            [eor,ndep] = deal(obj_.eor,obj_.ndep);
            % ndim = 1 + ndep*(eor+1); % dimension of n'th jet space
            [Pmat,Plen] = obj_.Pmat_Plen;

            ntheta = size(Rn_encoding{1},2)

            nobs = enc_specs.nobs;
            nset = enc_specs.nset;
            kor = enc_specs.kor_obs;
            ndim = enc_specs.ndim_obs;

            % R1_tns = enc_specs.R1_tns;

            function [r_,s_,V_] = svd_unpack(A_)
                [~,S_,V_] = svd(A_,'econ');
                r_ = rank(A_);
                s_ = diag(S_);
            end

            SVD_Rk_cell = cell(3,kor,nset);
            RRk_full_ttns = nan(ntheta,ntheta,kor,nset);
            for i = 1:nset
                for k = 1:kor

                    % concatenated matrix of Rk constraints, all dependent variables together
                    Rkmat_i = reshape( permute(Rn_encoding{k,i},[2 1 3]) ,ntheta,[])';

                    % [~,SVD_Rk_cell{2,k,i},SVD_Rk_cell{3,k,i}] = svd(Rkmat_i,'econ','vector');
                    % SVD_Rk_cell{1,k,i} = rank(Rkmat_i);
                    [SVD_Rk_cell{1,k,i},SVD_Rk_cell{2,k,i},SVD_Rk_cell{3,k,i}] = svd_unpack(Rkmat_i);

                    % save row space representation of this matrix for full null space problem
                    RRk_full_ttns(:,:,k,i) = SVD_Rk_cell{3,k,i}.*reshape( SVD_Rk_cell{2,k,i},1,[] );

                end
            end

            % Evaluate aggregate R1 svd, always exists, always relevant.
            R1_full_mat = ( reshape(RRk_full_ttns(:,:,1,:),ntheta,[]) )';
            % [~,s_R1_full,V_R1_full] = svd( R1_full_mat ,'econ','vector');
            % rank_R1_full = rank(R1_full_mat)
            [rank_R1_full,s_R1_full,V_R1_full] = svd_unpack(R1_full_mat);

            W1_mat = V_R1_full.*(1.0./(s_R1_full/s_R1_full(end)))';

            fspc.print_short_polynomial_theta(W1_mat(:,end),obj_.P_mat);
            fspc.print_short_polynomial_theta(W1_mat(:,end-1),obj_.P_mat);

            % split_Vxu_mat = @(V_) deal( V_(1:Plen,:) , V_((Plen+1):end,:) );
            [W1x,W1u] = obj_.split_Vxu_mat(W1_mat);

            % svd of W1x matrix
            % [~,s_W1x,V_W1x] = svd(W1x','econ','vector');
            % rank_W1x = rank(W1x)
            [rank_W1x,s_W1x,V_W1x] = svd_unpack(W1x');

            s_W1x_row = s_W1x'

            WW1_x = V_W1x.*(s_W1x'); % these are the scaled singular vectors of the Vx submatrix.
            rank_WW1_x = rank(WW1_x)

            % these are the scaled singular vectors of the Vx submatrix.
            fspc.print_short_polynomial_theta_z(WW1_x(:,1),obj_.P_mat);
            fspc.print_short_polynomial_theta_z(WW1_x(:,2),obj_.P_mat);


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

            [dx_l_S,gxu_l_S,L_S] = obj_.dx_l(obj_,xumat,d1xumat);
            l_S = (reshape(prod(L_S,1),Plen,nobs))';

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
            end

            R1_tns = nan(ndep,ntheta,nobs);
            for i = 1:nobs
                % check1 = (-d1xumat(:,i)*l_S(i,:)) % also works
                R1_tns(:,:,i) = [ (-d1xumat(:,i).*l_S(i,:)) , immerse_lambda(l_S(i,:)) ];
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
            'L_S', L_S ...
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
