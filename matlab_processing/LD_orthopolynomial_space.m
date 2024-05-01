classdef LD_orthopolynomial_space < LD_power_space
    properties (Constant)
        h_range = [-1, 1];
    end
    properties
        poly_coeffs;
        icoeff_mat;

        fmap_m;
        fmap_b;
        bmap_m;
        bmap_b;
    end
    methods
        function obj = LD_orthopolynomial_space(bor_, meta_)
            icoeff_mat = zeros(bor_+1);
            for j = 0:bor_
                jj = j + 1;
                icoeff_mat(jj,1) = 1;
                jpow = double(j);
                for i = 1:j
                    ii = i+1;
                    icoeff_mat(jj,ii) = icoeff_mat(jj,ii-1)*jpow;
                    jpow = jpow-1;
                end
            end

            obj@LD_power_space(bor_,meta_);

            obj.poly_coeffs = zeros(bor_+1);
            obj.icoeff_mat = icoeff_mat;

            obj.fmap_m = ones(1,1+obj.ndep);
            obj.fmap_b = zeros(1,1+obj.ndep);
            obj.bmap_m = ones(1,1+obj.ndep);
            obj.bmap_b = zeros(1,1+obj.ndep);
        end
        function obj_out = read_domain_configuration(obj,name_)
            [ndep nvar combo_len] = deal(obj.ndep, obj.ndep+1, obj.combo_len);
            fprintf('(orthopolynomial_basis::read_domain_configuration) reading %s\n', name_);
            file = fopen(name_);

            hlen = fread(file,1,'int=>int');
            header_rest = fread(file,hlen,'int=>int');
            [ord_mat_len map_len len_sym] = deal(header_rest(1), header_rest(2), header_rest(3));
            ord_mat_flat = fread(file,ord_mat_len,'int=>int');
            map_params_flat = fread(file,map_len,'double=>double');
            poly_coeffs_flat = fread(file,len_sym,'double=>double');
            fclose(file);

            %% ignore other data for now, only need the map params
            % ord_mat = double(reshape(ord_mat_flat,nvar,combo_len)');
            mparams_mat = reshape(map_params_flat,nvar,4)';

            obj_out = obj;
            % obj_out.order_mat = ord_mat;
            obj_out.fmap_m = mparams_mat(1,:);
            obj_out.fmap_b = mparams_mat(2,:);
            obj_out.bmap_m = mparams_mat(3,:);
            obj_out.bmap_b = mparams_mat(4,:);
        end
        function write_domain_configuration(obj,name_)
            hlen = 3;
            nvar = 1+obj.ndep;
            ord_mat_len = obj.combo_len*nvar;
            map_len = 4*nvar;
            borp1 = obj.bor + 1;
            len_sym = (borp1*(borp1+1))/2;
            coeff_vec_out = zeros(len_sym,1);
            icoeff = 1;
            for i = 0:(obj.bor)
                for j = 0:i
                    coeff_vec_out(icoeff) = obj.poly_coeffs(i+1,j+1);
                    icoeff = icoeff + 1;
                end
            end
            ord_vec = reshape((obj.order_mat)',[],1);
            mparams_mat = [obj.fmap_m'; obj.fmap_b'; obj.bmap_m'; obj.bmap_b'];
            mparams_vec = mparams_mat(:);

            header = [hlen,ord_mat_len,map_len,len_sym]';
            file = fopen(name_,'w+');
            fwrite(file, header, 'int');
            fwrite(file, ord_vec, 'int');
            fwrite(file, mparams_vec, 'double');
            fwrite(file, coeff_vec_out, 'double');
            fclose(file);
            fprintf('(orthopolynomial_basis::write_ode_experiment_data) wrote file: %s \n', name_);
        end
        function De_out = dp_de_ratio_n1_nfull_stabilized_fast(obj,xu_,H_)
            [eor,ndep,nvar,ndim,npts] = deal(obj.eor,obj.ndep,1 + obj.ndep, obj.ndim, size(xu_,2));
            [De_n1,params_n1] = obj.dp_de_ratio_n1_stabilized_fast_itemized(xu_,H_);
            De_out = [De_n1; nan(ndim-nvar,npts)];
            for ipts = 1:npts
                dnxu_i = reshape(De_out(2:end,ipts),ndep,eor+1);
                for k = 1:eor
                    vik = obj.dp_de_Hmat([xu_(:,ipts); reshape(dnxu_i(:,1:k),[],1)],params_n1(:,ipts),k);
                    dnxu_i(:,k+1) = vik((end-ndep+1):end);
                end
                De_out(2:end,ipts) = reshape(dnxu_i,[],1);
            end
        end
        function De_out = dp_de_ratio_n1_stabilized_fast(obj,p_,params_)
            [ndep,nvar,npts,nparams,combo_len,bor,poly_coeffs,order_mat] = deal(obj.ndep,1+obj.ndep,size(p_,2),size(params_,2),obj.combo_len,double(obj.bor),obj.poly_coeffs,obj.order_mat);
            xu_raw_mat = p_(1:nvar,:);
            xu_mat = ((obj.fmap_m.*xu_raw_mat') + obj.fmap_b)';
            Theta_x_mat = (params_(1:combo_len,:))';
            Theta_u_mat_T = params_((combo_len+1):end,:);
            Lxu_inds = order_mat + [0, (bor+1)*(1:ndep)] + 1;
            De_u_out = nan(ndep,npts);
            for ipts = 1:npts
                xu_pow = poly_coeffs*((xu_mat(:,ipts).^(0:bor))');
                Lvec_i = prod(xu_pow(Lxu_inds),2);
                vx_i_vec = Theta_x_mat*Lvec_i;
                De_u_out(:,ipts) = ((Lvec_i')*reshape((Theta_u_mat_T*vx_i_vec)/(vx_i_vec'*vx_i_vec),combo_len,ndep))';
            end
            De_out = [ones(1,npts) ; De_u_out];
        end
        function [De_out,params_out] = dp_de_ratio_n1_stabilized_fast_itemized(obj,p_,params_)
            [ndep,nvar,npts,nparams,combo_len,bor,poly_coeffs,order_mat] = deal(obj.ndep,1+obj.ndep,size(p_,2),size(params_,2),obj.combo_len,double(obj.bor),obj.poly_coeffs,obj.order_mat);
            xu_raw_mat = p_(1:nvar,:);
            xu_mat = ((obj.fmap_m.*xu_raw_mat') + obj.fmap_b)';
            Theta_x_mat = (params_(1:combo_len,:))';
            Lxu_inds = order_mat + [0, (bor+1)*(1:ndep)] + 1;
            De_u_out = nan(ndep,npts);
            params_out = nan(obj.ndof_full,npts);
            for ipts = 1:npts
                xu_pow = poly_coeffs*((xu_mat(:,ipts).^(0:bor))');
                Lvec_i = prod(xu_pow(Lxu_inds),2);
                vx_i_vec = Theta_x_mat*Lvec_i;
                params_out(:,ipts) = (params_*vx_i_vec)/(vx_i_vec'*vx_i_vec);
                De_u_out(:,ipts) = ((Lvec_i')*reshape(params_out((combo_len+1):end,ipts),combo_len,ndep))';
            end
            De_out = [ones(1,npts) ; De_u_out];
        end
        function [De_out,dnxu_update] = dp_de_ratio_stabilized_update(obj,p_,params_)
            [eor,ndep,nvar,ndim,npts] = deal(obj.eor,obj.ndep,1 + obj.ndep, obj.ndim, size(p_,2));
            xu = p_(1:nvar,:);
            dnxu_old = p_((nvar+1):end,:);
            [De_n1,params_n1] = obj.dp_de_ratio_n1_stabilized_fast_itemized(xu,params_);
            De_out = nan(size(p_));
            De_out(1:nvar,:) = De_n1;
            dnxu_update = nan(obj.ndim-nvar,npts);
            dnxu_update(1:ndep,:) = 0.5*(dnxu_old(1:ndep,:)+De_n1(2:nvar,:));
            for ipts = 1:npts
                for k = 2:eor
                    vik = obj.dp_de_Hmat([xu(:,ipts) ; dnxu_update(1:(ndep*(k-1)),ipts)],params_n1(:,ipts),k-1);
                    i_update = (ndep*(k-1) + 2):(ndep*k + 1);
                    De_out(i_update,ipts) = vik(i_update);
                    dnxu_update(((ndep*(k-1))+1):(ndep*k),ipts) = 0.5*(dnxu_old((ndep*(k-1)+1):(ndep*k),ipts) + vik(i_update));
                end
            end
        end
        function [De_out,dnxu_update] = dp_de_ratio_stabilized_brgman(obj,p_,p_og_,params_)
            [eor,ndep,nvar,ndim,npts] = deal(obj.eor,obj.ndep,1 + obj.ndep, obj.ndim, size(p_,2));
            xu = p_(1:nvar,:);
            dnxu_old = p_((nvar+1):end,:);
            dnxu_og = p_og_((nvar+1):end,:);
            [De_n1,params_n1] = obj.dp_de_ratio_n1_stabilized_fast_itemized(xu,params_);
            De_out = nan(size(p_));
            De_out(1:nvar,:) = De_n1;
            dnxu_update = nan(obj.ndim-nvar,npts);
            % dnxu_update(1:ndep,:) = dnxu_og(1:ndep,:)+(dnxu_old(1:ndep,:)-De_n1(2:nvar,:));
            dnxu_update(1:ndep,:) = De_n1(2:nvar,:);
            for ipts = 1:npts
                for k = 2:eor
                    % vik = obj.dp_de_Hmat([xu(:,ipts) ; dnxu_old(1:(ndep*(k-1)),ipts)],params_n1(:,ipts),k-1);
                    % vik = obj.dp_de_Hmat([xu(:,ipts) ; dnxu_update(1:(ndep*(k-1)),ipts)],params_n1(:,ipts),k-1);
                    % vik = obj.dp_de_Hmat([xu(:,ipts) ; dnxu_update(1:(ndep*(k-1)),ipts)],params_n1(:,ipts),k-1);
                    vik = obj.dp_de_Hmat([xu(:,ipts) ; 0.5*(dnxu_old(1:(ndep*(k-1)),ipts)+dnxu_update(1:(ndep*(k-1)),ipts))],params_n1(:,ipts),k-1);

                    i_update = (ndep*(k-1) + 2):(ndep*k + 1);
                    De_out(i_update,ipts) = vik(i_update);
                    % dnxu_update(((ndep*(k-1))+1):(ndep*k),ipts) = dnxu_og((ndep*(k-1)+1):(ndep*k),ipts) + (dnxu_old((ndep*(k-1)+1):(ndep*k),ipts) - vik(i_update));
                    dnxu_update(((ndep*(k-1))+1):(ndep*k),ipts) = De_out(i_update,ipts);
                end
            end
        end
        function [De_out,v_out,r_out,w_out,params_out] = dp_de_ratio_n1_stabilized_itemized(obj,p_,K_)
            nvar = 1+obj.ndep;
            [param_len,Kdim] = size(K_);
            npts = size(p_,2);
            De_out = nan(nvar,npts);
            [r_out,v_out] = deal(nan(nvar,Kdim,npts));
            w_out = nan(Kdim,npts);
            params_out = nan(param_len,npts);
            for i = 1:npts
                v_out(:,:,i) = obj.dp_de_Hmat(p_(1:nvar,i),K_,1);
                r_out(:,:,i) = v_out(:,:,i)./v_out(1,:,i);
                w_out(:,i) = ( v_out(1,:,i)/( sum( (v_out(1,:,i)).^2 ) ) )';
                params_out(:,i) = K_*w_out(:,i);
                De_out(:,i) = v_out(:,:,i)*w_out(:,i);
            end
        end
        function [De_out,v_out,r_out,w_out,params_out] = dp_de_ratio_stabilized_itemized(obj,p_,K_)
            [param_len,Kdim] = size(K_);
            [ndof_ODE,npts] = size(p_);
            De_out = nan(ndof_ODE,npts);
            [r_out,v_out] = deal(nan(ndof_ODE,Kdim,npts));
            w_out = nan(Kdim,npts);
            params_out = nan(param_len,npts);
            for i = 1:npts
                v_out(:,:,i) = obj.dp_de_Hmat(p_(:,i),K_,obj.eor);
                r_out(:,:,i) = v_out(:,:,i)./v_out(1,:,i);
                w_out(:,i) = ( v_out(1,:,i)/( sum( (v_out(1,:,i)).^2 ) ) )';
                params_out(:,i) = K_*w_out(:,i);
                De_out(:,i) = v_out(:,:,i)*w_out(:,i);
            end
        end
        function [dls_out v_forw_eval V_mat w_vec J_V_tns] = dls_eval(obj,p_,H_)
        % function [dls_out v_forw_eval V_mat w_vec J_V_tns J_V_ttns] = dls_eval(obj,p_,H_)
            % [V_mat J_V_tns J_V_ttns] = obj.dp_de_Hmat_Jac(p_,H_,obj.eor);
            [V_mat J_V_tns] = obj.dp_de_Hmat_Jac(p_,H_,obj.eor);
            vt_vec = V_mat(1,:)';
            vdnm1_mat = V_mat((obj.ndim-(2*obj.ndep)+1):(obj.ndim-obj.ndep),:)';

            del_vt_mat = reshape(J_V_tns(1,:,:),obj.ndim,[])';
            del_vdnm1_tns = permute(reshape(J_V_tns((obj.ndim-(2*obj.ndep)+1):(obj.ndim-obj.ndep),:,:),obj.ndep,obj.ndim,[]),[3 2 1]);

            w_vec = vt_vec/sum(sum(vt_vec.*vt_vec));
            v_forw_eval = V_mat*w_vec;
            % v_forw_eval = [1; (V_mat(2:end,:)*w_vec)];
            % v_d1x_check = [(reshape(J_V_tns(2,:,:),obj.ndim,[])' - p_(3)*(reshape(J_V_tns(1,:,:),obj.ndim,[])')) * [1; p_(3); p_(4); 0], V_mat(3,:)'];
            % v_d2x_check = [(reshape(J_V_tns(3,:,:),obj.ndim,[])' - p_(4)*(reshape(J_V_tns(1,:,:),obj.ndim,[])')) * [1; p_(3); p_(4); 0], V_mat(4,:)'];

            del_denom = (vt_vec'*vt_vec)*(vt_vec'*vt_vec);

            dls_out = ones(obj.ndim,obj.ndep);
            for idep = 1:obj.ndep
                for idim = 1:(obj.ndim-obj.ndep)
                    del_numer = vt_vec'*((vt_vec'*del_vdnm1_tns(:,idim,idep) + vdnm1_mat(:,idep)'*del_vt_mat(:,idim))-2.0*(del_vt_mat(:,idim)*vdnm1_mat(:,idep)'))*vt_vec;
                    dls_out(idim,idep) = -1.0*del_numer/del_denom;
                end
            end
        end
        function De_out = dp_de_ratio_nk(obj,p_,h_,k_)
            npts = size(p_,2);
            if (npts>1)
                De_out = nan(size(p_));
                for i = 1:npts
                    v = obj.dp_de(p_(:,i),h_,k_);
                    De_out(:,i) = [1.0; p_((obj.ndep+2):end,i); v((length(v)-obj.ndep+1):end)/v(1)];
                end
            else
                v = obj.dp_de(p_,h_,k_);
                De_out = [1.0; p_((obj.ndep+2):end); v((length(v)-obj.ndep+1):end)/v(1)];
            end
        end
        function De_out = dp_de_proj(obj,p_,h_)
            npts = size(p_,2);
            if (npts>1)
                De_out = nan(size(p_));
                for i = 1:npts
                    v = obj.dp_de(p_(:,i),h_);
                    v_short = v(1:(obj.ndim-obj.ndep));
                    d_short = [1.0; p_((obj.ndep+2):end,i)];
                    % De_out(:,i) = v*((v_short/norm(v_short))'*(d_short/norm(d_short)));
                    De_out(:,i) = v*((v_short/norm(v_short))'*(d_short));
                end
            else
                v = obj.dp_de(p_,h_);
                v_short = v(1:(obj.ndim-obj.ndep));
                d_short = [1.0; p_((obj.ndep+2):end)];
                % De_out = v*((v_short/norm(v_short))'*(d_short/norm(d_short)));
                De_out = v*((v_short/norm(v_short))'*(d_short));
            end
        end
        function De_out = dp_de_ratio_n1_nfull_stabilized(obj,p_,H_)
            vn1 = obj.dp_de_ratio_n1_stabilized(p_,H_);
            [eor ndep ndim] = deal(obj.eor,obj.ndep,obj.ndim);
            npts = size(p_,2);
            De_out = nan(ndim, npts);
            De_out(1:(1+ndep),:) = vn1;
            if (npts>1)
                for i = 1:npts
                    for k = 2:(eor+1)
                        vk = obj.dp_de_Hmat([ p_(:,i); De_out(2:(ndep*(k-1) + 1),i)],H_,k-1);

                        wvec = (vk(1,:).*vk(1,:))'; %% redundant calc, probably want to recycle
                        vdkt_use = (vk((size(vk,1)-ndep+1):end,:)./vk(1,:))*(wvec/sum(wvec));

                        De_out((ndep*(k-1) + 2):(ndep*k + 1),i) = vdkt_use;
                    end
                end
            else
                for k = 2:(eor+1)
                    vk = obj.dp_de_Hmat([ p_; De_out(2:(ndep*(k-1) + 1))],H_,k-1);

                    wvec = (vk(1,:).*vk(1,:))'; %% redundant calc, probably want to recycle
                    vdkt_use = (vk((size(vk,1)-ndep+1):end,:)./vk(1,:))*(wvec/sum(wvec));

                    De_out((ndep*(k-1) + 2):(ndep*k + 1)) = vdkt_use;
                end
            end
        end
        function De_out = dp_de_ratio_n1_nfull(obj,p_,h_)
            vn1 = obj.dp_de_ratio_n1(p_,h_);
            [eor ndep ndim] = deal(obj.eor,obj.ndep,obj.ndim);
            npts = size(p_,2);
            De_out = nan(ndim, npts);
            De_out(1:(1+ndep),:) = vn1;
            if (npts>1)
                for i = 1:npts
                    for k = 2:(eor+1)
                        vk = obj.dp_de([ p_(:,i); De_out(2:(ndep*(k-1) + 1),i)],h_,k-1);
                        De_out((ndep*(k-1) + 2):(ndep*k + 1),i) = vk((length(vk)-ndep+1):end)/vk(1);
                    end
                end
            else
                for k = 2:(eor+1)
                    vk = obj.dp_de([ p_; De_out(2:(ndep*(k-1) + 1))],h_,k-1);
                    De_out((ndep*(k-1) + 2):(ndep*k + 1)) = vk((length(vk)-ndep+1):end)/vk(1);
                end
            end
        end
        function De_out = dp_de_ratio_stabilized(obj,p_,H_)
            [ndof_ODE npts] = size(p_);
            [eor_set,inds_use] = deal(max([obj.eor-1 1]), (2 + obj.ndep*(obj.eor-1) ):(1 + obj.ndep*obj.eor) );
            if (npts>1)
                De_out = nan(ndof_ODE,npts);
                if (ndof_ODE == obj.ndim)
                    dn_inds = (ndof_ODE-obj.ndep+1):ndof_ODE;
                    for i = 1:npts
                        v = obj.dp_de_Hmat(p_(:,i),H_,obj.eor);
                        wvec = (v(1,:).*v(1,:))';
                        vdtn_use = (v(dn_inds,:)./v(1,:))*(wvec/sum(wvec));
                        De_out(:,i) = [1.0; p_((obj.ndep+2):end,i); vdtn_use];
                    end
                else
                    for i = 1:npts
                        v = obj.dp_de_Hmat(p_(:,i),H_,eor_set);
                        wvec = (v(1,:).*v(1,:))';
                        vdtn_use = (v(inds_use,:)./v(1,:))*(wvec/sum(wvec));
                        De_out(:,i) = [1.0; p_((obj.ndep+2):inds_use(end),i); vdtn_use];
                    end
                end
            else
                v = obj.dp_de_Hmat(p_,H_,eor_set);
                wvec = (v(1,:).*v(1,:))';
                vdtn_use = (v(inds_use,:)./v(1,:))*(wvec/sum(wvec));
                De_out = [1.0; p_((obj.ndep+2):inds_use(end)); vdtn_use];
            end
        end
        function De_out = dp_de_ratio_n1_stabilized(obj,p_,H_)
            [ndof_ODE npts] = size(p_);
            if (npts>1)
                De_out = nan(ndof_ODE,npts);
                for i = 1:npts
                    v = obj.dp_de_Hmat([p_(:,i); ones(obj.ndep,1)],H_,1); % dummy values for n = 1 derivative
                    wvec = (v(1,:).*v(1,:))';
                    vdt_use = (v(2:(1+obj.ndep),:)./v(1,:))*(wvec/sum(wvec));
                    De_out(:,i) = [1.0; vdt_use];
                end
            else
                v = obj.dp_de_Hmat([p_; ones(obj.ndep,1)],H_,1); % dummy values for n = 1 derivative
                wvec = (v(1,:).*v(1,:))';
                vdt_use = (v(2:(1+obj.ndep),:)./v(1,:))*(wvec/sum(wvec));
                De_out = [1.0; vdt_use];
            end
        end
        function De_out = dp_de_ratio_n1(obj,p_,h_)
            npts = size(p_,2);
            if (npts>1)
                De_out = nan(size(p_));
                for i = 1:npts
                    v = obj.dp_de([p_(:,i); ones(obj.ndep,1)],h_,1); % dummy values for n = 1 derivative
                    De_out(:,i) = [1.0; v(2:(1+obj.ndep))/v(1)]; % evaluate ratio for first derivative
                end
            else
                v = obj.dp_de([p_; ones(obj.ndep,1)],h_,1); % dummy values for n = 1 derivative
                De_out = [1.0; v(2:(1+obj.ndep))/v(1)]; % evaluate ratio for first derivative
            end
        end
        function De_out = dp_de_ratio(obj,p_,h_)
            npts = size(p_,2);
            eor_set = max([obj.eor-1 1]);
            if (npts>1)
                De_out = nan(size(p_));
                for i = 1:npts
                    v = obj.dp_de(p_(:,i),h_,eor_set);
                    De_out(:,i) = [1.0; p_((obj.ndep+2):end,i); v((length(v)-obj.ndep+1):end)/v(1)];
                end
            else
                v = obj.dp_de(p_,h_,eor_set);
                De_out = [1.0; p_((obj.ndep+2):end); v((length(v)-obj.ndep+1):end)/v(1)];
            end
        end
        function [nc_out np_out] = get_def_nc_np(obj)
            [nc_out np_out] = get_def_nc_np@tangent_space_basis(obj);
        end
        function coeff_out = get_coeff(obj,k_,j_,coeffs_in_)
            if (nargin==4)
                coeffs_in = coeffs_in_;
            else
                coeffs_in = obj.poly_coeffs;
            end

            if ((j_>=0)&&(j_<=k_))
                coeff_out = coeffs_in(k_+1,j_+1);
            else
                coeff_out = 0.0;
            end
        end
        function obj_out = configure_self(obj,ldd_)
            [ndep nvar combo_len] = deal(obj.ndep, obj.ndep+1, obj.combo_len);
            input_name = [ldd_.base_dir_name ldd_.basis_dir obj.name num2str(obj.bor) '_domain_config.' LD_demo_data.LD_data_extension];
            fprintf('(orthopolynomial_basis::configure_self) reading %s\n', input_name);
            file = fopen(input_name);

            hlen = fread(file,1,'int=>int');
            header_rest = fread(file,hlen,'int=>int');
            [ord_mat_len map_len len_sym] = deal(header_rest(1), header_rest(2), header_rest(3));
            ord_mat_flat = fread(file,ord_mat_len,'int=>int');
            map_params_flat = fread(file,map_len,'double=>double');
            poly_coeffs_flat = fread(file,len_sym,'double=>double');
            fclose(file);

            %% ignore other data for now, only need the map params
            % ord_mat = double(reshape(ord_mat_flat,nvar,combo_len)');
            mparams_mat = reshape(map_params_flat,nvar,4)';

            obj_out = obj;
            % obj_out.order_mat = ord_mat;
            obj_out.fmap_m = mparams_mat(1,:);
            obj_out.fmap_b = mparams_mat(2,:);
            obj_out.bmap_m = mparams_mat(3,:);
            obj_out.bmap_b = mparams_mat(4,:);
        end
        function De_out = dp_de_time(obj,p_,Hmat_)
            De_H = obj.dp_de_Hmat(p_,Hmat_);

            [ptall nbse nvar ndim] = deal(reshape(p_,[],1), size(Hmat_,2), obj.ndep+1, obj.ndim);
            [p_tx p_xp] = deal(ptall(1:nvar), ptall((nvar+1):end));

            Emat = De_H(1:nbse,:);
            tau = [1; p_xp];
            beta = linsolve(Emat,tau);

            %% we have a few choices now for how to compute De_out

            %% option 1: take the raw linear combination of the computed epsilon vectors
            % De_out = De_H*beta;

            %% option 2: use the input values accordingly, then use linear combination of final dimensions
            De_out = nan(ndim,1);
            De_out(1) = 1.0;
            De_out(2:nbse) = p_xp;
            De_out((nbse+1):end) = De_H((nbse+1):end,:)*beta;
        end
        function De_out_pts = dp_de_ptsmat(obj,ptsmat_,Hmat_)
            [ndim npts] = size(ptsmat_);
            nH = size(Hmat_,2);
            if (nH==1)
                De_out_pts = nan(ndim,npts);
                for i = 1:npts
                    De_out_pts(:,i) = obj.dp_de(ptsmat_(:,i),Hmat_);
                end
            else
                De_out_pts = nan(ndim,nH,npts);
                for i = 1:npts
                    De_out_pts(:,:,i) = obj.dp_de_Hmat(ptsmat_(:,i),Hmat_);
                end
            end
        end
        function [De_out_H Jac_v_out Jac_L_tx_h_out] = dp_de_Hmat_Jac(obj,p_,Hmat_,kcap_)
            if (nargin == 3)
                kcap = obj.eor;
            else
                kcap = kcap_;
            end
            [eor ndep ndim ndof combo_len bor] = deal(kcap, obj.ndep, 1+(obj.ndep*(kcap+1)), obj.ndof, obj.combo_len, double(obj.bor));
            nvar = ndep + 1;
            pflat = reshape(p_,1,[]);
            Htall = reshape(Hmat_,ndof,[]); %% ensuring we have the right dimensions
            nbse = size(Htall,2); % number of solutions to check

            H_all_tens = reshape(Htall,combo_len,nvar,nbse);
            ht_mat = reshape(H_all_tens(:,1,:), combo_len, nbse);
            Hx_tens = reshape(H_all_tens(:,2:end,:), combo_len, ndep, nbse);
            fmm = obj.fmap_m;
            fmb = obj.fmap_b;

            ptx_raw = pflat(1:(ndep+1));
            ptx = fmm.*ptx_raw + fmb;
            pxp = pflat((ndep+2):end);
            [order_mat, poly_coeffs, icoeff_mat] = deal(obj.order_mat, obj.poly_coeffs, obj.icoeff_mat);

            txP = ptx.^((0:bor)');
            Ltx = comp_Lt_Lx(combo_len,bor,txP,poly_coeffs,order_mat);

            De_tx_H = nan(nvar,nbse);
            for ibse = 1:nbse
                De_tx_H(:,ibse) = (H_all_tens(:,:,ibse)')*prod(Ltx,2); %% compute indep. and dep. var terms
            end

            dtxc_init = zeros(1,ndep+1);
            De_xp_dxde_mat = zeros(eor,combo_len);
            coeffs_init = ones(combo_len,1);
            Ptx_init = order_mat;
            Otx_init = order_mat;
            Ltxh_init = Ltx;
            Pxp_init = zeros(1,eor*ndep);

            specs = struct( 'eor',eor, ...
                            'ndep',ndep, ...
                            'bor', bor, ...
                            'p',pflat, ...
                            'txP',txP, ...
                            'fmm',fmm, ...
                            'pc',poly_coeffs, ...
                            'ic',icoeff_mat);

            Jac_vx = zeros(ndim,ndim,nbse);
            [De_xp_dxde_mat Jac_vx] = compute_dxde_scalar_Jac(specs,1,dtxc_init,De_xp_dxde_mat,Jac_vx,H_all_tens,coeffs_init,Otx_init,Ptx_init,Ltxh_init,Pxp_init);

            De_xp_dxde_H = nan(eor*ndep,nbse);
            for ibse = 1:nbse
                De_xp_dxde_H(:,ibse) = reshape((De_xp_dxde_mat*Hx_tens(:,:,ibse))',[],1);
            end

            De_xp_dtde_mat = zeros(eor*ndep,combo_len);
            Jac_vt = zeros(ndim,ndim,nbse);
            Jac_L_tx_h_out = zeros(1+ndep,combo_len);

            bor_p1 = bor + 1;

            nDzi = obj.not_D_zero_inds;
            idt = nDzi(:,1);
            if (sum(idt))
                dtc_p1 = dtxc_init(1) + 1;
                idt_p1 = dtc_p1 + 1;
                idt_front = 1:(bor_p1-dtc_p1);
                idt_back = idt_p1:bor_p1;

                Otx_dt = Otx_init(idt,:);
                c_dt = (-1.0*coeffs_init(idt))*fmm(1);

                Ptx_dt = [Ptx_init(idt,1)-1, Ptx_init(idt,2:end)];

                dtLtxh = Ltxh_init(idt,:);
                dtLtxh(:,1) = comp_dnL(dtc_p1,txP(idt_front,1),poly_coeffs(idt_back,idt_back),icoeff_mat(idt_back,idt_p1),Otx_dt(:,1));

                dt_dtxc = dtxc_init;
                dt_dtxc(1) = dtc_p1;

                Jac_t_t = ((c_dt).*prod(dtLtxh,2))';
                Jac_t_t_h = -1.0*Jac_t_t;
                Jac_vt(1,1,:) = pagemtimes(Jac_t_t_h,H_all_tens(idt,1,:));
                Jac_L_tx_h_out(1,idt) = Jac_t_t_h;

                for kder = 1:eor
                    kder_inds = ((ndep*(kder-1))+1):(ndep*kder);
                    De_xp_dtde_mat(kder_inds,idt) = De_xp_dtde_mat(kder_inds,idt) + (pxp(kder_inds)')*Jac_t_t;
                    [De_xp_dtde_mat(:,idt) Jac_vt] = compute_dxde_vector_Jac(specs,kder+1,kder,dt_dtxc,De_xp_dtde_mat(:,idt),Jac_vt,H_all_tens(idt,:,:),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_init);
                end
            end

            for ix = 1:ndep
                idxi = nDzi(:,ix+1);
                if (sum(idxi))
                    dxic_p1 = dtxc_init(ix+1) + 1; % increase the xi derivative count by one
                    idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                    idxi_front = 1:(bor_p1-dxic_p1);
                    idxi_back = idxi_p1:bor_p1;

                    Otx_dxi = Otx_init(idxi,:);
                    c_dxi = (-1.0*coeffs_init(idxi))*fmm(ix+1); % update coefficients

                    Ptx_dxi = Ptx_init(idxi,:); % initialize power records to selected previous values
                    Ptx_dxi(:,ix+1) = Ptx_dxi(:,ix+1)-1; % initialize powers of xi

                    dxiLtxh = Ltxh_init(idxi,:); % initialize L functional history to selected previous values
                    % update L functional history to reflect derivative in t
                    dxiLtxh(:,ix+1) = comp_dnL(dxic_p1,txP(idxi_front,ix+1),poly_coeffs(idxi_back,idxi_back),icoeff_mat(idxi_back,idxi_p1),Otx_dxi(:,ix+1));

                    dxi_dtxc = dtxc_init;
                    dxi_dtxc(ix+1) = dxic_p1;

                    Jac_t_xi = ((c_dxi.*prod(dxiLtxh,2))');
                    Jac_t_xi_h = -1.0*Jac_t_xi;
                    Jac_vt(1,ix+1,:) = pagemtimes(Jac_t_xi_h,H_all_tens(idxi,1,:));
                    Jac_L_tx_h_out(ix+1,idxi) = Jac_t_xi_h;

                    Pxp_dxi = Pxp_init; % initialize xp power records to previous values
                    Pxp_dxi(ix) = Pxp_init(ix) + 1; % update xp power records to reflect derivative in xi
                    for kder = 1:eor
                        kder_inds = ((ndep*(kder-1))+1):(ndep*kder);
                        De_xp_dtde_mat(kder_inds,idxi) = De_xp_dtde_mat(kder_inds,idxi) + ((pxp(kder_inds)')*Jac_t_xi)*pxp(ix); % add contribution of freshly computed t derivative term
                        [De_xp_dtde_mat(:,idxi) Jac_vt] = compute_dxde_vector_Jac(specs,kder+1,kder,dxi_dtxc,De_xp_dtde_mat(:,idxi),Jac_vt,H_all_tens(idxi,:,:),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
                    end
                end
            end

            Jac_v_out = Jac_vt + Jac_vx;

            De_xp_dtde_H = De_xp_dtde_mat*ht_mat;

            De_xp_H = De_xp_dxde_H + De_xp_dtde_H;

            De_out_H = [De_tx_H; De_xp_H];
        end
        function De_out_H = dp_de_Hmat(obj,p_,Hmat_,kcap_)
            if (nargin == 3)
                kcap = obj.eor;
            else
                kcap = kcap_;
            end
            [eor ndep ndof combo_len bor] = deal(kcap, obj.ndep, obj.ndof, obj.combo_len, double(obj.bor));
            nvar = ndep + 1;
            pflat = reshape(p_,1,[]);
            Htall = reshape(Hmat_,ndof,[]); %% ensuring we have the right dimensions
            nbse = size(Htall,2); % number of solutions to check

            H_all_tens = reshape(Htall,combo_len,nvar,nbse);
            ht_mat = reshape(H_all_tens(:,1,:), combo_len, nbse);
            Hx_tens = reshape(H_all_tens(:,2:end,:), combo_len, ndep, nbse);
            fmm = obj.fmap_m;
            fmb = obj.fmap_b;

            ptx_raw = pflat(1:(ndep+1));
            ptx = fmm.*ptx_raw + fmb;
            pxp = pflat((ndep+2):end);
            [order_mat, poly_coeffs, icoeff_mat] = deal(obj.order_mat, obj.poly_coeffs, obj.icoeff_mat);

            txP = ptx.^((0:bor)');
            Ltx = comp_Lt_Lx(combo_len,bor,txP,poly_coeffs,order_mat);

            De_tx_H = nan(nvar,nbse);
            for ibse = 1:nbse
                De_tx_H(:,ibse) = (H_all_tens(:,:,ibse)')*prod(Ltx,2); %% compute indep. and dep. var terms
            end

            if (~prod(size(pxp)))
                De_out_H = De_tx_H;
            else
                dtxc_init = zeros(1,ndep+1);
                De_xp_dxde_mat = zeros(eor,combo_len);
                coeffs_init = ones(combo_len,1);
                Ptx_init = order_mat;
                Otx_init = order_mat;
                Ltxh_init = Ltx;
                Pxp_init = zeros(1,eor*ndep);

                specs = struct( 'eor',eor, ...
                                'ndep',ndep, ...
                                'bor', bor, ...
                                'p',pflat, ...
                                'txP',txP, ...
                                'fmm',fmm, ...
                                'pc',poly_coeffs, ...
                                'ic',icoeff_mat);

                De_xp_dxde_mat = compute_dxde_scalar(specs,1,dtxc_init,De_xp_dxde_mat,coeffs_init,Otx_init,Ptx_init,Ltxh_init,Pxp_init);

                De_xp_dxde_H = nan(eor*ndep,nbse);
                for ibse = 1:nbse
                    De_xp_dxde_H(:,ibse) = reshape((De_xp_dxde_mat*Hx_tens(:,:,ibse))',[],1);
                end

                De_xp_dtde_mat = zeros(eor*ndep,combo_len);

                bor_p1 = bor + 1;

                nDzi = obj.not_D_zero_inds;
                idt = nDzi(:,1);
                if (sum(idt))
                    dtc_p1 = dtxc_init(1) + 1;
                    idt_p1 = dtc_p1 + 1;
                    idt_front = 1:(bor_p1-dtc_p1);
                    idt_back = idt_p1:bor_p1;

                    Otx_dt = Otx_init(idt,:);
                    c_dt = (-1.0*coeffs_init(idt))*fmm(1);

                    Ptx_dt = [Ptx_init(idt,1)-1, Ptx_init(idt,2:end)];

                    dtLtxh = Ltxh_init(idt,:);
                    dtLtxh(:,1) = comp_dnL(dtc_p1,txP(idt_front,1),poly_coeffs(idt_back,idt_back),icoeff_mat(idt_back,idt_p1),Otx_dt(:,1));

                    dt_dtxc = dtxc_init;
                    dt_dtxc(1) = dtc_p1;

                    for kder = 1:eor
                        kder_inds = ((ndep*(kder-1))+1):(ndep*kder);
                        De_xp_dtde_mat(kder_inds,idt) = De_xp_dtde_mat(kder_inds,idt) + (pxp(kder_inds)')*((c_dt.*prod(dtLtxh,2))');
                        De_xp_dtde_mat(:,idt) = compute_dxde_vector(specs,kder+1,kder,dt_dtxc,De_xp_dtde_mat(:,idt),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_init);
                    end
                end

                for ix = 1:ndep
                    idxi = nDzi(:,ix+1);
                    if (sum(idxi))
                        dxic_p1 = dtxc_init(ix+1) + 1; % increase the xi derivative count by one
                        idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                        idxi_front = 1:(bor_p1-dxic_p1);
                        idxi_back = idxi_p1:bor_p1;

                        Otx_dxi = Otx_init(idxi,:);
                        c_dxi = (-1.0*coeffs_init(idxi))*fmm(ix+1); % update coefficients

                        Ptx_dxi = Ptx_init(idxi,:); % initialize power records to selected previous values
                        Ptx_dxi(:,ix+1) = Ptx_dxi(:,ix+1)-1; % initialize powers of xi

                        dxiLtxh = Ltxh_init(idxi,:); % initialize L functional history to selected previous values
                        % update L functional history to reflect derivative in t
                        dxiLtxh(:,ix+1) = comp_dnL(dxic_p1,txP(idxi_front,ix+1),poly_coeffs(idxi_back,idxi_back),icoeff_mat(idxi_back,idxi_p1),Otx_dxi(:,ix+1));

                        dxi_dtxc = dtxc_init;
                        dxi_dtxc(ix+1) = dxic_p1;

                        Pxp_dxi = Pxp_init; % initialize xp power records to previous values
                        Pxp_dxi(ix) = Pxp_init(ix) + 1; % update xp power records to reflect derivative in xi
                        for kder = 1:eor
                            kder_inds = ((ndep*(kder-1))+1):(ndep*kder);
                            De_xp_dtde_mat(kder_inds,idxi) = De_xp_dtde_mat(kder_inds,idxi) + ((pxp(kder_inds)')*((c_dxi.*prod(dxiLtxh,2))'))*pxp(ix); % add contribution of freshly computed t derivative term
                            De_xp_dtde_mat(:,idxi) = compute_dxde_vector(specs,kder+1,kder,dxi_dtxc,De_xp_dtde_mat(:,idxi),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
                        end
                    end
                end

                De_xp_dtde_H = De_xp_dtde_mat*ht_mat;

                De_xp_H = De_xp_dxde_H + De_xp_dtde_H;

                De_out_H = [De_tx_H; De_xp_H];
            end
        end
        function De_out = dp_de(obj,p_,h_,kcap_)
            if (nargin == 3)
                kcap = obj.eor;
            else
                kcap = kcap_;
            end
            pflat = reshape(p_,1,[]);
            htall = reshape(h_,[],1); %% ensuring we have the right dimensions
            [eor ndep ndof combo_len bor] = deal(kcap, obj.ndep, obj.ndof, obj.combo_len, double(obj.bor));

            H_all = reshape(htall,combo_len,1+ndep);
            ht = htall(1:combo_len);
            Hx = reshape(htall((combo_len+1):end), combo_len, ndep);
            fmm = obj.fmap_m;
            fmb = obj.fmap_b;

            ptx_raw = pflat(1:(ndep+1));
            ptx = fmm.*ptx_raw + fmb;
            pxp = pflat((ndep+2):end);
            [order_mat, poly_coeffs, icoeff_mat] = deal(obj.order_mat, obj.poly_coeffs, obj.icoeff_mat);

            txP = ptx.^((0:bor)');
            Ltx = comp_Lt_Lx(combo_len,bor,txP,poly_coeffs,order_mat);
            De_tx = (H_all')*prod(Ltx,2); %% compute indep. and dep. var terms

            dtxc_init = zeros(1,ndep+1);
            De_xp_dxde_mat = zeros(eor,combo_len);
            coeffs_init = ones(combo_len,1);
            Ptx_init = order_mat;
            Otx_init = order_mat;
            Ltxh_init = Ltx;
            Pxp_init = zeros(1,eor*ndep);

            specs = struct( 'eor',eor, ...
                            'ndep',ndep, ...
                            'bor', bor, ...
                            'p',pflat, ...
                            'txP',txP, ...
                            'fmm',fmm, ...
                            'pc',poly_coeffs, ...
                            'ic',icoeff_mat);

            De_xp_dxde_mat = compute_dxde_scalar(specs,1,dtxc_init,De_xp_dxde_mat,coeffs_init,Otx_init,Ptx_init,Ltxh_init,Pxp_init);
            De_xp_dxde = reshape((De_xp_dxde_mat*Hx)',[],1);

            De_xp_dtde_mat = zeros(eor*ndep,combo_len);

            bor_p1 = bor + 1;

            nDzi = obj.not_D_zero_inds;
            idt = nDzi(:,1);
            if (sum(idt))
                dtc_p1 = dtxc_init(1) + 1;
                idt_p1 = dtc_p1 + 1;
                idt_front = 1:(bor_p1-dtc_p1);
                idt_back = idt_p1:bor_p1;

                Otx_dt = Otx_init(idt,:);
                c_dt = (-1.0*coeffs_init(idt))*fmm(1);

                Ptx_dt = [Ptx_init(idt,1)-1, Ptx_init(idt,2:end)];

                dtLtxh = Ltxh_init(idt,:);
                dtLtxh(:,1) = comp_dnL(dtc_p1,txP(idt_front,1),poly_coeffs(idt_back,idt_back),icoeff_mat(idt_back,idt_p1),Otx_dt(:,1));

                dt_dtxc = dtxc_init;
                dt_dtxc(1) = dtc_p1;

                for kder = 1:eor
                    kder_inds = ((ndep*(kder-1))+1):(ndep*kder);
                    De_xp_dtde_mat(kder_inds,idt) = De_xp_dtde_mat(kder_inds,idt) + (pxp(kder_inds)')*((c_dt.*prod(dtLtxh,2))');
                    De_xp_dtde_mat(:,idt) = compute_dxde_vector(specs,kder+1,kder,dt_dtxc,De_xp_dtde_mat(:,idt),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_init);
                end
            end

            for ix = 1:ndep
                idxi = nDzi(:,ix+1);
                if (sum(idxi))
                    dxic_p1 = dtxc_init(ix+1) + 1; % increase the xi derivative count by one
                    idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                    idxi_front = 1:(bor_p1-dxic_p1);
                    idxi_back = idxi_p1:bor_p1;

                    Otx_dxi = Otx_init(idxi,:);
                    c_dxi = (-1.0*coeffs_init(idxi))*fmm(ix+1); % update coefficients

                    Ptx_dxi = Ptx_init(idxi,:); % initialize power records to selected previous values
                    Ptx_dxi(:,ix+1) = Ptx_dxi(:,ix+1)-1; % initialize powers of xi

                    dxiLtxh = Ltxh_init(idxi,:); % initialize L functional history to selected previous values
                    % update L functional history to reflect derivative in t
                    dxiLtxh(:,ix+1) = comp_dnL(dxic_p1,txP(idxi_front,ix+1),poly_coeffs(idxi_back,idxi_back),icoeff_mat(idxi_back,idxi_p1),Otx_dxi(:,ix+1));

                    dxi_dtxc = dtxc_init;
                    dxi_dtxc(ix+1) = dxic_p1;

                    Pxp_dxi = Pxp_init; % initialize xp power records to previous values
                    Pxp_dxi(ix) = Pxp_init(ix) + 1; % update xp power records to reflect derivative in xi
                    for kder = 1:eor
                        kder_inds = ((ndep*(kder-1))+1):(ndep*kder);
                        De_xp_dtde_mat(kder_inds,idxi) = De_xp_dtde_mat(kder_inds,idxi) + ((pxp(kder_inds)')*((c_dxi.*prod(dxiLtxh,2))'))*pxp(ix); % add contribution of freshly computed t derivative term
                        De_xp_dtde_mat(:,idxi) = compute_dxde_vector(specs,kder+1,kder,dxi_dtxc,De_xp_dtde_mat(:,idxi),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
                    end
                end
            end

            De_out = [De_tx; (De_xp_dtde_mat*ht) + De_xp_dxde];
            % De_out = De_out/norm(De_out);

            % De_xp_dtde = De_xp_dtde_mat*ht;
            % De_xp = De_xp_dxde + De_xp_dtde;
            % De_out = [De_tx; De_xp];
            % De_out = De_out/norm(De_out);
        end
        function debugging_description(obj)
            fprintf('%s function space (basis id: %d), order %d, %d degrees of freedom\n', ...
            obj.name, ...
            obj.bid, ...
            obj.bor, ...
            obj.ndof);

            fprintf('\nt  powers: ');
            for i = 1:obj.combo_len
                fprintf('%d ', obj.order_mat(i,1));
            end
            for l = 1:obj.ndep
                fprintf('\nx%d powers: ',l);
                for i = 1:obj.combo_len
                    fprintf('%d ', obj.order_mat(i,l+1));
                end
                fprintf('\n');
            end


            fprintf('orthogonal domain: (%.2e, %.2e)', orthopolynomial_basis.h_range(1), orthopolynomial_basis.h_range(2));

            % fprintf('\n%.2e < t < %.2e,', obj.sve(1,1), obj.sve(2,1));
            % for ix = 1:(obj.ndep)
            %     fprintf('\n%.2e < x%d < %.2e', obj.sve(1,ix+1), ix, obj.sve(2,ix+1));
            % end

            fprintf('\n\nforward mappings: t_hat = %.2e t + %.2e,', obj.fmap_m(1), obj.fmap_b(1));
            for ix = 1:(obj.ndep)
                fprintf('x%d_hat = %.2e x%d + %.2e\n', ix, obj.fmap_m(ix+1), ix, obj.fmap_b(ix+1));
            end
            fprintf('\ninverse mappings: t = %.2e t_hat + %.2e,', obj.bmap_m(1), obj.bmap_b(1));
            for ix = 1:(obj.ndep)
                fprintf('x%d = %.2e x%d_hat + %.2e\n', ix, obj.bmap_m(ix+1), ix, obj.bmap_b(ix+1));
            end

            fprintf('\nt  powers: ');
            for i = 1:obj.combo_len
                fprintf('%d ', obj.order_mat(i,1));
            end
            for l = 1:obj.ndep
                fprintf('\nx%d powers: ',l);
                for i = 1:obj.combo_len
                    fprintf('%d ', obj.order_mat(i,l+1));
                end
                fprintf('\n');
            end
            fprintf('\n%s polynomials:\n', obj.name);
            for i = 0:obj.bor
                fprintf('L%d: ',i);
                for j = 0:i
                    fprintf(' + (%.2e)u^%d', obj.poly_coeffs(i+1,j+1), j);
                end
                fprintf('\n');
            end
        end
    end
end

function [eor_ ndep_ bor_ p_ txP_ fmm_ pc_ ic_] = unpack_struct(str)
    [eor_ ndep_ bor_ p_ txP_ fmm_ pc_ ic_] = deal(str.eor, str.ndep, str.bor, str.p, str.txP, str.fmm, str.pc, str.ic);
end

function [De_xp_out Jac_out] = compute_dxde_vector_Jac(specs_,k_,source_,dtxc_,De_xp_in_,Jac_,H_,c_,Otx_,Ptx_,Ltxh_,Pxp_)
    Jac_v_k = zeros(size(Jac_));
    [eor_ ndep_ bor_ p_ txP_ fmm_ pc_ ic_] = unpack_struct(specs_);
    [bor_p1 nvar] = deal(bor_+1,ndep_+1);
    [tx_ xp_] = deal(p_(1:nvar), p_((nvar+1):end));

    iJ_k = ((k_-1)*ndep_+2):((k_*ndep_)+1);

    kder_inds = ((ndep_*(k_-1))+1):(ndep_*k_);
    source_inds = ((ndep_*(source_-1))+1):(ndep_*source_);

    id_all = ~(Ptx_==0);
    idt = id_all(:,1);
    if (k_==(specs_.eor+1))
        De_xp_out = De_xp_in_;
        if (sum(idt))
            dtc_p1 = dtxc_(1) + 1; % increase the t derivative count by one
            idt_p1 = dtc_p1 + 1; % index of relevant indices for derivative evaluation
            idt_front = 1:(bor_p1-dtc_p1);
            idt_back = idt_p1:bor_p1;

            Otx_dt = Otx_(idt,:);
            c_dt = c_(idt)*fmm_(1); % update coefficients

            Ptx_dt = [Ptx_(idt,1)-1, Ptx_(idt,2:end)]; % update powers

            dtLtxh = Ltxh_(idt,:); % initialize L functional history to selected previous terms
            % update L functional history to reflect derivative in t
            dtLtxh(:,1) = comp_dnL(dtc_p1,txP_(idt_front,1),pc_(idt_back,idt_back),ic_(idt_back,idt_p1),Otx_dt(:,1));

            dt_dtxc = dtxc_;
            dt_dtxc(1) = dtc_p1;

            Jac_dkm1_t = (xp_(source_inds)')*(c_dt.*prod(dtLtxh,2))';
            Jac_v_k(iJ_k,1,:) = pagemtimes(Jac_dkm1_t*prod(xp_.^Pxp_,2), H_(idt,1,:));
        end

        for ix = 2:nvar
            idxi = id_all(:,ix);
            if (sum(idxi))
                dxic_p1 = dtxc_(ix) + 1; % increase the xi derivative count by one
                idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                idxi_front = 1:(bor_p1-dxic_p1);
                idxi_back = idxi_p1:bor_p1;

                Otx_dxi = Otx_(idxi,:);
                c_dxi = c_(idxi)*fmm_(ix); % update coefficients

                Ptx_dxi = Ptx_(idxi,:); % initialize power records to selected previous values
                Ptx_dxi(:,ix) = Ptx_dxi(:,ix)-1; % initialize powers of xi

                dxiLtxh = Ltxh_(idxi,:); % initialize L functional history to selected previous values
                % update L functional history to reflect derivative in t
                dxiLtxh(:,ix) = comp_dnL(dxic_p1,txP_(idxi_front,ix),pc_(idxi_back,idxi_back),ic_(idxi_back,idxi_p1),Otx_dxi(:,ix));

                dxi_dtxc = dtxc_;
                dxi_dtxc(ix) = dxic_p1;

                Jac_dkm1_xi = (xp_(source_inds)')*((c_dxi.*prod(dxiLtxh,2))');

                Jac_v_k(iJ_k,ix,:) = pagemtimes(Jac_dkm1_xi*prod(xp_.^Pxp_,2),H_(idxi,1,:));
            end
        end

        idim = 0;
        idim_full = 1 + ndep_;
        Ltxh_eval = prod(Ltxh_,2);
        for kder = 1:(k_-1)
            for idep = 1:ndep_
                idim = idim + 1;
                idim_full = idim_full + 1;
                if (Pxp_(idim))
                    c_dxpi = c_*Pxp_(idim);
                    Pxp_dxpi = Pxp_;

                    Jac_dkm1_dxi = (xp_(source_inds)')*((c_dxpi.*Ltxh_eval)');

                    Pxp_dxpi(idim) = Pxp_dxpi(idim)-1;
                    Jac_v_k(iJ_k,idim_full,:) = pagemtimes(Jac_dkm1_dxi*prod(xp_.^Pxp_dxpi),H_(:,1,:));
                end
            end
        end
        Jac_dkm1_dxi_vec = ((c_.*Ltxh_eval)')*prod(xp_.^Pxp_,2);
        Jac_v_k(iJ_k,source_inds+ndep_+1,:) = Jac_v_k(iJ_k,source_inds+ndep_+1,:) + pagemtimes(Jac_dkm1_dxi_vec,H_(:,1,:));

        Jac_out = Jac_ + Jac_v_k;
    else
        De_xp_k = zeros(size(De_xp_in_));
        if (sum(idt))
            dtc_p1 = dtxc_(1) + 1; % increase the t derivative count by one
            idt_p1 = dtc_p1 + 1; % index of relevant indices for derivative evaluation
            idt_front = 1:(bor_p1-dtc_p1);
            idt_back = idt_p1:bor_p1;

            Otx_dt = Otx_(idt,:);
            c_dt = c_(idt)*fmm_(1); % update coefficients

            Ptx_dt = [Ptx_(idt,1)-1, Ptx_(idt,2:end)]; % update powers

            dtLtxh = Ltxh_(idt,:); % initialize L functional history to selected previous terms
            % update L functional history to reflect derivative in t
            dtLtxh(:,1) = comp_dnL(dtc_p1,txP_(idt_front,1),pc_(idt_back,idt_back),ic_(idt_back,idt_p1),Otx_dt(:,1));

            dt_dtxc = dtxc_;
            dt_dtxc(1) = dtc_p1;

            Jac_dkm1_t = (xp_(source_inds)')*(c_dt.*prod(dtLtxh,2))';

            De_xp_k(kder_inds,idt) = Jac_dkm1_t*prod(xp_.^Pxp_,2); % add contribution of freshly computed t derivative term
            Jac_v_k(iJ_k,1,:) = pagemtimes(Jac_dkm1_t*prod(xp_.^Pxp_,2), H_(idt,1,:));

            [De_xp_k(:,idt) Jac_v_k] = compute_dxde_vector_Jac(specs_,k_+1,source_,dt_dtxc,De_xp_k(:,idt),Jac_v_k,H_(idt,:,:),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_);
        end

        for ix = 2:nvar
            idxi = id_all(:,ix);
            if (sum(idxi))
                dxic_p1 = dtxc_(ix) + 1; % increase the xi derivative count by one
                idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                idxi_front = 1:(bor_p1-dxic_p1);
                idxi_back = idxi_p1:bor_p1;

                Otx_dxi = Otx_(idxi,:);
                c_dxi = c_(idxi)*fmm_(ix); % update coefficients

                Ptx_dxi = Ptx_(idxi,:); % initialize power records to selected previous values
                Ptx_dxi(:,ix) = Ptx_dxi(:,ix)-1; % initialize powers of xi

                dxiLtxh = Ltxh_(idxi,:); % initialize L functional history to selected previous values
                % update L functional history to reflect derivative in t
                dxiLtxh(:,ix) = comp_dnL(dxic_p1,txP_(idxi_front,ix),pc_(idxi_back,idxi_back),ic_(idxi_back,idxi_p1),Otx_dxi(:,ix));

                dxi_dtxc = dtxc_;
                dxi_dtxc(ix) = dxic_p1;

                Jac_dkm1_xi = (xp_(source_inds)')*((c_dxi.*prod(dxiLtxh,2))');

                Pxp_dxi = Pxp_; % initialize xp power records to previous values
                Pxp_dxi(ix-1) = Pxp_dxi(ix-1) + 1; % update xp power records to reflect derivative in xi

                De_xp_k(kder_inds,idxi) = De_xp_k(kder_inds,idxi) + Jac_dkm1_xi*prod(xp_.^Pxp_dxi,2); % add contribution of freshly computed t derivative term
                Jac_v_k(iJ_k,ix,:) = pagemtimes(Jac_dkm1_xi*prod(xp_.^Pxp_,2),H_(idxi,1,:));

                [De_xp_k(:,idxi) Jac_v_k] = compute_dxde_vector_Jac(specs_,k_+1,source_,dxi_dtxc,De_xp_k(:,idxi),Jac_v_k,H_(idxi,:,:),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
            end
        end

        idim = 0;
        idim_full = 1 + ndep_;
        Ltxh_eval = prod(Ltxh_,2);
        for kder = 1:(k_-1)
            for idep = 1:ndep_
                idim = idim + 1;
                idim_full = idim_full + 1;
                if (Pxp_(idim))
                    c_dxpi = c_*Pxp_(idim);
                    Pxp_dxpi = Pxp_;

                    Jac_dkm1_dxi = (xp_(source_inds)')*((c_dxpi.*Ltxh_eval)');

                    Pxp_dxpi(idim) = Pxp_dxpi(idim)-1;
                    Jac_v_k(iJ_k,idim_full,:) = pagemtimes(Jac_dkm1_dxi*prod(xp_.^Pxp_dxpi),H_(:,1,:));

                    Pxp_dxpi(idim+ndep_) = Pxp_dxpi(idim+ndep_) + 1;
                    De_xp_k(kder_inds,:) = De_xp_k(kder_inds,:) + Jac_dkm1_dxi*prod(xp_.^Pxp_dxpi,2);

                    [De_xp_k Jac_v_k] = compute_dxde_vector_Jac(specs_,k_+1,source_,dtxc_,De_xp_k,Jac_v_k,H_,c_dxpi,Otx_,Ptx_,Ltxh_,Pxp_dxpi);
                end
            end
        end

        Jac_dkm1_dxi_vec = ((c_.*Ltxh_eval)')*prod(xp_.^Pxp_,2);

        De_xp_k(kder_inds,:) = De_xp_k(kder_inds,:) + (xp_(source_inds+ndep_)')*Jac_dkm1_dxi_vec;
        Jac_v_k(iJ_k,source_inds+ndep_+1,:) = Jac_v_k(iJ_k,source_inds+ndep_+1,:) + pagemtimes(Jac_dkm1_dxi_vec,H_(:,1,:));

        [De_xp_k Jac_v_k] = compute_dxde_vector_Jac(specs_,k_+1,k_,dtxc_,De_xp_k,Jac_v_k,H_,c_,Otx_,Ptx_,Ltxh_,Pxp_);

        De_xp_out = De_xp_in_ + De_xp_k;
        Jac_out = Jac_ + Jac_v_k;
    end
end

function [De_xp_out Jac_out] = compute_dxde_scalar_Jac(specs_,k_,dtxc_,De_xp_in_,Jac_,H_,c_,Otx_,Ptx_,Ltxh_,Pxp_)
    Jac_v_k = zeros(size(Jac_));
    [eor_ ndep_ bor_ p_ txP_ fmm_ pc_ ic_] = unpack_struct(specs_);
    [bor_p1 nvar] = deal(bor_+1,ndep_+1);
    [tx_ xp_] = deal(p_(1:nvar), p_((nvar+1):end));
    iJ_k = ((k_-1)*ndep_+2):((k_*ndep_)+1);

    id_all = ~(Ptx_==0);
    idt = id_all(:,1);
    if (k_==(specs_.eor + 1))
        De_xp_out = De_xp_in_;
        if (sum(idt)) %% if there are derivatives in t we may determine the Jacobian of the previous order
            dtc_p1 = dtxc_(1) + 1; % increase the t derivative count by one
            idt_p1 = dtc_p1 + 1; % index of relevant indices for derivative evaluation
            idt_front = 1:(bor_p1-dtc_p1);
            idt_back = idt_p1:bor_p1;

            Otx_dt = Otx_(idt,:);
            c_dt = c_(idt)*fmm_(1); % update coefficients

            Ptx_dt = [Ptx_(idt,1)-1, Ptx_(idt,2:end)]; % update powers

            dtLtxh = Ltxh_(idt,:); % initialize L functional history to selected previous terms
            % update L functional history to reflect derivative in t
            dtLtxh(:,1) = comp_dnL(dtc_p1,txP_(idt_front,1),pc_(idt_back,idt_back),ic_(idt_back,idt_p1),Otx_dt(:,1));

            dt_dtxc = dtxc_;
            dt_dtxc(1) = dtc_p1;

            Jac_dkm1_t = (c_dt.*prod(dtLtxh,2))';
            Jac_v_k(iJ_k,1,:) = pagemtimes(Jac_dkm1_t*prod(xp_.^Pxp_,2),H_(idt,2:end,:));
        end

        for ix = 2:nvar
            idxi = id_all(:,ix);
            if (sum(idxi))
                dxic_p1 = dtxc_(ix) + 1; % increase the xi derivative count by one
                idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                idxi_front = 1:(bor_p1-dxic_p1);
                idxi_back = idxi_p1:bor_p1;

                Otx_dxi = Otx_(idxi,:);
                c_dxi = c_(idxi)*fmm_(ix); % update coefficients

                Ptx_dxi = Ptx_(idxi,:); % initialize power records to selected previous values
                Ptx_dxi(:,ix) = Ptx_dxi(:,ix)-1; % initialize powers of xi

                dxiLtxh = Ltxh_(idxi,:); % initialize L functional history to selected previous values
                % update L functional history to reflect derivative in t
                dxiLtxh(:,ix) = comp_dnL(dxic_p1,txP_(idxi_front,ix),pc_(idxi_back,idxi_back),ic_(idxi_back,idxi_p1),Otx_dxi(:,ix));

                dxi_dtxc = dtxc_;
                dxi_dtxc(ix) = dxic_p1;

                Jac_dkm1_xi = (c_dxi.*prod(dxiLtxh,2))';
                Jac_v_k(iJ_k,ix,:) = pagemtimes(Jac_dkm1_xi*prod(xp_.^Pxp_,2),H_(idxi,2:end,:));
            end
        end

        idim = 0;
        idim_full = 1 + ndep_;
        Ltxh_eval = prod(Ltxh_,2);
        for kder = 1:(k_-1)
            for idep = 1:ndep_
                idim = idim + 1;
                idim_full = idim_full + 1;
                if (Pxp_(idim))
                    c_dxpi = c_*Pxp_(idim);
                    Pxp_dxpi = Pxp_;

                    Jac_dkm1_dxi = ((c_dxpi.*Ltxh_eval)');

                    Pxp_dxpi(idim) = Pxp_dxpi(idim)-1;
                    Jac_v_k(iJ_k,idim_full,:) = pagemtimes(Jac_dkm1_dxi*prod(xp_.^Pxp_dxpi),H_(:,2:end,:));
                end
            end
        end
        Jac_out = Jac_ + Jac_v_k;
    else
        De_xp_k = zeros(size(De_xp_in_));
        if (sum(idt)) %% if there are derivatives in t we may determine the Jacobian of the previous order
            dtc_p1 = dtxc_(1) + 1; % increase the t derivative count by one
            idt_p1 = dtc_p1 + 1; % index of relevant indices for derivative evaluation
            idt_front = 1:(bor_p1-dtc_p1);
            idt_back = idt_p1:bor_p1;

            Otx_dt = Otx_(idt,:);
            c_dt = c_(idt)*fmm_(1); % update coefficients

            Ptx_dt = [Ptx_(idt,1)-1, Ptx_(idt,2:end)]; % update powers

            dtLtxh = Ltxh_(idt,:); % initialize L functional history to selected previous terms
            % update L functional history to reflect derivative in t
            dtLtxh(:,1) = comp_dnL(dtc_p1,txP_(idt_front,1),pc_(idt_back,idt_back),ic_(idt_back,idt_p1),Otx_dt(:,1));

            dt_dtxc = dtxc_;
            dt_dtxc(1) = dtc_p1;

            Jac_dkm1_t = (c_dt.*prod(dtLtxh,2))';

            De_xp_k(k_,idt) = Jac_dkm1_t*prod(xp_.^Pxp_,2); % add contribution of freshly computed t derivative term
            Jac_v_k(iJ_k,1,:) = pagemtimes(Jac_dkm1_t*prod(xp_.^Pxp_,2),H_(idt,2:end,:));

            [De_xp_k(:,idt) Jac_v_k] = compute_dxde_scalar_Jac(specs_,k_+1,dt_dtxc,De_xp_k(:,idt),Jac_v_k,H_(idt,:,:),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_);
        end

        for ix = 2:nvar
            idxi = id_all(:,ix);
            if (sum(idxi))
                dxic_p1 = dtxc_(ix) + 1; % increase the xi derivative count by one
                idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                idxi_front = 1:(bor_p1-dxic_p1);
                idxi_back = idxi_p1:bor_p1;

                Otx_dxi = Otx_(idxi,:);
                c_dxi = c_(idxi)*fmm_(ix); % update coefficients

                Ptx_dxi = Ptx_(idxi,:); % initialize power records to selected previous values
                Ptx_dxi(:,ix) = Ptx_dxi(:,ix)-1; % initialize powers of xi

                dxiLtxh = Ltxh_(idxi,:); % initialize L functional history to selected previous values
                % update L functional history to reflect derivative in t
                dxiLtxh(:,ix) = comp_dnL(dxic_p1,txP_(idxi_front,ix),pc_(idxi_back,idxi_back),ic_(idxi_back,idxi_p1),Otx_dxi(:,ix));

                dxi_dtxc = dtxc_;
                dxi_dtxc(ix) = dxic_p1;

                Jac_dkm1_xi = (c_dxi.*prod(dxiLtxh,2))';

                Pxp_dxi = Pxp_; % initialize xp power records to previous values
                Pxp_dxi(ix-1) = Pxp_dxi(ix-1) + 1; % update xp power records to reflect derivative in xi

                De_xp_k(k_,idxi) = De_xp_k(k_,idxi) + Jac_dkm1_xi*prod(xp_.^Pxp_dxi,2); % add contribution of freshly computed t derivative term
                Jac_v_k(iJ_k,ix,:) = pagemtimes(Jac_dkm1_xi*prod(xp_.^Pxp_,2),H_(idxi,2:end,:));

                [De_xp_k(:,idxi) Jac_v_k] = compute_dxde_scalar_Jac(specs_,k_+1,dxi_dtxc,De_xp_k(:,idxi),Jac_v_k,H_(idxi,:,:),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
            end
        end

        idim = 0;
        idim_full = 1 + ndep_;
        Ltxh_eval = prod(Ltxh_,2);
        for kder = 1:(k_-1)
            for idep = 1:ndep_
                idim = idim + 1;
                idim_full = idim_full + 1;
                if (Pxp_(idim))
                    c_dxpi = c_*Pxp_(idim);
                    Pxp_dxpi = Pxp_;

                    Jac_dkm1_dxi = ((c_dxpi.*Ltxh_eval)');

                    Pxp_dxpi(idim) = Pxp_dxpi(idim)-1;
                    Jac_v_k(iJ_k,idim_full,:) = pagemtimes(Jac_dkm1_dxi*prod(xp_.^Pxp_dxpi),H_(:,2:end,:));

                    Pxp_dxpi(idim+ndep_) = Pxp_dxpi(idim+ndep_) + 1;
                    De_xp_k(k_,:) = De_xp_k(k_,:) + Jac_dkm1_dxi*prod(xp_.^Pxp_dxpi,2);

                    [De_xp_k Jac_v_k] = compute_dxde_scalar_Jac(specs_,k_+1,dtxc_,De_xp_k,Jac_v_k,H_,c_dxpi,Otx_,Ptx_,Ltxh_,Pxp_dxpi);
                end
            end
        end
        De_xp_out = De_xp_in_ + De_xp_k;
        Jac_out = Jac_ + Jac_v_k;
    end
end

function De_xp_out = compute_dxde_vector(specs_,k_,source_,dtxc_,De_xp_in_,c_,Otx_,Ptx_,Ltxh_,Pxp_)
    if (k_>(specs_.eor))
        De_xp_out = De_xp_in_;
    else
        De_xp_k = zeros(size(De_xp_in_));

        [eor_ ndep_ bor_ p_ txP_ fmm_ pc_ ic_] = unpack_struct(specs_);
        [bor_p1 nvar] = deal(bor_+1,ndep_+1);
        [tx_ xp_] = deal(p_(1:nvar), p_((nvar+1):end));

        kder_inds = ((ndep_*(k_-1))+1):(ndep_*k_);
        source_inds = ((ndep_*(source_-1))+1):(ndep_*source_);

        id_all = ~(Ptx_==0);
        idt = id_all(:,1);
        if (sum(idt))
            dtc_p1 = dtxc_(1) + 1; % increase the t derivative count by one
            idt_p1 = dtc_p1 + 1; % index of relevant indices for derivative evaluation
            idt_front = 1:(bor_p1-dtc_p1);
            idt_back = idt_p1:bor_p1;

            Otx_dt = Otx_(idt,:);
            c_dt = c_(idt)*fmm_(1); % update coefficients

            Ptx_dt = [Ptx_(idt,1)-1, Ptx_(idt,2:end)]; % update powers

            dtLtxh = Ltxh_(idt,:); % initialize L functional history to selected previous terms
            % update L functional history to reflect derivative in t
            dtLtxh(:,1) = comp_dnL(dtc_p1,txP_(idt_front,1),pc_(idt_back,idt_back),ic_(idt_back,idt_p1),Otx_dt(:,1));

            dt_dtxc = dtxc_;
            dt_dtxc(1) = dtc_p1;

            De_xp_k(kder_inds,idt) = (xp_(source_inds)')*(c_dt.*prod(dtLtxh,2))'*prod(xp_.^Pxp_,2); % add contribution of freshly computed t derivative term
            De_xp_k(:,idt) = compute_dxde_vector(specs_,k_+1,source_,dt_dtxc,De_xp_k(:,idt),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_);
        end

        for ix = 2:nvar
            idxi = id_all(:,ix);
            if (sum(idxi))
                dxic_p1 = dtxc_(ix) + 1; % increase the xi derivative count by one
                idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                idxi_front = 1:(bor_p1-dxic_p1);
                idxi_back = idxi_p1:bor_p1;

                Otx_dxi = Otx_(idxi,:);
                c_dxi = c_(idxi)*fmm_(ix); % update coefficients

                Ptx_dxi = Ptx_(idxi,:); % initialize power records to selected previous values
                Ptx_dxi(:,ix) = Ptx_dxi(:,ix)-1; % initialize powers of xi

                dxiLtxh = Ltxh_(idxi,:); % initialize L functional history to selected previous values
                % update L functional history to reflect derivative in t
                dxiLtxh(:,ix) = comp_dnL(dxic_p1,txP_(idxi_front,ix),pc_(idxi_back,idxi_back),ic_(idxi_back,idxi_p1),Otx_dxi(:,ix));

                dxi_dtxc = dtxc_;
                dxi_dtxc(ix) = dxic_p1;

                Pxp_dxi = Pxp_; % initialize xp power records to previous values
                Pxp_dxi(ix-1) = Pxp_dxi(ix-1) + 1; % update xp power records to reflect derivative in xi

                De_xp_k(kder_inds,idxi) = De_xp_k(kder_inds,idxi) + (xp_(source_inds)')*((c_dxi.*prod(dxiLtxh,2))')*prod(xp_.^Pxp_dxi,2); % add contribution of freshly computed t derivative term
                De_xp_k(:,idxi) = compute_dxde_vector(specs_,k_+1,source_,dxi_dtxc,De_xp_k(:,idxi),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
            end
        end

        idim = 0;
        Ltxh_eval = prod(Ltxh_,2);
        for kder = 1:(k_-1)
            for idep = 1:ndep_
                idim = idim + 1;
                if (Pxp_(idim))
                    c_dxpi = c_*Pxp_(idim);
                    Pxp_dxpi = Pxp_;
                    Pxp_dxpi(idim) = Pxp_dxpi(idim)-1;
                    Pxp_dxpi(idim+ndep_) = Pxp_dxpi(idim+ndep_) + 1;

                    De_xp_k(kder_inds,:) = De_xp_k(kder_inds,:) + (xp_(source_inds)')*((c_dxpi.*Ltxh_eval)')*prod(xp_.^Pxp_dxpi,2);
                    De_xp_k = compute_dxde_vector(specs_,k_+1,source_,dtxc_,De_xp_k,c_dxpi,Otx_,Ptx_,Ltxh_,Pxp_dxpi);
                end
            end
        end

        De_xp_k(kder_inds,:) = De_xp_k(kder_inds,:) + (xp_(source_inds+ndep_)')*((c_.*Ltxh_eval)')*prod(xp_.^Pxp_,2);
        De_xp_k = compute_dxde_vector(specs_,k_+1,k_,dtxc_,De_xp_k,c_,Otx_,Ptx_,Ltxh_,Pxp_);
        De_xp_out = De_xp_in_ + De_xp_k;
    end
end

function De_xp_out = compute_dxde_scalar(specs_,k_,dtxc_,De_xp_in_,c_,Otx_,Ptx_,Ltxh_,Pxp_)
    if (k_>(specs_.eor))
        De_xp_out = De_xp_in_;
    else
        De_xp_k = zeros(size(De_xp_in_));
        [eor_ ndep_ bor_ p_ txP_ fmm_ pc_ ic_] = unpack_struct(specs_);
        [bor_p1 nvar] = deal(bor_+1,ndep_+1);
        [tx_ xp_] = deal(p_(1:nvar), p_((nvar+1):end));

        id_all = ~(Ptx_==0);
        idt = id_all(:,1);
        if (sum(idt))
            dtc_p1 = dtxc_(1) + 1; % increase the t derivative count by one
            idt_p1 = dtc_p1 + 1; % index of relevant indices for derivative evaluation
            idt_front = 1:(bor_p1-dtc_p1);
            idt_back = idt_p1:bor_p1;

            Otx_dt = Otx_(idt,:);
            c_dt = c_(idt)*fmm_(1); % update coefficients

            Ptx_dt = [Ptx_(idt,1)-1, Ptx_(idt,2:end)]; % update powers

            dtLtxh = Ltxh_(idt,:); % initialize L functional history to selected previous terms
            % update L functional history to reflect derivative in t
            dtLtxh(:,1) = comp_dnL(dtc_p1,txP_(idt_front,1),pc_(idt_back,idt_back),ic_(idt_back,idt_p1),Otx_dt(:,1));

            dt_dtxc = dtxc_;
            dt_dtxc(1) = dtc_p1;

            De_xp_k(k_,idt) = (c_dt.*prod(dtLtxh,2))'*prod(xp_.^Pxp_,2); % add contribution of freshly computed t derivative term
            De_xp_k(:,idt) = compute_dxde_scalar(specs_,k_+1,dt_dtxc,De_xp_k(:,idt),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_);
        end

        for ix = 2:nvar
            idxi = id_all(:,ix);
            if (sum(idxi))
                dxic_p1 = dtxc_(ix) + 1; % increase the xi derivative count by one
                idxi_p1 = dxic_p1 + 1; % index of relevant indices for derivative evaluation
                idxi_front = 1:(bor_p1-dxic_p1);
                idxi_back = idxi_p1:bor_p1;

                Otx_dxi = Otx_(idxi,:);
                c_dxi = c_(idxi)*fmm_(ix); % update coefficients

                Ptx_dxi = Ptx_(idxi,:); % initialize power records to selected previous values
                Ptx_dxi(:,ix) = Ptx_dxi(:,ix)-1; % initialize powers of xi

                dxiLtxh = Ltxh_(idxi,:); % initialize L functional history to selected previous values
                % update L functional history to reflect derivative in t
                dxiLtxh(:,ix) = comp_dnL(dxic_p1,txP_(idxi_front,ix),pc_(idxi_back,idxi_back),ic_(idxi_back,idxi_p1),Otx_dxi(:,ix));

                dxi_dtxc = dtxc_;
                dxi_dtxc(ix) = dxic_p1;

                Pxp_dxi = Pxp_; % initialize xp power records to previous values
                Pxp_dxi(ix-1) = Pxp_dxi(ix-1) + 1; % update xp power records to reflect derivative in xi

                De_xp_k(k_,idxi) = De_xp_k(k_,idxi) + (c_dxi.*prod(dxiLtxh,2))'*prod(xp_.^Pxp_dxi,2); % add contribution of freshly computed t derivative term
                De_xp_k(:,idxi) = compute_dxde_scalar(specs_,k_+1,dxi_dtxc,De_xp_k(:,idxi),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
            end
        end

        idim = 0;
        Ltxh_eval = prod(Ltxh_,2);
        for kder = 1:(k_-1)
            for idep = 1:ndep_
                idim = idim + 1;
                if (Pxp_(idim))
                    c_dxpi = c_*Pxp_(idim);
                    Pxp_dxpi = Pxp_;
                    Pxp_dxpi(idim) = Pxp_dxpi(idim)-1;
                    Pxp_dxpi(idim+ndep_) = Pxp_dxpi(idim+ndep_) + 1;

                    De_xp_k(k_,:) = De_xp_k(k_,:) + (((c_dxpi.*Ltxh_eval)')*prod(xp_.^Pxp_dxpi,2));
                    De_xp_k = compute_dxde_scalar(specs_,k_+1,dtxc_,De_xp_k,c_dxpi,Otx_,Ptx_,Ltxh_,Pxp_dxpi);
                end
            end
        end
        De_xp_out = De_xp_in_ + De_xp_k;
    end
end

function Ltx_out = comp_Lt_Lx(combo_len_, bor_, txP_, poly_coeffs_, poly_ords_)
    nvar = size(txP_,2);
    Ltx_out = nan(combo_len_,nvar);
    for i = 0:bor_
        ii = i+1;
        Lpq = poly_coeffs_(ii,1:ii)*txP_(1:(ii),:);
        ipq = poly_ords_==i;
        for j = 1:nvar
            Ltx_out(ipq(:,j),j) = Lpq(j);
        end
    end
end

function dL_out = comp_dnL(dc_,vp_,pc_,ic_,O_)
    dL_out = zeros(size(O_));
    ishift = 1;
    for i = dc_:max(O_)
        ii = i+1;

        dL_O = (pc_(ishift,1:ishift))*(ic_(1:ishift).*vp_(1:ishift));

        dL_out(O_==i) = dL_O;
        ishift = ishift+1;
    end
end
