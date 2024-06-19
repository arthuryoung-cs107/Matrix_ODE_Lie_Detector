classdef LD_orthopolynomial_space < LD_power_space
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
            fprintf('(LD_orthopolynomial_space::read_domain_configuration) reading %s\n', name_);
            file = fopen(name_);

            hlen = fread(file,1,'int=>int');
            header_rest = fread(file,hlen,'int=>int');
            [bor_in ord_mat_len map_len len_sym] = deal(header_rest(1), header_rest(2), header_rest(3), header_rest(4));
            ord_mat_flat = fread(file,ord_mat_len,'int=>int');
            map_params_flat = fread(file,map_len,'double=>double');
            poly_coeffs_flat = fread(file,len_sym,'double=>double');
            fclose(file);

            nvar = map_len/4;
            ndep = nvar - 1;
            ndof_full = ord_mat_len;
            perm_len = ndof_full/nvar;

            order_mat = double(reshape(ord_mat_flat,nvar,[]))';
            dim0_lens = zeros(bor_in+1,1);
            not_D_zero_inds = ~(order_mat==0);

            poly_coeffs = obj.poly_coeffs;
            [coeff_len,icoeff] = deal(1);
            for ipow = 0:bor_in
                dim0_lens(ipow+1) = sum(double(order_mat(:,1)==ipow));
                poly_coeffs(ipow+1,1:coeff_len) = poly_coeffs_flat(icoeff:(icoeff+coeff_len-1));
                icoeff = icoeff + coeff_len;
                coeff_len = coeff_len + 1;
            end

            mparams_mat = reshape(map_params_flat,nvar,4)';

            obj_out = obj;
            obj_out.bor = bor_in;
            obj_out.perm_len = perm_len;
            obj_out.ndep = ndep;
            obj_out.ndim = 1 + ndep*(obj.eor+1); % assume input order is the same, does not affect actual parameters
            obj_out.order_mat = order_mat;
            obj_out.dim0_lens = dim0_lens;
            obj_out.not_D_zero_inds = not_D_zero_inds;
            obj_out.poly_coeffs = poly_coeffs;

            obj_out.fmap_m = mparams_mat(1,:);
            obj_out.fmap_b = mparams_mat(2,:);
            obj_out.bmap_m = mparams_mat(3,:);
            obj_out.bmap_b = mparams_mat(4,:);
        end
        function write_domain_configuration(obj,name_)
            hlen = 4;
            nvar = 1+obj.ndep;
            bor = obj.bor;
            ord_mat_len = obj.perm_len*nvar;
            map_len = 4*nvar;
            borp1 = bor + 1;
            len_sym = (borp1*(borp1+1))/2;
            coeff_vec_out = zeros(len_sym,1);
            icoeff = 1;
            for i = 0:bor
                for j = 0:i
                    coeff_vec_out(icoeff) = obj.poly_coeffs(i+1,j+1);
                    icoeff = icoeff + 1;
                end
            end
            ord_vec = reshape((obj.order_mat)',[],1);
            mparams_mat = [obj.fmap_m'; obj.fmap_b'; obj.bmap_m'; obj.bmap_b'];
            mparams_vec = mparams_mat(:);

            header = [hlen,bor,ord_mat_len,map_len,len_sym]';
            file = fopen(name_,'w+');
            fwrite(file, header, 'int');
            fwrite(file, ord_vec, 'int');
            fwrite(file, mparams_vec, 'double');
            fwrite(file, coeff_vec_out, 'double');
            fclose(file);
            fprintf('(LD_orthopolynomial_space::write_ode_experiment_data) wrote file: %s \n', name_);
        end
        function v_mat_out = ds_de_ratio_ptsmat_stabilized_fast(obj,s_,params_)
            [~,params_mat] = ds_de_ratio_n1_stabilized_fast_itemized(obj,s_,params_);
            [ndim,npts] = size(s_);
            v_mat_out = nan(ndim,npts);
            for i = 1:npts
                v_mat_out(:,i) = obj.ds_de_thetamat_Jac(s_(:,i),params_mat(:,i));
            end
        end
        function v_out = ds_de_ratio_n1_stabilized_fast(obj,s_,params_)
            [ndep,nvar,npts,nparams,perm_len,bor,poly_coeffs,order_mat] = deal(obj.ndep,1+obj.ndep, ...
                size(s_,2),size(params_,2),obj.perm_len,double(obj.bor),obj.poly_coeffs,obj.order_mat);
            xu_raw_mat = s_(1:nvar,:);
            xu_mat = ((obj.fmap_m.*xu_raw_mat') + obj.fmap_b)';
            Theta_x_mat = (params_(1:perm_len,:))';
            Theta_u_mat_T = params_((perm_len+1):end,:);
            Lxu_inds = order_mat + [0, (bor+1)*(1:ndep)] + 1;
            v_u_out = nan(ndep,npts);
            for ipts = 1:npts
                xu_pow = poly_coeffs*((xu_mat(:,ipts).^(0:bor))');
                Lvec_i = prod(xu_pow(Lxu_inds),2);
                vx_i_vec = Theta_x_mat*Lvec_i;
                v_u_out(:,ipts) = ((Lvec_i')*reshape((Theta_u_mat_T*vx_i_vec)/(vx_i_vec'*vx_i_vec),perm_len,ndep))';
            end
            v_out = [ones(1,npts) ; v_u_out];
        end
        function [v_out,params_out] = ds_de_ratio_n1_stabilized_fast_itemized(obj,s_,params_)
            [ndep,nvar,npts,nparams,perm_len,bor,poly_coeffs,order_mat] = deal(obj.ndep,1+obj.ndep, ...
                size(s_,2),size(params_,2),obj.perm_len,double(obj.bor),obj.poly_coeffs,obj.order_mat);
            xu_raw_mat = s_(1:nvar,:);
            xu_mat = ((obj.fmap_m.*xu_raw_mat') + obj.fmap_b)';
            Theta_x_mat = (params_(1:perm_len,:))';
            Lxu_inds = order_mat + [0, (bor+1)*(1:ndep)] + 1;
            v_u_out = nan(ndep,npts);
            params_out = nan(obj.ndof_full,npts);
            for ipts = 1:npts
                xu_pow = poly_coeffs*((xu_mat(:,ipts).^(0:bor))');
                Lvec_i = prod(xu_pow(Lxu_inds),2);
                vx_i_vec = Theta_x_mat*Lvec_i;
                params_out(:,ipts) = (params_*vx_i_vec)/(vx_i_vec'*vx_i_vec);
                v_u_out(:,ipts) = ((Lvec_i')*reshape(params_out((perm_len+1):end,ipts),perm_len,ndep))';
            end
            v_out = [ones(1,npts) ; v_u_out];
        end
        function v_mat_out = ds_de_ptsmat(obj,s_,theta_)
            [ndim,npts] = size(s_);
            v_mat_out = nan(ndim,npts);
            for i = 1:npts
                v_mat_out(:,i) = obj.ds_de_thetamat_Jac(s_(:,i),theta_);
            end
        end
        function [v_out_theta_mat Jac_v_out Jac_L_xu_theta_out] = ds_de_thetamat_Jac(obj,s_,theta_mat_,kcap_)
            if (nargin == 3)
                kcap = obj.eor;
            else
                kcap = kcap_;
            end
            [eor ndep ndim ndof perm_len bor] = deal(   kcap, obj.ndep, 1+(obj.ndep*(kcap+1)), obj.ndof_full, ...
                                                        obj.perm_len, double(obj.bor));
            [order_mat, poly_coeffs, icoeff_mat] = deal(obj.order_mat, obj.poly_coeffs, obj.icoeff_mat);
            [fmm,fmb] = deal(obj.fmap_m,obj.fmap_b);
            nvar = ndep + 1;
            sflat = reshape(s_,1,[]);
            theta_mat_tall = reshape(theta_mat_,ndof,[]); %% ensuring we have the right dimensions
            nbse = size(theta_mat_tall,2); % number of solutions to check

            theta_all_tens = reshape(theta_mat_tall,perm_len,nvar,nbse);
            theta_x_mat = reshape(theta_all_tens(:,1,:), perm_len, nbse);
            theta_u_tens = reshape(theta_all_tens(:,2:end,:), perm_len, ndep, nbse);

            xu_raw = sflat(1:(nvar));
            xu = fmm.*xu_raw + fmb;
            pdxu = sflat((ndep+2):end);

            xuP = xu.^((0:bor)');
            Lxu = comp_Lx_Lu(perm_len,bor,xuP,poly_coeffs,order_mat);

            v_xu_theta = nan(nvar,nbse);
            for ibse = 1:nbse
                v_xu_theta(:,ibse) = (theta_all_tens(:,:,ibse)')*prod(Lxu,2); %% compute indep. and dep. var terms
            end

            dxuc_init = zeros(1,nvar);
            v_dxu_dude_theta_mat = zeros(eor,perm_len);
            coeffs_init = ones(perm_len,1);
            [Pxu_init,Oxu_init] = deal(order_mat);
            Lxuh_init = Lxu;
            Pdxu_init = zeros(1,eor*ndep);

            specs = struct( 'eor',eor, ...
                            'ndep',ndep, ...
                            'bor', bor, ...
                            's',sflat, ...
                            'xuP',xuP, ...
                            'fmm',fmm, ...
                            'pc',poly_coeffs, ...
                            'ic',icoeff_mat);

            Jac_vu = zeros(ndim,ndim,nbse);
            [v_dxu_dude_theta_mat Jac_vu] = compute_dude_scalar_Jac(specs,1,dxuc_init,v_dxu_dude_theta_mat, ...
                Jac_vu,theta_all_tens,coeffs_init,Oxu_init,Pxu_init,Lxuh_init,Pdxu_init);

            v_dxu_dude_theta = nan(eor*ndep,nbse);
            for ibse = 1:nbse
                v_dxu_dude_theta(:,ibse) = reshape((v_dxu_dude_theta_mat*theta_u_tens(:,:,ibse))',[],1);
            end

            v_dxu_dxde_mat = zeros(eor*ndep,perm_len);
            Jac_vx = zeros(ndim,ndim,nbse);
            Jac_L_xu_theta_out = zeros(nvar,perm_len);

            bor_p1 = bor + 1;

            nDzi = obj.not_D_zero_inds;
            idx = nDzi(:,1);
            if (sum(idx))
                dxc_p1 = dxuc_init(1) + 1;
                idx_p1 = dxc_p1 + 1;
                idx_front = 1:(bor_p1-dxc_p1);
                idx_back = idx_p1:bor_p1;

                Oxu_dx = Oxu_init(idx,:);
                c_dx = (-1.0*coeffs_init(idx))*fmm(1);

                Pxu_dx = [Pxu_init(idx,1)-1, Pxu_init(idx,2:end)];

                dxLxuh = Lxuh_init(idx,:);
                dxLxuh(:,1) = comp_dnL(    dxc_p1,xuP(idx_front,1),poly_coeffs(idx_back,idx_back), ...
                                                icoeff_mat(idx_back,idx_p1),Oxu_dx(:,1));

                dx_dxuc = dxuc_init;
                dx_dxuc(1) = dxc_p1;

                Jac_x_x = ((c_dx).*prod(dxLxuh,2))';
                Jac_x_x_theta = -1.0*Jac_x_x;
                Jac_vx(1,1,:) = pagemtimes(Jac_x_x_theta,theta_all_tens(idx,1,:));
                Jac_L_xu_theta_out(1,idx) = Jac_x_x_theta;

                for kder = 1:eor
                    kder_inds = ((ndep*(kder-1))+1):(ndep*kder);
                    v_dxu_dxde_mat(kder_inds,idx) = v_dxu_dxde_mat(kder_inds,idx) + (pdxu(kder_inds)')*Jac_x_x;
                    [v_dxu_dxde_mat(:,idx) Jac_vx] = compute_dude_vector_Jac(specs,kder+1,kder,dx_dxuc, ...
                        v_dxu_dxde_mat(:,idx),Jac_vx,theta_all_tens(idx,:,:),c_dx,Oxu_dx,Pxu_dx,dxLxuh,Pdxu_init);
                end
            end

            for iu = 1:ndep
                idui = nDzi(:,iu+1);
                if (sum(idui))
                    duic_p1 = dxuc_init(iu+1) + 1; % increase the ui derivative count by one
                    idui_p1 = duic_p1 + 1; % index of relevant indices for derivative evaluation
                    idui_front = 1:(bor_p1-duic_p1);
                    idui_back = idui_p1:bor_p1;

                    Oxu_dui = Oxu_init(idui,:);
                    c_dui = (-1.0*coeffs_init(idui))*fmm(iu+1); % update coefficients

                    Pxu_dui = Pxu_init(idui,:); % initialize power records to selected previous values
                    Pxu_dui(:,iu+1) = Pxu_dui(:,iu+1)-1; % initialize powers of ui

                    duiLxuh = Lxuh_init(idui,:); % initialize L functional history to selected previous values
                    % update L functional history to reflect derivative in x
                    duiLxuh(:,iu+1) = comp_dnL( duic_p1,xuP(idui_front,iu+1),poly_coeffs(idui_back,idui_back), ...
                                                icoeff_mat(idui_back,idui_p1),Oxu_dui(:,iu+1));

                    dui_dxuc = dxuc_init;
                    dui_dxuc(iu+1) = duic_p1;

                    Jac_x_ui = ((c_dui.*prod(duiLxuh,2))');
                    Jac_x_ui_h = -1.0*Jac_x_ui;
                    Jac_vx(1,iu+1,:) = pagemtimes(Jac_x_ui_h,theta_all_tens(idui,1,:));
                    Jac_L_xu_h_out(iu+1,idui) = Jac_x_ui_h;

                    Pdxu_dui = Pdxu_init; % initialize dxu power records to previous values
                    Pdxu_dui(iu) = Pdxu_init(iu) + 1; % update dxu power records to reflect derivative in ui
                    for kder = 1:eor
                        kder_inds = ((ndep*(kder-1))+1):(ndep*kder);
                        v_dxu_dxde_mat(kder_inds,idui) = v_dxu_dxde_mat(kder_inds,idui) + ((pdxu(kder_inds)')*Jac_x_ui)*pdxu(iu); % add contribution of freshly computed t derivative term
                        [De_dxu_dxde_mat(:,idui) Jac_vx] = compute_dude_vector_Jac(specs,kder+1,kder,dui_dxuc, ...
                            v_dxu_dxde_mat(:,idui),Jac_vx,theta_all_tens(idui,:,:),c_dui, ...
                            Oxu_dui,Pxu_dui,duiLxuh,Pdxu_dui);
                    end
                end
            end

            Jac_v_out = Jac_vx + Jac_vu;

            v_dxu_dxde_theta = v_dxu_dxde_mat*theta_x_mat;

            v_dxu_theta = v_dxu_dude_theta + v_dxu_dxde_theta;

            v_out_theta_mat = [v_xu_theta; v_dxu_theta];
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    [eor_ ndep_ bor_ p_ txP_ fmm_ pc_ ic_] = deal(str.eor, str.ndep, str.bor, str.s, str.xuP, str.fmm, str.pc, str.ic);
end

%{
    NOTE: methods below use outdated variable naming convention -
        t := indep var, x := dep var, h/H := parameters / parameter matrix, De := vector field value, v
%}

function Ltx_out = comp_Lx_Lu(combo_len_, bor_, txP_, poly_coeffs_, poly_ords_)
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

function [De_xp_out Jac_out] = compute_dude_vector_Jac(specs_,k_,source_,dtxc_,De_xp_in_,Jac_,H_,c_,Otx_,Ptx_,Ltxh_,Pxp_)
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

            [De_xp_k(:,idt) Jac_v_k] = compute_dude_vector_Jac(specs_,k_+1,source_,dt_dtxc,De_xp_k(:,idt),Jac_v_k,H_(idt,:,:),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_);
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

                [De_xp_k(:,idxi) Jac_v_k] = compute_dude_vector_Jac(specs_,k_+1,source_,dxi_dtxc,De_xp_k(:,idxi),Jac_v_k,H_(idxi,:,:),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
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

                    [De_xp_k Jac_v_k] = compute_dude_vector_Jac(specs_,k_+1,source_,dtxc_,De_xp_k,Jac_v_k,H_,c_dxpi,Otx_,Ptx_,Ltxh_,Pxp_dxpi);
                end
            end
        end

        Jac_dkm1_dxi_vec = ((c_.*Ltxh_eval)')*prod(xp_.^Pxp_,2);

        De_xp_k(kder_inds,:) = De_xp_k(kder_inds,:) + (xp_(source_inds+ndep_)')*Jac_dkm1_dxi_vec;
        Jac_v_k(iJ_k,source_inds+ndep_+1,:) = Jac_v_k(iJ_k,source_inds+ndep_+1,:) + pagemtimes(Jac_dkm1_dxi_vec,H_(:,1,:));

        [De_xp_k Jac_v_k] = compute_dude_vector_Jac(specs_,k_+1,k_,dtxc_,De_xp_k,Jac_v_k,H_,c_,Otx_,Ptx_,Ltxh_,Pxp_);

        De_xp_out = De_xp_in_ + De_xp_k;
        Jac_out = Jac_ + Jac_v_k;
    end
end

function [De_xp_out Jac_out] = compute_dude_scalar_Jac(specs_,k_,dtxc_,De_xp_in_,Jac_,H_,c_,Otx_,Ptx_,Ltxh_,Pxp_)
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

            [De_xp_k(:,idt) Jac_v_k] = compute_dude_scalar_Jac(specs_,k_+1,dt_dtxc,De_xp_k(:,idt),Jac_v_k,H_(idt,:,:),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_);
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

                [De_xp_k(:,idxi) Jac_v_k] = compute_dude_scalar_Jac(specs_,k_+1,dxi_dtxc,De_xp_k(:,idxi),Jac_v_k,H_(idxi,:,:),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
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

                    [De_xp_k Jac_v_k] = compute_dude_scalar_Jac(specs_,k_+1,dtxc_,De_xp_k,Jac_v_k,H_,c_dxpi,Otx_,Ptx_,Ltxh_,Pxp_dxpi);
                end
            end
        end
        De_xp_out = De_xp_in_ + De_xp_k;
        Jac_out = Jac_ + Jac_v_k;
    end
end

function De_xp_out = compute_dude_vector(specs_,k_,source_,dtxc_,De_xp_in_,c_,Otx_,Ptx_,Ltxh_,Pxp_)
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
            De_xp_k(:,idt) = compute_dude_vector(specs_,k_+1,source_,dt_dtxc,De_xp_k(:,idt),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_);
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
                De_xp_k(:,idxi) = compute_dude_vector(specs_,k_+1,source_,dxi_dtxc,De_xp_k(:,idxi),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
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
                    De_xp_k = compute_dude_vector(specs_,k_+1,source_,dtxc_,De_xp_k,c_dxpi,Otx_,Ptx_,Ltxh_,Pxp_dxpi);
                end
            end
        end

        De_xp_k(kder_inds,:) = De_xp_k(kder_inds,:) + (xp_(source_inds+ndep_)')*((c_.*Ltxh_eval)')*prod(xp_.^Pxp_,2);
        De_xp_k = compute_dude_vector(specs_,k_+1,k_,dtxc_,De_xp_k,c_,Otx_,Ptx_,Ltxh_,Pxp_);
        De_xp_out = De_xp_in_ + De_xp_k;
    end
end

function De_xp_out = compute_dude_scalar(specs_,k_,dtxc_,De_xp_in_,c_,Otx_,Ptx_,Ltxh_,Pxp_)
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
            De_xp_k(:,idt) = compute_dude_scalar(specs_,k_+1,dt_dtxc,De_xp_k(:,idt),c_dt,Otx_dt,Ptx_dt,dtLtxh,Pxp_);
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
                De_xp_k(:,idxi) = compute_dude_scalar(specs_,k_+1,dxi_dtxc,De_xp_k(:,idxi),c_dxi,Otx_dxi,Ptx_dxi,dxiLtxh,Pxp_dxi);
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
                    De_xp_k = compute_dude_scalar(specs_,k_+1,dtxc_,De_xp_k,c_dxpi,Otx_,Ptx_,Ltxh_,Pxp_dxpi);
                end
            end
        end
        De_xp_out = De_xp_in_ + De_xp_k;
    end
end
