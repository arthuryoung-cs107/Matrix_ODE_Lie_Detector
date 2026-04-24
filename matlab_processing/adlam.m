% classdef adlam < adobj
classdef adlam

    properties (Constant)

        %% quick handles for simply processed quantities, inline fast
        Omap = @(o_) o_.Jac * o_.xu + o_.b;
        vevl = @(o_,v_) o_.Jac * v_(:);

        ndep = @(o_) length( o_.xu )-1;
        kor = @(o_) size( o_.dxu , 2 );
        nvar = @(o_) length( o_.xu );
        ndim = @(o_) 1 + ( length(o_.xu)-1 )*( size( o_.dxu , 2 )+1 );

        bspc_dims = @(o_) deal( adlam.ndep(o_), adlam.nvar(o_) );

        Pmat = @(o_) o_.lam.Pmat;
        Plen = @(o_) size(adlam.Pmat(o_),2);

        %% handles for processed quantities of moderate complexity, still inline fast
        Lmat_full = @(o_) [ ones( adlam.nvar(o_), 1) , o_.Lmat];
        dkL_mat_full = @(o_,k_) [ zeros( adlam.nvar(o_), 1) , o_.dL_tns(:,:,k_) ];

        iLmatP = @(sz_,Pl_,Pm_) sub2ind(sz_,reshape((1:sz_(1))' .* ones(sz_(1),Pl_),1,[]), Pm_(:)'+1 ) ;
        LmPf_vec = @(o_,Lmf_) Lmf_(adlam.iLmatP(size(Lmf_),adlam.Plen(o_),adlam.Pmat(o_)));

        LmatP = @(o_) reshape( adlam.LmPf_vec(o_,adlam.Lmat_full(o_)), adlam.nvar(o_), adlam.Plen(o_) );
        dkLmatP = @(o_,k_) reshape( adlam.LmPf_vec(o_,adlam.dkL_mat_full(o_,k_)), adlam.nvar(o_), adlam.Plen(o_) );

        lrowP = @(o_) prod(adlam.LmatP(o_),1);
        lvec = @(o_) adlam.lrowP(o_)';

    end

    properties

        xu;
        lam;

        A;
        b;

        s0;

        Lmat;

        dL_tns;

        pr0;

        dxu;
        pr;

    end

    methods
        function obj = adlam(lam_,xu_,dxu_)
            if ( isfield(lam_,'Omap_A') )
                A_ = lam_.Omap_A;
            else
                A_ = eye(length(xu_(:)));
            end
            if ( isfield(lam_,'Omap_b') )
                b_ = lam_.Omap_b;
            else
                b_ = zeros(length(xu_(:)),1);
            end
            if ( isfield(lam_,'Pmat') )
                Pmat = lam_.Pmat;
                bor = max(Pmat(:));
            elseif ( isfield(lam_,'bor') )
                bor = lam_.bor;
                [~, Pmat, ~] = ldaux.count_set_P_len(bor, length(xu_(:)) );
            else
                bor = 3; % default to cubic
                [~, Pmat, ~] = ldaux.count_set_P_len(bor, length(xu_(:)) );
            end
            if ( nargin==3 )
                dxu = reshape(dxu_, length(xu_(:))-1, []);
                kor = size(dxu,2);
            else
                dxu = [];
                if ( isfield(lam_,'kor') )
                    kor = max([ 1 lam_.kor ]) ;
                else
                    kor = 1;
                end
            end

            xu = xu_(:);
            lam = lam_;
            b =  b_(:);
            A =  reshape( A_, length(b), length(xu) );
            s0 = adobj( A*xu + b, A );

            %% stage prolongation workspace
            % bor = double(lam_.bor);
            % Pmat = double(lam_.Pmat);
            Plen = size(Pmat,2);

            nvar = size(Pmat,1);
            ndep = nvar - 1;
            Oa = diag(A);
            xu_O = s0.val;

            Lmat = nan(nvar,bor); % matrix of values for each L polynomial
            c_tns = nan(nvar,bor,kor); % aggregate chain coefficient
            dL_tns = nan(nvar,bor,kor); % derivatives of ordered single variate monomials

            Lmat(:,1) = xu_O; % image of xu under diagonal linear map A * [x ; u] + b into O
            c_tns(:,1,1) = Oa; % first derivative chain coefficient
            dL_tns(:,1,1) = c_tns(:,1,1); % first order derivatives of linear map are just chain coefficients
            c_tns(:,1,2:end) = 0; % onward chain coefficients are zero
            dL_tns(:,1,2:end) = 0; % onward dL value is zero
            for ib = 2:bor
                Lmat(:,ib) = Lmat(:,ib-1) .* xu_O;
                c_tns(:,ib,1) = ib*Oa;
                dL_tns(:,ib,1) = c_tns(:,ib,1) .* Lmat(:,ib-1);
                for ik = 2:kor
                    if (ik<=ib)
                        c_tns(:,ib,ik) = c_tns(:,ib,ik-1) .* ((ib-ik+1)*Oa);
                        if (ik<ib)
                            dL_tns(:,ib,ik) = c_tns(:,ib,ik) .* Lmat(:,ib-ik);
                        else
                            dL_tns(:,ib,ik) = c_tns(:,ib,ik);
                        end
                    else
                        c_tns(:,ib,ik) = 0;
                        dL_tns(:,ib,ik) = 0;
                    end
                end
            end

            %% core initializations
            obj.xu = xu;
            obj.lam = lam;

            obj.A =  A;
            obj.b =  b;

            obj.s0 = s0;

            obj.Lmat = Lmat;
            % obj.c_tns = c_tns;
            obj.dL_tns = dL_tns;

            %% prologation initializations
            LmatP = adlam.LmatP(obj);
            d1LmatP = adlam.dkLmatP(obj,1);

            % (transposed) Jacobian of l over the base space
            JlT = zeros( nvar,Plen );
            for iv = 1:nvar
                JlT(iv,:) = 1;
                for iiv = 1:(iv-1)
                    JlT(iv,:) = JlT(iv,:) .* LmatP(iiv,:);
                end
                JlT(iv,:) = JlT(iv,:) .* d1LmatP(iv,:);
                for iiv = (iv+1):nvar
                    JlT(iv,:) = JlT(iv,:) .* LmatP(iiv,:);
                end
            end

            % the base space is fully characterized by the lambda vector and its Jacobian
            pr0 = struct( ...
                'l', adobj( prod( LmatP, 1 )', JlT' ), ...
                'v', @(o_,th_) (reshape( th_, length( o_.l.val ), [] ))' * o_.l.val ...
            );
            obj.pr0 = pr0;

            pr0_data = struct( ...
                            'ndep', ndep, ...
                            'Plen', Plen, ...
                            'Pmat', Pmat, ...
                            'inds', 1:Plen ... % all indices contribute to base space image
                            );
            %% if derivatives are provided, prolong observed jet space
            if ( size(dxu,2) )
                obj = adlam.prolong_mvpolynomial(obj,dxu,pr0_data);
            else
                obj.dxu = [];
                obj.pr = [];
            end
        end

    end

    methods (Static)

        function obj_out = prolong_mvpolynomial(obj,dxu_,pr0d_)
            [ndep, nvar] = adlam.bspc_dims(obj);
            Pmat = adlam.Pmat(obj);
            Plen = size(Pmat,2);

            dxu = reshape( dxu_,ndep,[] );
            kor = size(dxu,2);
            ndim = 1+ndep*(kor+1);

            LmatP = adlam.LmatP(obj);
            dkLtnsP = nan(nvar,Plen,kor);
            dkLtnsP_full = nan(nvar,Plen,kor+1);
            dkLtnsP_full(:,:,1) = LmatP;
            for k = 1:kor
                % dkLtnsP(:,:,k) = adlam.dkLmatP(obj,k);
                [dkLtnsP(:,:,k),dkLtnsP_full(:,:,k+1)] = deal(adlam.dkLmatP(obj,k));
            end

            pr_data = pr0d_;

            pr_data.kor = kor;
            pr_data.ndim = 1+ndep*(kor+1);
            pr_data.LmatP = LmatP;
            pr_data.dkLtnsP = dkLtnsP;
            pr_data.nder = [ ones(nvar,1) ; zeros(ndim-nvar,1) ]; % initialize derivative block

            pr_data.dxu = dxu;
            pr_data.tau_prk = @(o_,k_) [1 ; reshape( o_.dxu(:,1:k_) , [], 1 )];

            dkLP = @(iv_,k_) dkLtnsP_full(iv_,:,k_);


            dkxl = zeros(kor,Plen);
            Jdkxl = zeros(ndim,Plen,kor);

            dkxl(1,:) = ( pr_data.tau_prk(pr_data,1)' )*(obj.pr0.l.Jac');
            % compute Jacobian of resultant dkxl term, pass to second prolongation




            for iv = 1:nvar % compute gradient of each Jacobian row
                Jdkxl(iv,:,1) = 1;
                for iiv = 1:(iv-1)
                    Jdkxl(iv,:,1) = Jdkxl(iv,:,1) .* dkLtnsP(iiv,:,1);
                end
                Jdkxl(iv,:,1) = Jdkxl(iv,:,1) .* dkLtnsP(iv,:,2);
                for iiv = (iv+1):nvar
                    Jdkxl(iv,:,1) = Jdkxl(iv,:,1) .* LmatP(iiv,:);
                end
            end


            pause

            % propogate base space partial x term

            %


            lkx_out = zeros(prd_.ndep,prd_.Plen,prd_.kor);
            Jlkx_out = zeros(prd_.ndep,prd_.Plen,prd_.kor,prd_.ndim);

            % [dkxl,Jdkxl,lkx,Jlkx] = adlam.prolong_vxu_mvpolynomial(obj,pr_data);

            pr = 0;
            % pause

            obj_out = obj;
            % obj_out.dxu = dxu;
            obj_out.pr = pr;
        end
        function [dkxl_out,Jdkxl_out] = prolong_vu_mvpolynomial(obj,k_,prd_)
            dkxl_out = zeros(prd_.kor,prd_.Plen);
            Jdkxl_out = zeros(prd_.ndim,prd_.Plen,prd_.kor);

        end
        function vec_out = coordgrads2Jac(mat_)
            [M,nobs] = size(mat_);
            N = size( mat_(1).Jac,2 );
            for iobs = 1:nobs
                val_i = nan(M,1);
                Jac_i = nan(M,N);
                for i = 1:M
                    val_i(i) = mat_(i,iobs).val;
                    Jac_i(i,:) = mat_(i,iobs).Jac;
                end
                vec_out(iobs) = adobj(val_i,Jac_i);
            end
        end
        function obj_o = trunc_to_adobj(obj)
            obj_o = adobj(obj.val,obj.Jac);
        end

    end
end


            % obj.pr = [];


            %% recursive prolongation from the base space
            % Lmat_full = [ ones(nvar,1) , Lmat ];
            % LmatP = reshape(Lmat_full(sub2ind([nvar,bor+1],reshape((1:nvar)'.*ones(nvar,Plen),1,[]),Pmat(:)'+1)),nvar,Plen)
            % lrowP = prod(LmatP,1)
            % LmatP = adlam.LmatP(obj)
            % lrowP = adlam.lrowP(obj)

            % pr0_init = struct( ...
            %     'k', 0, ...
            %     'Pmat', Pmat, ...
            %     'nder', zeros(ndim,1), ...
            %     'LmatP', [ LmatP ; ones(ndim-nvar,Plen) ], ...
            %     'dkxlP', [ lrowP ; zeros(kor,Plen) ] ...
            % );
            % obj.pr = adlam.prolong_mvpolynomial(obj,pr0_init);

            % obj.pr = adlam.prolong_mvpolynomial(obj);

            %% relegated
            % Lmat_full = [ ones(nvar,1) , Lmat ];

            % tau_k = [1 ; dxu(:) ];
            % tau_0 = tau_k(1:nvar);
            %
            % tau_k_ad = adobj(tau_k, [ zeros(1,ndim) ; zeros(ndim-nvar,nvar), eye(ndim-nvar) ] );
            % tau_0_ad = tau_k_ad.qdim(1:nvar)

            % iiv = (1:nvar)';
            % iiD_base = sub2ind([nvar,nvar],1:nvar,1:nvar);
            % LmatP = ones(nvar,Plen);
            % lvecP = ones(Plen,1);
            % J_dkx_ltnsP = zeros( Plen,ndim, kor+1 );
            % dkx_lmatP = zeros( Plen,kor );
            %
            % dkx_l_mat = zeros( Plen,kor );
            % Jdkx_l_tns = zeros( Plen,kor,ndim );
            %
            % for ip = 1:Plen
            %
            %     iip = Pmat( :,ip )
            %     iip_nz = iip~=0
            %     if (sum(iip_nz))
            %         iiL = sub2ind([nvar,bor+1], iiv, iip+1 ) % linear indices of LmatP_i values in Lmat_full
            %         LmatP(:,ip) = Lmat_full( iiL ) %
            %         lvecP(ip) = prod( LmatP(:,ip) )
            %
            %
            %
            %         dL_mat_i_full = zeros(nvar,bor+1);
            %         dL_mat_i_full(:,2:end) = dL_tns(:,:,1);
            %
            %         gLi = LmatP(:,ip)*ones(1,nvar);
            %         gLi( iiD_base ) = dL_mat_i_full( iiL );
            %         gli = prod(gLi,1);
            %         J_dkx_ltnsP(ip,1:nvar,1) = gli;
            %
            %         dkx_lmatP(ip,1) = J_dkx_ltnsP(ip,1:nvar,1) * tau_0;
            %
            %         li_ad = adobj( LmatP(1,ip), [ dL_mat_i_full(iiL(1)) , zeros(1, ndim-1) ] );
            %         for iv = 2:nvar
            %             li_ad = li_ad * adobj( LmatP(iv,ip), [ zeros(1,iv-1) , dL_mat_i_full(iiL(iv)), zeros(1,ndim-iv) ] );
            %         end
            %
            %         % compute grad (d1xl) = grad ( grad l * tau_u )
            %
            %         pause
            %
            %         li_ad.val,li_ad.Jac
            %         tau_k_ad.val, tau_k_ad.Jac
            %
            %         % gli_ad = adobj( li_ad.Jac(:),prod(,1) );
            %         % d1xli_ad = .qdim()
            %         %
            %         % for iv = 2:nvar
            %         %
            %         % end
            %
            %         % d1xli_ad.val,d1xli_ad.Jac
            %
            %         % J_dkx_ltnsP(ip,:,2) =
            %
            %         for ik = 2:kor
            %
            %             dkx_lmatP(ip,ik) = J_dkx_ltnsP(ip,:,ik) * tau_k;
            %
            %
            %         end
            %
            %         pause
            %
            %     end

                % J_dkx_ltnsP

                % iip_nz = iip~=0;
                % if (sum(iip_nz))
                %     iiL = sub2ind([nvar,bor], iiv(iip_nz) , iip(iip_nz)' ) % linear indices of non unitary LmatP_i values in
                %     LmatP(:,ip) = Lmat( iiL );
                %     lvecP(ip) = prod( LmatP(:,ip) );
                %
                %     dL_mat_i_full = zeros(nvar,bor+1);
                %     dL_mat_i_full(:,2:end) = dL_tns(:,:,1)
                %
                %     dL_mat_i_full(iiL)
                %
                %     gli = LmatP(:,ip)*ones(1,nvar)
                %     gli( iiD_base(iip_nz) ) = dL_mat_i_full( iiL )
                %
                %     pause
                %
                % end
            % end

            % LmatP
            % lvecP



            % function prn = prolong_mvpolynomial(s0_,dxu_,bor_)
            %     bor = double(bor_);
            %     nvar = length(s0_.xu);
            %     ndep = nvar-1;
            %     dxu = reshape(dxu_,ndep,[]);
            %     kor = size(dxu,2);
            %     kp1 = kor+1;
            %
            %     xu = s0_.val;
            %     Oa = diag(s0_.Jac);
            %
            %     c_tns = nan(nvar,bor,kor+1);
            %     L_cell = cell(nvar,bor,kor+1);
            %     for iv = 1:nvar
            %         ai = Oa(iv);
            %         L_cell{iv,1,1} = adobj( xu(iv) , ai );
            %         for ik = 2:kp1
            %             L_cell{iv,1,ik} = adobj(L_cell{iv,1,ik-1}.Jac, 0.0);
            %         end
            %
            %         for ib = 2:bor
            %             % compute L_ib and its first partial
            %             L_cell{iv,ib,1} = L_cell{iv,1,1}.*L_cell{iv,ib-1,1};
            %
            %             for ik = 2:kp1
            %                 if (ik<ib)
            %                     % L_cell{iv,ib,ik} = ...
            %                     %     adobj( ...
            %                     %     L_cell{iv,ib,ik-1}.Jac, ...
            %                     %
            %                     %      );
            %                 elseif (ik==ib)
            %
            %                 else
            %
            %                 end
            %             end
            %
            %             L_cell{iv,ib,ik} = adobj( pLs(iv,ib,ik-1), pLs(iv,ib,ik) );
            %
            %
            %             % second partial to ib-1
            %             for ik = 2:(ib-1)
            %                 pLs(iv,ib,ik) = ...
            %                     pLs(iv,ib,ik-1) ...
            %                     *( ai*(ib-ik+1) ) ...
            %                     *( L_cell{iv,ib-ik,1}.val );
            %                 L_cell{iv,ib,ik} = adobj( pLs(iv,ib,ik-1), pLs(iv,ib,ik) );
            %             end
            %
            %             if (ik)
            %             % partial of same order
            %             ik = ib;
            %             pLs(iv,ib,ik) = ...
            %                 pLs(iv,ib,ik-1) ...
            %                 *( ai*(ib-ik+1) ) ;
            %             L_cell{iv,ib,ik} = adobj( pLs(iv,ib,ik-1), pLs(iv,ib,ik) );
            %
            %             for ik = (ib+1):kp1
            %                 L_cell{iv,ib,ik} = adobj( L_cell{iv,ib,ik-1}.Jac, 0.0);
            %             end
            %         end
            %     end
            %
            %
            %     Lk_cell{1,1} = adobj(sh.val,eye(nvar));
            %     for ik = 1:kor
            %         % iik = ik+1;
            %         % Lk_cell{1,}
            %     end
            %     for ib = 2:bor_
            %         Lk_cell{1,ib} = Lk_cell{1,1}.*Lk_cell{1,ib-1};
            %         for ik = 1:kor
            %
            %         end
            %     end
            % end

            % function d1xl_out = d1xl(l_,d1xu_)
            %     ndep = length(l_.xu)-1;
            %     vo = adlam.vevl(l_,d1xu_);
            %     Jo = [,];
            %     d1xl_out = adobj(vo,Jo);
            % end
