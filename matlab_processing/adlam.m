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
        % pr;
        dkxl;
        Jdkxl;
        lkx;
        Jlkx;

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

            korp1 = kor+1; % compute L partials up to n+1 for the computation of n'th prolongation Jac
            c_tns = nan(nvar,bor,korp1); % aggregate chain coefficient
            dL_tns = nan(nvar,bor,korp1); % derivatives of ordered single variate monomials

            Lmat(:,1) = xu_O; % image of xu under diagonal linear map A * [x ; u] + b into O
            c_tns(:,1,1) = Oa; % first derivative chain coefficient
            dL_tns(:,1,1) = c_tns(:,1,1); % first order derivatives of linear map are just chain coefficients
            c_tns(:,1,2:end) = 0; % onward chain coefficients are zero
            dL_tns(:,1,2:end) = 0; % onward dL value is zero
            for ib = 2:bor
                Lmat(:,ib) = Lmat(:,ib-1) .* xu_O;
                c_tns(:,ib,1) = ib*Oa;
                dL_tns(:,ib,1) = c_tns(:,ib,1) .* Lmat(:,ib-1);
                for ik = 2:korp1
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

            %% if derivatives are provided, prolong observed jet space
            if ( size(dxu,2) )

                % obj = adlam.prolong_mvpolynomial(obj,dxu,pr0_data);
                obj.dxu = dxu;
                kor = size(dxu,2);
                ndim = 1+ndep*(kor+1);

                dkLtnsP = nan(nvar,Plen,korp1+1);
                dkLtnsP = nan(nvar,Plen,korp1+1);
                dkLtnsP(:,:,1) = LmatP;
                for k = 1:korp1
                    dkLtnsP(:,:,k+1) = adlam.dkLmatP(obj,k);
                end
                dkLP = @(iv_,k_) dkLtnsP(iv_,:,k_+1);

                pr0_data = struct( ...
                'ndep', ndep, ...
                'nvar', nvar, ...
                'Plen', Plen, ...
                'Pmat', Pmat, ...
                's', [xu ; dxu(:) ], ...
                'kor', kor, ...
                'ndim', ndim, ...
                'LmatP', LmatP, ...
                'dkLtnsP', dkLtnsP, ...
                'nder', zeros(ndim,1), ...
                'c_net', 1, ... % power rule coefficient for partials of L_dkxu
                'Pmat_dkxu', zeros(ndep,kor), ... % no derivatives in 0'th prolongation
                'inds', 1:Plen, ... % all indices contribute to 0'th prolongation
                'Ldkxu', @(o_,dxu_) o_.c_net*prod( reshape((dxu_).^(o_.Pmat_dkxu),[],1) ), ...
                'dkxu_inds', @(o_,k_) o_.nvar+sub2ind([o_.ndep,o_.kor],1:o_.ndep,k_*ones(1,o_.ndep) ) ...
                );
                obj.dkxl = zeros(kor,Plen);
                obj.Jdkxl = zeros(ndim,Plen,kor);
                obj.lkx = zeros(ndep,Plen,kor);
                obj.Jlkx = zeros(ndep,Plen,kor,ndim);
                for iv = 1:nvar
                    obj = adlam.prolong_vu_mvpolynomial(obj,1,iv,JlT(iv,:),pr0_data);
                    obj = adlam.prolong_vx_mvpolynomials(obj,iv,JlT(iv,:),pr0_data);
                end
            else
                obj.dxu = [];
                obj.dkxl = [];
                obj.Jdkxl = [];
                obj.lkx = [];
                obj.Jlkx = [];
            end
        end

    end

    methods (Static)
        function obj_out = prolong_vx_mvpolynomials(obj,ivp_,pvl_,prd_)
            obj_out = obj;
            prd_k = prd_;
            kor = prd_.kor;
            nvar = prd_k.nvar;

            %% first step: accumulate input partial derivative term's contribution to lkx. Handle product rule.
            % step 1a: accumulate input partial derivative term's contribution to lkx
            if (ivp_>1) % if the input partial is with respect to a dependent variable
                prd_k.Pmat_dkxu(ivp_-1) = prd_k.Pmat_dkxu(ivp_-1) + 1; % increment polynomial derivative power
            end
            % evaluate input chained dkxu polynomial
            Ldkxu_k0 = prd_k_k.Ldkxu(prd_k,obj.dxu); % product of derivative polynomials. Scalar
            Ldkxu_kderset = obj.dxu * Ldkxu_k0;
            for k = 1:kor
                prolong_vx_mvpolynomial(obj,k,ivp_,pvl_,prd_)
            end

            obj_out.lkx(:,prd_k.inds,k_) = obj_out.lkx(:,prd_k.inds,k_) + Ldkxu_k*pvl_; % accumulate total der

            Ldkxu_k = obj.dxu(:,src_)*Ldkxu_k0; % aggregate product of derivative polynomials times src derivatives. Vector
            obj_out.lkx(:,prd_k.inds,k_) = obj_out.lkx(:,prd_k.inds,k_) + Ldkxu_k*pvl_; % accumulate total der
            prd_k.nder(ivp_) = prd_k.nder(ivp_) + 1; % increment count of partial derivatives of ivp_ variable

            % step 1b: accumulate first half of product rule to lkx Jacobian,
            Ldkxu_k05 = Ldkxu_k0*pvl_;
            inds_src = prd_.dkxu_inds(prd_,src_);
            obj_out.Jlkx( :, prd_.inds, k_, inds_src ) = obj_out.Jlkx( :, prd_.inds, k_, inds_src ) ...
                                                            + ones(prd_.ndep,1)*Ldkxu_k05;
            if (k_<prd_.kor)
                obj_out = prolong_vx_mvpolynomial(obj_out,k_+1,src_+1,ivp_,pvl_,prd_k);
            end

        end
        function obj_out = prolong_vx_mvpolynomial(obj,k_,src_,ivp_,pvl_,prd_)
            obj_out = obj;
            prd_k = prd_;
            nvar = prd_k.nvar;

            %% first step: accumulate input partial derivative term's contribution to lkx.
            % step 1a: accumulate input partial derivative term's contribution to lkx
            if (ivp_>1) % if the input partial is with respect to a dependent variable
                prd_k.Pmat_dkxu(ivp_-1) = prd_k.Pmat_dkxu(ivp_-1) + 1; % increment polynomial derivative power
            end
            % evaluate input chained dkxu polynomial
            Ldkxu_k0 = prd_k_k.Ldkxu(prd_k,obj.dxu); % product of derivative polynomials. Scalar
            Ldkxu_k = obj.dxu(:,src_)*Ldkxu_k0; % aggregate product of derivative polynomials times src derivatives. Vector
            obj_out.lkx(:,prd_k.inds,k_) = obj_out.lkx(:,prd_k.inds,k_) + Ldkxu_k*pvl_; % accumulate total der
            prd_k.nder(ivp_) = prd_k.nder(ivp_) + 1; % increment count of partial derivatives of ivp_ variable

            % step 1b: accumulate first half of product rule to lkx Jacobian,
            Ldkxu_k05 = Ldkxu_k0*pvl_;
            inds_src = prd_.dkxu_inds(prd_,src_);
            obj_out.Jlkx( :, prd_.inds, k_, inds_src ) = obj_out.Jlkx( :, prd_.inds, k_, inds_src ) ...
                                                            + ones(prd_.ndep,1)*Ldkxu_k05;
            if (k_<prd_.kor)
                obj_out = prolong_vx_mvpolynomial(obj_out,k_+1,src_+1,ivp_,pvl_,prd_k);
            end

            %% second step: accumulate input's contribution to the Jacobian of lkx
            % prd_.dkLtnsP
            dkLP = @(iv_,ii_,k_) prd_k.dkLtnsP(iv_,ii_,k_+1);
            % step 2a: accumulate this term's base space partial derivatives
            % identify base space lambda indices which have nonzero contribution to this order
            iipnz = ( prd_k.Pmat(:,prd_k.inds) - prd_k.nder(1:nvar) ) > 0;
            % accumulate partial derivatives over base space variables
            for iv = 1:nvar
                if ( iipnz(iv,:) )
                    % accumulate partials over the base space
                    inds_iv = prd_k.inds( iipnz(iv,:) );
                    pvl_iv = ones( 1,length(inds_iv) );
                    for iiv = 1:(iv-1)
                        pvl_iv = pvl_iv .* dkLP(iiv, inds_iv, prd_k.nder(iiv));
                    end
                    pvl_iv = pvl_iv .* dkLP(iv, inds_iv, prd_k.nder(iv) + 1);
                    for iiv = (iv+1):nvar
                        pvl_iv = pvl_iv .* dkLP(iiv, inds_iv, prd_k.nder(iiv));
                    end
                    % pvl_iv

                    % accumulate partial derivative of dkxl wrt to base space iv
                    % obj_out.Jdkxl(iv,inds_iv,k_) = obj_out.Jdkxl(iv,inds_iv,k_) + Ldkxu_k*pvl_iv;
                    obj_out.Jlkx(:,inds_iv,k_,iv) = obj_out.Jlkx(:,inds_iv,k_,iv) + Ldkxu_k*pvl_iv;

                    if (k_<prd_.kor)
                        prd_k_piv.inds = inds_iv;
                        obj_out = prolong_vx_mvpolynomial(obj_out,k_+1,src_,iv,pvl_iv,prd_k_piv);
                    end
                end
            end

            % step 2b: accumulate this term's base space partial derivatives
            % identify jet space lambda indices which have nonzero contribution to this order
            iipdxu_nz = prd_k.Pmat_dkxu > 0
            ndep = nvar-1;
            iv = nvar;
            for k = 1:k_
                for idep = 1:ndep
                    iv = iv + 1;
                    if ( iipdxu_nz(idep,k) )
                        prd_k_pdkxui = prd_k;

                            % accumulate power rule coefficient
                        prd_k_pdkxui.c_net = prd_k_pdkxui.c_net * prd_k_pdkxui.Pmat_dkxu(idep,k);
                            % decrement derivative power
                        prd_k_pdkxui.Pmat_dkxu(idep,k) = prd_k_pdkxui.Pmat_dkxu(idep,k) - 1;
                            % update derivative polynomial product
                        pdkxui_Ldkxu_k = prd_k_pdkxui.Ldkxu(prd_k_pdkxui,obj.dxu);

                        % accumulate partial derivative of dkxl wrt to jet space iv
                        obj_out.Jdkxl(iv,prd_k.inds,k_) = obj_out.Jdkxl(iv,prd_k.inds,k_) + pdkxui_Ldkxu_k*pvl_;

                        if (k_<prd_.kor)
                            obj_out = prolong_vu_mvpolynomial(obj_out,k_+1,iv,pvl_,prd_k_pdkxui);
                        end
                    end
                end
            end

        end
        function obj_out = prolong_vu_mvpolynomial(obj,k_,ivp_,pvl_,prd_)
            obj_out = obj;
            prd_k = prd_;
            nvar = prd_k.nvar;

            %% first step: accumulate this input partial derivative term's contribution to dkxl
            if (ivp_>1) % if the input partial is with respect to a dependent variable
                prd_k.Pmat_dkxu(ivp_-1) = prd_k.Pmat_dkxu(ivp_-1) + 1; % increment polynomial derivative power
            end
            % evaluate input chained dkxu polynomial
            Ldkxu_k = prd_k.Ldkxu(prd_k,obj.dxu); % compute the product of derivative polynomials. Scalar valued
            obj_out.dkxl(k_,prd_k.inds) = obj_out.dkxl(k_,prd_k.inds) + Ldkxu_k*pvl_; % accumulate total der
            prd_k.nder(ivp_) = prd_k.nder(ivp_) + 1; % increment count of partial derivatives of ivp_ variable

            %% second step: accumulate this term's contribution to the Jacobian of dkxl
            % prd_.dkLtnsP
            dkLP = @(iv_,ii_,k_) prd_k.dkLtnsP(iv_,ii_,k_+1);

            % step 2a: accumulate this term's base space partial derivatives
            % identify base space lambda indices which have nonzero contribution to this order
            iipnz = ( prd_k.Pmat(:,prd_k.inds) - prd_k.nder(1:nvar) ) > 0;
            % accumulate partial derivatives over base space variables
            for iv = 1:nvar
                if ( iipnz(iv,:) )
                    % accumulate partials over the base space
                    inds_iv = prd_k.inds( iipnz(iv,:) );
                    pvl_iv = ones( 1,length(inds_iv) );
                    for iiv = 1:(iv-1)
                        pvl_iv = pvl_iv .* dkLP(iiv, inds_iv, prd_k.nder(iiv));
                    end
                    pvl_iv = pvl_iv .* dkLP(iv, inds_iv, prd_k.nder(iv) + 1);
                    for iiv = (iv+1):nvar
                        pvl_iv = pvl_iv .* dkLP(iiv, inds_iv, prd_k.nder(iiv));
                    end
                    % pvl_iv

                    % accumulate partial derivative of dkxl wrt to base space iv
                    obj_out.Jdkxl(iv,inds_iv,k_) = obj_out.Jdkxl(iv,inds_iv,k_) + Ldkxu_k*pvl_iv;

                    if (k_<prd_.kor)
                        prd_k_piv.inds = inds_iv;
                        obj_out = prolong_vu_mvpolynomial(obj_out,k_+1,iv,pvl_iv,prd_k_piv);
                    end
                end
            end

            % step 2b: accumulate this term's base space partial derivatives
            % identify jet space lambda indices which have nonzero contribution to this order
            iipdxu_nz = prd_k.Pmat_dkxu > 0
            ndep = nvar-1;
            iv = nvar;
            for k = 1:k_
                for idep = 1:ndep
                    iv = iv + 1;
                    if ( iipdxu_nz(idep,k) )
                        prd_k_pdkxui = prd_k;

                            % accumulate power rule coefficient
                        prd_k_pdkxui.c_net = prd_k_pdkxui.c_net * prd_k_pdkxui.Pmat_dkxu(idep,k);
                            % decrement derivative power
                        prd_k_pdkxui.Pmat_dkxu(idep,k) = prd_k_pdkxui.Pmat_dkxu(idep,k) - 1;
                            % update derivative polynomial product
                        pdkxui_Ldkxu_k = prd_k_pdkxui.Ldkxu(prd_k_pdkxui,obj.dxu);

                        % accumulate partial derivative of dkxl wrt to jet space iv
                        obj_out.Jdkxl(iv,prd_k.inds,k_) = obj_out.Jdkxl(iv,prd_k.inds,k_) + pdkxui_Ldkxu_k*pvl_;

                        if (k_<prd_.kor)
                            obj_out = prolong_vu_mvpolynomial(obj_out,k_+1,iv,pvl_,prd_k_pdkxui);
                        end
                    end
                end
            end
        end

        function obj_out = prolong_mvpolynomial(obj,dxu_,pr0d_)
            [ndep, nvar] = adlam.bspc_dims(obj);
            Pmat = adlam.Pmat(obj);
            Plen = size(Pmat,2);

            dxu = reshape( dxu_,ndep,[] );
            kor = size(dxu,2);
            ndim = 1+ndep*(kor+1);

            LmatP = adlam.LmatP(obj);

            dkLtnsP = nan(nvar,Plen,kor+1);
            dkLtnsP(:,:,1) = LmatP;
            for k = 1:kor
                % dkLtnsP(:,:,k) = adlam.dkLmatP(obj,k);
                % [dkLtnsP(:,:,k),dkLtnsP_full(:,:,k+1)] = deal(adlam.dkLmatP(obj,k));
                dkLtnsP(:,:,k+1) = adlam.dkLmatP(obj,k);
            end
            % dkLtnsP_short = dkLtnsP(:,:,2:end);

            pr_data = pr0d_;

            pr_data.inds = 1:Plen; % all indices contributed to zeroeth prolongation
            pr_data.kor = kor;
            pr_data.ndim = 1+ndep*(kor+1);
            pr_data.LmatP = LmatP;
            pr_data.dkLtnsP = dkLtnsP;
            pr_data.nder = zeros(ndim,1) ; % initialize derivative block

            pr_data.dxu = dxu;
            pr_data.tau_prk = @(o_,k_) [1 ; reshape( o_.dxu(:,1:k_) , [], 1 )];

            dkLP = @(iv_,k_) dkLtnsP(iv_,:,k_+1);

            pr = struct( ...
                'dkxl', zeros(kor,Plen), ...
                'Jdkxl', zeros(ndim,Plen,kor), ...
                'lkx', zeros(ndep,Plen,kor), ...
                'Jlkx', zeros(ndep,Plen,kor,prd_.ndim) ...
            );
            iv = 1;
            pr.Jdkxl(iv,:,1) = dkLP(iv,1);
            for iiv = 2:nvar


                pr.Jdkxl(iv,:,1) = pr.Jdkxl(iv,:,1) .* dkLP(iiv,0);
            end


            for iv = 2:nvar % compute gradient of each Jacobian row for flat l
                pr.Jdkxl(iv,:,1) = 1;
                for iiv = 1:(iv-1)
                    pr.Jdkxl(iv,:,1) = pr.Jdkxl(iv,:,1) .* dkLP(iiv,0);
                end
                pr.Jdkxl(iv,:,1) = pr.Jdkxl(iv,:,1) .* dkLP(iv,1);
                for iiv = (iv+1):nvar
                    pr.Jdkxl(iv,:,1) = pr.Jdkxl(iv,:,1) .* dkLP(iiv,0);
                end
            end

            % pass partal x coordinate to first prolongation
            iv = 1;
            pr = adlam.prolong_vu_mvpolynomial()

            for iv = 2:nvar % compute gradient of each Jacobian row

            end

            % dkxl(1,:) = ( pr_data.tau_prk(pr_data,1)' )*(obj.pr0.l.Jac');
            % compute Jacobian of resultant dkxl term, pass to second prolongation







            pause

            % propogate base space partial x term

            %

            % [dkxl,Jdkxl,lkx,Jlkx] = adlam.prolong_vxu_mvpolynomial(obj,pr_data);

            % pr = 0;
            % pause

            obj_out = obj;
            % obj_out.dxu = dxu;
            obj_out.pr = pr;
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
