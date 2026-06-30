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

        ntheta = @(o_) adlam.nvar(o_) * adlam.Plen(o_);

        %% handles for processed quantities of moderate complexity, still inline fast
        Lmat_full = @(o_) [ ones( adlam.nvar(o_), 1) , o_.Lmat];
        dkL_mat_full = @(o_,k_) [ zeros( adlam.nvar(o_), 1) , o_.dL_tns(:,:,k_) ];

        % functions for indexing into (nvar)x(bor+1) Lmat, yielding permutation Lmat values
        iLmatP = @(sz_,Pl_,Pm_) sub2ind(sz_,reshape((1:sz_(1))' .* ones(sz_(1),Pl_),1,[]), Pm_(:)'+1 ) ;
        LmPf_vec = @(o_,Lmf_) Lmf_(adlam.iLmatP(size(Lmf_),adlam.Plen(o_),adlam.Pmat(o_)));

        % permutation Lmat, multiplicatively accumulated down columns for lambda rows
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

        lrow_vals;
        Jl;
        pr0;

        dxu;
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
            Plen = size(Pmat,2);

            nvar = size(Pmat,1);
            ndep = nvar - 1;
            Oa = diag(A);
            xu_O = s0.val;

            [Lmat,Lmat_raw] = deal(nan(nvar,bor)); % matrix of values for each L univariate polynomial

            korp1 = kor+1; % compute L partials up to n+1 for the computation of n'th prolongation Jac
            [c_tns,dL_tns,dL_tns_raw] = deal(nan(nvar,bor,korp1)); % aggregate chain coefficient and derivatives of ordered univariate polynomial

            Lmat_raw(:,1) = xu_O; % image of xu under diagonal linear map A * [x ; u] + b into O
            c_tns(:,1,1) = Oa; % first derivative chain coefficient
            dL_tns_raw(:,1,1) = c_tns(:,1,1); % first order derivatives of linear map are just chain coefficients
            c_tns(:,1,2:end) = 0; % onward chain coefficients are zero
            dL_tns_raw(:,1,2:end) = 0; % onward dL value is zero
            for ib = 2:bor
                Lmat_raw(:,ib) = Lmat_raw(:,ib-1) .* xu_O;
                c_tns(:,ib,1) = ib*Oa;
                dL_tns_raw(:,ib,1) = c_tns(:,ib,1) .* Lmat_raw(:,ib-1);
                for ik = 2:korp1
                    if (ik<=ib)
                        c_tns(:,ib,ik) = c_tns(:,ib,ik-1) .* ((ib-ik+1)*Oa);
                        if (ik<ib)
                            dL_tns_raw(:,ib,ik) = c_tns(:,ib,ik) .* Lmat_raw(:,ib-ik);
                        else
                            dL_tns_raw(:,ib,ik) = c_tns(:,ib,ik);
                        end
                    else
                        c_tns(:,ib,ik) = 0;
                        dL_tns_raw(:,ib,ik) = 0;
                    end
                end
            end
            % use univariate values and derivatives to compute image in polynomial family space
            if (isfield(lam,'poly_coeffs'))
                [Lmat,dL_tns] = adlam.compute_Lmat_family(lam.poly_coeffs,Lmat_raw,dL_tns_raw);
            else
                [Lmat,dL_tns] = deal(Lmat_raw,dL_tns_raw);
            end


            %% core initializations
            obj.xu = xu;
            obj.lam = lam;

            obj.A =  A;
            obj.b =  b;

            obj.s0 = s0;

            obj.Lmat = Lmat;
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

            obj.lrow_vals = adlam.lrowP(obj);
            obj.Jl = JlT;
            % the base space is fully characterized by the lambda vector and its Jacobian
            pr0 = struct( ...
                'l', adobj( prod( LmatP, 1 )', JlT' ), ...
                'v', @(o_,th_) (reshape( th_, length( o_.l.val ), [] ))' * o_.l.val ...
            );
            obj.pr0 = pr0;

            %% if derivatives are provided, prolong observed jet space
            if ( size(dxu,2) )

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
                obj.dkxl = zeros(kor,Plen); % matrix of total derivatives of lambda row vector
                obj.Jdkxl = zeros(ndim,Plen,kor); % Jacobian of matrix of total derivatives of lambda row vector
                obj.lkx = zeros(ndep,Plen,kor);
                obj.Jlkx = zeros(ndep,Plen,kor,ndim);
                for iv = 1:nvar
                    obj = adlam.prolong_vu_mvpolynomial(obj,1,iv,JlT(iv,:),pr0_data);
                    for k = 1:kor
                        obj = adlam.prolong_vx_mvpolynomials(obj,k,k,iv,JlT(iv,:),pr0_data);
                    end
                end
                % obj = adlam.verify_prolongation(obj); % debugging
            else
                obj.dxu = [];
                obj.dkxl = [];
                obj.Jdkxl = [];
                obj.lkx = [];
                obj.Jlkx = [];
            end
        end
        function obj_out = prolong_jet_space(obj_in,dxu_) % for updating initially unprolonged instances
            xu = obj_in.xu(:);
            nvar = length(xu(:));
            ndep = nvar-1;

            lam = obj_in.lam;
            b =  obj_in.b;
            A =  obj_in.A;
            s0 = obj_in.s0;

            Pmat = lam.Pmat;
            bor = max(Pmat(:));
            Plen = size(Pmat,2);

            if ( (size(obj_in.dL_tns,3)-1) == length(dxu_(:))/ndep ) % nothing to do for N=1
                dL_tns = obj_in.dL_tns;
                korp1 = size(dL_tns,3);
                kor = korp1-1;
                dxu = reshape(dxu_,ndep,kor); % enforce shape to correspond with N'th jet space

                [Lmat,dL_tns] = deal(obj_in.Lmat,obj_in.dL_tns);

            else % recompute base space partial derivatives
                dxu = reshape(dxu_,ndep,[]);
                kor = size(dxu,2);
                korp1 = kor+1;

                Oa = diag(A);
                xu_O = s0.val;

                [Lmat,Lmat_raw] = nan(nvar,bor); % matrix of values for each L univariate polynomial
                [c_tns,dL_tns,dL_tns_raw] = deal(nan(nvar,bor,korp1)); % aggregate chain coefficient and derivatives of ordered univariate polynomial

                Lmat_raw(:,1) = xu_O; % image of xu under diagonal linear map A * [x ; u] + b into O
                c_tns(:,1,1) = Oa; % first derivative chain coefficient
                dL_tns_raw(:,1,1) = c_tns(:,1,1); % first order derivatives of linear map are just chain coefficients
                c_tns(:,1,2:end) = 0; % onward chain coefficients are zero
                dL_tns_raw(:,1,2:end) = 0; % onward dL value is zero
                for ib = 2:bor
                    Lmat_raw(:,ib) = Lmat_raw(:,ib-1) .* xu_O;
                    c_tns(:,ib,1) = ib*Oa;
                    dL_tns_raw(:,ib,1) = c_tns(:,ib,1) .* Lmat_raw(:,ib-1);
                    for ik = 2:korp1
                        if (ik<=ib)
                            c_tns(:,ib,ik) = c_tns(:,ib,ik-1) .* ((ib-ik+1)*Oa);
                            if (ik<ib)
                                dL_tns_raw(:,ib,ik) = c_tns(:,ib,ik) .* Lmat_raw(:,ib-ik);
                            else
                                dL_tns_raw(:,ib,ik) = c_tns(:,ib,ik);
                            end
                        else
                            c_tns(:,ib,ik) = 0;
                            dL_tns_raw(:,ib,ik) = 0;
                        end
                    end
                end
                % use univariate values and derivatives to compute image in polynomial family space
                if (isfield(lam,'poly_coeffs'))
                    [Lmat,dL_tns] = adlam.compute_Lmat_family(lam.poly_coeffs,Lmat_raw,dL_tns_raw);
                else
                    [Lmat,dL_tns] = deal(Lmat_raw,dL_tns_raw);
                end
            end

            %% core initializations
            obj_out = obj_in;

            obj_out.xu = xu;
            obj_out.lam = lam;

            obj_out.A =  A;
            obj_out.b =  b;

            obj_out.s0 = s0;

            obj_out.Lmat = Lmat;
            obj_out.dL_tns = dL_tns;

            %% prologation initializations
            LmatP = adlam.LmatP(obj_out);
            d1LmatP = adlam.dkLmatP(obj_out,1);

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

            obj_out.lrow_vals = adlam.lrowP(obj_out);
            obj_out.Jl = JlT;
            % the base space is fully characterized by the lambda vector and its Jacobian
            pr0 = struct( ...
                'l', adobj( prod( LmatP, 1 )', JlT' ), ...
                'v', @(o_,th_) (reshape( th_, length( o_.l.val ), [] ))' * o_.l.val ...
            );
            obj_out.pr0 = pr0;

            %% derivatives are provided, prolong observed jet space
            obj_out.dxu = dxu;
            ndim = 1+ndep*(kor+1);

            dkLtnsP = nan(nvar,Plen,korp1+1);
            dkLtnsP = nan(nvar,Plen,korp1+1);
            dkLtnsP(:,:,1) = LmatP;
            for k = 1:korp1
                dkLtnsP(:,:,k+1) = adlam.dkLmatP(obj_out,k);
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
            obj_out.dkxl = zeros(kor,Plen); % matrix of total derivatives of lambda row vector
            obj_out.Jdkxl = zeros(ndim,Plen,kor); % Jacobian of matrix of total derivatives of lambda row vector
            obj_out.lkx = zeros(ndep,Plen,kor);
            obj_out.Jlkx = zeros(ndep,Plen,kor,ndim);
            for iv = 1:nvar
                obj_out = adlam.prolong_vu_mvpolynomial(obj_out,1,iv,JlT(iv,:),pr0_data);
                for k = 1:kor
                    obj_out = adlam.prolong_vx_mvpolynomials(obj_out,k,k,iv,JlT(iv,:),pr0_data);
                end
            end
            % obj_out = adlam.verify_prolongation(obj_out); % debugging

        end
        function J_out = J_tau_dNm1xu(obj,th_,tdNm1xu_)
            tdNm1xu = tdNm1xu_(:);
            nvar = length(obj.xu(:));
            ndep = nvar-1;
            kor = max([1,size(obj.dxu,2)]);
            ndim = 1 + ndep*(kor+1);

            Th = reshape(th_,[],nvar);
            th_x = Th(:,1);
            Th_u = Th(:,2:end);

            if (kor==1)
                J_out = [(obj.Jl*( Th_u - th_x * tdNm1xu'))',zeros(ndep,ndim-nvar)];
            else
                Plen = length(th_x(:));
                J_out = ( Th_u' * obj.Jdkxl(:,:,kor-1)' ) ...
                    - (reshape(pagemtimes(reshape(obj.Jlkx(:,:,kor-1,:),ndep,Plen,ndim),th_x),ndep,ndim) ...
                        + tdNm1xu*([ obj.Jl ; zeros(ndim-nvar,Plen) ]*th_x)');
            end
        end
        function [J_out,tau_uN_out] = J_tau_uN(obj,th_)
            lv = obj.lrow_vals; % 1 by Plen
            dkxlv = obj.dkxl; % kor by Plen
            lkx = obj.lkx; % ndep by Plen by kor
            Jlv0 = obj.Jl; % nvar by Plen
            Jdkxlv = obj.Jdkxl; % ndim by Plen by kor
            Jlkx = obj.Jlkx; % ndep by Plen by kor by ndim

            [ndep,Plen,kor,ndim]  = size(Jlkx);
            nvar = ndep+1;
            korp1 = kor+1;

            Jlv = [ Jlv0 ; zeros(ndim-nvar,Plen) ];

            Th = reshape(th_,Plen,nvar);
            th_x = Th(:,1);
            Th_u = Th(:,2:end);

            J_out = zeros(ndep,ndim,korp1);
            tau_uN_out = nan(ndep,korp1); % free (and necessary) evaluation of tau_uN

            thx_comp = @(J_) reshape(pagemtimes(reshape(J_,ndep,Plen,ndim),th_x),ndep,ndim);

            % base space
            tau_uN_out(:,1) = ( lv * Th_u )'; % tau_ukm1 = dxu = tau_u
            J_out(:,1:nvar,1) = (Jlv0*( Th_u - th_x * tau_uN_out(:,1)' ))'; % J tau_u = J dxu
            gl_thx = (Jlv*th_x)';
            for k = 2:kor % first thru N-1'th jet space
                tau_uN_out(:,k) = ( dkxlv(k-1,:) * Th_u )' ...
                    - lkx(:,:,k-1)*th_x; % tau_ukm1 = dkxu
                J_out(:,:,k) = ( Th_u' * Jdkxlv(:,:,k-1)' ) ...
                    - ( thx_comp( Jlkx(:,:,k-1,:) ) + tau_uN_out(:,k)*gl_thx ); % J dkxu
            end
            % N'th jet space
            tau_uN_out(:,korp1) = ( dkxlv(kor,:) * Th_u )' ...
                - lkx(:,:,kor)*th_x; % tau_uN = dNp1xu
            J_out(:,:,korp1) = ( Th_u' * Jdkxlv(:,:,kor)' ) ...
                - ( thx_comp( Jlkx(:,:,kor,:) ) + tau_uN_out(:,korp1)*gl_thx ); % J dkxu
        end
        function Rtns_out = Renc_tns(obj,ords_)
            Plen = adlam.Plen(obj);
            [ndep,kor] = deal( adlam.ndep(obj),adlam.kor(obj) );
            %% initialize injection into Lambda column space
            i_imm = zeros(ndep*Plen,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            function l_imm = immerse_lambda(l_)
                l_imm = zeros((ndep*Plen)*ndep,1);
                l_imm(logical( i_imm(:) )) = reshape( l_(:) * ones(1,ndep), [], 1);
                l_imm = (reshape(l_imm,ndep*Plen,ndep))';
                % l_ -> [l_ 0 ... 0 ; 0 l_ ... 0 ; ...]
            end

            nords = length(ords_(:));
            ntheta = (1+ndep)*Plen;

            Rtns_out = nan(ndep,ntheta,nords);
            for iord = 1:nords
                if ( ords_(iord)==1 )
                    Rtns_out(:,:,iord) = [ -obj.dxu(:,1)*obj.lrow_vals, immerse_lambda(obj.lrow_vals) ];
                else
                    k_i = ords_(iord);
                    if ( k_i <= kor )
                        Rtns_out(:,:,iord) = [  ...
                        -obj.lkx(:,:,k_i-1)-obj.dxu(:,k_i)*obj.lrow_vals, ...
                        immerse_lambda(obj.dkxl(k_i-1,:)) ...
                        ];
                    else
                        Rtns_out(:,:,iord) = zeros(ndep,theta);
                    end
                end
            end
        end
        function obj_out = prn_dxu(obj,dxu_)
            obj_out = adlam(obj.lam,obj.xu,dxu_);
        end
        function i_imm_out = i_imm_data(obj)
            ndep = adlam.ndep(obj);
            Plen = adlam.Plen(obj);
            ntheta = (1+ndep)*Plen;
            ntheta_u = ntheta - Plen;
            i_imm = zeros(ntheta_u,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            i_imm_out = struct( ...
            'ndep', ndep, ...
            'Plen', ntheta, ...
            'ntheta_u', ntheta_u, ...
            'i_imm', i_imm ...
            );
        end
        function Lambda_uk_out = Lambda_u_k(obj,k_)
            % i_imm_dat = obj.i_imm_data();
            % ndep = i_imm_dat.ndep;
            % Plen = i_imm_dat.Plen;
            % ntheta = i_imm_dat.ntheta;
            % ntheta_u = i_imm_dat.ntheta_u;
            % i_imm = i_imm_dat.i_imm;

            ndep = adlam.ndep(obj);
            Plen = adlam.Plen(obj);
            %% initialize injection into Lambda column space
            i_imm = zeros(ndep*Plen,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            function l_imm = immerse_lambda(l_)
                l_imm = zeros((ndep*Plen)*ndep,1);
                % l_imm(logical( i_imm(:) )) = reshape( (ones(ndep,1)*l_)', [], 1);
                l_imm(logical( i_imm(:) )) = reshape( l_(:) * ones(1,ndep), [], 1);
                % l_imm = reshape(l_imm,ndep,ndep*Plen);
                l_imm = (reshape(l_imm,ndep*Plen,ndep))';
                % l_ -> [l_ 0 ... 0 ; 0 l_ ... 0 ; ...]
            end

            if (k==0)
                Lambda_uk_out = [ zeros(ndep,Plen) , l_imm( adlam.lrowP(obj) ) ];
            else
                Lambda_uk_out = [ obj.lkx(:,:,k_) , l_imm( obj.dkxl(k_,:) ) ];
            end
        end
    end

    methods (Static)
        function fspace_out = init_fspace_family(fspace_)
            nvar = size(fspace_.Pmat,1);
            bor = fspace_.bor;
            borp1 = bor+1;

            poly_coeffs = zeros(borp1,borp1);
            function c_out = get_coeff(j_,i_)
                if ( (0<=j_)&&(j_<=i_) )
                    c_out = poly_coeffs(j_+1,i_+1);
                else
                    c_out = 0.0;
                end
            end
            function comp_Legendre_coeffs()
                poly_coeffs(1,1) = 1.0;
                if (bor > 0)
                    [poly_coeffs(1,2),poly_coeffs(2,2)] = deal(0.0,1.0);
                    for ib = 2:bor
                        k = ib;
                        k2_m1 = 2*k - 1;
                        k_m1 = k - 1;
                        for jb = 0:ib
                            poly_coeffs(jb+1,ib+1) = ( k2_m1*(get_coeff(jb-1,ib-1)) - k_m1*(get_coeff(jb,ib-2)) )/k;
                        end
                    end
                end
            end
            function comp_Hermite1_coeffs() % physicist's Hermite polynomials
                poly_coeffs(1,1) = 1.0;
                if (bor > 0)
                    [poly_coeffs(1,2),poly_coeffs(2,2)] = deal(0.0,2.0);
                    for ib = 2:bor
                        poly_coeffs(1,ib+1) = -poly_coeffs(2,ib);
                        for jb = 1:ib
                            poly_coeffs(jb+1,ib+1) = 2.0*get_coeff(jb-1,ib-1) - (jb+1)*get_coeff(jb+1,ib-1);
                        end
                    end
                end
            end
            function comp_Hermite2_coeffs() % probabilist's Hermite polynomials
                poly_coeffs(1,1) = 1.0;
                if (bor > 0)
                    [poly_coeffs(1,2),poly_coeffs(2,2)] = deal(0.0,1.0);
                    for ib = 2:bor
                        poly_coeffs(1,ib+1) = -poly_coeffs(2,ib);
                        for jb = 1:ib
                            poly_coeffs(jb+1,ib+1) = get_coeff(jb-1,ib-1) - (jb+1)*get_coeff(jb+1,ib-1);
                        end
                    end
                end
            end
            function comp_Chebyshev_coeffs()
                for ib = 2:bor
                    poly_coeffs(1,ib+1) = -poly_coeffs(1,ib-1);
                    for jb = 1:ib
                        poly_coeffs(jb+1,ib+1) = 2.0*get_coeff(jb-1,ib-1) - get_coeff(jb,ib-2);
                    end
                end
            end
            function comp_Chebyshev1_coeffs()
                poly_coeffs(1,1) = 1.0;
                if (bor > 0)
                    [poly_coeffs(1,2),poly_coeffs(2,2)] = deal(0.0,1.0);
                    comp_Chebyshev_coeffs()
                end
            end
            function comp_Chebyshev2_coeffs()
                poly_coeffs(1,1) = 1.0;
                if (bor > 0)
                    [poly_coeffs(1,2),poly_coeffs(2,2)] = deal(0.0,2.0);
                    comp_Chebyshev_coeffs()
                end
            end

            fspace_out = fspace_;
            if (isfield(fspace_,'fam'))
                switch (fspace_.fam)
                    case 'Legendre'
                        comp_Legendre_coeffs();
                        fspace_out.poly_coeffs = poly_coeffs;
                    case 'Hermite1'
                        comp_Hermite1_coeffs();
                        fspace_out.poly_coeffs = poly_coeffs;
                    case 'Hermite2'
                        comp_Hermite2_coeffs();
                        fspace_out.poly_coeffs = poly_coeffs;
                    case 'Chebyshev1'
                        comp_Chebyshev1_coeffs();
                        fspace_out.poly_coeffs = poly_coeffs;
                    case 'Chebyshev2'
                        comp_Chebyshev2_coeffs();
                        fspace_out.poly_coeffs = poly_coeffs;
                end
            end
        end
        function [Lmat_fam, dL_tns_fam] = compute_Lmat_family(poly_coeffs_,Lmat_,dL_tns_)
            [nvar,bor] = size(Lmat_);
            korp1 = size(dL_tns_,3);
            borp1 = bor+1;
            Lmat_all = [ ones(nvar,1) , Lmat_ ];

            Lmat_fam = nan(size(Lmat_)); % only need to compute (1+ndep)x(bor) terms
            dL_tns_fam = nan(size(dL_tns_));

            Lmat_fam(:,1) = Lmat_all(:,1)*poly_coeffs_(1,1);
            dL_tns_fam(:,1,1) = dL_tns_(:,1,1);
            dL_tns_fam(:,1,2:end) = 0;
            for ib = 2:bor
                iib = 1:(ib+1);
                Lmat_fam(:,ib) = Lmat_all(:,iib)*poly_coeffs_(iib,ib);
                iiib = 1:ib;
                pcoeff_ib = poly_coeffs_(2:(ib+1),ib);
                for ik = 1:korp1
                    dL_tns_fam(:,ib,ik) = dL_tns_(:,iiib,ik)*pcoeff_ib;
                end
            end
        end
        function obj_out = verify_prolongation(obj)

            [ndep,kor] = deal( adlam.ndep(obj),adlam.kor(obj) );
            Plen = adlam.Plen(obj);

            nvar = ndep+1;
            ndim = nvar + kor*ndep;

            xu = obj.xu;
            dxu = obj.dxu;

            tkderset = @(k_) [ 1 ; reshape( dxu(:,1:(k_+1)), [],1 ) ];
            tk_N = @(tk_) [ tk_ ; zeros(ndim-length(tk_),1) ];
            tk_Nthspace = @(k_) tk_N(tkderset(k_));

            %% check adherence to 2 primary prolongation verification quantities

            dkp1x_l = nan(kor,Plen);
            lkp1x = nan(ndep,Plen,kor);

            % 1) verify dkp1x l = tau^(k) . grad^(k) dkx l
            dkp1x_l(1,:) = ( tkderset(0) )' * obj.Jl;
            % 2) verify lkp1x l = tau^(k) . grad^(k) lkp1x + dkp1x u . dxl
            lkp1x(:,:,1) = dxu(:,1) * obj.dkxl(1,:);
            for k = 2:kor
                iJk = nvar + ((k-1)*ndep);
                dkp1x_l(k,:) = ( tk_Nthspace(k-1) )' * obj.Jdkxl(:,:,k-1);
                for i = 1:ndep
                    lkp1x(i,:,k) = dxu(i,k) * obj.dkxl(1,:) ...
                                   + ( tk_Nthspace(kor-1) )' * reshape( obj.Jlkx(i,:,k-1,:), Plen, ndim )';
                end
            end

            dkp1x_l
            obj.dkxl

            lkp1x
            obj.lkx

            Pmat = adlam.Pmat(obj)
            d0xl = obj.lrow_vals
            LmatP = adlam.LmatP(obj)
            d1LmatP = adlam.dkLmatP(obj,1)
            Lmat = obj.Lmat
            dL_tns = obj.dL_tns

            function dkxli_out = dxli_evl(i_)
                Pvals0 = Pmat(1:2,i_);
                [Lx,Lu] = deal(LmatP(1,i_), LmatP(2,i_));
                if (Pvals0(1) > 1)
                    dxLx = Pvals0(1)*Lmat(1,Pvals0(1)-1);
                    if (Pvals0(1) > 2)
                        d2xLx = Pvals0(1)*(Pvals0(1)-1)*Lmat(1,Pvals0(1)-2);
                    else
                        d2xLx = Pvals0(1);
                    end
                else
                    dxLx = Pvals0(1);
                    d2xLx = 0.0;
                end
                if (Pvals0(2) > 1)
                    dxLu = Pvals0(2)*Lmat(2,Pvals0(2)-1)*dxu(1,1);
                    if (Pvals0(2) > 2)
                        d2xLu = Pvals0(2)*(Pvals0(2)-1)*Lmat(2,Pvals0(2)-2)*dxu(1,1)*dxu(1,1) ...
                                + Pvals0(2)*Lmat(2,Pvals0(2)-1)*dxu(1,2);
                    else
                        d2xLu = Pvals0(2)*(Pvals0(2)-1)*dxu(1,1)*dxu(1,1) ...
                                + Pvals0(2)*Lmat(2,Pvals0(2)-1)*dxu(1,2);
                    end
                else
                    dxLu = Pvals0(2)*dxu(1,1);
                    d2xLu = Pvals0(2)*dxu(1,2);
                end
                % d1xli_out = [ Pvals0(1)*LmatP(2,i_) , Pvals0(2)*LmatP(1,i_)*dxu(1,1) ] * [ iele1 ; iele2 ];
                dkxli_out = [ dxLx*Lu + Lx*dxLu ; d2xLx*Lu + 2*(dxLx*dxLu) + Lx*d2xLu ];
            end

            dkxli = nan(2,Plen);
            likx = nan(ndep,Plen,2);
            for iP = 1:Plen
                dxli(:,iP) = dxli_evl(iP);
                likx(:,iP,1) = dxu(:,1) * dxli(1,iP);
                likx(:,iP,2) = ( dxu(:,1) * dxli(2,iP) ) + 2*dxu(:,2) * dxli(1,iP);
            end
            check1 = [ ...
            Pmat ; ...
            dxli ; ...
            dkp1x_l; ...
            obj.dkxl ...
            ]

            check2 = [ ...
            reshape(permute(likx(1,:,1:2),[3 2 1]), 2, Plen) ; ...
            reshape(permute(lkp1x(1,:,1:2),[3 2 1]), 2, Plen) ; ...
            reshape(permute(obj.lkx(1,:,1:2),[3 2 1]), 2, Plen)  ...
            ]

            pause

            obj_out = obj;

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
            dkLP = @(iv_,ii_,k_) prd_k.dkLtnsP(iv_,ii_,k_+1);

            % step 2a: accumulate this term's base space partial derivatives
            % identify base space lambda indices which have nonzero contribution to this order
            iipnz = ( prd_k.Pmat(:,prd_k.inds) - prd_k.nder(1:nvar) ) > 0;
            % accumulate partial derivatives over base space variables
            for iv = 1:nvar
                if ( sum(iipnz(iv,:)) )
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

                    % accumulate partial derivative of dkxl wrt to base space iv
                    obj_out.Jdkxl(iv,inds_iv,k_) = obj_out.Jdkxl(iv,inds_iv,k_) + Ldkxu_k*pvl_iv;

                    if (k_<prd_.kor)
                        prd_k_piv = prd_k;
                        prd_k_piv.inds = inds_iv;
                        obj_out = adlam.prolong_vu_mvpolynomial(obj_out,k_+1,iv,pvl_iv,prd_k_piv);
                    end
                end
            end

            % step 2b: accumulate this term's jet space partial derivatives
            % identify jet space lambda indices which have nonzero contribution to this order
            iipdxu_nz = prd_k.Pmat_dkxu > 0;
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
                        obj_out.Jdkxl(iv,prd_.inds,k_) = obj_out.Jdkxl(iv,prd_.inds,k_) + pdkxui_Ldkxu_k*pvl_;

                        if (k_<prd_.kor)
                            obj_out = adlam.prolong_vu_mvpolynomial(obj_out,k_+1,iv,pvl_,prd_k_pdkxui);
                        end
                    end
                end
            end
        end
        function obj_out = prolong_vx_mvpolynomials(obj,k_,src_,ivp_,pvl_,prd_)
            obj_out = obj;
            prd_k = prd_;
            nvar = prd_k.nvar;

            %% first step: accumulate input partial derivative term's contribution to lkx.
            % step 1a: accumulate input partial derivative term's contribution to lkx
            if (ivp_>1) % if the input partial is with respect to a dependent variable
                prd_k.Pmat_dkxu(ivp_-1) = prd_k.Pmat_dkxu(ivp_-1) + 1; % increment polynomial derivative power
            end
            % evaluate input chained dkxu polynomial
            Ldkxu_k0 = prd_k.Ldkxu(prd_k,obj.dxu); % product of derivative polynomials. Scalar
            Ldkxu_k = obj.dxu(:,src_)*Ldkxu_k0; % net product of derivative polynomials times src derivatives. Vector.

            % acc tot der of lkm1x into lkx
            obj_out.lkx(:,prd_k.inds,k_) = obj_out.lkx(:,prd_k.inds,k_) + Ldkxu_k*pvl_;

            prd_k.nder(ivp_) = prd_k.nder(ivp_) + 1; % increment count of partial derivatives of ivp_ variable

            % step 1b: accumulate first half of product rule to lkx Jacobian
            inds_src = prd_.dkxu_inds(prd_,src_);
            % obj_out.Jlkx( :, prd_.inds, k_, inds_src ) = obj_out.Jlkx( :, prd_.inds, k_, inds_src ) ...
            %                                             + ones(prd_.ndep,1)*Ldkxu_k0*pvl_;
            ider_src = (nvar + ( 1:prd_.ndep )) + (src_-1)*prd_.ndep;
            for idep = 1:prd_.ndep
                obj_out.Jlkx(idep,prd_.inds,k_,ider_src(idep)) = obj_out.Jlkx(idep,prd_.inds,k_,ider_src(idep)) ...
                    + Ldkxu_k0*pvl_;
            end
            if ( k_ < prd_.kor ) % pass identical gradient information to next prolongation. Source incremented by one.
                obj_out = adlam.prolong_vx_mvpolynomials(obj_out,k_+1,src_+1,ivp_,pvl_,prd_);
            end

            %% second step: accumulate input's contribution to the Jacobian of lkx
            dkLP = @(iv_,ii_,k_) prd_k.dkLtnsP(iv_,ii_,k_+1);
            % step 2a: accumulate this term's base space partial derivatives
            % identify base space lambda indices which have nonzero contribution to this order
            iipnz = ( prd_k.Pmat(:,prd_k.inds) - prd_k.nder(1:nvar) ) > 0;
            % accumulate partial derivatives over base space variables
            for iv = 1:nvar
                if ( sum(iipnz(iv,:)) )
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
                    % accumulate partial derivative of dkxl wrt to base space iv
                    obj_out.Jlkx(:,inds_iv,k_,iv) = obj_out.Jlkx(:,inds_iv,k_,iv) + Ldkxu_k*pvl_iv;
                    if (k_<prd_.kor)
                        prd_k_piv = prd_k;
                        prd_k_piv.inds = inds_iv;
                        obj_out = adlam.prolong_vx_mvpolynomials(obj_out,k_+1,src_,iv,pvl_iv,prd_k_piv);
                    end
                end
            end

            % step 2b: accumulate this term's jet space partial derivatives
            % identify jet space lambda indices which have nonzero contribution to this order
            iipdxu_nz = prd_k.Pmat_dkxu > 0;
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
                        % accumulate partial derivative of lkx wrt to jet space (derivative) iv
                        obj_out.Jlkx(:,prd_.inds,k_,iv) = obj_out.Jlkx(:,prd_.inds,k_,iv) ...
                                                        + obj.dxu(:,src_)*prd_k_pdkxui.Ldkxu(prd_k_pdkxui,obj.dxu)*pvl_;
                        if (k_<prd_.kor)
                            obj_out = adlam.prolong_vx_mvpolynomials(obj_out,k_+1,src_,iv,pvl_,prd_k_pdkxui);
                        end
                    end
                end
            end
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
