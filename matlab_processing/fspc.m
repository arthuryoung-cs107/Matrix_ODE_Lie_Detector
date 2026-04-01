classdef fspc < jspc
    properties (Constant)

        bor_cap = @(obj_) (obj_.bor)*(obj_.ndep+1);

    end
    properties

        bor;

    end
    methods

        function obj = fspc(obs_)

            obj@jspc(obs_.ndep,obs_.eor);

        end

    end

    methods (Static)
        function [Bout_,dim_out,s_out] = append_to_basis(Bin_,Uin_,dim_scl_)
            if (nargin == 3)
                dim_scl = dim_scl_;
            else
                dim_scl = max(size(Uin_));
            end

            % the column space of the resultant matrix are the new vector fields
            [r_add,s_out,~,U_new] = ...
                fspc.rsV_unpack( Uin_ - (Bin_)*(Bin_')*(Uin_), dim_scl );

            r_add =  sum(double( s_out==1 ));
            Bout_ = [Bin_, U_new(:,1:r_add) ];
            dim_out = size(Bout_,2);
        end
        function str_out = svd_methods_package()
            str_out = struct( ...
            'k', @(svd_) svd_.ncol-svd_.r , ...
            'D', @(svd_) svd_.V(:,1:(svd_.r)) , ...
            'K', @(svd_) svd_.V(:,(svd_.r+1):end) , ...
            'PD', @(svd_) [ svd_.V(:,1:(svd_.r)) , zeros(svd_.ncol,svd_.ncol-svd_.r) ], ...
            'PK', @(svd_) [ zeros(svd_.ncol,svd_.r) , svd_.V(:,(svd_.r+1):end) ], ...
            'PK_k', @(svd_,k_) [ zeros(svd_.ncol,svd_.ncol-k_) , svd_.V(:,(end-(k_-1)):end) ], ...
            'PK_r', @(svd_,r_) [ zeros(svd_.ncol,r_) , svd_.V(:,(r_+1):end) ] ...
            );
        end
        function str_out = compute_svd_package(A_,Uflag_)
            if (nargin==2)
                Uflag = Uflag_;
            else
                Uflag = false;
            end

            if (Uflag)
                [W_A,r_A,s_A,V_A,Amat_,U_A] = fspc.safely_process_net_svd(A_);
            else
                [W_A,r_A,s_A,V_A,Amat_] = fspc.safely_process_net_svd(A_);
                U_A = [];
            end

            str_out = fspc.svd_methods_package();

            str_out.mat = Amat_;
            str_out.ncol = length( s_A );
            str_out.W = W_A;
            str_out.r = r_A;
            str_out.s = s_A;
            str_out.V = V_A;
            str_out.U = U_A;
            str_out.Y =  V_A.*( (s_A(:)/s_A(1))' );

            % A_comp = A_;
            % [W_A,r_A,s_A,V_A,Amat_] = fspc.safely_process_net_svd(A_comp);
            % str_out = struct( ...
            %     'mat', Amat_, ...
            %     'ncol', length( s_A ), ...
            %     'W', W_A, ...
            %     'r', r_A, ...
            %     's', s_A, ...
            %     'V', V_A, ...
            %     'Y',  V_A.*( (s_A(:)/s_A(1))' ),  ...
            %     'D', @(svd_) [ svd_.V(:,1:(svd_.r)) , zeros(svd_.ncol,svd_.ncol-svd_.r) ], ...
            %     'PK', @(svd_) [ zeros(svd_.ncol,svd_.r) , svd_.V(:,(svd_.r+1):end) ], ...
            %     'PK_k', @(svd_,k_) [ zeros(svd_.ncol,svd_.ncol-k_) , svd_.V(:,(end-(k_-1)):end) ], ...
            %     'PK_r', @(svd_,r_) [ zeros(svd_.ncol,r_) , svd_.V(:,(r_+1):end) ] ...
            % );
        end

        function [s_scaled_pseudo_inverse,s_scl] = s_spi(s_)

            % s_scaled_pseudo_inverse = reshape(1.0./(s_/s_(end)),1,[]);
            ss = s_;
            logc_zero = (s_ == 0.0);
            [s_scl,imin_nzero] = min( s_(~logc_zero) );
            ss((imin_nzero+1):end) = s_scl;

            s_scaled_pseudo_inverse = reshape(1.0./(ss/s_scl),1,[]);

        end

        function [WA_,rA_,sA_,VA_,Afull_,UA_] = safely_process_net_svd(Atns_)
            size_A = size(Atns_);
            if (length(size_A) > 2)
                Afull_ = ( reshape( Atns_, size(Atns_,1),[] ) )';
            else
                if ( size_A(1)<size_A(2) ) % flip if short matrix
                    Afull_ = Atns_';
                else
                    Afull_ = Atns_;
                end
            end

            if (nargout == 6)
                [rA_,sA_,VA_,UA_] = fspc.rsV_unpack( Afull_ );
            else
                [rA_,sA_,VA_] = fspc.rsV_unpack( Afull_ );
            end

            ssA = sA_;
            logc_zero = (sA_ == 0.0);
            [s_min_nzero,imin_nzero] = min( sA_(~logc_zero) );
            ssA((imin_nzero+1):end) = s_min_nzero;

            WA_ = VA_.*reshape(1.0./(ssA/s_min_nzero),1,[]);

        end

        function str_out = rsV_struct(A_)
            [r_,s_,V_] = fspc.rsV_unpack(A_);
            str_out = struct( ...
                'r', r_, ...
                's', s_, ...
                'V', V_ ...
            );
        end

        function [r_,s_,V_,U_] = rsV_unpack(A_, dim_scl_)
            if (nargin == 2)
                dim_scl = dim_scl_;
            else
                dim_scl = max(size(A_));
            end
            if (nargout == 4)
                [U_,S_,V_] = svd(A_,'econ');
            else
                [~,S_,V_] = svd(A_,'econ');
            end
            s_ = reshape(diag(S_),[],1); % output s_ is col vector
            r_ = sum(double(s_ > dim_scl*eps(S_(1)))); % default matlab tol
        end

        function Lam_out = comp_prk_Lambda(l_,Lx_,Lu_)
            [ndep,Plen,kor,nobs] = size(Lx_);
            nvar = ndep+1;
            ndim = 1 + ndep*(kor+1);
            ntheta = Plen*nvar;

            ncol_Lambda_u = ntheta-Plen;
            i_imm = zeros(ncol_Lambda_u,ndep);
            for i = 1:ndep
                idel = (i-1)*Plen;
                i_imm( (1+idel):(Plen+idel), i ) = 1;
            end
            i_imm_vec = logical( i_imm(:) );
            len_Lambda_u = prod(size(i_imm));
            l_imm_init = zeros(len_Lambda_u,1);
            ones_imm = ones(ndep,1);
            function l_imm = immerse_lambda(l_)
                l_imm = l_imm_init;
                l_imm(i_imm_vec) = reshape( (ones_imm.*(l_(:)'))', len_Lambda_u ,1);
                l_imm = reshape(l_imm,ndep,ncol_Lambda_u);
                % l_ -> [l_ 0 ... 0 ; 0 l_ ... 0 ; ...]
            end

            zero_Lx_block = zeros(ndep, Plen);
            zero_Lu_block = zeros(ndep, len_Lambda_u);

            Lam_out = zeros(ndim,ntheta,nobs);
            % set 0th and 1st jet space blocks of Lambda
            for i = 1:nobs
                Lam_out(1:(nvar+ndep),:,i) = [ ...
                l_(i,:) , zero_Lu_block ; ...
                zero_Lx_block , immerse_lambda(l_(i,:)) ; ...
                reshape(Lx_(:,:,1,i), ndep, Plen) , immerse_lambda(Lu_(:,1,i)) ; ...
                ];
            end
            if (kor>1)
                % To do
            end
        end

        function print_vshort_polynomial_theta_z(theta_z_,Pmat_,Zname_,preamble_)
            if (nargin==2)
                Zname = 'z';
                preamble = '\n';
            else
                Zname = Zname_;
                if (nargin==3)
                    preamble = '\n';
                else
                    preamble = preamble_;
                end
            end
            ndep_ = size(Pmat_,1)-1;
            small_tol = 1e-6;

            fprintf(preamble);
            fprintf('multivariate (1+Q = %d) polynomial of form\n',1+ndep_);

            Plen = length(theta_z_);

            ord_z = 0;
            theta_z = theta_z_;
            theta_z_maxmag = ( max(abs( theta_z )) );
            theta_z_print = theta_z/theta_z_maxmag;
            fprintf('   v_%s (x,u) = (%.1e)( ',Zname,theta_z_maxmag);
            for i = 1:Plen
                if ( abs(theta_z_print(i)) > small_tol )
                    ord_z = max([ord_z, sum(Pmat_(:,i))]);
                    fprintf('+ (%.2f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                    for idep = 1:ndep_
                        fprintf(',%d',Pmat_(idep+1,i) );
                    end
                    fprintf('] ');
                end
            end
            fprintf('), %d/%d largest terms (ord=%d) \n', ...
            sum(double( abs(theta_z_print)>small_tol )),Plen,ord_z);

            fprintf('(small_tol = %.1e)\n', small_tol);

        end
        function print_vshort_polynomial_theta(theta_,Pmat_,preamble_)
            if (nargin==2)
                preamble = '\n';
            else
                preamble = preamble_;
            end
            ndep_ = size(Pmat_,1)-1;
            small_tol = 1e-6;

            fprintf(preamble);
            fprintf('multivariate (1+Q = %d) polynomials of form\n',1+ndep_);

            theta_maxmag = ( max(abs(theta_(:))) );
            theta_print = theta_/theta_maxmag;

            ntheta = length(theta_);
            Plen = ntheta/(1+ndep_);

            theta_mat = reshape(theta_,[],1+ndep_);

            ord_z = 0;
            theta_z = theta_mat(:,1);
            theta_z_maxmag = ( max(abs( theta_z )) );
            theta_z_print = theta_z/theta_z_maxmag;
            fprintf('   v_x (x,u) = (%.1e)( ',theta_z_maxmag);
            for i = 1:Plen
                if ( abs(theta_z_print(i)) > small_tol )
                    ord_z = max([ord_z, sum(Pmat_(:,i))]);
                    fprintf('+ (%.2f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                    for idep = 1:ndep_
                        fprintf(',%d',Pmat_(idep+1,i) );
                    end
                    fprintf('] ');
                end
            end
            fprintf('), %d/%d largest terms (ord=%d) \n', ...
            sum(double( abs(theta_z_print)>small_tol )),Plen,ord_z);

            for jdep = 1:ndep_

                ord_z = 0;
                theta_z = theta_mat(:,jdep+1);
                theta_z_maxmag = ( max(abs( theta_z )) );
                theta_z_print = theta_z/theta_z_maxmag;

                fprintf('   v_u_%d (x,u) = (%.1e)( ',jdep,theta_z_maxmag);
                for i = 1:Plen
                    if ( abs(theta_z_print(i)) > small_tol )
                        ord_z = max([ord_z, sum(Pmat_(:,i))]);
                        fprintf('+ (%.2f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                        for idep = 1:ndep_
                            fprintf(',%d',Pmat_(idep+1,i) );
                        end
                        fprintf('] ');
                    end
                end
                fprintf('), %d/%d largest terms (ord=%d) \n', ...
                sum(double( abs(theta_z_print)>small_tol )),Plen,ord_z);
            end
            fprintf('(small_tol = %.1e)\n', small_tol);

        end
        function print_short_polynomial_theta_z(theta_z_,Pmat_,Zname_,preamble_)
            if (nargin==2)
                Zname = 'z';
                preamble = '\n';
            else
                Zname = Zname_;
                if (nargin==3)
                    preamble = '\n';
                else
                    preamble = preamble_;
                end
            end
            ndep_ = size(Pmat_,1)-1;
            small_tol = 1e-6;

            fprintf(preamble);
            fprintf('multivariate (1+Q = %d) polynomial of form\n',1+ndep_);
            fprintf('   v_%s(xu) = Sum_i theta(i,z)*l_i(xu) ',Zname);
            fprintf(' = Sum_i theta(i,z)*(xu)^[k(i,xu)] = Sum_i theta(i,z)*(x^k(i,x)');
            for i = 1:ndep_
                fprintf('*u_%d^k(i,u_%d)',i,i);
            end
            fprintf(')\n');
            fprintf('   where theta(i,z) is coeff of mv polynomial l_i(xu)=(xu)^[k(i,xu)], combining into v_z : z = x, u_1, ..., u_Q (Q=%d)\n',ndep_);

            Plen = length(theta_z_);

            theta_z = theta_z_;
            theta_z_maxmag = ( max(abs( theta_z )) );
            theta_z_print = theta_z/theta_z_maxmag;
            fprintf('   v_%s (x,u) = (%.1e)( ',Zname,theta_z_maxmag);
            for i = 1:Plen
                if ( abs(theta_z_print(i)) > small_tol )
                    fprintf('+ (%.2f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                    for idep = 1:ndep_
                        fprintf(',%d',Pmat_(idep+1,i) );
                    end
                    fprintf('] ');
                end
            end
            fprintf('), %d/%d largest terms\n', ...
            sum(double( abs(theta_z_print)>small_tol )),Plen );

            fprintf('(small_tol = %.1e)\n', small_tol);

        end

        function print_short_polynomial_theta(theta_,Pmat_,preamble_)
            if (nargin==2)
                preamble = '\n';
            else
                preamble = preamble_;
            end
            ndep_ = size(Pmat_,1)-1;
            small_tol = 1e-6;

            fprintf(preamble);
            fprintf('multivariate (1+Q = %d) polynomials of form\n',1+ndep_);
            fprintf('   v_z(xu) = Sum_i theta(i,z)*l_i(xu) ');
            fprintf(' = Sum_i theta(i,z)*(xu)^[k(i,xu)] = Sum_i theta(i,z)*(x^k(i,x)');
            for i = 1:ndep_
                fprintf('*u_%d^k(i,u_%d)',i,i);
            end
            fprintf(')\n');
            fprintf('   where theta(i,z) is coeff of mv polynomial l_i(xu)=(xu)^[k(i,xu)], combining into v_z : z = x, u_1, ..., u_Q (Q=%d)\n',ndep_);

            theta_maxmag = ( max(abs(theta_(:))) );
            theta_print = theta_/theta_maxmag;

            ntheta = length(theta_);
            Plen = ntheta/(1+ndep_);

            theta_mat = reshape(theta_,[],1+ndep_);

            theta_z = theta_mat(:,1);
            theta_z_maxmag = ( max(abs( theta_z )) );
            theta_z_print = theta_z/theta_z_maxmag;
            fprintf('   v_x (x,u) = (%.1e)( ',theta_z_maxmag);
            for i = 1:Plen
                if ( abs(theta_z_print(i)) > small_tol )
                    fprintf('+ (%.2f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                    for idep = 1:ndep_
                        fprintf(',%d',Pmat_(idep+1,i) );
                    end
                    fprintf('] ');
                end
            end
            % fprintf('), %d/%d largest terms (tol=%1.e) \n', ...
            % sum(double( abs(theta_z_print)>small_tol )),Plen,small_tol );
            fprintf('), %d/%d largest terms\n', ...
            sum(double( abs(theta_z_print)>small_tol )),Plen );

            for jdep = 1:ndep_

                theta_z = theta_mat(:,jdep+1);
                theta_z_maxmag = ( max(abs( theta_z )) );
                theta_z_print = theta_z/theta_z_maxmag;

                fprintf('   v_u_%d (x,u) = (%.1e)( ',jdep,theta_z_maxmag);
                for i = 1:Plen
                    if ( abs(theta_z_print(i)) > small_tol )
                        fprintf('+ (%.2f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                        for idep = 1:ndep_
                            fprintf(',%d',Pmat_(idep+1,i) );
                        end
                        fprintf('] ');
                    end
                end
                fprintf('), %d/%d largest terms\n', ...
                sum(double( abs(theta_z_print)>small_tol )),Plen );
            end
            fprintf('(small_tol = %.1e)\n', small_tol);

        end

    end

end
