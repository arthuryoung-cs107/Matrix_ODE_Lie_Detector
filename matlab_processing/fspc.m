classdef fspc < jspc

    properties

        bor;

    end

    methods

        function obj = fspc(obs_)

            obj@jspc(obs_.ndep,obs_.eor);

        end

    end

    methods (Static)

        function str_out = compute_svd_package(A_)
            A_comp = A_;
            [W_A,r_A,s_A,V_A,~] = fspc.safely_process_net_svd(A_comp);
            str_out = struct( ...
                'ncol', length( s_A ), ...
                'W', W_A, ...
                'r', r_A, ...
                's', s_A, ...
                'V', V_A, ...
                'Y',  V_A.*( (s_A(:)/s_A(1))' ) ...
            );
        end

        function [s_scaled_pseudo_inverse,s_scl] = s_spi(s_)

            % s_scaled_pseudo_inverse = reshape(1.0./(s_/s_(end)),1,[]);
            ss = s_;
            logc_zero = (s_ == 0.0);
            [s_scl,imin_nzero] = min( s_(~logc_zero) );
            ss((imin_nzero+1):end) = s_scl;

            s_scaled_pseudo_inverse = reshape(1.0./(ss/s_scl),1,[]);

        end

        function [WA_,rA_,sA_,VA_,Afull_] = safely_process_net_svd(Atns_)
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
            [rA_,sA_,VA_] = fspc.rsV_unpack( Afull_ );

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

        function [r_,s_,V_] = rsV_unpack(A_)
            [~,S_,V_] = svd(A_,'econ');
            s_ = reshape(diag(S_),[],1); % output s_ is col vector
            r_ = sum(double(s_ > max(size(A_))*eps(S_(1)))); % default matlab tol
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
