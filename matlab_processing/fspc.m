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

        function print_short_polynomial_theta_z(theta_z_,Pmat_)
            ndep_ = size(Pmat_,1)-1;
            small_tol = 1e-6;

            fprintf('multivariate (1+Q = %d) polynomials of form\n',1+ndep_);
            fprintf('   v_z(xu) = Sum_i theta(i,z)*l_i(xu) ');
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
            fprintf('   v_z (x,u) = (%.1e)( ',theta_z_maxmag);
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

        function print_short_polynomial_theta(theta_,Pmat_)
            ndep_ = size(Pmat_,1)-1;
            small_tol = 1e-6;

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
            fprintf('   v_z (x,u) = (%.1e)( ',theta_z_maxmag);
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
