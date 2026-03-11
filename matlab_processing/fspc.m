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
        function print_short_polynomial_theta(theta_,Pmat_)
            ndep_ = size(Pmat_,1)-1;
            small_tol = 1e-6;

            fprintf('multivariate (1+Q = %d) polynomials of form v_z(x,u) = Sum_i theta_z(i)',1+ndep_);
            fprintf('*x^k(i,x)');
            for i = 1:ndep_
                fprintf('*u_%d^k(i,u_%d)',i,i);
            end
            fprintf(',\n');
            fprintf('   where theta_z : z = x, u_1, ..., u_Q, and k(i,z) is an assigned integer value, the exponent of z at i.\n');

            theta_maxmag = ( max(abs(theta_(:))) );
            theta_print = theta_/theta_maxmag;

            ntheta = length(theta_);
            Plen = ntheta/(1+ndep_);

            theta_mat = reshape(theta_,[],1+ndep_);

            theta_z = theta_mat(:,1);
            theta_z_maxmag = ( max(abs( theta_z )) );
            theta_z_print = theta_z/theta_z_maxmag;
            fprintf('   In short,');
            fprintf('   v_x (x,u) = (%.1f)( ',theta_z_maxmag);
            for i = 1:Plen
                if ( abs(theta_z_print(i)) > small_tol )
                    fprintf('+ (%.1f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                    for idep = 1:ndep_
                        fprintf(',%d',Pmat_(idep+1,i) );
                    end
                    fprintf('] ');
                end
            end
            fprintf('), %d (of %d) largest terms \n', ...
            sum(double( abs(theta_z_print)>small_tol )),Plen );

            for jdep = 1:ndep_

                theta_z = theta_mat(:,jdep+1);
                theta_z_maxmag = ( max(abs( theta_z )) );
                theta_z_print = theta_z/theta_z_maxmag;

                fprintf('   In short,');
                fprintf('   v_u_%d (x,u) = (%.1f)(',jdep,theta_z_maxmag);
                for i = 1:Plen
                    if ( abs(theta_z_print(i)) > small_tol )
                        fprintf('+ (%.1f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                        for idep = 1:ndep_
                            fprintf(',%d',Pmat_(idep+1,i) );
                        end
                        fprintf('] ');
                    end
                end
                fprintf('), %d (of %d) largest terms \n', ...
                sum(double( abs(theta_z_print)>small_tol )),Plen );
            end
            fprintf('\n');

        end
        function print_polynomial_theta(theta_,Pmat_)
            ndep_ = size(Pmat_,1)-1;
            small_tol = 1e-6;

            fprintf('multivariate (1+Q = %d) polynomials of form v_z(x,u) = Sum_i theta_z(i)',1+ndep_);
            fprintf('*x^k(i,x)');
            for i = 1:ndep_
                fprintf('*u_%d^k(i,u_%d)',i,i);
            end
            fprintf(',\n');
            fprintf('   where theta_z : z = x, u_1, ..., u_Q,');
            fprintf('   and k(i,z) is an assigned integer value, the exponent of z at i.\n');

            theta_maxmag = ( max(abs(theta_(:))) );
            theta_print = theta_/theta_maxmag;

            % fprintf('   In full, theta = (%.2e)( ', theta_maxmag);
            % fprintf('%.1e ', theta_print);
            % fprintf(')\n');

            ntheta = length(theta_);
            Plen = ntheta/(1+ndep_);

            theta_mat = reshape(theta_,[],1+ndep_);


            theta_z = theta_mat(:,1);
            theta_z_maxmag = ( max(abs( theta_z )) );
            theta_z_print = theta_z/theta_z_maxmag;
            fprintf('   v_x (x,u) = (%.2e)(\n',theta_z_maxmag);
            for i = 1:Plen
                fprintf('(%.2e) x^%d',theta_z_print(i),Pmat_(1,i));
                for idep = 1:ndep_
                    fprintf(' u_%d^%d',idep,Pmat_(idep+1,i) );
                end
                fprintf('\n');
            end
            fprintf(')\n');

            fprintf('   In short,');
            fprintf('   v_x (x,u) = (%.1f)( ',theta_z_maxmag);
            for i = 1:Plen
                if ( abs(theta_z_print(i)) > small_tol )
                    fprintf('+ (%.1f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                    for idep = 1:ndep_
                        fprintf(',%d',Pmat_(idep+1,i) );
                    end
                    fprintf('] ');
                end
            end
            fprintf(')\n');

            for jdep = 1:ndep_

                theta_z = theta_mat(:,jdep+1);
                theta_z_maxmag = ( max(abs( theta_z )) );
                theta_z_print = theta_z/theta_z_maxmag;

                fprintf('   v_u_%d (x,u) = (%.2e)(\n',jdep,theta_z_maxmag);
                for i = 1:Plen
                    fprintf('(%.2e) x^%d',theta_z_print(i),Pmat_(1,i));
                    for idep = 1:ndep_
                        fprintf(' u_%d^%d',idep,Pmat_(idep+1,i) );
                    end
                    fprintf('\n');
                end
                fprintf(')\n');

                fprintf('   In short,');
                fprintf('   v_u_%d (x,u) = (%.1f)(',jdep,theta_z_maxmag);
                for i = 1:Plen
                    if ( abs(theta_z_print(i)) > small_tol )
                        fprintf('+ (%.1f)xu^[%d',theta_z_print(i),Pmat_(1,i));
                        for idep = 1:ndep_
                            fprintf(',%d',Pmat_(idep+1,i) );
                        end
                        fprintf('] ');
                    end
                end
                fprintf(')\n');
            end
            fprintf('\n');

        end

    end

end
