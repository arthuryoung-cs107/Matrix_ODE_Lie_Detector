classdef ldaux

    properties (Constant)
        %% true utility
        iidiag = @(N_) logical( reshape(eye(N_),[],1) );
        % iidiag = @(N_) sub2ind( [N_,N_],1:N_,1:N_ );

        %% jet space adjacent

        nvar = @(o_) 1+o_.ndep; % ambient dimension of embedding base space (0'th jet space)
        ndim = @(o_) 1+o_.ndep*(o_.eor+1); % ambient dimension of embedding N'th jet space
                               % [ndep_, eor_, nvar_, ndim_]
        jspc_dims = @(o_) deal(o_.ndep, o_.eor, ldaux.nvar(o_), ldaux.ndim(o_));

        Acell_2_Amat = @(Ac_) cell2mat( Ac_(:) );

    end

    properties

    end

    methods

        function obj = ldaux()

        end

    end

    methods (Static)
        function [Sobs,dat_out,JF,dNp1xu] = generate_double_oscillator_data()
            % equation specification
            fcn_name = 'generate_double_oscillator_data';
            eqn_name = 'double_oscillator';
            eor = 2;
            ndep = 2;
            ndim = 1 + ndep*(eor+1);
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);
            [mmass1,mmass2,kstiff1,kstiff2,kstiff3] = deal(0.5,1.25,0.75,1.25,1.75);
            cm = [ 1./mmass1 ; 1./mmass2 ].*[ ...
                -kstiff1 , kstiff2; ...
                -kstiff2 , -kstiff3 ...
            ];
            f_eqn = @(u1_,u2_,dxu1_,dxu2_) [ ...
                cm(1,1).*u1_ + cm(1,2).*(u2_-u1_) ; ...
                cm(2,1).*(u2_-u1_) + cm(2,2).*u2_ ;
            ];
            function [JF_out,dNp1xu_out] = JF_dxf_eqn(s_)
                npts_evl = size(s_,2);
                JF_out = zeros(ndep,ndim,npts_evl);
                dNp1xu_out = zeros(ndep,npts_evl);
                for i = 1:npts_evl
                    sol_i = adobj.seed_sol(s_(1:7,i));
                    % x_i = sol_i.qdim(1); % autonomous
                    u1_i = sol_i.qdim(2);
                    u2_i = sol_i.qdim(3);
                    dxu1_i = sol_i.qdim(4);
                    dxu2_i = sol_i.qdim(5);
                    d2xu1_i = sol_i.qdim(6);
                    d2xu2_i = sol_i.qdim(7);

                    udiff_ad = u2_i-u1_i;

                    f1_i = cm(1,1).*u1_i + cm(1,2).*udiff_ad ;
                    f2_i = cm(2,1).*udiff_ad + cm(2,2).*u2_i ;
                    F1_i = f1_i - d2xu1_i;
                    F2_i = f2_i - d2xu2_i;

                    [JF_out(1,:,i),JF_out(2,:,i)] = deal(F1_i.Jac,F2_i.Jac);

                    dNp1xu_out(:,i) = linsolve(JF_out(:,(end-ndep+1):end,i),-JF_out(:,1:(end-ndep),i)*[1;s_((2+ndep):end,i)]);
                end
            end
            Fode_sys_evl = @(e_,s_) [ ...
                ones(1,size(s_,2)) ; ...
                s_((ndep+2):end,:) ; ...
                f_eqn( s_(2,:) , s_(3,:), s_(4,:), s_(5,:) ) ...
            ];
            fode = struct( ...
            'name', eqn_name, ...
            'eor', eor, ...
            'ndep', ndep, ...
            'f', @(s_) f_eqn(s_(2,:),s_(3,:),s_(4,:),s_(5,:) ), ...
            'JF_dxf', @(s_) JF_dxf_eqn(s_) ...
            );
            % 'RelTol', 1e-3, ...
            % 'AbsTol', 1e-5, ...

            % ncrv = 28;
            ncrv = 10;
            seed = 6; % Shay's choice
            seed0 = rng(seed);

            % trajectory specification
            x0 = 0.0;
            % ef = 10.0; % epsilon varies from e0 = 0 to ef > 0
            ef = 9.5; % epsilon varies from e0 = 0 to ef > 0, but good to have >= 1 for default randomization

            u00_vals = 2*(rand(4,ncrv)-0.5);
            u0_vals = u00_vals .* ([ pi/4 ; pi/4 ; pi/8 ; pi/8 ]);
            % u0_vals = u00_vals .* ([ pi/4 ; pi/4 ; 0.0 ; 0.0 ]);

            %{
                Uniform count of observed points per curve (M)
                Induces MNQ R matrix constraints per curve,
                so we'd like M >= ((Q+1)*P)/(NQ) observations
                for a well defined P = (O+1)^(Q+1) local curve model,
                where O = 3 (cubic permutation) by default
                Does not actually need to be this many, but convenient.
            %}
            Odef = 3; % cubic permutation model (default)
            Pdef = (Odef+1)^(ndep+1);
            % nevl = ceil( (ndep+1)*Pdef/(eor*ndep) ) + 1;
            nevl = 3*(ceil(1* (ndep+1)*Pdef/(eor*ndep) ) + 1 ); % number of observed solutions per observed integral curve

            % s0 = [ x0*ones(1,ncrv); u0_vals ];
            s0 = [ x0 + rand(1,ncrv) ; u0_vals ];
            % epsf =  ef;
            epsf =  ef + rand(ncrv,1);
            % epsevl = linspace(0.0,ef,nevl);
            epsevl = cell([ncrv,1]);
            del_edge = rand(ncrv,2);
            for icrv = 1:ncrv
                epsevl{icrv,1} = linspace(del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
            end
            trj_specs = struct( ...
                's0', s0, ...
                'epsf',  epsf ...
            );
            trj_specs.epsevl = epsevl;

            [Sobs,trjs,JF,dNp1xu] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs);

            fprintf(' nobs=%d, ncrv=%d \n', ...
                ncrv*nevl, ncrv );

            dat_out = fode;

        end
        function [Sobs,dat_out,JF,dNp1xu] = generate_double_pendulum_data()
            % equation specification
            fcn_name = 'generate_double_pendulum_data';
            eqn_name = 'double_pendulum';
            eor = 2;
            ndep = 2;
            ndim = 1 + ndep*(eor+1);
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);
            [ggrav,llen] = deal(9.81,2.0);
            cm = [ ...
                4.0/3.0 , 0.5 , 0.5 , 1.5 ; ...
                0.5 , 1.0/3.0 , -0.5 , 0.5  ...
            ];
            cm(:,1:3) = llen*cm(:,1:3);
            cm(:,4) = ggrav*cm(:,4);
            function d2xu_out = f_eqn(u_,dxu_)
                udiff_F = u_(1)-u_(2);
                cos_udiff_F = cos(udiff_F);
                A_F = [ cm(1,1), cm(1,2)*cos_udiff_F  ; ...
                        cm(2,1)*cos_udiff_F, cm(2,2) ...
                ];

                sin_udiff_F = sin(udiff_F);
                
                b_F = -[ cm(1,3)*sin_udiff_F.*dxu_(2).*dxu_(2) + cm(1,4)*sin(u_(1)) ; ...
                        cm(2,3)*sin_udiff_F.*dxu_(1).*dxu_(1) + cm(2,4)*sin(u_(2)) ...
                ];
                d2xu_out = linsolve(A_F,b_F);
            end
            function [JF_out,dNp1xu_out] = JF_dxf_eqn(s_)
                npts_evl = size(s_,2);
                JF_out = zeros(ndep,ndim,npts_evl);
                dNp1xu_out = zeros(ndep,npts_evl);
                for i = 1:npts_evl
                    sol_i = adobj.seed_sol(s_(1:ndim,i));
                    % x_i = sol_i.qdim(1); % autonomous
                    u1_i = sol_i.qdim(2);
                    u2_i = sol_i.qdim(3);
                    dxu1_i = sol_i.qdim(4);
                    dxu2_i = sol_i.qdim(5);
                    d2xu1_i = sol_i.qdim(6);
                    d2xu2_i = sol_i.qdim(7);

                    udiff_ad = u1_i-u2_i;
                    cos_udiff_ad = adobj.cosine(udiff_ad);
                    sin_udiff_ad = adobj.sine(udiff_ad);

                    F1_i = cm(1,1).*d2xu1_i + cm(1,2).*d2xu2_i.*cos_udiff_ad ...
                            + cm(1,3).*dxu2_i.*dxu2_i.*sin_udiff_ad ...
                            + cm(1,4).*adobj.sine(u1_i);
                    F2_i = cm(2,1).*d2xu1_i.*cos_udiff_ad + cm(2,2).*d2xu2_i ...
                            + cm(2,3).*dxu1_i.*dxu1_i.*sin_udiff_ad ...
                            + cm(2,4).*adobj.sine(u2_i);
                    [JF_out(1,:,i),JF_out(2,:,i)] = deal(F1_i.Jac,F2_i.Jac);

                    dNp1xu_out(:,i) = linsolve(JF_out(:,(end-ndep+1):end,i),-JF_out(:,1:(end-ndep),i)*[1;s_((2+ndep):end,i)]);
                end
            end
            Fode_sys_evl = @(e_,s_) [   1 ; ...
                                        s_(4:end,1) ; ...
                                        f_eqn( s_(2:3,1) , s_(3:4,1) ) ...
            ];
            function d2xu_out = f_evl(s_)
                npts_evl = size(s_,2);
                d2xu_out = nan(2,npts_evl);
                for i = 1:npts_evl
                    d2xu_out(:,i) = f_eqn(s_(2:3,i),s_(4:5,i) );
                end
            end
            fode = struct( ...
            'name', eqn_name, ...
            'eor', eor, ...
            'ndep', ndep, ...
            'llen', llen, ...
            'f', @(s_) f_evl(s_), ...
            'JF_dxf', @(s_) JF_dxf_eqn(s_) ...
            );
            % 'RelTol', 1e-3, ...
            % 'AbsTol', 1e-5, ...

            % ncrv = 28;
            ncrv = 10;
            seed = 6; % Shay's choice
            seed0 = rng(seed);

            % trajectory specification
            x0 = 0.0;
            % ef = 10.0; % epsilon varies from e0 = 0 to ef > 0
            ef = 6.0; % epsilon varies from e0 = 0 to ef > 0, but good to have >= 1

            u00_vals = 2*(rand(4,ncrv)-0.5);
            u0_vals = u00_vals .* ([ pi/4 ; pi/4 ; pi/8 ; pi/8 ]);
            % u0_vals = u00_vals .* ([ pi/4 ; pi/4 ; 0.0 ; 0.0 ]);

            %{
                Uniform count of observed points per curve (M)
                Induces MNQ R matrix constraints per curve,
                so we'd like M >= ((Q+1)*P)/(NQ) observations
                for a well defined P = (O+1)^(Q+1) local curve model,
                where O = 3 (cubic permutation) by default
                Does not actually need to be this many, but convenient.
            %}
            Odef = 3; % cubic permutation model (default)
            Pdef = (Odef+1)^(ndep+1);
            % nevl = ceil( (ndep+1)*Pdef/(eor*ndep) ) + 1;
            nevl = 7*(ceil(1* (ndep+1)*Pdef/(eor*ndep) ) + 1 ); % number of observed solutions per observed integral curve

            % s0 = [ x0*ones(1,ncrv); u0_vals ];
            s0 = [ x0 + rand(1,ncrv) ; u0_vals ];
            % epsf =  ef;
            epsf =  ef + rand(ncrv,1);
            % epsevl = linspace(0.0,ef,nevl);
            epsevl = cell([ncrv,1]);
            del_edge = rand(ncrv,2);
            for icrv = 1:ncrv
                epsevl{icrv,1} = linspace(del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
            end
            trj_specs = struct( ...
                's0', s0, ...
                'epsf',  epsf ...
            );
            % 'epsevl', epsevl, ...
            trj_specs.epsevl = epsevl;
            % trj_specs = struct( ...
            %     's0', [ x0*ones(1,ncrv) ; u0_vals ], ...
            %     'epsevl', linspace(0.0,ef,nevl), ...
            %     'epsf',  ef ...
            % );

            [Sobs,trjs,JF,dNp1xu] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs);

            fprintf(' nobs=%d, ncrv=%d \n', ...
                ncrv*nevl, ncrv );

            dat_out = fode;

        end
        % function [Sobs,dat_out] = generate_oscillator_cart_data()
        function [Sobs,dat_out,JF,dNp1xu] = generate_oscillator_cart_data()
            fcn_name = 'generate_oscillator_cart_data';
            eqn_name = 'oscillator_cart';
            eor = 2;
            ndep = 2;
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);

            % [Sobs_1D,fode_1D] = ldaux.generate_oscillator_polr_data();
            [Sobs_1D,fode_1D,JF,dNp1xu] = ldaux.generate_oscillator_polr_data();

            fode_out = fode_1D;
            fode_out.name = eqn_name;
            fode_out.eor = eor;
            fode_out.ndep = ndep;

            rrad = fode_1D.rrad;
            ncrv = length(Sobs_1D);
            Sobs = cell([ncrv,1]);
            for icrv = 1:ncrv
                npts_i = size(Sobs_1D{icrv},2);
                Sobs_1D_i = Sobs_1D{icrv};
                sin_th_i = sin( Sobs_1D_i(2,:) );
                cos_th_i = cos( Sobs_1D_i(2,:) );

                Sobs{icrv} =  [ ...
                    Sobs_1D_i(1,:) ; ...
                    rrad*cos_th_i; ...
                    rrad*sin_th_i ; ...
                    -rrad*sin_th_i .* Sobs_1D_i(3,:) ; ...
                    rrad*cos_th_i .* Sobs_1D_i(3,:) ; ...
                    rrad*( -cos_th_i.*Sobs_1D_i(3,:).*Sobs_1D_i(3,:) - sin_th_i.*Sobs_1D_i(4,:) ) ; ...
                    rrad*( -sin_th_i.*Sobs_1D_i(3,:).*Sobs_1D_i(3,:) + cos_th_i.*Sobs_1D_i(4,:) )  ...
                ];
            end

            % fode_out.f = @(u1_,u2_,dxu1_,dxu2_) rrad*[-1 0; 0 1]*[]

            dat_out = fode_out;
        end
        function [Sobs,dat_out,JF,dNp1xu] = generate_oscillator_polr_data()
            % equation specification
            fcn_name = 'generate_oscillator_polr_data';
            eqn_name = 'oscillator_polr';
            eor = 2;
            ndep = 1;
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);
            % can be viewed as linear or torsional spring
            % [kstiff,cdamp,mmass,rrad] = deal(1.0,0.1,1.25,2.0);
            [kstiff,cdamp,mmass,rrad] = deal(1.0,0.0,1.25,2.0);
            inv_Iinrt = (rrad*rrad)/mmass;
            f_eqn = @(u_,dxu_) ( kstiff*u_ + cdamp*dxu_ )*(-inv_Iinrt);
            % Jf_eqn = @(u_,dxu_) ( kstiff*[0 1 0] + cdamp*[0 0 1] )*(-inv_Iinrt); % constant
            function Jf_out = Jf_eqn(s_)
                Jf_out = ((kstiff*[0 1 0] + cdamp*[0 0 1])*(-inv_Iinrt)) ...
                            .* ones(1,3,size(s_,2));
            end
            tau_Nm1_tns = @(s_) reshape( ...
                [ ones(1,size(s_,2)) ; s_(3,:) ; f_eqn(s_(2,:),s_(3,:)) ], ...
                3, 1, size(s_,2) ...
            );
            dxf_eqn = @(J_,s_) reshape(pagemtimes(J_,tau_Nm1_tns(s_)), 1,size(s_,2));
            Jf2JF = @(Jf_) permute(reshape([ Jf_ ; -ones(1,size(Jf_,2)) ],4,1,size(Jf_,2)),[2 1 3]);
            JF_dxf_eqn = @(J_,s_) deal( Jf2JF(reshape(J_,3,[])), dxf_eqn(J_,s_) );
            Fode_sys_evl = @(e_,s_) [ ...
                ones(1,size(s_,2)) ; ...
                s_((ndep+2):end,:) ; ...
                f_eqn( s_(2,:) , s_(3,:) ) ...
            ];
            fode = struct( ...
            'name', eqn_name, ...
            'eor', eor, ...
            'ndep', ndep, ...
            'rrad', rrad, ...
            'f', @(s_) f_eqn( s_(2,:),s_(3,:) ), ...
            'JF_dxf', @(s_) JF_dxf_eqn( Jf_eqn(s_), s_(1:3,:) ) ....
            );
            % 'Jf', @(s_) Jf_eqn( s_(1:3,:) ), ...
            % 'dxf', @(s_) dxf_eqn( Jf_eqn(s_), s_(1:3,:) ), ...
            % 'f_ad', @(s_) f_eqn( s_.qdim(2), s_.qdim(3) ) ...

            % ncrv = 28;
            ncrv = 10;
            seed = 6; % Shay's choice
            seed0 = rng(seed);

            % trajectory specification
            x0 = 0.0;
            % ef = 10.0; % epsilon varies from e0 = 0 to ef > 0
            ef = 9.0; % epsilon varies from e0 = 0 to ef > 0

            u00_vals = 2*(rand(2,ncrv)-0.5);
            u0_vals = u00_vals .* [ 3*pi/4 ; pi/2 ];

            %{
                Uniform count of observed points per curve (M)
                Induces MNQ R matrix constraints per curve,
                so we'd like M >= ((Q+1)*P)/(NQ) observations
                for a well defined P = (O+1)^(Q+1) local curve model,
                where O = 3 (cubic permutation) by default
                Does not actually need to be this many, but convenient.
            %}
            Odef = 3; % cubic permutation model (default)
            Pdef = (Odef+1)^(ndep+1);
            % nevl = ceil( (ndep+1)*Pdef/(eor*ndep) ) + 1;
            nevl = 10*( ceil(1* (ndep+1)*Pdef/(eor*ndep) ) + 1 ); % number of observed solutions per observed integral curve

            % s0 = [ x0*ones(1,ncrv); u0_vals ];
            s0 = [ x0 + rand(1,ncrv) ; u0_vals ];
            % epsf =  ef;
            epsf =  ef + rand(ncrv,1);
            % epsevl = linspace(0.0,ef,nevl);
            epsevl = cell([ncrv,1]);
            del_edge = rand(ncrv,2);
            for icrv = 1:ncrv
                % epsevl{icrv,1} = linspace(x0 + del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
                epsevl{icrv,1} = linspace(del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
            end
            trj_specs = struct( ...
                's0', s0, ...
                'epsf',  epsf ...
            );
            % 'epsevl', epsevl, ...
            trj_specs.epsevl = epsevl;
            % trj_specs = struct( ...
            %     's0', [ x0*ones(1,ncrv) ; u0_vals ], ...
            %     'epsevl', linspace(0.0,ef,nevl), ...
            %     'epsf',  ef ...
            % );

            [Sobs,trjs,JF,dNp1xu] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs);

            fprintf(' nobs=%d, ncrv=%d \n', ...
                ncrv*nevl, ncrv );

            dat_out = fode;

        end
        function [Sobs,dat_out] = generate_pendulum_cart_data()
            fcn_name = 'generate_pendulum_cart_data';
            eqn_name = 'pendulum_cart';
            eor = 2;
            ndep = 2;
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);

            [Sobs_1D,fode_1D] = ldaux.generate_pendulum_polr_data();
            rrad = fode_1D.rrad;

            fode_out = fode_1D;
            fode_out.name = eqn_name;
            fode_out.eor = eor;
            fode_out.ndep = ndep;
            % fode_out.f =
            fode_out.upolr_2_ucart_ad = @(u_) rrad*[ adobj.sine(u_) -adobj.cosine(u_) ];


            ncrv = length(Sobs_1D);
            Sobs = cell([ncrv,1]);
            for icrv = 1:ncrv
                npts_i = size(Sobs_1D{icrv},2);
                Sobs_1D_i = Sobs_1D{icrv};
                sin_th_i = sin( Sobs_1D_i(2,:) );
                cos_th_i = cos( Sobs_1D_i(2,:) );

                Sobs{icrv} =  [ ...
                    Sobs_1D_i(1,:) ; ...
                    rrad*sin_th_i; ...
                    -rrad*cos_th_i ; ...
                    rrad*cos_th_i .* Sobs_1D_i(3,:) ; ...
                    rrad*sin_th_i .* Sobs_1D_i(3,:) ; ...
                    rrad*( -sin_th_i.*(Sobs_1D_i(3,:).*Sobs_1D_i(3,:)) + cos_th_i.*Sobs_1D_i(4,:) ) ; ...
                    rrad*( cos_th_i.*(Sobs_1D_i(3,:).*Sobs_1D_i(3,:)) + sin_th_i.*Sobs_1D_i(4,:) )  ...
                ];
            end
            dat_out = fode_out;
        end
        function [Sobs,dat_out,JF,dNp1xu] = generate_pendulum_polr_data()
            % equation specification
            fcn_name = 'generate_pendulum_polr_data';
            eqn_name = 'pendulum_polr';
            eor = 2;
            ndep = 1;
            ndim = 1 + ndep*(eor+1);
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);
            mass = 1.0;
            % [rrad,mmass,cdamp,ggrav] = deal(2.0,1.0,1.0,9.81);
            [rrad,mmass,cdamp,ggrav] = deal(2.0,1.0,0.0,9.81);
            [c1,c2] = deal(ggrav/rrad,cdamp/(mmass*rrad*rrad));
            f_eqn = @(u_,dxu_) -c1*sin(u_) - c2*dxu_ ;
            Fode_sys_evl = @(e_,s_) [ ...
                ones(1,size(s_,2)) ; ...
                s_((ndep+2):end,:) ; ...
                f_eqn( s_(2,:) , s_(3,:) ) ...
            ];
            function [JF_out,dNp1xu_out] = JF_dxf_eqn(s_)
                npts_evl = size(s_,2);
                JF_out = zeros(1,ndim,npts_evl);
                dNp1xu_out = zeros(1,npts_evl);
                for i = 1:npts_evl
                    sol_i = adobj.seed_sol(s_(1:ndim,i));
                    % x_i = sol_i.qdim(1); % autonomous
                    u1_i = sol_i.qdim(2);
                    dxu1_i = sol_i.qdim(3);
                    d2xu1_i = sol_i.qdim(4);

                    f1_i = -c1*adobj.sine(u1_i) - c2*dxu1_i ; % forward/evolutionary model
                    F1_i = f1_i - d2xu1_i;

                    JF_out(1,:,i) = F1_i.Jac;
                end
            end
            fode = struct( ...
            'name', eqn_name, ...
            'eor', eor, ...
            'ndep', ndep, ...
            'rrad', rrad, ...
            'f', @(s_) f_eqn( s_(2,:),s_(3,:) ), ...
            'f_ad', @(s_) f_eqn( adobj(s_(2),[0 1 0]), adobj(s_(2),[0 0 1]) ), ...
            'JF_dxf', @(s_) JF_dxf_eqn(s_) ...
            );

            ncrv = 10;
            seed = 6; % Shay's choice
            seed0 = rng(seed);

            % trajectory specification
            x0 = 0.0;
            % ef = 10.0; % epsilon varies from e0 = 0 to ef > 0
            ef = 9.0; % epsilon varies from e0 = 0 to ef > 0

            u00_vals = 2*(rand(2,ncrv)-0.5);
            u0_vals = u00_vals .* [ 0.25*pi/2 ; 0.2*pi/2 ];
            % u0_vals = u00_vals .* [ 0.5*pi/2 ; 0.0 ];
            % u0_vals = u00_vals .* [ 2 ; 0 ];

            %{
                Uniform count of observed points per curve (M)
                Induces MNQ R matrix constraints per curve,
                so we'd like M >= ((Q+1)*P)/(NQ) observations
                for a well defined P = (O+1)^(Q+1) local curve model,
                where O = 3 (cubic permutation) by default
                Does not actually need to be this many, but convenient.
            %}
            Odef = 3; % cubic permutation model (default)
            Pdef = (Odef+1)^(ndep+1);
            nevl = 7*( ceil(1* (ndep+1)*Pdef/(eor*ndep) ) + 1 ); % number of observed solutions per observed integral curve

            % s0 = [ x0*ones(1,ncrv); u0_vals ];
            s0 = [ x0 + rand(1,ncrv) ; u0_vals ];
            % epsf =  ef;
            epsf =  ef + rand(ncrv,1);
            % epsevl = linspace(0.0,ef,nevl);
            epsevl = cell([ncrv,1]);
            del_edge = rand(ncrv,2);
            for icrv = 1:ncrv
                % epsevl{icrv,1} = linspace(x0 + del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
                epsevl{icrv,1} = linspace(del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
            end
            trj_specs = struct( ...
                's0', s0, ...
                'epsf',  epsf ...
            );
            % 'epsevl', epsevl, ...
            trj_specs.epsevl = epsevl;

            % [Sobs,trjs] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs);
            [Sobs,trjs,JF,dNp1xu] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs);

            fprintf(' nobs=%d, ncrv=%d (nevl=%d) \n', ...
                ncrv*nevl, ncrv, nevl);

            dat_out = fode;

        end
        function [Sobs,dat_out,JF,dNp1xu] = generate_Van_der_Pol_data()
            % equation specification
            fcn_name = 'generate_Van_der_Pol_data';
            eqn_name = 'VanderPol';
            eor = 2;
            ndep = 1;
            ndim = 1 + ndep*(eor+1);
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n ' , ...
                fcn_name,eqn_name,eor,ndep);
            mu = 1.75;
            f_eqn = @(u_,dxu_) mu.*( (1.0-( u_ .* u_ )).*dxu_ ) - u_;
            Fode_sys_evl = @(e_,s_) [ ...
                ones(1,size(s_,2)) ; ...
                s_((ndep+2):end,:) ; ...
                f_eqn( s_(2,:) , s_(3,:) ) ...
            ];
            function [JF_out,dNp1xu_out] = JF_dxf_eqn(s_)
                npts_evl = size(s_,2);
                JF_out = zeros(1,ndim,npts_evl);
                dNp1xu_out = zeros(1,npts_evl);
                for i = 1:npts_evl
                    sol_i = adobj.seed_sol(s_(1:ndim,i));
                    % x_i = sol_i.qdim(1); % autonomous
                    u1_i = sol_i.qdim(2);
                    dxu1_i = sol_i.qdim(3);
                    d2xu1_i = sol_i.qdim(4);

                    f1_i = mu.*( (1.0-( u1_i .* u1_i )).*dxu1_i ) - u1_i; % forward/evolutionary model
                    F1_i = f1_i - d2xu1_i;

                    JF_out(1,:,i) = F1_i.Jac;
                    dNp1xu_out(:,i) = linsolve(JF_out(:,(end-ndep+1):end,i),-JF_out(:,1:(end-ndep),i)*[1;s_((2+ndep):end,i)]);
                end
            end
            fode = struct( ...
            'name', eqn_name, ...
            'eor', eor, ...
            'ndep', ndep, ...
            'f', @(s_) f_eqn( s_(2,:),s_(3,:) ), ...
            'f_ad', @(s_) f_eqn( adobj(s_(2),[0 1 0]), adobj(s_(2),[0 0 1]) ), ...
            'JF_dxf', @(s_) JF_dxf_eqn(s_) ...
            );

            ncrv = 10;
            seed = 6; % Shay's choice
            seed0 = rng(seed);

            % trajectory specification
            x0 = 0.0;
            % ef = 30.0; % epsilon varies from e0 = 0 to ef > 0
            ef = 9.0; % epsilon varies from e0 = 0 to ef > 0

            u00_vals = 2*(rand(2,ncrv)-0.5);
            u0_vals = u00_vals .* [ 2 ; 2 ];
            % u0_vals = u00_vals .* [ 2 ; 0 ];

            %{
                Uniform count of observed points per curve (M)
                Induces MNQ R matrix constraints per curve,
                so we'd like M >= ((Q+1)*P)/(NQ) observations
                for a well defined P = (O+1)^(Q+1) local curve model,
                where O = 3 (cubic permutation) by default
                Does not actually need to be this many, but convenient.
            %}
            Odef = 3; % cubic permutation model (default)
            Pdef = (Odef+1)^(ndep+1);
            nevl = 10*(ceil(1* (ndep+1)*Pdef/(eor*ndep) ) + 1);

            s0 = [ x0 + rand(1,ncrv) ; u0_vals ];
            epsf =  ef + rand(ncrv,1);
            epsevl = cell([ncrv,1]);
            del_edge = rand(ncrv,2);
            for icrv = 1:ncrv
                % epsevl{icrv,1} = linspace(x0 + del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
                epsevl{icrv,1} = linspace(del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
            end
            trj_specs = struct( ...
                's0', s0, ...
                'epsf',  epsf ...
            );
            trj_specs.epsevl = epsevl;
            % trj_specs = struct( ...
            %     's0', [ x0*ones(1,ncrv) ; u0_vals ], ...
            %     'epsevl', linspace(0.0,ef,nevl), ...
            %     'epsf',  ef ...
            % );

            % [Sobs,trjs] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs);
            [Sobs,trjs,JF,dNp1xu] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs);

            fprintf(' nobs=%d, ncrv=%d \n', ...
                ncrv*nevl, ncrv );

            dat_out = fode;

        end

        function [Sobs,dat_out,JF,dNp1xu] = generate_Brusselator_data()
            % equation specification
            fcn_name = 'generate_Brusselator_data';
            eqn_name = 'Brusselator';
            eor = 1;
            ndep = 2;
            ndim = 1 + ndep*(eor+1);
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);
            cm = [ ...
                1.0,-4.0,1.0; ...
                0.0,3.0,-1.0 ...
            ];
            f_eqn = @(u1_,u2_) cm(:,1) + cm(:,2).*u1_ + cm(:,3).*( u1_.*u1_.*u2_ );
            function [JF_out,dNp1xu_out] = JF_dxf_eqn(s_)
                npts_evl = size(s_,2);
                JF_out = zeros(ndep,ndim,npts_evl);
                dNp1xu_out = zeros(ndep,npts_evl);
                for i = 1:npts_evl
                    sol_i = adobj.seed_sol(s_(1:ndim,i));
                    % x_i = sol_i.qdim(1); % autonomous
                    u1_i = sol_i.qdim(2);
                    u2_i = sol_i.qdim(3);
                    dxu1_i = sol_i.qdim(4);
                    dxu2_i = sol_i.qdim(5);

                    u1u1u2_i = u1_i.*u1_i.*u2_i;

                    f1_i = cm(1,1) + cm(1,2).*u1_i + cm(1,3).*( u1u1u2_i );
                    f2_i = cm(2,1) + cm(2,2).*u1_i + cm(2,3).*( u1u1u2_i );
                    F1_i = f1_i - dxu1_i;
                    F2_i = f2_i - dxu2_i;

                    [JF_out(1,:,i),JF_out(2,:,i)] = deal(F1_i.Jac,F2_i.Jac);

                    dNp1xu_out(:,i) = linsolve(JF_out(:,(end-ndep+1):end,i),-JF_out(:,1:(end-ndep),i)*[1;s_((2+ndep):end,i)]);
                end
            end
            Fode_sys_evl = @(e_,s_) [ ...
                ones(1,size(s_,2)) ; ...
                f_eqn( s_(2,:) , s_(3,:) ) ...
            ];
            fode = struct( ...
            'name', eqn_name, ...
            'eor', eor, ...
            'ndep', ndep, ...
            'f', @(s_) f_eqn( s_(2,:),s_(3,:) ), ...
            'JF_dxf', @(s_) JF_dxf_eqn(s_) ...
            );
            x0 = 0.0;
            ef = 9.5; % epsilon varies from e0 = 0 to ef > 0

            ncrv = 10;
            seed = 34; % nice to look at
            seed0 = rng(seed);
            u00_vals = 2*(rand(2,ncrv)-0.5);
            % figure
            % scatter(u00_vals(1,:), u00_vals(2,:))
            % pause
            u0_vals = u00_vals .* (0.75*[ 1.5 ; 3.0 ]) + [ 1.5 ; 3.0 ];

            %{
                Uniform count of observed points per curve (M)
                Induces MNQ R matrix constraints per curve,
                so we'd like M >= ((Q+1)*P)/(NQ) observations
                for a well defined P = (O+1)^(Q+1) local curve model,
                where O = 3 (cubic permutation) by default
                Does not actually need to be this many, but convenient.
            %}
            Odef = 3; % cubic permutation model (default)
            Pdef = (Odef+1)^(ndep+1);
            nevl = 3*ceil( (ndep+1)*Pdef/(eor*ndep) ) + 1;

            trj_specs = struct( ...
                's0', [ x0*ones(1,ncrv) ; u0_vals ], ...
                'epsevl', linspace(0.0,ef,nevl), ...
                'epsf',  ef ...
            );

            [Sobs,trjs,JF,dNp1xu] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs)

            fprintf(' nobs=%d, ncrv=%d \n', ...
                ncrv*nevl, ncrv );

            dat_out = fode;

        end
        function [Sobs,dat_out,JF,dNp1xu] = generate_Riccati_data_2()
            % equation specification
            fcn_name = 'generate_Riccati_data';
            eqn_name = 'Riccati';
            eor = 1;
            ndep = 1;
            ndim = 1 + ndep*(eor+1);
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);
            cm = [ ...
                2.0,-1.0 ...
            ];
            f_eqn = @(x_,u_) cm(1,1).*( u_ ./ x_ ) + cm(1,2).*(x_.*x_).*(u_.*u_) ;
            function [JF_out,dNp1xu_out] = JF_dxf_eqn(s_)
                npts_evl = size(s_,2);
                JF_out = zeros(ndep,ndim,npts_evl);
                dNp1xu_out = zeros(ndep,npts_evl);
                for i = 1:npts_evl
                    sol_i = adobj.seed_sol(s_(1:ndim,i));
                    x_i = sol_i.qdim(1);
                    u1_i = sol_i.qdim(2);
                    dxu1_i = sol_i.qdim(3);

                    f1_i = cm(1,1).*(u1_i./x_i) + cm(1,2).*x_i.*x_i.*u1_i.*u1_i;
                    F1_i = f1_i - dxu1_i;

                    JF_out(1,:,i) = F1_i.Jac;

                    dNp1xu_out(:,i) = linsolve(JF_out(:,(end-ndep+1):end,i),-JF_out(:,1:(end-ndep),i)*[1;s_((2+ndep):end,i)]);
                end
            end
            Fode_sys_evl = @(e_,s_) [ ...
                ones(1,size(s_,2)) ; ...
                f_eqn( s_(1,:) , s_(2,:) ) ...
            ];
            fode = struct( ...
            'name', eqn_name, ...
            'eor', eor, ...
            'ndep', ndep, ...
            'f', @(s_) f_eqn(s_(1,:),s_(2,:)), ...
            'JF_dxf', @(s_) JF_dxf_eqn(s_) ...
            );

            x0 = 1e-1;
            ef = 10;

            ncrv = 10;
            % seed = 34; % nice to look at
            % seed0 = rng(seed);
            % u00_vals = rand(1,ncrv);
            % u0_vals = u00_vals .* (1e1);
            u0_vals = logspace(-1,1,ncrv);

            % nevl = 33; % one more than the cubic Rmat curve matrix
            Odef = 3; % cubic permutation model (default)
            Pdef = (Odef+1)^(ndep+1);
            nevl = 2*(ceil(1* (ndep+1)*Pdef/(eor*ndep) ) + 1);

            s0 = [ x0 + rand(1,ncrv) ; u0_vals ];
            epsf =  ef + rand(ncrv,1);
            epsevl = cell([ncrv,1]);
            del_edge = rand(ncrv,2);
            for icrv = 1:ncrv
                epsevl{icrv,1} = linspace(del_edge(icrv,1), ef - del_edge(icrv,1), nevl);
            end
            trj_specs = struct( ...
                's0', s0, ...
                'epsf',  epsf ...
            );
            trj_specs.epsevl = epsevl;

            [Sobs,trjs,JF,dNp1xu] = ldaux.evaluate_trajectories(Fode_sys_evl,fode,trj_specs);

            fprintf(' nobs=%d, ncrv=%d \n', ...
                ncrv*nevl, ncrv );

            % dat_out = struct( ...
            %     'name', eqn_name, ...
            %     'eor', eor, ...
            %     'ndep', ndep ...
            % );
            dat_out = fode;

        end
        function [Sobs,dat_out] = generate_Riccati_data()
            % equation specification
            fcn_name = 'generate_Riccati_data';
            eqn_name = 'Riccati';
            eor = 1;
            ndep = 1;
            fprintf('(ldaux::%s) Generating %s observations (ndep=%d, eor=%d) : \n' , ...
                fcn_name,eqn_name,eor,ndep);
            f_eqn = @(x_,u_) 2.0*( u_ ./ x_ ) - (x_.*x_).*(u_.*u_) ;
            gradf_eqn = @(x_,u_) [ ...
                -2.0*(u_ ./ (x_.*x_)) - 2.0*(x_).*(u_.*u_) ; ...
                2.0*(1.0./x_) - 2.0*(x_.*x_).*u_ ...
            ];
            dxf_eqn = @(x_,u_) sum( ...
                gradf_eqn( x_, u_(1,:) )' ...
                .* [ ones(size(x_))' , u_(2,:)' ] ...
                , 2  )';

            fode = struct( ...
            'name', eqn_name, ...
            'eor', eor, ...
            'ndep', ndep, ...
            'f', @(x_,u_) f_eqn(x_,u_), ...
            'gradf', @(x_,u_) gradf_eqn(x_,u_), ...
            'dxf', @(x_,u_) dxf_eqn(x_,u_) ...
            );

            xu_check = [ 1e-1 ; 1e1 ];
            xu_ad = adobj(xu_check,eye(length(xu_check)));
            x_ad = xu_ad.qdim(1);
            u_ad = xu_ad.qdim(2);
            f_check = fode.f(xu_check(1),xu_check(2));
            g_check = fode.gradf(xu_check(1),xu_check(2));
            f_ad = 2.0.*( u_ad ./ x_ad ) - (x_ad.*x_ad).*(u_ad.*u_ad);
            fmap_ad = @(x_,u_) 2.0*( u_ ./ x_ ) - (x_.*x_).*(u_.*u_);
            fm_ad = fmap_ad( xu_ad.qdim(1),xu_ad.qdim(2) );
            if ( fm_ad.Jac(:) ~= g_check(:) )
                fprintf('autodiff is broken\n');
                pause
            end
            Fode_sys_evl = @(e_,s_) [ ones(1,size(s_,2)) ; f_eqn( s_(1,:),s_(2:end,:) ) ];

            x0 = 1e-1;
            ncrv = 10;

            u0_vals = logspace(-1,1,ncrv);
            ef = 2.0; % epsilon varies from e0 = 0 to ef > 0
            xf = x0 + ef;

            nevl = 33; % one more than the cubic Rmat curve matrix

            epsevl = linspace(0.0,ef,nevl);
            x_evl = epsevl+x0;
            u_evl = nan(length(u0_vals),nevl );
            Sobs = cell(ncrv,1);
            for i = 1:length(u0_vals)
                Fode_sys = ode45(Fode_sys_evl,[0.0,ef],[x0 ; u0_vals(i)]);
                phi_xu_i = deval(Fode_sys,epsevl);

                u_evl(i,:) = phi_xu_i(2:end,:);

                Sobs{i} = [ phi_xu_i ; fode.f( phi_xu_i(1,:),phi_xu_i(2:end,:) ) ]; % sols are the graph of f
            end

            fprintf(' nobs=%d, ncrv=%d \n', ...
                ncrv*nevl, ncrv );

            % dat_out = struct( ...
            %     'name', eqn_name, ...
            %     'eor', eor, ...
            %     'ndep', ndep ...
            % );
            dat_out = fode;

        end

        function [Sobs_out, trjs,JF_out,dNp1xu_out] = evaluate_trajectories(Fevl_,fode_,trjs_)
            trjs = trjs_;
            ncrv = size(trjs.s0,2);
            sol_opts = odeset( ...
            'RelTol', 1e-12, ...
            'AbsTol', 1e-10 ...
            );
            if (isfield(fode_,'RelTol'))
                sol_opts.RelTol = fode_.RelTol;
            end
            if (isfield(fode_,'AbsTol'))
                sol_opts.AbsTol = fode_.AbsTol;
            end
            if (ncrv==1)
                if ( iscell(trjs.epsevl) )
                    epsevl = trjs.epsevl{1};
                else
                    epsevl = reshape(trjs.epsevl,1,[]);
                end
                phi_s_i = deval( ode45(Fevl_,[0.0 trjs.epsf(1)], trjs.s0(:,1), sol_opts) , epsevl );

                Sobs_out = [ ...
                phi_s_i ; ...
                fode_.f(phi_s_i) ...
                ]; % sols are the graph of f
            else
                if ( prod(size(trjs.epsf)) == 1 ) % provided as a scalar
                    epsf = trjs.epsf*ones(1,ncrv);
                else
                    epsf = trjs.epsf;
                end
                if ( iscell(trjs.epsevl) )
                    epsevl = trjs.epsevl;
                elseif ( prod(size(trjs.epsevl))==max(size(trjs.epsevl)) ) % if a vector
                    epsevl = num2cell( ones(ncrv,1)*reshape(trjs.epsevl,1,length(trjs.epsevl)) , 2 );
                else % otherwise, it should be a matrix
                    epsevl = num2cell(trjs.epsevl,2);
                end
                Sobs_out = cell([ncrv,1]);
                for i = 1:ncrv
                    % trjs.s0(:,i)
                    phi_s_i = deval( ode45(Fevl_,[0.0 epsf(i)],trjs.s0(:,i), sol_opts) , epsevl{i} );
                    Sobs_out{i} = [ ...
                    phi_s_i ; ...
                    fode_.f(phi_s_i) ...
                    ]; % sols are the graph of f
                end
            end
            if (isfield(fode_,'JF_dxf'))
                [JF_out,dNp1xu_out] = deal(cell([ncrv,1]));
                for i = 1:ncrv
                    [JF_out{i},dNp1xu_out{i}] = fode_.JF_dxf(Sobs_out{i});
                end
            else
                [JF_out,dNp1xu_out] = deal([]);
            end
        end

        %% simple jet space models

        function mod_out = first_order_ld_ad_cubic_model(S_,obs_)
            if ( isfield(obs_,'Omap_A') )
                Omap_A = obs_.Omap_A;
            elseif (isfield(obs_,'Omap_a'))
                Omap_A = (obs_.Omap_a(:)) .* eye(length(obs_.Omap_a(:))) ;
            else
                Omap_A = eye(1+obs_.ndep);
            end
            if ( isfield(obs_,'Omap_b') )
                Omap_b = obs_.Omap_b(:);
            else
                Omap_b = zeros(1+obs_.ndep,1);
            end

            bor = 3;

            [ndep,eor,nvar,ndim] = ldaux.jspc_dims(obs_);
            [Smat,nobs,nset,kor_obs,ndim_obs] = ldaux.unpack_Scell(S_,ndep);
            [Plen, Pmat, dim0_lens] = ldaux.count_set_P_len(bor,nvar);

            %{
                s = (x,u,dxu, ..., dnx u)
                xu = (x,u)
                un = (u,dxu, ..., dnx u) \in R^{Q \times (N+1)}
                dxu = (dxu, ..., dnx u)
                d1xu = ( dxu )
            %}

            S1mat = Smat(1:(nvar+ndep),:);
            xumat = Smat(1:nvar,:);
            untns = reshape( Smat(2:end,:), ndep,eor+1,nobs );
            dxuntns = reshape(untns(:,2:end,:),ndep,eor,nobs);
            d1xumat = reshape(dxuntns(:,1,:),ndep,nobs);

            poly_lam_spc = struct( ...
                'Omap_A', Omap_A, ...
                'Omap_b', Omap_b, ...
                'bor', bor, ...
                'Pmat', Pmat ...
            );
            % poly_lam_spc = 0;
            % poly_lam_spc = struct( ...
            %     'Pmat', Pmat ...
            % );

            Lval_cell = cell(3,nobs);
            ilv = 1;
            for i = 1:nobs
                % s0(i) = adlam(xumat(:,i),Omap_A,Omap_b);
                % s0(i) = adlam(xumat(:,i),poly_lam_spc);
                % lambdas(i) = adlam(xumat(:,i),Smat((nvar+1):end,i),poly_lam_spc);

                lambdas(i) = adlam(poly_lam_spc, xumat(:,i), Smat((nvar+1):end,i));

                s0(i) = lambdas(i).s0;

                % prl(i) = adlam.prolong_mvpolynomial(s0(i),dxuntns(:,:,i),bor);

                s0i = s0(i);

                Lval_cell{1,i} = s0i;
                Lval_cell{2,i} = Lval_cell{1,i}.*Lval_cell{1,i};
                Lval_cell{3,i} = Lval_cell{1,i}.*Lval_cell{2,i};
                Lvals(:,i) = [Lval_cell{1,i} ; Lval_cell{2,i} ; Lval_cell{3,i}];

                Lvals_full_i = [ adobj.ad_constant(ones(nvar,1),nvar); Lvals(:,i) ];

                for p = 1:Plen
                    lvip = Lvals_full_i( Pmat(1,p)+1 ).qdim(1);
                    for idep = 1:ndep
                        lvip = lvip.*( Lvals_full_i( Pmat(1+idep,p)+1 ).qdim(1+idep) );
                    end
                    lvals(ilv) = lvip;
                    % pr1l(ilv) = adlam.d1xl(lvip,d1xumat(:,i));

                    ilv = ilv + 1;
                end
            end
            lvals = reshape(lvals,Plen,nobs)
            levls = adlam.coordgrads2Jac(lvals)

            icheck = 21
            l_check = [levls(icheck).val, adlam.lvec(lambdas(icheck))]
            % lambdas(icheck).dL_tns
            lambdas(icheck)

            pause



            dxlmat = nan(Plen,nobs);
            for i = 1:nobs
                dxlmat(:,i) =  adlam.vevl( levls(i), [1 ; d1xumat(:,i)] );

                % for p = 1:Plen
                %
                % end
            end

            mod_out = obs_;
            mod_out.Omap_A = Omap_A;
            mod_out.Omap_b = Omap_b;
            mod_out.lambdas = lambdas;
            mod_out.s0 = s0;
            mod_out.Lvals = Lvals;
            mod_out.lvals = lvals;
            mod_out.levls = levls;
            mod_out.dxlmat = dxlmat;
        end

        %% true utility functions


        %% jet space adjacent


        function [Smat,nobs,nset,kor_obs,ndim_obs,npts_per_crv,ipts_crv] = unpack_Scell(Scell_,ndep_)

            nset = prod(size(Scell_));
            ndim_obs = size(Scell_{1},1);
            kor_obs = (ndim_obs-1)/ndep_ - 1;
            Smat = ldaux.Scell_2_Smat( Scell_,ndim_obs );
            nobs = size(Smat,2);
            if (nargout>=6)
                npts_per_crv = vertcat(cellfun( @(C_) size(C_,2),Scell_(:)));
                if (nargout>6)
                    csum_npts_per_crv = cumsum(npts_per_crv);
                    ipts_crv = [ [ 1 , (csum_npts_per_crv(1:(end-1))'+1) ] ; csum_npts_per_crv' ];
                end
            end
        end
        function Smat_out = Scell_2_Smat(Scell_,ndim_)
            if (nargin==1)
                ndim = size(Scell_{1},2);
            else
                ndim = ndim_;
            end
            Smat_out = reshape(cell2mat(cellfun( @(s_) s_(:) ,Scell_,'UniformOutput',false )),ndim,[]);
        end
        function [xu_out,nobs_out] = xumat_nobs(obj,xu_)
            xu_out = reshape(xu_,ldaux.nvar(obj),[]);
            nobs_out = size(xu_out,2);
        end
        % function [P_len_out pow_mat_out dim0_lens_out] = count_set_P_len(bor_,ndep_)
        function [P_len_out pow_mat_out dim0_lens_out] = count_set_P_len(bor_,d_)
            % d_ = ndep_ + 1;1
            P_len_out = (double(bor_+1))^(double(d_));

            dim0_lens_out = nan(bor_+1,1);
            pow_mat_out = nan(P_len_out,d_);

            k_perm = 0;
            for i = 0:bor_
                [delk pow_mat_out] = set_powmat_full_recursive(pow_mat_out,1,k_perm);
                dim0_lens_out(i+1) = delk;
                for k = k_perm:(k_perm+delk-1)
                    pow_mat_out(k+1,1) = i;
                end
                k_perm = k_perm + delk;
            end

            pow_mat_out = pow_mat_out'; % lambda map as a row is preferred now

            function [dk_out pm_out] = set_powmat_full_recursive(mat_, ilevel_, k_perm_)
                pm_out = mat_;
                if (ilevel_>=d_)
                    dk_out = 0;
                else
                    dk_level = 0;
                    if (ilevel_ == (d_-1))
                        for ib = 0:bor_
                            pm_out(k_perm_+dk_level+1,ilevel_+1) = ib;
                            dk_level = dk_level + 1;
                        end
                    else
                        for ib = 0:bor_
                            [dk pm_out] = set_powmat_full_recursive(pm_out,ilevel_+1,dk_level + k_perm_);
                            for k = 0:(dk-1)
                                pm_out(k_perm_+dk_level+1,ilevel_+1) = ib;
                                dk_level = dk_level + 1;
                            end
                        end
                    end
                    dk_out = dk_level;
                end
            end


        end
    end
end
