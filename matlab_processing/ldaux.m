classdef ldaux

    properties (Constant)
        %% true utility
        iidiag = @(N_) logical( reshape(eye(N_),[],1) );

        %% jet space adjacent

        nvar = @(o_) 1+o_.ndep; % ambient dimension of embedding base space (0'th jet space)
        ndim = @(o_) 1+o_.ndep*(o_.eor+1); % ambient dimension of embedding N'th jet space
                               % [ndep_, eor_, nvar_, ndim_]
        jspc_dims = @(o_) deal(o_.ndep, o_.eor, ldaux.nvar(o_), ldaux.ndim(o_));

    end

    properties

    end

    methods

        function obj = ldaux()

        end

    end

    methods (Static)

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
            [Plen, Pmat, dim0_lens] = ldaux.count_set_P_len(bor,ndep);
            Pmat = Pmat';

            S1mat = Smat(1:(nvar+ndep),:);
            xumat = Smat(1:nvar,:);
            untns = reshape( Smat(2:end,:), ndep,eor+1,nobs );
            dxuntns = reshape(untns(:,2:end,:),ndep,eor,nobs);
            d1xumat = reshape(dxuntns(:,1,:),ndep,nobs);

            Lval_cell = cell(3,nobs);
            ilv = 1;
            for i = 1:nobs
                s0(i) = adlam(xumat(:,i),Omap_A,Omap_b);
                prl(i) = adlam.prolong_mvpolynomial(s0(i),dxuntns(:,:,i),bor);

                s0i = s0(i);

                Lval_cell{1,i} = adlam.trunc_to_adobj(s0i);
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
            mod_out.s0 = s0;
            mod_out.Lvals = Lvals;
            mod_out.lvals = lvals;
            mod_out.levls = levls;
            mod_out.dxlmat = dxlmat;
        end

        %% true utility functions


        %% jet space adjacent


        function [Smat,nobs,nset,kor_obs,ndim_obs] = unpack_Scell(Scell_,ndep_)

            nset = prod(size(Scell_));
            ndim_obs = size(Scell_{1},1);
            kor_obs = (ndim_obs-1)/ndep_ - 1;
            Smat = jspc.Scell_2_Smat( Scell_,ndim_obs );
            nobs = size(Smat,2);

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
            xu_out = reshape(xu_,jspc.nvar(obj),[]);
            nobs_out = size(xu_out,2);
        end
        function [P_len_out pow_mat_out dim0_lens_out] = count_set_P_len(bor_,ndep_)
            d_ = ndep_ + 1;
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
