% classdef adlam < adobj
classdef adlam

    properties (Constant)

        Omap = @(o_) o_.Jac * o_.xu + o_.b;
        vevl = @(o_,v_) o_.Jac * v_(:);

    end

    properties

        xu;
        dxu;
        lam;

        A;
        b;

        s0;
        Lmat;
        c_tns;
        dL_tns;



        % pr1;
        % pr;

    end

    methods
        % function obj = adlam(xu_,A_,b_)
        %     obj@adobj( A_*(xu_(:)) + b_(:) , reshape(A_,length(b_(:)),length(xu_(:))) );
        %     obj.xu = xu_(:);
        %     obj.b =  b_(:);
        % end
        function obj = adlam(xu_,dxu_,lam_)
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

            % obj@adobj( A_*(xu_(:)) + b_(:) , reshape( A_, length(b_(:)), length(xu_(:)) ) );

            xu = xu_(:);
            dxu = reshape(dxu_,length(xu_(:))-1,[]);
            lam = lam_;

            b =  b_(:);
            A =  reshape( A_, length(b), length(xu) );

            s0 = adobj( A*xu + b, A );

            %% stage prolongation workspace
            bor = double(lam_.bor);
            Pmat = double(lam_.Pmat);
            Plen = size(Pmat,2);

            % nvar = length(xu);
            nvar = size(Pmat,1);
            ndep = nvar - 1;
            kor = size(dxu,2);
            kp1 = kor + 1;
            ndim = 1+(ndep*kp1);
            Oa = diag(A);
            xu_O = s0.val;

            Lmat = nan(nvar,bor); % matrix of values for each L polynomial
            c_tns = nan(nvar,bor,kor);
            dL_tns = nan(nvar,bor,kor);

            Lmat(:,1) = xu_O; % image of xu under diagonal linear map A * [x ; u] + b
            c_tns(:,1,1) = Oa; % first derivative chain coefficient
            dL_tns(:,1,1) = c_tns(:,1,1); % first order derivatives are chain coefficients

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

            Lmat_full = [ ones(nvar,1) , Lmat ];

            tau_k = [1 ; dxu(:) ];
            tau_0 = tau_k(1:nvar);

            tau_k_ad = adobj(tau_k, [ zeros(1,ndim) ; zeros(ndim-nvar,nvar), eye(ndim-nvar) ] );
            tau_0_ad = tau_k_ad.qdim(1:nvar)

            iiv = (1:nvar)';
            iiD_base = sub2ind([nvar,nvar],1:nvar,1:nvar);
            LmatP = ones(nvar,Plen);
            lvecP = ones(Plen,1);
            J_dkx_ltnsP = zeros( Plen,ndim, kor+1 );
            dkx_lmatP = zeros( Plen,kor );
            for ip = 1:Plen

                iip = Pmat( :,ip )
                iip_nz = iip~=0

                if (sum(iip_nz))
                    iiL = sub2ind([nvar,bor+1], iiv, iip+1 ) % linear indices of LmatP_i values in Lmat_full
                    LmatP(:,ip) = Lmat_full( iiL )
                    lvecP(ip) = prod( LmatP(:,ip) )

                    dL_mat_i_full = zeros(nvar,bor+1);
                    dL_mat_i_full(:,2:end) = dL_tns(:,:,1);

                    gLi = LmatP(:,ip)*ones(1,nvar);
                    gLi( iiD_base ) = dL_mat_i_full( iiL );
                    gli = prod(gLi,1);
                    J_dkx_ltnsP(ip,1:nvar,1) = gli;

                    dkx_lmatP(ip,1) = J_dkx_ltnsP(ip,1:nvar,1) * tau_0;

                    li_ad = adobj( LmatP(1,ip), [ dL_mat_i_full(iiL(1)) , zeros(1, ndim-1) ] );
                    for iv = 2:nvar
                        li_ad = li_ad * adobj( LmatP(iv,ip), [ zeros(1,iv-1) , dL_mat_i_full(iiL(iv)), zeros(1,ndim-iv) ] );
                    end

                    gLi
                    % HLi = 
                    gLi .* ones(nvar,ndim,nvar)

                    pause

                    li_ad.val,li_ad.Jac
                    tau_k_ad.val, tau_k_ad.Jac

                    % gli_ad = adobj( li_ad.Jac(:),prod(,1) );
                    % d1xli_ad = .qdim()
                    %
                    % for iv = 2:nvar
                    %
                    % end

                    % d1xli_ad.val,d1xli_ad.Jac

                    % J_dkx_ltnsP(ip,:,2) =

                    for ik = 2:kor

                        dkx_lmatP(ip,ik) = J_dkx_ltnsP(ip,:,ik) * tau_k;


                    end

                    pause

                end

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
            end

            % LmatP
            % lvecP


            %% initializations

            obj.xu = xu;
            obj.lam = lam;

            obj.A =  A;
            obj.b =  b;

            obj.s0 = s0;

            obj.Lmat = Lmat;
            obj.c_tns = c_tns;
            obj.dL_tns = dL_tns;
        end

    end

    methods (Static)

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
