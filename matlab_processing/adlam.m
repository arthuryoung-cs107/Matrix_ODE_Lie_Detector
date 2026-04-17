classdef adlam < adobj

    properties (Constant)

        A = @(o_) o_.Jac ;
        Omap = @(o_) o_.Jac * o_.xu + o_.b;
        vevl = @(o_,v_) o_.Jac * v_(:);

    end

    properties

        xu;
        b;

    end

    methods
        function obj = adlam(xu_,A_,b_)
            obj@adobj( A_*(xu_(:)) + b_(:) , reshape(A_,length(b_(:)),length(xu_(:))) );
            obj.xu = xu_(:);
            obj.b =  b_(:);
        end

    end

    methods (Static)

        function prn = prolong_mvpolynomial(s0_,dxu_,bor_)
            bor = double(bor_);
            nvar = length(s0_.xu);
            ndep = nvar-1;
            dxu = reshape(dxu_,ndep,[]);
            kor = size(dxu,2);
            kp1 = kor+1;

            xu = s0_.val;
            Oa = diag(s0_.Jac);

            c_tns = nan(nvar,bor,kor+1);
            L_cell = cell(nvar,bor,kor+1);
            for iv = 1:nvar
                ai = Oa(iv);
                L_cell{iv,1,1} = adobj( xu(iv) , ai );
                for ik = 2:kp1
                    L_cell{iv,1,ik} = adobj(L_cell{iv,1,ik-1}.Jac, 0.0);
                end

                for ib = 2:bor
                    % compute L_ib and its first partial
                    L_cell{iv,ib,1} = L_cell{iv,1,1}.*L_cell{iv,ib-1,1};

                    for ik = 2:kp1
                        if (ik<ib)
                            L_cell{iv,ib,ik} = ...
                                adobj( ...
                                L_cell{iv,ib,ik-1}.Jac, ...
                                
                                 );
                        elseif (ik==ib)

                        else

                        end
                    end

                    L_cell{iv,ib,ik} = adobj( pLs(iv,ib,ik-1), pLs(iv,ib,ik) );


                    % second partial to ib-1
                    for ik = 2:(ib-1)
                        pLs(iv,ib,ik) = ...
                            pLs(iv,ib,ik-1) ...
                            *( ai*(ib-ik+1) ) ...
                            *( L_cell{iv,ib-ik,1}.val );
                        L_cell{iv,ib,ik} = adobj( pLs(iv,ib,ik-1), pLs(iv,ib,ik) );
                    end

                    if (ik)
                    % partial of same order
                    ik = ib;
                    pLs(iv,ib,ik) = ...
                        pLs(iv,ib,ik-1) ...
                        *( ai*(ib-ik+1) ) ;
                    L_cell{iv,ib,ik} = adobj( pLs(iv,ib,ik-1), pLs(iv,ib,ik) );

                    for ik = (ib+1):kp1
                        L_cell{iv,ib,ik} = adobj( L_cell{iv,ib,ik-1}.Jac, 0.0);
                    end

                end
            end


            Lk_cell{1,1} = adobj(sh.val,eye(nvar));
            for ik = 1:kor
                iik = ik+1;
                Lk_cell{1,}
            end
            for ib = 2:bor_
                Lk_cell{1,ib} = Lk_cell{1,1}.*Lk_cell{1,ib-1};
                for ik = 1:kor

                end
            end

        end

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
