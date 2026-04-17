classdef adobj

    properties (Constant)

    end

    properties

        val;
        Jac;

    end

    methods
        function obj = adobj(val_,Jac_)
            [obj.val,obj.Jac] = deal(val_,Jac_);
        end

        function obj_o = qdim(obj_,q_)
            [vi,Ji] = adobj.unpack_valJac(obj_);
            obj_o = adobj(vi(q_),Ji(q_,:));
        end

        function obj_o = plus(obj_l_,obj_r_)
            [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
            obj_o = adobj(vl+vr,Jl+Jr);
        end
        function obj_o = minus(obj_l_,obj_r_)
            [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
            obj_o = adobj(vl-vr,Jl-Jr);
        end

        function obj_o = uminus(obj_)
            [vi,Ji] = adobj.unpack_valJac(obj_);
            obj_o = adobj(-vi,-Ji);
        end
        function obj_o = uplus(obj_)
            [vi,Ji] = adobj.unpack_valJac(obj_);
            obj_o = adobj(vi,Ji);
        end

        function obj_o = times(obj_l_,obj_r_)
            [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
            obj_o = adobj( vl.*vr , Jl.*vr + vl.*Jr );
        end
        function obj_o = mtimes(obj_l_,obj_r_)
            [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
            obj_o = adobj( vl*vr , Jl*vr + vl*Jr );
        end

        function obj_o = rdivide(obj_l_,obj_r_)
            [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
            obj_o = adobj( vl./vr , ( Jl.*vr - vl.*Jr )./( vr.^2 ) );
        end
        function obj_o = mrdivide(obj_l_,obj_r_)
            [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
            obj_o = adobj( vl/vr , ( Jl*vr - vl*Jr )/( vr^2 ) );
        end

        function obj_o = power(obj_l_,obj_r_)
            % [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
            % vo = vl.^vr;
            % obj_o = adobj( vo , vo.*( log(vl).*Jr + ( Jl./vl ).*vr ) );

            if ( isnumeric(obj_l_) ) % base is constant
                if ( isnumeric(obj_r_) ) % exponent is also constant. Jac = 0
                    obj_o = adobj( obj_l_.^obj_r_ , 0 );
                else % exponent is variable -> exponential
                    [vr,Jr] = deal(obj_r_.val,obj_r_.Jac);
                    vo = obj_l_.^vr;
                    obj_o = adobj( vo , vo.*log(obj_l_).*Jr  );
                end
            else
                [vl,Jl] = deal(obj_l_.val,obj_l_.Jac);
                if ( isnumeric(obj_r_) ) % exponent is a constant -> true power rule
                    vo = vl.^obj_r_;
                    obj_o = adobj( vo , vo.*Jl./vl.*obj_r_ );
                else % exponent is also variable -> generalized exponential
                    [vr,Jr] = deal(obj_r_.val,obj_r_.Jac);
                    vo = vl.^vr;
                    obj_o = adobj( vo , vo.*( log(vl).*Jr + Jl./vl.*vr ) );
                end
            end
        end
        function obj_o = mpower(obj_l_,obj_r_)
            % [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
            % vo = vl^vr;
            % obj_o = adobj( vl^vr , vo*( log(vl)*Jr + ( Jl/vl )*vr ) );

            if ( isnumeric(obj_l_) ) % base is constant
                if ( isnumeric(obj_r_) ) % exponent is also constant. Jac = 0
                    obj_o = adobj( obj_l_^obj_r_ , 0 );
                else % exponent is variable -> exponential
                    [vr,Jr] = deal(obj_r_.val,obj_r_.Jac);
                    vo = obj_l_^vr;
                    obj_o = adobj( vo , vo*log(obj_l_)*Jr  );
                end
            else
                [vl,Jl] = deal(obj_l_.val,obj_l_.Jac);
                if ( isnumeric(obj_r_) ) % exponent is a constant -> true power rule
                    vo = vl^obj_r_;
                    obj_o = adobj( vo , vo*Jl/vl*vr );
                else % exponent is also variable -> generalized exponential
                    [vr,Jr] = deal(obj_r_.val,obj_r_.Jac);
                    vo = vl^vr;
                    obj_o = adobj( vo , vo*( log(vl)*Jr + Jl/vl*vr ) );
                end
            end
        end

        function obj_o = ctranspose(obj_)
            [vi,Ji] = adobj.unpack_valJac(obj_);
            obj_o = adobj(vi',Ji');
        end
        function obj_o = transpose(obj_)
            [vi,Ji] = adobj.unpack_valJac(obj_);
            obj_o = adobj(vi.',Ji.');
        end

        %% glitchy, not clear if this is actually useful. Lazy list inits more useful
        % function out = subsref(obj_,s_)
        %     s_
        %
        %     switch s_.type
        %         case '()'
        %             [vi,Ji] = adobj.unpack_valJac(obj_);
        %             out = adobj( vi(s_.subs{1}),Ji(s_.subs{1},:) );
        %         case '.'
        %             out = obj_.(s_.subs);
        %         otherwise
        %             out = obj_.(s_.subs);
        %     end
        % end
        % function obj_o = subsasgn(obj_l_,s_,obj_r_)
        %     [vl,Jl,vr,Jr] = adobj.unpack_valJac_pair(obj_l_,obj_r_);
        %     [vo,Jo] = deal(vl,Jl);
        %     if (~isempty(s_))
        %         vo(s_.subs) = vr;
        %         Jo(s_.subs,:) = Jr;
        %     end
        %     obj_o = adobj( vo,Jo );
        % end

    end

    methods (Static)

        % "ad constant": constant value point in coordinates
        function obj_o = ad_constant(val_,nvar_)
            if (nargin == 2) % fix derivatives to 0 mat
                obj_o = adobj(val_,zeros(length(val_(:)),nvar_));
            else % fix derivatives to 0, generic
                obj_o = adobj(val_(:),0);
            end
        end
        % "ad seed": use this to seed automatic differentiation in generic coordinates
        function obj_o = ad_seed(val_)
            obj_o = adobj(val_(:),1.0);
        end
        % "ad identity": use this to seed automatic differentiation, enforcing native coords
        function obj_o = ad_identity(val_)
            obj_o = adobj(val_(:),eye(length(val_(:))));
        end




        function [vl,Jl,vr,Jr] = unpack_valJac_pair(obj_l_,obj_r_)
            if ( isnumeric(obj_l_) )
                [vl,Jl] = deal(obj_l_,0);
            else
                [vl,Jl] = deal(obj_l_.val,obj_l_.Jac);
            end
            if ( isnumeric(obj_r_) )
                [vr,Jr] = deal(obj_r_,0);
            else
                [vr,Jr] = deal(obj_r_.val,obj_r_.Jac);
            end
        end
        function [vo,Jo] = unpack_valJac(obj_)
            if ( isnumeric(obj_) )
                [vo,Jo] = deal(obj_,0);
            else
                [vo,Jo] = deal(obj_.val,obj_.Jac);
            end
        end
    end

end
