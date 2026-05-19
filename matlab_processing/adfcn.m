classdef adfcn

    properties (Constant)

    end

    properties

        val; % scalar
        grad; % tall vector

    end

    methods
        function obj = adfcn(val_,grad_)
            [obj.val,obj.grad] = deal(val_(1),grad_(:));
        end

        %% overloaded operators

        function obj_o = plus(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn(vl+vr,gl+gr);
        end
        function obj_o = minus(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn(vl-vr,gl-gr);
        end

        function obj_o = uminus(obj_)
            [vi,Ji] = adfcn.unpack_valgrad(obj_);
            obj_o = adfcn(-vi,-Ji);
        end
        function obj_o = uplus(obj_)
            [vi,Ji] = adfcn.unpack_valgrad(obj_);
            obj_o = adfcn(vi,Ji);
        end

        function obj_o = times(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn( vl.*vr , gl.*vr + vl.*gr );
        end
        function obj_o = mtimes(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn( vl*vr , gl*vr + vl*gr );
        end

        function obj_o = rdivide(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn( vl./vr , ( gl.*vr - vl.*gr )./( vr.^2 ) );
        end
        function obj_o = mrdivide(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn( vl/vr , ( gl*vr - vl*gr )/( vr^2 ) );
        end

        function obj_o = power(obj_l_,obj_r_)

            if ( isnumeric(obj_l_) ) % base is constant
                if ( isnumeric(obj_r_) ) % exponent is also constant. grad = 0
                    obj_o = adfcn( obj_l_.^obj_r_ , 0 );
                else % exponent is variable -> exponential
                    [vr,gr] = deal(obj_r_.val,obj_r_.grad);
                    vo = obj_l_.^vr;
                    obj_o = adfcn( vo , vo.*log(obj_l_).*gr  );
                end
            else % base is variable
                [vl,gl] = deal(obj_l_.val,obj_l_.grad);
                if ( isnumeric(obj_r_) ) % exponent is a constant -> true power rule
                    vo = vl.^obj_r_;
                    obj_o = adfcn( vo , vo.*gl./vl.*obj_r_ );
                else % exponent is also variable -> generalized exponential
                    [vr,gr] = deal(obj_r_.val,obj_r_.grad);
                    vo = vl.^vr;
                    obj_o = adfcn( vo , vo.*( log(vl).*gr + gl./vl.*vr ) );
                end
            end
        end
        function obj_o = mpower(obj_l_,obj_r_)
            if ( isnumeric(obj_l_) ) % base is constant
                if ( isnumeric(obj_r_) ) % exponent is also constant. grad = 0
                    obj_o = adfcn( obj_l_^obj_r_ , 0 );
                else % exponent is variable -> exponential
                    [vr,gr] = deal(obj_r_.val,obj_r_.grad);
                    vo = obj_l_^vr;
                    obj_o = adfcn( vo , vo*log(obj_l_)*gr  );
                end
            else
                [vl,gl] = deal(obj_l_.val,obj_l_.grad);
                if ( isnumeric(obj_r_) ) % exponent is a constant -> true power rule
                    vo = vl^obj_r_;
                    obj_o = adfcn( vo , vo*gl/vl*vr );
                else % exponent is also variable -> generalized exponential
                    [vr,gr] = deal(obj_r_.val,obj_r_.grad);
                    vo = vl^vr;
                    obj_o = adfcn( vo , vo*( log(vl)*gr + gl/vl*vr ) );
                end
            end
        end

        function obj_o = ctranspose(obj_)
            [vi,gi] = adfcn.unpack_valgrad(obj_);
            obj_o = adfcn( vi', gi' );
        end
        function obj_o = transpose(obj_)
            [vi,gi] = adfcn.unpack_valgrad(obj_);
            obj_o = adfcn( vi.', gi.' );
        end

    end

    methods (Static)

        % "ad seed": use this to seed automatic differentiation in generic coordinates
        function obj_o = ad_seed(val_)
            obj_o = adfcn(val_(:),1.0);
        end

        function [vl,gl,vr,gr] = unpack_valgrad_pair(obj_l_,obj_r_)
            if ( isnumeric(obj_l_) )
                [vl,gl] = deal(obj_l_,0);
            else
                [vl,gl] = deal(obj_l_.val,obj_l_.grad);
            end
            if ( isnumeric(obj_r_) )
                [vr,gr] = deal(obj_r_,0);
            else
                [vr,gr] = deal(obj_r_.val,obj_r_.grad);
            end
        end
        function [vo,go] = unpack_valgrad(obj_)
            if ( isnumeric(obj_) )
                [vo,go] = deal(obj_,0);
            else
                [vo,go] = deal(obj_.val,obj_.grad);
            end
        end
    end

end
