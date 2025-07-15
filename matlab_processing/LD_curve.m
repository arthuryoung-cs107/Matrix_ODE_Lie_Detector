classdef LD_curve
    properties (Constant)

    end

    properties

        eor;
        ndep;

        pts_in;
        dnp1xu_in;
        JFs_in;

        npts;

        jets;
        rjets;
    end
    methods
        function obj = LD_curve(eor_,ndep_,pts_,dnp1xu_,JFs_)
            if (nargin == 3)
                dnp1xu_in = [];
                JFs_in = [];
            else
                dnp1xu_in = dnp1xu_(:);
                JFs_in = JFs_(:);
            end
            pts_in = pts_(:);

            [eor,ndep] = deal(eor_,ndep_);
            meta = LD_observations_set.make_meta_data(eor,ndep);
            npts = length(pts_in)/(meta.ndim);

            if (isempty(dnp1xu_in))
                dnp1xu_in = zeros(ndep*npts,1);
            end

            %% initializations

            obj.eor = eor;
            obj.ndep = ndep;

            obj.pts_in = pts_in;
            obj.dnp1xu_in = dnp1xu_in;
            obj.JFs_in = JFs_in;

            obj.npts = npts;
            obj = obj.load_jets;
        end
        function obj_out = load_jets(obj)
            obj_out = obj;
            obj_out.jets = obj_out.make_jets;
        end
        function jets_out = make_jets(obj,kor_)
            if (nargin == 1)
                if (prod(double(obj.dnp1xu_in == 0))) % if using zeros as placeholders
                    kor = obj.eor;
                else
                    kor = obj.eor+1;
                end
            else
                kor = kor_;
            end
            [eor,ndep,ndim,npts] = obj.get_meta_dims;
            pts_mat = reshape(obj.pts_in,ndim,obj.npts);

            if (kor>eor)
                pts_mat = [pts_mat ; reshape(obj.dnp1xu_in,ndep,[])];
            else
                pts_mat = pts_mat(1:(1+ndep*(kor+1)),:);
            end

            jets_out = LD_jets(pts_mat,ndep);
        end
        function u_hat_out = u_hat(obj)
            [eor,ndep,ndim] = obj.get_meta_dims;
            np1 = double(eor+1);
            o_range = (0:np1)';

            pts_mat = reshape(obj.pts_in,ndim,[]);
            nptsm1 = size(pts_mat,2)-1;

            x_vec = pts_mat(1,:);
            u_np1_mat = [pts_mat(2:end,:) ; reshape(obj.dnp1xu_in,ndep,[])];

            sgn_mat = ones(ndep,np1+1);
            sgn_mat(:,2:2:end) = -1;
            sgn_vec = sgn_mat(:);

            U_n_vec = 0.5*(u_np1_mat(:,1:(end-1)) + (sgn_vec).*u_np1_mat(:,2:end));
            U_n_tns = reshape(U_n_vec,ndep,np1+1,nptsm1);
            h_vec = 0.5*(x_vec(2:end) - x_vec(1:(end-1)));
            h_mat = ([1; cumprod(o_range(2:end))]).*(h_vec.^(o_range));

            u_hat_out = nan(2*ndep,nptsm1);
            for i = 1:nptsm1
                u_hat_out(:,i) = U_n_tns(:,:,i)*h_mat(:,i);
            end
        end
        function tau_mat_out = tau_mat(obj)
            [eor,ndep,ndim] = obj.get_meta_dims;
            pts_mat = reshape(obj.pts_in,ndim,[]);
            dnp1xu_mat = reshape(obj.dnp1xu_in,ndep,[]);
            tau_mat_out = [ ones(1,obj.npts) ; pts_mat((ndep+2):end,:) ; dnp1xu_mat ];
        end
        function [eor,ndep,ndim,npts] = get_meta_dims(obj)
            [eor,ndep] = deal(obj.eor,obj.ndep);
            ndim = 1 + ndep*(eor+1);
            if (nargout == 4)
                npts = obj.npts;
            end
        end
    end
    methods (Static)

    end
end
