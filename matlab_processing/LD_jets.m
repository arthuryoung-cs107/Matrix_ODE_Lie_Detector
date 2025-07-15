classdef LD_jets
    properties
        pts_mat;
        ndep;

        pts;
        pts0;
        pts1;

        kor; % order of input pt. ders,  which induce kor + 1 constraints per dim. per pt.
        jor; % order of induced jet polynomial, generally equal to (2*(kor + 1) - 1)*m

        xh_vec; % origin of midpoint jets
        hx_vec; % radius of midpoint jets
        a_tns; % jet coefficients at colocation points

        %% legacy
        xh; % origin of midpoint jet
        hx; % radius of midpoint jet

        xh0;
        xh1;

        Umat;
        Vmnd;
        a_mat;

    end
    methods
        function obj = LD_jets(pts_,ndep_)
            % pts = pts_(:,1:2);
            pts = [pts_(:,1), pts_(:,end)];

            [pts0,pts1,ndep] = deal(pts(:,1),pts(:,end),ndep_);

            korp1 = (size(pts,1) - 1)/ndep;
            kor = korp1-1;
            krange = double((0:kor));

            jorp1 = 2*korp1; % given 2*(kor+1) constraints in each dimension,
            jor = jorp1 - 1; % determines coefficients of order 2*(kor+1)-1 jet in each dim.
            jrange = double((0:jor));

            %% basic data management
            pts_mat = pts_;
            npts = size(pts_mat,2);
            nptsm1 = npts-1;
            pts_mat_0 = pts_mat(:,1:(end-1));
            pts_mat_1 = pts_mat(:,2:end);
            x_vec = pts_mat(1,:);

            xh_vec = 0.5*(x_vec(2:end)+x_vec(1:(end-1))); % origin of midpoint jet
            hx_vec = 0.5*(x_vec(2:end)-x_vec(1:(end-1))); % radius of midpoint jet
            hx_0_P = (-hx_vec).^(jrange');
            hx_1_P = (hx_vec).^(jrange');

            U_tns = reshape(pts_mat(2:end,:),ndep,korp1,npts);
            U0_tns = U_tns(:,:,1:(end-1));
            U1_tns = U_tns(:,:,2:end);

            cVmnd = [ones(1,jorp1) ; zeros(kor,jorp1)];
            for iord = 1:kor
                cVmnd(iord+1,(iord+1):end) = jrange(2:(end-iord+1)).*cVmnd(iord,(iord+1):end);
            end
            % cVmnd

            a_tns = nan(ndep,jorp1,nptsm1);
            for ipts = 1:nptsm1 % solve linear system for each curve point solution and dim
                % dVmnd_i = [ cVmnd.*(hx_0_P(:,ipts)') ; cVmnd.*(hx_1_P(:,ipts)') ];
                [dVmnd_i0,dVmnd_i1] = deal(cVmnd);
                dVmnd_i0(1,:) = dVmnd_i0(1,:).*hx_0_P(:,ipts)';
                dVmnd_i1(1,:) = dVmnd_i1(1,:).*hx_1_P(:,ipts)';
                for iord = 2:korp1
                    dVmnd_i0(iord,iord:end) = dVmnd_i0(iord,iord:end).*hx_0_P(1:(end-iord+1),ipts)';
                    dVmnd_i1(iord,iord:end) = dVmnd_i1(iord,iord:end).*hx_1_P(1:(end-iord+1),ipts)';
                end
                dVmnd_i = [dVmnd_i0 ; dVmnd_i1];
                for idep = 1:ndep
                    a_tns(idep,:,ipts) = dVmnd_i\[ U0_tns(idep,:,ipts)' ; U1_tns(idep,:,ipts)' ];
                end
            end

            obj.pts_mat = pts_mat;
            obj.ndep = ndep;
            obj.xh_vec = xh_vec;
            obj.hx_vec = hx_vec;
            obj.a_tns = a_tns;
            %% Hermite polynomial generation
            xh = 0.5*(pts1(1)+pts0(1)); % origin of midpoint jet
            hx = 0.5*(pts1(1)-pts0(1)); % radius of midpoint jet

            xh0 = (-hx).^(jrange);
            xh1 = (hx).^(jrange);

            U0mat = reshape(pts0(2:end,:),ndep,korp1)';
            U1mat = reshape(pts1(2:end,:),ndep,korp1)';
            Umat = [U0mat ; U1mat];

            Vmnd1 = zeros(2,jorp1,korp1);
            coeffs = ones(size(jrange));
            Vmnd1(:,:,1) = [xh0 ; xh1];
            for ior = 1:kor % walking through derivative constraints
                coeffs = coeffs(2:end).*jrange(2:(end-ior+1)); % update coeffs, skip constant.
                Vmnd1(:,(ior+1):end,ior+1) = coeffs.*[xh0(1:(end-ior)) ; xh1(1:(end-ior))];
            end
            Vmnd1_01 = permute(Vmnd1,[3 2 1]);
            Vmnd = [Vmnd1_01(:,:,1) ; Vmnd1_01(:,:,2)];

            a_mat = nan(jorp1,ndep);
            for idep = 1:ndep
                a_mat(:,idep) = Vmnd\Umat(:,idep);
            end


            %% initializations

            obj.pts = pts;
            obj.pts0 = pts0;
            obj.pts1 = pts1;

            obj.kor = kor;
            obj.jor = jor;

            obj.xh = xh;
            obj.hx = hx;

            obj.xh0 = xh0;
            obj.xh1 = xh1;

            obj.Umat = Umat;
            obj.Vmnd = Vmnd;
            obj.a_mat = a_mat;
        end
        function [u_out,dxu_out] = u_hat(obj, x_, kor_)
            if (nargin==3)
                kor = kor_;
            else
                kor = obj.kor;
            end
            xvec = reshape(x_,[],1);
            [ndep,npts,xh_vec,jor] = deal(obj.ndep,length(xvec),obj.xh_vec,obj.jor);
            [~, idx_h] = min(abs( xvec - xh_vec ),[],2);

            sh_mat_i = obj.pts_mat(:,idx_h);
            xh_vec_i = xh_vec(idx_h);
            a_tns_i = obj.a_tns(:,:,idx_h);
            h_mat_i_P = (xvec'-xh_vec(idx_h)).^double((0:jor)');

            if (nargout==1)
                u_out = reshape(pagemtimes(a_tns_i,reshape(h_mat_i_P,[],1,npts)),ndep,npts);
            else
                korp1 = kor+1;
                jor = obj.jor;
                [jorp1,jrange] = deal(jor+1, double(0:jor));
                cVmnd = [ones(1,jorp1) ; zeros(kor,jorp1)];
                for iord = 1:kor
                    cVmnd(iord+1,(iord+1):end) = jrange(2:(end-iord+1)).*cVmnd(iord,(iord+1):end);
                end

                u_out = nan(ndep,npts);
                dxu_out = nan(ndep,kor,npts);
                for ipts = 1:npts
                    dVmnd_i = cVmnd;
                    dVmnd_i(1,:) = dVmnd_i(1,:).*h_mat_i_P(:,ipts)';

                    for iord = 2:korp1
                        dVmnd_i(iord,iord:end) = dVmnd_i(iord,iord:end).*h_mat_i_P(1:(end-iord+1),ipts)';
                    end

                    uk_i = (dVmnd_i)*(a_tns_i(:,:,ipts)');
                    u_out(:,ipts) = (uk_i(1,:))';
                    dxu_out(:,:,ipts) = uk_i(2:end,:)';

                end
                dxu_out = reshape(dxu_out,[],npts);
            end
        end
    end
    methods (Static)
        function jet_out = regularize_jet(jet_,sigma_,lam_)
            if (nargin==1)
                sigma = ones(jet_.kor+1,jet_.ndep);
                lam = ones(jet_.ndep,1);
            else
                if (nargin==2)
                    sigma = sigma_;
                    lam = ones(jet_.ndep,1);
                else
                    [sigma,lam] = deal(sigma_,lam_);
                end
            end
            ndep = jet_.ndep;
            kor = jet_.kor; % order of observed derivatives
            kp1 = kor+1; % order of applied regularization
            kp1_range = double(0:kp1);
            M = jet_.jor+1; % order of regularized jet, one greater than exact jet
            Mp1 = M+1; % coefficient count of regularized jet, one greater than order (due to k = 0)

            M_range = double(0:M);
            Fkp1_mat_kp1 = [ ones(1,Mp1) ; zeros(kp1,Mp1) ];
            for i = 1:kp1 % for evaluation of derivatives up to order k+1
                Fkp1_mat_kp1(i+1,(i+1):end) = Fkp1_mat_kp1(i,(i+1):end) .* M_range(2:(end-i+1));
            end
            Fkp1_mat_k = Fkp1_mat_kp1(1:(end-1),:);
            Fkp1_vec = Fkp1_mat_kp1(end, (kp1+1):end);

            [xh_vec,hx_vec] = deal(jet_.xh_vec,jet_.hx_vec);
            hx_mat_Mp1 = hx_vec.^( [M_range, double(Mp1)]' );
            sgn0_hx_Mp1 =  ((-1).^[M_range, double(Mp1)]');
            [hx_mat0_Mp1,hx_mat1_Mp1] = deal(sgn0_hx_Mp1.*hx_mat_Mp1,hx_mat_Mp1);
            [hx_mat0,hx_mat1] = deal(hx_mat0_Mp1(1:(end-1),:),hx_mat1_Mp1(1:(end-1),:));
            nh = length(xh_vec);

            U_0 = reshape(jet_.pts_mat(2:end,1:(end-1)),ndep,kp1,nh);
            U_1 = reshape(jet_.pts_mat(2:end,2:end),ndep,kp1,nh);
            ismat = sigma.^(-2);
            istns = reshape(ismat,kp1,1,ndep);

            a_tns = zeros(ndep,Mp1,nh);
            for ipts = 1:nh
                [U0_i,U1_i] = deal(U_0(:,:,ipts)',U_1(:,:,ipts)');
                [dVmnd_i0,dVmnd_i1] = deal(Fkp1_mat_k);
                [dVmnd_i0(1,:),dVmnd_i1(1,:)] =  deal((hx_mat0(:,ipts)'),(hx_mat1(:,ipts)'));
                [Dtns0_i,Dtns1_i] = deal(zeros(Mp1,Mp1,ndep));
                [Dtns0_i(1,:,:),Dtns1_i(1,:,:)] = deal( istns(1,:,:).*dVmnd_i0(1,:), ...
                                                        istns(1,:,:).*dVmnd_i1(1,:) );
                [bmat0_i,bmat1_i] = deal(zeros(Mp1,ndep));
                [bmat0_i(1,:),bmat1_i(1,:)] = deal( ismat(1,:).*(U_0(:,1,ipts)'), ...
                                                    ismat(1,:).*(U_1(:,1,ipts)') );
                for iord = 2:kp1
                    dVmnd_i0(iord,iord:end) = dVmnd_i0(iord,iord:end).*(hx_mat0(1:(end-iord+1),ipts)');
                    dVmnd_i1(iord,iord:end) = dVmnd_i1(iord,iord:end).*(hx_mat1(1:(end-iord+1),ipts)');
                    [Dtns0_i(iord,:,:),Dtns1_i(iord,:,:)] = deal(   istns(1,:,:).*dVmnd_i0(1,iord)*dVmnd_i0(1,:), ...
                                                                    istns(1,:,:).*dVmnd_i1(1,iord)*dVmnd_i1(1,:) );
                    [bmat0_i(iord,:),bmat1_i(iord,:)] = deal(   ismat(1,:).*dVmnd_i0(1,iord)*(U_0(:,1,ipts)'), ...
                                                                ismat(1,:).*dVmnd_i1(1,iord)*(U_1(:,1,ipts)') );
                    for l = 2:iord
                         Dtns0_i(iord,:,:) = Dtns0_i(iord,:,:) + istns(l,:,:).*dVmnd_i0(l,iord)*dVmnd_i0(l,:);
                         Dtns1_i(iord,:,:) = Dtns1_i(iord,:,:) + istns(l,:,:).*dVmnd_i1(l,iord)*dVmnd_i1(l,:);
                         bmat0_i(iord,:) = bmat0_i(iord,:) + ismat(l,:).*dVmnd_i0(l,iord)*(U_0(:,l,ipts)');
                         bmat1_i(iord,:) = bmat1_i(iord,:) + ismat(l,:).*dVmnd_i1(l,iord)*(U_1(:,l,ipts)');
                    end
                end
                Imat_i = zeros(Mp1,Mp1);
                for iord = (kp1+1):Mp1
                    [Dtns0_i(iord,:,:),Dtns1_i(iord,:,:)] = deal(   istns(1,:,:).*dVmnd_i0(1,iord)*dVmnd_i0(1,:), ...
                                                                    istns(1,:,:).*dVmnd_i1(1,iord)*dVmnd_i1(1,:) );
                    [bmat0_i(iord,:),bmat1_i(iord,:)] = deal(   ismat(1,:).*dVmnd_i0(1,iord)*(U_0(:,1,ipts)'), ...
                                                                ismat(1,:).*dVmnd_i1(1,iord)*(U_1(:,1,ipts)') );
                    for l = 2:kp1
                         Dtns0_i(iord,:,:) = Dtns0_i(iord,:,:) + istns(l,:,:).*dVmnd_i0(l,iord)*dVmnd_i0(l,:);
                         Dtns1_i(iord,:,:) = Dtns1_i(iord,:,:) + istns(l,:,:).*dVmnd_i1(l,iord)*dVmnd_i1(l,:);
                         bmat0_i(iord,:) = bmat0_i(iord,:) + ismat(l,:).*dVmnd_i0(l,iord)*(U_0(:,l,ipts)');
                         bmat1_i(iord,:) = bmat1_i(iord,:) + ismat(l,:).*dVmnd_i1(l,iord)*(U_1(:,l,ipts)');
                    end
                    i_range = double(iord-kp1) + kp1_range;
                    Imat_i(iord,(kp1+1):end) = Fkp1_vec(iord-kp1)*((Fkp1_vec.*(hx_mat1_Mp1(i_range+1,ipts)-hx_mat0_Mp1(i_range+1,ipts))')./i_range);
                end

                Dtns_i = Dtns0_i+Dtns1_i;
                bmat_i = bmat0_i+bmat1_i;

                % Dtns_i
                % pause

                for idep = 1:ndep
                    a_tns(idep,:,ipts) = [Dtns_i(:,:,idep)+lam(idep)*Imat_i] \ bmat_i(:,idep);
                end
            end
            jet_out = jet_;
            jet_out.jor = M;
            jet_out.a_tns = a_tns;

        end
    end
end
