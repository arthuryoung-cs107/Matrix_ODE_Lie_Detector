%{
    ----------------------------
    matlab preliminaries
    ----------------------------
%}
clear;
close all;

%{
    ----------------------------
    plotting contexts
    ----------------------------
%}
sys_screens = apv_plots.get_sys_screens();
screen_specs = sys_screens(1,:);
[o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1)
scrn_id = 1;

%{
    ----------------------------
    load/generate observations
    ----------------------------
%}

%{
    ----------------------------
    load/generate observations
    ----------------------------
%}
% % generate a data set
tic0 = tic;
[Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Linden_bouyancy_data(); % N = 2, Q = 2, linear homogenous
toc1 = toc(tic0);
fprintf('generated jet space data in %.3f seconds \n', toc1);

ncrv = length(Sobs(:));
npts_per_crv = zeros(ncrv,1);
for i = 1:ncrv
    npts_per_crv(i) = size(Sobs{i},2);
end
icrv_check = 1;
isol_check = 10;
iisol_check = sum(npts_per_crv(1:(icrv_check-1))) + isol_check;

dat = struct('ndep',dat_true.ndep,'eor',dat_true.eor);

ndep = dat.ndep;
eor = dat.eor;
nvar = 1+ndep;
ndim = 1+ndep*(eor + 1);
nvar_N1 = ndim - ndep; % = 1+Q(N+1) - Q = 1 + QN

[Smat,nobs,ncrv,kor,ndim,npts_per_crv,ipts_crv] = ldaux.unpack_Scell(Sobs,ndep);
Xspan = max(Smat(1,:))-min(Smat(1,:));
%{
    ----------------------------
    plot observations
    ----------------------------
%}
plt0 = apv_plots('Linden_buoyancy_jetspace', ...
                [5 6],...
                [2 6],[3 1], ...
                scrn_id);
plt0 = plt0.init_tiles_safe(1,4);
axs = plt0.axs;
axs_mat = plt0.axs_mat;
hold(axs, 'on');
box(axs,'on');
dat_plt0 = dat;
dat_plt0.LineStyle = '-';
dat_plt0.MarkerSize = 4;
dat_plt0.Color = apv_plots.green4;
plt0 = apv_plots.plot_Sobs(plt0,Sobs,dat_plt0);
dat_plt_i = dat_plt0;
dat_plt_i.LineStyle = 'none';
dat_plt_i.Color = apv_plots.blue1;
dat_plt_i.MarkerSize = 6;
plt0 = apv_plots.plot_Sobs(plt0,Sobs{1},dat_plt_i);

plot3( axs(4), ...
    Smat(1,:), Smat(2,:) , Smat(3,:) , ...
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'Color', dat_plt0.Color, ...
    'Marker', 'o', ...
    'MarkerSize', dat_plt0.MarkerSize, ...
    'MarkerFaceColor', dat_plt0.Color, ...
    'MarkerEdgeColor', dat_plt0.Color ...
    );
for i = 1:ncrv
    plot3( axs(4), ...
        Sobs{i}(1,:), Sobs{i}(2,:) , Sobs{i}(3,:) , ...
        'LineStyle', '-', ...
        'LineWidth', 0.5, ...
        'Color', dat_plt0.Color, ...
        'Marker', 'none' );
end
plot3( axs(4), ...
    Sobs{1}(1,:), Sobs{1}(2,:) , Sobs{1}(3,:) , ...
    'Color', dat_plt_i.Color, ...
    'Marker', 'o', ...
    'MarkerSize', dat_plt_i.MarkerSize, ...
    'MarkerFaceColor', dat_plt_i.Color, ...
    'MarkerEdgeColor', dat_plt_i.Color ...
    );
xlabel(plt0.axs(1:3), '$$ t \textrm{ [s]} $$', 'Interpreter','Latex','FontSize',16);
ylabel(plt0.axs(1), '$$ h \textrm{ [m]} $$', 'Interpreter','Latex','FontSize',16);
ylabel(plt0.axs(2), '$$ \dot h \textrm{ [m/s]} $$', 'Interpreter','Latex','FontSize',16);
ylabel(plt0.axs(3), '$$ \ddot h \textrm{ [m/s} ^2 \textrm{]} $$', 'Interpreter','Latex','FontSize',16);

xlabel(plt0.axs(4), '$$ t \textrm{ [s]} $$', 'Interpreter','Latex','FontSize',16);
ylabel(plt0.axs(4), '$$ h \textrm{ [m]} $$', 'Interpreter','Latex','FontSize',16);
zlabel(plt0.axs(4), '$$ \dot h \textrm{ [m/s]} $$', 'Interpreter','Latex','FontSize',16);
% zlabel(plt0.axs(4), '$$ \ddot h \textrm{ [m/s} ^2 \textrm{]} $$', 'Interpreter','Latex','FontSize',16);
view(axs(4), apv_plots.view_mat(6, :));
plt0.show_toolbar

%{
    ----------------------------
    specify function space
    ----------------------------
%}
fmap0 = struct( ...
'Omap_a', ones(nvar,1), ...
'Omap_b', zeros(nvar,1) ...
);
fmap = fmap0;

%{
    ----------------------------
    model solution space
    ----------------------------
%}
tic0 = tic;
mod = LDsol.model_solspace(Sobs,dat,fmap);
toc1 = toc(tic0);
fprintf('built jet space model in %.3f seconds \n', toc1);

mod
sO = mod.s_O

Jf_RN1_check = mod.J_tau_u_RN1(:,:,1,iisol_check)
Jdxf_RN1_check = mod.J_tau_u_RN1(:,:,2,iisol_check)
JF_tru_check = JF_obs{icrv_check}( :,:,isol_check)
JF_N1_check = mod.JF_N1( :,:, iisol_check )

dNp1xu_tru_check = dNp1xu_obs{icrv_check}( :,isol_check )
dNp1xu_N1mod_check = mod.tau_uN_RN1_net((end-ndep+1):end,2,iisol_check)

function out = plot_tvector(axi_,sx_,su_,vx_,vu_,color_,LS_,mrkr_)
    plot(axi_, ...
    [sx_ , sx_+vx_], [su_ , su_+vu_], ...
    'LineStyle', LS_, ...
    'LineWidth', 1, ...
    'Color', color_, ...
    'Marker', mrkr_, ...
    'MarkerSize', 3, ...
    'MarkerFaceColor', color_, ...
    'MarkerEdgeColor', color_ ...
    );
end
plt_3D = @(axi_,s_,mrkr_,c_,LS_) plot3(axi_, ...
    s_(1,:), s_(2,:), s_(3,:), ...
    'LineStyle', LS_, ...
    'LineWidth', 1, ...
    'Color', c_, ...
    'Marker', mrkr_, ...
    'MarkerSize', 3, ...
    'MarkerFaceColor', c_, ...
    'MarkerEdgeColor', c_ ...
);

tvf_scl = 0.25*Xspan;
[xO,uO] = deal( mod.s_O(1),reshape(mod.sNp1_O(2:end),ndep,kor+2) );
tvfO_x = tvf_scl;
tvfO_u = tvf_scl*uO(:,2:end);
tvfO_s = [tvfO_x ; tvfO_u(:) ];
Vspc_O = mod.GN1_sO_basis.Vspc_sO;
Vspc_O = Vspc_O*( norm(tvfO_s) / sqrt(max(sum(Vspc_O.*Vspc_O,1))) ); % rescale wrt plotting tvf magnitude
Vspc_O_x = Vspc_O(1,:);
Vspc_O_u = reshape( Vspc_O(2:end,:), ndep,kor+1,[] );
cmat_i = hsv(nvar_N1);
for k = 1:(kor+1)
    plot(axs(k), ...
    xO, uO(1,k), ...
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'Color', [1 1 1], ...
    'Marker', 'o', ...
    'MarkerSize', 3, ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [1 1 1] ...
    );
    plot_tvector(axs(k),xO,uO(1,k),tvfO_x,tvfO_u(1,k),[1 1 1],'-','none');
    for iq = 1:nvar_N1
        plot_tvector(axs(k),xO,uO(1,k),Vspc_O_x(iq),Vspc_O_u(1,k,iq),cmat_i(iq,:),'-','none');
    end
end
plt_3D(axs(4), mod.s_O(1:3),'o',[1 1 1],'none')
plt_3D(axs(4),[ mod.s_O(1:3), mod.s_O(1:3)+tvfO_s(1:3) ],'none',[1 1 1],'-')
for iq = 1:nvar_N1
    plt_3D(axs(4),[ mod.s_O(1:3), mod.s_O(1:3)+Vspc_O(1:3,iq) ],'none',cmat_i(iq,:),'-')
end
% plt0.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])

%% visualize solution space in canonical coordinates
plt1 = apv_plots('Linden_buoyancy_canon_coord_space', ...
                [5 6],...
                [2 6],[3 1], ...
                scrn_id);
plt1 = plt1.init_tiles_safe(1,4);
axs = plt1.axs;
axs_mat = plt1.axs_mat;
hold(axs, 'on');
box(axs,'on');
dat_plt1 = dat;
dat_plt1.eor = 0;
dat_plt1.ndep = nvar_N1;
dat_plt1.LineStyle = '-';
dat_plt1.MarkerSize = 4;
dat_plt1.Color = apv_plots.green4;
Sobs_Xi = mod.dXi_S_sO_cell;
for i = 1:ncrv
    Sobs_Xi{i} = [Sobs{i}(1,:) ; Sobs_Xi{i}];
end
[plt1,dat_plt1] = apv_plots.plot_Sobs(plt1,Sobs_Xi,dat_plt1);
Smat_Xi = mod.GN1_sO_basis.Xi_S - mod.GN1_sO_basis.Xi_sO;
plot3( axs(4), ...
    Smat_Xi(1,:), Smat_Xi(2,:) , Smat_Xi(3,:) , ...
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'Color', dat_plt1.Color, ...
    'Marker', 'o', ...
    'MarkerSize', dat_plt1.MarkerSize, ...
    'MarkerFaceColor', dat_plt1.Color, ...
    'MarkerEdgeColor', dat_plt1.Color ...
    );
xlabel(plt1.axs(1:3), '$$ t \textrm{ [s]} $$', 'Interpreter','Latex','FontSize',16);
ylabel(plt1.axs(1), '$$ \xi_1 $$', 'Interpreter','Latex','FontSize',16);
ylabel(plt1.axs(2), '$$ \xi_2 $$', 'Interpreter','Latex','FontSize',16);
ylabel(plt1.axs(3), '$$ \xi_3 $$', 'Interpreter','Latex','FontSize',16);

xlabel(plt1.axs(4), '$$ \xi_1 $$', 'Interpreter','Latex','FontSize',16);
ylabel(plt1.axs(4), '$$ \xi_2 $$', 'Interpreter','Latex','FontSize',16);
zlabel(plt1.axs(4), '$$ \xi_3 $$', 'Interpreter','Latex','FontSize',16);

view(axs(4), apv_plots.view_mat(6, :));
% plt1.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])
plt1.show_toolbar

function plot_labelled_svds(axi_,svds_,lbls_,leg_loc_)
    nsvd = length(svds_(:));
    cmat = hsv(nsvd);
    if (nsvd==1)
        leg_i(1) = plot(axi_, ...
        1:length(svds_.s), svds_.s, ...
        'LineStyle', 'none', ...
        'LineWidth', 0.5, ...
        'Color', cmat(1,:), ...
        'Marker', 's', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', cmat(1,:), ...
        'MarkerEdgeColor', cmat(1,:), ...
        'DisplayName', lbls_{1} ...
        );
    else
        for i = 1:nsvd
            leg_i(i) = plot(axi_, ...
            1:length(svds_{i}.s), svds_{i}.s, ...
            'LineStyle', 'none', ...
            'LineWidth', 0.5, ...
            'Color', cmat(i,:), ...
            'Marker', 's', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', cmat(i,:), ...
            'MarkerEdgeColor', cmat(i,:), ...
            'DisplayName', lbls_{i} ...
            );
        end
    end
    set(axi_, ...
        'YScale', 'log' );
    ylabel(axi_, '$ \sigma_i $', 'Interpreter','Latex','FontSize',16);
    legend(axi_, leg_i(1:nsvd),'Location', leg_loc_, 'Interpreter', 'Latex', 'NumColumns',1,'FontSize',14);
end
rlbl = @(s_) [ ', \textrm{rank} =' num2str(s_.r) '/' num2str(s_.dim) '\textrm{ columns}'];
nlbl = @(s_) [ ', \| \cdot \|_* =' num2str(sum(s_.s/s_.s(1)),'%.1f') '/' num2str(s_.dim) ];
clbl = @(s_) [ ', \sigma_{\textrm{max}} / \sigma_{\textrm{min}} =' num2str(s_.s(1)/s_.s(end),'%.1e') ];

dat_plt2 = dat_plt0;
plt2 = apv_plots('Linden_buoyancy_model_summary', ...
                [5 6],...
                [2 6],[5 1], ...
                scrn_id);
plt2 = plt2.init_tiles_safe(1,2);
hold(plt2.axs, 'on');
box(plt2.axs,'on');
axs = plt2.axs;
axs_mat = plt2.axs_mat;
set(axs, ...
    'XScale', 'linear',  ...
    'TickLabelInterpreter','Latex', ...
    'FontSize',16 );
xlabel(axs, '$$ i $$', 'Interpreter','Latex','FontSize',16);
plt2.show_toolbar

axi = axs_mat(1,1);
svds_i = { ...
    mod.Rsvd_N1_net, '$R'; ...
    mod.Gsvd_N1_net,  '$G'; ...
};
svds = cell([ size(svds_i,1),1 ]);
labels = cell([length(svds),1]);
for isvd = 1:size(svds_i,1)
    labels{isvd} = [ svds_i{isvd,2} rlbl(svds_i{isvd,1}) clbl(svds_i{isvd,1}) '$'];
    svds{isvd} = svds_i{isvd,1};
end
plot_labelled_svds(axi,svds,labels,'NorthEast');

axi = axs_mat(1,2);
svds_i = { ...
mod.GN1_sO_basis , '$\Lambda^{(0)} |_{s_O} W_G'; ...
mod.GN1_sNO_basis, '$\Lambda^{(1)} |_{s_O} W_G'; ...
};
svds = cell([ size(svds_i,1),1 ]);
labels = cell([length(svds),1]);
for isvd = 1:size(svds_i,1)
    labels{isvd} = [ svds_i{isvd,2} rlbl(svds_i{isvd,1}) clbl(svds_i{isvd,1}) '$'];
    svds{isvd} = svds_i{isvd,1};
end
plot_labelled_svds(axi,svds,labels,'SouthWest');
% plt2.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])
% set(axs(2),'YScale', 'linear')
