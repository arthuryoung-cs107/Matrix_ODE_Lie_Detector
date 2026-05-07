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
% sys_screens = get(groot,'MonitorPositions');

screen_specs = sys_screens(1,:);
[o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1)

% scrn_id = 1;
scrn_id = 3;

%{
    ----------------------------
    load/generate observations
    ----------------------------
%}
% load a data set

% generate a data set
tic0 = tic;
% [Sobs,dat_true] = ldaux.generate_Riccati_data();
[Sobs,dat_true] = ldaux.generate_Brusselator_data();
toc1 = toc(tic0);
fprintf('generated jet space data in %.3f seconds \n', toc1);

% dat = dat_true; % oracle information, debugging purposes
dat = struct('ndep',dat_true.ndep,'eor',dat_true.eor);

%{
    ----------------------------
    plot observations
    ----------------------------
%}
dat_plt0 = dat;
dat_plt0.LineStyle = '-';
dat_plt0.Color = apv_plots.green4;
plt0 = apv_plots('jetspace', ...
                [3 4],...
                [1 3],[1 1], ...
                scrn_id);
plt0 = apv_plots.plot_Sobs(plt0,Sobs,dat_plt0);

%{
    ----------------------------
    build jet space model
    ----------------------------
%}
fspc = struct( ...
'bor', 3 ...
);
tic0 = tic;
mod = mvp_jspc_model(Sobs,dat,fspc)
toc1 = toc(tic0);
fprintf('built jet space model in %.3f seconds \n', toc1);

mod = mvp_jspc_model.verify(mod);

dat_plt1 = dat_plt0;
% dat_plt1.LineStyle = '-';
% dat_plt1.Color = apv_plots.green4;
plt1 = apv_plots('svds', ...
                [3 4],...
                [1 3],[3 1], ...
                scrn_id);
% plt1 = apv_plots.plot_SVDs(plt1,mod.Rsvd_cell,dat_plt1);
plt1 = apv_plots.plot_SVDs(plt1,vertcat(mod.Rsvd_cell(:),mod.Rksvds(:),mod.Rk_sub_svds(:)),dat_plt1);
