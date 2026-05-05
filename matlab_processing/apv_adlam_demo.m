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

scrn_id = 1;
% scrn_id = 3;

% load a data set
[Sobs,dat_true] = ldaux.generate_Riccati_data();
% dat = struct('ndep',dat_true.ndep,'eor',dat_true.eor);
dat = dat_true;

dat_plt0 = dat;
dat_plt0.LineStyle = '-';
dat_plt0.Color = apv_plots.green4;
plt0 = apv_plots('S', ...
                [3 4],...
                [1 3],[1 1], ...
                1);
plt0 = apv_plots.plot_Sobs(plt0,Sobs,dat_plt0);
