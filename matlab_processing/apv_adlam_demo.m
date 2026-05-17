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

%{
    ----------------------------
    load/generate observations
    ----------------------------
%}
% load a data set


% generate a data set
tic0 = tic;
% [Sobs,dat_true] = ldaux.generate_Riccati_data(); % N = 1, Q = 1
[Sobs,dat_true] = ldaux.generate_Brusselator_data(); % N = 1, Q = 2
% [Sobs,dat_true] = ldaux.generate_Van_der_Pol_data(); % N = 2, Q = 1
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
% fspc = struct( ...
% 'bor', 3 ...
% );
fspc = struct( ...
'bor', 3 ...
);
tic0 = tic;
mod = mvp_jspc_model(Sobs,dat,fspc)
toc1 = toc(tic0);
fprintf('built jet space model in %.3f seconds \n', toc1);

%{
    ----------------------------
    review jet space model
    ----------------------------
%}
mod = mvp_jspc_model.verify(mod);
obj = mod;
ndep = obj.Sdat.ndep;
[ndim_obs,nobs] = size(obj.Smat);
kor_obs = (ndim_obs-1)/ndep - 1;

Smat = obj.Smat;
xvec = Smat(1,:);
unmat = Smat(2:end,:);
untns = reshape( unmat, ndep, kor_obs+1, nobs );
umat = reshape( untns(:,1,:), ndep, nobs );
dxutns = reshape( untns(:,2:end,:), ndep, kor_obs, nobs );
% minmag_dxutns = min(abs( dxutns(:) ));
% maxmag_dxutns = max(abs( dxutns(:) ));

tau_uk_glb = obj.tau_uk_glb(:,1:(end-1),:);
tau_uk_glb_sub = obj.tau_uk_glb_sub(:,1:(end-1),:,:);
% minmag_tau_uk_glb = min(abs( tau_uk_glb(:) ));
% maxmag_tau_uk_glb = max(abs( tau_uk_glb(:) ));
tau_uN_crv = obj.tau_uN_crv(:,1:(end-1),:);
tau_uN_crv_sub = obj.tau_uN_crv_sub(:,1:(end-1),:,:);

tmod_glb = struct('t_uk', tau_uk_glb, 't_uk_sub', tau_uk_glb_sub);
tmod_crv = struct('t_uk', tau_uN_crv, 't_uk_sub', tau_uN_crv_sub);
tmodels = [tmod_glb tmod_crv];

abserr_tol = 1e-2;
for imod = 1:length(tmodels)
    tmod_i = tmodels(imod);
    [tau_uk_i,tau_uk_i_sub] = deal(tmod_i.t_uk,tmod_i.t_uk_sub);

    err_uk = tau_uk_i-dxutns;
    abserr_uk = abs(err_uk);
    [abserr_uk_sorted,i_abserr_uk_sorted] = sort(abserr_uk(:));
    min_abserr_uk = abserr_uk_sorted(1);
    avg_abserr_uk = mean(abserr_uk_sorted);
    med_abserr_uk = median(abserr_uk_sorted);
    max_abserr_uk = abserr_uk_sorted(end);
    fprintf('(b=%d,full) abserr [min,avg,med,max]=[%.1e,%.1e,%.1e,%.1e]. Success = %d \n', ...
        fspc.bor, min_abserr_uk,avg_abserr_uk,med_abserr_uk,max_abserr_uk, ...
        max_abserr_uk<abserr_tol ...
    );
    err_uk_sub = nan(size(tau_uk_i_sub));
    for b = 1:fspc.bor
        err_uk_sub(:,:,:,b) = tau_uk_i_sub(:,:,:,b)-dxutns;

        abserr_uk_bsub = abs(err_uk_sub(:,:,:,b));
        [abserr_uk_sorted_bsub,i_abserr_uk_sorted_bsub] = sort(abserr_uk_bsub(:));
        min_abserr_uk_bsub = abserr_uk_sorted_bsub(1);
        avg_abserr_uk_bsub = mean(abserr_uk_sorted_bsub);
        med_abserr_uk_bsub = median(abserr_uk_sorted_bsub);
        max_abserr_uk_bsub = abserr_uk_sorted_bsub(end);
        fprintf('   (b=%d,sub) abserr [min,avg,med,max]=[%.1e,%.1e,%.1e,%.1e]. Success = %d \n', ...
            b, min_abserr_uk_bsub,avg_abserr_uk_bsub,med_abserr_uk_bsub,max_abserr_uk_bsub, ...
            max_abserr_uk_bsub<abserr_tol ...
        );
    end
end

%{
    ----------------------------
    plot jet space model data
    ----------------------------
%}
dat_plt1 = dat_plt0;
% dat_plt1.LineStyle = '-';
% dat_plt1.Color = apv_plots.green4;
plt1 = apv_plots('svds', ...
                [3 4],...
                [1 3],[3 1], ...
                scrn_id);
plt1 = apv_plots.plot_SVDs(plt1, mvp_jspc_model.vertcat_glb_SVDs(mod), dat_plt1);
