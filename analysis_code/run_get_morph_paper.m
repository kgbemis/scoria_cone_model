% script to run get_volcano_morph

dirname='../data_for_DistVolcPaper/';
% build_ver code filename
file1='httrials_volcano_vel30_NppP100000_NP200_gs0.1_annulus0dir360std55dir10_vL0_rv3d_01-03-21_194350.mat';
code1='httrials_194350';
file2='httrials_volcano_vel40_NppP100000_NP200_gs0.1_annulus0dir360std55dir10_vL0_rv3d_10-04-21_181259.mat';
code2='httrials_181259';
file3='httrials_volcano_vel50_NppP100000_NP200_gs0.1_annulus0dir360std55dir10_vL0_rv3d_22-02-21_174838.mat';
code3='httrials_174838';
% duplicate runs used in original plot for shape with growth
%file4='httrials_volcano_vel40_NppP100000_NP200_gs0.1_annulus0dir360std55dir10_vL0_rv3d_26-03-21_165117.mat';
%code4='httrials_165117';
%file5='httrials_volcano_vel30_NppP100000_NP200_gs0.1_annulus0dir360std55dir10_vL0_rv3d_09-04-21_170216.mat';
%code5='httrials_170216';

load(fullfile(dirname,file1))
filecode=code1;

stage=200;
%zz=stages(:,:,stage);
zz=stages;
[m,n,~]=size(stages);
xx=1:m;  % these might be opposite if n~=m
yy=1:n;
get_volcano_morph_auto_from_fig_paper(xx,yy,zz,stage,10,filecode)

%%

openfig('code1_morph_profiles.fig')
cumfall=load('../diagnostics_code/pro_200_fall.mat','pro_200_fall');
hold on
    plot(1:1000,pro_200_fall,'-m','LineWidth',2)
    daspect([1 1 1])
hold off
axis([0 1000 0 300])
title('')
%print('cummulative_fall_combined_cone.png','-dpng','-r300')
