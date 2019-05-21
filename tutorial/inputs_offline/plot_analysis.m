% Script to plot the analysis output file from the PDAF tutorial

dim_x = 36
dim_y = 18
dim_ens = 9

% Load and plot analysis state

load state_ana.txt

field_plot=zeros(dim_y+1, dim_x+1);
field_plot(1:dim_y,1:dim_x) = state_ana;
figure
pcolor(field_plot)
set(gca,'fontsize',16)
cb=colorbar
set(cb,'fontsize',16)
title('Analysis state','fontsize',18)
