load ens_01_step04_for.txt
load ens_02_step04_for.txt
load ens_03_step04_for.txt
load ens_04_step04_for.txt

load ens_01_step04_ana.txt
load ens_02_step04_ana.txt
load ens_03_step04_ana.txt
load ens_04_step04_ana.txt

load state_step04_for.txt
load state_step04_ana.txt
load ../inputs_online/true_step4.txt
load ../inputs_online/obs_step4.txt
load ../inputs_online/true_initial.txt
load ../inputs_online/state_ini.txt


f1=figure;
subplot(2,2,1)
pcolor(state_step04_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'clim',[-1 1])
title('Forecast','fontsize',20)
subplot(2,2,2)
pcolor(state_step04_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Analysis','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,4)
pcolor(true_step4)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('True field','fontsize',20)
set(gca,'clim',[-1 1])
set(gcf,'position',[1,91,946, 878])

for i=1:12
    for j=1:12
        if obs_step4(i,j)==-999.0
            obs_step4(i,j) = NaN;
        end
    end
end
%f4=figure;
subplot(2,2,3)
pcolor(obs_step4)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Observations','fontsize',20)
set(gca,'clim',[-1 1])
rmse_for=0;
rmse_ana=0;
for i=1:12
    for j=1:12
    rmse_for = rmse_for + (state_step04_for(i,j) - true_step4(i,j))^2;
    rmse_ana = rmse_ana + (state_step04_ana(i,j) - true_step4(i,j))^2;
    end    
end
rmse_for = sqrt(rmse_for/144)
rmse_ana = sqrt(rmse_ana/144)

f2=figure;
subplot(2,2,1)
pcolor(ens_01_step04_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 1','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,2)
pcolor(ens_02_step04_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 2','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,3)
pcolor(ens_03_step04_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 3','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,4)
pcolor(ens_04_step04_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 4','fontsize',20)
set(gca,'clim',[-1 1])
set(gcf,'position',[300,91,946, 878])


f3=figure;
subplot(2,2,1)
pcolor(ens_01_step04_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 1','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,2)
pcolor(ens_02_step04_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 2','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,3)
pcolor(ens_03_step04_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 3','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,4)
pcolor(ens_04_step04_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 4','fontsize',20)
set(gca,'clim',[-1 1])
set(gcf,'position',[600,91,946, 878])

f4=figure;
pcolor(state_ini)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Initial ensemble mean','fontsize',20)
set(gca,'clim',[-1 1])
set(gcf,'position',[200,100,408,369])

f5=figure;
pcolor(true_initial)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('True initial state','fontsize',20)
set(gca,'clim',[-1 1])
set(gcf,'position',[1000,100,408, 369])

f6=figure;
pcolor(true_initial)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('True initial state','fontsize',20)
set(gca,'fontsize',20)
set(gca,'clim',[-1 1])
colorbar
set(gcf,'position',[500,1,500,500])
