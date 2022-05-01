load ens_01_step02_for.txt
load ens_02_step02_for.txt
load ens_03_step02_for.txt
load ens_04_step02_for.txt


load ens_01_step02_for.txt
load ens_02_step02_for.txt
load ens_03_step02_for.txt
load ens_04_step02_for.txt

load ens_01_step02_ana.txt
load ens_02_step02_ana.txt
load ens_03_step02_ana.txt
load ens_04_step02_ana.txt

load state_step02_for.txt
load state_step02_ana.txt
load ../inputs_online/true_step2.txt
load ../inputs_online/obs_step2.txt
load ../inputs_online/true_initial.txt
load ../inputs_online/state_ini.txt


f1=figure;
subplot(2,2,1)
pcolor(state_step02_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'clim',[-1 1])
title('Forecast','fontsize',20)
subplot(2,2,2)
pcolor(state_step02_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Analysis','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,4)
pcolor(true_step2)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('True field','fontsize',20)
set(gca,'clim',[-1 1])
set(gcf,'position',[1,91,946, 878])

for i=1:12
    for j=1:12
        if obs_step2(i,j)==-999.0
            obs_step2(i,j) = NaN;
        end
    end
end
%f4=figure;
subplot(2,2,3)
pcolor(obs_step2)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Observations','fontsize',20)
set(gca,'clim',[-1 1])
rmse_for=0;
rmse_ana=0;
for i=1:12
    for j=1:12
    rmse_for = rmse_for + (state_step02_for(i,j) - true_step2(i,j))^2;
    rmse_ana = rmse_ana + (state_step02_ana(i,j) - true_step2(i,j))^2;
    end    
end
rmse_for = sqrt(rmse_for/144)
rmse_ana = sqrt(rmse_ana/144)

f2=figure;
subplot(2,2,1)
pcolor(ens_01_step02_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 1','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,2)
pcolor(ens_02_step02_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 2','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,3)
pcolor(ens_03_step02_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 3','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,4)
pcolor(ens_04_step02_for)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 4','fontsize',20)
set(gca,'clim',[-1 1])
set(gcf,'position',[300,91,946, 878])


f3=figure;
subplot(2,2,1)
pcolor(ens_01_step02_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 1','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,2)
pcolor(ens_02_step02_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 2','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,3)
pcolor(ens_03_step02_ana)
set(gca,'ytick',[])
set(gca,'xtick',[])
title('Ensemble state 3','fontsize',20)
set(gca,'clim',[-1 1])
subplot(2,2,4)
pcolor(ens_04_step02_ana)
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



% Compute rmse for an ensemble member
rmse_ens =0;
for i=1:12
    for j=1:12
    rmse_ens = rmse_ens + (ens_04_step02_ana(i,j) - true_step2(i,j))^2;
    end    
end
%rmse_ens = sqrt(rmse_ens/144)

var_for=zeros(12);
var_ana=zeros(12);
for i=1:12
    for j=1:12
            var_for(i,j) = var_for(i,j) ...
                + (ens_01_step02_for(i,j) - state_step02_for(i,j))^2 ...
                + (ens_02_step02_for(i,j) - state_step02_for(i,j))^2 ...
                + (ens_03_step02_for(i,j) - state_step02_for(i,j))^2 ...
                + (ens_04_step02_for(i,j) - state_step02_for(i,j))^2;
            var_ana(i,j) = var_ana(i,j) ...
                + (ens_01_step02_ana(i,j) - state_step02_ana(i,j))^2 ...
                + (ens_02_step02_ana(i,j) - state_step02_ana(i,j))^2 ...
                + (ens_03_step02_ana(i,j) - state_step02_ana(i,j))^2 ...
                + (ens_04_step02_ana(i,j) - state_step02_ana(i,j))^2;
    end
end
var_for = 0.2*var_for;
var_ana = 0.2*var_ana;

figure
subplot(1,2,1)
pcolor(var_for)
colormap(inferno)
set(gca,'clim',[0, 0.7],'fontsize',20)
title('estimated forecast variances')
set(gca,'ytick',[])
set(gca,'xtick',[])
colorbar
subplot(1,2,2)
pcolor(var_ana)
title('estimated analysis variances')
set(gca,'ytick',[])
set(gca,'xtick',[])
colormap(inferno)
set(gca,'clim',[0, 0.1],'fontsize',20)
colorbar
set(gcf,'position',[1000,100,946, 369])

ensf=zeros(144,4);
ensa=zeros(144,4);
cnt=1;
for j=1:12
    for i=1:12
    ensf(cnt,1) = ens_01_step02_for(i,j)-state_step02_for(i,j);
    ensf(cnt,2) = ens_02_step02_for(i,j)-state_step02_for(i,j);
    ensf(cnt,3) = ens_03_step02_for(i,j)-state_step02_for(i,j);
    ensf(cnt,4) = ens_04_step02_for(i,j)-state_step02_for(i,j);
    ensa(cnt,1) = ens_01_step02_ana(i,j)-state_step02_ana(i,j);
    ensa(cnt,2) = ens_02_step02_ana(i,j)-state_step02_ana(i,j);
    ensa(cnt,3) = ens_03_step02_ana(i,j)-state_step02_ana(i,j);
    ensa(cnt,4) = ens_04_step02_ana(i,j)-state_step02_ana(i,j);
    cnt = cnt + 1;
    end
end
    
cov_for = 0.2*ensf*ensf';
cov_ana = 0.2*ensa*ensa';
figure
subplot(1,2,1)
pcolor(flipud(cov_for))
set(gca,'ytick',[])
set(gca,'xtick',[])
colormap(cm_nogreenc)
shading flat
set(gca,'fontsize',20)
title('forecast covariance matrix')
colorbar
subplot(1,2,2)
pcolor(flipud(cov_ana))
title('analysis covariance matrix')
set(gca,'ytick',[])
set(gca,'xtick',[])
colormap(cm_nogreenc)
set(gca,'fontsize',20)
colorbar
shading flat
set(gcf,'position',[100,100,946, 369])

for i=1:144
    for j=1:144
        corf(i,j) = cov_for(i,j)/sqrt(cov_for(i,i))/sqrt(cov_for(j,j));
        cora(i,j) = cov_ana(i,j)/sqrt(cov_ana(i,i))/sqrt(cov_ana(j,j));
    end
end
figure
subplot(1,2,1)
pcolor(flipud(corf))
shading flat
title('forecast correlation matrix')
set(gca,'fontsize',20)
set(gca,'ytick',[])
set(gca,'xtick',[])
colorbar
subplot(1,2,2)
pcolor(flipud(cora))
shading flat
title('analysis correlation matrix')
set(gca,'fontsize',20)
set(gca,'ytick',[])
set(gca,'xtick',[])
colorbar
set(gcf,'position',[100,550,946, 369])

