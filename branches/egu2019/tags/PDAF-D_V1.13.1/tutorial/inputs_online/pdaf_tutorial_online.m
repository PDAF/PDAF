% Script to generate files for the PDAF online tutorial

dim_x = 36
dim_y = 18
dim_ens = 9
dim_step = 18
stddev_obs = 0.5
dxobs = 5;
dyobs = 4;


% True field
for j=1:dim_x
    for i=1:dim_y
        field(i,j,1) = sin(2*pi*(i/18+j/36));
    end
end

for step=2:dim_step+1
    for j=1:dim_x
        for i=1:dim_y-1
            field(i+1,j,step) = field(i,j,step-1);
        end
    end
    field(1,:,step) = field(dim_y,:,step-1);
end

for step=1:dim_step+1
    field_plot=zeros(dim_y+1, dim_x+1);
    field_plot(1:dim_y,1:dim_x) = field(:,:,step);
    figure
    pcolor(field_plot)
    set(gca,'fontsize',16)
    cb=colorbar
    set(cb,'fontsize',16)
    title(['True field, step ' num2str(step)],'fontsize',18)
end


% Ensemble states
for k=1:dim_ens
    for j=1:dim_x
        for i=1:dim_y
            ens(i,j,k) = sin(2*pi*(i/18+j/36)+2*0.5*pi*(k+5)/dim_ens);
        end
    end
    figure
    field_plot=zeros(dim_y+1, dim_x+1);
    field_plot(1:dim_y,1:dim_x) = ens(:,:,k);
    pcolor(field_plot)
    set(gca,'fontsize',16)
    cb=colorbar
    set(cb,'fontsize',16)
    title(['Ensemble member ' num2str(k)],'fontsize',18)
end

figure
field_plot=zeros(dim_y+1, dim_x+1);
ensmean = mean(ens,3);
field_plot(1:dim_y,1:dim_x) = ensmean;
pcolor(field_plot)
set(gca,'fontsize',16)
cb=colorbar
set(cb,'fontsize',16)
title('Initial estimate (ensemble mean)','fontsize',18)

% Observations
obs_error = stddev_obs * randn(dim_y, dim_x, dim_step+1);
full_obs = field + obs_error;

for step=1:dim_step+1
    field_plot=zeros(dim_y+1, dim_x+1);
    field_plot(1:dim_y,1:dim_x) = full_obs(:,:,step);
    figure
    pcolor(field_plot)
    set(gca,'fontsize',16)
    cb=colorbar
    set(cb,'fontsize',16)
    title(['Disturbed true state, step ' num2str(step-1)],'fontsize',18)
end

obs = zeros(dim_y, dim_x, dim_step+1)-999;
for step=1:dim_step+1
    for j=dxobs:dxobs:dim_x
        for i=dyobs:dyobs:dim_y
            obs(i,j,step) = full_obs(i,j,step);
        end
    end
end

for step=1:dim_step+1
    field_plot=zeros(dim_y+1, dim_x+1);
    field_plot(1:dim_y,1:dim_x) = obs(:,:,step);
    figure
    pcolor(field_plot)
    set(gca,'fontsize',16)
    cb=colorbar
    set(cb,'fontsize',16)
    title(['28 Observations used for analysis, step ' num2str(step-1)],'fontsize',18)
    set(gca,'clim',[-3 3])
end


% Write files

% True field
fid = fopen(['true_initial.txt'],'w');
for i=1:dim_y
    fprintf(fid,'%14.8f',field(i,:,1));
    fprintf(fid,'\n')
end
fclose(fid)
for step=2:dim_step+1
    fid = fopen(['true_step' num2str(step-1) '.txt'],'w');
    for i=1:dim_y
        fprintf(fid,'%14.8f',field(i,:,step));
        fprintf(fid,'\n')
    end
    fclose(fid)
end

% Observations
for step=2:dim_step+1
    fid = fopen(['obs_step' num2str(step-1) '.txt'],'w');
    for i=1:dim_y
        fprintf(fid,'%14.6f',obs(i,:,step));
        fprintf(fid,'\n')
    end
    fclose(fid)
end

% Ensemble
for k=1:dim_ens
    fid = fopen(['ens_' num2str(k) '.txt'],'w');
    for i=1:dim_y
        fprintf(fid,'%14.8f',ens(i,:,k));
        fprintf(fid,'\n')
    end
    fclose(fid)
end


% Ensemble mean
fid = fopen('state_ini.txt','w');
for i=1:dim_y
    fprintf(fid,'%14.6f',ensmean(i,:));
    fprintf(fid,'\n')
end
fclose(fid)
