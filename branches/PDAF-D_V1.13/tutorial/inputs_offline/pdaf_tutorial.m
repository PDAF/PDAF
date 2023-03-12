% Script to generate files for the PDAF offline tutorial

dim_x = 36
dim_y = 18
dim_ens = 9
stddev_obs = 0.5
dxobs = 5;
dyobs = 4;


% True field
for j=1:dim_x
    for i=1:dim_y
        field(i,j) = sin(pi*(i/18+j/36));
    end
end


field_plot=zeros(dim_y+1, dim_x+1);
field_plot(1:dim_y,1:dim_x) = field;
figure
pcolor(field_plot)
set(gca,'fontsize',16)
cb=colorbar
set(cb,'fontsize',16)
title('True field','fontsize',18)


% Ensemble states
for k=1:dim_ens
    for j=1:dim_x
        for i=1:dim_y
            ens(i,j,k) = sin(pi*(i/18+j/36)+0.5*pi*(k+5)/dim_ens);
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
obs_error = stddev_obs * randn(dim_y, dim_x);
full_obs = field + obs_error;

field_plot=zeros(dim_y+1, dim_x+1);
field_plot(1:dim_y,1:dim_x) = full_obs;
figure
pcolor(field_plot)
set(gca,'fontsize',16)
cb=colorbar
set(cb,'fontsize',16)
title('Disturbed true state','fontsize',18)

obs = zeros(dim_y, dim_x)-999;
for j=dxobs:dxobs:dim_x
    for i=dyobs:dyobs:dim_y
        obs(i,j) = full_obs(i,j);
    end
end

field_plot=zeros(dim_y+1, dim_x+1);
field_plot(1:dim_y,1:dim_x) = obs;
figure
pcolor(field_plot)
set(gca,'fontsize',16)
cb=colorbar
set(cb,'fontsize',16)
title(['28 Observations used for analysis'],'fontsize',18)
set(gca,'clim',[-3 3])


% Write files

% True field
fid = fopen('true.txt','w');
for i=1:dim_y
    fprintf(fid,'%14.8f',field(i,:));
    fprintf(fid,'\n')
end
fclose(fid)


% Observations
fid = fopen('obs.txt','w');
for i=1:dim_y
    fprintf(fid,'%14.6f',obs(i,:));
    fprintf(fid,'\n')
end
fclose(fid)

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
