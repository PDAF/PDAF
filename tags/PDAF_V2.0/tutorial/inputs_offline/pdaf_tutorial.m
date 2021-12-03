% Script to generate files for the PDAF offline tutorial

dim_x = 36;         % Grid dimension in x-direction
dim_y = 18;         % Grid dimension in y-direction
dim_ens = 9;        % Maximum ensemble size
stddev_obs = 0.5;   % Observation error standard deviation
dxobs = 5;          % x-Grid spacing for observations type A
dyobs = 4;          % y-Grid spacing for observations type A
dxobsB = 11;        % x-Grid spacing for observations type B
dyobsB = 8;         % y-Grid spacing for observations type B
obsB_offsetx = -4;  % x-offset in position of observations type B
obsB_offsety = -2;  % y-offset in position of observations type B

dowrite = 1;        % 1 to write files

% Locations of observations not placed at grid points (x, y) - type C
obs_interp = [3.0 2.1; ...
     3.4 6.8; ...
     6.1 6.8; ...
     8.9 7.6; ...
     8.9 14.9; ...
     20.0 6.4; ...
     20.4 16.1; ...
     14.1 10.2; ...
     31.0 5.2; ...
     31.2 11.9; ...
     28.9 14.9];

% Reset random number generation
rng('default')

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
set(gca,'fontsize',20)
cb=colorbar
set(cb,'fontsize',20)
title('True field','fontsize',24)


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
    set(gca,'fontsize',20)
    cb=colorbar
    set(cb,'fontsize',20)
    title(['Ensemble member ' num2str(k)],'fontsize',24)
end

figure
field_plot=zeros(dim_y+1, dim_x+1);
ensmean = mean(ens,3);
field_plot(1:dim_y,1:dim_x) = ensmean;
pcolor(field_plot)
set(gca,'fontsize',20)
cb=colorbar
set(cb,'fontsize',20)
title('Initial estimate (ensemble mean)','fontsize',24)

% Observations
obs_error = stddev_obs * randn(dim_y, dim_x);
full_obs = field + obs_error;

field_plot=zeros(dim_y+1, dim_x+1);
field_plot(1:dim_y,1:dim_x) = full_obs;
figure
pcolor(field_plot)
set(gca,'fontsize',20)
cb=colorbar
set(cb,'fontsize',20)
title('Perturbed true state','fontsize',24)

obs = zeros(dim_y, dim_x)-999;
for j=dxobs:dxobs:dim_x
    for i=dyobs:dyobs:dim_y
        obs(i,j) = full_obs(i,j);
    end
end

obsB = zeros(dim_y, dim_x)-999;
for j=dxobsB+obsB_offsetx:dxobsB:dim_x
    for i=dyobsB+obsB_offsety:dyobsB:dim_y
        obsB(i,j) = full_obs(i,j);
    end
end

field_plot=zeros(dim_y+1, dim_x+1);
field_plot(1:dim_y,1:dim_x) = obs;
figure
pcolor(field_plot)
set(gca,'fontsize',20)
cb=colorbar
set(cb,'fontsize',20)
title(['Type A: 28 Observations'],'fontsize',24)
set(gca,'clim',[-3 3])

field_plot=zeros(dim_y+1, dim_x+1);
field_plot(1:dim_y,1:dim_x) = obsB;
figure
pcolor(field_plot)
set(gca,'fontsize',20)
cb=colorbar
set(cb,'fontsize',20)
title(['Type B: 6 Observations'],'fontsize',24)
set(gca,'clim',[-3 3])

% Interpolated observations
iobs_error = stddev_obs * randn(length(obs_interp));
for i=1:length(obs_interp)
    % Get closted grid points
    gx(1) = floor(obs_interp(i,1));
    gx(2) = ceil(obs_interp(i,1));
    if gx(2)==gx(1)
        gx(2) = gx(2)+1;
    end
    gy(1) = floor(obs_interp(i,2));
    gy(2) = ceil(obs_interp(i,2));
    if gy(2)==gy(1)
        gy(2) = gy(2)+1;
    end

    % Compute interpolation coefficients
    denum = (gx(2)-gx(1))*(gy(2)-gy(1));
    icoeff(1) = (gx(2)-obs_interp(i,1))*(gy(2)-obs_interp(i,2))/denum;
    icoeff(2) = (obs_interp(i,1)-gx(1))*(gy(2)-obs_interp(i,2))/denum;
    icoeff(3) = (gx(2)-obs_interp(i,1))*(obs_interp(i,2)-gy(1))/denum;
    icoeff(4) = (obs_interp(i,1)-gx(1))*(obs_interp(i,2)-gy(1))/denum;

    % Interpolate
    iobs(i,1) = icoeff(1)*field(gy(1),gx(1))+icoeff(2)*field(gy(1),gx(2))+ ...
        icoeff(3)*field(gy(2),gx(1))+icoeff(4)*field(gy(2),gx(2));

    % Add error
    iobs(i,1) = iobs(i,1) + iobs_error(i);

    % Augment with coordinates
    iobs(i,2:3) = obs_interp(i,1:2);
end 

field_plot=zeros(dim_y+1, dim_x+1) - 999;
for i=1:length(obs_interp)
    field_plot(floor(obs_interp(i,2)),floor(obs_interp(i,1))) = iobs(i,1);
end
figure
pcolor(field_plot)
set(gca,'fontsize',20)
cb=colorbar;
set(cb,'fontsize',20)
title([num2str(length(obs_interp)) ' Obs. with bi-linear interpolation'],'fontsize',24)
set(gca,'clim',[-3 3])


%%%%%%%%%%%%%%%%%%% Write files

if dowrite == 1

    % True field
    fid = fopen('true.txt','w');
    for i=1:dim_y
        fprintf(fid,'%14.8f',field(i,:));
        fprintf(fid,'\n')
    end
    fclose(fid)


    % Observations - type A
    fid = fopen('obs.txt','w');
    for i=1:dim_y
        fprintf(fid,'%14.6f',obs(i,:));
        fprintf(fid,'\n')
    end
    fclose(fid)

    % Observations - type B
    fid = fopen('obsB.txt','w');
    for i=1:dim_y
        fprintf(fid,'%14.6f',obsB(i,:));
        fprintf(fid,'\n')
    end
    fclose(fid)

    % Interpolated observations
    fid = fopen(['obsC.txt'],'w');
    fprintf(fid,'%5i',length(obs_interp));
    fprintf(fid,'\n');
    for i=1:length(obs_interp)
        fprintf(fid,'%14.6f',iobs(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);

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
end