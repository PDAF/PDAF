% Script to generate files for the PDAF online tutorial

dim_x = 36;         % Grid dimension in x-direction
dim_y = 18;         % Grid dimension in y-direction
dim_ens = 9;        % Maximum ensemble size
dim_step = 18;      % Number of time steps
stddev_obs = 0.5;   % Observation error standard deviation
dxobs = 5;          % x-Grid spacing for observations type A
dyobs = 4;          % y-Grid spacing for observations type A
dxobsB = 11;        % x-Grid spacing for observations type B
dyobsB = 8;         % y-Grid spacing for observations type B
obsB_offsetx = -4;  % x-offset in position of observations type B
obsB_offsety = -2;  % y-offset in position of observations type B

dowrite = 1;        % 1 to write files

% Locations of observations not placed at grid points (x, y)
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
    set(gca,'fontsize',20)
    cb=colorbar;
    set(cb,'fontsize',20)
    title(['True field, step ' num2str(step)],'fontsize',24)
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
    set(gca,'fontsize',20)
    cb=colorbar;
    set(cb,'fontsize',20)
    title(['Ensemble member ' num2str(k)],'fontsize',24)
end

figure
field_plot=zeros(dim_y+1, dim_x+1);
ensmean = mean(ens,3);
field_plot(1:dim_y,1:dim_x) = ensmean;
pcolor(field_plot)
set(gca,'fontsize',20)
cb=colorbar;
set(cb,'fontsize',20)
title('Initial estimate (ensemble mean)','fontsize',24)


% Ensemble states - inverted
for k=1:dim_ens
    for j=1:dim_x
        for i=1:dim_y
            ensB(i,j,k) = sin(2*pi*(i/18-j/36)+2*0.5*pi*(k+5)/dim_ens);
        end
    end
    figure
    field_plot=zeros(dim_y+1, dim_x+1);
    field_plot(1:dim_y,1:dim_x) = ensB(:,:,k);
    pcolor(field_plot)
    set(gca,'fontsize',20)
    cb=colorbar;
    set(cb,'fontsize',20)
    title(['B-Ensemble member ' num2str(k)],'fontsize',24)
end

figure
field_plot=zeros(dim_y+1, dim_x+1);
ensmeanB = mean(ensB,3);
field_plot(1:dim_y,1:dim_x) = ensmeanB;
pcolor(field_plot)
set(gca,'fontsize',20)
cb=colorbar;
set(cb,'fontsize',20)
title('B: Initial estimate (ensemble mean)','fontsize',24)

% Observations
obs_error = stddev_obs * randn(dim_y, dim_x, dim_step+1);
full_obs = field + obs_error;

for step=1:dim_step+1
    field_plot=zeros(dim_y+1, dim_x+1);
    field_plot(1:dim_y,1:dim_x) = full_obs(:,:,step);
    figure
    pcolor(field_plot)
    set(gca,'fontsize',20)
    cb=colorbar;
    set(cb,'fontsize',20)
    title(['Perturbed true state, step ' num2str(step-1)],'fontsize',24)
end

obs = zeros(dim_y, dim_x, dim_step+1)-999;
for step=1:dim_step+1
    for j=dxobs:dxobs:dim_x
        for i=dyobs:dyobs:dim_y
            obs(i,j,step) = full_obs(i,j,step);
        end
    end
end

obsB = zeros(dim_y, dim_x, dim_step+1)-999;
for step=1:dim_step+1
    for j=dxobsB+obsB_offsetx:dxobsB:dim_x
        for i=dyobsB+obsB_offsety:dyobsB:dim_y
            obsB(i,j,step) = full_obs(i,j,step);
        end
    end
end

for step=1:dim_step+1
    field_plot=zeros(dim_y+1, dim_x+1);
    field_plot(1:dim_y,1:dim_x) = obs(:,:,step);
    figure
    pcolor(field_plot)
    set(gca,'fontsize',20)
    cb=colorbar;
    set(cb,'fontsize',20)
    title(['Type A: 28 Observations, step ' num2str(step-1)],'fontsize',24)
    set(gca,'clim',[-3 3])
end

for step=1:dim_step+1
    field_plot=zeros(dim_y+1, dim_x+1);
    field_plot(1:dim_y,1:dim_x) = obsB(:,:,step);
    figure
    pcolor(field_plot)
    set(gca,'fontsize',20)
    cb=colorbar;
    set(cb,'fontsize',20)
    title(['Type B: 6 Observations, step ' num2str(step-1)],'fontsize',24)
    set(gca,'clim',[-3 3])
end

% Interpolated observations
iobs_error = stddev_obs * randn(length(obs_interp), dim_step+1);
for step=1:dim_step+1
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
        iobs(i,1,step) = icoeff(1)*field(gy(1),gx(1),step)+icoeff(2)*field(gy(1),gx(2),step)+ ...
            icoeff(3)*field(gy(2),gx(1),step)+icoeff(4)*field(gy(2),gx(2),step);

        % Add error
        iobs(i,1,step) = iobs(i,1,step) + iobs_error(i,step);

        % Augment with coordinates
        iobs(i,2:3,step) = obs_interp(i,1:2);
    end 
    
    field_plot=zeros(dim_y+1, dim_x+1) - 999;
    for i=1:length(obs_interp)
        field_plot(floor(obs_interp(i,2)),floor(obs_interp(i,1))) = iobs(i,1,step);
    end
    figure
    pcolor(field_plot)
    set(gca,'fontsize',20)
    cb=colorbar;
    set(cb,'fontsize',20)
    title([num2str(length(obs_interp)) ' Obs. with bi-linear interpolation, step ' num2str(step-1)],'fontsize',24)
    set(gca,'clim',[-3 3])
end
    

%%%%%%%%%%%%%%%%%%% Write files

if dowrite == 1
    
    % True field
    fid = fopen(['true_initial.txt'],'w');
    for i=1:dim_y
        fprintf(fid,'%14.8f',field(i,:,1));
        fprintf(fid,'\n');
    end
    fclose(fid);
    for step=2:dim_step+1
        fid = fopen(['true_step' num2str(step-1) '.txt'],'w');
        for i=1:dim_y
            fprintf(fid,'%14.8f',field(i,:,step));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end

    % Observations - Type A
    for step=2:dim_step+1
        fid = fopen(['obs_step' num2str(step-1) '.txt'],'w');
        for i=1:dim_y
            fprintf(fid,'%14.6f',obs(i,:,step));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end

    % Observations - Type B
    for step=2:dim_step+1
        fid = fopen(['obsB_step' num2str(step-1) '.txt'],'w');
        for i=1:dim_y
            fprintf(fid,'%14.6f',obsB(i,:,step));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end


    % Interpolated observations
    for step=2:dim_step+1
        fid = fopen(['iobs_step' num2str(step-1) '.txt'],'w');
        fprintf(fid,'%5i',length(obs_interp));
        fprintf(fid,'\n');
        for i=1:length(obs_interp)
            fprintf(fid,'%14.6f',iobs(i,:,step));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end

    % Ensemble
    for k=1:dim_ens
        fid = fopen(['ens_' num2str(k) '.txt'],'w');
        for i=1:dim_y
            fprintf(fid,'%14.8f',ens(i,:,k));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end

    % Ensemble
    for k=1:dim_ens
        fid = fopen(['ensB_' num2str(k) '.txt'],'w');
        for i=1:dim_y
            fprintf(fid,'%14.8f',ensB(i,:,k));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end


    % Ensemble mean
    fid = fopen('state_ini.txt','w');
    for i=1:dim_y
        fprintf(fid,'%14.6f',ensmean(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end