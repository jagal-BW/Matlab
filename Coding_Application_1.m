%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%            PART 1            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Team Members:  Jason Galvan + Emily Hatton 
%%%%%% Contribution:  Worked on all problems together

clear;  close all;

% Specify location of NetCDF data source

opendap = {'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.ACCESS1-0/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.BNU-ESM/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.CanESM2/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.CCSM4/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.CESM1-CAM5/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.CMCC-CESM/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.CNRM-CM5/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.CSIRO-Mk3-6-0/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.EC-EARTH/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.FGOALS-g2/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.FIO-ESM/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.GFDL-CM3/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.GISS-E2-R/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.HadGEM2-AO/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.inmcm4/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.IPSL-CM5A-LR/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.MIROC5/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.MPI-ESM-LR/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.MRI-CGCM3/.r1i1p1/.pr/dods',...
    'http://strega.ldeo.columbia.edu:81/CMIP5/.monthly/.byScenario/.rcp85/.atmos/.mon/.pr/.NorESM1-M/.r1i1p1/.pr/dods'};


for i = 1:20
    
    % latitudes in degrees N
    lats = ncread(opendap{i},'lat');

    % longitudes in degrees E
    lons = ncread(opendap{i},'lon');

    % Read data
    pr = ncread(opendap{i},'pr',[1 1 1],[length(lons) length(lats) 1140]);
    
    % Temporal mean (annual mean)
    for j = 1:95
        
       tempmean(:,:,j) = mean(pr(:,:,(j-1)*12+1:j*12),3); 
        
    end
    
    % Add Longitude mean to Temporal mean
    longmean = mean(tempmean,1);
    longmean = permute(longmean,[2 3 1]);
    
    % Add Latitude mean to Longitude and Temporal mean
    for j = 1:95
        
        annualmean(j) = sum(longmean(:,j).*cosd(lats))/sum(cosd(lats));
    
    end
    
    % Transpose and convert to mm/day
    modelmean(i,:) = 86400*annualmean';
    
    clearvars lats lons pr tempmean longmean annualmean j
    
end

multimodelmean = mean(modelmean);

% Plotting

% Configuring plot style
set(gca,'linestyleorder',{'-',':','-.'},...
'colororder',get(gca,'colororder'),'nextplot','add')

for i = 1:20
   
   plot(2005+[1:95],modelmean(i,:))
   
   if i == 1
       hold on
   end
    
end

% Add multimodel to plot
plot(2005+[1:95],multimodelmean)

hold off

legend('ACCESS1-0','BNU-ESM','CanESM2','CCSM4','CESM1-CAM5','CMCC-CESM',...
    'CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-g2','FIO-ESM','GFDL-CM3',...
    'GISS-E2-R','HadGEM2-AO','inmcm4','IPSL-CM5A-LR','MIROC5','MPI-ESM-LR',...
    'MRI-CGCM3','NorESM1-M','Multimodel Mean')
title('Annual mean precipitation for different models')
xlabel('Year')
ylabel('Precipitation [mm/day]')
xlim([2006 2100])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%            PART 2            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% geo map of multimodel mean precip change between 2006-2025 and 2081-2100

for i = 1:20
    
    % latitudes in degrees N
    lats = ncread(opendap{i},'lat');

    % longitudes in degrees E
    lons = ncread(opendap{i},'lon');

    % Read data
    pr = ncread(opendap{i},'pr',[1 1 1],[length(lons) length(lats) 1140]);
    
    % Mean change for ith model
    meanchange = mean(pr(:,:,901:1140),3) - mean(pr(:,:,1:240),3);
    
    % Regrid
    [~,~,meanchanger(:,:,i)] = regrid_cmip5(lons,lats,meanchange);
    
    clearvars lats lons pr meanchange
    
end

multimodelmeanchange = 86400*mean(meanchanger,3);

% Get matrix of lat/lon values
latm=double(repmat([-89.5:89.5]',1,360));
lonm=double(repmat([0.5:359.5],180,1));

% Change longitude range
lonm = lonm + 180;

figure

% Define geo axes
gx = axesm('MapProjection','eckert4');

% Display temperature data
geoshow(gx, latm, lonm, multimodelmeanchange', 'DisplayType', 'surface');

% Add frame
framem;

% Add colorbar
cb = colorbar('southoutside');

% Add label to colorbar and adjust font sizes
cb.Label.String = 'Precipitation change [mm/day]';
cb.Label.FontSize = 14;
cb.FontSize = 12;

% Add title
title('Multimodel mean change in precipitations (2081-2100;2006-2025)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%            PART 3            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% map where at least 90% of models agree on precipitation 

countsign = zeros();
for i = 1:20
    
   countsign = countsign + sign(meanchanger(:,:,i));
    
end

ind = (countsign > 16) + (countsign < -16);

s = 0;
for i = 1:360
    for j = 1:180
        if ind(i,j) == 1
            s = s + 1;
            lonvect(s) = 179.5 + i;
            latvect(s) = j - 90.5;
        end
    end
end


% Get matrix of lat/lon values
latm=double(repmat([-89.5:89.5]',1,360));
lonm=double(repmat([0.5:359.5],180,1));

% Change longitude range
lonm = lonm + 180;

figure

% Define geo axes
gx = axesm('MapProjection','eckert4');

% Display temperature data
geoshow(gx, latm, lonm, multimodelmeanchange', 'DisplayType', 'surface');

% Add frame
framem;

% Add colorbar
cb = colorbar('southoutside');

% Add label to colorbar and adjust font sizes
cb.Label.String = 'Precipitation change [mm/day]';
cb.Label.FontSize = 14;
cb.FontSize = 12;

% Add title
title('Multimodel mean change in precipitations (2081-2100;2006-2025)');
hold on

geoshow(latvect, lonvect, 'DisplayType', 'point', 'Marker', '+');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%            PART 4            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Creating NetCdf file with necessary data to reproduce plots

nccreate('part4.nc','lat','Dimensions',{'lat',180});
nccreate('part4.nc','lon','Dimensions',{'lon',360});
nccreate('part4.nc','multimodelmeanchange','Dimensions',{'lon',360,'lat',180});
nccreate('part4.nc','ind','Dimensions',{'lon',360,'lat',180});

ncwrite('part4.nc','lat',-89.5:89.5);
ncwrite('part4.nc','lon',0.5:359.5);
ncwrite('part4.nc','multimodelmeanchange',multimodelmeanchange);
ncwrite('part4.nc','ind',ind);

ncwriteatt('part4.nc','lat','units','degrees N')
ncwriteatt('part4.nc','lon','units','degrees E')
ncwriteatt('part4.nc','multimodelmeanchange','units','mm/day')
ncwriteatt('part4.nc','ind','units','N/A')
