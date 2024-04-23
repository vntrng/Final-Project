
load('C:\Users\woods\EnvDataExp\PartnerInfo\Final Project\GLODAPv2.2023_Merged_Master_File.mat')

lat = G2latitude(:,:);
lon = G2longitude(:,:);

year = G2year(:,:);
hour = G2hour(:,:);
minutes = G2minute(:,:);
month = G2month(:,:);
day = G2day(:,:);

pH = G2phts25p0f(:,:);

temperature = G2temperature(:,:);

%variables = table('temperature', 'year', 'month', 'minutues')


%Attempting to make the time in a format that matlab can read
dates = {};

for i = 1:length(year)
    
    dates{i} = datetime(year(i,:),month(i,:),day(i,:),hour(i,:),minutes(i,:),0,0);
    
end

%%
addpath("C:\Users\woods\Downloads")

filename1 = 'hawaii_d90f_20ee_c4cb_e8d8_79ca_4468.nc';
%1a. Use the function "ncdisp" to display information about the data contained in this file
%-->
ncdisp(filename1)

%% Global_Coral_Bleaching Database
%addpath('/Users/amyz/Documents/Data Exploration/Final-Project/')
filename = 'Global_Coral_Bleaching_Database.csv';
stationdata = readtable(filename);

lat = table2array(stationdata(:,5));
lon = table2array(stationdata(:,6));
comments = table2array(stationdata(:,16));

%row_has_NA = any(strcmp(comments(:,1), 'N/A'));
rows = find(contains(comments,'N/A'));

goodIndices = setdiff(1:35053, rows);
cleanedcomments = comments(goodIndices, :);
cleanedlat = lat(goodIndices,:);
cleanedlon = lon(goodIndices,:);
figure(2); clf
worldmap world
plotm(cleanedlat,cleanedlon,'m.','MarkerSize',10);
geoshow('landareas.shp','FaceColor','white')
title('Location for Coral Reef Bleaching')

%%
filename3 = 'NOAA_DHW_monthly_1a19_b145_a3b0_U1713885190139.nc'

target_lat = 42.171949, % Replace with your desired latitude
target_lon = -68.261707; % Replace with your desired longitude


ncdisp(filename3)
lat2 = ncread(filename3, 'latitude');
lon2 = ncread(filename3, 'longitude');
lat2 = double(lat2);
lon2 = double(lon2);
mask2 = ncread(filename3, 'mask'); % Read the land-sea mask
mask_slice = mask2(:,:,1); 

time2 = ncread(filename3, 'time');
time3 = datetime(time2, 'ConvertFrom', 'posixtime');

sst2 = ncread(filename3, 'sea_surface_temperature');
ssta2 = ncread(filename3, 'sea_surface_temperature_anomaly');

% Calculate mean SST anomaly across time
mean_ssta = mean(ssta2, 3);

% Create a plot for global mean SST anomaly
figure;
worldmap('World'); % Set map to world
load coastlines; % Load coastlines data for plotting
plotm(coastlat, coastlon, 'k'); % Plot coastlines
pcolorm(lat2, lon2, mean_ssta'); % Plot mean SST anomaly
cb = colorbar; % Create colorbar
cb.Label.String = 'Temperature (°C)'; % Add units to the colorbar label
title('Global Mean Sea Surface Temperature Anomaly 1985-Present');
hold on
plotm(cleanedlat,cleanedlon,'m.','MarkerSize',10);
hold off

% Calculate mean SST across time
mean_sst = mean(sst2, 3);

% Create a plot for global mean SST
figure;
worldmap('World'); % Set map to world
plotm(coastlat, coastlon, 'k'); % Plot coastlines
pcolorm(lat2, lon2, mean_sst'); % Plot mean SST
cb = colorbar; % Create colorbar
cb.Label.String = 'Temperature (°C)'; % Add units to the colorbar label
title('Average Sea Surface Temperature 1985-Present');
hold on
plotm(cleanedlat,cleanedlon,'m.','MarkerSize',10);
hold off


[lon_size, lat_size, ~] = size(sst2);
rate_of_change = zeros(lon_size, lat_size);

% Fit a linear model to get the rate of change at each grid point
for i = 1:lon_size
    for j = 1:lat_size
        if mask_slice(i, j) == 0 % Check if the point is water
            % Use polyfit to get slope (rate of change)
            p = polyfit(time2, squeeze(sst2(i, j, :)), 1);
            rate_of_change(i, j) = p(1); % Slope
        end
    end
end

% Replace zero values with NaN
rate_of_change(rate_of_change == 0) = NaN;

% Create a plot for global rate of change of SST
figure;
worldmap('World'); % Set map to world
plotm(coastlat, coastlon, 'k'); % Plot coastlines
pcolorm(lat2, lon2, rate_of_change'); % Plot rate of change
colorbar; % Add colorbar
title('Rate of Change of Sea Surface Temperature');
hold on
plotm(cleanedlat,cleanedlon,'m.','MarkerSize',10);
hold off

% Find the nearest grid point
[~,lat_idx] = min(abs(lat2 - target_lat));
[~,lon_idx] = min(abs(lon2 - target_lon));

sst_target_data = sst2(lon_idx, lat_idx,:)
sst_target_data = squeeze(sst_target_data)
ssta_target_data = ssta2(lon_idx, lat_idx,:)
ssta_target_data = squeeze(ssta_target_data)

figure;
plot(time3, sst_target_data, 'b-', 'DisplayName', 'SST'); % Blue line for SST

x = datenum(time3)

p = polyfit(x,sst_target_data,1)

trendlinesst = polyval(p,x)

hold on;
plot(time3, trendlinesst, '-', 'DisplayName', 'SST Trendline'); % Plot the trendline
slope_SST = p(1) % Slope of the line

hold on;
plot(time3, ssta_target_data, 'r--', 'DisplayName', 'SSTA'); % Red dashed line for SSTA

p2 = polyfit(x,ssta_target_data,1);

trendlinessta = polyval(p2,x);

hold on;
plot(time3, trendlinessta, '-', 'DisplayName', 'SST Anomaly Trendline'); % Plot the trendline
slope_SSTA = p(1) % Slope of the line

xlabel('Time');
ylabel('Temperature (°C)');
title('Time Series of Sea Surface Temperature and Anomaly 1985 - Present');
legend() % Add a legend to identify the lines'
hold off;



%% Global_Coral_Bleaching Database
%addpath('/Users/amyz/Documents/Data Exploration/Final-Project/')
addpath("C:\Users\woods\EnvDataExp\PartnerInfo\Final-Project\")
filename = 'Global_Coral_Bleaching_Database.csv';
stationdata = readtable(filename);

lat = table2array(stationdata(:,5));
lon = table2array(stationdata(:,6));



%% Global_Coral_Bleaching Database
%addpath('/Users/amyz/Documents/Data Exploration/Final-Project/')
filename = 'Global_Coral_Bleaching_Database.csv';
stationdata = readtable(filename);

lat = table2array(stationdata(:,5));
lon = table2array(stationdata(:,6));
comments = table2array(stationdata(:,16));

%row_has_NA = any(strcmp(comments(:,1), 'N/A'));
rows = find(contains(comments,'N/A'));

goodIndices = setdiff(1:35053, rows);
cleanedcomments = comments(goodIndices, :);
cleanedlat = lat(goodIndices,:);
cleanedlon = lon(goodIndices,:);
figure(2); clf
worldmap world
plotm(cleanedlat,cleanedlon,'m.','MarkerSize',10);
geoshow('landareas.shp','FaceColor','white')
title('Location for Coral Reef Bleaching')
