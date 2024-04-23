%spots = readtable('spots.csv');

addpath('C:\Users\woods\EnvDataExp\PartnerInfo\Final-Project\')





file = 'spots.csv'
data = readtable(file)




site = table2array(data(:,1));
ph = table2array(data(:,62));
lat = table2array(data(:,8));
lon = table2array(data(:,9));
pressure = table2array(data(:,10))

surfacepressureind = find(pressure<= 10)


ind = find(ph ~= -999)

pH_insitu = CO2SYS(data.PH_TOT(1:86706,:), data.TCARBN(1:86706,:), 3, 2, data.CTDSAL(1:86706,:), data.PH_TMP(1:86706,:), data.CTDTMP(1:86706,:), 0, data.CTDPRS(1:86706,:), 0, 0, 1, 10, 1);
%8.744402, -172.731079

%30.119914, -142.660725

latitudes = find(lat>=20 & lat <= 25);
longitudes = find(lon>= -160 & lon <= -150);

%inds = pH_insitu(latitudes,:)

surfacepressure = []


for i = 1:length(surfacepressureind)
    if surfacepressureind(i,1) <= 86706;
        surfacepressure(i) = surfacepressureind(i);
       

    end
end

inds = pH_insitu(surfacepressure',:);

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
