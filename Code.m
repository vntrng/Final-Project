%spots = readtable('spots.csv');

addpath('C:\Users\woods\EnvDataExp\PartnerInfo\Final-Project\')

<<<<<<< HEAD




file = 'spots.csv'
data = readtable(file);


=======
file = 'spots.csv';
data = readtable(file);
>>>>>>> 647c469ed6e9cf9da6f31ca4c736731562810915


site = table2array(data(:,1));
ph = table2array(data(:,62));
lat = table2array(data(:,8));
lon = table2array(data(:,9));
pressure = table2array(data(:,10));
date = table2array(data(:,6));
surfacepressureind = find(pressure<= 10);


%pH_insitu = CO2SYS(data.PH_TOT(1:86706,:), data.TCARBN(1:86706,:), 3, 2, data.CTDSAL(1:86706,:), data.PH_TMP(1:86706,:), data.CTDTMP(1:86706,:), 0, data.CTDPRS(1:86706,:), 0, 0, 1, 10, 1);
%8.744402, -172.731079

%30.119914, -142.660725

latitudes = find(lat>=20 & lat <= 25);
longitudes = find(lon>= -160 & lon <= -150);

%latitudes = find(lat>=-50 & lat <= -40);
%longitudes = find(lon>= 165 & lon <= 175);
%inds = pH_insitu(latitudes,:)

surfacepressure = [];


for i = 1:length(surfacepressureind);
    %if surfacepressureind(i,1) <= 104915 ;
    if surfacepressureind(i,1) >= 86707;
        %if surfacepressureind(i,1) >= 104595;
        surfacepressure(i) = surfacepressureind(i);
<<<<<<< HEAD
        %end

=======
       
>>>>>>> 647c469ed6e9cf9da6f31ca4c736731562810915
    end
end
inds = surfacepressure';

%inds = intersect(latitudes(:,1),surfacepressure(:,1))

%inds = pH_insitu(surfacepressure',:);

ind = find(ph ~= -999);

inds2 = intersect(ind(:,1),inds(:,1));

%pH_measures = ph(inds2,1)

pH_insitu = CO2SYS(data.PH_TOT(inds2,1), data.TCARBN(inds2,1), 3, 2, data.CTDSAL(inds2,1), data.PH_TMP(inds2,1), data.CTDTMP(inds2,1), 0, data.CTDPRS(inds2,1), 0, 0, 1, 10, 1);

ph_output = pH_insitu(:,38);
x = date(inds2,1);
p = polyfit(x,ph_output,1);
<<<<<<< HEAD
slope = p(1) % Slope of the line
=======
slope = p(1); % Slope of the line
>>>>>>> 647c469ed6e9cf9da6f31ca4c736731562810915

trendline = polyval(p,x);

figure;

x = arrayfun(@num2str, x, 'UniformOutput', false);
x = datetime(x, 'InputFormat', 'yyyyMMdd');

plot(x, ph_output,'r', 'DisplayName', 'Hawaii pH');
%plot(x, ph_output,'r', 'DisplayName', 'Japan pH');
hold on
plot(x,trendline,'r--','DisplayName','Hawaii Trendline')
%plot(x,trendline,'r--','DisplayName','Japan Trendline')

hold on

%%
latitudes2 = find(lat>=0 & lat <= 20);
longitudes2 = find(lon>= -70 & lon <= -60);

%inds = pH_insitu(latitudes,:)

surfacepressure2 = [];


for i = 1:length(surfacepressureind)
    if surfacepressureind(i,1) <= 90982 
        if surfacepressureind(i,1) >= 86707;
        surfacepressure(i) = surfacepressureind(i);
        end
       

    end
end
inds3 = surfacepressure';
inds4 = intersect(ind(:,1),inds3(:,1));
pH_insitu = CO2SYS(data.PH_TOT(inds4,1), data.TCARBN(inds4,1), 3, 2, data.CTDSAL(inds4,1), data.PH_TMP(inds4,1), data.CTDTMP(inds4,1), 0, data.CTDPRS(inds4,1), 0, 0, 1, 10, 1);
ph_output1 = pH_insitu(:,38);
ph_output1 = real(ph_output1);
x2 = date(inds4,1);
p2 = polyfit(x2,ph_output1,1);
slope = p2(1) % Slope of the line

trendline2 = polyval(p2,x2);

x2 = arrayfun(@num2str, x2, 'UniformOutput', false);
x2 = datetime(x2, 'InputFormat', 'yyyyMMdd');

plot(x2, ph_output1,'b', 'DisplayName', 'Carribean pH');
hold on
plot(x2,trendline2,'b--','DisplayName', 'Carribean Trendline')
hold on
xlabel('Date');
ylabel('pH');
title('pH Levels Over Time');
legend('show'); % show legend with DisplayName labels
grid on; % enable grid
hold off



%%
addpath("C:\Users\woods\Downloads")

filename1 = 'hawaii_d90f_20ee_c4cb_e8d8_79ca_4468.nc';
%1a. Use the function "ncdisp" to display information about the data contained in this file
%-->
ncdisp(filename1)


%%
filename3 = 'NOAA_DHW_monthly_1a19_b145_a3b0_U1713885190139.nc';

%target_lat = 8.6167, % Replace with your desired latitude
%target_lon = -79.050; % Replace with your desired longitude
%target_lat = 12.6244
target_lon = -61.3506;
target_lat = 24.2282; % Replace with your desired latitude
%target_lon = 123.9218; 
cleanedlat = [8.6167,24.2282,12.6244];
cleanedlon = [-79.050,123.9218,-61.3506];

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
plotm(cleanedlat,cleanedlon,'m.','MarkerSize',20);
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
plotm(cleanedlat,cleanedlon,'m.','MarkerSize',20);
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
[~,lat_idx] = min(abs(target_lat-lat2));
[~,lon_idx] = min(abs(target_lon-lon2));

sst_target_data = sst2(lon_idx, lat_idx,:);
sst_target_data = squeeze(sst_target_data);
ssta_target_data = ssta2(lon_idx, lat_idx,:);
ssta_target_data = squeeze(ssta_target_data);

figure;
plot(time3, sst_target_data, 'b-', 'DisplayName', 'SST'); % Blue line for SST

x = datenum(time3);

p = polyfit(x,sst_target_data,1);

trendlinesst = polyval(p,x);

hold on;
plot(time3, trendlinesst, '-', 'DisplayName', 'SST Trendline'); % Plot the trendline
slope_SST = p(1); % Slope of the line

hold on;
plot(time3, ssta_target_data, 'r--', 'DisplayName', 'SSTA'); % Red dashed line for SSTA

p2 = polyfit(x,ssta_target_data,1);

trendlinessta = polyval(p2,x);

hold on;
plot(time3, trendlinessta, '-', 'DisplayName', 'SST Anomaly Trendline'); % Plot the trendline
slope_SSTA = p(1); % Slope of the line

xlabel('Time');
ylabel('Temperature (°C)');
title('Time Series of Sea Surface Temperature and Anomaly 1985 - Present');
legend() % Add a legend to identify the lines'
hold off;


%% Global_Coral_Bleaching_Database
addpath('/Users/amyz/Documents/Data Exploration/Final-Project/')
%addpath('C:\Users\woods\EnvDataExp\PartnerInfo\Final-Project\')
filename = 'Global_Coral_Bleaching_Database.csv';
stationdata = readtable(filename);

lat = table2array(stationdata(:,5));
lon = table2array(stationdata(:,6));
comments = table2array(stationdata(:,16));
precent_bleached = table2array(stationdata(:,11));
%% Raw Data
figure(1);clf
worldmap World
load coastlines
plotm(lat,lon,'m.','MarkerSize',10);
geoshow('landareas.shp','FaceColor','white')
title('Location for Coral Reef Bleaching')
%% Cleaned Data
%row_has_NA = any(strcmp(comments(:,1), 'N/A'));
rows = find(contains(precent_bleached,'N/A')|contains(precent_bleached,'>')|contains(precent_bleached,'<')|contains(precent_bleached,'%')|contains(precent_bleached,'0')|contains(precent_bleached,'-'));

goodIndices = setdiff(1:35053, rows);
%cleanedcomments = comments(goodIndices, :);
cleaned_precent = precent_bleached(goodIndices, :);

clean_precent_array = NaN(1,length(cleaned_precent));
for i = 1:length(cleaned_precent)
    clean_precent_array(i) = str2double(cell2mat(cleaned_precent(i)));
end

%array_precent = str2double(cell2mat(cleaned_precent));
cleanedlat = lat(goodIndices,:);
cleanedlon = lon(goodIndices,:);

figure(2); clf
worldmap world
load coastlines
plotm(coastlat,coastlon)

%plotm(cleanedlat,cleanedlon,cleaned_precent(:,1),'MarkerSize',10);
scatterm(cleanedlat,cleanedlon,50,clean_precent_array','filled','MarkerEdgeColor','k')

%geoshow('landareas.shp','FaceColor','white')
c = colorbar;
c.Title.String = '% bleached';
title('Location for Coral Reef Bleaching')

%% pie chart
rows_25 = find(clean_precent_array < 25);
rows_50 = find(clean_precent_array > 25 & clean_precent_array <= 50);
rows_75 = find(clean_precent_array > 50 & clean_precent_array <= 75);
rows_100 = find(clean_precent_array > 75 & clean_precent_array <= 100);
size_pie = [length(rows_25)/length(clean_precent_array)*100, length(rows_50)/length(clean_precent_array)*100,length(rows_75)/length(clean_precent_array)*100,length(rows_100)/length(clean_precent_array)*100];
labels_pie = {'25%', '50%', '75%', '100%'};
figure(4),clf
pie(size_pie, labels_pie);
title('% of Coral Bleech Pie Chart');

%% Location of most and least coral bleaching
max_region = max(clean_precent_array);
index_max = find(clean_precent_array == max_region);
lan_max = cleanedlat(index_max);
lon_max = cleanedlon(index_max);


min_region = min(clean_precent_array);
index_min = find(clean_precent_array == min_region);
%lower_lon = find(lon_min(index_min) < 0);

lan_min = cleanedlat(1597);
lon_min = cleanedlon(1597);
%1597

index_mid = find(min(abs(clean_precent_array - 50)));
lan_mid = cleanedlat(index_mid);
lon_mid = cleanedlon(index_mid);

index_aus = find(clean_precent_array == 51.52824717);
lan_aus = cleanedlat(index_aus);
lon_aus = cleanedlon(index_aus);

examine_regions = [max_region; min_region; clean_precent_array(index_mid)];
examine_lan = [lan_max;lan_min;lan_mid];
examine_lon = [lon_max;lon_min;lon_mid];

figure(3); clf
worldmap world
load coastlines
plotm(coastlat,coastlon)
scatterm(examine_lan,examine_lon,100,examine_regions,'filled','MarkerEdgeColor','k')

c = colorbar;
c.Title.String = 'precent_bleached(%)';
title('Location for Coral Reef Bleaching')








