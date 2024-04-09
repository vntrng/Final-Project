
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

filename = 'hawaii_d90f_20ee_c4cb_e8d8_79ca_4468.nc';
%1a. Use the function "ncdisp" to display information about the data contained in this file
%-->
ncdisp(filename)




%%
lat = ncread(filename,'latitude');
lon = ncread(filename,'longitude');

temperature = ncread(filename,'temp');

figure()
worldmap('World')

contourfm(lat,lon,temperature')
scatterm()
colorbar()

%% Global_Coral_Bleaching Database
%addpath('/Users/amyz/Documents/Data Exploration/Final-Project/')
filename = 'Global_Coral_Bleaching_Database.csv';
stationdata = readtable(filename);

lat = table2array(stationdata(:,5));
lon = table2array(stationdata(:,6));
comments = table2array(stationdata(:,16));
%%
figure(1);clf
worldmap World
load coastlines
plotm(lat,lon,'m.','MarkerSize',10);
geoshow('landareas.shp','FaceColor','white')
title('Location for Coral Reef Bleaching')
%%
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

