
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

