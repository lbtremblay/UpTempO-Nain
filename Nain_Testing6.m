
clear
%% Global Constants
plotcolor = ['r'; 'y'; 'g'; 'b'; 'm'; 'k'];
depthprofile = [5, 20, 40, 60, 90, 100];

%% Import Data

%Import UpTempO and PAWS data into timetables
filename = 'SVP-I-UT UpTempo (6740)-300234065546740-20221206T222323UTC.csv';
TT = readtimetable(filename);
TT = removevars(TT, [1 2 3 4 5 12 13 14 15 16 17 18]);
TT = sortrows(TT, 'DataDate_UTC_', 'ascend');

filename = 'PAWS (6430)-300234065006430-20221206T223237UTC.csv';
WT = readtimetable(filename);
WT = removevars(WT, [1 2 3 4 14 15 16 17 18 19 20]);
WT = sortrows(WT, 'DataDate_UTC_', 'ascend');

%Import tide data into timetables
filename = 'tide.csv';
tide = readtimetable(filename);
filename = 'tide2.csv';
tide2 = readtimetable(filename);
filename = 'tide3.csv';
tide3 = readtimetable(filename);

%Import ERA5 data and isolate region around Nain
time = ncread('weatherdata.nc', 'time'); time = datetime(1900, 1, 1) + hours(time);
u = ncread('weatherdata.nc', 'u10'); u = reshape(u(2,2,:), [8664 1]);
v = ncread('weatherdata.nc', 'v10'); v = reshape(v(2,2,:), [8664 1]);
sw = ncread('weatherdata.nc', 'ssrd'); sw = reshape(sw(2,2,:), [8664 1]);
lw = ncread('weatherdata.nc', 'strd'); lw = reshape(lw(2,2,:), [8664 1]);
snow = ncread('snow.nc', 'sd'); snow = reshape(snow(2,2,:), [8664 1]);

%% Clean/Format Data

%Define date range for timeseries
range = datetime({'2020-02-21', '2020-05-22'});
timescale = timerange(range(1), range(2));
[~, i] = withinrange(TT, timescale);
t20 = cleantable(TT(i,:));
%Resample PAWS data and combine with UpTempO timetable
[~, i] = withinrange(WT, timescale);
tmp = WT(i,:);
i = removerep(tmp);
tmp(i,:) = [];
t20w = retime(tmp, 'hourly', 'spline');
t20 = synchronize(t20, t20w );

range = datetime({'2018-04-18', '2018-05-19'});
timescale = timerange(range(1), range(2));
[~, i] = withinrange(TT, timescale);
t18 = cleantable(TT(i,:));
[~, i] = withinrange(WT, timescale);
tmp = WT(i,:);
i = removerep(tmp);
tmp(i,:) = [];
t18w = retime(tmp, 'hourly', 'spline');
t18 = synchronize(t18, t18w );

%Add ERA5 data to timetable
i = find(time == t18.DataDate_UTC_(1) | time == t18.DataDate_UTC_(end));
t18 = addvars(t18, u(i(1):i(2)), v(i(1):i(2)), sw(i(1):i(2)), lw(i(1):i(2)), ...
    snow(i(1):i(2)), 'NewVariableNames', {'u', 'v', 'sw', 'lw', 'sd'});

i = find(time == t20.DataDate_UTC_(1) | time == t20.DataDate_UTC_(end));
t20 = addvars(t20, u(i(1):i(2)), v(i(1):i(2)), sw(i(1):i(2)), lw(i(1):i(2)), ...
    snow(i(1):i(2)), 'NewVariableNames', {'u', 'v', 'sw', 'lw', 'sd'});

t20d = retime(t20,'daily','mean');
t18d = retime(t18,'daily','mean');
%% Timeseries Plot
%Plotting temperature time series for each depth and year

%Define figure and axes 
figure()
temp_fig = tiledlayout(2,1);
ylabel(temp_fig,'Temperature (C)')
ax1=nexttile; ax2=nexttile; hold ([ax1 ax2], 'on')
temp_fig.TileSpacing = 'tight'; temp_fig.Padding = 'compact';
%Create colormap for plotting
C = linspecer(6, 'sequential'); C = flip(C);

%Plot each depth
for i = 1:6
    plot(ax1, t18.DataDate_UTC_, t18.(i), 'Color', C(i,:), 'LineWidth', 1)
    plot(ax2, t20.DataDate_UTC_, t20.(i), 'Color', C(i,:), 'LineWidth', 1)
end

%Define axis limits
xlim(ax1, [datetime('2018-02-21') datetime('2018-05-22')])
xlim(ax2, [datetime('2020-02-21') datetime('2020-05-22')])
ylim([ax1 ax2], [min(t20{:,1:6}, [], 'all') max(t20{:,1:6}, [], 'all')])

%Formatting
lg = legend(ax1, '5m', '20m', '40m', '60m', '90m', '100m');
lg.Location = 'northwest'; fontsize(lg, 20, 'points')
fontsize(temp_fig, 20, 'points')

saveas(gcf, 'figures/timeseries.png')

%% Daily Temperature Profiles
% Plotting vertical temperature profile of each day of a year

%Define figure and axes
figure()
profile_fig = (tiledlayout(1,2));
xlabel(profile_fig, 'Daily Mean Temp (C)'), ylabel(profile_fig, 'Depth (m)')
ax1=nexttile; ax2=nexttile; yticks([ax1 ax2], depthprofile)
hold ([ax1 ax2], 'on'), set([ax1 ax2], 'YDir','reverse')
profile_fig.Padding = 'compact'; profile_fig.TileSpacing = 'compact';

%Plot
years = {t18d, t20d}; ax = {ax1, ax2};
for i = 1:length(years)
    yr = years{i};
    C = linspecer(height(yr.(1)), 'sequential');
    for j = 1:height(yr.(1))
        plot(ax{i}, yr{j,1:6}, depthprofile, 'Color', C(j,:))
    end
end

%Formatting
ylim([ax1 ax2], [5 100])
fontsize(12, 'points')
saveas(gcf, 'figures/tempprofile.png')

%% Depth Correlation Scatter Plots
% Plotting correlation scatter plots between each layer
% Variables:
yr = t20;
%

figure(), corr_fig = tiledlayout(5,5);
corr_fig.Padding = 'compact'; corr_fig.TileSpacing = 'none';
%title(corr_fig, ['Scatterplot of Temperature Profiles ', int2str(year(yr.DataDate_UTC_(1)))])
xlabel(corr_fig, 'Temp (C)'), ylabel(corr_fig, 'Temp (C)')
in_use = [1,6,11,16,21,7,12,17,22,13,18,23,19,24,25];
index = 1;
for i = 1:6
    for j = i:6
        if i == j
            continue
        end
        nexttile(in_use(index)), hold on
        C = linspecer(length(yr.(i)), 'sequential');
        scatter(yr.(i), yr.(j), 15, C)
        f = polyfit(yr.(i), yr.(j),1);
        x = linspace(min(yr.(i)), max(yr.(i)), length(yr.(i)));
        y = polyval(f, x);
        plot(x,y, 'k')
        xlim([-1.8 -1.6]), ylim([-1.8 -1.6])
        xticks([-1.77 -1.65]), yticks([-1.77 -1.65])
        text(min(xlim)+0.01, max(ylim)-0.01, [int2str(depthprofile(i)), ' v. ', int2str(depthprofile(j))])
        if ismember(index,1:4)
            h = gca; h.XAxis.Visible = 'off';
        elseif ismember(index,[9,12,14,15])
            h = gca; h.YAxis.Visible = 'off';
        elseif index ~= 5
            h = gca; h.XAxis.Visible = 'off';
            h = gca; h.YAxis.Visible = 'off';
        end
        index = index+1;
    end
end
fontsize(15,"points")
text()
saveas(gcf, 'figures/depthscatter.png')

%% Fourier Analysis
% Fourier analysis of each depth time series
% Variables:
yr = t20;
%

%Detrend each temp timeseries
for i = 1:6
    yr.(i) = yr.(i)-movmean(yr.(i), 5*24);
end

%Fourier transform
Ts = 3600; %Period: 1 hour
fs = 1/Ts; %Frequency 1/h
fftTemp = fft(yr.Variables);
n = length(fftTemp);
fshift = (-n/2:n/2-1)*(fs/n);

%Plotting
figure(); 
fft_plot = tiledlayout(1,1); ax1 = nexttile;
hold(ax1, 'on'), xlim(ax1, [2 30]), ylim(ax1, [0 3])
set(ax1, 'XScale', 'log');
fft_plot.Padding = 'compact'; fft_plot.TileSpacing = 'none';
xlabel('Period (hours)'), ylabel('Magnitude'), xticks(0:3:30)

for i = 1:6
    semilogx(1./(Ts*fshift),abs(fftshift(fftTemp(:,i))), plotcolor(i), 'LineWidth', 1)
end

lg = legend('5m', '20m', '40m', '60m', '90m', '100m'); lg.Location = 'northwest';
fontsize(20,"points")
saveas(gcf, 'figures/fft.png')

%% Wavelet Analysis
% Wavelet analysis of each depth of one year
% Variables:
yr = t20;
%

figure(), wavelet_fig = tiledlayout(2,3);
wavelet_fig.Padding = 'compact'; wavelet_fig.TileSpacing = 'tight';
xlabel(wavelet_fig,'Hours'), ylabel(wavelet_fig,'Hours/cycle')
for i = 1:6
    [cfs, frq] = cwt(yr.(i), 'bump');
    nexttile
    n = 0:length(yr.(i))-1;
    pcolor(n,frq,abs(cfs)), shading interp
    xt = get(gca, 'YTick'); xp = round(1./(xt), 1);
    set(gca, 'YTick',xt, 'YTickLabel',xp)
    title([int2str(depthprofile(:,i)), 'm'])
end
fontsize(20,"points")
saveas(gcf, 'figures/wavelet.png')

%% Tide Plots

figure()
tide_fig = tiledlayout(2,2);
ax1=nexttile; ax3=nexttile; ax2=nexttile; ax4=nexttile;
tide_fig.TileSpacing = 'tight'; tide_fig.Padding = 'compact';

plot(ax1, t18.DataDate_UTC_, t18.(6)-movmean(t18.(6), 5*24), 'LineWidth', 1)
xlim(ax1, [datetime('2018-02-21') datetime('2018-05-22')])

ylim(ax1, [-0.02 0.02])
ylabel(ax1 ,'Temp Anomaly (C)')

yyaxis(ax1, 'right')
plot(ax1, tide2.date, tide2.(1), 'LineWidth', 1)
ylim(ax1, [-8.5 10.5])
ylabel(ax1,'Tide Height (m)')

plot(ax2, t20.DataDate_UTC_, t20.(6)-movmean(t20.(6), 5*24),'LineWidth', 1)
xlim(ax2, [datetime('2020-02-21') datetime('2020-05-22')])

ylim(ax2, [-0.02 0.02])
ylabel(ax2 ,'Temp Anomaly (C)')

yyaxis(ax2, 'right')
plot(ax2, tide.date, tide.(1), 'LineWidth', 1)
ylim(ax2, [-8.5 10.5])
ylabel(ax2,'Tide Height (m)')

x = t20.(6)(5:2165);
dx = x-movmean(x, 24*5);
y = tide.(1);

[c, lags] = xcorr(dx,y,12);
stem(ax4,lags, c)
xlabel(ax4, 'Lag (hr)'), ylabel(ax4, 'Correlation')

x = t18.(6)(5:725);
dx = x-movmean(x, 24*5);
y = tide2.(1);

[c, lags] = xcorr(dx,y,12);
stem(ax3,lags, c)
xlabel(ax3, 'Lag (hr)'), ylabel(ax3, 'Correlation')

grid([ax1 ax2 ax3 ax4], 'on')
fontsize(20,'points')
saveas(gcf, 'figures/lags.png')

%% Shortwave Down Plot

figure()
plot(t20.DataDate_UTC_, t20.(6))
ylabel('Temperature (C)')
yyaxis right
plot(t20d.DataDate_UTC_, t20d.(18)/(60*60), 'LineWidth', 1.5)
ylabel('SW Down (W/m^2)')

fontsize(20, 'points')

%% Snow Depth Plot

figure()
plot(t20.DataDate_UTC_, t20.(6))
ylabel('Temperature (C)')

yyaxis right
plot(t20.DataDate_UTC_, t20.(20))
ylabel('Snow Depth (m)')

title('Temperature (100m) and Snow Depth')
fontsize(13, 'points')

%% Hovmoller Diagram

figure()
hov_fig = tiledlayout(2,1);
ax1 = nexttile; ax2 = nexttile;
x = 1:length(t20d.DataDate_UTC_);
y = depthprofile';
z = t20d{:,1:6}';

contourcolor = [-1.778:0.001:-1.76, -1.76:0.01:-1.6];
contourcolor = [contourcolor, flip(contourcolor, 2)];
contourf(ax1,x,y,z,contourcolor)
set(ax1, 'YDir','reverse')
shading flat
title(ax1, 'Daily Average Temperature in Nain, 2020')
ylabel(ax1, 'Depth (m)')

cMap = turbo(256);
dataMax = max(z, [], 'all');
dataMin = min(z, [], 'all');

centerPoint = -1.77;
scalingIntensity = 5;

x = 1:length(cMap);
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity*x/max(abs(x));

x = sign(x).*exp(abs(x));
x = x - min(x);
x = x*511/max(x)+1;
newMap = interp1(x, cMap, 1:512);
colormap(newMap);

start = datetime(2020,4,18);%datetime(2020,3,19);
stop = datetime(2020,5,18);%datetime(2020,3,27);
num1 = days(datetime(start) - datetime(2020,2,20));
num2 = days(datetime(stop) - datetime(2020,2,20));
xlim(ax1,[num1 num2])
xt = num1:2:num2; xp = datetime(2020,2,21) + days(xt);
datelabel = string(datetime(xp, 'InputFormat', 'yyyy-MM-dd HH:mm', 'Format', 'dd-MM-yy'));
set(ax1, 'XTick',xt, 'XTickLabel', datelabel)
colorbar

x = 1:length(t18d.DataDate_UTC_);
y = depthprofile';
z = t18d{:,1:6}';
%figure()
contourf(ax2,x,y,z,contourcolor)
set(gca, 'YDir','reverse')
shading flat
title('Daily Average Temperature in Nain, 2018')
ylabel('Depth (m)')
clim([-1.7755   -1.6200])
xt = 1:2:31; xp = datetime(2018,4,18) + days(xt);
datelabel = string(datetime(xp, 'InputFormat', 'yyyy-MM-dd HH:mm', 'Format', 'dd-MM-yy'));
set(ax2, 'XTick',xt, 'XTickLabel', datelabel)
%colormap(newMap);

%% Air Temperature Plot
figure(), hold on
C = linspecer(6, 'sequential'); C = flip(C);

for i = 1:6
    plot(t20.DataDate_UTC_, t20.(i), 'Color', C(i,:), 'LineWidth', 1)
end
ylabel('Ocean Temp (C)')
yyaxis right, ylabel('Air Temp (C)')
plot(t20.DataDate_UTC_, t20.(10), 'k')
title('Ocean and Air Temperature Nain, 2020')
lg = legend('5m', '20m', '40m', '60m', '90m', '100m');
lg.Location = 'northwest';

%% NSIDC Sea Ice Concentration Plot
lati = ncread('seaice_conc_daily_nh_2020_v04r00.nc', 'latitude');
long = ncread('seaice_conc_daily_nh_2020_v04r00.nc', 'longitude');
conc = ncread('seaice_conc_daily_nh_2020_v04r00.nc', 'cdr_seaice_conc');
nsidctime = ncread('seaice_conc_daily_nh_2020_v04r00.nc', 'time');
[c, index] = min(sqrt((lati-56.4779).^2+(long-(-61.1298)).^2), [], 'all');

avgSIC = [];
num1 = days(datetime(2020,02,21) - datetime(2020,01,01));
num2 = days(datetime(2020,05,21) - datetime(2020,01,01));
conc = conc(:,:,num1:num2);

for q = 1:num2-num1+1
    layer = conc(:,:,q);
    avgSIC(end+1) = conc(113,378,q);
end

figure()
plot(t20.DataDate_UTC_, t20.(6))
ylabel('Temperature (C)')

yyaxis right
plot(t20d.DataDate_UTC_, avgSIC)
ylabel('SIC (%)')

title('Temperature (100m) and Sea Ice Concentration')
fontsize(13, 'points')

%% Functions

function outtable = cleantable(table)
    %Takes in a timetable containing hourly temperature data at 6 depths.
    %Returns a timetable after removing duplicate entries, resampling at
    %hour intervals, removing errors, and interpolating gaps in the data
    i = removerep(table);
    table(i,:) = [];
    outtable = retime(table, 'hourly', 'spline');
    maxlimit = -1.5; minlimit = -1.8;
    if year(table.DataDate_UTC_) == 2018
        maxlimit = -1.7; minlimit = -1.78;
    end
    for i = 1:6
        rowsToChange = (outtable.(i) >= maxlimit | outtable.(i) <= minlimit);
        if ~isempty(rowsToChange)
            outtable.(i)(rowsToChange) = nan;
        end
    end
    if year(table.DataDate_UTC_) ~= 2019
        outtable = fillmissing(outtable, 'linear');
    end
end

function index = removerep(table)
    %Takes in a timetable, returns the indices of every duplicate entry
    time = table.DataDate_UTC_;
    index = [];
    for i = 2:length(time)
        if time(i) == time(i-1)
            index(end+1) = i;
        end
    end
end

function [ax1, ax2] = correlate(table, var)
    %Takes in a timetable and an integer representing the position of a
    %variable to be compared
    %Plots a correlation plot between the 6th temperature timeseries and
    %the variable at the position entered. Returns the axes of the
    %subplots.
    table.Variables = table.Variables-movmean(table.Variables, 24*5);
    avgTemp = retime(table, 'daily', 'mean');
    avgTemp{:,1:6} = 1000*avgTemp{:,1:6};

    figure()
    weather_fig = tiledlayout(2, 1);
    weather_fig.Padding = 'compact'; weather_fig.TileSpacing = 'tight';
    ax1 = nexttile; ax2 = nexttile; hold([ax1 ax2], 'on'), grid(ax2, 'on')

    plot(ax1, avgTemp.DataDate_UTC_, avgTemp.(6), 'b')
    yyaxis(ax1, 'right'), plot(ax1, avgTemp.DataDate_UTC_, avgTemp.(var))

    scatter(avgTemp.(var), avgTemp.(6))
    r = corrcoef(avgTemp.(var), avgTemp.(6));
    r = round(r(2), 2);
    text('Units', 'Normalized', 'Position', [0.05, 0.85], 'string', ['r = ', num2str(r)], 'FontSize', 13 )


    f = polyfit(avgTemp.(var), avgTemp.(6), 1);
    x = linspace(min(avgTemp.(var)), max(avgTemp.(var)), length(avgTemp.(var)));
    y = polyval(f, x);
    plot(x,y, 'r')   
end
