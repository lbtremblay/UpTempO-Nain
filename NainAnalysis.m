clear
%% Global Constants
depthprofile = [5, 20, 40, 60, 90, 100];
colors = linspecer(6, 'sequential'); colors = flip(colors);

%% Import Data
%Import UpTempO
filename = 'SVP-I-UT UpTempo (6740)-300234065546740-20221206T222323UTC.csv';
TT = readtimetable(filename);
TT = removevars(TT, [1:5, 12:18]);
TT = sortrows(TT, 'DataDate_UTC_', 'ascend');

%Import PAWS
filename = 'PAWS (6430)-300234065006430-20221206T223237UTC.csv';
WT = readtimetable(filename);
WT = removevars(WT, [1:6, 8:20]);
WT = sortrows(WT, 'DataDate_UTC_', 'ascend');

%Import tide data
filename = 'tideheight.csv';
tide = readtimetable(filename);

%Import ERA5 data and isolate region around Nain
time = ncread('era5.nc', 'time'); time = datetime(1900, 1, 1) + hours(time);
u = ncread('era5.nc', 'u10'); u = reshape(u(2,2,:), [8664 1]);
v = ncread('era5.nc', 'v10'); v = reshape(v(2,2,:), [8664 1]);
sw = ncread('era5.nc', 'ssrd'); sw = reshape(sw(2,2,:), [8664 1]);
lw = ncread('era5.nc', 'strd'); lw = reshape(lw(2,2,:), [8664 1]);
snow = ncread('era5.nc', 'sd'); snow = reshape(snow(2,2,:), [8664 1]);

%Import NSIDC SIC data
conc20 = ncread('seaice_conc_daily_nh_2020_v04r00.nc', 'cdr_seaice_conc');
conc18 = ncread('seaice_conc_daily_nh_2018_v04r00.nc', 'cdr_seaice_conc');

%% Clean/Format Data

%Define date range for timeseries
range = datetime({'2020-02-21', '2020-05-22'});
timescale = timerange(range(1), range(2));
[~, i] = withinrange(TT, timescale);
t20 = cleantable(TT(i,:));

%Resample PAWS/tide data and combine with UpTempO timetable
[~, i] = withinrange(WT, timescale);
tmp = WT(i,:);
i = removerep(tmp);
tmp(i,:) = [];
temptable = retime(tmp, 'hourly', 'spline');
t20 = synchronize(t20, temptable);

[~, i] = withinrange(tide, timescale);
t20 = synchronize(t20, tide(i,:));

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

[~, i] = withinrange(tide, timescale);
t18 = synchronize(t18, tide(i,:));

%Add ERA5 data to timetable
i = find(time == t18.DataDate_UTC_(1) | time == t18.DataDate_UTC_(end));
t18 = addvars(t18, u(i(1):i(2)), v(i(1):i(2)), sw(i(1):i(2)), lw(i(1):i(2)), ...
    snow(i(1):i(2)), 'NewVariableNames', {'u', 'v', 'sw', 'lw', 'sd'});

i = find(time == t20.DataDate_UTC_(1) | time == t20.DataDate_UTC_(end));
t20 = addvars(t20, u(i(1):i(2)), v(i(1):i(2)), sw(i(1):i(2)), lw(i(1):i(2)), ...
    snow(i(1):i(2)), 'NewVariableNames', {'u', 'v', 'sw', 'lw', 'sd'});

%Create daily timetables
t20d = retime(t20,'daily','mean');
t18d = retime(t18,'daily','mean');

%Add NSIDC SIC data to daily timetable
sic = zeros(1,height(t20d));
num1 = days(datetime(2020,02,21) - datetime(2020,01,01));
num2 = days(datetime(2020,05,21) - datetime(2020,01,01));

for q = num1:num2
    layer = conc20(:,:,q);
    sic(num2+1-q) = conc20(113,378,q);
end
t20d = addvars(t20d, flipud(sic'), 'NewVariableNames', 'sic');

sic = zeros(1,height(t18d));
num1 = days(datetime(2018,04,18) - datetime(2018,01,01));
num2 = days(datetime(2018,05,18) - datetime(2018,01,01));

for q = num1:num2
    layer = conc18(:,:,q);
    sic(num2+1-q) = conc18(113,378,q);
end
t18d = addvars(t18d, flipud(sic'), 'NewVariableNames', 'sic');

%% Timeseries Plot
%Plotting temperature time series for each depth and year

%Define figure and axes 
figure()
temp_fig = tiledlayout(2,1);
ylabel(temp_fig,'Temperature (C)')
ax1=nexttile; ax2=nexttile; hold ([ax1 ax2], 'on')
temp_fig.TileSpacing = 'tight'; temp_fig.Padding = 'compact';
title(temp_fig, 'Temperature with Depth, 2018/2020')

%Plot
for i = 1:6
    plot(ax1, t18.DataDate_UTC_, t18.(i), 'Color', colors(i,:), 'LineWidth', 1)
    plot(ax2, t20.DataDate_UTC_, t20.(i), 'Color', colors(i,:), 'LineWidth', 1)
end

%Define axis limits
xlim(ax1, [datetime('2018-02-21') datetime('2018-05-22')])
xlim(ax2, [datetime('2020-02-21') datetime('2020-05-22')])
ylim([ax1 ax2], [min(t20{:,1:6}, [], 'all') max(t20{:,1:6}, [], 'all')])

%Formatting
lg = legend(ax1, '5m', '20m', '40m', '60m', '90m', '100m');
fontsize(temp_fig, 15, 'points')
lg.Location = 'northwest'; fontsize(lg, 10, 'points')

saveas(gcf, 'figures/timeseries.png')

%% Daily Temperature Profiles
% Plotting vertical temperature profile of each day of a year

%Define figure and axes
figure()
profile_fig = (tiledlayout(1,2));
xlabel(profile_fig, 'Daily Mean Temp (C)'), ylabel(profile_fig, 'Depth (m)')
ax1=nexttile; ax2=nexttile; yticks([ax1 ax2], depthprofile)
hold ([ax1 ax2], 'on'), set([ax1 ax2], 'YDir','reverse')
title(profile_fig, 'Daily Vertical Temp Profiles, 2018/2020')
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

%Format
ylim([ax1 ax2], [5 100])
fontsize(13, 'points')
saveas(gcf, 'figures/tempprofile.png')

%% Depth Correlation Scatter Plots
% Plotting correlation scatter plots between each layer
% Variables:
yr = t20;
%

%Figure/Axes
figure(), corr_fig = tiledlayout(5,5);
corr_fig.Padding = 'compact'; corr_fig.TileSpacing = 'none';
title(corr_fig, ['Scatterplot of Temperature Profiles, ', int2str(year(yr.DataDate_UTC_(1)))])
xlabel(corr_fig, 'Temp (C)'), ylabel(corr_fig, 'Temp (C)')

%Plot
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

fontsize(11,"points")
saveas(gcf, 'figures/depthscatter.png')

%% Fourier Analysis
% Fourier analysis of each depth time series
% Year:
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

%Figure/Axes
figure(); 
fft_plot = tiledlayout(1,1); ax1 = nexttile;
hold(ax1, 'on'), xlim(ax1, [2 30]), ylim(ax1, [0 3])
set(ax1, 'XScale', 'log');
fft_plot.Padding = 'compact'; fft_plot.TileSpacing = 'none';
xlabel('Period (hours)'), ylabel('Magnitude'), xticks(0:3:24)
title('Fourier of Temp Timeseries, 2020')

%Plot
for i = 1:6
    semilogx(1./(Ts*fshift),abs(fftshift(fftTemp(:,i))), 'Color', colors(i,:), 'LineWidth', 1)
end

%Format
lg = legend('5m', '20m', '40m', '60m', '90m', '100m'); lg.Location = 'northwest';
fontsize(12,"points")
saveas(gcf, 'figures/fft.png')

%% Wavelet Analysis
% Wavelet analysis of each depth of one year
% Year:
yr = t20;
%

%Figure/axes
figure(), wavelet_fig = tiledlayout(2,3);
wavelet_fig.Padding = 'compact'; wavelet_fig.TileSpacing = 'tight';
xlabel(wavelet_fig,'Hours'), ylabel(wavelet_fig,'Hours/cycle')
title(wavelet_fig, 'Wavelet Analysis, 2020')

%Plot
for i = 1:6
    [cfs, frq] = cwt(yr.(i), 'bump');
    nexttile
    n = 0:length(yr.(i))-1;
    pcolor(n,frq,abs(cfs)), shading interp
    xt = get(gca, 'YTick'); xp = round(1./(xt), 1);
    set(gca, 'YTick',xt, 'YTickLabel',xp)
    title([int2str(depthprofile(:,i)), 'm'])
end

%Define Colormap
cMap = turbo(256);
dataMax = max(abs(cfs), [], 'all');
dataMin = min(abs(cfs), [], 'all');

centerPoint = 3e-3;
scalingIntensity = 4;

x = 1:length(cMap);
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity*x/max(abs(x));

x = sign(x).*exp(abs(x));
x = x - min(x);
x = x*511/max(x)+1;
newMap = interp1(x, cMap, 1:512);
colormap(newMap);

%Format
fontsize(15,"points")
saveas(gcf, 'figures/wavelet.png')

%% Tide Plots

%Figure/axes
figure()
tide_fig = tiledlayout(2,2);
colororder([0 0.4470 0.7410;0.6350 0.0780 0.1840])
title(tide_fig, 'Detrended Temp vs Tide Height, 2018/2020')
ax1=nexttile; ax3=nexttile; ax2=nexttile; ax4=nexttile;
tide_fig.TileSpacing = 'tight'; tide_fig.Padding = 'compact';

%Plot/format
plot(ax1, t18.DataDate_UTC_, t18.(6)-movmean(t18.(6), 5*24), 'LineWidth', 1)
xlim(ax1, [datetime('2018-02-21') datetime('2018-05-22')])

ylim(ax1, [-0.02 0.02])
ylabel(ax1 ,'Temp Anomaly (C)')

yyaxis(ax1, 'right')
plot(ax1, t18.DataDate_UTC_, t18.(8), 'LineWidth', 1)
ylim(ax1, [-8.5 10.5])
ylabel(ax1,'Tide Height (m)')

plot(ax2, t20.DataDate_UTC_, t20.(6)-movmean(t20.(6), 5*24),'LineWidth', 1)
xlim(ax2, [datetime('2020-02-21') datetime('2020-05-22')])

ylim(ax2, [-0.02 0.02])
ylabel(ax2 ,'Temp Anomaly (C)')

yyaxis(ax2, 'right')
plot(ax2, t20.DataDate_UTC_, t20.(8), 'LineWidth', 1)
ylim(ax2, [-8.5 10.5])
ylabel(ax2,'Tide Height (m)')

x = t20.(6)(5:2165);
dx = x-movmean(x, 24*5);
y = t20.(8);

[c, lags] = xcorr(dx,y,12);
stem(ax4,lags, c)
xlabel(ax4, 'Lag (hr)'), ylabel(ax4, 'Correlation')

x = t18.(6)(5:725);
dx = x-movmean(x, 24*5);
y = t18.(8);

[c, lags] = xcorr(dx,y,12);
stem(ax3,lags, c)
xlabel(ax3, 'Lag (hr)'), ylabel(ax3, 'Correlation')

grid([ax1 ax2 ax3 ax4], 'on')
fontsize(13,'points')
saveas(gcf, 'figures/lags.png')

%% Shortwave Down Plot

figure(), swd_fig = tiledlayout(2,1);
ax1 = nexttile; ax2 = nexttile;
swd_fig.TileSpacing = 'tight'; swd_fig.Padding = 'compact';
colororder([0 0.4470 0.7410;0.6350 0.0780 0.1840])

yyaxis(ax1, 'left')
plot(ax1, t20.DataDate_UTC_, t20.(6))
ylabel(ax1, 'Temperature (C)')

yyaxis(ax1, 'right')
plot(ax1, t20d.DataDate_UTC_, t20d.(11)/(60*60), 'LineWidth', 1.5)
ylabel(ax1, 'SW Down (W/m^2)')

yyaxis(ax2, 'left')
plot(ax2, t18.DataDate_UTC_, t18.(6))
ylabel(ax2, 'Temperature (C)')

yyaxis(ax2, 'right')
plot(ax2, t18d.DataDate_UTC_, t18d.(11)/(60*60), 'LineWidth', 1.5)
ylabel(ax2, 'SW Down (W/m^2)')

%Define axis limits
xlim(ax1, [datetime('2020-02-21') datetime('2020-05-22')])
xlim(ax2, [datetime('2018-02-21') datetime('2018-05-22')])
yyaxis(ax1, 'left')
yyaxis(ax2, 'left')
ylim([ax1 ax2], [min(t20{:,6}, [], 'all') max(t20{:,6}, [], 'all')])

title(swd_fig, 'Temperature (100m) and SWD')
fontsize(12, 'points')

%% Snow Depth Plot

figure(), sd_fig = tiledlayout(2,1);
ax1 = nexttile; ax2 = nexttile;
sd_fig.TileSpacing = 'tight'; sd_fig.Padding = 'compact';
colororder([0 0.4470 0.7410;0.6350 0.0780 0.1840])

yyaxis(ax1, 'left')
plot(ax1, t20.DataDate_UTC_, t20.(6))
ylabel(ax1, 'Temperature (C)')

yyaxis(ax1, 'right')
plot(ax1, t20.DataDate_UTC_, t20.(13), 'LineWidth', 1.5)
ylabel(ax1, 'Snow Depth (m)')

yyaxis(ax2, 'left')
plot(ax2, t18.DataDate_UTC_, t18.(6))
ylabel(ax2, 'Temperature (C)')

yyaxis(ax2, 'right')
plot(ax2, t18.DataDate_UTC_, t18.(13), 'LineWidth', 1.5)
ylabel(ax2, 'Snow Depth (m)')

%Define axis limits
xlim(ax1, [datetime('2020-02-21') datetime('2020-05-22')])
xlim(ax2, [datetime('2018-02-21') datetime('2018-05-22')])
yyaxis(ax1, 'right'), yyaxis(ax2, 'right')
ylim([ax1 ax2], [min(t20{:,13}, [], 'all') max(t18{:,13}, [], 'all')])
yyaxis(ax1, 'left'), yyaxis(ax2, 'left')
ylim([ax1 ax2], [min(t20{:,6}, [], 'all') max(t20{:,6}, [], 'all')])

title(sd_fig, 'Temperature (100m) and Snow Depth')
fontsize(12, 'points')

%% Hovmoller Diagram

figure()
hov_fig = tiledlayout(2,1);
ax1 = nexttile; ax2 = nexttile;

%Plot 1st Hovmoller
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
hov_fig.TileSpacing = 'tight'; hov_fig.Padding = 'compact';

%Define Colormap
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

%Plot 2nd Hovmoller
x = 1:length(t20d.DataDate_UTC_);
y = depthprofile';
z = zeros(6,91);
z(z==0) = NaN;
z(1:6, 58:88) = t18d{:,1:6}';
contourf(ax2,x,y,z,contourcolor)
set(gca, 'YDir','reverse')
shading flat

%Axes
set(ax1,'xtick',[])
start = datetime(2018,2,21);
stop = datetime(2018,5,21);
num1 = days(datetime(start) - datetime(2018,2,20));
num2 = days(datetime(stop) - datetime(2018,2,20));
xlim(ax1,[num1 num2])
xt = num1:14:num2; xp = datetime(2018,2,20) + days(xt);
datelabel = string(datetime(xp, 'InputFormat', 'yyyy-MM-dd HH:mm', 'Format', 'MMM dd'));
set(ax2, 'XTick',xt, 'XTickLabel', datelabel)

%Figure/format
title('Daily Average Temperature in Nain 2018')
ylabel('Depth (m)')
clim([-1.7755   -1.6200])

fontsize(13, 'points')

%% Air Temperature Plot

figure(), hold on
colororder([0.6350 0.0780 0.1840;0 0.4470 0.7410])

%Plot
for i = 1:6
    plot(t20.DataDate_UTC_, t20.(i), 'Color', colors(i,:), 'LineWidth', 1)
end

%Format
ylabel('Ocean Temp (C)')
yyaxis right, ylabel('Air Temp (C)')
plot(t20.DataDate_UTC_, t20.(7), '--')
title('Ocean and SST Nain, 2020')
lg = legend('5m', '20m', '40m', '60m', '90m', '100m');
lg.Location = 'northwest';

fontsize(14, 'points')

%% NSIDC Sea Ice Concentration Plot

figure(), sic_fig = tiledlayout(2,1);
ax1 = nexttile; ax2 = nexttile;
sic_fig.TileSpacing = 'tight'; sic_fig.Padding = 'compact';
colororder([0 0.4470 0.7410;0.6350 0.0780 0.1840])

yyaxis(ax1, 'left')
plot(ax1, t20.DataDate_UTC_, t20.(6))
ylabel(ax1, 'Temperature (C)')

yyaxis(ax1, 'right')
plot(ax1, t20d.DataDate_UTC_, t20d.(14), 'LineWidth', 1.5)
ylabel(ax1, 'SIC (%)')

yyaxis(ax2, 'left')
plot(ax2, t18.DataDate_UTC_, t18.(6))
ylabel(ax2, 'Temperature (C)')

yyaxis(ax2, 'right')
plot(ax2, t18d.DataDate_UTC_, t18d.(14), 'LineWidth', 1.5)
ylabel(ax2, 'SIC (%)')

%Define axis limits
xlim(ax1, [datetime('2020-02-21') datetime('2020-05-22')])
xlim(ax2, [datetime('2018-02-21') datetime('2018-05-22')])
yyaxis(ax1, 'right'), yyaxis(ax2, 'right')
ylim([ax1 ax2], [min(t18d{:,14}, [], 'all') max(t18d{:,14}, [], 'all')])
yyaxis(ax1, 'left'), yyaxis(ax2, 'left')
ylim([ax1 ax2], [min(t20{:,6}, [], 'all') max(t20{:,6}, [], 'all')])

title(sic_fig, 'Temperature (100m) and Sea Ice Concentration')
fontsize(12, 'points')

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