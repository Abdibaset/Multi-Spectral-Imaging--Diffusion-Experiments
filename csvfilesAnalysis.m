clear; clc; close all;
parentDIR = fullfile("..", "MaestroData_fall22", "11.11.22-IR800");
acqfolders = dir(fullfile(parentDIR, "acq*"));

%% setting up a timetamp folder and a dairy to save outcome 
today = datetime('today');
date_string = datestr(today, 'mm.dd.yyyy');
date_string = strcat(date_string, ".fiji_mat");

csv_dir = fullfile(parentDIR, date_string);
if ~exist(csv_dir, 'dir')
    mkdir(csv_dir);
end
% diary(fullfile(csv_dir, strcat(date_string, ".txt"))); diary on; 

%% sort folders in chronological order
fnum = zeros(1, numel(acqfolders));
for i = 1:numel(acqfolders)
    str_folder = strsplit(acqfolders(i).name, "ACQ");
    fnum(i) = str2double(str_folder{2}); 
end
[~, ind] = sort(fnum);
acqfolders = acqfolders(ind);

%% prompt - analyze by 60mins or 3hrs
c = input('Duration to analyze? Select:\n(1)for 60mins\n(2)for 3hrs\n');
if c == 1
    xlsFile = fullfile(parentDIR, "time_60mins_run2.xlsx");
elseif c == 2
    xlsFile = fullfile(parentDIR, "time_pearl_3hrs.xlsx");
end
time = (readtable(xlsFile));
time = (table2array(time))';
time_vect = zeros(1, length(time));
time_vect = time;
radii_vect = zeros(1, length(time_vect));
count = 1;
if c == 1
    numfolders = (length(time_vect)/4) + 3;
elseif c == 2
    numfolders = length(acqfolders);
end 
count = 1;

%% autofluorescence data
%acquisition 1 - contains autofluorescent data
noiseFolder = fullfile(parentDIR, acqfolders(1).name);

%vertical line ROI data from autofluorescent data
noise_vertical = dir(fullfile(noiseFolder, "vertical*"));
noise_vertFile = fullfile(noiseFolder, noise_vertical(1).name);
vertical_auto = readmatrix(noise_vertFile);

%horizontal line ROI data from autofluorescence data
noise_horizontal = dir(fullfile(noiseFolder, "center*"));
noise_horiFile = fullfile(noiseFolder, noise_horizontal(1).name);
horizontal_auto = readmatrix(noise_horiFile);   

%left diagonal line ROI data from autofluorescence data
noise_lDiagonal = dir(fullfile(noiseFolder, "left*"));
noise_ldFile = fullfile(noiseFolder, noise_lDiagonal(1).name);
lDiagonal_auto = readmatrix(noise_ldFile);

%right diagonal line ROI data from autofluorescence dat
noise_rDiagonal = dir(fullfile(noiseFolder, "right*"));
noise_rdFile = fullfile(noiseFolder, noise_rDiagonal(1).name);
rDiagonal_auto = readmatrix(noise_rdFile);

figure(1);
plot(vertical_auto(:, 1), vertical_auto(:, 2), 'Color','red', 'LineWidth', 1.5);
hold on
plot(horizontal_auto(:, 1), horizontal_auto(:, 2), 'Color','blue', 'LineWidth', 1.5);
plot(lDiagonal_auto(:, 1), lDiagonal_auto(:, 2), 'Color','green', 'LineWidth', 1.5);
plot(rDiagonal_auto(:, 1), rDiagonal_auto(:, 2), 'Color','black', 'LineWidth', 1.5);
hold off
xlim([0 3]);
ylim([0, 30]);
xlabel("Distance(cm)");
ylabel("Pixel Intensity");
title("IRDye800 Autofluorescence");
legend({'Vertical', 'horizontal', 'ldiagonal', 'rdiagonal'}, 'Location','northeast');

acq_dir = fullfile(csv_dir, strcat("acq", num2str(1)));
if ~exist(acq_dir,  'dir')
    mkdir(acq_dir);
end 
saveas(gcf, fullfile(acq_dir, sprintf('autofluorescence_csv')), 'fig');


%% loop over all folders
for i = 4: numfolders
    floFolder = fullfile(parentDIR, acqfolders(i).name);
    acq_dir = fullfile(csv_dir, strcat('acq', num2str(i)));
    if ~exist(acq_dir,  'dir')
        new_acq_dir = mkdir(acq_dir);
    end 
    %vertical line ROI
    floVert = dir(fullfile(floFolder, "vertical*"));
    floVert_File = fullfile(floFolder, floVert(1).name);
    floVert_table = readmatrix(floVert_File);
    floVert_table(:, 2) = floVert_table(:, 2) - vertical_auto(:, 2);
    floVert_table = smoothdata(floVert_table, "gaussian", 15);

    %center line ROI
    floCenter = dir(fullfile(floFolder, "center*"));
    floCenter_File = fullfile(floFolder, floCenter(1).name);
    floCenter_table = readmatrix(floCenter_File);
    floCenter_table(:, 2) = floCenter_table(:, 2) - horizontal_auto(:, 2);
    floCenter_table = smoothdata(floCenter_table, "gaussian", 15);

    %left diagonal ROI
    flo_lDiagonal = dir(fullfile(floFolder, "left*"));
    flo_lDiagonal_File = fullfile(floFolder, flo_lDiagonal(1).name);
    flo_lDiagonal_table = readmatrix(flo_lDiagonal_File);
    flo_lDiagonal_table(:, 2) = flo_lDiagonal_table(:, 2) - lDiagonal_auto(:, 2);
    flo_lDiagonal_table = smoothdata(flo_lDiagonal_table, "gaussian", 15);

    %right diagonal ROI
    flo_rDiagonal = dir(fullfile(floFolder, "right*"));
    flo_rDiagonal_File = fullfile(floFolder, flo_rDiagonal(1).name);
    flo_rDiagonal_table = readmatrix(flo_rDiagonal_File);
    flo_rDiagonal_table(:, 2) = flo_rDiagonal_table(:, 2) - rDiagonal_auto(:, 2);
    flo_rDiagonal_table = smoothdata(flo_rDiagonal_table, "gaussian", 15);
    
    [leftslope_c, rightslope_c, radius_c] = generate_slopes(floCenter_table);
    [leftslope_v, rightslope_v, radius_v] = generate_slopes(floVert_table);
    [leftslope_ld, rightslope_ld, radius_ld] = generate_slopes(flo_lDiagonal_table);
    [leftslope_rd, rightslope_rd, radius_rd] = generate_slopes(flo_rDiagonal_table);
   
    figure(i);
    plot(floVert_table(:, 1), floVert_table(:, 2), 'LineWidth',1.5);
    hold on
    plot(leftslope_v(:, 1), leftslope_v(:, 2), 'LineWidth', 2);
    plot(rightslope_v(:, 1), rightslope_v(:, 2), 'LineWidth', 2);
    plot(floCenter_table(:, 1), floCenter_table(:, 2), 'LineWidth',1.5);
    plot(leftslope_c(:, 1), leftslope_c(:, 2), 'LineWidth', 2);     plot(rightslope_c(:, 1), rightslope_c(:, 2), 'LineWidth', 2);
    plot(flo_lDiagonal_table(:, 1), flo_lDiagonal_table(:, 2), 'LineWidth',1.5);
    plot(leftslope_ld(:, 1), leftslope_ld(:, 2), 'LineWidth', 2);
    plot(rightslope_ld(:, 1), rightslope_ld(:, 2), 'LineWidth', 2);
    plot(flo_rDiagonal_table(:, 1), flo_rDiagonal_table(:, 2), 'LineWidth', 1.5);
    plot(leftslope_rd(:, 1), leftslope_rd(:, 2), 'LineWidth', 2);
    plot(rightslope_rd(:, 1), rightslope_rd(:, 2), 'LineWidth', 2);
    hold off
    xlabel("Distance(cm)");
    ylabel("Pixel Instensity")
    title("IRDye 800 Fluorescence");
    legend({'Vertical', 'horizontal', 'lDiagonal', 'rDiagonal'}, 'Location','northeast');
    
    saveas(gcf, fullfile(acq_dir, sprintf('fluorescence')), 'fig');
    if contains(floFolder, 'IR800') || contains(floFolder, 'AF750')
       ylim([0 350]);                            %y-lim can changed depending on the strength of signal of dye
    elseif contains(floFolder, 'AF700')
       ylim([0 2000]);
    else
       ylim([0 5000]);
    end

%     radii_vect(count) = radius_c^2; 
    radii_vect(count) = radius_v^2; 
    radii_vect(count+1) = radius_rd^2;
    radii_vect(count+2) = radius_ld^2;
    count = count + 3;
    fprintf("\n==========================================================================\n");
    fprintf("\nRun: %d, corresponding to acq %d\n", i-1, i);
    fprintf("Center: %d(cm)\t\tC_r^2: %d(cm^2)\nLeft diagonal: %d(cm)\t\tld_r^2: %d(cm^2)\nVertical: %d(cm)\t\tv_r^2: %d(cm^2)\n" + ...
            "Right diagonal: %d(cm)\trd_r^2: %d(cm^2)\n", ...
            radius_c,radius_c^2,radius_ld,radius_ld^2, radius_v, radius_v^2, radius_rd, radius_rd^2);

%     fprintf("Left diagonal: %d(cm)\t\tld_r^2: %d(cm^2)\nVertical: %d(cm)\t\tv_r^2: %d(cm^2)\n" + ...
%             "Right diagonal: %d(cm)\trd_r^2: %d(cm^2)\n", ...
%             radius_ld,radius_ld^2, radius_v, radius_v^2, radius_rd, radius_rd^2);
end 


%% plot for the brownian motion equation - x^2 = 4Dt 
final_vec = zeros(length(time_vect), 2);
final_vec(:, 1) = time_vect;
final_vec(:, 2) = radii_vect;

figure(6);
coefficients = polyfit(final_vec(:, 1), final_vec(:, 2), 1);        %return 2*1 matrix- slope &y-intercept respectively
xfit = linspace(min(final_vec(:, 1)), max(final_vec(:, 1)), 1000);  %x-values to plot against
yfit = polyval(coefficients, xfit);                                 %evaluates polynomial p at each point in x. returns array of y-values                          


scatter(final_vec(:, 1), final_vec(:, 2));                          %scatter plots for all 4 cross-sections for visualization
hold on
plot(xfit, yfit, 'LineWidth',2);                                    %best line of fit plot
hold off
xlabel('time(min)'); ylabel('r^2(cm^2)'); 
title('Graph of radius^2(cm^2) vs time(min)');
if c == 1
    xlim([0 80]); ylim([0 0.01]);
    saveas(gcf, fullfile(acq_dir, sprintf('radii_time_60mins.fig')), 'fig');
else
    xlim([0 200]); ylim([0 1]);
    saveas(gcf, fullfile(acq_dir, sprintf('radii_time_3hrs.fig')), 'fig');
end 
close(figure(6));



diff_cofficient = coefficients(1)*(1/60) * 0.25;
fprintf('\nDiffusion co-efficient: %d (cm^2 per second)\n', diff_cofficient)
diary off   