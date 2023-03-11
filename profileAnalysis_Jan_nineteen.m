%% open directory for maestro
close all; clear; clc;
%diary 12.7.22-AF700.jan31.22; diary on;                                %saving the command window output to file

%providing the path
parentdir = fullfile("..", "MaestroData_fall22", "11.11.22-IR800");     %'folder_name' -> folder name can change to one required
acqfolders = dir(fullfile(parentdir, "acq*"));




%% sorting folders in an ascending order  
fnum = zeros(1,numel(acqfolders));
for a = 1:numel(acqfolders)
    str_fld = strsplit(acqfolders(a).name,'ACQ');
    fnum(a) = str2double(str_fld{2});
end
[~,ind] = sort(fnum);
acqfolders = acqfolders(ind);



%% prompt - analyze by 60mins or 3hrs
c = input('Select duration to analyze for:\n(1)for 60mins\n(2)for 3hrs\n');
if c == 1
    xlsFile = fullfile(parentdir, "time_60mins_run2.xlsx");
elseif c == 2
    xlsFile = fullfile(parentdir, "time_af700_3hrs_run3.xlsx");
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



%% select one well to create sample roi to reuse
imagedir = fullfile(parentdir, acqfolders(2).name);
acq_folder = dir(fullfile(imagedir, "*IR800*"));
createSample_profileIM = fullfile(imagedir, acq_folder(1).name);
fprintf("Used file: %s for sample ROI\n", createSample_profileIM)       %printing the acqusition file used



%% loop over all folders and apply RO
for i = 1: numfolders
    %% open folder
    imdirectory = fullfile(parentdir, acqfolders(i).name);
    im_toRead_dir = dir(fullfile(imdirectory, "*IR800*"));
    im_toRead = fullfile(imdirectory, im_toRead_dir(1).name);


    %% creating and saving sample for selected well
    for j = 1: 4
        ROI = fullfile(imagedir, sprintf("Sample_roi_%d.mat",j));
        if ~exist(ROI, 'file')
            f1 = figure;
            imshow(imread(createSample_profileIM)); clim([0 10]); axis image; axis off;
            [cx1, cy1, p1, ix1, iy1] = improfile;
            save(ROI, 'cx1', 'cy1', 'p1', 'ix1', 'iy1');
            close(f1);
        else
            load(ROI, 'cx1', 'cy1', 'p1', 'ix1', 'iy1');               %load if sample roi exists 
        end
    
    %% using the sample ROI to create a ROI on other acquisitions
    ROI_saved = fullfile(imdirectory, sprintf("roi_saved_%d.mat", j));
    if ~exist(ROI_saved, 'file')
        [cx2, cy2, p2, ix2, iy2] = improfile(imread(im_toRead), ix1, iy1);
        save(ROI_saved, 'cx2', 'cy2', 'p2', 'ix2', 'iy2');
    else 
        load(ROI_saved);
    end 

   


    %% calculating distance using euclidean equation
    distance_vect = zeros(1, numel(cx2));
    midp_index = round(length(cx2)/2);
    midX = cx2(midp_index); midY = cy2(midp_index);

    %calculating distance using Euclidean equation
    for n = 1: numel(cx2)
        distance = ((cx2(n) - midX)^2 + (cy2(n) - midY)^2)^0.5;         %distance between two points
        if n < midp_index
            distance = -distance;
        end
        distance_vect(n) = distance/1415;                             %converting from pixel to cm - 1cm -> 142.003
    end 


      
    %% autofluorescence & fluorescence data - raw and smoothened data 
    if i == 1
        %profile through the horizontal center
        if j == 1
            center_autoData = zeros(length(cx2), 2);
            center_autoData(:, 1) = distance_vect;                      %column 1-distance(cm)
            center_autoData(:, 2) = p2;                                 %column 2-pixels
            center_autoData = smoothdata(center_autoData, "gaussian", 5);
            
            center_rawdata_Auto = zeros(length(cx2), 2);
            center_rawdata_Auto(:, 1) = distance_vect;
            center_rawdata_Auto(:, 2) = p2;

        %profile through the left diagonal
        elseif j == 2
            leftDiagonal_autoData = zeros(length(cx2), 2);
            leftDiagonal_autoData(:, 1) = distance_vect;
            leftDiagonal_autoData(:, 2) = p2;

            leftDiagonal_rawdata_Auto = zeros(length(cx2), 2);
            leftDiagonal_rawdata_Auto(:, 1) = distance_vect;
            leftDiagonal_rawdata_Auto(:, 2) = p2;

        %profile through the vertical center
        elseif j == 3
            vertical_autoData = zeros(length(cx2), 2);
            vertical_autoData(:, 1) = distance_vect;
            vertical_autoData(:, 2) = p2;

            vertical_rawdata_Auto = zeros(length(cx2), 2);
            vertical_rawdata_Auto(:, 1) = distance_vect;
            vertical_rawdata_Auto(:, 2) = p2;
        %profile through the right diagonal of the well
        else
            rightDiagonal_autoData = zeros(length(cx2), 2);
            rightDiagonal_autoData(:, 1) = distance_vect;
            rightDiagonal_autoData(:, 2) = p2;

            rightDiagonal_rawdata_Auto = zeros(length(cx2), 2);
            rightDiagonal_rawdata_Auto(:, 1) = distance_vect;
            rightDiagonal_rawdata_Auto(:, 2) = p2;

        end

    else
        %profile through horizontal center
        if j == 1
            center_floData = zeros(length(cx2), 2);
            center_floData(:, 1) = distance_vect;
            center_floData(:, 2) = p2;
            center_floData(:, 2) = center_floData(:, 2) - center_autoData(:, 2);
            center_floData = smoothdata(center_floData, "gaussian", 15);

            center_rawdata_flo = zeros(length(cx2), 2);
            center_rawdata_flo(:, 1) = distance_vect;
            center_rawdata_flo(:, 2) = p2 - center_rawdata_Auto(:, 2);
        %profile through the left diagonal
        elseif j == 2
            leftDiagonal_floData = zeros(length(cx2), 2);
            leftDiagonal_floData(:, 1) = distance_vect;
            leftDiagonal_floData(:, 2) = p2;
            leftDiagonal_floData(:, 2) = leftDiagonal_floData(:, 2) - leftDiagonal_autoData(:, 2);
            leftDiagonal_floData = smoothdata(leftDiagonal_floData, "gaussian", 15);

            leftDiagonal_rawdata_flo = zeros(length(cx2), 2);
            leftDiagonal_rawdata_flo(:, 1) = distance_vect;
            leftDiagonal_rawdata_flo(:, 2) = p2 - leftDiagonal_rawdata_Auto(:, 2);
        %profile through vertical center
        elseif j == 3
            vertical_floData = zeros(length(cx2), 2);
            vertical_floData(:, 1) = distance_vect;
            vertical_floData(:, 2) = p2;
            vertical_floData(:, 2) = vertical_floData(:, 2) - vertical_autoData(:, 2);
            vertical_floData = smoothdata(vertical_floData, "gaussian", 15);

            vertical_rawdata_flo = zeros(length(cx2), 2);
            vertical_rawdata_flo(:, 1) = distance_vect;
            vertical_rawdata_flo(:, 2) = p2 - vertical_rawdata_Auto(:, 2);
        %profile through the right diagonal
        else
            rightDiagonal_floData = zeros(length(cx2), 2);
            rightDiagonal_floData(:, 1) = distance_vect;
            rightDiagonal_floData(:, 2) = p2;
            rightDiagonal_floData(:, 2) = rightDiagonal_floData(:, 2) - rightDiagonal_autoData(:, 2);
            rightDiagonal_floData = smoothdata(rightDiagonal_floData, "gaussian", 15);

            rightDiagonal_rawdata_flo = zeros(length(cx2), 2);
            rightDiagonal_rawdata_flo(:, 1) = distance_vect;
            rightDiagonal_rawdata_flo(:, 2) = p2 - rightDiagonal_rawdata_Auto(:, 2);
        end
    end

     %plotting the profiles created to visualize 
      figure(j);      
      imshow(imread(im_toRead)); clim([0 50]); axis image; axis off;
      hold on
      plot( cx2, cy2, 'y--', 'LineWidth', 2.5);
      hold off
      profile = fullfile(imdirectory, sprintf("profile_%d.png", j));
      if ~exist(profile, 'file')
        saveas(gcf, fullfile(imdirectory, sprintf("profile_%d.png", j)), 'png');
      end
      close(figure(j))
    end  

    
    %% plotting the autofluorescence graphs 
    f4 = figure;
    if i == 1
        plot(center_autoData(:, 1), center_autoData(:, 2), 'LineWidth', 2);
        hold on
        plot(leftDiagonal_autoData(:, 1), leftDiagonal_autoData(:, 2), 'LineWidth', 2);
        plot(rightDiagonal_autoData(:, 1), rightDiagonal_autoData(:, 2), 'LineWidth', 2);
        plot(vertical_autoData(:, 1), vertical_autoData(:, 2), 'LineWidth',2);
        hold off
        xlim([-2 2]); ylim([0 50]);
        legend({'center', 'ldiagonal', 'rdiagonal', 'vertical'});
        xlabel('distance(cm)'); ylabel('GrayScale/Pixel Value');
        title('Graph of autofluorescence');
        saveas (gcf, fullfile(imdirectory, sprintf('autofluoresence_%d.fig', i)), 'fig');
    end 


    %% calculating the distance peaks - x-intercepts on either sides of the peak
    if i > 3
    
       [leftslope_c, rightslope_c, radius_c] = generate_slopes(center_floData);
       [leftslope_v, rightslope_v, radius_v] = generate_slopes(vertical_floData);
       [leftslope_ld, rightslope_ld, radius_ld] = generate_slopes(leftDiagonal_floData);
       [leftslope_rd, rightslope_rd, radius_rd] = generate_slopes(rightDiagonal_floData);
       

    %% plot fluoresecence data 
       plot(center_floData(:, 1), center_floData(:, 2), 'LineWidth', 2, 'Color','blue');
       hold on
%        plot(center_rawdata_flo(:, 1), center_rawdata_flo(:, 2), 'LineWidth', 2);
       plot(leftslope_c(:, 1), leftslope_c(:, 2), 'LineWidth', 2);
       plot(rightslope_c(:, 1), rightslope_c(:, 2), 'LineWidth', 2);
       plot(leftDiagonal_floData(:, 1), leftDiagonal_floData(:, 2), 'LineWidth', 2, 'Color', 'g');
       plot(leftslope_ld(:, 1), leftslope_ld(:, 2), 'LineWidth', 2);
       plot(rightslope_ld(:, 1), rightslope_ld(:, 2), 'LineWidth', 2);
       plot(vertical_floData(:, 1), vertical_floData(:, 2), 'LineWidth', 2, 'Color', 'r');
       plot(leftslope_v(:, 1), leftslope_v(:, 2), 'LineWidth', 2);
       plot(rightslope_v(:, 1), rightslope_v(:, 2), 'LineWidth', 2);
       plot(rightDiagonal_floData(:, 1), rightDiagonal_floData(:, 2), 'LineWidth', 2, 'Color', 'c');
       plot(leftslope_rd(:, 1), leftslope_rd(:, 2), 'LineWidth', 2);
       plot(rightslope_rd(:, 1), rightslope_rd(:, 2), 'LineWidth', 2);
       hold off
      
       if contains(im_toRead, 'IR800') || contains(im_toRead, 'AF750')
            ylim([0 350]);                            %y-lim can changed depending on the strength of signal of dye
       elseif contains(im_toRead, 'AF700')
            ylim([0 2000]);
       else

           ylim([0 5000]);
       end
       xlim([-1.5 1.5]); 
       legend({'c', 'cr', 'ld', 'ldr', 'v', 'vr', 'rd', 'rdr'});
        xlabel('distance(cm)'); ylabel('GrayScale/Pixel');
       title('Graph of Fluorescence');
       saveas(gcf, fullfile(imdirectory, sprintf("fluorescence_fig_slopes%d",i)), 'fig');
       radii_vect(count) = radius_c^2; radii_vect(count+1) = radius_v^2; radii_vect(count+2) = radius_rd^2;
       radii_vect(count+3) = radius_ld^2;
       count = count + 4;
       fprintf("\n=====================================================================\n");
       fprintf("Run: %d\n", i-1);
        fprintf("Center: %d(cm)\t\tC_r^2: %d(cm^2)\nLeft diagonal: %d(cm)\t\tld_r^2: %d(cm^2)\nVertical: %d(cm)\t\tv_r^2: %d(cm^2)\n" + ...
            "Right diagonal: %d(cm)\trd_r^2: %d(cm^2)\n", ...
            radius_c,radius_c^2,radius_ld,radius_ld^2, radius_v, radius_v^2, radius_rd, radius_rd^2);
    end

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
    xlim([0 80]); ylim([0 0.2]);
    saveas(gcf, fullfile(imdirectory, sprintf('radii_time_60mins.fig')), 'fig');
else
    xlim([0 200]); ylim([0 1]);
    saveas(gcf, fullfile(imdirectory, sprintf('radii_time_3hrs.fig')), 'fig');
end 
close(figure(6));



diff_cofficient = (coefficients(1)/4)/60;
fprintf('\nDiffusion co-efficient: %d (cm^2 per second)\n', diff_cofficient)
diary off   