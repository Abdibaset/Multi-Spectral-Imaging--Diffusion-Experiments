close all; clear; clc;
diary name_of_choice; diary on;                                    %saving the command window output to file

%providing the path
parentdir = fullfile("..", "MaestroData_fall22", folder_name);     %'folder_name' -> folder name can change to one required
acqfolders = dir(fullfile(parentdir, "acq*"));

% sorting folders in an ascending order  
fnum = zeros(1,numel(acqfolders));
for a = 1:numel(acqfolders)
    str_fld = strsplit(acqfolders(a).name,'ACQ');
    fnum(a) = str2double(str_fld{2});
end
[~,ind] = sort(fnum);
acqfolders = acqfolders(ind);

c = input('Select duration to analyze for:\n(1)for 60mins\n(2)for 3hrs\n');
t1 = 4;
% len = length(time_vect);
if c == 1
    %array for radius from all cross-sections -> 17*4
    radii_vect = zeros(1, (length(acqfolders)-5)*4);
    count = 0;                                                          %iterator to keep track of current column in radii_vec
    numfolders = numel(acqfolders)-4;
    increment_time = 0;                                                 %keeps track of time
    time_vect = zeros(1, (length(acqfolders)-5)*4);                     %4 slots for each time point
    
    %populating the vector with time
    while 1 
        increment_time = increment_time + 5;                           %increment current time by 5 for first 13 acquistio
        for t2 = 1:4
            time_vect(t1+t2) = increment_time;
        end 
        t1 = t1 + 4;
        if t1 >= (length(time_vect))
            break
        end 
    end
elseif c == 2
    numfolders = numel(acqfolders);
    %creating time vector 
    radii_vect = zeros(1, (length(acqfolders)-1)*4);
    count = 0;                                                          %iterator to keep track of current column in radii_ve
    increment_time = 0;                                                 %keeps track of time
    time_vect = zeros(1, (length(acqfolders)-1)*4);                     %4 slots for each time point
    
    %populating the vector with time
    while 1 
     if  t1 > (12*4)
           increment_time = increment_time + 30;                     %increment current time by 30mins
      else
         increment_time = increment_time + 5;                           %increment current time by 5 for first 13 acquistions
      end
    
        for t2 = 1:4
            time_vect(t1+t2) = increment_time;
        end 
        t1 = t1 + 4;
        if t1 == (length(time_vect))
            break
        end 
    end
end


%median arrays 
median_radii = zeros(1, (length(radii_vect)/4));
median_time = zeros(1, length(median_radii));
count_median = 0;
time_passed = 0;
t = 1;
while 1
    if t > 12
        time_passed = time_passed + 30;
    else 
        time_passed = time_passed + 5;
    end 
    t = t + 1;
    median_time(t) = time_passed;
    if t >= (length(median_time))
        break
    end

end
%target directory to create sample roi to reuse
imagedir = fullfile(parentdir, acqfolders(2).name);
acq_folder = dir(fullfile(imagedir, dye_name));
createSample_profileIM = fullfile(imagedir, acq_folder(1).name);
fprintf("Used file: %s for sample ROI\n", createSample_profileIM)       %printing the acqusition file used

%loop over all folders
for i = 1: numfolders
    %open folder
    imdirectory = fullfile(parentdir, acqfolders(i).name);
    im_toRead_dir = dir(fullfile(imdirectory, dye_name));
    im_toRead = fullfile(imdirectory, im_toRead_dir(1).name);
 
    
    %creating and saving rios
    for j = 1: 4
        %creating the sample roi is doesn't exist in target folder for
        %creating sample roi
        ROIs = fullfile(imagedir, sprintf("Sample_roi_%d.mat",j));
        if ~exist(ROIs, 'file')
            f1 = figure;
            imshow(createSample_profileIM); clim([0 50]); axis image; axis off;
            [cx1, cy1, p1, ix1, iy1] = improfile;
            save(ROIs, 'cx1', 'cy1', 'p1', 'ix1', 'iy1');
            close(f1);
        else
            load(ROIs, 'cx1', 'cy1', 'p1', 'ix1', 'iy1');               %load if sample roi exist 
        end
    
        %using the sample ROI to create a profile on other acquistions
        ROI_saved = fullfile(imdirectory, sprintf("roi_saved_%d.mat", j));
        if ~exist(ROI_saved, 'file')
            [cx2, cy2, p2, ix2, iy2] = improfile(imread(im_toRead), ix1, iy1);
            save(ROI_saved, 'cx2', 'cy2', 'p2', 'ix2', 'iy2');
        else 
            load(ROI_saved);
        end 

    %storing the data from every cross section of the well
    distance_vect = zeros(1, numel(cx2));
    midp_index = round(length(cx2)/2);
    midX = cx2(midp_index); midY = cy2(midp_index);

    %calculating distance using Euclidean equation
    for n = 1: numel(cx2)
        distance = ((cx2(n) - midX)^2 + (cy2(n) - midY)^2)^0.5;         %distance between two points
        if n < midp_index
            distance = -distance;
        end
        distance_vect(n) = distance/142.03;                             %converting from pixel to cm - 1cm -> 142.003
    end 
      
    %autofluorescence data
    if i == 1
        %profile through the horizontal center
        if j == 1
            center_autoData = zeros(length(cx2), 2);
            center_autoData(:, 1) = distance_vect;                      %column 1 - distance (cm)
            center_autoData(:, 2) = p2;                                 %column 2 - pixels 
            center_autoData = smoothdata(center_autoData);
        %profile through the left diagonal
        elseif j == 2
            leftDiagonal_autoData = zeros(length(cx2), 2);
            leftDiagonal_autoData(:, 1) = distance_vect;
            leftDiagonal_autoData(:, 2) = p2;
            leftDiagonal_autoData = smoothdata(leftDiagonal_autoData);
        %profile through the vertical center
        elseif j == 3
            vertical_autoData = zeros(length(cx2), 2);
            vertical_autoData(:, 1) = distance_vect;
            vertical_autoData(:, 2) = p2;
            vertical_autoData = smoothdata(vertical_autoData);
        %profile through the right diagonal of plate
        else
            rightDiagonal_autoData = zeros(length(cx2), 2);
            rightDiagonal_autoData(:, 1) = distance_vect;
            rightDiagonal_autoData(:, 2) = p2;
            rightDiagonal_autoData = smoothdata(rightDiagonal_autoData);
        end
    %fluorescence data 
    else
        %profile through horizontal center
        if j == 1
            center_floData = zeros(length(cx2), 2);
            center_floData(:, 1) = distance_vect;
            center_floData(:, 2) = p2;
            center_floData = smoothdata(center_floData);
            center_floData(:, 2) = center_floData(:, 2) - center_autoData(:, 2);
        %profile through the left diagonal
        elseif j == 2
            leftDiagonal_floData = zeros(length(cx2), 2);
            leftDiagonal_floData(:, 1) = distance_vect;
            leftDiagonal_floData(:, 2) = p2;
            leftDiagonal_floData = smoothdata(leftDiagonal_floData);
            leftDiagonal_floData(:, 2) = leftDiagonal_floData(:, 2) - leftDiagonal_autoData(:, 2);
        %profile through vertical center
        elseif j == 3
            vertical_floData = zeros(length(cx2), 2);
            vertical_floData(:, 1) = distance_vect;
            vertical_floData(:, 2) = p2;
            vertical_floData = smoothdata(vertical_floData);
            vertical_floData(:, 2) = vertical_floData(:, 2) - vertical_autoData(:, 2);
        %profile through the right diagonal
        else
            rightDiagonal_floData = zeros(length(cx2), 2);
            rightDiagonal_floData(:, 1) = distance_vect;
            rightDiagonal_floData(:, 2) = p2;
            rightDiagonal_floData = smoothdata(rightDiagonal_floData);
            rightDiagonal_floData(:, 2) = rightDiagonal_floData(:, 2) - rightDiagonal_autoData(:, 2);
        end
    end

     %plotting the profiles created to visualize 
      f3 = figure;      
      imshow(im_toRead); clim([0 50]); axis image; axis off;
      hold on
      plot( cx2, cy2, 'y--', 'LineWidth', 2.5);
      hold off
      profile = fullfile(imdirectory, sprintf("profile_%d.png", j));
      if ~exist(profile, 'file')
        saveas(gcf, fullfile(imdirectory, sprintf("profile_%d.png", j)), 'png');
      end 
    end
   close(f3);
    
    %plotting the autofluorescence graphs 
    f4 = figure;
    if i == 1
        plot(center_autoData(:, 1), center_autoData(:, 2), 'LineWidth', 2);
        hold on
        plot(leftDiagonal_autoData(:, 1), leftDiagonal_autoData(:, 2), 'LineWidth', 2);
        plot(vertical_autoData(:, 1), vertical_autoData(:, 2), 'LineWidth',1.5);
        plot(rightDiagonal_autoData(:, 1), rightDiagonal_autoData(:, 2), 'LineWidth', 2);
        hold off
        xlim([-2 2]); ylim([0 50]);
        legend({'c', 'ld', 'v', 'rd'});
        xlabel('distance(cm)'); ylabel('GrayScale/Pixel');
        title('Graph of autofluorescence');
    %plotting aggregate -> fluorescence - autofluorescence 
    elseif i > 1
       plot(center_floData(:, 1), center_floData(:, 2), 'LineWidth', 2);
       hold on
       plot(leftDiagonal_floData(:, 1), leftDiagonal_floData(:, 2), 'LineWidth', 2);
       plot(vertical_floData(:, 1), vertical_floData(:, 2), 'LineWidth', 2);
       plot(rightDiagonal_floData(:, 1), rightDiagonal_floData(:, 2), 'LineWidth', 2);
       hold off
       if contains(im_toRead, 'IR800') || contains(im_toRead, 'AF750')
            ylim([0 350]);                            %y-lim can changed depending on the strength of signal of dye
       else
           xlim([0 2500]);
       end
       xlim([-2 2]); 
       legend({'c', 'ld', 'v', 'rd'});
       xlabel('distance(cm)'); ylabel('GrayScale/Pixel');
       title('Graph of Fluorescence');
    end

    %saving the plots as a graph
    if i > 2
        graph = fullfile(imdirectory, sprintf('fluoresence_%d.fig', i));
        if ~exist(graph, 'file')
            saveas (gcf, fullfile(imdirectory, sprintf('fluoresence_%d.fig', i)), 'fig');
        end
    else
        autograph = fullfile(imdirectory, sprintf('autofluoresence_%d.fig', 1));
        if ~exist(autograph, 'file')
            saveas (gcf, fullfile(imdirectory, sprintf('autofluoresence.fig')), 'fig');
        end
    end 
    close(f4);

    %calculating the distance peaks - x-intercepts on either sides of the
    %peak
    if i > 1
        %per cross-section

        %arrays with derivates greater than zero for each array with
        %fluorescent data
        c_ind = find(abs(diff(center_floData(:, 2))) > 0);      
        ld_ind = find(abs(diff(leftDiagonal_floData(:, 2))) > 0);
        rd_ind = find(abs(diff(rightDiagonal_floData(:, 2))) > 0);
        v_ind = find(abs(diff(vertical_floData(:, 2))) > 0);
        
        cl_intercept = center_floData(c_ind(1), 1);
        cr_intercept = center_floData(c_ind(length(c_ind))+1, 1);   
        dist_center = cr_intercept -cl_intercept;                   %diameter        
        center_radiusSquared = (dist_center/2)^2;

        ldl_intercept = leftDiagonal_floData(ld_ind(1), 1);
        ldr_intercept = leftDiagonal_floData(ld_ind(length(ld_ind))+1, 1);
        dist_ld = ldr_intercept - ldl_intercept;
        ld_radiusSquared = (dist_ld/2)^2;

        rdl_intercept = rightDiagonal_floData(rd_ind(1), 1);
        rdr_intercept = rightDiagonal_floData(rd_ind(length(rd_ind))+1, 1);
        dist_rd = rdr_intercept - rdl_intercept;
        rd_radiusSquared = (dist_rd/2)^2;
        
        vl_intercept = vertical_floData(v_ind(1), 1);
        vr_intercept = vertical_floData(v_ind(length(v_ind))+1, 1);
        dist_v = vr_intercept - vl_intercept;
        v_radiusSquared = (dist_v/2)^2;
        ave_dist = (dist_v + dist_rd + dist_ld + dist_center)/4;
        
        %printing the diameter for storing the output 
        fprintf("\nThe diameter, radius squared from different cross-sections run: %d\n", i-1);
        fprintf("Center: %d(cm)\t\t\tC_r^2: %d(cm^2)\nLeft diagonal: %d(cm)\t\t\tld_r^2: %d(cm^2)\nVertical: %d(cm)\t\t\tv_r^2: %d(cm^2)\n" + ...
            "Right diagonal: %d(cm)\t\trd_r^2: %d(cm^2)\nAverage of the four: %d(cm)\n", ...
            dist_center,center_radiusSquared,dist_ld,ld_radiusSquared, dist_v, v_radiusSquared, dist_rd, rd_radiusSquared, ave_dist);
        
        %populating the radii_vect with radius squared from each
        %cross-section
        radii_vect(1, count+1) = center_radiusSquared;
        radii_vect(1, count+2) = ld_radiusSquared;
        radii_vect(1, count+3) = v_radiusSquared;
        radii_vect(1, count+4) = rd_radiusSquared;

        all_radii = zeros(1, 4);
        all_radii(1,1) = center_radiusSquared; all_radii(1, 2) = ld_radiusSquared;
        all_radii(1, 3) = v_radiusSquared; all_radii(1, 4) = rd_radiusSquared;
        
        median_radii(1, count_median+1) = median(all_radii);
        count_median = count_median + 1;
        count = count + 4; %updating the count 
    end 
 end 

%conslidates the radii and time vectors to one vector
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

if c == 1
    xlim([0 80]); ylim([0 1]);
    saveas(gcf, fullfile(imdirectory, sprintf('radii_time_60mins.fig')), 'fig');
else
    xlim([0 200]); ylim([0 1]);
    saveas(gcf, fullfile(imdirectory, sprintf('radii_time_3hrs.fig')), 'fig');
end 
title('Graph of radius^2(cm^2) vs time(min)');
close(figure(6));


figure(7);
median_coefficients = polyfit(median_time, median_radii, 1);
median_xfit = linspace(min(median_time), max(median_time), 1000);
median_yfit = polyval(median_coefficients, median_xfit);
scatter(median_time, median_radii);
hold on
plot(median_xfit, median_yfit, 'LineWidth', 2);
hold off 
if c == 1
    xlim([0, 80]); ylim([0, 1]);
    saveas(gcf, fullfile(imdirectory, sprintf('Median_radii_time_60mins.fig')), 'fig');
else 
    xlim([0, 200]); ylim([0, 1]);
    saveas(gcf, fullfile(imdirectory, sprintf('Median_radii_time.fig')), 'fig');
end 
xlabel('time(min)'); ylabel('r^2(cm^2)'); 
title('Median r^2 against time(mins)');
close(figure(7));


diff_cofficient = coefficients(1)*(1/60) * 0.25;
fprintf('\nDiffusion co-efficient: %d (cm^2 per second)\n', diff_cofficient)
medianDiff_cofficient = median_coefficients(1)*(1/60) * 0.25;                   %Einstein's equation
fprintf('Median radii Diffusion co-efficient: %d (cm^2 per second)\n', medianDiff_cofficient)
diary off   