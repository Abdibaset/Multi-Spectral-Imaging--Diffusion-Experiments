%% open directory for maestr
close all; clear; clc;
parentdir = fullfile("..", "PearlData_Divya", "06.15.23-AF750+PD-L1", "every_15mins");
acqfolders = dir(fullfile(parentdir, "00*"));

diary (strcat(parentdir,  filesep(), 'diary')) 
%% imaging interval
fprintf('Remember 60mins imaging time means 5mins intervals\n');
interval = input('Enter imaging interval:');
fprintf('Interval choosen: %f\n', interval);

%% sorting folders by the order images were taken
fnum = zeros(1,numel(acqfolders));
for a = 1:numel(acqfolders)
    str_fld = strsplit(acqfolders(a).name,'_');
    fnum(a) = str2double(str_fld{1});
end
[~,ind] = sort(fnum);
acqfolders = acqfolders(ind);

%% debug
for ind=1: numel(acqfolders)
    fprintf("dir: %s\n", acqfolders(ind).name);
end 

%% create dir to save data 
recordsDate = datetime('today');
foldername = parentdir.split("/"); f = foldername(length(foldername));
datafld = fullfile(parentdir, "analysisResults/");
folder = strcat(datafld, "/", datestr(recordsDate));
image_fld = strcat(fullfile(parentdir, 'images/'), '/', datestr(recordsDate));
fprintf("%s\n", folder)
if ~exist(folder, 'dir')
    mkdir(folder);
end 

if ~exist(image_fld, 'dir')
    mkdir(image_fld);
end

numfiles = length(acqfolders);
allImages = zeros(482, 650, numfiles);
%% loop over all folders and apply ROI
for i = 1: numfiles
    %% open folder
    imagedir = fullfile(parentdir, acqfolders(i).name);
    imFiles = dir(fullfile(imagedir, "*800*"));           %specify the dye
    image = bfopen(convertStringsToChars(fullfile(imagedir, imFiles(1).name)));
    im_firstElement = image{1};
    im_secondElement = im_firstElement{1};
    allImages(:,:,i) = im_secondElement;
    f1 = figure;
    imshow(imadjust(im_secondElement, [0, 100/255])); colorbar; 
    saveas(gca, fullfile(image_fld, sprintf('iteration_%d.png', i)), 'png');
    close(f1);
 end 

%% halfMax full width implementation
pixels = squeeze(sum(allImages(:,:,:),1));
distance = 0:1:size(allImages,2)-1;
distance = distance.*0.01692307692;      % converting into cm - 1cm --> 142pixels.

for i = 1:numfiles
    pdata = pixels(:,i)-pixels(:,1);             %background substraction of pixels
    halfMax = (min(pdata) + max(pdata))/2;       %Find the half max value.
    index1 = find(pdata >= halfMax, 1, 'first'); %Find where the data first drops below half the max.
    index2 = find(pdata >= halfMax, 1, 'last');  %Find where the data last rises above half the max.
    fprintf("left:%d, right:%d\n", index1, index2);
    fullwidthHalfMax(i) = index2-index1 + 1;     %FWHM in indexes.
    pixels(:, i) = pdata;

    f1=figure;
    plot(distance, pdata, 'LineWidth', 2);
    hold on
    plot(distance(index1), pdata(index1), "*");
    plot(distance(index2), pdata(index2), "*");
    hold off 
    fullWidthHalfMaxDistance(i) = fullwidthHalfMax(i).*0.01692307692; % FWHM in cm
    legend({'curve', 'fwhm1', 'fwhm2'});
%     xlim([0 5]); ylim([0 6000]);
    saveas(gca, fullfile(folder, sprintf('acq%d.fig', i)), 'fig');
    title(['difference = ',num2str(fullWidthHalfMaxDistance(i)),' cm']);
    close(f1);
end

f2 = figure;
plot(distance, pixels, 'LineWidth',2);
xlabel('Distance (cm)')
ylabel('PixelIntensity')
saveas(gca, fullfile(folder, sprintf('allGraphs.fig')), 'fig');
title('All curves');
close(f2);
fprintf('diameter ');
disp(fullWidthHalfMaxDistance);
time = linspace(0,interval*60*(numfiles-1),numfiles);
radiusSquared = (fullWidthHalfMaxDistance/2).^2;
fprintf('radius squared ');
disp(radiusSquared)

linefit = polyfit(time(1:length(time)-1),radiusSquared(2:end),1);
pixelfit = linefit(1)*time(1:length(time)-1)+linefit(2);
f1=figure;
scatter(time(1:length(time)-1),radiusSquared(2:end));hold on
plot(time(1:length(time)-1),pixelfit,'blue');
sl = (0.106^2)/(4*(linefit(1)))/60;
fprintf('slope: %f\n', sl);
timeTaken = (0.106^2)/(4*(linefit(1)/4))/60;  %0.106 - depth for dye to diffuse through
title (['D-Coeff = ',num2str(linefit(1)/4),' cm^2/s', '  D-Time = ', num2str(timeTaken), ' minutes']);
saveas(gca, fullfile(folder, sprintf('radius^{2}_vs_time.fig')), 'fig');
for i = 1:length(pixelfit)
    fitannotation = sprintf('fity = %f', pixelfit(i));
    originalannotation = sprintf('initY = %f', radiusSquared(i+1));
    text(time(i), pixelfit(i), fitannotation, 'Color', 'red', 'FontSize', 12);
    text(time(i), radiusSquared(i+1), originalannotation, 'Color', 'black', 'FontSize', 12);
end
hold off
