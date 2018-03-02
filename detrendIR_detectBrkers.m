
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to detrend Thermal Infrared images and isolate breaking waves.
%
% Created on 2 March 2018 by Roxanne J Carini at Applied Physics
% Laboratory, University of Washington.
% Contact: rjcarini@uw.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data

% Load IR images into array I (#rows x #cols x #frames):
% This file available in matlab_SampleData repository.
path = '/Users/rjcarini/Documents/GitHub/matlab_TestData/';
load([path,'sampleImageSeries.mat'])


%% Calculate trend or drift of pixel intensity over time

% Choose reference region centered on your area of interest for use in
% thresholding algorithm:
ref3d = I(160:189,60:89,:);
ref2d = reshape(ref3d,[size(ref3d,1)*size(ref3d,2),size(ref3d,3)]);

% Use detrend the data through two median-smoothing operations:
reftrend = nanmedian(ref2d);
smoothparam = round(numel(reftrend)/10);
reftrendsmooth = movmedian(reftrend,smoothparam);
reftrendsmooth = movmedian(reftrendsmooth,smoothparam); 
reftrendsmooth = reshape(reftrendsmooth,1,1,length(reftrendsmooth));
refdetrend = ref3d-repmat(reftrendsmooth,size(ref3d,1),size(ref3d,2));
reftrendmat = repmat(reftrendsmooth,size(I,1),size(I,2));

% Visualize trends:
figure('position',[50 50 600 400])
hold on
for i=1:length(timevec)
    p0 = plot(repmat(timevec(i),size(ref2d,1),1),ref2d(:,i),'.k');
end
p1 = plot(timevec,reftrend,'-r');
p2 = plot(timevec,squeeze(reftrendsmooth),'-c');
axis([timevec(1) timevec(end) 2e4 2.25e4])
grid on
ylabel('Pixel Intensity')
xlabel('Time')
datetick('x','keepticks','keeplimits')
legend([p0,p1,p2],'All reference data','Median at each time step','Smoothed reference level')


%% Run detection algorithm on raw and detrended pixel intensities

% Apply thresholding algorithm to raw images:
[Iraw,brkMask] = detectBreakers_singlecam(I,ref2d);

% Apply thresholding algorithm to detrended images:
Idetrend = I-reftrendmat;
[Idetrend,brkMask] = detectBreakers_singlecam(Idetrend,refdetrend);

%% Play IR imagery with mask overlay

figure('position',[50 50 600 600])
for i=1:length(timevec)
    clf
    h = imagesc(imrotate(I(:,:,i),180));
    colormap gray
    hold on
    RGB = cat(3,brkMask(:,:,i),zeros(size(brkMask,1),size(brkMask,2)));
    RGB = cat(3,RGB,zeros(size(brkMask,1),size(brkMask,2)));
    M = image(imrotate(RGB,180));
    alphadat = imrotate(RGB(:,:,1),180).*0.3;
    set(M,'AlphaData',alphadat);
    set(gca,'XTickLabel',[],'YTickLabel',[])
    axis equal
    axis tight
    title('TOWER IR')
    drawnow
end
