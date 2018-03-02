% Load IR images into array I (rows x cols x frames):
load('sampleImageSeries.mat')

% Choose reference region centered on your area of interest for use in
% thresholding algorithm:
ref = I(310:339,210:239,:);
ref = reshape(ref,[size(ref,1)*size(ref,2),size(ref,3)]);

% Use detrend the data through two median-smoothing operations:
reftrend = nanmedian(ref);
smoothparam = round(numel(reftrend)/10);
reftrendsmooth = movmedian(reftrend,smoothparam);
reftrendsmooth = movmedian(reftrendsmooth,smoothparam); 
reftrendmat = repmat(reftrendsmooth,size(I,1),size(I,2));
refdetrend = ref-reftrendmat;

% Visualize trends:
figure('position',[50 50 800 300])
hold on
for i=1:length(timevec)
    p0 = plot(repmat(timevec(i),size(ref,1),1),ref(:,i),'.k');
end
p1 = plot(timevec,reftrend,'-r');
p2 = plot(timevec,reftrendsmooth,'-c');
axis([timevec(1) timevec(end) 2e4 2.25e4])
grid on
ylabel('Pixel Intensity')
xlabel('Time')
datetick('x','keepticks','keeplimits')
legend([p0,p1,p2],'All reference data','Median at each time step','Smoothed reference level')

% Apply thresholding algorithm to raw images:
[Iraw,brkMask] = detectBreakers_singlecam(I,ref);


% Apply thresholding algorithm to detrended images:
Idetrend = I-reftrendmat;
[Idetrend,brkMask] = detectBreakers_singlecam(Idetrend,refdetrend);

