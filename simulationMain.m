close all
clear all
%% generate default image
xPix = 512; % dimension of image in x
yPix = 512; % dimension of image in y
pixS = 0.1; % pixel size in nm
centerX = xPix / 2;
centerY = yPix / 2;

%% protein structre
r = 20; % radius nanopore in nm
N = 6; % number of subunits (labelled)

%% donut properties
lambda = 488; % wavelength [nm]
NA = 1.4; % numerical aperture
I0 = 100; % laser intensity at max of Gauss
doN = 360; % how many points for donut circle movement

%% generate protein structure
protImg = zeros(xPix, yPix); 
% angles for N points
angles = linspace(0, 2*pi, N+1);
angles(end) = []; % remove last point (= first point)
% set points on circle
for i = 1:N
    x = round(centerX + r/pixS * cos(angles(i)));
    y = round(centerY + r/pixS * sin(angles(i)));
    if x > 0 && x <= xPix && y > 0 && y <= yPix
        protImg(y, x) = 1;
    end
end


%% generate donuts with different positions
[X, Y] = meshgrid(1:xPix, 1:yPix); % meshgrid
w0 = lambda/(pi*NA)/pixS; % Gauss
doR = r/pixS;
doAngles = linspace(0, 2*pi, doN+1);
doAngles(end) = []; % remove last point (= first point)
resDo = zeros(xPix, yPix, doN);
for i = 1:doN
    doX = round(centerX + doR * cos(doAngles(i)));
    doY = round(centerY + doR * sin(doAngles(i)));
    doRR = sqrt((X-doX).^2 + (Y-doY).^2);  
    resDo(:,:,i) = I0 * (doRR.^2 / w0^2) .* exp(-2 * doRR.^2 / w0^2);
end

%% calc resulting intensity profiles
resProt = protImg.*resDo; 
resSumInt = squeeze(sum(resProt,[1 2]));
resDoSat = resDo;
maxDo = max(resDo, [], 'all');
thDo = 0.50; % prozent of saturation
resDoSat(resDoSat>thDo*maxDo) = thDo*maxDo;
resProtSat = protImg.*resDoSat;
resSumIntSat = squeeze(sum(resProtSat,[1 2]));
%resProtBlurr = imfilter(resProt, fspecial('gaussian', [50 50], 10), 'replicate');
%% plot stuff
close all
figure
plot(resSumInt)
figure
plot(resSumIntSat)
figure
sliceViewer(resDoSat, 'DisplayRange', [min(resDoSat,[],'all') max(resDoSat,[],'all')], 'SliceNumber',1)
figure
sliceViewer(resProtSat, 'DisplayRange', [min(resProtSat,[],'all') max(resProtSat,[],'all')], 'SliceNumber',1)
% figure
% montage(resProt, 'DisplayRange',[])
% figure
% montage(resProtSat, 'DisplayRange',[])
% figure
% montage(resDo, 'DisplayRange',[])
% figure
% montage(resDoSat, 'DisplayRange',[])