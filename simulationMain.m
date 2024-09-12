%close all
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
doN = 12; % how many points for donut circle movement
w0 = lambda/(pi*NA)/pixS; % Gauss width
thDo = 0.50; % prozent of saturation for superRes

%% generate images 
% protein structure
protImg = genProtImg(xPix, yPix, pixS, r, N);

%generate donuts with different positions
resDo = genDonImgs(xPix, yPix, w0, 50, pixS, doN, I0);

%% calc resulting intensity profiles
resProt = protImg.*resDo; 
resSumInt = squeeze(sum(resProt,[1 2]));
resMod = (max(resSumInt)-min(resSumInt))/(max(resSumInt)+min(resSumInt));

resDoSat = resDo;
maxDo = max(resDo, [], 'all');
resDoSat(resDoSat>thDo*maxDo) = thDo*maxDo;
resProtSat = protImg.*resDoSat;
resSumIntSat = squeeze(sum(resProtSat,[1 2]));
resModSat = (max(resSumIntSat)-min(resSumIntSat))/(max(resSumIntSat)+min(resSumIntSat));
%% plot stuff
% figure
% plot(resSumInt)
figure
plot(resSumIntSat)
% figure
% sliceViewer(resDoSat, 'DisplayRange', [min(resDoSat,[],'all') max(resDoSat,[],'all')], 'SliceNumber',1)
% figure
% sliceViewer(resProtSat, 'DisplayRange', [min(resProtSat,[],'all') max(resProtSat,[],'all')], 'SliceNumber',1)
% % figure
% montage(resProt, 'DisplayRange',[])
% figure
% montage(resProtSat, 'DisplayRange',[])
% figure
% montage(resDo, 'DisplayRange',[])
% figure
% montage(resDoSat, 'DisplayRange',[])
%% plot
% figure
% for j = 1:max(size(doRR))
%     %Legend{j} = doRR(j)
%     plot(loopSumIntSat(j,:),'DisplayName',sprintf('%.0f',doRR(j))); hold on
% end
% title('superResolved')
% legend
% figure
% for j = 1:max(size(doRR))
%     plot(loopSumInt(j,:), 'DisplayName',sprintf('%.0f',doRR(j))); hold on
% end
% title('dif limited')