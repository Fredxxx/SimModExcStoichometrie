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
doN = 36; % how many points for donut circle movement
w0 = lambda/(pi*NA)/pixS; % Gauss width
thDo = 0.05; % prozent of saturation for superRes

%% generate images 
% protein structure
protImg = genProtImg(xPix, yPix, pixS, r, N);

minR = 5;
maxR = 50;
nR = 46;
R = linspace(minR,maxR,nR);

for j = 1:nR % j = 1
    %R(j) = 20;
    %generate donuts with different positions
    resDo = genDonImgs(xPix, yPix, w0, R(j), pixS, doN, I0);
    
    %% calc resulting intensity profiles
    resProt = protImg.*resDo; 
    resSumInt(j,:) = squeeze(sum(resProt,[1 2]));
    resMod(j) = (max(resSumInt(j,:))-min(resSumInt(j,:)))/(max(resSumInt(j,:))+min(resSumInt(j,:)));
    %resModSim(j) = (max(resSumInt(j,:))-min(resSumInt(j,:)));
    
    resDoSat = resDo;
    maxDo = max(resDo, [], 'all');
    resDoSat(resDoSat>thDo*maxDo) = thDo*maxDo;
    resProtSat = protImg.*resDoSat;
    resSumIntSat(j,:) = squeeze(sum(resProtSat,[1 2]));
    resModSat(j) = (max(resSumIntSat(j,:))-min(resSumIntSat(j,:)))/(max(resSumIntSat(j,:))+min(resSumIntSat(j,:)));
    %resModSatSim(j) = (max(resSumIntSat(j,:))-min(resSumIntSat(j,:)));
end

%% plot stuff
% figure
% plot(resSumInt)
% figure
% plot(resSumIntSat)
% figure
% sliceViewer(resDoSat, 'DisplayRange', [min(resDoSat,[],'all') max(resDoSat,[],'all')], 'SliceNumber',1)
% figure
% sliceViewer(resProtSat, 'DisplayRange', [min(resProtSat,[],'all') max(resProtSat,[],'all')], 'SliceNumber',1)
%% plot rdependency
figure
plot(R,resModSat,R,resMod)

%% plot all traces of summed intensities r-dependency
% figure
% for j = 1:nR
%     %Legend{j} = doRR(j)
%     plot(resSumIntSat(j,:),'DisplayName',sprintf('%.0f',R(j))); hold on
% end
% title('superResolved')
% legend
% figure
% for j = 1:nR
%     plot(resSumInt(j,:), 'DisplayName',sprintf('%.0f',R(j))); hold on
% end
% title('dif limited')
%% montage
% figure
% montage(resProt, 'DisplayRange',[])
% figure
% montage(resProtSat, 'DisplayRange',[])
% figure
% montage(resDo, 'DisplayRange',[])
% figure
% montage(resDoSat, 'DisplayRange',[])