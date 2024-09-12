function resDo = genDonImgs(xPix, yPix, w0, r, pixS, doN, I0)
    [X, Y] = meshgrid(1:xPix, 1:yPix); % meshgrid
    doR = r/pixS;
    doAngles = linspace(0, 2*pi, doN+1);
    doAngles(end) = []; % remove last point (= first point)
    resDo = zeros(xPix, yPix, doN);
    for i = 1:doN
        doX = round(xPix/2 + doR * cos(doAngles(i)));
        doY = round(yPix/2 + doR * sin(doAngles(i)));
        doRRR = sqrt((X-doX).^2 + (Y-doY).^2);  
        resDo(:,:,i) = I0 * (doRRR.^2 / w0^2) .* exp(-2 * doRRR.^2 / w0^2);
    end
end