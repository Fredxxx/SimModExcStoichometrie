function protImg = genProtImg(xPix, yPix, pixS, r, N)
    protImg = zeros(xPix, yPix);
    % angles for N points
    angles = linspace(0, 2*pi, N+1);
    angles(end) = []; % remove last point (= first point)
    % set points on circle
    for i = 1:N
        x = round(xPix/2 + r/pixS * cos(angles(i)));
        y = round(yPix/2 + r/pixS * sin(angles(i)));
        if x > 0 && x <= xPix && y > 0 && y <= yPix
            protImg(y, x) = 1;
        end
    end
end