function [sep_red_acf,sep_green_acf,ratio] = calc_polar(Cmax,dr,R_min,Rmax,stats,scale)
% calculate red ratio and separation width
% we are processing the image in the R_min, R_max range, with dr length
% divisions
% scale: um/pixel
% every (length-type) input parameter is in pixel!
% stats: properties of the ellipse
% ratio: calculated red ratio
% sep_red_acf: red separation from acf
% sep_green_acf: green separation from acf

nlag = 1000; % number of lags for acf calculation
corr_lim = 0.2;

% coordinates of the original picture
[Y,X] = ndgrid(1:size(Cmax,1),1:size(Cmax,2));
X = X-stats.Centroid(1);
Y = Y-stats.Centroid(2);

% transform according to the ellipse
a = stats.MajorAxisLength/2;
b = stats.MinorAxisLength/2;
r0 = (a+b)/2; % equivalent circle radius

% transform back the points
rotation = -stats.Orientation*pi/180;
rotate_inv = [cos(rotation), sin(rotation);...
  -sin(rotation), cos(rotation)];
stretch_inv = [r0/a,0;0,r0/b];

Minv = stretch_inv*rotate_inv;

X_trans = X;
Y_trans = Y;
for i = 1:size(X,1)
    for j = 1:size(X,2)
        trans = Minv*[X(i,j);Y(i,j)];
        X_trans(i,j) = trans(1);
        Y_trans(i,j) = trans(2);
    end
end
% X = X+stats.Centroid(1);
% Y = Y+stats.Centroid(2);

% transform to polar coordinate system
[theta,rho] = cart2pol(X_trans,Y_trans);
clear X_trans Y_trans X Y

% wee use the data only in the given radius range
r_range = rho(:)>=R_min & rho(:)<=Rmax;

rho = rho(r_range);
theta = theta(r_range);

% change theta for the arc length calculations with ellipsoids
% theta originally in [-pi,pi] => [0,2pi]
theta(theta<0) = theta(theta<0)+2*pi;

[rho,idx] = sort(rho(:));
theta = theta(idx);

% average red ratio
Rx = Cmax(:,:,1);
Gx = Cmax(:,:,2);
Rx = double(Rx(r_range));
Gx = double(Gx(r_range));
ratio = mean(Rx./(Rx+Gx));

% for acf calculation
Rx = Rx(idx);
Gx = Gx(idx);

r = (R_min:dr:Rmax).';
r_avr = (r(1:end-1)+r(2:end))/2;
nr = length(r)-1;

% memory allocation
sep_red_acf = zeros(nr,1);
sep_green_acf = zeros(nr,1);

angles = linspace(0,2*pi,1e3);
arc_length = zeros(length(angles),1);
for i = 1:length(angles)
    arc_length(i) = integral(@(fi) 1/r0*sqrt(a^2*sin(fi).^2+b^2*cos(fi).^2),0,angles(i));
end

for ri = 1:nr

    % average circle/arc radius
    rho_avr = r_avr(ri);

    start_pos = find(rho>=r(ri),1,'first');
    stop_pos = find(rho<r(ri+1),1,'last');
    positions = (start_pos:stop_pos).';

    % order them according to the polar angle
    [theta_s,idx] = sort(theta(positions));
    positions = positions(idx);


    % https://www.johndcook.com/blog/2022/11/02/elliptic-arc-length/
    % arc lenght computed from the origo without r
    s = rho_avr.*interp1(angles,arc_length,theta_s,'pchip');

    % acf calculation
    s_uni = linspace(0,s(end),length(s));

    % wee need a uniform distribution for acf
    R_uni = interp1(s,Rx(positions),s_uni);
    G_uni = interp1(s,Gx(positions),s_uni);

    % acf calculation
    acf_R = autocorr(R_uni,min(nlag,length(s_uni)-2));
    acf_G = autocorr(G_uni,min(nlag,length(s_uni)-2));
    tmp = s_uni(find(acf_R<=corr_lim,1,'first'));
    if isempty(tmp)
        sep_red_acf(ri) = nan;
    else
        sep_red_acf(ri) = tmp;
    end
    tmp = s_uni(find(acf_G<=corr_lim,1,'first'));
    if isempty(tmp)
        sep_green_acf(ri) = nan;
    else
        sep_green_acf(ri) = tmp;
    end

end

% pixel => um
sep_red_acf = scale*sep_red_acf;
sep_green_acf = scale*sep_green_acf;

end

