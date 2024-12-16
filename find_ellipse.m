function [stats,rmax] = find_ellipse(A_in)
% find the parameters of the ellipse fitted on the colony
% rmax / pixel: equivalent circle radius
% stats: output parameters from regionprops

% the area where we have cells
A = A_in(:,:,1)==2 | A_in(:,:,2)==2;

% fill up holes
A = imfill(A,"holes");

% find ellipses
stats = regionprops(A, 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Centroid', 'Perimeter');

% the center should be around the center of the image
image_center = size(A)/2;
notgoodenough = sum(((reshape([stats.Centroid],2,length(stats))).'-image_center).^2,2)>(size(A,1)/10).^2;
stats(notgoodenough) = [];
% radius should be large enough
notgoodenough = [stats.MinorAxisLength]<size(A,1)/4;
stats(notgoodenough) = [];
% find the maximal perimeter
stats = stats(find([stats.Perimeter]==min([stats.Perimeter]),1,'first'));

a = stats.MinorAxisLength/2;
b = stats.MajorAxisLength/2;
rmax = (a+b)/2;

end

