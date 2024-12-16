function [Cmax,scale] = draw_max(name,serie)
% from a czi figure it is making a composite picture with the maximal
% values
% name: the name and path of the czi image (eg: 'Experiment-1268.czi')
% serie: the number of the series in the czi (eg: 1)
% Cmax: the composit matrix from the maximal values
% scale: scale in um


disp('Read image...')
% Construct an empty Bio-Formats reader
r = bfGetReader();
% Decorate the reader with the Memoizer wrapper
r = loci.formats.Memoizer(r);
% Initialize the reader with an input file
% If the call is longer than a minimal time, the initialized reader will
% be cached in a file under the same directory as the initial file
% name .large_file.bfmemo
r.setId(name);

% Perform work using the reader

% set the serie, C-type numbering, (0:end-1)
r.setSeries(serie - 1);

C = uint16(zeros(r.getSizeY(),r.getSizeX(),3));
Cmax = C;

for plane = 1:2:r.getImageCount()
    C(:,:,1) = bfGetPlane(r, plane);
    C(:,:,2) = bfGetPlane(r, plane+1);
    Cmax(C>Cmax) = C(C>Cmax);
end

scale = double(r.getMetadataStore().getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER));

% Close the reader
r.close()

end

