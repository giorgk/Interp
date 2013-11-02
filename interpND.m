function v0 = interpND(p, x, y, z, v, no_data, method)
%interpND: interpolates 1-2-3D gridded data. Missing values are
%interpolated using inverse distance weight. Elevation in layers is also
%supported
if size(p,2) > 3
    error('interpolation greater than 3D is not supported')
end

if nargin ~= 7
    error('The number of input arguments must be 6')
end


if exist('OCTAVE_VERSION')
    v0 = interpNd_oct(p, x, y, z, v, no_data, method);
else
    v0 = interpNd_mat(p, x, y, z, v, no_data, method);
end

