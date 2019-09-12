function data = getProbMap(data)
% Create Foreground Probablities
% data = getProbMap(data)

data = double(data);

for z=1:size(data, 3)
    
    img = data(:, :, z);
    imgvec = img(:);
    
    % Calculate foreground probabilities
    probImg = normcdf(double(img), mean(imgvec), std(double(imgvec)));
    
    data(:, :, z) = probImg;
    
end

end
