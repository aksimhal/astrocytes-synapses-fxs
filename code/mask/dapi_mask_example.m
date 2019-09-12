%% DAPI Mask Example
% Contact Anish with any questions  
% '2ss_DAPIAligned.tif' can be downloaded here: 
% https://figshare.com/articles/2ss-F000/8019605

% Location of the data 
data_location = '2ss_DAPIAligned.tif'; 

% Get information about the image
info = imfinfo(data_location);
num_images = numel(info);

% Create empty mask 
mask = zeros(info(1).Height, info(1).Width, numel(info)); 

% Iterate over each slice of the 
for n=1:num_images
    DAPI = imread(data_location, n); 
    probdapi = getProbMap(DAPI);
    bwimg = probdapi > 0.6;

    J = imclose(bwimg, strel('disk', 4, 4));
    J_morph = bwmorph(J,'bridge'); 
    J2 = imfill(J_morph, 'holes');
    J3 = bwareaopen(J2, 50); 

    mask(:, :, n) = J3;
end 

output_folder = 'output_masks';

if isfolder(output_folder) == false
    mkdir(output_folder)
end

output_fn = strcat(output_folder, filesep, 'DAPI-mask.tiff');

writeTIFFStacks(mask, output_fn);
disp(output_fn);

