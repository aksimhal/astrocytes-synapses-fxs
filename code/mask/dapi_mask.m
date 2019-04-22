%% Make DAPI Masks for 2ss, 3ss, 4ss, 5ss, 6ss, 7ss

% Location of the data 
base_dir = '/Users/anish/Documents/yi_mice/'; 

section_list = {'F000', 'F001', 'F002', 'F003'}; 
folder_str = 'ss_stacks';
fname = 'ss_DAPIAligned.tif';

for mouse_num=2:1:7
    for section_num=1:length(section_list)
        datalocation = strcat(base_dir, num2str(mouse_num), folder_str,...
            filesep, section_list{section_num}, filesep, num2str(mouse_num), fname); 
        disp(datalocation) 
        info = imfinfo(datalocation);
        
        num_images = numel(info);
        mask = zeros(info(1).Height, info(1).Width, numel(info)); 
 
        for n=1:num_images

            DAPI = imread(datalocation, n); 
            probdapi = getProbMap(DAPI, -1);
            bwimg = probdapi > 0.6;

            J = imclose(bwimg, strel('disk', 4, 4));
            J_morph = bwmorph(J,'bridge'); 
            J2 = imfill(J_morph, 'holes');
            J3 = bwareaopen(J2, 50); 

            mask(:, :, n) = J3;
        end 
        
        output_folder = strcat(base_dir, filesep, 'masks', filesep, 'DAPI',...
            filesep, num2str(mouse_num), folder_str, filesep,...
            section_list{section_num}, filesep);
        
        if isfolder(output_folder) == false
            mkdir(output_folder)
        end
        
        output_fn = strcat(base_dir, filesep, 'masks', filesep, 'DAPI',...
            filesep, num2str(mouse_num), folder_str, filesep,...
            section_list{section_num}, filesep, num2str(mouse_num), fname(1:2), '-DAPI-mask.tiff');
        
        writeTIFFStacks(mask, output_fn);
        disp(output_fn);
    end 
end


