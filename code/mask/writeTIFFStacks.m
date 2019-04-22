function writeTIFFStacks(data, path)


imwrite(data(:, :, 1), path, 'compression', 'none')

for n=2:size(data, 3)
    
    imwrite(data(:, :, n), path, 'WriteMode', 'append', 'compression', 'none')
    
end
