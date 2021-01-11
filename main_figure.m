function main_figure(unmixed, M0, lambda, channel, filename)
    [fPath, fName, ~] = fileparts(filename);
    [~, maxdim] = max(M0);
    maxL = lambda(maxdim);
    map = struct;
    unmix = struct;
    colors = [];
    for i=1:size(M0,2)-1
        sRGB = spectrumRGB(maxL(i));
        if contains(channel{i}, 'Lipo')
            sRGB = [1,1,1];
        end
        unmix.(channel{i}) = squeeze(unmixed(:,:,i,:));
        figure(i)
        MIP = max(unmix.(channel{i})/max(unmix.(channel{i}),[],'all'), [], 3);
        map.(channel{i}) = [linspace(0,sRGB(1),256)', linspace(0,sRGB(2),256)', linspace(0,sRGB(3),256)'];
        imagesc(imadjust(MIP)); colormap(map.(channel{i}));
        MIPTiff = ind2rgb(im2uint8(MIP), map.(channel{i}));
        imgPath = fullfile(fPath,'Figures', [fName, '_MIP_', channel{i}, '.tif']);
        imwrite(MIPTiff, imgPath);
    end
end