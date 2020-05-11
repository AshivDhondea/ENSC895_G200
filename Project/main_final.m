%% ENSC 895 G200 Project
% Segment lungs in coronal CT scans
% Note:
% I have commented out code which creates a number of figures. You can
% uncomment them to see that they produce the exact same results as in the
% output_images folder.
% If you want to test on a different folder of images, replace the input
% directory with the correct name.
%% Clear memory, figures and clean up
clear;close all;clc;
%% Results file naming
% Find the name of this script. Used to name all output created by this
% script.
script_name = mfilename;
%% Load utilities
addpath(genpath('Utilities'));
%% Load functions
addpath(genpath('Functions'));
%% Select folder for testing.
test_img_folder_str = 'ENSC478_2020_ProjectTestSet';
fprintf('Classifiying the CT scans in the test set \n')
fnCreateFileList(script_name,test_img_folder_str);

testset_str = strcat(strcat(script_name,'_'),test_img_folder_str);
files_table = readtable(strcat(testset_str,'.xlsx'));
%% Process the images
% array preallocation
test_set_total_nnz = zeros(size(files_table,1),1);
test_set_threshed2_nnz = zeros(size(files_table,1),1);
test_set_threshed2_nnz_ratio = zeros(size(files_table,1),1);
num_pix_touching_border = zeros(size(files_table,1),1);
diagnostic = false(size(files_table,1),1);

% Loop through all the test images
for index = 1:size(files_table,1)
    %% Read image file
    file_str = strcat(test_img_folder_str,'\');
    img_str = strcat(file_str,string(files_table{index,1}));
    %% Otsu thresholding
    % Convert to grayscale and threshold using Otsu's method.
    J_grayscale = fnReadImage(img_str);
    J_rgb = fnReadImage_rgb(img_str);
    otsu_level = graythresh(J_grayscale);
    otsu_bw = imbinarize(J_grayscale, otsu_level);
    % Display the Otsu thresholded image
    f = figure;
    title('Binary image from Otsu thresholding','Interpreter','latex');
    hold on;
    imshow(otsu_bw);
    hold off;
    figfile = strcat(script_name,strcat('_otsu_',string(files_table{index,1})));
    print(gcf,'-dpng',strcat('output_images\',figfile));
    %%
    % Get rid of stuff touching the border (the area around the patient's
    % body)
    binaryImage = imclearborder(~otsu_bw);
    % Extract only the two largest blobs. (because we want both lungs)
    binaryImage = bwareafilt(binaryImage, 2);
    % Fill holes in the blobs to make them solid.
    binaryImage = imfill(binaryImage, 'holes');

    % Mask all pixels outside the regions identified as lungs to 0. 
     J_grayscale_masked = J_grayscale;
     J_grayscale_masked(~binaryImage) = 0;
     
%     f = figure;
%     title('Masked test image','Interpreter','latex');
%     hold on;
%     imshow(im2uint8(mat2gray(J_grayscale_masked)));
%     hold off;
%     figfile = strcat(script_name,strcat('_masked_',string(files_table{index,1})));
%     print(gcf,'-dpng',strcat('output_images\',figfile));

    %% Morphological processing to segment the lungs
    % Use two linear structuring elements of length 3 and angle 0 and 90.
    struc_el_90 = strel('line', 3, 90);
    struc_el_0 = strel('line', 3, 0);
    % Dilate image
    edges_dilated = imdilate(binaryImage, [struc_el_90 struc_el_0]);
    % fill interior gaps and remove connected objects on the image border
    edges_dilated_filled = imfill(edges_dilated, 'holes');
    edges_dilated_clear = imclearborder(edges_dilated_filled, 4);
    % Use a diamond structuring element to erode the image to smoothen the
    % lungs
    struc_el_dia = strel('diamond',1); % 1
    img_smoothened_lungs = imerode(edges_dilated_clear,struc_el_dia);
    segmented_lungs = imerode(img_smoothened_lungs,struc_el_dia);

    % Outline the perimeter of the lungs in the original image
    ggo_patches = bwperim(segmented_lungs);
    % draw lines to indicate the border of the lungs
    thickOutlines = imdilate(ggo_patches, true(3));
    % mask all pixels outside the lungs to zero.
    segmented_lungs_outlined = J_grayscale_masked;
    segmented_lungs_outlined(thickOutlines) = 255;
    
    test_set_total_nnz(index) = nnz(J_grayscale_masked);
    
    % Label the lungs
    Labels = bwlabel(segmented_lungs);
    info = regionprops(Labels,'Centroid'); % get the centroids to label the
    % lungs in the plots
    %% Plot the segmented lungs: using the perimeters as boundary
    % Plot the segmented and labeled lungs
    f = figure;
    imshow(segmented_lungs_outlined);
    hold on
    for k = 1 : length(info)
        c = info(k).Centroid;    
        text(c(1), c(2), sprintf('%d', k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','color','r','BackgroundColor','w');
    end
    title('Segmented lungs','Interpreter','latex');
    hold off;
    figfile = strcat(script_name,strcat('_seg_',string(files_table{index,1})));
    print(gcf,'-dpng',strcat('output_images\',figfile));
    
    %% Plot the segmented lungs: using bounding boxes
    Labels = bwlabel(segmented_lungs_outlined);
    info = regionprops(Labels,'Boundingbox');
%     f = figure;
%     imshow(segmented_lungs_outlined);
%     hold on
     
    lung_coords = zeros(4,2);
    for k = 1 : length(info)
        BB = info(k).BoundingBox;
%         rectangle('Position', [BB(1),BB(2),BB(3),BB(4)],'EdgeColor','r','LineWidth',2) ;    
%         text((BB(1)+round(0.5*BB(3))), BB(2)+BB(4), sprintf('%d', k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','color','r','BackgroundColor','w');
        
        lung_coords(:,k) = [BB(1),BB(2),BB(3),BB(4)];
    end
%     title('Segmented lungs: labeled','Interpreter','latex');

%     overall_bbox = [min(lung_coords(1,:)),min(lung_coords(2,:)),max(lung_coords(1,1)+lung_coords(3,1),lung_coords(1,2)+lung_coords(3,2)),max(lung_coords(2,1)+lung_coords(4,1),lung_coords(2,2)+lung_coords(4,2))];
%     rectangle('Position', [overall_bbox(1),overall_bbox(2),overall_bbox(3)-overall_bbox(1),overall_bbox(4)-overall_bbox(2)],'EdgeColor','g','LineWidth',2) ; 
%     
%     hold off;
%     figfile = strcat(script_name,strcat('_bbox_',string(files_table{index,1})));
%     print(gcf,'-dpng',strcat('output_images\',figfile));
    
    % Use a median filter of kernel size 7
    median_filter_size = 7;
    J_grayscale_masked_medfilt = medfilt2(J_grayscale_masked,[median_filter_size median_filter_size]); % with zero padding
    
    f = figure;
    imshow(im2uint8(mat2gray(J_grayscale_masked_medfilt)));
    hold on;
    title('Segmented lungs: filtered','Interpreter','latex');
    overall_bbox = [min(lung_coords(1,:)),min(lung_coords(2,:)),max(lung_coords(1,1)+lung_coords(3,1),lung_coords(1,2)+lung_coords(3,2)),max(lung_coords(2,1)+lung_coords(4,1),lung_coords(2,2)+lung_coords(4,2))];
    rectangle('Position', [overall_bbox(1),overall_bbox(2),overall_bbox(3)-overall_bbox(1),overall_bbox(4)-overall_bbox(2)],'EdgeColor','g','LineWidth',2) ; 
    hold off;
    figfile = strcat(script_name,strcat('_bbox_medfilt_',string(files_table{index,1})));
    print(gcf,'-dpng',strcat('output_images\',figfile));
    %% Multiple level thresholding
    numlevels = 5; % 5 levels of thresholding
    thresh = multithresh(J_grayscale_masked_medfilt,numlevels);
    seg_I = imquantize(J_grayscale_masked_medfilt,thresh);
        
    RGB = label2rgb(seg_I); % label according to threshold	 
    figure;
    imshow(RGB)
    colorbar('Ticks',thresh,'TickLabels',1:numlevels)
    axis off
    title(sprintf('Multiple level thresholding with %d levels', numlevels) ,'Interpreter','latex');
    
    figfile = strcat(script_name,strcat('_multithresh_',string(files_table{index,1})));
    print(gcf,'-dpng',strcat('output_images\',figfile));
    %%  Rejecting pixels belonging to classes which are not relevant to us
    below = 3; above = 5;
    threshold_chosen = thresh(below);
    
    % filter the image
    J_grayscale_masked_threshed = J_grayscale_masked_medfilt;
    J_grayscale_masked_threshed(J_grayscale_masked_medfilt < threshold_chosen) = 0;
    threshold_chosen_above = thresh(above);
    J_grayscale_masked_threshed(J_grayscale_masked_medfilt > threshold_chosen_above) = 0;
    J_grayscale_masked_threshed_rescaled = mat2gray(J_grayscale_masked_threshed);
    
%     figure;
%     imshow(im2uint8(J_grayscale_masked_threshed_rescaled));
%     title(sprintf('Multiple level thresholding: chosen levels [%d, %d]', below,above) ,'Interpreter','latex');
%     figfile = strcat(script_name,strcat('_multithresh_35_',string(files_table{index,1})));
%     print(gcf,'-dpng',strcat('output_images\',figfile));

    %%  Otsu thresholding to reject pixels belong to diffuse white patches
    otsu_level = graythresh(J_grayscale_masked_threshed);
    otsu_bw = imbinarize(J_grayscale_masked_threshed, otsu_level);

    % Extract only the largest blobs.
    big_blobs = bwareafilt(otsu_bw, 3);
    % Fill holes in the blobs to make them solid.
    big_blobs = imfill(big_blobs, 'holes');
    J_grayscale_masked_threshed(~big_blobs) = 0;
    
    %% Count the number of pixels touching the border
    num_pix_touching_border(index) = nnz(big_blobs & thickOutlines);
    
    if num_pix_touching_border(index) == 0
        fprintf('Patient %d tests negative for GGO.\n',index);
        diagnostic(index) = false;
        title_str = strcat(strcat('Patient',num2str(index)),'tests negative for GGO');
    else
        fprintf('Patient %d tests positive for GGO.\n',index);
        diagnostic(index) = true;
        
        % Count number of non-zero pixels in the image
        test_set_threshed2_nnz(index) = nnz(J_grayscale_masked_threshed);
        test_set_threshed2_nnz_ratio(index) = test_set_threshed2_nnz(index)/test_set_total_nnz(index); % this indicates the severity of the infection
        fprintf('The infection is %.2f %% \n' ,100*test_set_threshed2_nnz_ratio(index));
        title_str = strcat('Patient',{' '},num2str(index),' tests positive for GGO');
        
        %% Indicate where the GGOs occur.
        % Outline the perimeter of the lungs in the original image
        ggo_patches = bwperim(big_blobs);
        segmented_ggo_outlined = J_rgb;

        [r,c] = find(ggo_patches);
        for idx = 1:size(r,1)
            segmented_ggo_outlined(r(idx),c(idx),:) = [255,0,0];
        end
        f = figure;
        imshow(im2uint8(segmented_ggo_outlined));
        title(title_str,'Interpreter','latex');
        figfile = strcat(script_name,strcat('_labeled_',string(files_table{index,1})));
        print(gcf,'-dpng',strcat('output_images\',figfile));
    end
    %% show the lungs' outline
%     J_grayscale_masked_threshed(thickOutlines) = 255;
%     
%     f = figure;
%     imshow(im2uint8(J_grayscale_masked_threshed));
%     title(title_str,'Interpreter','latex');
%     figfile = strcat(script_name,strcat('_blobs_',string(files_table{index,1})));
%     print(gcf,'-dpng',strcat('output_images\',figfile));
%     
end

% close all;
