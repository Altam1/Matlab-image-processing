%////1//////////////////////////////////////////////////
clc
clear all
close all

%% First
% Read an image
my_image = imread('my_image3.png');

% Display the original image
figure; imshow(my_image); title('Original Image');

% Add Gaussian noise to the image
noisy_image = imnoise(my_image, 'gaussian', 0.05);

% Display the noisy image
figure; imshow(noisy_image); title('Noisy Image');

% Apply spatial domain smoothing filter
kernel = fspecial('average', [5 5]); % 5x5 mean filter
smooth_image = imfilter(noisy_image, kernel);
imshow(smooth_image);

% Apply different threshold smoothing
threshold1 = 10;
threshold2 = 50;
threshold3 = 100;

smooth_image1 = smooth_image .* uint8(smooth_image > threshold1);
smooth_image2 = smooth_image .* uint8(smooth_image > threshold2);
smooth_image3 = smooth_image .* uint8(smooth_image > threshold3);

% Display the smoothed images with different thresholds
figure;
subplot(2,2,1); imshow(smooth_image); title('Smooth Image');
subplot(2,2,2); imshow(smooth_image1); title('Smooth Image (Threshold 10)');
subplot(2,2,3); imshow(smooth_image2); title('Smooth Image (Threshold 50)');
subplot(2,2,4); imshow(smooth_image3); title('Smooth Image (Threshold 100)');
figure, imhist(my_image), title('Histogram of the Original Image')
%% Part 2

% Convert the input image to a binary image using a threshold value of 0.4
binary_image = im2bw(my_image, 0.4);

% Display the binary image
figure, imshow(binary_image), title('Binary Image');

% Resize the input image to a size of 200x200 and display
resized_image = imresize(my_image, [200 200]);
figure, imshow(resized_image), title('Resized Image');

% Adjust the intensity levels of the input image using the imadjust function
%P = imadjust(my_image,[0.3 0.7],[0 1]);
adjusted_image = imadjust(my_image,[.1 .1 0; .9 .9 1],[]); % Adjust the contrast of an RGB image by specifying contrast limits
figure, imshow(adjusted_image), title('Intensity Adjusted Image');
figure, imhist(adjusted_image), title('Histogram of the Intensity Adjusted Image');

%% Use Laplacian operator to enhance the image
%laplacian = [-1 0 1];
%laplacian = [1 0 0; 0 0 0; 0 0 0];
%laplacian = [0 -1 0; -1 5 -1; 0 -1 0]; % Laplacian operator
laplacian = [1 0 0; -1 1 0; 0 0 0];
sharpened_image = imfilter(my_image, laplacian);

% Display the sharpened image
figure; imshow(sharpened_image); title('Sharpened Image');

%% Use different operators to process the image
sobel_x = [-1 0 1; -2 0 2; -1 0 1]; % Sobel X operator
sobel_y = [-1 -2 -1; 0 0 0; 1 2 1]; % Sobel Y operator
prewitt_x = [-1 0 1; -1 0 1; -1 0 1]; % Prewitt X operator
prewitt_y = [-1 -1 -1; 0 0 0; 1 1 1]; % Prewitt Y operator

% Apply the operators to the input image using the imfilter function
sobel_x_image = imfilter(my_image, sobel_x);
sobel_y_image = imfilter(my_image, sobel_y);
prewitt_x_image = imfilter(my_image, prewitt_x);
prewitt_y_image = imfilter(my_image, prewitt_y);

% Display the resulting images
figure;
subplot(2,2,1); imshow(sobel_x_image); title('Sobel X Image');
subplot(2,2,2); imshow(sobel_y_image); title('Sobel Y Image');
subplot(2,2,3); imshow(prewitt_x_image); title('Prewitt X Image');
subplot(2,2,4); imshow(prewitt_y_image); title('Prewitt Y Image');
%% Part 3
% Applying different low-pass filters
gaussian_filter = fspecial('gaussian', [5 5], 2); % Gaussian filter
average_filter = fspecial('average', [5 5]); % Mean filter

gaussian_image = imfilter(my_image, gaussian_filter);
average_image = imfilter(my_image, average_filter);

% Applying different high-pass filters
laplacian_of_gaussian = fspecial('log', [9 9], 9); % LoG filter
sobel = [-1 0 1; -2 0 2; -1 0 1]; % Sobel operator

log_image = imfilter(my_image, laplacian_of_gaussian);
sobel_image = imfilter(my_image, sobel);
% Displaying the images with different filters
figure;
subplot(2,2,1); imshow(gaussian_image); title('Gaussian Filtered Image');
subplot(2,2,2); imshow(average_image); title('Average Filtered Image');
subplot(2,2,3); imshow(log_image); title('Laplacian of Gaussian Filtered Image');
subplot(2,2,4); imshow(sobel_image); title('Sobel Filtered Image');