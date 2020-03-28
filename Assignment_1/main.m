%% ENSC 895 G200 Assignment 1
% Due: Friday 17 January 2020
% -----
% @author: Ashiv Hans Dhondea
% Course: ENSC 895 G200
% Institution: Simon Fraser University - Engineering Science
% email: hdhondea@sfu.ca
% ----
% This is the base script which calls functions which implement each
% question in this assignment.
% ----
% Version history:
% 13/01/20: Created
% 15/01/20: Edited: did question 2,

%% Clear memory, figures and clean up
clear;close all;clc;
%% Question 1
% Plot a 2D sinc function over a 2D cartesian grid where x ranges from 100 
% to 355 and y ranges from 0 to 255 and the sinc is centered at somewhere 
% not midway in the grid. First plot the grid, and then plot the function 
% on the grid (both as grayscale and as a surface plot).
fprintf('Question 1\n');
% Add path to function file.
addpath('Codes and Images\Question 1');

% Define the start and end of the x and y dimensions
xrange = [100,355];
yrange = [0,255];
% Define the centre of the sinc function
centre = [300,200];

% Function call to create the 2D sinc surface and plot the results
fprintf('Now generating 2D sinc surface\n');
fnCreate2DSinc( xrange,yrange,centre);
%% Question 2
% Use your smartphone to take a selfie. Load this image in Matlab. Display 
% the image, and display the red, green and blue channels separately.

fprintf('Question 2\n');
% Add path to function file.
addpath('Codes and Images\Question 2');
% Selfie input file. This can be changed to a different file name.
input_file = ('Codes and Images\Question 2\Images\IMG_20190227_101447.jpg');
% Read in image, convert to double and scale to (0,1)
I = imread(input_file);
J_double = double(I);
J = mat2gray(J_double);

% Find factor to resize image to less than 400 pixels wide.
% It is assumed that input image is larger than 400 x 400.
J_ratio = size(J)/400;
chosen_factor = 1./max(J_ratio);
% Image is resized using the nearest neighbour interpolation method.
J_resized = imresize(J,chosen_factor,'nearest');
% Confirm that resized image is smaller than 400x400.
fprintf('The resized image is %d by %d.\n',size(J_resized,1),size(J_resized,2));

% Function call for Question 2
[R,G,B] = fnRGBchannels( J_resized );
%% Question 3
% Convert to grayscale (average the three channels I = (R+G+B)/3) and display the grayscale image. 

fprintf('Question 3\n');
% Add path to function file.
addpath('Codes and Images\Question 3');

[ J_average ] = fnGrayscale( R,G,B );
%% Question 4
% Display the grayscale image as a height-map (surface plot).

fprintf('Question 4\n');
% Add path to function file.
addpath('Codes and Images\Question 4');

fnImageSurf( J_average )
