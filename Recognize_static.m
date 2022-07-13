close all; clc; clear all;
%% Load reference image, and compute surf features
ref_img = imread('Example_2.jpg');
[BW,maskedRGBImage] = createMask(ref_img);
ref_img_gray = rgb2gray(maskedRGBImage);
ref_pts = detectSURFFeatures(ref_img_gray);
[ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
figure; imshow(ref_img);
ref_pts_max=ref_pts.selectStrongest(50);
hold on; plot(ref_pts.selectStrongest(50));

%% Visual 25 SURF features
figure;
subplot(5,5,3); title('First 25 Features');
for i=1:25
    scale = ref_pts(i).Scale;
    image = imcrop(ref_img,[ref_pts(i).Location-10*scale 20*scale 20*scale]);
    subplot(5,5,i);
    imshow(image);
    hold on;
    rectangle('Position',[5*scale 5*scale 10*scale 10*scale],'Curvature',1,'EdgeColor','g');
end
%% Masking the video frame and computing surf features
image = imread('3.jpg');
[BW_2,maskedRGBImage_2] = createMask(image);
figure;
montage ({image, maskedRGBImage_2, BW_2}, 'Size', [1 3], 'BorderSize', [2 2], 'BackgroundColor', 'yellow')
I = rgb2gray(maskedRGBImage_2);
% Detect features
I_pts = detectSURFFeatures(I);
[I_features, I_validPts] = extractFeatures(I, I_pts);
figure;imshow(image);
hold on; plot(I_pts.selectStrongest(50));
%% Compare ref image to video frame
index_pairs = matchFeatures(ref_features, I_features);
ref_matched_pts = ref_validPts(index_pairs(:,1)).Location;
I_matched_pts = I_validPts(index_pairs(:,2)).Location;
figure, showMatchedFeatures(image, ref_img, I_matched_pts, ref_matched_pts, 'montage');
title('Showing all matches');

%% Convex hull reference image
CH = bwconvhull(BW);
s  = regionprops(CH, 'centroid');
centroids = cat(1, s.Centroid);
area_1=regionprops(CH, 'Area');
figure
imshow(CH)
hold on
plot(centroids(:,1), centroids(:,2), 'b*')
hold off
%% Convex hull video frame
CH_2 = bwconvhull(BW_2);
s_2  = regionprops(CH_2, 'centroid');
area_2=regionprops(CH_2, 'Area');
centroids_2 = cat(1, s_2.Centroid);
figure
imshow(CH_2)
hold on
plot(centroids_2(:,1), centroids_2(:,2), 'b*')
hold off
%% Skin effect plotting
u0          = 4*pi*1e-7;   % Permeability constant in [Vs/Am]
ur          = 1;           % relative Permeability of Material   e.g. Copper
rho         = 1.72e-8;     % Resistivity in [Ohm*m],             e.g. Copper
D_wire_m    = 1e-3;        % Wire-Diameter 0.001m
f           = [1e-3, 1e0, 16+2/3, 50 , 60,  1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]';       % frequency vector

% Calculate additonial parameter
R_m         = D_wire_m/2;  % Wire-Radius in [m²]
area        = R_m^2*pi;    % Circle_area in [m²]
sigma       = 1/rho;       % Conductivity, [S/m]

R_DC = 1/(sigma*area);     % DC-Resistance load per unit length in [Ohm/m]
omega  =  2*pi*f;                                  % angular frequency in [Hz]
delta  = 1./sqrt(omega.*sigma.*u0.*ur./2);         % skindepth in [m]

% calculate impedance of this wire;
z_math = Z_wire (omega, R_m,  sigma, u0, ur);

% R_DC_Ratio by an tubemodel-calculation (DC-Resistance in depence of skindepth and wire-radius)
R_dc_tube_ratio = tube_model_func(delta ,R_m ,sigma);


% *Collect the data for the table*
tab_header = {'frequency', 'skindepth', 'L lpul', 'R_DC lpul', 'R_DC-ratio ', 'R_DC-tubemodel_ratio'; ...
              'Hz',        'm',          'H/m',   'Ohm/m' ,    '[1]',          '[1]' };
tab2       = num2cell([ f, ...                                              % 1. column frequency
                        delta, ...                                          % 2. column skindepth
                        imag(z_math)./omega,  ...                           % 3. column impedance  (H/m)
                        real(z_math), ...                                   % 4. column resistance (Ohm/m)
                        real(z_math)/R_DC,...                               % 5. column resistance-Ratio [1]
                        R_dc_tube_ratio   ]);                               % 6. column resistance tube_model ratio [1]
tab       = [tab_header; cellfun(@num2str,(tab2), 'UniformOutput',0)];
% lpul: load per unit length
% PART II - Plot currentdensity on wire for 50 Hz
f_circle_Hz = 50;     % frequency in Hz
current_A   = 1;      % current in wire in Amps
dr          = 200;    % number of steps (dicretization

% create a grid for that circle
[X,Y] = meshgrid(-R_m: R_m/(dr-1) :R_m , -R_m: R_m/(dr-1) :R_m  );

% Calculate radius for each Point:
r = sqrt(X.^2+Y.^2);

% calculation of currentdensity:
A      =  sqrt(-1j*2*pi*f_circle_Hz*sigma*u0*ur);       % Nomalized argument for Bessel function
[J0]   =  besselj ( 0 , A.* r );                        % first bessel-function of zeroth order
[J1]   =  besselj ( 1 , A.*R_m );                       % first bessel-function of first order
J_vec  =  A.*current_A.*1./(2.*pi.*R_m).*J0./J1;

% Recalc only the wire with radius
J_vec((X.^2+Y.^2)>=R_m^2)   = NaN;
X(isnan(J_vec))             = NaN;
Y(isnan(J_vec))             = NaN;

% Plot results
figure;

set(gcf, 'color', 'none');    
set(gca, 'color', 'none');
% Plot the Current density
Surface = surf(X,Y,real(J_vec),'LineStyle','none');
set(gcf, 'color', 'none');    
set(gca, 'color', 'none');
% colorbar;
axis square
view([0 0 90]);
axis off;

set(gcf, 'color', 'none');    
set(gca, 'color', 'none');
saveas(Surface, 'Surface', 'png');
surf_img = imread('Surface.png');
[BW_Surf, Masked_Surf]=createMaskSurf(surf_img);
CH_surf = bwconvhull(BW_Surf);
area_surf=regionprops(CH_surf, 'Area');
% title (['Currentdensity in wire (R=', num2str(R_m), ' m, I=',num2str(current_A),' A, f=', num2str(f_circle_Hz), ' Hz)'  ]);

length_m    = 1000;  % lenghth of wire;

E_Vpm       = current_A .* A./(2 .* pi.* R_m.* sigma)  .* besselj ( 0 , A.* R_m )./ besselj ( 1 , A.*R_m )   % electrical field_strength in Volt per m
deltaU_V    = E_Vpm.*length_m                              % complex voltage drop in Volt
absU_V      = abs(deltaU_V)                                % absolute voltage drop in Volt
R_Ohm       = real(deltaU_V/current_A)                     % Resistance in Ohm
Li_H        = imag(deltaU_V/(2*pi*f_circle_Hz*current_A))  % inductance in Henry
set(gcf, 'color', 'none');    
set(gca, 'color', 'none');

%% Superimposing the skin effect chart and ref image
clearvars -except s s_2 area_1 area_2 area_surf
s_mat=cell2mat(struct2cell(s));
s_2_mat=cell2mat(struct2cell(s_2));

[RGB, MAP] = imread('Example_2.jpg');
I=imread('Surface.png');
[X_bg,cmap] = rgb2ind(RGB,256);
bg = ind2rgb(X_bg,cmap);
[RGB_surf, MAP_surf] = imread('Surface.png'); 
[X_surf,cmap_surf] = rgb2ind(RGB_surf,256);
im = ind2rgb(X_surf,cmap_surf);

scale=sqrt(area_1.Area/area_surf.Area)
t = 0:.1:2*pi;
x = round(266*cos(t)+454);
y = round(266*sin(t)+317);
imCircleMask = roipoly(im,x,y);
imCircleMask_resized = imresize(imCircleMask,scale);
im_resized = imresize(im,scale);

image(bg)
hold on;

[rows, columns]=size(im_resized);
[rows_res columns_res]=size(im_resized);
image(im_resized, 'alphadata',imCircleMask_resized,'XData',s_mat(1,1)-454*scale , 'YData', s_mat(1,2)-317*scale);
axis off
hold off
