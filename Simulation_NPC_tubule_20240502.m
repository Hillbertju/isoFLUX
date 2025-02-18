% This code is to simulation isoFLUX raw data for reconstruction
% Author: Hongfei Zhu
% Note: if you want to use the code in your manuscript, please cite our paper

close all
clear
clc
addpath('shared\')


mode = 'saddle';
% mode = 'astig';
name = ['exp_',mode,'_'];
tmp = ['coeff_',mode,'_20240415.mat'];
load(tmp);

% data = readmatrix('simulation csv data\NPC0905_2.csv');  %==========================REVISE==============================================
% save_name = 'simulation results\NPC';   %============================REVISE==============================================
% x_real = data(:,1);  % x/y/z_real
% y_real = data(:,2);
% z_real = data(:,3);

data = readmatrix('simulation csv data\data_EPFL.csv');  %==========================REVISE==============================================
save_name = 'simulation results\EPFL';   %============================REVISE==============================================
x_real = data(:,3);  
y_real = data(:,4);
z_real = data(:,5);


num = 1;   % Repeat times for each localization 

x_pixelsize = 100;
z_pixelsize = 10;
off = 7;

[x, y, z , ~] = size(coeff);
zz =  (size(coeff,3))/2;

xsize = x+1;
ROI = xsize - 2*off;
r_ROI = (ROI-1)/2;

%% 0. system parameters
temp_xsize = 301; 
ysize = temp_xsize;
mid = (temp_xsize+1)/2;
wavelength = 640;
lambda = wavelength;
pixelsize = x_pixelsize;
k = 2*pi/wavelength;
NA = 1.49;  % 1.45
kmin = 2*pi/pixelsize/(temp_xsize-1);

sigma_est = lambda/4/NA/pixelsize;   
n = 1.518;   % sample medium 
swcycle = lambda/2/n;   
swcycle_pixel = swcycle/pixelsize;   
k_pixel = 2*pi/swcycle_pixel;   


Nx = (2500)/3;
Ny = (2500)/3;
mx = 0.95;
my = 0.95;
bglist = [5, 5, 5, 5, 5, 5];
sita_x = (0.1)/180*pi;
sita_y = (90.1)/180*pi;
phase_x = (60)/180*pi;
phase_y = (60)/180*pi;
sita_z = (40)/180*pi;


%% 1. import raw data
x_final1 = x_real./x_pixelsize;
y_final1 = y_real./x_pixelsize;
z_final1 = z_real./z_pixelsize;

x_final = repmat( x_final1,num,1)+(r_ROI+3);
y_final = repmat( y_final1,num,1)+(r_ROI+3);
z_final = repmat( z_final1,num,1);


%% 2. generate raw images
x_rough = round(x_final);   
y_rough = round(y_final);   
x_loc = x_final - x_rough;   
y_loc = y_final - y_rough;   
z_loc = z_final;

locs = size(x_loc,1);
image_stack  = zeros(ROI,ROI,locs,6);

tic
parfor index = 1:locs
            xc = x_loc(index);
            yc = y_loc(index);
            zc = zz+z_loc(index);

            delta_x = -1*xc;
            xstart = floor(delta_x);
            delta_x = delta_x - xstart;

            delta_y = -1*yc;
            ystart = floor(delta_y);
            delta_y = delta_y - ystart;
            
            delta_z = zc - floor(zc);
            [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2((delta_x),(delta_y),(delta_z));

            
            z_index = floor(zc);

            % 注意z/2这一层表示的是 zc属于0~1
            if z_index <=0
                z_index = 1;
            end
            if z_index >z
                z_index = z;
            end
            
            test1 = zeros(ROI, ROI);
            for ii = 0:ROI-1
                for jj = 0:ROI-1
                    temp = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,z,delta_f,coeff);
                    test1(ii+1,jj+1)=temp;
                end
            end

            test1 = test1./sum(test1(:));  % 光子数归一化

            xfinal = x_final(index);
            yfinal = y_final(index);
            zfinal = z_final(index)*z_pixelsize/x_pixelsize;
            N_x1 = Nx * ( 1+mx*sin( k_pixel*cos(sita_x)*cos(sita_z)*xfinal+k_pixel*sin(sita_x)*cos(sita_z)*yfinal + k_pixel*sin(sita_z)*zfinal + phase_x -2/3*pi ) ); 
            N_x2 = Nx * ( 1+mx*sin( k_pixel*cos(sita_x)*cos(sita_z)*xfinal+k_pixel*sin(sita_x)*cos(sita_z)*yfinal + k_pixel*sin(sita_z)*zfinal + phase_x ) );
            N_x3 = Nx * ( 1+mx*sin( k_pixel*cos(sita_x)*cos(sita_z)*xfinal+k_pixel*sin(sita_x)*cos(sita_z)*yfinal + k_pixel*sin(sita_z)*zfinal + phase_x +2/3*pi) );
            N_y1 = Ny * ( 1+my*sin( k_pixel*cos(sita_y)*cos(sita_z)*xfinal+k_pixel*sin(sita_y)*cos(sita_z)*yfinal + k_pixel*sin(sita_z)*zfinal + phase_y -2/3*pi ) ); 
            N_y2 = Ny * ( 1+my*sin( k_pixel*cos(sita_y)*cos(sita_z)*xfinal+k_pixel*sin(sita_y)*cos(sita_z)*yfinal + k_pixel*sin(sita_z)*zfinal + phase_y ) ); 
            N_y3 = Ny * ( 1+my*sin( k_pixel*cos(sita_y)*cos(sita_z)*xfinal+k_pixel*sin(sita_y)*cos(sita_z)*yfinal + k_pixel*sin(sita_z)*zfinal + phase_y +2/3*pi ) );
            Nlist = [N_x1, N_x2, N_x3, N_y1, N_y2, N_y3];
            for j = 1:1:6
                    image = Nlist(j)*test1;
        %             image = image(row_start:row_end,  col_start:col_end);
                    image = uint16(image)+bglist(j);   % add background noise                 
                    image = imnoise(image, 'poisson');   % add shot noise (poisson noise)
                    image_stack ( :, :, index , j) = double(image);
            end
    
end
toc
















