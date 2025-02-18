% This code is to calculate theoretical CRLB for ModLoc and conventional
% 3D-SMLM (6bg),using 3D patterns
% Author: Hongfei Zhu
% Note: if you want to use our code in your manuscript, please cite the paper!!!



close all;clear;clc;
addpath('helper function\')
% close all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase_list = [0 15 30 45 60 75 90 105]/180*pi;    
phase_list = [60]/180*pi;    
thetaz = (20)/180*pi;  %*************************************************************************************************************************************************************
% the angle between the side beam and transverse plane (focal plane)
x0 = 0.00001;
y0 = 0.00001;
xc = x0;
yc = y0;


% load('bead_astig_3dcal.mat')  % PSF_XYZ_3D 
% load('DH-Exp_3Dcorr.mat')
% coeff = cspline.coeff;
% psf = cspline_all.PSFsmooth;


mode = 'saddle';
% mode = 'astig';
name = ['exp_',mode,'_'];
tmp = ['coeff_',mode,'_20240415.mat'];
load(tmp);


% load('saddle_20240424.mat')
% name = 'exp_saddle_';
% psf  = PSF_XYZ(:,:,1:end-2);
% clear PSF_XYZ
% coeff = Spline3D_interpv1(psf);

off = 3;


%% 1.preprocess
% [xsize,ysize,zsize] = size(psf);
% [xsize,~,zsize] = size(psf);
[x,y,zzz,~] = size(coeff);
xsize = x+1;
zsize = zzz+1;
% psf = psf/max(psf(:));
zz = (zsize+1)/2;
ROI = xsize- 2*off;

%% 1.1 if needed, can flip the PSF
% for i = 1:zsize
%     temp = psf(:, :, i);
%     temp = rot90(temp,-1);
%     temp = fliplr(temp);
%     psf(:,:,i) = temp;
% end


%% 2. cubic spline and generate coeffcient
% coeff = Spline3D_interp(psf); 




%% Loop start
for global_phase = phase_list


%%1.parameters of the blinking event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_set = 3000;    %**************************************************************** total photon number
m = 0.95;    %**************************************************************** modulation depth
bg = 5;   %**************************************************************** background photon number for each sub image (6 sub images)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dd = 5;  
zrange = dd:dd:600;   %*****************************************************
zrange = [-1*fliplr(zrange) 0.00001 zrange];



x_pixelsize = 100;    %***************************************************** lateral pixelsize
z_pixelsize = 10;    %***************************************************** axial pixelsize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. system parameters 
lambda = 640;
n = 1.518;
k0 = 2*pi/lambda;
kn = n*k0;
swcycle = lambda/2/n;   % minimum period of the pattern (non-evanescent wave )(nm)
swcycle_pixel = swcycle/x_pixelsize;   % minimum period of the pattern (non-evanescent wave )(pixel)


mid = (zsize+1)/2;
zrange1 = zrange;
zrange = zrange/z_pixelsize;
zrange = zrange + mid;
num = size(zrange,2);
CRLBx_m = zeros(1,num);
CRLBy_m = zeros(1,num);
CRLBz_m = zeros(1,num);
index = 1;






for z = zrange
    % reset photon to normalize (align with N_set)
    z0 = zrange1(index);
    Nx = 1300;
    Ny = 1300;
    zc = (z);

    delta_x = -1*xc;
    xstart = floor(delta_x);
    delta_x = delta_x - xstart;

    delta_y = -1*yc;
    ystart = floor(delta_y);
    delta_y = delta_y - ystart;

    delta_z = zc - floor(zc);
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2((delta_x),(delta_y),(delta_z));

    z_index = floor(zc);
    if z_index <=0
        z_index = 1;
    end

    if z_index >zzz
        z_index = zzz;
    end

    testf = zeros(ROI, ROI);
    testx = zeros(ROI, ROI);
    testy = zeros(ROI, ROI);
    testz = zeros(ROI, ROI);
    for ii = 0:ROI-1
        for jj = 0:ROI-1
            tempf = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,zzz,delta_f,coeff);
            tempx = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,zzz,delta_dxf,coeff);
            tempy = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,zzz,delta_dyf,coeff);
            tempz = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,zzz,delta_dzf,coeff);
            testf(ii+1,jj+1)=tempf;
            testx(ii+1,jj+1)=tempx;
            testy(ii+1,jj+1)=tempy;
            testz(ii+1,jj+1)=tempz;
        end
    end

    total_Nx = 3*Nx*sum(testf(:));
    total_Ny = 3*Ny*sum(testf(:));
    total_N = total_Nx +total_Ny;
    %     fac = N_set/total_N;
    %     tempx = tempx*fac;
    %     tempy = tempy*fac;
    %     tempz = tempz*fac;
    Nx = Nx*(N_set/total_N);
    Ny = Ny*(N_set/total_N);
    total_Nx = 3*Nx*sum(testf(:));
    total_Ny = 3*Ny*sum(testf(:));
    total_N = total_Nx +total_Ny;


    % [Nx, Ny, mx, my, x0, y0, z0, bg1~bg6]  13 elements in total
    d_matrix = zeros(ROI, ROI, 13, 6);
    temp_M = zeros(13,13);
    for channel = 0:1:2
        d_matrix(:,:,1,channel+1) = ( 1+m*sin( kn*x0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) )*testf;
        d_matrix(:,:,3,channel+1) = Nx*testf*sin( kn*x0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase);
        d_matrix(:,:,5,channel+1) = -Nx*testx/x_pixelsize *( 1+m*sin( kn*x0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) )+Nx*testf*m*kn*cos(thetaz)*cos( kn*x0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase );
        d_matrix(:,:,6,channel+1) = -Nx*testy/x_pixelsize *( 1+m*sin( kn*x0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) );
        d_matrix(:,:,7,channel+1) = Nx*testz/z_pixelsize *( 1+m*sin( kn*x0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) )+Nx*testf*m*kn*(1-sin(thetaz))*cos( kn*x0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase );
    end
    for channel = 0:1:2
        d_matrix(:,:,2,channel+4) = ( 1+m*sin( kn*y0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) )*testf;
        d_matrix(:,:,4,channel+4) = Ny*testf*sin( kn*y0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase);
        d_matrix(:,:,5,channel+4) = -Ny*testx/x_pixelsize *( 1+m*sin( kn*y0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) );
        d_matrix(:,:,6,channel+4) = -Ny*testy/x_pixelsize *( 1+m*sin( kn*y0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) )+Ny*testf*m*kn*cos(thetaz)*cos( kn*y0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase );
        d_matrix(:,:,7,channel+4) = Ny*testz/z_pixelsize *( 1+m*sin( kn*y0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) )+Ny*testf*m*kn*(1-sin(thetaz))*cos( kn*y0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase );
    end
    
    d_matrix(:,:,8,1) = ones(ROI,ROI);
    d_matrix(:,:,9,2) = ones(ROI,ROI);
    d_matrix(:,:,10,3) = ones(ROI,ROI);
    d_matrix(:,:,11,4) = ones(ROI,ROI);
    d_matrix(:,:,12,5) = ones(ROI,ROI);
    d_matrix(:,:,13,6) = ones(ROI,ROI);
    for iii = 1:1:13
        for jjj = 1:1:iii
            for channel = 0:1:2
                temp_M(iii, jjj) = temp_M(iii, jjj) + sum(sum(d_matrix(:,:,iii,channel+1).*d_matrix(:,:,jjj,channel+1)./...
                    (bg+Nx*( 1+m*sin( kn*x0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) )*testf )));
                temp_M(iii, jjj) = temp_M(iii, jjj) + sum(sum(d_matrix(:,:,iii,channel+4).*d_matrix(:,:,jjj,channel+4)./...
                    (bg+Ny*( 1+m*sin( kn*y0*cos(thetaz)+kn*z0*(1-sin(thetaz)) +2*pi/3*(channel-1)+global_phase) )*testf )));
            end
        end
    end

    tmp_M = temp_M';
    for iii = 1:1:13
        tmp_M(iii,iii) = 0;
    end



    M = tmp_M + temp_M;
    M1 = inv(M);
    aCRLBx = (M1(5,5))^(0.5);
    aCRLBy = (M1(6,6))^(0.5);
    aCRLBz = (M1(7,7))^(0.5);
    CRLBx_m(index) = aCRLBx;
    CRLBy_m(index) = aCRLBy;
    CRLBz_m(index) = aCRLBz;
    index = index +1;
    %

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bg = 6*bg;
CRLBx_w = zeros(1,size(zrange,2));
CRLBy_w = zeros(1,size(zrange,2));
CRLBz_w = zeros(1,size(zrange,2));


index = 1;
for z = zrange
     % reset photon to normalize (align with N_set)
    NN = 1000;
    zc = (z);

    delta_x = -1*xc;
    xstart = floor(delta_x);
    delta_x = delta_x - xstart;

    delta_y = -1*yc;
    ystart = floor(delta_y);
    delta_y = delta_y - ystart;

    delta_z = zc - floor(zc);
    [delta_f,delta_dxf,delta_ddxf,delta_dyf,delta_ddyf,delta_dzf,delta_ddzf]=computeDelta3Dj_v2((delta_x),(delta_y),(delta_z));

    z_index = floor(zc);
    if z_index <=0
        z_index = 1;
    end

    if z_index >zzz
        z_index = zzz;
    end

    testf = zeros(ROI, ROI);
    testx = zeros(ROI, ROI);
    testy = zeros(ROI, ROI);
    testz = zeros(ROI, ROI);
    for ii = 0:ROI-1
        for jj = 0:ROI-1
            tempf = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,zzz,delta_f,coeff);
            tempx = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,zzz,delta_dxf,coeff);
            tempy = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,zzz,delta_dyf,coeff);
            tempz = fAt3Dj_v2(ii+xstart+off,jj+ystart+off,z_index,x,y,zzz,delta_dzf,coeff);
            testf(ii+1,jj+1)=tempf;
            testx(ii+1,jj+1)=tempx;
            testy(ii+1,jj+1)=tempy;
            testz(ii+1,jj+1)=tempz;
        end
    end

    total_N = NN*sum(testf(:));
    NN = NN*(N_set/total_N);


    % N, x, y, z, bg
    d_matrix = zeros(ROI, ROI, 5);
    temp_M = zeros(5,5);
    d_matrix(:, :, 1) = testf;
    d_matrix(:, :, 2) = -1*(NN)*testx/x_pixelsize;
    d_matrix(:, :, 3) = -1*(NN)*testy/x_pixelsize;
    d_matrix(:, :, 4) = 1*(NN)*testz/z_pixelsize;
    d_matrix(:, :, 5) = ones(ROI,ROI);
    for iii = 1:1:5
        for jjj = 1:1:iii
            temp_M(iii, jjj) = sum(sum(d_matrix(:,:,iii).*d_matrix(:,:,jjj)./(bg+1*(NN)*testf)));  
        end
    end
    tmp_M = temp_M';
    for iii = 1:1:5
        tmp_M(iii,iii) = 0;
    end

    M = tmp_M + temp_M;
    M1 = inv(M);
    aCRLBx = (M1(2,2))^(0.5);
    aCRLBy = (M1(3,3))^(0.5);
    aCRLBz = (M1(4,4))^(0.5);
    CRLBx_w(index) = aCRLBx; 
    CRLBy_w(index) = aCRLBy; 
    CRLBz_w(index) = aCRLBz; 
    index = index +1;

end
% 
% figure;
% plot(zrange1+dd, CRLBx_w);hold on 
% plot(zrange1+dd, CRLBy_w);hold on 
% legend('x 6bg','y 6bg')
% xlim([-750 750]);










imp_xy = mean(CRLBx_w)/mean(CRLBx_m)
imp_z = mean(CRLBz_w)/mean(CRLBz_m)
%
close all
zrange2 = zrange1 ;
% zrange1 = zrange1+40;
figure;
plot(zrange1, CRLBx_m);hold on 
plot(zrange1, CRLBy_m);hold on 
plot(zrange1, CRLBz_m);hold on 
plot(zrange1, CRLBx_w);hold on 
plot(zrange1, CRLBy_w);hold on 
plot(zrange1, CRLBz_w);hold on 
legend('x pre','y pre','z pre','x 6bg','y 6bg','z 6g')
% xlim([-750 750]);
% title('+-750 nm(precision)')
% xlim([-1000 1000]);
% title('+-1000 nm(precision)')
xlim([-600 600]);
title('+-600 nm(precision)')

figure;
plot(zrange1, CRLBx_w./CRLBx_m);hold on 
plot(zrange1, CRLBy_w./CRLBy_m);hold on 
plot(zrange1, CRLBz_w./CRLBz_m);hold on 
legend('x imp','y imp','z imp')
% xlim([-750 750]);
% title('+-750 nm(improvement)')
% xlim([-1000 1000]);
% title('+-1000 nm(precision)')
xlim([-600 600]);
title('+-600 nm(precision)')
x_imp_mean =  mean(CRLBx_w)./mean(CRLBx_m)
y_imp_mean =  mean(CRLBy_w)./mean(CRLBy_m)
z_imp_mean =  mean(CRLBz_w)./mean(CRLBz_m)

% figure;
% plot(zrange1, CRLBx_m);hold on 
% plot(zrange1, CRLBy_m);hold on 
% legend('x pre','y pre')
% xlim([-750 750]);
% title('+-750 nm(precision, m)')

% figure;
% plot(zrange1, CRLBx_m);hold on 
% plot(zrange1, CRLBy_m);hold on 
% plot(zrange1, CRLBx_w);hold on 
% plot(zrange1, CRLBy_w);hold on 
% legend('x pre','y pre','x 6bg','y 6bg')
% title('+-1000 nm(precision)')
% 
% figure;
% plot(zrange1, CRLBx_w./CRLBx_m);hold on 
% plot(zrange1, CRLBy_w./CRLBy_m);hold on 
% legend('x imp','y imp')
% title('+-1000 nm(improvement)')


% CRLBx_m = CRLBx_m + 0.103;  % for experiment astigmatism
% CRLBy_m = CRLBy_m + 0.103;  % for experiment astigmatism








% 
% temp = [num2str(thetaz/pi*180),'_',num2str(global_phase/pi*180),'degree_CRLB6bg.mat'];
% output = [name, temp];
% save(output, 'CRLBx_m',  'CRLBx_w','CRLBy_m',  'CRLBy_w','CRLBz_m',  'CRLBz_w', 'N_set', 'm', 'bg',  'zrange1','xc','yc','thetaz');



end







%%
% temp = [num2str(thetaz/pi*180),'_',num2str(global_phase/pi*180),'degree_CRLB6bg','_',num2str(N_set),'photons','.mat'];
% output = [name, temp];
% save(output, 'CRLBx_m',  'CRLBx_w','CRLBy_m',  'CRLBy_w','CRLBz_m',  'CRLBz_w', 'N_set', 'm', 'bg',  'zrange1','xc','yc','thetaz');


