function pd = fAt3Dj_v2(xc, yc, zc,xsize,ysize,zsize,delta_f,coeff)
% fAt3Dj_v2(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
% 这里xsize和ysize都是30，zsize是197，delta_f就是之前computeDelta3Dj_v2计算的

xc = max(xc,0);
xc = min(xc,xsize-1);

yc = max(yc,0);
yc = min(yc,ysize-1);

zc = max(zc,0);
zc = min(zc,zsize-1);

temp = coeff(xc+1,yc+1,zc+1,:);
pd=sum(delta_f.*(temp(:)));





