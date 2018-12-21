clc,clear all,close all,
load('..\Data\before_rotation_2.mat')
dx = 0.2; dy=0.05;%range and cross range resolution (cm)
                 %scaling factor
%F3 = F3 (end:-1:1,:); % Mirror the image vertically
[n3 n4]=size(F3_norm); % F3 is the raw reconstructed image
extX=5000; extY=1; % extend both dimensions of F3 to make sure the corrected image fits properly after rotating the pixels
A=(10^-100)*ones(extX,extY+n4);
B=(10^-100)*ones(n3,extY);
F33=cat(2,B,F3_norm); % Extend the range
ReconImag_i=cat(1,F33,A); % Extend the cross-range
ZeroMat = ones(size(ReconImag_i));
[n3 n4]=size(ReconImag_i);

%ReconImag=ReconImag_i(end:-1:1,:);
ReconImag=ReconImag_i;

[s1 s2]=size(ReconImag);
n3=s1;
n4=s2
ys=(1:s1)*dy;
xs=(1:s2)*dx;
[Xs Ys]=meshgrid(xs,ys);
% surf(Xs,Ys,abs(ReconImag(end:-1:1,500:1300)));
figure(1);
surf(Xs,Ys,20*(abs(ReconImag)));
pbaspect([1 s1*dy/(s2*dx) 1])
shading interp
axis tight
view([0 90])
fprintf('please check the graph to input the extreme points on the line:\n')
fprintf('now enter point 1:\n');
l11=input('enter point 1 array:')
l12=input('\nenter point 2 array:')
%% selection of vertices

% %% finding the slope of the lines

sm1=(l12(2)-l11(2))/(l12(1)-l11(1)); %line 1 slope
b1 = -sm1*l11(1)+l11(2);
xa = (1:n3)*dx;  % range coordinates with respect to indices
ya = (1:n4)*dy; % dy cross-range coordinates with respect to indices
y1a = sm1*xa + b1; %equation for line 1

c1=zeros(1,n3-extX); % finds the indices on the raw image where the mirror exists
for ii=1:n4  % indices for range
    
    for jj=1:(n3-extX) %indices for cross-range  

%find the difference between the scanning point and the line
sk1=abs(y1a(jj)-ya(ii)); %%%%%% 
%sk2=abs(y1a(ii)-ya(jj));
sk11(ii,jj)=sk1;
    if sk1<(5*dy) %dy
        c1(ii)=jj; %%%%%  
    end
    

    end
end


%% rotation with respect to line 1
td1=2*atan(sm1);
%td1=degtorad(td1);
checkpoint=1;
for ii=1:(n3-extX)  % scan the lines
%for ii=n3-extX:n3  
    for jj=(c1(ii)):n4 % scan the columns 
    
y = 0; x= (jj-c1(ii))*dx; %temporary coordinates  ############ also will not go all possible values of c1 why y=0 ?

x_pr = cos(td1)*x - sin(td1)*y; y_pr = sin(td1)*x + cos(td1)*y ;% local rotated point
x_2pr = x_pr + c1(ii)*dx; y_2pr = y_pr +ii*dy; %shifted point
yin = abs(round(x_2pr/dx)); xin = abs(round(y_2pr/(dy))); % convert to indices

if  xin>n3 || yin>n4 || xin<1 || yin<1

  ReconImag(ii,jj)=10^-100;
 
else
%   tic 
% ReconImag(xin,yin) = ReconImag(xin,yin)+7*ReconImag(ii,jj);
ReconImag(xin,yin) = ReconImag(xin,yin)+ReconImag(ii,jj);
% toc
ZeroMat(ii,jj) = (10^-100);

end

    end  
end

%    ReconImag = ReconImag.*ZeroMat;
 ReconImag = ReconImag.*ZeroMat;


ReconImagabs=20*log10(abs(ReconImag));
figure(2);
[s1 s2]=size(ReconImag);
ys=(1:s1)*dy;
% xs=(1:801)*dx;
xs=(1:s2)*dx;
[Xs Ys]=meshgrid(xs,ys);
% surf(Xs,Ys,abs(ReconImag(end:-1:1,500:1300)));
surf(Xs,Ys,abs(ReconImag(end:-1:1,:)));

pbaspect([1 s1*dy/(s2*dx) 1])
shading interp
axis tight
view([0 90])

figure(3)
[s11 s22]=size(ReconImagabs);
ys=(1:s11)*dy;
% xs=(1:801)*dx;
xs=(1:s22)*dx;
[Xs Ys]=meshgrid(xs,ys);
surf(Xs,Ys,ReconImagabs);
surf(Xs,Ys,ReconImagabs(end:-1:1,:));

% pbaspect([1 s1*dy/(801*dx) 1])
 pbaspect([1 s11*dy/(s22*dx) 1])
shading interp
axis tight
view([0 90])


