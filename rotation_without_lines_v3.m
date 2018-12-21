clc,clear all,close all,
load('main_data_range1_box.mat')
dx = 0.2; dy=0.1;dz=0.2;%range and cross range resolution (cm)
k=dx/dy;                 %scaling factor
[n3 n4]=size(F3);
extX=500; extY=500;
A=zeros(extX,extX+n4);
B=zeros(n3,extY);
F3=cat(2,B,F3);
ReconImag=cat(1,F3,A);
ZeroMat = ones(size(ReconImag));
[n3 n4]=size(ReconImag);
% MaxDim = max([N3 N4]);
% MaxDim = max([n3*dx n4*dy]);
% % MaxD = max([dx dy]);
% % MaxDim = MaxDim/MaxD;
% IncFc = 3;
% ReconImag = zeros(IncFc*MaxDim);
% ReconImag(((IncFc/2)*MaxDim-N3/2+1):((IncFc/2)*MaxDim+N3/2),...
%     ((IncFc/2)*MaxDim-N4/2+1):((IncFc/2)*MaxDim+N4/2))=F3;
% [n3 n4]=size(ReconImag);

figure(1);
[s1 s2]=size(F3);
ys=(1:s1)*dy;
xs=(1:801)*dx;
[Xs Ys]=meshgrid(xs,ys);
surf(Xs,Ys,abs(F3(:,500:1300)));
pbaspect([1 s1*dy/(801*dx) 1])
shading interp
axis tight
view([0 90])

%% selection of vertices
% xd=input('enter each range point cm:');
% yd=input('enter each crossrange point mm:');
% sm1=input('give the slope of the first wall:');
% sm2=input('give the slope of the second wall:');

% %% finding the slope of the lines
l11=[(911)*dx 500*dy];
l12=[(902)*dx 338*dy];
l21=l12;
l22=[825*dx 189*dy];

sm1=round((l12(2)-l11(2))/(l12(1)-l11(1))); %line 1 slope
sm2=round((l22(2)-l21(2))/(l22(1)-l21(1))); %line 2 slope 
b1 = -sm1*l11(1)+l11(2);
b2 = -sm2*l21(1)+l21(2);
xa = (1:n3)*dx;
ya = (1:n4)*dy;
y1a = sm1*xa + b1; %equation for line 1
y2a = sm2*xa + b2; %equation for line 2
c1=zeros(1,n4); c2=zeros(1,n4);
for ii=1:n3
    
    for jj=1:n4

%find the difference between the scanning point and the line
sk1=abs(y1a(ii)-ya(jj));
sk2=abs(y2a(ii)-ya(jj));

    if sk1<(10*dy)
        c1(jj)=ii;
    end
    
    if sk2<(10*dy)
       c2(jj)=ii;
   end
%    if kk==1
%        break
%    end
    end
end
% Find the lines common point
for ii=1:n4
    
    df=abs(c1(ii)-c2(ii));    
   if df<2 && c1(ii)>0 && c2(ii)>0
       s=0;
       cx=c1(ii); % range coordinate
       cy=ii; % cross-range coordinate
   end
   
end
%% rotation with respect to line 1
td1=2*atan(sm1);
% td1=0.95*pi;

for ii=cy:(n3-extX)  %lines
%     for jj=((IncFc/2)*MaxDim-N4/2+1):(((IncFc/2)*MaxDim+N4/2)-1) %columns

    for jj=c1(ii):n4 %columns
      
y = 0; x= (jj-c1(ii))*dx;%temporary coordinates
%         if jj>c1(ii)
            %s=0
x_pr = cos(td1)*x - sin(td1)*y; y_pr = sin(td1)*x + cos(td1)*y ;% local rotated point
x_2pr = x_pr + c1(ii)*dx; y_2pr = y_pr +ii*dy; %shifted point
yin = round(x_2pr/dx); xin = round(y_2pr/(dy)); % convert to indices

if xin<1 || yin<1 || xin>n3 || yin>n4
 
  ReconImag(ii,jj)=0;
 
else
%   tic 
ReconImag(xin,yin) = ReconImag(xin,yin)+2*ReconImag(ii,jj);
% toc
ZeroMat(ii,jj) = 0;


%         end
end

    end  
   
end
%% rotation with respect to line 2
[sd1, sd2]=size(ReconImag);
td2=2*atan(sm2);
% td2 = pi;
for ii=1:(cy)  %lines
    for jj=c2(ii):n4 %columns
    
           
y = 0; x = (jj-c2(ii))*dx;%temporary coordinates
x_pr = cos(td2)*x - sin(td2)*y; y_pr = sin(td2)*x + cos(td2)*y ;% local rotated point
x_2pr = x_pr + c2(ii)*dx; y_2pr = y_pr +ii*dy; % shifted point
yin = round(x_2pr/dx); xin = round(y_2pr/(dy)); % convert to indices
if xin<1 || yin<1 || xin>n3 || yin>n4
  ReconImag(ii,jj) = 0;
else
ReconImag(xin,yin) = ReconImag(xin,yin)+2*ReconImag(ii,jj);
ZeroMat(ii,jj) = 0;
%         end
end
    end 
        

end
   ReconImag = ReconImag.*ZeroMat;  
figure(2)

% [Xs Ys]=meshgrid(xs,ys);
surf(Xs,Ys,abs(ReconImag));
pbaspect([1 s1*dy/(s2*dx) 1])
shading interp
view([0 90])
axis tight
