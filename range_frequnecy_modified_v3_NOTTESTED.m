clc;
clear all;
cd 'C:\Users\Public\Documents\Labview\SAi_test\5_23_2018\case_1_1_5000_1target_nlos';
s=load('point_data.mat');
s_u=getfield(s,'sort_data');
% x=input('please give the calibration value:');
% s_u=s_u-x;
freq_p=load('frequency_points.mat');
freq_q=getfield(freq_p,'fr');
fno=length(freq_q);
 B=zeros(250,fno);
 s_u=cat(1,B,s_u,B);
[m1 m2]=size(s_u);
for f=1:m2
   S_ku(:,f)=fftshift(fft(fftshift(s_u(:,f))));
 

end
w=2*pi*freq_q;
k=w/(3*10^8);
du_x=1*10^-3;
dku_x=2*pi/(m1*du_x);
ku_x=(dku_x)*(-m1/2:m1/2-1);
for g=1:m2
    for h=1:m1
        ku(h,g)=4*k(g)^2-ku_x(h)^2;
%         s_c(h,g)=exp(1i*k(g)*((m1/2)*du_x));
    end
end
ku1=sqrt(ku.*(ku>0));
rn=input('Range in cm:');
kk=1;
dz = 0.2; % range step in cm
z_start=1; % starting range value in cm
for x=z_start:dz:rn
S_0=exp(1i.*ku1*x*10^-2);
F=S_ku.*(S_0);
F2=zeros(m1,1);
sum=zeros(m1,1);
for p=1:m2
    sum=F(:,p);
    F2(:,1)=sum+F2(:,1);
end
F1=fftshift(ifftn(fftshift(F2),[m1 1]));
F3(:,kk)=F1;
kk=kk+1;
end
%  x=linspace(0,15,301);
%  y=(0:150*(3e8/(2*(-freq_q(1,1)+freq_q(2000,1))))*100:150);
%  [X Y]=meshgrid(x,y);
% F3=F3.';
figure(1);
[Zc,Xc] = meshgrid(z_start:dz:rn,(1:m1)*(du_x*100));
surf(Zc,Xc,abs(F3));
xlabel('Range cm');
ylabel('Cross Range mm');
title('reconstructed 3d');
% pbaspect([10 1 1]);
shading interp;
axis tight;
view([0 90]);
pbaspect([(rn-z_start)/(m1*du_x*100) 1 1]);
figure(2)
surf(Zc,Xc,20*log10(abs(F3)));
shading interp
pbaspect([(rn-z_start)/(m1*du_x*100) 1 1]);
S=abs(F3);
axis tight
view([0 90]);
% save('data.mat','S','Zc','Xc','du_x','m1','m2');
save('test.mat');
colorbar
caxis([-60 10])