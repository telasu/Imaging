clc,clear all,close all,
load('..\polar_image3.mat')
%% Loading the image of the walls
load('..\main_data_range1_box.mat')
dx = 0.2; dy=0.05;%range and cross range resolution (cm)
                  %scaling factor
[n3 n4]=size(F3_norm); % F3 is the raw reconstructed image
F3_norm=abs(F3_norm);
F5 = zeros(n3,n4);
zmax=((n3-1)*dy*10^-2)/2;

%%% rectangular co-ordinates
xx=(0:n3)*dy*10^-2-zmax;%%offset zmax
yy=(0:n4)*dx*10^-2;
hh=1;
for cc=1:n3   %%cross-range scanning
    for rr=1:n4  %%range scanning
        [t r]= cart2pol(xx(cc),yy(rr));%% polar point conversion
        %%indices formation
        Ri=round(r*10^2/DR)+1; 
        Ri1(hh)=Ri;
        Ti=round(sample*t/pi)+1;
%         for qwer=1:1:length(RR)
%            for zft=1:1:Ntheta
               if Ri<length(RR)
                   F5(cc,rr)=abs(5e4*Fn(Ri,Ti))+F5(cc,rr);
                   
               end
               hh=hh+1;
    end
end
figure(2)
dz = 0.2; % range step in cm
z_start=1; % starting range value in cm
rn=300;
m1=n3;
du_x=0.5*10^-3;
[Zc,Xc] = meshgrid(z_start:dz:rn,(1:m1)*(du_x*100));
% surf(Zc,Xc,(F3_norm));
surf(Zc,Xc,(F5));

shading interp
pbaspect([(rn-z_start)/(m1*du_x*100) 1 1]);
axis tight
view([0 90]);

figure(1)
surf(Z,Y,(abs(Fn)));
shading interp
% pbaspect([1/3 1 1])
xlim([-0.5 0.5])
ylim([0 3])
axis tight
view([0 90]);