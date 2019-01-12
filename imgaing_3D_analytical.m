%%% 3D imaging analytical example using point sources as targets based on
%%% Sheen's 3D imaging paper and Soumkeh theroy 
%%% written by Mohammed Aladsani 2018 email : maladsan@asu.edu 

clc;
clear all;


%%%%%%%%%%%%%%%%%%%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=20*10^9:10^9:45*10^9 ; % bandwidth of freq. to be used
nf=length(f); % range of Freq. how many frequecny points 

%X=cell(nf,1);

Lamb=(3*10^8)./f(nf); %wavelength assuming highest freq. taking from somkeh book and literature to be used later on 

LambK=(3*10^8)./f; % for range of k

k=(2*pi)./LambK; %wavenumbers

dX=Lamb/4; %sampling distance in X assuming highest freq

dY=Lamb/4; %sampling distance in Y assuming highest freq

Ltx=20*Lamb; %lenght of traget space in X assuming highest freq
Lty=20*Lamb; %length of traget space in Y assuming highest freq


%z=15*Lamb; %range or distance from target

Nx=ceil(Ltx/dX); %get total number of smaples = length/sampling distance in X assuming highest freq
dKx=(2*pi)./(Nx.*dX); %spatial wavenumber in x per unit assuming highest freq
ScaleX=(-Nx/2:(Nx/2)-1); %just a scale for later assuming highest freq

X=dX.*ScaleX;% total length in x assuming highest freq
Kx=dKx.*ScaleX;% total spatial wavenumber in x assuming highest freq

Ny=ceil(Lty/dY); %get total number of smaples = length/sampling distance in Y assuming highest freq
dKy=(2*pi)./(Ny*dY); %spatial wavenumber in y per unit assuming highest freq
ScaleY=(-Ny/2:(Ny/2)-1); %just a scale for later assuming highest freq

Y=dY.*ScaleY;% total length in y assuming highest freq
Ky=dKy.*ScaleY;% total spatial wavenumber in y assuming highest freq


[Xs,Ys] = meshgrid(X,Y); %assuming highest freq for plotting
[Kxs,Kys]=meshgrid(Kx,Ky); %assuming highest freq for plotting



% get kz for each freq. points 
for pp=1:nf
for hh=1:Nx
    for gg=1:Ny
       Kzs(hh,gg,pp)=4*k(pp)^2-Kx(hh)^2-Ky(gg)^2; % get Kz for each layer aka freq.

        KzsA(hh,gg)=4*k(pp)^2-Kx(hh)^2-Ky(gg)^2; % get Kz for each layer aka freq.
        
        
    end
   
end
kza{pp}=KzsA;
end
Kzs2=sqrt(Kzs.*(Kzs>0)); % get real value of Kz since Kz is real

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% system or target creation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%XYn=[0;0]; % traget location each row traget in X Y  
Xn=[0;0.04;-0.06]; %traget location in X
Yn=[0;0.04;-0.05];%traget location in Y
Zn=[0.15;0.1;0.4];%traget location in Z
fn=[1;1;1]; % traget reflectivty function 

m=length(Xn);
n=length(Yn); % m should = n 

%sys=zeros(Nx,Ny);
for mm=1:nf
for p=1:Nx %each ant pos. at x
    for q=1:Ny %%each ant pos. at y
        for a=1:m % traget in x or y 
        wave(a)=fn(a).*exp(-1i*2*k(mm)*(sqrt((-Zn(a)).^2+((Xs(p,q)-Xn(a)).^2)+((Ys(p,q)-Yn(a)).^2)))); % wave at each pos of (x,y,z)
      %  sys(p,q)=wave(a)+sys(p,q);

        
        end
        sys(p,q,mm)=sum(sum(wave)); % the sum of wavs for each pos of ant. (X,Y,0)
        
        
    end
      
end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for ooo=1:nf
sysku=fftshift(fft2(fftshift((sys)))); % 2D-FTT of system function
%end
%expos=exp(1i.*Kzs{1}.*z(1)); %foucsing function



%%%%%%%%%%%%%%%%%%%%%%%%% generating 3D-IFDT here it basiclly scans each layer of (x,y,freq.) and foucs for each range points of interest the more samples in range used here the more
% accurate you get but more computational heavy but since here in this example we know what range the target just make sure the points in traget target range are within
% this search or they will not be seen  
i=1;
for zt=0.1:0.05:0.4 % searching in range 
expos=exp(1i.*Kzs2*zt); % foucsing functions 
%%%%reconstruction%%%
F=sysku.*(expos);
F2=zeros(Nx,Ny);
sum=zeros(Nx,Ny); % for compute perpouses 
for p=1:nf % for each freq. points
    sum=F(:,:,p);
    F2=sum+F2;
end
%F2=F2./nf; %scaling factor
F1=fftshift(ifft2(fftshift(F2)));,%[Nx Ny])); %take 2DIFDT of each slice
Fn(:,:,i)=F1; %saving into matrix
i=i+1;
end





%Fn=sysku.*expos;

%F=ifft2(Fn); % reconstructed image



%%%%%%%%%%%%%%%%%%%%%%%% ploting %%%%%%%%%%%%%%%%%%%%%%

RF=real(Fn);
AF=abs(Fn);

Z=0.1:0.05:0.4;
[XX YY ZZ]=meshgrid(X,Y,Z);
slice(XX,YY,ZZ,AF,[],[],Z);
alpha(0.2)
shading interp
colorbar
colormap jet
xlabel('x');
ylabel('y');
zlabel('z');
title('reconstructed signal');

