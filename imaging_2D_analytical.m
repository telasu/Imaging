%%% 2D imaging analytical example using point sources as targets based on
%%% Sheen's 3D imaging paper and Soumkeh theroy 
%%% written by Mohammed Aladsani 2018 email : maladsan@asu.edu 


clc;
clear all;


%%%%%%%%%%%%%%%%%%%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=60*10^9 ; % frequency
Lamb=(3*10^8)/f; %wavelength
k=(2*pi)/Lamb; %wavenumber 

dX=Lamb/4; %sampling distance in X

dY=Lamb/4; %sampling distance in Y

Ltx=1000*Lamb; %lenght of traget space in X or cross range X
Lty=1000*Lamb; %length of traget space in Y or cross range Y


z=5*Lamb; %range or distance from target

Nx=ceil(Ltx/dX); %get total number of smaples = length/sampling distance in X
dKx=(2*pi)/(Nx*dX); %spatial wavenumber in x per unit 
ScaleX=(-Nx/2:(Nx/2)-1); %just a scale for later
X=dX*ScaleX;% total length in x
Kx=dKx*ScaleX;% total spatial wavenumber in x


Ny=ceil(Lty/dY); %get total number of smaples = length/sampling distance in Y
dKy=(2*pi)/(Ny*dY); %spatial wavenumber in y per unit 
ScaleY=(-Ny/2:(Ny/2)-1); %just a scale for later
Y=dY*ScaleY;% total length in y
Ky=dKy*ScaleY;% total spatial wavenumber in y


%Kz=sqrt(((2*k)^2)-((Kx).^2)-((Ky).^2)); %% total spatial wavenumber in z


[Xs,Ys] = meshgrid(X,Y); % to use later for plotting 
[Kxs,Kys]=meshgrid(Kx,Ky); % to use later plotting

%Kzs=sqrt(((2*k)^2)-((Kxs).^2)-((Kys).^2));
%Kzs=sqrt(Kzs2.*(Kzs2>0));

% get Kz for each point (x,y)
for h=1:Nx
    for g=1:Ny
        Kzu(h,g)=4*k^2-Kx(h)^2-Ky(g)^2;
    end
end
Kzs=sqrt(Kzu.*(Kzu>0)); % make sure it's real





%Xs2=6.*Xs;
%Ys2=6.*Ys;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% system or traget creation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%XYn=[0;0]; % traget location each row traget in X Y  
Xn=[0.2;0.25;0.3;-3;-3.1;-3.15;2.5]; %tragert location in x
Yn=[0.2;0.25;0.3;4;4.1;4.15;2.5];%target location in y

fn=[1;1;1;1;1;1;2]; % traget reflectivty function for each target

m=length(Xn); % represent number of targets with should be the same as n with respect of X
n=length(Yn); % represent number of targets with should be the same as m with respect of Y

sys=zeros(Nx,Ny); % just initialization for faster computation 

for p=1:Nx %each ant pos. at x
    for q=1:Ny %%each ant pos. at y
        for a=1:m % traget in x or y since they should be the same 
        wave(a)=fn(a)*exp(-1i*2*k*(sqrt(z^2+((Xs(p,q)-Xn(a))^2)+((Ys(p,q)-Yn(a))^2))));
      %  sys(p,q)=wave(a)+sys(p,q);

        
        end
            
        sys(p,q)=sum(sum(wave)); % system responce for each (x,y) simpling distance 
        
        end
end
    



    % the following in not important just to see the real part of system
    % reaponce 
Rsys=real(sys); 

subplot(2,1,1);
surf(Xs,Ys,Rsys);
colormap jet
xlabel('cross range x');
ylabel('cross range y');
zlabel('real part');
title('received signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sysku=fftshift(fft2((sys))); % 2D-FTT of system function

expos=exp(1i*Kzs*z); %foucsing function

Fn=sysku.*expos;

F=ifft2(Fn); % reconstructed image

RF=real(F);
AF=abs(F); % magnitude of the reflectivy function 

subplot(2,1,2);
surf(Xs,Ys,AF);
colormap jet
xlabel('cross range');
ylabel('magnitude');
title('reconstructed signal');

