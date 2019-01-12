
%%% 1D imaging analytical example using point sources as targets based on
%%% Sheen's 3D imaging paper and Soumkeh theroy 
%%% written by Mohammed Aladsani 2018 email : maladsan@asu.edu 



clc;
clear all;


%%%%%%%%%%%%%%%%%%%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=3*10^9 ; % frequency
Lamb=(3*10^8)/f; %wavelength
k=(2*pi)/Lamb; %wavenumber 

dX=Lamb/50; %sampling distance

Lt=20*Lamb; %length of tagret space or cross range or scaning arputre 
z=5*Lamb; %range or distance from target
N=ceil(Lt/dX); %get total number of smaples = length/sampling distance
dKx=(2*pi)/(N*dX); %spatial wavenumber in x per unit

Scale=(-N/2:(N/2)-1); %just a scale for later

X=dX*Scale;% total lenght in x
Kx=dKx*Scale;% total spatial wavenumber in x

Kz=sqrt(((2*k)^2)-((Kx).^2)); %% total spatial wavenumber in z or range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% system responce or target creation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xn=[0.75;-0.5]; % traget location in cross range 
fn=[1;1]; % traget reflectivty function 

[m,n]=size(Xn); % m is rows and n is colums

for p=1:N %each ant pos.
    for q=1:m %each traget 
        wave(q)=fn(q)*exp(-1i*2*k*sqrt(z^2+(X(p)-Xn(q))^2));
        
 %       w=fn(q)*exp(-1i*2*k*sqrt(z^2+(X(p)-Xn(q))^2));
%sys(p)=w+sys(p);
        
        
    end
    sys(p)=sum(sum(wave)); % the system response at each sampling distance 
end

subplot(2,1,1);
plot(X./Lamb,real(sys)); %just a plot of the real part of the signal
xlabel('cross range');
ylabel('real part');
title('received signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sysku=fftshift(fft((sys))); % FFT of the system Responce

expos=exp(i*Kz*z); % the "foucsing function "

Fn=sysku.*expos; 

F=ifft(Fn); % IFFT to get the target function 

%RF=real(F);
AF=abs(F); % get the magnitude of the target function 

subplot(2,1,2);
plot(X,AF);
xlabel('cross range');
ylabel('magnitude');
title('reconstructed signal');

