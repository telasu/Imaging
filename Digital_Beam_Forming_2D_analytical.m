
%%%%% Digital Beam Forming to find AoA " Angle of Arival " and range
%%%%% example of analytical model using point charge as tragets
%%%% written by Mohammed Aladsani 2018 email : maladsan@asu.edu 




clear;
clc;
tic;


% setup
N = 150; % Number of Elements in the imaging system\comms.
n = 0:(N-1); % Antenna Array Vector
f=200.03*10^9:3*10^7:350*10^9 ; % frequency of operation 
nf=length(f); % range of Freq.
c=3*10^8; %speed of light

lamd=c/f(nf); %wavelenght
kd=(2*pi)/lamd; %wavenumber
d=0.5*lamd; %spacing between each element in array

%sampling_array=(1*10^-3)/lamd


%%%%%%%%%%%%%% scene\traget definiton %%%%%%%%

r1=20*10^-2; %raidal distance from first traget till the array in meters
r2=12*10^-2; %raidal distance from 2nd traget till the array in meters
theta1=45; %angle of arrival of the first target " results in the end should match this one"
theta2=120; %angle of arrival of the 2nd target " results in the end should match this one"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j=1:1:nf %%%% for each freq. point 
    
    ff=f(j);
lam=c/ff; %wavelenght
k=(2*pi)/lam; %wavenumber
kk(j)=k;







E1=1*exp(i*k*r1); %plane wave going into -ve raidal direction , for simplcity ignore 1/r factor as it doesnt imapct the results for first traget  
E2=1*exp(i*k*r2); %plane wave going into -ve raidal direction , for simplcity ignore 1/r factor as it doesnt imapct the results for 2nd traget 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



sai1 = k*d*cos(theta1*pi/180); %phase factor of the plane wave accounting for AoA 
sai2= k*d*cos(theta2*pi/180);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%phase_induced_upon_the_signal per array element or you could call it
%normalized array responce vector

PE1 = exp(-1*i*(n - (N-1)/2)*sai1);
PE2 = exp(-1*i*(n - (N-1)/2)*sai2);





%%%%%%%%%%%%%% total responce of the array %%%%%%%%%%%%%%%%%%%
RE1 = E1*PE1; % traget 1
RE2= E2*PE2; % target 2


RE=RE1+RE2; %% all signals 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%% weighting for each ULA element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sweep Angle Theta
sample = 1000; % how many angle samples you want
theta_d = 0:(180/sample):180;
sai_d = k*d*cos(theta_d*pi/180); % phase for each sample
Ntheta=length(theta_d);
% Weight Coefficients
w = (1/N)*exp(i*(n.'-(N-1)/2)*sai_d); % weight factor matrix for each element with all possible sample angles









% sweeping the Beampattern to find AoA 

for ii=1:length(theta_d) % for all samples 

 w_sai=w(:,ii); %for all elements in array for a search angle 
B(ii) = RE*w_sai; %responce after multiplying if angle is AoA this should be maximum in magnitude 

end

BF(j,:)=B; %% saving into new matrix

end






%%%%%%% for ranging similar to imaging 

for ooo=1:nf
sysku(ooo,:)=fftshift(fft(fftshift((BF(ooo,:))))); % 2D-FTT of system function
end


iii=1;
for R=0:0.5*10^-2:23*10^-2 %%% possible ranges search 
expos=exp(1i.*kk*R)';

%%%%reconstruction%%%
F0=sysku.*(expos);




F2=zeros(1,Ntheta);

for p=1:nf
 %   for p=1:10000
   sum=F0(p,:);
  
    F2(1,:)=sum+F2(1,:);
end
%F2=F2./nf; %scaling factor
F1=fftshift(ifft(fftshift(F2)));
Fn(iii,:)=F1; %saving into matrix
iii=iii+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% ploting %%%%%%%%%%%%%%

RR=0:0.5*10^-2:23*10^-2;

[ThetaDD RRR] = meshgrid(theta_d,RR);

surf(RRR,ThetaDD,20*log(abs(Fn)))
shading interp

xlabel('Radius');
ylabel('Theta ');
title(' Analytical data FFT - Reconstructed image ');

shading interp;
axis tight;
view([0 90]);
pbaspect([(23*10^-2)/(150*d) 1 1]);
%theta_d(index)
toc

%figure(1)

%figure(2)
%plot(theta_d,20*log10(abs(B)));
