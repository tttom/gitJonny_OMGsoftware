close all; clear all;

%initialise variables (all values in SI units)
N=300; % sampling number
L=50*10^-3; % sample size
dx=L/N; % step size
lambda=633*10^-9;
k0=(2*pi)/lambda;
w0=10*10^-3; %input waist
n=1;
zR = (pi*n*w0^2)/lambda; % Rayleigh range
f=20*10^-3; 
z=f;

for n= 1:N+1
    for m=1:N+1
       %Space axis
        x(m)=(m-1)*dx-L/2;
        y(n)=(n-1)*dx-L/2;
        %Frequency axis
        Kx(m)=(2*pi)/L*(m-1) - pi/dx;
        Ky(n)=(2*pi)/L*(n-1)- pi/dx;
    end
end

[X,Y]=meshgrid(x,y);
[KX, KY]=meshgrid(Kx,Kx);

% Gaussian Beam in space domain
gauss=exp(-(X.^2+Y.^2)./(w0^2));
%gauss = w0^2/2 *exp(-(w0^2)/4 *(KX.^2 + KY.^2));

% lens transfer function (real space and reciprocal space eqns)
%t=exp( ((1i*k0)/(2*f)) *(X.^2+Y.^2));
t=(1i*f)/k0 *exp((1i*f)/(2*k0) *(KX.^2 + KY.^2));

% free space transfer function
H=exp(1i/(2*k0)*z*(KX.^2+KY.^2));

%setting radius of curvature of gaussian = inf (const. phase)
Fgauss = fftshift(fft2(gauss));
FgaussR = real(Fgauss);
FgaussI = repmat(0.002, [N+1 N+1]);
Fgauss = complex(FgaussR,FgaussI);

prop1=Fgauss.*t;%after lens
prop2=prop1.*H;%after propagation 

figure(1)
colormap(hot)
subplot(1,2,1)
% surf(x,y,abs(ifft2(prop1)))
imagesc(x,y,abs(ifft2(prop1)));axis image;
title('After lens')
shading flat
xlabel('x');ylabel('y');zlabel('Amplitude');
ylim([-0.025 0.025]); xlim([-0.025 0.025])

subplot(1,2,2)
% surf(x,y,abs(ifft2(prop2)))
imagesc(x,y,abs(ifft2(prop2)));axis image;
title('After propagation')
shading flat
xlabel('x');ylabel('y');zlabel('Amplitude');
ylim([-0.025 0.025]); xlim([-0.025 0.025])

%calculating FWHM of final beam
halfMax =max(max(abs(prop2)))/2;
index1 = find(abs(prop2(:,151))>=halfMax,1, 'first');
index2 = find(abs(prop2(:,151))>=halfMax,1, 'last');
FWHM = (x(index2)-x(index1))/sqrt(2*log(2));
FWHM = FWHM*10^6;
sprintf('Waist after lens [microns], sim. %0.2f',FWHM)

%theoretical waist of final beam
zf=(f*zR^2)/(f^2+zR^2); %position of wf, should be equal to f
wf=sqrt( (f^2*w0^2)/(f^2+zR^2) );
wf=wf*10^6;
sprintf('Waist after lens [microns], theoret. %0.2f',2*wf)
