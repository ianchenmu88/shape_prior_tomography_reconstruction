%This programe is for a set of parameters, registrate using the Old methods

arraytheta = [9, 30, 17, 62, 51];
arrayb = [10, -8; -15, 21; -30, 20; 16, -27; 11, -25];

ADDNOISE = false;

I = imread('T1.tif');
Irv = imread('T1dd.tif');

for i=1:3


I = mat2gray(I);
Irv = mat2gray(Irv);
[I1 I2] = myimt(I, 0, 1, [0, 0]);
disp('theta is:');
arraytheta(i)
disp('vector b is')
arrayb(i,:)

[Irv1 Irv2] = myimt(Irv,arraytheta(i), 1, arrayb(i,:));

%This programe can find the parameter of translate using the music method.
disp('The size of the input image is: ');
size(I)
disp('The size of the image wrap the image I is ');
size(I1)

I2 = Irv2;   %Very Important
%imshow(I1); figure; imshow(I2);figure;


%Compute the Radon transform for image I1 and I2

tic;

theta = 0:179;

R1 = radon(I1,theta);
R2 = radon(I2,theta);


if ADDNOISE == true   % IF need add noise to image, then process into this segment.
R1 = (R1-min(R1(:)))/(max(R1(:))-min(R1(:)));   %normalized for display
R2 = (R2-min(R2(:)))/(max(R2(:))-min(R2(:)));   
psnr = 25;
sigma = 1/power(10, psnr/10);
R1 = imnoise(R1, 'gaussian', 0, sigma);
R2 = imnoise(R2, 'gaussian', 0, sigma);
end
R1 = (R1-min(R1(:)))/(max(R1(:))-min(R1(:)));   %normalized for display
R2 = (R2-min(R2(:)))/(max(R2(:))-min(R2(:)));   






disp('The size of the Radon image');
size(R1)

R1b = R1;                  %R1 and R2 are stored for find parameter b
R2b = R2;

for i=1:size(R1, 2)
    R1(:, i) = abs(fft(R1(:, i)));
    R2(:, i) = abs(fft(R2(:, i)));
end
  

m = size(R1,1);
n = size(R1,2);

f = 10;
v1 = R1(f,:);
v2 = R2(f,:);

e_theta = MCPT(v1, v2);
e_theta

%Shif the image by the Theta computed from e_theta

ae_theta  = floor(e_theta+0.5);  %Convert the e_theta to integer.
R2b = circshift(R2b, [0 -ae_theta]);   % shift at x direction to inverse direction by etheta
if ae_theta > 0                        
    for j=1:ae_theta
        t = R2b(:, 180-j+1);
        R2b(:, 180-j+1) = t(end:-1:1);
    end

else
    for j = ae_theta:-1
        t = R2b(:, j-ae_theta+1);
        R2b(:,j-ae_theta+1) = t(end:-1:1);
    end
        
end

R1b = (R1b-min(R1b(:)))/(max(R1b(:))-min(R1b(:)));   %normalized for display
R2b = (R2b-min(R2b(:)))/(max(R2b(:))-min(R2b(:)));   
%figure; imshow(R1b); figure; imshow(R2b); figure;
%Compute translation in each colum.
tphi = zeros(size(R1b,2),1);   %  Mtheta*e_b = tphi find e_b using least square
l = size(R1b,1);
for i=1:size(R1b,2)            % iterate for each colum (angle) in R
    d = CPT(R1b(:, i), R2b(:, i));
    tphi(i) = d;
end


Mtheta = zeros(length(tphi),2);
for i=1:length(tphi)
    Mtheta(i,1) = cos(-(i-1)/180.0*pi);
    Mtheta(i,2) = sin(-(i-1)/180.0*pi);
end

e_b = lscov(Mtheta, tphi);             %Mtheta * e_b = tphi  => e_b = Mtheta/tphi

e_b

toc;

end

