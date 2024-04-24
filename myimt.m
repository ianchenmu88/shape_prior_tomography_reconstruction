function [I1 I2]=myimt(I, theta, s, b)
%This function set two image for registration, 
% I, is an input image.
% theta, is the angle to be rotated.
% s, is the scale parameter.
% b, is the translate parameter, and b = [b1, b2].
% Put the I in the meddle of image I1;
% Put the translated I in the image I2;


%imshow(I) 
theta = theta*pi/180;

Mr = [ cos(theta)  -sin(theta)  0
       sin(theta)  cos(theta)  0
           0              0        1 ];
[m,n] = size(I);
Mb = [1   0   0
      0   1   0
      -b(1) -b(2)  1];
 
Ms = [s  0  0
       0  s  0
       0  0  1];
% set the original at the middle of the image using Mt.

Mt = [1                      0        0
      0                      1        0
      -floor((n+1)/2) -floor((m+1)/2) 1];   
   
M = Mt*Ms*Mb*Mr;

tformI1 = maketform('affine', Mt);
tformI2 = maketform('affine', M);

marginal = 50;
I1 = imtransform(I, tformI1,...
                'XData', [-n/2-marginal n/2+marginal],...
                'YData', [-m/2-marginal m/2+marginal], 'XYScale',1);

            
I2 = imtransform(I, tformI2,...
                'XData', [-n/2-marginal n/2+marginal],...
                'YData', [-m/2-marginal m/2+marginal], 'XYScale',1);
      
            
            
%figure, imshow(I1), figure, imshow(I2)


