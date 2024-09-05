function ellipsePixels = ar_1(i)

%Adjust for ellipse%
Scaling_Factor_1 = 1;
Scaling_Factor_2 = 1/i;

Image=zeros(344,344);

%Importing image%
% Image2 = imread('np3_2.jpg');
% Image2 = rgb2gray(Image2);
% Image2 = imresize(Image2,[344,344]);
% 
% Image(1:114,1:114) = Image2(1:114,1:114);
% Image(1:114,115:228) = Image2(1:114,1:114);
% Image(1:114,229:344) = Image2(1:114,1:116);
% 
% Image(115:228,1:114) = Image2(1:114,1:114);
% Image(115:228,115:228) = Image2(1:114,1:114);
% Image(115:228,229:344) = Image2(1:114,1:116);
% 
% 
% Image(229:344,1:114) = Image2(1:116,1:114);
% Image(229:344,115:228) = Image2(1:116,1:114);
% Image(229:344,229:344) = Image2(1:116,1:116);

%The radius of the filled circle%
Radius = 50;

%Grabbing the dimensions of the image%
[Image_Height,Image_Width] = size(Image);

%Evaluating the midpoint of the image%
Image_Midpoint = [round(Image_Height/2), round(Image_Width/2)];

%Scanning through the pixels of the image%
for Row_Scanner = 1: +1: Image_Height
   for Column_Scanner = 1: +1: Image_Width 
      
      
      
%The pixel coordinate%
Pixel_Coordinate = [Row_Scanner Column_Scanner];
Magnitude_From_Centre = sqrt((abs(Pixel_Coordinate(1,1) - Image_Midpoint(1,1))^2)/Scaling_Factor_1 + abs(Pixel_Coordinate(1,2) - Image_Midpoint(1,2))^2/Scaling_Factor_2);



%If the magnitude from the centre is smaller than the set radius set
%intensity to 0%
if(Magnitude_From_Centre <= Radius) 
    Image(Row_Scanner,Column_Scanner) = 1; 
end 

   end 
end

Image = ~Image;
%imshow(Image,[]);

ellipsePixels = Image;
%ellipsePixels = imresize(ellipsePixels,[172,172]);
end
