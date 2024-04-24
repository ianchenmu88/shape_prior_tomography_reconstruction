clc
clear all
close all


angle = [25:154];

gamma = 0.54
s = 0.98
theta = 1*180/pi      %radian*180/pi
a = 155.47
b = 176.66
g_snr = 5              %5:31, 15:22,  25:17

ori = myimt(~ar_t_man(s), theta, gamma, [a-122, 222-b]);

ori = imresize(ori,[344,344]);

ori_proj = radon(ori, angle) ;

ori_proj = ori_proj +  normrnd(0,g_snr,[491,130]);

%imshow(ori_proj,[])

%%  SIRT-TV-finite support reconstruction

P_num = 344;
Nx = P_num;
Ny = P_num;

Img_size = 344;   % in mm
Pixel_size = Img_size/P_num;


lamda_ART = 1;


N_iter = 10;


Ngrad = 3;      
geoOffset = 0;      
alpha_Grad = 0.08;  

Part0 = zeros(Ny, Nx);
Part1 = Part0;


theta = [25:154]; %projection angle 
F=double(ori);

R = ori_proj ;

At = iradon(R,theta,'linear', 'none', 1,344); %reconstruct noisy alien

n_sirt = 50;%iterations

Fk = At;%Matrix Fk is our solution at the k-th step, now it is our initial guess
At=(At-min(At(:)))/(max(At(:))-min(At(:)));


for n = 1:N_iter
   
    n
    



for  k=1:n_sirt
    
    t = iradon(radon(Fk,theta),theta, 'linear', 'none', 1,344);% reconstruct alien using Fk unfiltered sinogram
    t=(t-min(t(:)))/(max(t(:))-min(t(:)));
    
    Fk = Fk + At - t;
    
    Fk=(Fk-min(Fk(:)))/(max(Fk(:))-min(Fk(:)));
  
end

    level = graythresh(Fk)
    BW = imbinarize(Fk,level);

    Part1 = BW;
    
        % ----------------- Positivity constraint ----------------- %
    negtive_pixel_ind = find(Part1<0);
    if ~isempty(negtive_pixel_ind)
        Part1(negtive_pixel_ind) = 0;
    end    
    
    %     % -------- TV gradient descent ----------- %
    dA = sqrt(sum(sum((Part0-Part1).^2)));
    Ptv_grad = Part1;
    for m = 1: Ngrad
        V_vect =  gradient(Ptv_grad);
        V_vect_norm = V_vect/norm(V_vect);

        Ptv_grad = Ptv_grad - alpha_Grad*dA*V_vect_norm;
    end
    Part1 = Ptv_grad;
    Part0 = Ptv_grad;


    Fk = Part1;


end


imshow(Fk,[])

Fk = (Fk-min(Fk(:)))/(max(Fk(:))-min(Fk(:)));
ori = (ori-min(ori(:)))/(max(ori(:))-min(ori(:)));

mse = mean((Fk(:)-ori(:)).^2)

level = graythresh(Fk);
BW = imbinarize(Fk,level);
mask_ = boundarymask(BW);
level2 = graythresh(ori);
BW2 = imbinarize(ori,level2);
mask = boundarymask(BW2);
HausdorffDist = HausdorffDist(mask,mask_)


