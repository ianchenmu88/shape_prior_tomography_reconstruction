function recPixels = ar_r(i)

ar_true = i;

ar = (1-ar_true)/2;

a = zeros(344,344);

a(32:313,(32+floor(ar*261)):(313-floor(ar*261))) = 1;

b = zeros(344*3,344*3);

b(345:(345+343),345:(345+343))=a;
recPixels = ~b;

recPixels = imresize(recPixels,[344,344]);



%imshow(recPixels,[])
end
