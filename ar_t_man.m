function trianglePixels = ar_t_man(i)

a = zeros(150,150);
a(floor((150-i*100)/2),75) = 1;
%a(80*i+20,50) = 1;

b = 1:i*100;

for j = 1:size(b,2)
a(floor((150-i*100)/2)+b(j),floor(75-(50/size(b,2)*j)):floor(75+(50/size(b,2)*j))) = 1;
end

a=a(21:129,21:129);

a=imresize(a,[344,344]);


b = zeros(344*3,344*3);

b(345:(345+343),345:(345+343))=a;

b=imresize(b,[344,344]);

trianglePixels = ~b;
%imshow(a, [])
end


