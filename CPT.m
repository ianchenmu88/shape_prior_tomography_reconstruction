function t = CPT(t1, t2)
% This function is for estimating the translation between two 1D signals
% using phase correlation method. 

n = length(t1);

Ft2 = fft(t1, 2*n);
Ft1 = fft(t2, 2*n);

CPS = (Ft2./abs(Ft2))./(Ft1./abs(Ft1));
 [~, t] = max(real(ifft(CPS)));
 
 %real(ifft(CPS))
if t > n
    t =  t-2*n-1; 
end
