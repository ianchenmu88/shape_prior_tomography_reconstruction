function t = MCPT(t1, t2)
% This function is for estimating the translation between two 1D signals.
% The method is the modified phase correlation. 

n = length(t1);

Ft1 = fft(t1);
Ft2 = fft(t2);

MCPS = (Ft2.*conj(Ft1))./abs(Ft1);
 [~, t] = max(real(ifft(MCPS)));
 
 figure; plot(real(ifft(MCPS)),'LineWidth',2);
 set(gca, 'FontSize', 20 );
 print -depsc vvvv.eps ;
 
 %real(ifft(CPS))
if t > n
    t =  t-2*n-1; 
end

disp('In MCPT');
