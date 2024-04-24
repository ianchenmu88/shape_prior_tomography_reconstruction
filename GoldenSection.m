%%--- Golden Section algorithm --- %%
%%----- Made by Keya Ghonasgi -----%%
% Algorithm to minimize a given single variable function. Requires a
% function handle, initial guess, step size, and epsilon for stopping
% criteria. Default step size is 1, and eps is 0.01.
%  [x_st,f_st,iter] = GoldenSection(@func_name,x0,s,eps) where func_name is
%  the function to be minimized saved as a separate .m file.

function [x_st,f_st,iter] = GoldenSection(fun,x0,s,eps)

if nargin < 4
    eps = 0.01;
end
if nargin < 3
    s = 1;
end
if nargin < 2
    error('Not enough arguments!')
end

iter = 1;   % Initializing the loop
tau = 0.618;

x1(iter) = x0; x2(iter) = x1(iter) +  s;    % Initializing x1 and x2
f1 = fun(x1(iter));
f2 = fun(x2(iter));

% Check if the search needs to change direction
if f2 > f1
    x2 = x1 + x2; f2 = f1 + f2;
    x1 = x2 - x1; f1 = f2 - f1;
    x2 = x2 - x1; f2 = f2 - f1;
    s = -s;
end

s = s/tau;
x4(iter) = x2(iter) + s;    % Intitializing x4
f4 = fun(x4(iter));

if f4 < f2
    x1(iter) = x2(iter); f1 = f2; x2(iter) = x4(iter); f2 = f4;
end

diff = eps + 1; % Checking if limit is reached

while diff > eps    % stopping criteria
    iter = iter + 1;    % Iteration counter
    x3(iter) = tau*x4(iter - 1) + (1 - tau)*x1(iter - 1);
    f3 = fun(x3(iter));
    if f2 < f3      % if f2<f3, 1 -> 4, 3 -> 1
        x4(iter) = x1(iter - 1); f4 = f1;
        x1(iter) = x3(iter); f1 = f3;
        x2(iter) = x2(iter - 1);
    else            % else, 2 -> 1, 3 -> 2
        x1(iter) = x2(iter - 1); f1 = f2;
        x2(iter) = x3(iter); f2 = f3;
        x4(iter) = x4(iter - 1);
    end
    diff = abs(x1(iter) - x4(iter));    % Check stopping criteria
end

x_st = x1(iter);
% x_st = [x1;x2;x3;x4];     % Uncomment this to see how the iterations run through values
f_st = f1;

end