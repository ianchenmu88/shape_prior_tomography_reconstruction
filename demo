clear all;
close all;

theta = 25:154;
sigma = 15;

for S = 1:50
    % Specify following values
    arrayscale = [1]; 
    arrayar = [1];
    arraytheta = [1]; 
    arrayb = [1, 1];

    filename_noise = ['noise_fixed_15', '_', num2str(S), '.mat'];
    load(filename_noise, 'noise');

    Irv = ar_1(arrayar);
    Irv = ~mat2gray(Irv);
    Pixels = Irv;

    Irv2 = myimt(Irv, arraytheta, arrayscale, arrayb);
    R_radon = radon(Irv2, theta); 
    R_radon = R_radon + noise;         

    % Determine rotation bounds
    if 0 < arraytheta && arraytheta <= 45
        a3_0 = 0;
        b3_0 = 45;
    elseif 45 < arraytheta && arraytheta <= 90
        a3_0 = 45;
        b3_0 = 90;
    elseif 90 < arraytheta && arraytheta <= 135
        a3_0 = 90;
        b3_0 = 135;
    else
        a3_0 = 135;
        b3_0 = 180;
    end

    % Determine translation bounds for arrayb
    if -100 <= arrayb(1) && arrayb(1) <= 0
        a1_0 = -100;
        b1_0 = 0;
    else 
        a1_0 = 0;
        b1_0 = 100;
    end

    if -100 <= arrayb(2) && arrayb(2) <= 0
        a2_0 = -100;
        b2_0 = 0;
    else 
        a2_0 = 0;
        b2_0 = 100;
    end

    % Part 1
    g = @(x, y, z) -(log(1/3) - 1/2 * (sigma)^(-1) * size(Pixels, 1) - ...
        1/2 * (sigma)^(-1) * norm(R_radon - (radon(myimt(~mat2gray(ar_1(y)), z, x, [0, 0]), theta)), 'fro')^2 - ...
        size(Pixels, 1) / 2 * log(2 * pi));

    a1 = 0; % Lower bound for variable x (scaling)
    b1 = 2.5; % Upper bound for variable x
    a2 = 0; % Lower bound for variable y (aspect ratio)
    b2 = 1; % Upper bound for variable y
    a3 = a3_0;
    b3 = b3_0;
    epsilon = 0.1; % Termination criteria 
    tau = double((sqrt(5) - 1) / 2); % Golden number 
    k = 0; % Number of iterations 

    % Initial guesses
    x1 = a1 + (1 - tau) * (b1 - a1); 
    x2 = a1 + tau * (b1 - a1);
    y1 = a2 + (1 - tau) * (b2 - a2); 
    y2 = a2 + tau * (b2 - a2); 
    z1 = a3 + (1 - tau) * (b3 - a3); 
    z2 = a3 + tau * (b3 - a3);

    % Points A, B, C, D, E, F, G, H
    ek = [x1, y1, z1]; 
    fk = [x1, y2, z1]; 
    gk = [x2, y1, z1]; 
    hk = [x1, y2, z2]; 
    ik = [x2, y1, z1]; 
    jk = [x2, y2, z1]; 
    kk = [x2, y1, z2]; 
    lk = [x2, y2, z2]; 

    % Function values at points
    gek = g(x1, y1, z1); 
    gfk = g(x1, y2, z1);
    ggk = g(x1, y1, z2); 
    ghk = g(x1, y2, z2); 
    gik = g(x2, y1, z1); 
    gjk = g(x2, y2, z1);
    gkk = g(x2, y1, z2); 
    glk = g(x2, y2, z2); 

    % Optimization loop
    while sqrt((b1 - a1)^2 + (b2 - a2)^2 + (b3 - a3)^2) > epsilon
        k = k + 1;
        min1 = min([gek, ghk, gfk, ggk, gik, gjk, gkk, glk]); 

        if min1 == gek 
            b1 = x2; b2 = y2; b3 = z2; 
        elseif min1 == gfk 
            b1 = x2; a2 = y1; b3 = z2; 
        elseif min1 == ggk 
            b1 = x2; b2 = y2; a3 = z1; 
        elseif min1 == ghk 
            b1 = x2; a2 = y1; a3 = z1; 
        elseif min1 == gik 
            a1 = x1; b2 = y2; b3 = z2; 
        elseif min1 == gjk 
            a1 = x1; a2 = y1; b3 = z2; 
        elseif min1 == gkk 
            a1 = x1; b2 = y2; a3 = z1; 
        elseif min1 == glk 
            a1 = x1; a2 = y1; a3 = z1; 
        end

        % Update guesses
        x1 = a1 + (1 - tau) * (b1 - a1); 
        x2 = a1 + tau * (b1 - a1);
        y1 = a2 + (1 - tau) * (b2 - a2); 
        y2 = a2 + tau * (b2 - a2); 
        z1 = a3 + (1 - tau) * (b3 - a3); 
        z2 = a3 + tau * (b3 - a3); 

        % Function evaluations
        gek = g(x1, y1, z1); 
        gfk = g(x1, y2, z1); 
        ggk = g(x1, y1, z2); 
        ghk = g(x1, y2, z2); 
        gik = g(x2, y1, z1); 
        gjk = g(x2, y2, z1); 
        gkk = g(x2, y1, z2); 
        glk = g(x2, y2, z2); 

        min1 = min([gek, ghk, gfk, ggk, gik, gjk, gkk, glk]); 
    end 

    % Determine the minimum point
    if min1 == gek
        a = x1; b = y1; c = z1;
    elseif min1 == gfk
        a = x1; b = y2; c = z1;
    elseif min1 == ggk
        a = x1; b = y1; c = z2;
    elseif min1 == ghk
        a = x1; b = y2; c = z2;
    elseif min1 == gik
        a = x2; b = y1; c = z1;
    elseif min1 == gjk
        a = x2; b = y2; c = z1;
    elseif min1 == gkk
        a = x2; b = y1; c = z2;
    elseif min1 == glk
        a = x2; b = y2; c = z2;
    end

    % Part 2
    f = @(x, y) norm(R_radon - (radon(myimt(~mat2gray(ar_1(b)), c, a, [x, y]), theta)), 'fro')^2; 
    a1 = a1_0; % Lower bound for variable x 
    b1 = b1_0; % Upper bound for variable x 
    a2 = a2_0; % Lower bound for variable y 
    b2 = b2_0; % Upper bound for variable y 
    epsilon = 0.1; % Termination criteria 
    tau = double((sqrt(5) - 1) / 2); % Golden number

    k = 0; % Number of iterations
    x1 = a1 + (1 - tau) * (b1 - a1); 
    x2 = a1 + tau * (b1 - a1);
    y1 = a2 + (1 - tau) * (b2 - a2); 
    y2 = a2 + tau * (b2 - a2);

    % Points A, B, C, D
    ek = [x1, y1]; 
    fk = [x1, y2]; 
    hk = [x2, y1]; 
    gk = [x2, y2]; 

    % Function values at points
    fek = f(x1, y1);
    ffk = f(x1, y2);
    fhk = f(x2, y1);
    fgk = f(x2, y2);

    while sqrt((b1 - a1)^2 + (b2 - a2)^2) > epsilon
        k = k + 1;
        min1 = min([fek, fhk, ffk, fgk]); 

        if min1 == fek
            b1 = x2; 
            b2 = y2;
        elseif min1 == ffk
            b1 = x2; 
            a2 = y1;
        elseif min1 == fgk
            a1 = x1; 
            a2 = y1;
        elseif min1 == fhk
            a1 = x1; 
            b2 = y2;
        end

        % Update guesses
        x1 = a1 + (1 - tau) * (b1 - a1); 
        x2 = a1 + tau * (b1 - a1);
        y1 = a2 + (1 - tau) * (b2 - a2); 
        y2 = a2 + tau * (b2 - a2);

        % Function evaluations
        fek = f(x1, y1);
        ffk = f(x1, y2);
        fhk = f(x2, y1);
        fgk = f(x2, y2);
        
        min1 = min([fek, fhk, ffk, fgk]); 
    end 

    % Determine the minimum point
    if min1 == fek
        m = x1; n = y1;
    elseif min1 == ffk
        m = x1; n = y2;
    elseif min1 == fhk
        m = x2; n = y1;
    elseif min1 == fgk
        m = x2; n = y2;
    end 

    % Output results
    fprintf('scale: %f\n', arrayscale);
    fprintf('scale estimated: %f\n', a);
    fprintf('aspect ratio: %f\n', arrayar);
    fprintf('aspect ratio estimated: %f\n', b);
    fprintf('rotation: %f\n', arraytheta);
    fprintf('rotation estimated: %f\n', c);
    fprintf('translation x: %f\n', arrayb(1));
    fprintf('translation x estimated: %f\n', m);
    fprintf('translation y: %f\n', arrayb(2));
    fprintf('translation y estimated: %f\n', n);
end
