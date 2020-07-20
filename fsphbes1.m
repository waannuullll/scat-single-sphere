function j = fsphbes1(i,z)
    % i : indeks Bessel/Hankel
    % z : argumen Bessel/Hankel
    j = sqrt(pi/(2*z))*besselj(i+0.5,z);
end