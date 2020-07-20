function y = fsphbes2(i,z)
    % i : indeks Bessel/Hankel
    % z : argumen Bessel/Hankel
    y = sqrt(pi/(2*z))*bessely(i+0.5,z);
end