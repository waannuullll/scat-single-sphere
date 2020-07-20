function h1 = fsphhan1(i,z)
    % i : indeks Bessel/Hankel
    % z : argumen Bessel/Hankel
    h1 = fsphbes1(i,z) + fsphbes2(i,z)*1j;
end