function db1 = fdsphbes1(i,z)
    % i : indeks Bessel/Hankel
    % z : argumen Bessel/Hankel
    db1 = fsphbes1(i-1,z) - (i+1)/z*fsphbes1(i,z);
end