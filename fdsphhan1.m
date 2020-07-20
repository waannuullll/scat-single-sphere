function dh1 = fdsphhan1 (i,z)
    % i : indeks Bessel/Hankel
    % z : argumen Bessel/Hankel
    dh1 = 0.5*(fsphhan1(i-1,z)-(fsphhan1(i,z)+z*fsphhan1(i+1,z))/z);
end