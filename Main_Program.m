%{
    IM - 10216079
    Menghitung Cross-Section untuk 1 bola dengan J_C
    Semua satuan dalam SI
    Dibuat 13 Okt 19

    Benchmark: El-Sayed

    Catatan bimbingan:
    - 15 Okt 19 (Rev 1): Udah sesuai untuk pola, cuma besarnya masih 2x lipat
    - 22 Okt 19 (Rev 2): Masih lebih dikit tapi dimaklumi, salah lupa kali
    N di bilangan gelombang
%}
clc
clear all
%% Pemilihan Material
% 2-3 Copper, 4-5 Silver, 6-7 Gold
load J_C;                                   % Load file
[baris,kolom] = size(J_C);                  % Ukuran J_C
benda = 6;                                  % Pilih benda
N1 = J_C(:,benda) + J_C(:,benda+1)*1i;      % Indeks bias bola
N = 1.33;                                   % Indeks bias bg
m = N1/N;                                   % Relative refractive index
%% Definisi Variabel
c = 299792458;                              % Kecepatan cahaya
h = 6.626e-34;                              % Konstanta Planck
el = 1.602e-19;                             % Muatan elektron
En = J_C(:,1);                              % Energi (eV)
lambda = h*c./(En*el);                      % Panjang gelombang
k = 2*pi*N./lambda;                         % Bilangan gelombang
a = 40E-9;                                  % Jejari bola
x = 2*pi*N*a./lambda;                       % Size parameter
error = 1E-12;                              % Konvergen
atas = 800;                                 % Batas atas plot lambda
bawah = 400;                                % Batas bawah plot lambda
%% Definisi Fungsi
% i : indeks Bessel/Hankel
% j : indeks x(j)
% z : argumen Bessel/Hankel
psi = @(i,z) (z*fsphbes1(i,z));
xi = @(i,z) (z*fsphhan1(i,z));
dpsi = @(i,z) (fsphbes1(i,z)+z*fdsphbes1(i,z));
dxi = @(i,z) ((fsphhan1(i,z))+(z*fdsphhan1(i,z)));
an = @(i,j) ((m(j)*psi(i,m(j)*x(j))*dpsi(i,x(j))-psi(i,x(j))*dpsi(i,m(j)*x(j)))/...
    (m(j)*psi(i,m(j)*x(j))*dxi(i,x(j))-xi(i,x(j))*dpsi(i,m(j)*x(j))));
bn = @(i,j) ((psi(i,m(j)*x(j))*dpsi(i,x(j))-m(j)*psi(i,x(j))*dpsi(i,m(j)*x(j)))/...
    (psi(i,m(j)*x(j))*dxi(i,x(j))-m(j)*xi(i,x(j))*dpsi(i,m(j)*x(j))));
%% Perhitungan
for iter1 = 1:baris
    total1ext = 0;
    total1sca = 0;
    total2ext = 1E3;
    total2sca = 1E3;
    iter2 = 1;
    while (or((abs(total1ext-total2ext) >= error),abs(total1sca-total2sca) >= error))
        aan(iter2) = an(iter2,iter1);
        % C_ext
        total2ext = total1ext;
        cariext = (2*iter2+1)*real(an(iter2,iter1)+bn(iter2,iter1));
        total1ext = total1ext + cariext;
        % C_sca
        total2sca = total1sca;
        carisca = (2*iter2+1)*(an(iter2,iter1)*conj(an(iter2,iter1))+bn(iter2,iter1)*conj(bn(iter2,iter1)));
        total1sca = total1sca + carisca;
        % Iterasi
        iter2 = iter2 + 1;
    end
    Ce(iter1) = 2*pi*total1ext/(k(iter1))^2;
    Cs(iter1) = 2*pi*total1sca/(k(iter1))^2;
    Ca(iter1) = Ce(iter1) - Cs(iter1);
end
%% Plotting
% C_ext
figure(1);
plot(lambda*1E9,Ce/(pi*a^2),lambda*1E9,Cs/(pi*a^2),lambda*1E9,Ca/(pi*a^2));
legend ('Extinction','Scattering','Absorption');
xlim([bawah atas]);
title (['Cross-Section Gold untuk r = ' num2str(a*1E9) ' nm']);
xlabel ('Panjang gelombang (nm)');
ylabel ('Efisiensi');