%% Produz os atrasos de grupo do filtro de Bessel
close all
clear
clc

%% Parâmetros
N = 1:4;
w = logspace(-2, 2, 100).';
taug = zeros(length(w), length(N));

%% Iterações
for n=1:length(N)
	%% Filtro
	[z, p, k] = besselap(N(n));
	[b, a] = zp2tf(z, p, k);
	H = freqs(b, a, w);

	%% Atraso de grupo
	phaH = unwrap(atan2(imag(H), real(H)));
	taug(2:end, n) = -diff(phaH)./diff(w);
	taug(1, n) = taug(2, n);
	%% Renormalização do besselap()
	taug(:, n) = taug(:, n)/taug(1, n);
end

%% Plotagem
semilogx(w, taug(:, 1))
hold on
semilogx(w, taug(:, 2))
semilogx(w, taug(:, 3))
semilogx(w, taug(:, 4))
hold off
grid on

%% Exportar dados para Gnuplot
writematrix([w, taug], 'taugbessel.dat', ...
	'Delimiter', 'tab')