function [f,F]=getSpectrum(t,dt,u)
% Determine discrete frequency to use
M = length(t); % Length of time series
fN = 1./(2.*dt); % Nyquist frequency (Hz)
df = fN./floor(M./2); % Discrete frequency step (Hz)
% Compute spectrum
Y = fft(u)./M; % Discrete Fourier transform (m/s)
A = abs([Y(1) 2.*Y(2:floor(M/2)) Y(floor(M/2)+1)]); % Harmonic amplitudes (m/s)
f = [0:df:fN]; % Vector of positive frequencies (Hz)
F = 0.25.*A.^2./df; % Energy density (m^2/s^2/Hz = m^2/s)
% The factor 0.25 is used, since if u'=A*sin(t), k_u=0.5.*mean(u.^2)=0.25*A.^2
% Check against the turbulent kinetic energy from the velocity component (these should match!)
%~ ku_ts = 0.5.*mean(u.^2) % Calculate directly from time series (m^2/s^2)
%~ ku_sp = trapz(f,F) % Integral of energy density spectrum (m^2/s^2)
end
