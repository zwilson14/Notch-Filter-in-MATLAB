function [X,f]=fft330(x,Fs)
% fft330(x,Fs) produces the Fast Fourier Transform of x(t), given the 
% sampling rate Fs. 
% X(f) is returned in its -fmax to fmax (double-sided) form, and is
% converted to "true" units with the multiplication by dt = 1/Fs, 
% which assumes an energy signal.  
%
% [X,f]=fft330(x,Fs) returns both the X(f) and the frequency vector f.  

% Fourier transform - scaled by dt = 1/Fs, shifted to show +/- freq terms
X = fftshift(fft(x)/Fs);

if nargout == 2
% Generate double-sided freq vector for even or odd N
N = length(x); 
f = ([0:N-1]-floor(N/2))*Fs/N;

% Converting f to column vector if x is a column vector
[m,n]=size(x);
if m > n, f = f'; end 

end

end %function