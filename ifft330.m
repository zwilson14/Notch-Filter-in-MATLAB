function [x,t]=ifft330(X,Fs)
% ifft330(X,Fs) produces the inverse Fast Fourier Transform of X(f), given
% the sampling rate Fs.  
% x(t) is converted to "true" units with the multiplication by N*df = Fs 
% (the ifft() already includes a 1/N). This assumes an energy signal. 
% Only the real part is returned to avoid imaginary terms due to round-off. 
%
% [x,t]=ifft330(X,Fs) returns both the x(t) and the time vector t.  

% Inverse Fourier transform - unshifting X before applying ifft()
x = real(ifft(ifftshift(X))*Fs);

if nargout == 2
% Generate time vector
N = length(X); 
t = [0:N-1]/Fs;

% Converting t to column vector if X is a column vector
[m,n]=size(X);
if m > n, t = t'; end 

end

end %function