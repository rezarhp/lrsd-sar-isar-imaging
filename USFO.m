function y = USFO(x,mask,mode)
%   Define forward/adjoint Fourier and Sparse Sampling Operator
%
%   Undersampled Fourier transform operator (USFO) for 2D sequences
%
%   x: input sequence [nx ny]
%   mask: sampling pattern [nx ny] 
%
%   H.R. Hashempour
%   hrhashempour@shirazu.ac.ir
%   Nov 2020

[nx,ny] = size(mask);

if mode == 1      %  Forward operator
    x = reshape(x,[nx ny]);
    u = fft2(x);
    u = fftshift(fftshift(u,1),2);    
    u = u/sqrt(nx*ny);
    y = u(mask>0);
elseif mode == 2  %  Adjoint operator
    u = zeros(nx,ny);
    u(mask>0) = x;
    u = ifftshift(ifftshift(u,1),2);    
    y = sqrt(nx*ny)*ifft2(u);
end
y = y(:);
