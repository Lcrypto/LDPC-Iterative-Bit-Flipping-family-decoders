function [out] = ldpc_encode(in,G,q) 
% function [out] = ldpc_encode(in,H,q) 
% encodes data from "in" using G over GFq 
% q = 2,4,8,16,32,64,128 or 256 
% Please, make sure that the data you use is valid! 
% Requires matlab communication toolbox for GFq operations. 
% Basically, this function performs 'in*G' over GFq.  
% this function is slow, will write a C version some time 
 
%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu 
%   $Revision: 1.0 $  $Date: 1999/11/23 $ 
 
[k,n] = size(G);  
if q==2 % binary case 
   out = mod((in')*G,2); 
else 
   M=log2(q); % GFq exponent 
   [tuple power] = gftuple([-1:2^M-2]', M, 2);  
   alpha = tuple * 2.^[0 : M - 1]'; 
   beta(alpha + 1) = 0 : 2^M - 1; 
   ll = ones(1,n)*(-Inf);      % will store results here, initialize with all zeros (exp form) 
   kk=zeros(1,n);
   for i=1:k % multiply each row of G by the input symbol in GFq 
      ii = power(beta(in(i)+1)+1); % get expon. representation of in(i) 
      jj = power(beta(G(i,:)+1)+1);% same for the row of G 
        for j=1:n
       %parfor j=1:n
        kk(j) = gfmul(ii,jj(j,1),tuple); % this is exponential representation of the product 
       end
      ll = gfadd(ll,kk,tuple); 
   end 
    
   out=zeros(size(ll)); 
   nzindx = find(isfinite(ll)); 
   out(nzindx) = alpha(ll(nzindx)+2); 
   out = out(:); 
end 

