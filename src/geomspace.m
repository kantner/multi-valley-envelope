function x = geomspace(a,b,N,alpha)

if alpha == 1.0
    
  x = linspace(a,b,N);
    
else
    
  S  = ( 1-alpha^(N-1) )/(1-alpha);
  h1 = (b-a)/S;
  
  h  = h1 * alpha.^(0 : N-3);
  
  x(1)     = a;
  x(2:N-1) = a + cumsum(h);
  x(N)     = b;

end