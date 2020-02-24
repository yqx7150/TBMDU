%--------------------------------------------------------------
% soft for both real and  complex numbers
%--------------------------------------------------------------
function y = soft(x,T)
%y = sign(x).*max(abs(x)-tau,0);
y = max(abs(x) - T, 0);
y = y./(y+T) .* x;