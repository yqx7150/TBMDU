function Dnorm = get_operator_norm_Liu(D,disp_flag)

%       Dnorm = get_operator_norm(P,disp_flag)
%
% returns the norm of an operator defined in a SPARCO problem P
%
% if  disp_flag ~= 0, the norm is displayed at each iteration
% 

nbIterMin=30;
nbIterMax=3000;
acc = 1e-10;
n = size(D,2);

it=1;
Dnorm_old= 0;
v=randn(n,1);
Dnorm = sqrt(norm(v));
while  ( ((abs ( Dnorm_old - Dnorm ) > acc) || ( it<nbIterMin ) ) && it <nbIterMax )    ;
     v=v/norm(v);
     u = D*v; % u=D*v;
     v = D'*u; % v=Dt*u;
     if disp_flag ~= 0 ; 
         disp( [it , sqrt(norm(v))] ); 
     end;
     it=it+1;
     Dnorm_old = Dnorm;
     Dnorm=sqrt(norm(v));
end;

Dnorm=sqrt(norm(v)); % approximation of the operator norm
