function A = init_dictionary(bb,RR)
atomnum=RR*bb^2; % number of atoms in the dictionary
Pn=ceil(sqrt(atomnum));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);
A = DCT(:,1:atomnum );  clear DCT  