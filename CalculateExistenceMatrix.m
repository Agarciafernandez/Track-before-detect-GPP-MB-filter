function A=CalculateExistenceMatrix(n)

%Te calcula todas las posibles combinaciones de exitence para n (numero
%maximo de blancos)

%All possible combinations of existence variables for a maximum target of
%targets n


A=ones(2^n,n);
v=(0:2^n-1).';
for ii=n:-1:1
    A(:,ii)=bitget(v,ii);
end

A=A';