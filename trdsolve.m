function b=trdsolve(A1, B1, C1, b)
%int:i,n,m
%matrix:A,B,C,A1,B1,C1,b

n=size(A1, 3);
m=size(A1, 1);

% en matlab la copia de A1, B1, C1 -> A, B, C
% no hace falta, pero si en c++, ya que lutrb modifica A, B, C
A=zeros(m, m, n);
B=zeros(m, m, n);
C=zeros(m, m, n);
A=A1;
B=B1;
C=C1;
[A, B, C]=lutrb(A, B, C);

b(:,1)=fwdsolve(A(:, :, 1), b(:,1));
for i=2:n
  b(:, i)=maxpy(B(:,:,i-1), b(:,i-1), b(:,i));
  b(:,i)=fwdsolve(A(:, :, i), b(:,i));
end

b(:,n)=bwdsolve(A(:, :, n), b(:,n));

for i=n-1:-1:1
  b(:, i)=maxpy(C(:,:,i), b(:,i+1), b(:,i));
  b(:,i)=bwdsolve(A(:, :, i), b(:,i));
end

function [A, B, C]=lutrb(A, B, C)
%int:i,n
%matrix:A,B,C

%fprintf('lutrb\n');
n=size(A, 3);

for i=1:n-1
  A(:,:,i)=nplu(A(:,:,i));
  B(:,:,i)=utrsolve(A(:,:,i), B(:,:,i));
  C(:,:,i)=fwdsolve(A(:,:,i), C(:,:,i));
  A(:,:,i+1)=maxpy(B(:,:,i), C(:,:,i),A(:,:,i+1));
end

A(:,:,n)=nplu(A(:,:,n));


function A=nplu(A)
%int:i,j,k,n
%matrix:A

n=size(A, 1);
for k=1:n
  if ( abs(A(k,k)) < 1e-15 )
    error('Pivot too small');
  end
  
  for i=k+1:n
    A(i, k)=A(i, k)/A(k, k);
    for j=k+1:n
      A(i, j)=A(i, j)-A(i, k)*A(k, j);
    end
  end
end

function B=fwdsolve(A, B)
%int:i,j,k,s,n
%matrix:A,B

s=size(B, 2);
n=size(A, 1);
for k=1:s
  for i=1:n
    for j=1:i-1
      B(i, k)=B(i, k)-A(i, j)*B(j, k);
    end
  end
end


function B=bwdsolve(A, B)
%int:i,j,k,s,n
%matrix:A,B

s=size(B, 2);
n=size(A, 1);
for k=1:s
  for i=n:-1:1
    for j=i+1:n
      B(i, k)=B(i, k)-A(i, j)*B(j, k);
    end
    B(i, k)=B(i, k)/A(i, i);
  end
end


function B=utrsolve(A, B)
%int:i,j,k,s,n
%matrix:A,B

s=size(B, 2);
n=size(A, 1);
for k=1:s
  for i=1:n
    for j=1:i-1
      B(k, i)=B(k, i)-A(j, i)*B(k, j);
    end
    B(k, i)=B(k, i)/A(i, i);
  end
end


function y=maxpy(A, x, y)
%int:i,j,k,m,l,n
%matrix:A,x,y

m=size(y, 1);
l=size(y, 2);
n=size(A, 2);

for k=1:l
  for i=1:m
    for j=1:n
      y(i, k)=y(i, k)-A(i, j)*x(j, k);
    end
  end
end



