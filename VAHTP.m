function [X] = VAHTP(Y,A, MaxIter, tol)
%------------------------------------------------------------------------
%%Input:
% Y:          observed data
% A:          measurement matrix
% MaxIter:  maximum number of iteration
% tol:        the stop tolerance

%%Output:
% X: recovered solution
%--------------- ---------------------------------------------------------

[MM,N]=size(A);
L=size(Y,2);

X = zeros(N,L);
res = Y-A*X;
NbIter=1;

while (NbIter <= MaxIter) && (norm(res) > tol)
  U = X+A'*res;
  V_v=sum(U.^2,2)/L-(sum(U,2)/L).^2;
  [~,sorted_idx]=sort(V_v,'descend');
  index=min(NbIter,MM);
  Snew=sort(sorted_idx(1:index));
  X=zeros(N,L);
  X(Snew,:)=A(:,Snew)\Y;
  res = Y-A*X;
  NbIter = NbIter+1;
end

end

