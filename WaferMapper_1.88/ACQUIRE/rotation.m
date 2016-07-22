function [f,df]=rotation(R,A,X)
% MSE for rotation matrix, R
% X=R*A
% X: actual x,y stigmation values (2xn)
% R: rotation matrix, [a b c d] (1x4) => (2x2) [a b; c d]
% A: MAPFoSt estimated aon,adiag (2xn)
    %b=R(5);
    R=[R(1) R(2);R(3) R(4)]; %convert vector to 2x2 matrix
    f=sum(sum((X-(R*A)).^2));
%     df(1)=sum(sum(-2*X(1,:).*A(1,:)+2*R(1,1)*A(1,:).^2+2*R(1,2)*A(1,:).*A(2,:)));
%     df(2)=sum(sum(-2*X(1,:).*A(2,:)+2*R(1,2)*A(2,:).^2+2*R(1,1)*A(1,:).*A(2,:)));
%     df(3)=sum(sum(-2*X(2,:).*A(1,:)+2*R(2,1)*A(1,:).^2+2*R(2,2)*A(1,:).*A(2,:)));
%     df(4)=sum(sum(-2*X(2,:).*A(2,:)+2*R(2,2)*A(2,:).^2+2*R(2,1)*A(1,:).*A(2,:)));
    df=-2*(X-(R*A))*A';
    df=[df(1,:) df(2,:)]; %convert matrix to vector
%     temp=-2*(X-(R*A+b))*ones(size(A,2),2);
%     df=[df temp(2)];
end