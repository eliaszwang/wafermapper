
Z=3;

num_inliers=4*(ones(Z,Z)-eye(Z,Z));

P1=[0 0;
    0 1;
    1 1;
    1 0]';

p2angle=45*pi/180;
p2mat=[cos(p2angle) sin(p2angle);-sin(p2angle) cos(p2angle)];
P2=p2mat*P1+.5;
p3angle=115*pi/180;
p3mat=[cos(p2angle) sin(p2angle);-sin(p2angle) cos(p2angle)];
P3=p3mat*P1-.25;

P=cat(3,P1,P2,P3);

figure(3);
clf;
plot([P1(2,:);P2(2,:)],[P1(1,:);P2(1,:)],'x-');
hold on;
plot([P1(2,:);P3(2,:)],[P1(1,:);P3(1,:)],'o-');
axis image;

Pos1s=cell(Z,Z);
Pos2s=cell(Z,Z);

for i=1:Z
    for j=1:Z
        if i~=j
            
            Pos1s{i,j}=squeeze(P(:,:,i))';
            Pos2s{i,j}=squeeze(P(:,:,j))';
            
        end
    end
end

[dTdpX,dTdpY]=MakeConstraintMatrix(num_inliers,Pos1s,Pos2s,2);
        
%%
       
        lambda=100000;
        Z=size(num_inliers,1);
        C=zeros(3*Z,1);
        C(1:3:3*size(num_inliers,1))=lambda;
        M=(dTdpX + lambda*eye(3*Z,3*Z));
        
        P=C/lambda;
        
       
        delt=.00000001;
        for i=1:100
             dE=M*P-C;
             P=P-delt*dE;
             figure(3);
             clf;
             plot(P(1:3:end));
             pause(.1);
        end
        X=inv(M)*C