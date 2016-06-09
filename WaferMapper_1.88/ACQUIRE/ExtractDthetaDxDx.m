function [dtheta,dx,dy]=ExtractDthetaDxDx(transforms)

Z=size(transforms,1);

dtheta=zeros(Z,Z);
dx=zeros(Z,Z);
dy=zeros(Z,Z);

for i=1:Z
    for j=1:Z
        if ~isempty(transforms{i,j})
            trans=transforms{i,j};
            M=trans.tdata.T;
            dtheta(i,j)=atan2(M(1,2),M(1,1));
            dx(i,j)=M(2,2);
            dy(i,j)=M(3,2);

        end
    end
end

