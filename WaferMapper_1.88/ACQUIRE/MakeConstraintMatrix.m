function [dTdpX,dTdpY]=MakeConstraintMatrix(num_inliers,Pos1s,Pos2s,min_inliers)

Z=size(num_inliers,1);

dTdpX = zeros(3*Z,3*Z);
dTdpY = zeros(3*Z,3*Z);

for i=1:Z
    for j=1:Z
        if num_inliers(i,j)>min_inliers
            a1i_ind=(i-1)*3+1;
            a2i_ind=(i-1)*3+2;
            ti_ind=(i-1)*3+3;
            a1j_ind=(j-1)*3+1;
            a2j_ind=(j-1)*3+2;
            
            tj_ind=(j-1)*3+3;
            
            Pos1=Pos1s{i,j};
            Pos2=Pos2s{i,j};
            numK=size(Pos1,1);
            for k=1:numK
                
                v=zeros(1,3*Z);
                
                v(a1i_ind)=Pos1(k,2);
                v(a2i_ind)=Pos1(k,1);
                v(ti_ind)=1;
                v(a1j_ind)=-Pos2(k,2);
                v(a2j_ind)=-Pos2(k,1);
                v(tj_ind)=-1;
                
                dTdpX(a1i_ind,:) = dTdpX(a1i_ind,:) + v*Pos1(k,2);
                dTdpX(a2i_ind,:) = dTdpX(a2i_ind,:) + v*Pos1(k,1);
                dTdpX(ti_ind,:)  = dTdpX(ti_ind,:)  + v;
                
                dTdpX(a1j_ind,:) = dTdpX(a1j_ind,:) - v*Pos2(k,2);
                dTdpX(a2j_ind,:) = dTdpX(a2j_ind,:) - v*Pos2(k,1);
                dTdpX(tj_ind,:)  = dTdpX(tj_ind,:)  - v;
                
                dTdpY(a1i_ind,:) = dTdpY(a1i_ind,:) + v*Pos1(k,2);
                dTdpY(a2i_ind,:) = dTdpY(a2i_ind,:) + v*Pos1(k,1);
                dTdpY(ti_ind,:)  = dTdpY(ti_ind,:)  + v;
                
                dTdpY(a1j_ind,:) = dTdpY(a1j_ind,:) - v*Pos2(k,2);
                dTdpY(a2j_ind,:) = dTdpY(a2j_ind,:) - v*Pos2(k,1);
                dTdpY(tj_ind,:)  = dTdpY(tj_ind,:)  - v;
            end
        end
    end
end


% 
% for i=1:Z
%     for j=1:Z
%         if num_inliers(i,j)>min_inliers
%             a1i_ind=(i-1)*3+1;
%             a2i_ind=(i-1)*3+2;
%             ti_ind=(i-1)*3+3;
%             a1j_ind=(j-1)*3+1;
%             a2j_ind=(j-1)*3+2;
%             
%             tj_ind=(j-1)*3+3;
%             
%             Pos1=Pos1s{i,j};
%             Pos2=Pos2s{i,j};
%             numK=size(Pos1,1);
%             
%             v=zeros(3*Z,numK);
%             
%             v(a1i_ind,:) = Pos1(:,2)';
%             v(a2i_ind,:) = Pos1(:,1)';
%             v(ti_ind,:)  = ones(1,numK);
%             v(a1j_ind,:) = -Pos2(:,2)';
%             v(a2j_ind,:) = -Pos2(:,1)';
%             v(tj_ind,:)  = -ones(1,numK);
%             
%             
%             dTdpX(a1i_ind,:) = dTdpX(a1i_ind,:) + sum(v.*repmat(Pos1(:,2)',3*Z,1),2)';
%             dTdpX(a2i_ind,:) = dTdpX(a2i_ind,:) + sum(v.*repmat(Pos1(:,1)',3*Z,1),2)';
%             dTdpX(ti_ind,:)  = dTdpX(ti_ind,:)  + sum(v,2)';
%             
%             dTdpX(a1j_ind,:) = dTdpX(a1j_ind,:) - sum(v.*repmat(Pos2(:,2)',3*Z,1),2)';
%             dTdpX(a2j_ind,:) = dTdpX(a2j_ind,:) - sum(v.*repmat(Pos2(:,1)',3*Z,1),2)';
%             dTdpX(tj_ind,:)  = dTdpX(tj_ind,:)  - sum(v,2)';
%             
%             dTdpY(a1i_ind,:) = dTdpY(a1i_ind,:) + sum(v.*repmat(Pos1(:,2)',3*Z,1),2)';
%             dTdpY(a2i_ind,:) = dTdpY(a2i_ind,:) + sum(v.*repmat(Pos1(:,1)',3*Z,1),2)';
%             dTdpY(ti_ind,:)  = dTdpY(ti_ind,:)  + sum(v,2)';
%             
%             dTdpY(a1j_ind,:) = dTdpY(a1j_ind,:) - sum(v.*repmat(Pos2(:,2)',3*Z,1),2)';
%             dTdpY(a2j_ind,:) = dTdpY(a2j_ind,:) - sum(v.*repmat(Pos2(:,1)',3*Z,1),2)';
%             dTdpY(tj_ind,:)  = dTdpY(tj_ind,:)  - sum(v,2)';
%             
%         end
%     end
% end


figure(3);
clf;
imagesc(dTdpX);




