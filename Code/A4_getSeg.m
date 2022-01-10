% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% Get the segments
% If we want to get the original segments, the input should be "cp"
% for non-iso-segment, we will make the segment into 2D

function path = A4_getSeg(cp_varA4)

load('A2_parameters.mat','knots','cp','step_xy')
load('A2_TupleSeg.mat','TupleSeg')
surface = nrbmak(cp_varA4,knots);

max_delta_cpX = max(cp(1,end,:) - cp(1,1,:));
max_delta_cpY = max(cp(2,:,end) - cp(2,:,1));
delta_u = (TupleSeg(2,1) - TupleSeg(1,1));
delta_v = (TupleSeg(2,2) - TupleSeg(1,2));
numP_P = ceil(sqrt(((max_delta_cpX*delta_u/step_xy).^2 + (max_delta_cpY*delta_v/step_xy).^2)));

% segType        segment type,"1" for constant v, "2" for constant u, "3" for non-iso...
if delta_v==0
        segType=1;
elseif delta_u==0
        segType=2;
else
        segType=3;
end

% ut_path             discrete u for the path
% vt_path             discrete v for the path

% we can devided the path generation into 3 cases:
% case 1: iso-parametric segment with constant v
% case 2: iso-parametric segment with constant u
% case 3: non-iso-parametric segment

if segType==1       % u-segment with constant v
        ut_path = linspace(TupleSeg(1,1),TupleSeg(2,1),numP_P);
        vt_path = zeros(1,numP_P) + TupleSeg(1,2);  
        temp1 = nrbeval(surface,{ut_path,vt_path});
        path = temp1(:,:,1);
        
elseif segType==2       % v-segment with constant u
        ut_path = zeros(1,numP_P) + TupleSeg(1,1);
        vt_path = linspace(TupleSeg(1,2),TupleSeg(2,2),numP_P);
        temp2 = nrbeval(surface,{ut_path,vt_path});
        temp3 = permute(temp2,[1 3 2]);
        path = temp3(:,:,1);

else % segType==3, non-iso-parametric segment
       ut_path = linspace(TupleSeg(1,1),TupleSeg(2,1),numP_P);
       vt_path = linspace(TupleSeg(1,2),TupleSeg(2,2),numP_P);
       temp4 = zeros(3,numP_P);
       for j = 1:numP_P
           temp4(:,j) = nrbeval(surface,{ut_path(j),vt_path(j)});
       end
       
       temp5 = temp4(:,:,1);
       
       % for non-iso case, we want the path to be in (xx,yy)-coordinate
       % That is to say, we convert 3D-points into 2D
       
       xx = sqrt((temp5(1,:)-temp5(1,1)).^2 + (temp5(2,:)-temp5(2,1)).^2);
       yy = temp5(3,:)-temp5(3,1);
       
       path = [xx;yy];      
end
% end of the function main body
end