% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% This function is aimed to change 2D segment back into 3D segment

function seg3D = A8_2Dinto3D(cp_2D)
load('A2_parameters.mat','cp','knots','step_xy','knot')
load('A2_TupleSeg.mat','TupleSeg')

% 1. get the fitted segment
% FCur             fitted NURBS-curve
FCur = nrbmak(cp_2D,knot);

max_delta_cpX = max(cp(1,end,:) - cp(1,1,:));
max_delta_cpY = max(cp(2,:,end) - cp(2,:,1));
delta_u = (TupleSeg(2,1) - TupleSeg(1,1));
delta_v = (TupleSeg(2,2) - TupleSeg(1,2));
numP_P = ceil(sqrt(((max_delta_cpX*delta_u/step_xy).^2 + (max_delta_cpY*delta_v/step_xy).^2)));
ut_seg = linspace(0,1,numP_P);

temp1 = nrbeval(FCur,ut_seg);
FSeg = temp1(1:2,:);


% 2. get xx and yy
xx = FSeg(1,:)-FSeg(1,1);
yy = FSeg(2,:)-FSeg(2,1);


% 3. get real segment
surface = nrbmak(cp,knots);
ut_path = linspace(TupleSeg(1,1),TupleSeg(2,1),numP_P);
vt_path = linspace(TupleSeg(1,2),TupleSeg(2,2),numP_P);
temp2 = zeros(3,numP_P);
for j = 1:numP_P
    temp2(:,j) = nrbeval(surface,{ut_path(j),vt_path(j)});
end
realSeg = temp2(:,:,1);

% information of coordinate transformation
deltaXseg = realSeg(1,end) - realSeg(1,1);
deltaYseg = realSeg(2,end) - realSeg(2,1);


% 4. segment type
% segType        segment type,"1" for constant v, "2" for constant u, "3" for non-iso...
if delta_v==0
        segType=1;
elseif delta_u==0
        segType=2;
else
        segType=3;
end


% 5. 2D to 3D
if segType == 1
    x = xx + realSeg(1,1);
    y = zeros(1,numP_P) + realSeg(2,1);
    z = yy + realSeg(3,1);
    
    seg3D = [x;y;z];
elseif segType == 2
    x = zeros(1,numP_P) + realSeg(1,1);
    y = xx + realSeg(2,1);
    z = yy + realSeg(3,1);
    
    seg3D = [x;y;z];
else
    cos_theta = deltaXseg/sqrt(deltaXseg^2+deltaYseg^2);
    sin_theta = deltaYseg/sqrt(deltaXseg^2+deltaYseg^2);
    
    x = xx*cos_theta + realSeg(1,1);
    y = xx*sin_theta + realSeg(2,1);
    z = yy+ realSeg(3,1);
    
    seg3D = [x;y;z];
end


end