% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% non-linear condition for NURBS-curve-fitting


function [c,ceq] = A14_nonlcon1(var)
load('A2_parameters.mat','cp','tb_r','knot')

% 1. get real segment in 2D
[fittedCP_init,~,xx,yy] = A6_fit(cp);
realSeg2D = [xx ; yy];


% 2. get fitted segment in 2D
cp_var = [fittedCP_init(:,1)  var  fittedCP_init(:,end)];
curve = nrbmak(cp_var,knot);
numP_P = size(realSeg2D,2);
ut_seg = linspace(0,1,numP_P);
temp1 = nrbeval(curve,ut_seg);
temp2 = temp1(:,:,1);
seg2D = temp2(1:2,:);


% 3. build the equality constraints
temp3 = 0;
for i = 1:numP_P
    temp4 = sqrt((seg2D(1,i)-realSeg2D(1,i)).^2+(seg2D(2,i)-realSeg2D(2,i)).^2);
    if  temp4 <= tb_r
        temp3 = temp3 + 1;
    end
end
c = [];
ceq = temp3 - size(seg2D,2);


end