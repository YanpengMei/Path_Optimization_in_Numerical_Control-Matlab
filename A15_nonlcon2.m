% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% non-linear condition for the optimization of the segment in 2D

function [c,ceq] = A15_nonlcon2(var)
load('A2_parameters.mat','tb_r','cp','step_xy','knot')
load('A2_TupleSeg.mat','TupleSeg')
load('A2_fittedCP_best.mat','fittedCP_best')

% 1. get the fitted segment
% to make sure the programm works, we will only optimize "y" of the control points
cp_var = [fittedCP_best(1,:) ; var];

% FCur             fitted NURBS-curve
FCur = nrbmak(cp_var,knot);  


% 2. get the discreted segment
max_delta_cpX = max(cp(1,end,:) - cp(1,1,:));
max_delta_cpY = max(cp(2,:,end) - cp(2,:,1));
delta_u = (TupleSeg(2,1) - TupleSeg(1,1));
delta_v = (TupleSeg(2,2) - TupleSeg(1,2));
numP_P = ceil(sqrt(((max_delta_cpX*delta_u/step_xy).^2 + (max_delta_cpY*delta_v/step_xy).^2)));
ut_seg = linspace(0,1,numP_P);

% Fseg             discreted fitted segment
temp1 = nrbeval(FCur,ut_seg);
Fseg = temp1(1:2,:) - temp1(1:2,1);


% 3.get the tolerance band
[~,~,xx,yy] = A6_fit(cp);
realSeg = [xx ; yy];


% 4. build the equality constraint
temp3 = 0;
for i = 1:numP_P
    temp4 = sqrt((Fseg(1,i)-realSeg(1,i))^2 + (Fseg(2,i)-realSeg(2,i))^2);
    if  temp4 <= tb_r
        temp3 = temp3 + 1;
    end
end
c = [];
ceq = temp3 - numP_P;
end