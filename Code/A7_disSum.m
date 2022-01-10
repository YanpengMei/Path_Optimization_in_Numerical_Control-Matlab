% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% calculate the sum of the distance between the fitted NURBS-curve and discreted points


function disSum = A7_disSum(var)
load('A2_parameters.mat','cp','step_xy','knot')
[fittedCP_init,~,xx,yy] = A6_fit(cp);

load('A2_TupleSeg.mat','TupleSeg')

%% 1. make a initial fitted NURBS-curve
% FCur             fitted NURBS-curve
cp_var = [fittedCP_init(:,1)  var  fittedCP_init(:,end)];
FCur = nrbmak(cp_var,knot);        



%% 2. get "ut_seg"
max_delta_cpX = max(cp(1,end,:) - cp(1,1,:));
max_delta_cpY = max(cp(2,:,end) - cp(2,:,1));
delta_u = (TupleSeg(2,1) - TupleSeg(1,1));
delta_v = (TupleSeg(2,2) - TupleSeg(1,2));
numP_P = ceil(sqrt((max_delta_cpX*delta_u/step_xy).^2 + (max_delta_cpY*delta_v/step_xy).^2));
ut_seg = linspace(0,1,numP_P);



%% 3. get the sum of the distance 
% Fseg             discreted fitted segment
temp5 = nrbeval(FCur,ut_seg);
Fseg = temp5(1:2,:);

% disSum             sum of the distance
disSum = sum(sqrt((xx(1,:)-Fseg(1,:)).^2+(yy(1,:)-Fseg(2,:)).^2));



end