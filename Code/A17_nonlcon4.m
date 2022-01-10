% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% non-linear condition for the optimization with linearly defined tolerance band


function [c,ceq] = A17_nonlcon4(var)
load('A2_parameters.mat','cp','knot')
load('A2_fittedCP_best.mat','fittedCP_best')
cp_varAC2 = [fittedCP_best(1,:) ; var];

% 1. get 2D real path
[~,~,xx,yy] = A6_fit(cp);
realPath2D = [xx ; yy];

numP_P = size(realPath2D,2);


% 2. get 2D optimized path
curve = nrbmak(cp_varAC2,knot);
ut_path = linspace(0,1,numP_P);
temp1 = nrbeval(curve,ut_path);
temp2 = temp1(:,:,1);
path_opt = temp2(1:2,:) - temp2(1:2,1);


% 3. get the information from "A13_linearTB"
[disLTB2RP] = A13_linearTB(cp);


% 4. build the equality constraints
temp3 = 0;
for i = 1:numP_P
    
    if i ==1
        temp4 = sqrt((realPath2D(1,i) - path_opt(1,i))^2 + (realPath2D(2,i) - path_opt(2,i))^2);
        if temp4 <= disLTB2RP(i)
            temp3 = temp3 + 1;
        end
    elseif i >1 && i<numP_P
        temp5 = sqrt((realPath2D(1,i) - path_opt(1,i))^2 + (realPath2D(2,i) - path_opt(2,i))^2);
        if temp5 <= disLTB2RP(i-1) && temp5 <= disLTB2RP(i)
            temp3 = temp3 + 1;
        end
    else
        temp6 = sqrt((realPath2D(1,i) - path_opt(1,i))^2 + (realPath2D(2,i) - path_opt(2,i))^2);
        if temp6 <= disLTB2RP(i-1)
            temp3 = temp3 + 1;
        end
        
    end
end
c = [];
ceq = temp3 - numP_P;
end