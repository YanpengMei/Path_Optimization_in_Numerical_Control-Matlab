% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% nonlinear condition for the optimization of iso-parametric segment


function [c,ceq] = A16_nonlcon3(cp_var)

load('A2_parameters.mat','tb_r','cp')

% 1.get the optimized path and tolerance band
path_opt = A4_getSeg(cp_var);
origSeg = A4_getSeg(cp);

% 2. build the equality constraints
temp_1 = 0;
for i = 1:size(path_opt,2)
    temp_2 = sqrt((path_opt(1,i)-origSeg(1,i)).^2+(path_opt(2,i)-origSeg(2,i)).^2+(path_opt(3,i)-origSeg(3,i)).^2);
    if  temp_2 <= tb_r
        temp_1 = temp_1 + 1;
    end
end
c = [];
ceq = temp_1 - size(path_opt,2);
end