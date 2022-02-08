clear
close all
%% 0.Informations
% Author: Yanpeng Mei
% Last updated in: 20. Okt. 2021
% Matlab Version: R2020b

% This is the main file of my thesis "Path optimization when interpolating in Numerical Controls"
% In my thesis, the whole tool path will be divided into many segments based on the Tuple
% There are two kinds of segments: iso-parametric segment(when 'u' or 'v' is constant) and non-iso-parametric segment
% In short: iso-segment and non-iso-segment
% The toolpath I am considering is a combination of many NURBS-curves(segments)
% 'Optimization' means minimizing the curvature of each NURBS-segment
% For non-iso-segment, I will first translate the path in 3D space in 2D plane
% For iso-segment, I can also translate the path in 3D space in 2D plane, or optimize the curve in 3D space directly

% List of functions:
% 'A3_origPath' get the original tool path
% 'A4_getSeg' get the segment with given Tuple
% 'A5_Seg2Path' combines the segments into a path
% 'A6_fit' fits the discreted points into a NURBS-curve
% 'A7_disSum' calculates the sum of the distance between the fitted NURBS-curve and discreted points
% 'A8_2Dinto3D' changes 2D segment into 3D
% 'A9_curIso' calculates the curvature of an iso-segment at each sampling point
% 'A10_sumCurIso' calculates the sum of the curvature of an iso-segment at all sampling points
% 'A11_cur2D' calculates the curvature of a 2D-segment at each sampling point
% 'A12_sumCur2D' calculates the sum of the curvature of a 2D-segment at all sampling points
% 'A13_linearTB' generates a linear tolerance band
% 'A14_nonlcon1' is the non-linear condition for NURBS-curve-fitting
% 'A15_nonlcon2' is the non-linear condition for the optimization of the 2D-segment
% 'A16_nonlcon3' is the non-linear condition for the optimization of the iso-segment
% 'A17_nonlcon4' is the non-linear condition for the optimization of the segment with linear tolerance band



%% 1.Create a NURBS-surface as the test surface
% in the following, "test surface/path" will be also called as "original surface/path" or "real surface/path"
% 1.1 define control points
% allocate multi-dimensional array of control points for NURBS-surface
% cp               control points
cp = zeros(3,6,3);
cp(:,:,1) = [ 0  50  150  250  350  400 ;     % x-axis coordinate
              0  0   0    0    0    0   ;     % y-axis coordinate
             -5  10 -20   20  -10   5   ];    % z-axis coordinate

cp(:,:,2) = [ 0   50   150  250  350  400;
              50  50   50   50   50   50 ;
             -50  100 -200  200 -100  50];

cp(:,:,3) = [ 0   50   150  250  350 400 ;
              100 100  100  100  100 100 ;
             -5   10  -20   20  -10  5 ];


% 1.2 define knots
knots{1} = [0 0 0 1/4 2/4 3/4 1 1 1]; % knots along u
knots{2} = [0 0 0 1 1 1]; % knots along v


% 1.3 generate and plot the test NURBS-surface
% surf          test NURBS-surface
figure
set(gcf,'outerposition',get(0,'screensize'));
surf = nrbmak(cp,knots);
nrbplot(surf,[100,50]);

xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
zlabel('z','Interpreter','latex')
hold on



%% 2.Input Tuple(u,v), set parameters and get the test path 
% 2.1 input Tuple
% Tuple               Tuple is a n*2 matrix [u_i  v_i], n is the number of points in Tuple
% input Tuple with right order

% 1st, test with the simplest iso-parametric u-segment
% Tuple = [0         0.5;
%          1         0.5];

% 2nd, test with the simplest iso-parametric v-segment
% Tuple = [0.625         0;
%          0.625         1];

% 3rd, test with the simplest non-iso-parametric segment
% Tuple = [0       0;
%          1       1];

% finally, test with commen case, for example, 2 iso and 1 non-iso segments
% Tuple = [0         0.1;
%          0.625     0.1;
%          0.625     0.9;
%          0         0.5 ];

% All the above cases have been tested. For the sake of demonstration, the first case will be shown next.
% For the demonstration of other cases, section 5 should be changed
Tuple = [0         0.5;
         1         0.5];

     
% 2.2 find out the type of each segment
% numSeg              number of segments in the whole path
numSeg = size(Tuple,1)-1;

% delta_uv            difference of the Tuple in u- and v-direction
delta_uv = [(Tuple(2:end,1) - Tuple(1:(end-1),1))  (Tuple(2:end,2) - Tuple(1:(end-1),2))];

% segType        segment type,"1" for constant v, "2" for constant u, "3" for non-iso...
segType = zeros(numSeg,2);
segType(:,1) = 1:1:numSeg;
for i=1:numSeg
    if delta_uv(i,2)==0
        segType(i,2)=1;
    elseif delta_uv(i,1)==0
        segType(i,2)=2;
    else
        segType(i,2)=3;
    end
end


% 2.3 define the radius of the rolerance band and step size when generating the path
% tb_r               radius of the rolerance band
% step_xy            step size in 3D-coordinate system when generating the path
tb_r = 4;
step_xy = 2;


% 2.4 define the knot for fitting the discreted points 2D NURBS-curve
knot = [0 0 0 0 1/11 2/11 3/11 4/11 5/11 6/11 7/11 8/11 9/11 10/11 1 1 1 1];


% 2.5 save the parameters to a ".mat"
save A2_parameters.mat


% 2.6 generate and plot the original path
path = A3_origPath(cp);     
set(gca,'LooseInset',[0,0,0,0]);
plot3(path(1,:),path(2,:),path(3,:),'r','LineWidth',5,'DisplayName','Testpfad')
hold off
legend('test surface','test path')
title('Test NURBS-surface and test tool path')



%% 3.Prepare for the minimization of the path-curvature
% each segment will be optimized separately
% 3.1 set the optimization options
options1 = optimoptions('fmincon','MaxFunctionEvaluations',1000000,'MaxIterations',10000,'OptimalityTolerance',1e-4,'StepTolerance',1e-5,'ConstraintTolerance',1e-5);


% 3.2 create empty cells and arrays
% for the original path
origCur = cell(1,numSeg);
origCurSum = zeros(1,numSeg);
origCurMax = zeros(1,numSeg);

% "1" stands for optimization method 1, which means the optimization will be done by adjusting the control points of the test surface
optSeg1 = cell(1,numSeg);
optCur1 = cell(1,numSeg);
optCurSum1 = zeros(1,numSeg);
optCurMax1 = zeros(1,numSeg);

% "2" stands for optimization method 2, which means the optimization will be done by adjusting the control points of the fitted 2D NURBS-curve
optSeg2 = cell(1,numSeg);
optCur2 = cell(1,numSeg);
optCurSum2 = zeros(1,numSeg);
optCurMax2 = zeros(1,numSeg);

% "fit" means that these are the results of the fitted NURBS-curve
fitSeg = cell(1,numSeg);
fitCur = cell(1,numSeg);
fitCurSum = zeros(1,numSeg);
fitCurMax = zeros(1,numSeg);

% "wlTB" means that the optimization is done by using the linearly defined tolerance band
optSeg_wlTB = cell(1,numSeg);
optCur_wlTB = cell(1,numSeg);
optCurSum_wlTB = zeros(1,numSeg);
optCurMax_wlTB = zeros(1,numSeg);



%% 4. Optimization
for i = 1:numSeg
    TupleSeg = Tuple(i:i+1,:);
    save('A2_TupleSeg','TupleSeg')
    
    if segType(i,2)== 1 || segType(i,2)== 2    % optimize type 1 and 2(iso-parametric segments)
        
        var0_iso = cp;
        cp_opt = fmincon(@A10_sumCurIso,var0_iso,[],[],[],[],[],[],@A16_nonlcon3,options1);
        
        % original
        origCur{i} = A9_curIso(cp);
        origCurSum(i) = A10_sumCurIso(cp);
        origCurMax(i) = max(A9_curIso(cp));
        
        % optimize with method 1
        optSeg1{i} = A4_getSeg(cp_opt);
        optCur1{i} = A9_curIso(cp_opt);
        optCurSum1(i) = A10_sumCurIso(cp_opt);
        optCurMax1 = max(A9_curIso(cp_opt));
        
        % fit the discreted points into a 2D NURBS-curve
        [fittedCP_init,numCP,~,~] = A6_fit(cp);
        fittedVAR_best = fmincon(@A7_disSum,fittedCP_init(:,2:(numCP-1)),[],[],[],[],[],[],@A14_nonlcon1,options1);
        fittedCP_best = [fittedCP_init(:,1)  fittedVAR_best  fittedCP_init(:,end)];
        save('A2_fittedCP_best.mat','fittedCP_best')
        
        % This part will only be used to prove the reasonableness of the curve fit
        fitSeg{i} = A8_2Dinto3D(fittedCP_best);
        fitCurSum(i) = A12_sumCur2D(fittedCP_best(2,:));
        fitCur{i} = A11_cur2D(fittedCP_best(2,:));
        fitCurMax = max(A11_cur2D(fittedCP_best(2,:)));

        % optimize with method 2
        var_best = fmincon(@A12_sumCur2D,fittedCP_best(2,:),[],[],[],[],[],[],@A15_nonlcon2,options1);
        cp_best = [fittedCP_best(1,:) ; var_best];
        optSeg2{i} = A8_2Dinto3D(cp_best);
        optCurSum2(i) = A12_sumCur2D(var_best);
        optCur2{i} = A11_cur2D(var_best);
        optCurMax2 = max(A11_cur2D(var_best));
               
        % optimize with linearly defined tolerance band
        VARopt_wlTB = fmincon(@A12_sumCur2D,fittedCP_best(2,:),[],[],[],[],[],[],@A17_nonlcon4,options1);
        CPopt_wlTB = [fittedCP_best(1,:) ; VARopt_wlTB];
        optSeg_wlTB{i} = A8_2Dinto3D(CPopt_wlTB);
        optCurSum_wlTB(i) = A12_sumCur2D(VARopt_wlTB);
        optCur_wlTB{i} = A11_cur2D(VARopt_wlTB);
        optCurMax_wlTB = max(A11_cur2D(VARopt_wlTB));
        
    else    % optimize type 3(non-iso-parametric segments)
        % fit the discreted points into a 2D NURBS-curve
        [cp0_Niso,numCP,~,~] = A6_fit(cp);
        fittedVAR_best = fmincon(@A7_disSum,cp0_Niso(:,2:(end-1)),[],[],[],[],[],[],@A14_nonlcon1,options1);
        fittedCP_best = [cp0_Niso(:,1)  fittedVAR_best  cp0_Niso(:,end)];
        save('A2_fittedCP_best.mat','fittedCP_best')
        % for non-iso-parametric segment, we name the optimized control points cp_best
        var_best = fmincon(@A12_sumCur2D,fittedCP_best(2,:),[],[],[],[],[],[],@A15_nonlcon2,options1);
        cp_best = [fittedCP_best(1,:) ; var_best];
        
        % original 
        origCurSum(i) = A12_sumCur2D(fittedCP_best(2,:));
        origCur{i} = A11_cur2D(fittedCP_best(2,:));
        origCurMax(i) = max(A11_cur2D(fittedCP_best(2,:)));
        
        % optimize with method 2
        optSeg2{i} = A8_2Dinto3D(cp_best);
        optCurSum2(i) = A12_sumCur2D(var_best);
        optCur2{i} = A11_cur2D(var_best);
        optCurMax2 = max(A11_cur2D(var_best));
        
        % optimize with linearly defined tolerance band
        VARopt_wlTB = fmincon(@A12_sumCur2D,fittedCP_best(2,:),[],[],[],[],[],[],@A17_nonlcon4,options1);
        CPopt_wlTB = [fittedCP_best(1,:) ; VARopt_wlTB];
        optSeg_wlTB{i} = A8_2Dinto3D(CPopt_wlTB);
        optCurSum_wlTB(i) = A12_sumCur2D(VARopt_wlTB);
        optCur_wlTB{i} = A11_cur2D(VARopt_wlTB);
        optCurMax_wlTB = max(A11_cur2D(VARopt_wlTB));
         
    end
end



%% 5. Plot
% this section works only for the Tuple we choose in section 2.1
% for other situations, the plot-programms should be changed

% 5.1 compare the results of "path", "path_opt1", "path_wlTB" and  "path_opt2"
% sum of the curvature
sumCur_orig = sum(origCurSum);
sumCur_opt1 = sum(optCurSum1);
sumCur_opt2 = sum(optCurSum2);
sumCur_wlTB = sum(optCurSum_wlTB);
sumCur_fit = sum(fitCurSum);

% maximum of the curvature
maxCur_orig = max(origCurMax);
maxCur_opt1 = max(optCurMax1);
maxCur_opt2 = max(optCurMax2);
maxCur_wlTB = max(optCurMax_wlTB);
maxCur_fit = max(fitCurMax);


% 5.2 first, plot all the curve in 2D
% original path
figure
set(gcf,'outerposition',get(0,'screensize'));
% set(gca,'LooseInset',[0,0,0,0]);
plot(path(1,:),path(3,:),'r','LineWidth',1)
hold on

% optimized path with method 1
path_opt1 = A5_Seg2Path(optSeg1);
plot(path_opt1(1,:),path_opt1(3,:),'g','LineWidth',1)
hold on

% optimized path with method 2
path_opt2 = A5_Seg2Path(optSeg2);
plot(path_opt2(1,:),path_opt2(3,:),'b','LineWidth',1)
hold on

% fitted path
path_fit = A5_Seg2Path(fitSeg);
plot(path_fit(1,:),path_fit(3,:),'c','LineWidth',1)
hold on

% optimized path with linearly defined tolerance band
path_wlTB = A5_Seg2Path(optSeg_wlTB);
plot(path_wlTB(1,:),path_wlTB(3,:),'k','LineWidth',1)
hold off

xlabel('x','Interpreter','latex') 
ylabel('z','Interpreter','latex') 
legend('path_{orig}','path_{opt1}','path_{opt2}','path_{fit}','path_{wlTB}')
title('Plot of all the paths')


% 5.3 second, show that the NURBS-curve is well fitted
% distance           distance between the original path and the fitted path at each sampling point   
distance = sign(path(3,:)-path_fit(3,:)).*sqrt((path(1,:)-path_fit(1,:)).^2+(path(3,:)-path_fit(3,:)).^2);

figure
set(gcf,'outerposition',get(0,'screensize'));
% set(gca,'LooseInset',[0,0,0,0]);
bar(path(1,:),distance,'r')
hold on

plot([path(1,1) path(1,end)],[tb_r tb_r],'g')
hold on
plot([path(1,1) path(1,end)],[-tb_r -tb_r],'g')
hold off

ylim([-(tb_r+1) (tb_r+1)])
legend('Distance','upper tolerance band','lower tolerance band')
xlabel('x','Interpreter','latex') 
ylabel('Distance','Interpreter','latex') 
title('Distance between path_{orig} and path_{fit} at each sampling point')


% 5.4 third, compare the curvature of path_opt2 and path_wlTB
figure
set(gcf,'outerposition',get(0,'screensize'));
% set(gca,'LooseInset',[0,0,0,0]);
plot(path_opt2(1,:),optCur2{1},'r')
hold on

plot(path_wlTB(1,:),optCur_wlTB{1},'g')
hold off

xlabel('x','Interpreter','latex') 
ylabel('Curvature','Interpreter','latex') 
legend('Curvature of path_{opt2}','Curvature of path_{wlTB}')
title('Curvature of optimized path with method 2 and with linearly defined tolerance band')


% 5.5 finally, show the curvature reduction rate(CRR) at different tolerance band radius r
r = 4:0.5:8;
% These data were obtained by gradually increasing the tolerance band radius r from 4 to 8 
% These data were obtained when "step_xy = 2"
% We want the result in percentage, so we first multiply 100 in front of each data
CRR_opt2 = 100*(3.4635 - [3.3152 3.2947 3.2750 3.2632 3.2478 3.2343 3.2192 3.2049 3.1925])/3.4635;
CRR_wlTB = 100*(3.4635 - [3.3642 3.3478 3.3303 3.3114 3.2906 3.2682 3.2419 3.2143 3.1934])/3.4635;
figure
set(gcf,'outerposition',get(0,'screensize'));

plot(r,CRR_opt2,'o-','LineWidth',2)
hold on
plot(r,CRR_wlTB,'r*--','LineWidth',2)
hold off

% for y-axis, we use '%'
ytickformat('percentage')

xlabel('Tolerance band radius r','Interpreter','latex')
ylabel('CRR','Interpreter','latex')
legend('CRR_{opt2}','CRR_{wlTB}')
title('Curvature reduction rate(CRR) at different tolerance band radius r')


