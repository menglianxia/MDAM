%% MDAM (an innovative multi-temporal digital elevation model (DEM) adjustment model)
%% By Menglian Xia
%% Calculate correction parameters for sub-DEM in X-, Y and Z directions (dX, dY, a₀, a₁)
%{
Input:
    1. obs matrix (formed by GCPs and TPs), 3m*n
    m is the total number of GCPs and TPs, n corresonding 14 attributes of observations;
    each row in represent a pair of GCP of TP points;
    columns 1-14 represents 14 fields: X1, Y1, Z1, t1, VX1, VY1, DEM_ID1, X2, Y2, Z2, t2, VX2, VY2, DEM_ID2；
    therefore, observations form a matrix of m*14.
    for GCPs (ICESat-2 points), DEM_id1 = -1 or DEM_id2 = -1;
    2. DEM_center matrix, N*3
    formed by three paramesters: X, Y (central coordinates of DEM), sita (the Angle between the Y-axis of 
    the local coordinate system and the Y-axis of the polar projection coordinate system)
    N is the number of sub DEMs
    3. P (the weight matrix, formed according the manuscirpt), m*1

    Since the calculation of X, Y and Z directions are independent, the observations are organized by X, Y and Z, respectively. 
    First, complete the observation equations in the X direction [obs(1:m, 1:n)], 
    then move on to the Y direction [obs(m+1:2m, n+1:2n)], 
    and finally to the Z direction [obs(2m+1:3m, 2n+1:3n)]. 

Output:
    1. X: the correction parameters of sub DEMs (unknowns in the observation equation), 4N*1
    [dX_1; dX_2; ...; dX_N; dY_1; dY_2; ...; dY_N; a₀_1; a₁_1; a₀_2; a₁_2; ...; a₀_N; a₁_N]
    2. eps_val (Ɛ): the residuals of observations, 3m*1
    [Ɛ_X1; Ɛ_X2; ...; Ɛ_Xm; Ɛ_Y1; Ɛ_Y2; ...; Ɛ_Ym; Ɛ_Z1; Ɛ_Z2; ...; Ɛ_Zm]
    3. pts_corrections (DEM_ID, dX, dY, a₀, a₁, dZ): the corrections and the errors for all observations on the DEM
    for all points on an individual DEM, dX is the same value, as well as dY, a₀ and a₁; dZ = a₀+a₁*Y_OL, which is related to location (Y_OL)
    4. sigma_pts_corrections: the errors corresponding to pts_corrections
    5. obs_adjusted: the corrected values of the observed values
%}

function [X,eps_val, pts_corrections, sigma_pts_corrections, obs_adjusted] = MDAM(obs, DEM_center, P)
% obtain the number of DEMs (N) and points
% found unique DEM_id
DEM_id_list = [obs(:,7);obs(:,14)];
DEM_id_list = unique(DEM_id_list);
DEM_id_list = DEM_id_list(DEM_id_list>0);
% obtain the number of DEMs and the number of points (GCPs + TPs)
DEM_count = size(DEM_id_list,1);
point_count = size(obs,1);

% 1. formulate the observation equation
% (a) the constant vector of the observation equation
L = zeros(point_count*3,1); 
% (b) the known parameters of the observation equation
B = zeros(point_count*3,1);
% (c) coefficient matrix of the observation equation
A = zeros(point_count*3,DEM_count*4);
% the coefficients of a₀ and a₁ (used to calculate the errors of a₀, a₁ and dz)
f_T = [];

for i = 1:point_count
    % the equation of the point i listed in the X, Y and Z directions
    L(i,1) = obs(i,8) - obs(i,1);  % in X direction
    L(point_count+i,1) = obs(i,9) - obs(i,2);  % in Y direction
    L(point_count*2+i,1) = obs(i,10) - obs(i,3);  % in Z direction

    % calculate the average velocity of a pair of points
    vx = (obs(i,5)+obs(i,12))/2;
    vy = (obs(i,6)+obs(i,13))/2;

    B(i,1) = vx*(obs(i,11)-obs(i,4));
    B(point_count+i,1) = vy*(obs(i,11)-obs(i,4));
    B(point_count*2+i,1) = 0;

    DEM_id1 = obs(i,7);
    DEM_id2 = obs(i,14);
    
    if DEM_id1 == -1 & DEM_id2 > 0
        % point (X1, Y1, Z1) is an ICESat-2 point, point (X2, Y2, Z2) is a point on the DEM (GCP)
        A(i,DEM_id2) = -1;
        A(point_count+i,DEM_count+DEM_id2) = -1;
        % Since the coefficient of a0 is 1 and differs too much in order of magnitude from the coefficient Y1 of a₁, 
        % it is easy to cause the matrix to be singular. Therefore, in the calculation, the coefficient of a₀ is set to 10000
        A(point_count*2+i,DEM_count*2+DEM_id2*2-1) = -10000;
        % calculate Y_OL
        sita = DEM_center(DEM_id2, 3);
        Y_OL = (obs(i,9)-DEM_center(DEM_id2,2))*cos(sita/180*pi)-(obs(i,8)-DEM_center(DEM_id2,1))*sin(sita/180*pi);
        A(point_count*2+i,DEM_count*2+DEM_id2*2) = -Y_OL;
        f_T = [f_T;DEM_id2,i,10000,Y_OL];
    end
    
    if DEM_id1 > 0 & DEM_id2 == -1
         % point (X1, Y1, Z1) is a point on the DEM, point (X2, Y2, Z2) is an ICESat-2 point(GCP)
        A(i,DEM_id1) = 1;
        A(point_count+i,DEM_count+DEM_id1) = 1;
        A(point_count*2+i,DEM_count*2+DEM_id1*2-1) = 10000;
        sita = DEM_center(DEM_id1, 3);
        Y_OL = (obs(i,2)-DEM_center(DEM_id1,2))*cos(sita/180*pi)-(obs(i,1)-DEM_center(DEM_id1,1))*sin(sita/180*pi);
        A(point_count*2+i,DEM_count*2+DEM_id1*2) = Y_OL;
        f_T = [f_T;DEM_id1,i,10000,Y_OL];
    end
    
    if DEM_id1 > 0 && DEM_id2 > 0
        % both point (X1, Y1, Z1) is and point (X2, Y2, Z2) are on DEMs (TP)
        A(i,DEM_id1) = 1;
        A(point_count+i,DEM_count+DEM_id1) = 1;
        A(point_count*2+i,DEM_count*2+DEM_id1*2-1) = 10000;
        sita = DEM_center(DEM_id1, 3);
        Y_OL = (obs(i,2)-DEM_center(DEM_id1,2))*cos(sita/180*pi)-(obs(i,1)-DEM_center(DEM_id1,1))*sin(sita/180*pi);
        A(point_count*2+i,DEM_count*2+DEM_id1*2) = Y_OL;
        f_T = [f_T;DEM_id1,i,10000,Y_OL];

        A(i,DEM_id2) = -1;
        A(point_count+i,DEM_count+DEM_id2) = -1;
        A(point_count*2+i,DEM_count*2+DEM_id2*2-1) = -10000;
        sita = DEM_center(DEM_id2, 3);
        Y_OL = (obs(i,9)-DEM_center(DEM_id2,2))*cos(sita/180*pi)-(obs(i,8)-DEM_center(DEM_id2,1))*sin(sita/180*pi);
        A(point_count*2+i,DEM_count*2+DEM_id2*2) = -Y_OL;
        f_T = [f_T;DEM_id2,i,10000,Y_OL];
    end
end

% Construct the weight matrix corresponding to the observation matrix
P_XYZ = [P;P;P];
P = diag(P_XYZ,0);

% 2. calculate unknowns and residuals
l = L-B;
NA = A'*P*A;
U = A'*P*l; 
X = NA\U; 
eps_val = A*X+B-L; 

% 3. accuracy assessment of the adjustment results
% cofactor matrix
Q = inv(NA);

% the minimum number of required observations(n)
% the number of observations (t)
% sigma0 = sqrt(eps_val'*eps_val/(n-t));
% The unit-weight standard deviations for the three directions are evaluated independently.
eps_val_X = eps_val(1:point_count,1);
P_X = P(1:point_count,1:point_count);
eps_val_Y = eps_val(point_count+1:point_count*2,1);
P_Y = P(point_count+1:point_count*2,point_count+1:point_count*2);
eps_val_Z = eps_val(point_count*2+1:point_count*3,1); 
P_Z = P(point_count*2+1:point_count*3,point_count*2+1:point_count*3);
% Each DEM in the X and Y directions has one unknown parameter. Each DEM in the Z direction has two unknown parameters
sigma_dX0 = sqrt(eps_val_X'*P_X*eps_val_X/(point_count-DEM_count));
sigma_dY0 = sqrt(eps_val_Y'*P_Y*eps_val_Y/(point_count-DEM_count));
sigma_dZ0 = sqrt(eps_val_Z'*P_Z*eps_val_Z/(point_count-DEM_count*2));

% Extract the diagonal elements of the Q matrix to compute the unit-weight standard deviation.
Q_diag = diag(Q);
%  the standard deviations of the parameters dx, dy, a₀, a₁
sigma_para(1:DEM_count,1) = sigma_dX0*sqrt(Q_diag(1:DEM_count,1));
sigma_para(DEM_count+1:DEM_count*2,1) = sigma_dY0*sqrt(Q_diag(DEM_count+1:DEM_count*2,1));
sigma_para(DEM_count*2+1:DEM_count*4,1) = sigma_dZ0*sqrt(Q_diag(DEM_count*2+1:DEM_count*4,1));
% error propagation: dZ = a₀+a₁*Y_OL; a=[a₀ a₁]'; f_T = [1 Y_OL],dz=f_T*a; Q_dZdZ=f_T*Q_aa*f_T'; sigma_dZ = sigma_dZ0*Q_dzdz
Q_aa = Q(DEM_count*2+1:end,DEM_count*2+1:end);

% calculate the corrections and the errors for all observations on the DEM
% sorted by the DEM sequence number (primary) and point number (secondary)
f_T = sortrows(f_T,[1 2]);
pts_corrections = [];
sigma_pts_corrections = [];
for i = 1:size(f_T,1)
    j = f_T(i,1); % DEM_id
    f_T_i = f_T(i,3:4);
    % unknowns (dx,dy,a₀,a₁) and the calculated dz
    dx_j = X(j); 
    dy_j = X(DEM_count+j); 
    a0_j = X(DEM_count*2+j*2-1)*10000; 
    a1_j = X(DEM_count*2+j*2);
    dz_i = a0_j + a1_j*f_T_i(2);
    pts_corrections = [pts_corrections;j, dx_j,dy_j,a0_j,a1_j,dz_i];
    % the errora of dx,dy,a₀,a₁ and dz (error propagation)
    sigma_dx_j = sigma_para(j);  
    sigma_dy_j = sigma_para(DEM_count+j);
    sigma_a0_j = sigma_para(DEM_count*2+j*2-1);
    sigma_a1_j = sigma_para(DEM_count*2+j*2);
    Q_aa_j = Q_aa(j*2-1:j*2,j*2-1:j*2);
    Q_ptsI_dz = f_T_i*Q_aa_j*f_T_i';
    sigma_pts_dz_i = sigma_dZ0*sqrt(Q_ptsI_dz);
    sigma_pts_corrections = [sigma_pts_corrections;j, sigma_dx_j,sigma_dy_j,sigma_a0_j*10000,sigma_a1_j,sigma_pts_dz_i];
end

% 4. Calculate the corrected values of the observed values
obs_adjusted = obs;
for i = 1:point_count
    DEM_id1 = obs_adjusted(i,7);
    DEM_id2 = obs_adjusted(i,14);
    if DEM_id1 ~= -1
        dX = X(DEM_id1);
        dY = X(DEM_count+DEM_id1);
        a0 = X(DEM_count*2+DEM_id1*2-1);
        a1 = X(DEM_count*2+DEM_id1*2);
        sita = DEM_center(DEM_id1, 3);
        Y_OL = (obs(i,9)-DEM_center(DEM_id1,2))*cos(sita/180*pi)-(obs(i,8)-DEM_center(DEM_id1,1))*sin(sita/180*pi);
        obs_adjusted(i,1) = obs_adjusted(i,1) + dX;
        obs_adjusted(i,2) = obs_adjusted(i,2) + dY;
        obs_adjusted(i,3) = obs_adjusted(i,3) + a0*10000 + a1*Y_OL;
    end
    if DEM_id2 ~= -1
        dX = X(DEM_id2);
        dY = X(DEM_count+DEM_id2);
        a0 = X(DEM_count*2+DEM_id2*2-1);
        a1 = X(DEM_count*2+DEM_id2*2);
        sita = DEM_center(DEM_id2, 3);
        Y_OL = (obs(i,9)-DEM_center(DEM_id2,2))*cos(sita/180*pi)-(obs(i,8)-DEM_center(DEM_id2,1))*sin(sita/180*pi);
        obs_adjusted(i,8) = obs_adjusted(i,8) + dX;
        obs_adjusted(i,9) = obs_adjusted(i,9) + dY;
        obs_adjusted(i,10) = obs_adjusted(i,10) + a0*10000 + a1*Y_OL;
    end
end
save('result','X','eps_val', 'pts_corrections', 'sigma_pts_corrections', 'obs_adjusted');
end