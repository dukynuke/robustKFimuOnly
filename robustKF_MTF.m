% Robust Attitude Estimation via Multiple Tuning Factors (MTFs) using IMU-Only Measurements 
% Batu Candan and Halil Ersin Soken, 2021

% Process Model
F = eye(3) - dt*skew(wx,wy,wz) ;
Q = -skew(x_c(1),x_c(2),x_c(3))*g_var*skew(x_c(1),x_c(2),x_c(3))*dt^2 ;
% wx,wy,wz -> Gyroscope Measurements (rad/s).
% skew is the skew-symmetric function transforming the vector to the matrix.
% Q is the process noise covariance matrix.
% g_var is the gyroscope noise variance.

% Measurement Model
H = norm(g)*eye(3) ; 
R = a_var ;
z = [ax;ay;az] - c_a*a_est ;
% H is the measurement matrix.
% g is the local gravity vector (m/s^2). 
% R is the measurement noise covariance matrix.
% ax,ay,az -> Accelerometer Measurements (m/s^2).
% a_var is the accelerometer variance.
% a_est is the estimation for external acceleration (m/s^2).
% z is the measurement matrix.
% c_a is a dimensionless, determinator constant specified for the cutoff  
% frequency and the value of this constant is varying between 0 and 1.

% Time Update (Prediction)
x_p = F*x_c ;    % State Prediction
P = F*P*F' + Q ; % Covariance Prediction
e = z - H*x_p ;  % Innovation in KF

% R-Adaptation via Novel Multiple Tuning Factor (MTF) Method
ek = e ;
S = ek*ek'-H*P*H'-R ;         % Multiple 
diags = diag(S) ;             % Tuning 
S = diag([max(diags(1),0)...  % Factor 
          max(diags(2),0)...  % Construction and
          max(diags(3),0)]) ; % Diagonalization.

% Measurement Update (Correction)
K = P*H'*inv(H*P*H'+S+R) ;              % Kalman Gain Evaulation
x_c = x_p + K*e ; x_c = x_c/norm(x_c) ; % State Correction and Normalization
P = (eye(3)-K*H)*P ;                    % Covariance Correction

% External Acc. Estimation (Gravity Compensation)
a_est = [ax;ay;az] - norm(g)*x_c ;

r_hat = atan2(x_c(2),x_c(3)) ;                   % Roll Estimate (rad)
p_hat = atan2(-x_c(1),sqrt(x_c(2)^2+x_c(3)^2)) ; % Pitch Estimate (rad)

ax_hat = a_est(1) ;        % 
ay_hat = a_est(2) ;        % External Acc. Estimates (m/s^2)
az_hat = a_est(3) - g(3) ; %