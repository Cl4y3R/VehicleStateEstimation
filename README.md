# VehicleStateEstimation
Online vehicle state estimation algorithm for vehicle motion control using Unscented Kalman Filter(UKF).

## UKF core
The UKF core is created as pure header file(hpp), thus using it would be simple and easy: include it in your source code and create ukf estimator by calling `UKF yourEstimator`.

The UKF core contains 4 functions: `void init(), void predict(Eigen::MatrixXd X_d), void measurement(Eigen::MatrixXd Z_in), void update()`. To use UKF you must call these 4 functions exactly in the given order.

## UKF with offline identification of cornering stiffness
First, identify the cornering stiffness(linear region) by using Least Square method.  
The measuring parameters are:  
-Derivative of yawrate $\ddot{\psi}$  
-Lateral acceleration $a_y$  
The used parameter are:  
-Vehicle slip angle $\beta$  
-Vehicle speed $v$  
-Yawrate $\dot{\psi}$  
The identified parameter are:  
-Conrering stiffness $C_f,C_r$  

$A=\begin{bmatrix}\frac{\delta_f-\beta-\frac{l_f\dot{\psi}}{v}}{m}&\frac{-\beta+\frac{l_r\dot{\psi}}{v}}{m}\\\frac{l_f(\delta_f-\beta-\frac{l_f\dot{\psi}}{v})}{J_z}&\frac{l_f(-\beta+\frac{l_r\dot{\psi}}{v})}{J_z}\end{bmatrix}$

> $E_k = Y_k - A_k\theta_k$
> 
> 
> $K_{k} = P_kA^T(\lambda I+A_kP_kA^T)^{-1}$
> 
> $P_{k+1} = \frac{1}{\lambda}(P_k - K_{k}A_kP_k)$
> 
> $\theta_{k+1} =\theta_k+K_{k}E_k$
> 

## DUKF with parallel estimation of cornering stiffness
