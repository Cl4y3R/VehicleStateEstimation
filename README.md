# VehicleStateEstimation
Online vehicle state estimation algorithm for vehicle motion control using Unscented Kalman Filter(UKF).  
**NON-COMMERCIAL USE!**

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

$\lambda \in[0.98,1]$  

Then use the UKF for vehicle slip angle online estimation.
The estimated parameters are:  
-Vehicle slip angle $\beta$  
-Vehicle longitudinal forces $F_{xf}, F_{xr}$  
The inputs are:  
-Front wheel steering angle $\delta_f$  
-Yawrate $\dot{\psi}$  
The measurements are:  
-Longitudinal acceleration $a_x$  
-Lateral acceleration $a_y$  
The state space:  
$\dot{\beta}=-\dot{\psi}+\frac{1}{mv}[F_{xf}sin(\delta_f-\beta)+F_{yf}cos(\delta_f-\beta)+F_{xr}sin(-\beta)+F_{yr}cos(-\beta)]$
$\dot{F_{xf}}=0$  
$\dot{F_{xr}}=0$  
The measurement equation:  
$y_1=a_x=\frac{1}{m}[-F_{yf}sin(\delta_f)+F_{xf}cos(\delta_f)+F_{xr}]$
$y_2=a_y=\frac{1}{m}[F_{yf}cos(\delta_f)+F_{yr}+F_{xf}sin(\delta_f))]$  

## DUKF with parallel estimation of cornering stiffness  
When offline identification of the cornering stiffness is not avaliable, you can use the parallel structure of the UKF(DUKF) to estimate the vehicle slip angle as well as the cornering stiffness by setting a approximated initial value of the cornering stiffness. The UKF will change the cornering stiffness value simutanoeusly.  

**Observer A**

The estimated parameters are:  
-Vehicle speed $v$  
-Yawrate $\dot{\psi}$  
-Front axle lateral force $F_{yf}$  
-Rear axle lateral force $F_{yr}$  
-Front axle longitudinal force $F_{xf}$  
The inputs are:  
-Front wheel steering angle $\delta_f$  
-Vehicle slip angle $\beta$ 
The measurements are:  
-Vehicle speed $v$ 
-Yawrate $\dot{\psi}$  
-Longitudinal acceleration $a_x$  
-Lateral acceleration $a_y$  
The state space:  
$\dot{v}=\frac{1}{m}[F_{xf}cos(\beta-\delta_f)+F_{yf}sin(\beta-\delta_f)+F_{yr}sin(\beta)]$  
$\ddot{\psi}=\frac{1}{Jz}[(F_{yf}cos(\delta_f)+F_{xf}sin(\delta_f))lf + F_{yr}l_r]$  
$\dot{F_x}=\dot{F}\_{yf}=\dot{F}\_{yr}=0$  
The measurements are:  
$y_1=v=v$  
$y_2=\dot{\psi}=\dot{\psi}$  
$y_3=a_x=\frac{1}{m}[-F_{yf}sin(\delta_f)+F_{xf}cos(\delta_f)]$  
$y_4=a_y=\frac{1}{m}[F_{yf}cos(\delta_f+F_{yr}+F_{xf}sin(\delta_f))]$  

**Observer B**

The estimated parameters are:  
-Vehicle slip angle $\beta$ 
-Front axle cornering stiffness changing part $\Delta c_f$  
-Rear axle cornering stiffness changing part $\Delta c_r$  
The inputs are: 
-Front wheel steering angle $\delta_f$  
-Yawrate $\dot{\psi}$  
-Front axle longitudinal force $F_{xf}$  
The measurements are:
-Front axle lateral force $F_{yf}$  
-Rear axle lateral force $F_{yr}$  
-Lateral acceleration $a_y$  
The state space:  
$\dot{\beta}=-\dot{\psi}+\frac{1}{mv}[F_{xf}sin(\delta_f-\beta)+F_{yf}cos(\delta_f-\beta)+F_{yr}cos(\beta)]$  
$F_{yf}=(c_f+\Delta c_f)(\delta_f-\beta-\frac{l_f\dot{\psi}}{v})$  
$F_{yr}=(c_r+\Delta c_r)(-\beta+\frac{l_r\dot{\psi}}{v})$  
$\Delta \dot{c}\_f = 0$  
$\Delta \dot{c}\_r = 0$  
The measurement equations:  
$y_1=F_{yf}=(c_f+\Delta c_f)(\delta_f-\beta-\frac{l_f\dot{\psi}}{v})$  
$y_2=F_{yr}=(c_r+\Delta c_r)(-\beta+\frac{l_r\dot{\psi}}{v})$  
$y_3=\frac{1}{m}[(c_f+\Delta c_f)(\delta_f-\beta-\frac{l_f\dot{\psi}}{v})cos(\delta_f)+(c_r+\Delta c_r)(-\beta+\frac{l_r\dot{\psi}}{v})+F_{xf}sin(\delta_f)]$  

