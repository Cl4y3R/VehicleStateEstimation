# VehicleStateEstimation
Online vehicle state estimation algorithm for vehicle motion control using Unscented Kalman Filter(UKF).

## UKF core
The UKF core is created as pure header file(hpp), thus using it would be simple and easy: include it in your source code and create ukf estimator by calling `UKF yourEstimator`.

The UKF core contains 4 functions: `void init(), void predict(Eigen::MatrixXd X_d), void measurement(Eigen::MatrixXd Z_in), void update()`. To use UKF you must call these 4 functions exactly in the given order.

## UKF with offline identification of cornering stiffness
First, identify the cornering stiffness(linear region) by using Least Square method. The measuring parameters are:  
-Derivative of yawrate `$ a $`  
-Lateral acceleration `$ a_y $`  
`$$A=1$$`  
横摆角加速度$\ddot{\psi}$（微分可得），横向加速度$a_y$  

参数：侧偏角$\beta$，车速$v$，横摆角速度$\dot{\psi}$

测量量$y$：横摆角加速度$\ddot{\psi}$（微分可得），横向加速度$a_y$

标定量$\theta$：前后轴侧偏刚度$C_f,C_r$

则有$y=A\theta$

## DUKF with parallel estimation of cornering stiffness
