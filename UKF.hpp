#pragma once
#ifndef UKF_HPP
#define UKF_HPP

#include <iostream>
#include <Eigen/Dense>

class UKF{
    public:
        UKF(int dim_state, int dim_mea, Eigen::MatrixXd Q_in, Eigen::MatrixXd R_in,double sampletime):is_inited(false),just_begin_filt(false),n_state(dim_state),n_mea(dim_mea),Q(Q_in),R(R_in),dt(sampletime){
            std::cout<<"*********Unscented Kalman Filter created!*********"<<std::endl;
        };
        ~UKF(){
            std::cout<<"*********Unscented Kalman Filter destroyed!*********"<<std::endl;
        };

        void init();
        void predict(Eigen::MatrixXd X_d);
        void measurement(Eigen::MatrixXd Z_in);
        void update();


        double alpha = 1e-3;
        double kappa = 0;
        double b = 2;
        double dt;
        double lambda;
        int n_state;
        int n_mea;
        bool is_inited;
        bool just_begin_filt;

        Eigen::MatrixXd P_pre;
        Eigen::MatrixXd Q;
        Eigen::MatrixXd R;
        Eigen::VectorXd Wm = Eigen::VectorXd::Zero(2 * n_state + 1);
        Eigen::VectorXd Wc = Eigen::VectorXd::Zero(2 * n_state + 1);
        Eigen::MatrixXd X_sig = Eigen::MatrixXd::Zero(n_state, 2 * n_state + 1);
        Eigen::MatrixXd X_sig_pre = Eigen::MatrixXd::Zero(n_state, 2 * n_state + 1);
        Eigen::VectorXd X_pre_mean = Eigen::VectorXd::Zero(n_state);
        Eigen::MatrixXd Z_sig_pre = Eigen::MatrixXd::Zero(n_mea, 2 * n_state + 1);
        Eigen::VectorXd Z_pre_mean = Eigen::VectorXd::Zero(n_mea);
        Eigen::MatrixXd P_mea = Eigen::MatrixXd::Zero(n_mea,n_mea);
        Eigen::MatrixXd P_cross = Eigen::MatrixXd::Zero(n_state,n_mea);
        Eigen::MatrixXd K;
        Eigen::VectorXd X_out;
        Eigen::VectorXd Z_hat;
};

void UKF::init(){
    if (!is_inited){
        X_out = Eigen::VectorXd::Zero(n_state);
        for(int i=0;i<n_state;i++){
            X_out(i,0) = 0.0001;
        }
        Z_hat = Eigen::VectorXd::Zero(n_mea);
        P_pre = Eigen::MatrixXd::Identity(n_state,n_state);
        lambda = alpha * alpha * (n_state + kappa) - n_state;
        //sigma points weights
        //common choice
        Wm(0) = lambda / (n_state + lambda);
        Wc(0) = Wm(0) + (1 - alpha * alpha + b);
        for(int i = 1; i < 2 * n_state + 1;i++){
		    Wm(i) = 1 / (2 * (n_state + lambda));
		    Wc(i) = 1 / (2 * (n_state + lambda));
        }

        //matlab choice
        /*Wm(0) = 1 - n_state / (alpha * alpha *(n_state + kappa));
        Wc(0) = (2 - alpha * alpha + b) - n_state / (alpha * alpha *(n_state + kappa));
        for(int i = 1; i < 2 * n_state + 1;i++){
	    	Wm(i) = 1 / (2 * alpha * alpha * (n_state + kappa));
	    	Wc(i) = 1 / (2 * alpha * alpha * (n_state + kappa));
        }*/

        //original choice
        /*Wm(0) = 0 / (n_state + 0);
        Wc(0) = Wm(0);
        for(int i = 1; i < 2 * n_state + 1;i++){
		    Wm(i) = 1 / (2 * (n_state + 0));
		    Wc(i) = 1 / (2 * (n_state + 0));
        }*/
        is_inited = true;
        just_begin_filt = true;

    }
    else just_begin_filt = false;

}

void UKF::predict(Eigen::MatrixXd X_d){
    if (just_begin_filt)
        return;
    //create sigma points
    Eigen::MatrixXd A = P_pre.ldlt().matrixL();//robust cholesky decomposition,AA^T=P
    X_sig.fill(0.0);
    X_sig.col(0) = X_out;
    for (int i = 0; i < n_state; ++i) {
        X_sig.col(i + 1) = X_out + sqrt(lambda + n_state) * A.col(i);
        X_sig.col(i + 1 + n_state) = X_out - sqrt(lambda + n_state) * A.col(i);
    }

    //nonlinear transformation of the sigma points
    for(int i=0;i < 2 * n_state + 1;i++){		
        for(int j=0;j<n_state;j++){
            X_sig_pre(j,i) = X_sig(j,i) + X_d(j,0) * dt;
        }	
    }

    //prediction mean and covariance
    P_pre.fill(0.0);
    X_pre_mean.fill(0.0);
    for(int i = 0;i < 2 * n_state + 1;i++){
		X_pre_mean += Wm(i) * X_sig_pre.col(i);
		P_pre += Wc(i) * (X_sig_pre.col(i) - X_pre_mean) * (X_sig_pre.col(i) - X_pre_mean).transpose();
    }
    P_pre += Q;
}

void UKF::measurement(Eigen::MatrixXd Z_in){
    if (just_begin_filt)
        return;
    /*!!!!!!!!!!!WARNING!!!!!!!!!!!!!!
        Measurement of sigma points must be calculated before calling UKF::measurement
    */
    //read measurement of sigma points
    Z_sig_pre = Z_in;
    //measurement mean and covariance
    P_mea.fill(0.0);
    Z_pre_mean.fill(0.0);
    for(int i = 0;i < 2 * n_state + 1;i++){
        Z_pre_mean += Wm(i) * Z_sig_pre.col(i);
        P_mea += Wc(i) * (Z_sig_pre.col(i) - Z_pre_mean)* (Z_sig_pre.col(i) - Z_pre_mean).transpose();
    }
    P_mea += R;
    
    //cross covariance
    P_cross.fill(0.0);
    for(int i = 0;i < 2 * n_state + 1;i++){
        P_cross += Wc(i) * (X_sig_pre.col(i) - X_pre_mean) * (Z_sig_pre.col(i) - Z_pre_mean).transpose();
    }

    //kalman gain
    K = P_cross * P_mea.inverse();
}

void UKF::update(){
    if (just_begin_filt)
        return;
    X_out = X_pre_mean + K * (Z_hat - Z_pre_mean);
    P_pre += - K * P_mea * K.transpose();
}

#endif