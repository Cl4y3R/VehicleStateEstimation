#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "UKF.hpp"

#define PI 3.1415926

class csvdata{
public:
    double beta;
    double vx;
    double ax;
    double ay;
    double r;
    double deltaf;
};

int main()
{
    //read data from csv
    std::vector<csvdata> incsv;
    csvdata data;
    FILE *fp;
    fp=fopen("/Users/xujunqing/Downloads/VehicleStateEstimation-main/data.csv","r");
    while(1){
        fscanf(fp,"%lf,%lf,%lf,%lf,%lf,%lf",&data.beta,&data.vx,&data.ax,&data.ay,&data.r,&data.deltaf);
        incsv.push_back(data);
        if (feof(fp))break;
    }
    fclose(fp);
    std::cout<<"*************Showing first ten lines**************"<<std::endl;
    for(int i=0;i<10;i++)
    {
        std::cout<<incsv[i].beta<<" "<<incsv[i].vx<<" "<<incsv[i].ax<<" "<<incsv[i].ay<<std::endl;

    }
    std::cout<<"*************Read complete**************"<<std::endl;
    std::cout<<"                                        "<<std::endl;
    //vehicle parameters
    double lf = 1.110;
    double lr = 1.756;
    double m = 1549.8;
    double iz = 2315.3;
    double cf = 100000;
    double cr = 100000;
    double cx = 127000;
    //calling ukf

    Eigen::MatrixXd X_pA = Eigen::MatrixXd::Zero(5,1);
    Eigen::MatrixXd X_pB = Eigen::MatrixXd::Zero(3,1);
    Eigen::MatrixXd Z_meaA = Eigen::MatrixXd::Zero(4, 2 * 5 + 1);
    Eigen::MatrixXd Z_meaB = Eigen::MatrixXd::Zero(3, 2 * 3 + 1);
    Eigen::MatrixXd Q_stateA = Eigen::MatrixXd::Identity(5,5);
    Eigen::MatrixXd R_meaA = Eigen::MatrixXd::Identity(4,4);
    Eigen::MatrixXd Q_stateB = Eigen::MatrixXd::Identity(3,3);
    Eigen::MatrixXd R_meaB = Eigen::MatrixXd::Identity(3,3);
    Q_stateA(0,0) = 1e-5;
    Q_stateA(1,1) = 1e-6;
    Q_stateA(2,2) = 1e3;
    Q_stateA(3,3) = 1e3;
    Q_stateA(4,4) = 1e3;
    R_meaA(0,0) = 1e-4;
    R_meaA(1,1) = 1e-2;
    R_meaA(2,2) = 1e-3;
    R_meaA(3,3) = 1e-3;

    Q_stateB(0,0) = 1e-13;
    Q_stateB(1,1) = 1e-10;
    Q_stateB(2,2) = 1e-10;
    R_meaB(0,0) = 1e3;
    R_meaB(1,1) = 1e3;
    R_meaB(2,2) = 1e-1;
//  Q_stateA(0,0) = 1;
//  Q_stateA(1,1) = 1;
//  Q_stateA(2,2) = 1;
//  Q_stateA(3,3) = 1;
//  Q_stateA(4,4) = 1;
//  R_meaA(0,0) = 1;
//  R_meaA(1,1) = 1;
//  R_meaA(2,2) = 1;
//  R_meaA(3,3) = 1;
//  
//  Q_stateB(0,0) = 1;
//  Q_stateB(1,1) = 1;
//  Q_stateB(2,2) = 1;
//  R_meaB(0,0) = 1;
//  R_meaB(1,1) = 1;
//  R_meaB(2,2) = 1;

    UKF ForceEstimator(5,4,Q_stateA,R_meaA,0.01);
    UKF BetaEstimator(3,3,Q_stateB,R_meaB,0.01);

    std::vector<double> beta_est_deg;
    double beta_last = 0;
    //ukf main loop
    for(int i=0;i<incsv.size();i++){
        //Force Estimation
        //init
        ForceEstimator.init();
        //create sigma points
        ForceEstimator.sigma();
        for(int k=0;k<ForceEstimator.Xd_sig.cols();k++){
            double vx = ForceEstimator.X_sig(0,k);
            double r = ForceEstimator.X_sig(1,k);
            double fyf = ForceEstimator.X_sig(2,k);
            double fyr = ForceEstimator.X_sig(3,k);
            double fx = ForceEstimator.X_sig(4,k);
            ForceEstimator.Xd_sig(0,k) = (fx*cos(beta_last - incsv[i].deltaf) + fyf*sin(beta_last - incsv[i].deltaf) + fyr*sin(beta_last))/m;
            ForceEstimator.Xd_sig(1,k) = (lf* (fyf*cos(incsv[i].deltaf) + fx*sin(incsv[i].deltaf)) - lr*fyr)/iz;
        }
        //predict
        ForceEstimator.predict();
        //measurement
        ForceEstimator.Z_hat(0,0) = incsv[i].vx;
        ForceEstimator.Z_hat(1,0) = incsv[i].r;
        ForceEstimator.Z_hat(2,0) = incsv[i].ax;
        ForceEstimator.Z_hat(3,0) = incsv[i].ay;
        for(int k=0;k<ForceEstimator.Xd_sig.cols();k++){
            double vx_pre = ForceEstimator.X_sig_pre(0,k);
            double r_pre = ForceEstimator.X_sig_pre(1,k);
            double fyf_pre = ForceEstimator.X_sig_pre(2,k);
            double fyr_pre = ForceEstimator.X_sig_pre(3,k);
            double fx_pre = ForceEstimator.X_sig_pre(4,k);
            ForceEstimator.Z_sig_pre(0,k) = vx_pre; //vx
            ForceEstimator.Z_sig_pre(1,k) = r_pre; //r
            ForceEstimator.Z_sig_pre(2,k) = (-fyf_pre*sin(incsv[i].deltaf) + fx_pre*cos(incsv[i].deltaf))/m; //ax
            ForceEstimator.Z_sig_pre(3,k) = (fyf_pre*cos(incsv[i].deltaf) + fyr_pre + fx_pre*sin(incsv[i].deltaf))/m; //ay
        }
        ForceEstimator.measurement();
        //update
        ForceEstimator.update();

        //Beta Estimation
        //init
        BetaEstimator.init();
        //create sigma points
        BetaEstimator.sigma();
        for(int k=0;k<BetaEstimator.Xd_sig.cols();k++){
            double beta = BetaEstimator.X_sig(0,k);
            double cfd = BetaEstimator.X_sig(1,k);
            double crd = BetaEstimator.X_sig(2,k);
            double alphaf = incsv[i].deltaf - beta - lf*incsv[i].r/incsv[i].vx;
            double alphar = - beta + lr*incsv[i].r/incsv[i].vx;
            double fx_est = ForceEstimator.X_sig(4,k);
            BetaEstimator.Xd_sig(0,k) = -incsv[i].r + (fx_est*sin(incsv[i].deltaf - beta) + (cf+cfd)*alphaf*cos(incsv[i].deltaf - beta) + (cr+crd)*alphar*cos(-beta))/(m*incsv[i].vx);
        }
        //predict
        BetaEstimator.predict();
        //measurement
        BetaEstimator.Z_hat(0,0) = ForceEstimator.X_out(2,0);
        BetaEstimator.Z_hat(1,0) = ForceEstimator.X_out(3,0);
        BetaEstimator.Z_hat(2,0) = incsv[i].ay;
        for(int k=0;k<BetaEstimator.Xd_sig.cols();k++){
            double beta_sig = BetaEstimator.X_sig_pre(0,k);
            double cfd_sig = BetaEstimator.X_sig_pre(1,k);
            double crd_sig = BetaEstimator.X_sig_pre(2,k);
            double alphaf_sig = incsv[i].deltaf - beta_sig - lf*incsv[i].r/incsv[i].vx;
            double alphar_sig = - beta_sig + lr*incsv[i].r/incsv[i].vx;
            double fx_est = ForceEstimator.X_sig(4,k);
            BetaEstimator.Z_sig_pre(0,k) = (cf+cfd_sig)*alphaf_sig; //Fyf
            BetaEstimator.Z_sig_pre(1,k) = (cr+crd_sig)*alphar_sig; //Fyr
            BetaEstimator.Z_sig_pre(2,k) = (BetaEstimator.Z_sig_pre(0,k)*cos(incsv[i].deltaf) + fx_est*sin(incsv[i].deltaf))/m;
        }
        BetaEstimator.measurement();
        //update
        BetaEstimator.update();
        
        //print the results
        std::cout<<"Estimated beta(deg): "<<BetaEstimator.X_out(0,0)*180/PI<<", real beta(deg): "<<incsv[i].beta<<std::endl;
        std::cout<<"Fyf: "<<ForceEstimator.X_out(2,0)<<std::endl;
        std::cout<<"K_00: "<<ForceEstimator.P_pre(0,0)<<std::endl;
        beta_last = BetaEstimator.X_out(0,0);
        //std::cout<<"-------------------------------------------"<<std::endl;
        beta_est_deg.push_back(BetaEstimator.X_out(0,0)*180/PI);
    }

//write data to csv
std::ofstream outFile;
outFile.open("data_est_dukf.csv", std::ios::out);
for(int i=0;i<beta_est_deg.size();i++){
    outFile << beta_est_deg[i] << std::endl;
}
outFile.close();
system("pause");
return 0;
}
