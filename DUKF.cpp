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
    fp=fopen("D:\\00test\\data.csv","r");
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
    Q_stateA(2,2) = 1e-1;
    Q_stateA(3,3) = 1e-1;
    Q_stateA(4,4) = 1e-1;
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

    UKF ForceEstimator(5,4,Q_stateA,R_meaA,0.01);
    UKF BetaEstimator(3,3,Q_stateB,R_meaB,0.01);

    std::vector<double> beta_est_deg;
    double beta_last = 0;
    //ukf main loop
    for(int i=0;i<incsv.size();i++){
        //Force Estimation
        //init
        ForceEstimator.init();
        double vx = ForceEstimator.X_out(0,0);
        double r = ForceEstimator.X_out(1,0);
        double fyf = ForceEstimator.X_out(2,0);
        double fyr = ForceEstimator.X_out(3,0);
        double fx = ForceEstimator.X_out(4,0);
        X_pA(0,0) = (fx*cos(beta_last - incsv[i].deltaf) + fyf*sin(beta_last - incsv[i].deltaf) + fyr*sin(beta_last))/m;
        X_pA(1,0) = (lf* (fyf*cos(incsv[i].deltaf) + fx*sin(incsv[i].deltaf)) - lr*fyr)/m;
        ForceEstimator.Z_hat(0,0) = incsv[i].vx;
        ForceEstimator.Z_hat(1,0) = incsv[i].r;
        ForceEstimator.Z_hat(2,0) = incsv[i].ax;
        ForceEstimator.Z_hat(3,0) = incsv[i].ay;

        //predict
        ForceEstimator.predict(X_pA);
        //measurement of sigma points
        for(int i=0;i < Z_meaA.cols();i++){
            double vx_sig = ForceEstimator.X_sig_pre(0,i);
            double r_sig = ForceEstimator.X_sig_pre(1,i);
            double fyf_sig = ForceEstimator.X_sig_pre(2,i);
            double fyr_sig = ForceEstimator.X_sig_pre(3,i);
            double fx_sig = ForceEstimator.X_sig_pre(4,i);
            Z_meaA(0,i) = vx_sig; //vx
            Z_meaA(1,i) = r_sig; //r
            Z_meaA(2,i) = (-fyf_sig*sin(incsv[i].deltaf) + fx_sig*cos(incsv[i].deltaf))/m; //ax
            Z_meaA(3,i) = (fyf_sig*cos(incsv[i].deltaf) + fyr_sig + fx_sig*sin(incsv[i].deltaf))/m; //ay
        }
        ForceEstimator.measurement(Z_meaA);
        ForceEstimator.update();

        //Beta Estimation
        //init
        BetaEstimator.init();
        double beta = BetaEstimator.X_out(0,0);
        double cfd = BetaEstimator.X_out(1,0);
        double crd = BetaEstimator.X_out(2,0);
        double alphaf = incsv[i].deltaf - beta - lf*incsv[i].r/incsv[i].vx;
        double alphar = - beta + lr*incsv[i].r/incsv[i].vx;
        double fx_est = ForceEstimator.X_out(4,0);
        X_pB(0,0) = -incsv[i].r + (fx_est*sin(incsv[i].deltaf - beta) + (cf+cfd)*alphaf*cos(incsv[i].deltaf - beta) + (cr+crd)*alphar*cos(-beta))/(m*incsv[i].vx);
        BetaEstimator.Z_hat(0,0) = ForceEstimator.X_out(2,0);
        BetaEstimator.Z_hat(1,0) = ForceEstimator.X_out(3,0);
        BetaEstimator.Z_hat(2,0) = incsv[i].ay;

        //predict
        BetaEstimator.predict(X_pB);
        //measurement of sigma points
        for(int i=0;i < Z_meaB.cols();i++){
            double beta_sig = BetaEstimator.X_sig_pre(0,i);
            double cfd_sig = BetaEstimator.X_sig_pre(1,i);
            double crd_sig = BetaEstimator.X_sig_pre(2,i);
            double alphaf_sig = incsv[i].deltaf - beta_sig - lf*incsv[i].r/incsv[i].vx;
            double alphar_sig = - beta_sig + lr*incsv[i].r/incsv[i].vx;
            Z_meaB(0,i) = (cf+cfd_sig)*alphaf_sig; //Fyf
            Z_meaB(1,i) = (cr+crd_sig)*alphar_sig; //Fyr
            Z_meaB(2,i) = (Z_meaB(0,i)*cos(incsv[i].deltaf) + fx_est*sin(incsv[i].deltaf))/m;
        }
        BetaEstimator.measurement(Z_meaB);
        BetaEstimator.update();
        std::cout<<"Estimated beta: "<<BetaEstimator.X_out(0,0)*180/PI<<", real beta: "<<incsv[i].beta<<std::endl;
        beta_last = BetaEstimator.X_out(0,0);
        //std::cout<<"-------------------------------------------"<<std::endl;
        beta_est_deg.push_back(BetaEstimator.X_out(0,0)*180/PI);
    }

//write data to csv
std::ofstream outFile;
outFile.open("D:\\00test\\data_est_dukf.csv", std::ios::out);
for(int i=0;i<beta_est_deg.size();i++){
    outFile << beta_est_deg[i] << std::endl;
}
outFile.close();
system("pause");
return 0;
}