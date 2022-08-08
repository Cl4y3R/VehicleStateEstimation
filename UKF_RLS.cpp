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
    double cf = 140000;
    double cr = 89000;
    double cx = 127000;
    //calling ukf

    Eigen::MatrixXd X_p = Eigen::MatrixXd::Zero(3,1);
    Eigen::MatrixXd Z_mea = Eigen::MatrixXd::Zero(2, 2 * 3 + 1);
    Eigen::MatrixXd Q_state = Eigen::MatrixXd::Identity(3,3);
    Eigen::MatrixXd R_mea = Eigen::MatrixXd::Identity(2,2);
    Q_state(0,0) = 1e-1;
    Q_state(1,1) = 1e5;
    Q_state(2,2) = 1e5;
    R_mea(0,0) = 1e-3;
    R_mea(1,1) = 1e-3;
    UKF BetaEstimator(3,2,Q_state,R_mea,0.001);

    std::vector<double> beta_est_deg;
    //ukf main loop
    for(int i=0;i<incsv.size();i++){
        //init
        BetaEstimator.init();
        double beta = BetaEstimator.X_out(0,0);
        double fxf = BetaEstimator.X_out(1,0);
        double fxr = BetaEstimator.X_out(2,0);
        double alphaf = incsv[i].deltaf - beta - lf*incsv[i].r/incsv[i].vx;
        double alphar = - beta + lr*incsv[i].r/incsv[i].vx;
        X_p(0,0) = -incsv[i].r + (fxf*sin(incsv[i].deltaf - beta) + cf*alphaf*cos(incsv[i].deltaf - beta) + cr*alphar*cos(-beta) + fxr*sin(-beta))/(m*incsv[i].vx);
        X_p(1,0) = 0;
        X_p(2,0) = 0;
        BetaEstimator.Z_hat(0,0) = incsv[i].ax;
        BetaEstimator.Z_hat(1,0) = incsv[i].ay;

        //predict
        BetaEstimator.predict(X_p);
        //measurement of sigma points
        for(int i=0;i < Z_mea.cols();i++){
            double beta_sig = BetaEstimator.X_sig_pre(0,i);
            double fxf_sig = BetaEstimator.X_sig_pre(1,i);
            double fxr_sig = BetaEstimator.X_sig_pre(2,i);
            double alphaf_sig = incsv[i].deltaf - beta_sig - lf*incsv[i].r/incsv[i].vx;
            double alphar_sig = - beta_sig + lr*incsv[i].r/incsv[i].vx;
            Z_mea(0,i) = (-cf*alphaf_sig*sin(incsv[i].deltaf) + fxf_sig*cos(incsv[i].deltaf)) / m; //ax
            Z_mea(1,i) = (cf*alphaf_sig*cos(incsv[i].deltaf) + cr*alphar_sig + fxf_sig*sin(incsv[i].deltaf) + fxr_sig) / m; //ay
            //std::cout<<"alpfaf_sig_"<<i<<" = "<<alphaf_sig<<", alpfar_sig_"<<i<<" = "<<alphar_sig<<", ax_"<<i<<" = "<<Z_mea(0,i)<<", ay_"<<i<<" ="<<Z_mea(1,i)<<std::endl;
        }
        //std::cout<<"ax = "<< BetaEstimator.Z_hat(0,0)<<", ay = "<< BetaEstimator.Z_hat(1,0)<<std::endl;
        BetaEstimator.measurement(Z_mea);
        BetaEstimator.update();
        std::cout<<"Estimated beta: "<<BetaEstimator.X_out(0,0)*180/PI<<", real beta: "<<incsv[i].beta<<std::endl;
        /*std::cout<<"Estimated Fxf: "<<BetaEstimator.X_out(1,0)<<std::endl;
        std::cout<<"Estimated Fxr: "<<BetaEstimator.X_out(2,0)<<std::endl;
        std::cout<<"P_pre= "<<BetaEstimator.P_pre(0,0)<<", "<<BetaEstimator.P_pre(0,1)<<", "<<BetaEstimator.P_pre(0,2)<<std::endl;
        std::cout<<"       "<<BetaEstimator.P_pre(1,0)<<", "<<BetaEstimator.P_pre(1,1)<<", "<<BetaEstimator.P_pre(1,2)<<std::endl;
        std::cout<<"       "<<BetaEstimator.P_pre(2,0)<<", "<<BetaEstimator.P_pre(2,1)<<", "<<BetaEstimator.P_pre(2,2)<<std::endl;
        std::cout<<"P_mea= "<<BetaEstimator.P_mea(0,0)<<", "<<BetaEstimator.P_mea(0,1)<<std::endl;
        std::cout<<"       "<<BetaEstimator.P_mea(1,0)<<", "<<BetaEstimator.P_mea(1,1)<<std::endl;
        std::cout<<"X_pre_mean = "<<BetaEstimator.X_pre_mean(0,0)<<", "<<BetaEstimator.X_pre_mean(1,0)<<", "<<BetaEstimator.X_pre_mean(2,0)<<std::endl;
        std::cout<<"Z_pre_mean = "<<BetaEstimator.Z_pre_mean(0,0)<<", "<<BetaEstimator.Z_pre_mean(1,0)<<std::endl;
        std::cout<<"Z_predicted_in = "<<Z_mea(0,0)<<", and "<<Z_mea(1,0)<<std::endl;
        std::cout<<"P= "<<BetaEstimator.P_pre(0,0)<<", "<<BetaEstimator.P_pre(0,1)<<", "<<BetaEstimator.P_pre(0,2)<<std::endl;
        std::cout<<"       "<<BetaEstimator.P_pre(1,0)<<", "<<BetaEstimator.P_pre(1,1)<<", "<<BetaEstimator.P_pre(1,2)<<std::endl;
        std::cout<<"       "<<BetaEstimator.P_pre(2,0)<<", "<<BetaEstimator.P_pre(2,1)<<", "<<BetaEstimator.P_pre(2,2)<<std::endl;*/
        
        //std::cout<<"-------------------------------------------"<<std::endl;
        beta_est_deg.push_back(BetaEstimator.X_out(0,0)*180/PI);
    }

//write data to csv
std::ofstream outFile;
outFile.open("D:\\00test\\data_est_rlsukf.csv", std::ios::out);
for(int i=0;i<beta_est_deg.size();i++){
    outFile << beta_est_deg[i] << std::endl;
}
outFile.close();
system("pause");
return 0;
}