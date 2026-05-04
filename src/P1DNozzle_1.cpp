//亚声速-超声速等熵喷管(非守恒型)
#include<iostream>
#include<cstdio>
#include<cmath>
#include<iomanip>
#include<string>
#include <fstream>

using namespace std;

constexpr double Gamma = 1.4;
constexpr double R = 287.0;
constexpr double cfl = 0.5;
constexpr double deltaX = 0.1;
constexpr int cellNum = 30;
constexpr int timesteps = 1400;

void exportToFile(const string& filename,double X[],double A[],double rho[],
                  double V[],double T[],double p[],double Ma[],
                  double mdot[],int cellNum) {
    ofstream fout(filename);
    if (!fout) {
        cerr << "无法创建文件: " << filename << endl;
        return;
    }

    fout << fixed << setprecision(3);

    fout << "No" << ',' << "X" << ',' << "A" << ',' << "rho" << ','
         << "V" << ',' << "T" << ',' << "p" << ',' 
         << "Ma" << ',' <<"mdot" <<endl;
    for (int i = 0; i < cellNum + 1; i++) {
        fout << i + 1 << ','
             << X[i] << ',' << A[i] << ',' << rho[i] << ','
             << V[i] << ',' << T[i] << ',' << p[i] << ',' 
             << Ma[i] <<',' <<mdot[i] << endl;
    }
}

int main(){
    double X[cellNum+1],A[cellNum+1],rho[cellNum+1],V[cellNum+1],T[cellNum+1],
    p[cellNum+1],Ma[cellNum+1],mdot[cellNum+1];
    double drhodt[cellNum+1],dVdt[cellNum+1],dTdt[cellNum+1];
    double rhop[cellNum+1],Vp[cellNum+1],Tp[cellNum+1];
    double drhopdt[cellNum+1],dVpdt[cellNum+1],dTpdt[cellNum+1];
    double drhodtav[cellNum+1],dVdtav[cellNum+1],dTdtav[cellNum+1];
    double deltaT = 1.00,pdeltaT = 1.00;

    //给定初值条件
    for(int i=0;i<cellNum+1;i++){
        X[i] = i * deltaX;
        A[i] = 1.0 + 2.2 * (X[i]-1.5)*(X[i]-1.5);
        rho[i] = 1.0 - 0.314 * X[i];
        T[i] = 1.0 - 0.2314 * X[i];
        V[i] = (0.1 + 1.09 * X[i])*sqrt(T[i]);
        Ma[i] = V[i]/sqrt(T[i]);
        mdot[i] = rho[i] * V[i] * A[i];
    }

    exportToFile("./output/Q1_initval.csv",X,A,rho,V,T,p,Ma,mdot,cellNum);

    for (int j=0;j<timesteps;j++){
    //Maccormack
    //确定时间步长
    for(int i = 0;i < cellNum+1; i++){
        pdeltaT = (cfl * deltaX)/(V[i]+sqrt(T[i]));
        if(deltaT>pdeltaT)
        deltaT=pdeltaT;
    }

    //预估步
    for (int i = 0; i < cellNum; i++){
	drhodt[i] = -V[i] * (rho[i + 1] - rho[i]) / deltaX
		- rho[i] * (V[i + 1] - V[i]) / deltaX
		- rho[i] * V[i] * (log(A[i + 1]) - log(A[i])) / deltaX;//7-51,p215
	dVdt[i] = -V[i] * (V[i + 1] - V[i]) / deltaX
		- (1 / Gamma) * ((T[i + 1] - T[i]) / deltaX + T[i] / rho[i] * (rho[i + 1] - rho[i]) / deltaX);//7-52
	dTdt[i] = -V[i] * (T[i + 1] - T[i]) / deltaX
		- (Gamma - 1) * T[i] * ((V[i + 1] - V[i]) / deltaX + V[i] * (log(A[i + 1]) - log(A[i])) / deltaX);//7-53
	rhop[i] = rho[i] + drhodt[i] * deltaT;//7-54,p215
	Vp[i] = V[i] + dVdt[i] * deltaT;//7-55
	Tp[i] = T[i] + dTdt[i] * deltaT;//7-56
    }

    //校正步
    for (int i = 1; i < cellNum;i++){
	drhopdt[i] = -Vp[i] * (rhop[i] - rhop[i-1]) / deltaX
		- rhop[i] * (Vp[i] - Vp[i - 1]) / deltaX
		- rhop[i] * Vp[i] * (log(A[i]) - log(A[i - 1])) / deltaX;//7-57,p216
	dVpdt[i] = -Vp[i] * (Vp[i] - Vp[i - 1]) / deltaX
		- (1 / Gamma) * ((Tp[i] - Tp[i - 1]) / deltaX + Tp[i] / rhop[i] * (rhop[i] - rhop[i - 1]) / deltaX);//7-58
	dTpdt[i] = -Vp[i] * (Tp[i] - Tp[i - 1]) / deltaX
		- (Gamma - 1) * Tp[i] * ((Vp[i] - Vp[i - 1]) / deltaX + Vp[i] * (log(A[i]) - log(A[i - 1])) / deltaX);//7-59
	drhodtav[i] = (drhodt[i] +drhopdt[i])/2;//7-60,p215
	dVdtav[i] = (dVdt[i] +dVpdt[i])/2;//7-61
	dTdtav[i] = (dTdt[i] +dTpdt[i])/2;//7-62
	rho[i] += drhodtav[i] * deltaT;//7-63
	V[i] += dVdtav[i] * deltaT;//7-64
	T[i] += dTdtav[i] * deltaT;//7-65
    }
    for(int i=1;i < cellNum; i++){
        p[i] = rho[i] * T[i];
        Ma[i] = V[i]/sqrt(T[i]);
        mdot[i] = rho[i] * V[i] * A[i];
    }

    //头尾部
    V[0] = 2 * V[1] - V[2];
    p[0] = rho[0] * T[0];
    Ma[0] = V[0]/sqrt(T[0]);
    mdot[0] = rho[0] * V[0] * A[0];
    V[cellNum] = 2 * V[cellNum-1] - V[cellNum-2];
    rho[cellNum] = 2 * rho[cellNum-1] - rho[cellNum-2];
    T[cellNum] = 2 * T[cellNum-1] - T[cellNum-2];
    p[cellNum] = rho[cellNum] * T[cellNum];
    Ma[cellNum] = V[cellNum]/sqrt(T[cellNum]);
    mdot[cellNum] = rho[cellNum] * V[cellNum] * A[cellNum];

    if(j==0){
       exportToFile("./output/Q1_iteration1.csv",X,A,rho,V,T,p,Ma,mdot,cellNum);
    }
    }

    //输出timesteps步后的数据
    exportToFile("./output/Q1_final.csv",X,A,rho,V,T,p,Ma,mdot,cellNum);

    return 0;
}
