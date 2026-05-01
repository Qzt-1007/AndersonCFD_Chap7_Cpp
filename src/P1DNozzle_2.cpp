//亚声速-超声速等熵喷管(守恒型)
#include<iostream>
#include<cstdio>
#include<cmath>
#include<iomanip>
#include<string>
#include <fstream>

using namespace std;

constexpr double gamma = 1.4;
constexpr double R = 287.0;
constexpr double cfl = 0.5;
constexpr double deltaX = 0.1;
constexpr int cellNum = 30;
constexpr int timesteps = 1400;

void exportToFile(const string& filename,double X[],double A[],double rho[],
                  double V[],double T[],double p[],double Ma[],
                  double mdot[],double U1[],double U2[],double U3[],int cellNum) {
    ofstream fout(filename);
    if (!fout) {
        cerr << "无法创建文件: " << filename << endl;
        return;
    }

    fout << fixed << setprecision(3);

    fout << "No" << '\t' << "X" << '\t' << "A" << '\t' << "rho" << '\t'
         << "V" << '\t' << "T" << '\t' << "p" << '\t'
         << "Ma" << '\t' <<"mdot" << '\t' << "U1" << '\t' << "U2" << '\t'
         << "U3" << endl;
    for (int i = 0; i < cellNum + 1; i++) {
        fout << i + 1 << '\t'
             << X[i] << '\t' << A[i] << '\t' << rho[i] << '\t'
             << V[i] << '\t' << T[i] << '\t' << p[i] << '\t' 
             << Ma[i] <<'\t' <<mdot[i] << '\t' << U1[i] << '\t' << U2[i] <<'\t'
             << U3[i] <<endl;
    }
}

double F1(double U2){
    return U2;
}

double F2(double U1,double U2,double U3){
    double res;
    res = ((U2 * U2)/U1) + ((gamma - 1.0)/gamma)*(U3-((0.5 * gamma * U2 * U2)/U1));
    return res;
}

double F3(double U1,double U2,double U3){
    double res;
    res = (gamma*U2*U3/U1)-((gamma*(gamma-1.0)*U2*U2*U2)/(2*U1*U1));
    return res;
}

double Tcal(double U1,double U2,double U3){
    return (gamma-1.0) * ((U3/U1)-(( 0.5 * gamma * U2 *U2)/(U1*U1)));
}

//rho = U1/A
//J2与A有关，此处不给出

int main(){
    double U1[cellNum+1],U2[cellNum+1],U3[cellNum+1],pU1[cellNum+1],pU2[cellNum+1],pU3[cellNum+1];
    double X[cellNum+1],A[cellNum+1],rho[cellNum+1],V[cellNum+1],T[cellNum+1],
    p[cellNum+1],Ma[cellNum+1],mdot[cellNum+1],prho[cellNum+1],pT[cellNum+1];
    double dU1dt[cellNum+1] ,dU2dt[cellNum+1] ,dU3dt[cellNum+1] ,
    pdU1dt[cellNum+1] ,pdU2dt[cellNum+1] ,pdU3dt[cellNum+1],
    dU1dtav[cellNum+1] ,dU2dtav[cellNum+1] ,dU3dtav[cellNum+1];
    double deltaT = 1e6 ,pdeltaT = 1.00;

    //给定初值条件
    for(int i=0;i<cellNum+1;i++){
        X[i] = i * deltaX;
        A[i] = 1.0 + 2.2 * (X[i]-1.5)*(X[i]-1.5);
        if(i<=5){
            rho[i] = 1.0;
            T[i] = 1.0;
        }
        else if(i>=6 && i<=15){
            rho[i] = 1.0 - 0.366 * (X[i] - 0.5);
            T[i] = 1.0 - 0.167 * (X[i] - 0.5);
        }
        else{
            rho[i] = 0.634 - 0.3879 * (X[i]-1.5);
            T[i] = 0.833 - 0.3507 * (X[i]-1.5);
        }
        V[i] = 0.59/(rho[i]*A[i]);
        Ma[i] = V[i]/sqrt(T[i]);
        mdot[i] = 0.59;
        p[i] = rho[i]*T[i];
        U1[i] = rho[i]*A[i];
        U2[i] = 0.59;
        U3[i] = rho[i]*A[i]*((T[i]/(gamma-1))+(gamma*V[i]*V[i]/2));
    }

    exportToFile("../output/Q2_initval.csv",X,A,rho,V,T,p,Ma,mdot,U1,U2,U3,cellNum);

    for (int j=0;j<timesteps;j++){
    //Maccormack
    //确定时间步长
    deltaT = 1e6;
    for(int i = 0;i < cellNum+1; i++){
        pdeltaT = (cfl * deltaX)/(V[i]+sqrt(T[i]));
        if(deltaT>pdeltaT)
        deltaT=pdeltaT;
    }

    //预估步
    for (int i = 0; i < cellNum; i++){
	dU1dt[i] = -(F1(U2[i+1])-F1(U2[i]))/deltaX;
    dU2dt[i] = -(F2(U1[i+1],U2[i+1],U3[i+1])-F2(U1[i],U2[i],U3[i]))/deltaX
    + rho[i] * T[i] * (A[i+1]-A[i])/(gamma * deltaX);
    dU3dt[i] = -(F3(U1[i+1],U2[i+1],U3[i+1])-F3(U1[i],U2[i],U3[i]))/deltaX;
    pU1[i] = U1[i] + dU1dt[i] * deltaT;
    pU2[i] = U2[i] + dU2dt[i] * deltaT;
    pU3[i] = U3[i] + dU3dt[i] * deltaT;
    prho[i] = pU1[i]/A[i];
    pT[i] = Tcal(pU1[i],pU2[i],pU3[i]);
    }

    //校正步
    for (int i = 1; i < cellNum;i++){
	pdU1dt[i] = -(F1(pU2[i])-F1(pU2[i-1]))/deltaX;
    pdU2dt[i] = -(F2(pU1[i],pU2[i],pU3[i])-F2(pU1[i-1],pU2[i-1],pU3[i-1]))/deltaX
    + prho[i] * pT[i] * (A[i]-A[i-1])/(gamma * deltaX);
    pdU3dt[i] = -(F3(pU1[i],pU2[i],pU3[i])-F3(pU1[i-1],pU2[i-1],pU3[i-1]))/deltaX;
    dU1dtav[i] = (dU1dt[i]+pdU1dt[i])/2;
    dU2dtav[i] = (dU2dt[i]+pdU2dt[i])/2;
    dU3dtav[i] = (dU3dt[i]+pdU3dt[i])/2;
    U1[i] = U1[i] + dU1dtav[i] * deltaT;
    U2[i] = U2[i] + dU2dtav[i] * deltaT;
    U3[i] = U3[i] + dU3dtav[i] * deltaT;
    }

    //其他变量计算
    for(int i=1;i < cellNum; i++){
        rho[i] = U1[i]/A[i];
        V[i] = U2[i]/U1[i];
        T[i] = Tcal(U1[i],U2[i],U3[i]);
        mdot[i] = U2[i];
        Ma[i] = V[i] / sqrt(T[i]);
        p[i] = rho[i]*T[i];
    }

    //头尾部
    U2[0] = 2 * U2[1] - U2[2];
    V[0] = U2[0]/U1[0];
    U3[0] = U1[0] * ((1/(gamma-1.0))+(gamma * V[0]*V[0]/2));
    mdot[0] = U2[0];
    p[0] = rho[0]*T[0];
    Ma[0] = V[0] / sqrt(T[0]);
    U1[cellNum] = 2 * U1[cellNum-1] - U1[cellNum-2];
    U2[cellNum] = 2 * U2[cellNum-1] - U2[cellNum-2];
    U3[cellNum] = 2 * U3[cellNum-1] - U3[cellNum-2];
    rho[cellNum] = U1[cellNum]/A[cellNum];
    V[cellNum] = U2[cellNum]/U1[cellNum];
    T[cellNum] = Tcal(U1[cellNum],U2[cellNum],U3[cellNum]);
    mdot[cellNum] = U2[cellNum];
    p[cellNum] = rho[cellNum]*T[cellNum];
    Ma[cellNum] = V[cellNum] / sqrt(T[cellNum]);

    if(j==0){
       exportToFile("../output/Q2_iteration1.csv",X,A,rho,V,T,p,Ma,mdot,U1,U2,U3,cellNum);
    }
    }

    //输出timesteps步后的数据
    exportToFile("../output/Q2_final.csv",X,A,rho,V,T,p,Ma,mdot,U1,U2,U3,cellNum);

    return 0;
}
