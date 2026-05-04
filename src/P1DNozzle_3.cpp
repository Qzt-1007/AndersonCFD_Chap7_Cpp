// 亚声速-超声速等熵喷管(守恒型,激波捕捉)
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

constexpr double Gamma = 1.400;
constexpr double cfl = 0.500;
constexpr double deltaX = 0.100;
constexpr int cellNum = 30;
constexpr int timesteps = 1600;
constexpr double Cx = 0.200;

void exportToFile(const string& filename,
                  const vector<double>& X,
                  const vector<double>& A,
                  const vector<double>& rho,
                  const vector<double>& V,
                  const vector<double>& T,
                  const vector<double>& p,
                  const vector<double>& Ma,
                  const vector<double>& mdot,
                  const vector<double>& U1,
                  const vector<double>& U2,
                  const vector<double>& U3) {
    ofstream fout(filename);
    if (!fout) {
        cerr << "无法创建文件: " << filename << endl;
        return;
    }

    fout << fixed << setprecision(3);
    fout << "No,X,A,rho,V,T,p,Ma,mdot,U1,U2,U3\n";
    for (size_t i = 0; i < X.size(); ++i) {
        fout << i + 1 << ',' << X[i] << ',' << A[i] << ',' << rho[i] << ','
             << V[i] << ',' << T[i] << ',' << p[i] << ',' << Ma[i] << ','
             << mdot[i] << ',' << U1[i] << ',' << U2[i] << ',' << U3[i] << '\n';
    }
}

inline double F1(double U2) { return U2; }

inline double F2(double U1, double U2, double U3) {
    return (U2 * U2) / U1 + (Gamma - 1.0) / Gamma *
           (U3 - (0.5 * Gamma * U2 * U2) / U1);
}

inline double F3(double U1, double U2, double U3) {
    return (Gamma * U2 * U3 / U1) -
           (Gamma * (Gamma - 1.0) * U2 * U2 * U2) / (2.0 * U1 * U1);
}

inline double Tcal(double U1, double U2, double U3) {
    return (Gamma - 1.0) * ((U3 / U1) - (0.5 * Gamma * U2 * U2) / (U1 * U1));
}

int main() {
    const int totalPoints = cellNum + 1;
    vector<double> X(totalPoints), A(totalPoints), rho(totalPoints), V(totalPoints),
                   T(totalPoints), p(totalPoints), Ma(totalPoints), mdot(totalPoints),
                   U1(totalPoints), U2(totalPoints), U3(totalPoints),
                   pU1(totalPoints), pU2(totalPoints), pU3(totalPoints),
                   prho(totalPoints), pT(totalPoints), pp(totalPoints),
                   dU1dt(totalPoints), dU2dt(totalPoints), dU3dt(totalPoints),
                   pdU1dt(totalPoints), pdU2dt(totalPoints), pdU3dt(totalPoints),
                   dU1dtav(totalPoints), dU2dtav(totalPoints), dU3dtav(totalPoints);

    // 初始化几何和流场
    for (int i = 0; i < totalPoints; ++i) {
        X[i] = i * deltaX;
        A[i] = 1.0 + 2.2 * (X[i] - 1.5) * (X[i] - 1.5);

        if (i <= 5) {
            rho[i] = 1.0;
            T[i] = 1.0;
        } else if (i <= 15) {
            rho[i] = 1.0 - 0.366 * (X[i] - 0.5);
            T[i] = 1.0 - 0.167 * (X[i] - 0.5);
        } else if (i <= 20) {
            rho[i] = 0.634 - 0.702 * (X[i] - 1.5);
            T[i] = 0.833 - 0.4908 * (X[i] - 1.5);
        } else {
            rho[i] = 0.5892 - 0.10228 * (X[i] - 2.1);
            T[i] = 0.93968 - 0.0622 * (X[i] - 2.1);
        }

        V[i] = 0.59 / (rho[i] * A[i]);
        Ma[i] = V[i] / sqrt(T[i]);
        mdot[i] = 0.59;
        p[i] = rho[i] * T[i];
        U1[i] = rho[i] * A[i];
        U2[i] = 0.59;
        U3[i] = rho[i] * A[i] * (T[i] / (Gamma - 1.0) + Gamma * V[i] * V[i] / 2.0);
    }

    exportToFile("./output/Q3_initval.csv", X, A, rho, V, T, p, Ma, mdot, U1, U2, U3);

    // 人工粘性计算 lambda
    auto artificialViscosity = [](double Cx, double p2, double p1, double p0,
                                  double u2, double u1, double u0) {
        double numerator = Cx * abs(p2 - 2.0 * p1 + p0) *
                           (u2 - 2.0 * u1 + u0);
        double denominator = p2 + 2.0 * p1 + p0;
        return numerator / denominator;
    };

    for (int iter = 0; iter < timesteps; ++iter) {
        // --- 计算时间步长 ---
        double deltaT = numeric_limits<double>::max();
        for (int i = 0; i < totalPoints; ++i) {
            double localDT = cfl * deltaX / (V[i] + sqrt(T[i]));
            if (localDT < deltaT) deltaT = localDT;
        }

        // --- 预估步 (i = 0 .. cellNum-1) ---
        for (int i = 0; i < cellNum; ++i) {
            dU1dt[i] = -(F1(U2[i + 1]) - F1(U2[i])) / deltaX;
            dU2dt[i] = -(F2(U1[i + 1], U2[i + 1], U3[i + 1]) -
                         F2(U1[i], U2[i], U3[i])) / deltaX +
                       rho[i] * T[i] * (A[i + 1] - A[i]) / (Gamma * deltaX);
            dU3dt[i] = -(F3(U1[i + 1], U2[i + 1], U3[i + 1]) -
                         F3(U1[i], U2[i], U3[i])) / deltaX;

            if (i == 0) {
                pU1[i] = U1[i] + dU1dt[i] * deltaT;
                pU2[i] = U2[i] + dU2dt[i] * deltaT;
                pU3[i] = U3[i] + dU3dt[i] * deltaT;
            } else {
                pU1[i] = U1[i] + dU1dt[i] * deltaT +
                         artificialViscosity(Cx, p[i + 1], p[i], p[i - 1],
                                             U1[i + 1], U1[i], U1[i - 1]);
                pU2[i] = U2[i] + dU2dt[i] * deltaT +
                         artificialViscosity(Cx, p[i + 1], p[i], p[i - 1],
                                             U2[i + 1], U2[i], U2[i - 1]);
                pU3[i] = U3[i] + dU3dt[i] * deltaT +
                         artificialViscosity(Cx, p[i + 1], p[i], p[i - 1],
                                             U3[i + 1], U3[i], U3[i - 1]);
            }

            prho[i] = pU1[i] / A[i];
            pT[i] = Tcal(pU1[i], pU2[i], pU3[i]);
            pp[i] = prho[i] * pT[i];
        }

        //预估步未计算最后一个节点，用外推补充
        pU1[cellNum] = 2.0 * pU1[cellNum - 1] - pU1[cellNum - 2];
        pU2[cellNum] = 2.0 * pU2[cellNum - 1] - pU2[cellNum - 2];
        pU3[cellNum] = 2.0 * pU3[cellNum - 1] - pU3[cellNum - 2];
        prho[cellNum] = pU1[cellNum] / A[cellNum];
        pT[cellNum] = Tcal(pU1[cellNum], pU2[cellNum], pU3[cellNum]);
        pp[cellNum] = prho[cellNum] * pT[cellNum];

        // --- 校正步 (i = 1 .. cellNum-1) ---
        for (int i = 1; i < cellNum; ++i) {
            pdU1dt[i] = -(F1(pU2[i]) - F1(pU2[i - 1])) / deltaX;
            pdU2dt[i] = -(F2(pU1[i], pU2[i], pU3[i]) -
                          F2(pU1[i - 1], pU2[i - 1], pU3[i - 1])) / deltaX +
                        prho[i] * pT[i] * (A[i] - A[i - 1]) / (Gamma * deltaX);
            pdU3dt[i] = -(F3(pU1[i], pU2[i], pU3[i]) -
                          F3(pU1[i - 1], pU2[i - 1], pU3[i - 1])) / deltaX;

            dU1dtav[i] = 0.5 * (dU1dt[i] + pdU1dt[i]);
            dU2dtav[i] = 0.5 * (dU2dt[i] + pdU2dt[i]);
            dU3dtav[i] = 0.5 * (dU3dt[i] + pdU3dt[i]);

            U1[i] = U1[i] + dU1dtav[i] * deltaT +
                    artificialViscosity(Cx, pp[i + 1], pp[i], pp[i - 1],
                                        pU1[i + 1], pU1[i], pU1[i - 1]);
            U2[i] = U2[i] + dU2dtav[i] * deltaT +
                    artificialViscosity(Cx, pp[i + 1], pp[i], pp[i - 1],
                                        pU2[i + 1], pU2[i], pU2[i - 1]);
            U3[i] = U3[i] + dU3dtav[i] * deltaT +
                    artificialViscosity(Cx, pp[i + 1], pp[i], pp[i - 1],
                                        pU3[i + 1], pU3[i], pU3[i - 1]);
        }

        // --- 更新原始变量 (i = 1 .. cellNum-1) ---
        for (int i = 1; i < cellNum; ++i) {
            rho[i] = U1[i] / A[i];
            V[i] = U2[i] / U1[i];
            T[i] = Tcal(U1[i], U2[i], U3[i]);
            mdot[i] = U2[i];
            Ma[i] = V[i] / sqrt(T[i]);
            p[i] = rho[i] * T[i];
        }

        // --- 边界条件处理 ---
        // 入口边界（保持 U1[0] 不变，外推 U2[0]）
        U2[0] = 2.0 * U2[1] - U2[2];
        V[0] = U2[0] / U1[0];
        U3[0] = U1[0] * (1.0 / (Gamma - 1.0) + Gamma * V[0] * V[0] / 2.0);
        mdot[0] = U2[0];
        p[0] = rho[0] * T[0];
        Ma[0] = V[0] / sqrt(T[0]);

        // 出口边界（外推所有守恒变量）
        U1[cellNum] = 2.0 * U1[cellNum - 1] - U1[cellNum - 2];
        U2[cellNum] = 2.0 * U2[cellNum - 1] - U2[cellNum - 2];
        U3[cellNum] = 2.0 * U3[cellNum - 1] - U3[cellNum - 2];
        rho[cellNum] = U1[cellNum] / A[cellNum];
        V[cellNum] = U2[cellNum] / U1[cellNum];
        T[cellNum] = Tcal(U1[cellNum], U2[cellNum], U3[cellNum]);
        mdot[cellNum] = U2[cellNum];
        p[cellNum] = rho[cellNum] * T[cellNum];
        Ma[cellNum] = V[cellNum] / sqrt(T[cellNum]);

        // 保存第一次迭代结果
        if (iter == 0) {
            exportToFile("./output/Q3_iteration1.csv", X, A, rho, V, T, p, Ma, mdot, U1, U2, U3);
        }
    }

    // 最终结果输出
    exportToFile("./output/Q3_final.csv", X, A, rho, V, T, p, Ma, mdot, U1, U2, U3);
    return 0;
}