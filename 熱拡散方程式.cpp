#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <omp.h>
#include "Eigen/Core"

// 各種定数を設定
const double dt = 0.2;  // 時間離散化の単位
const double d = 0.01;  // 空間離散化の単位

const double t_start = 0.;   // 開始時間
const double t_end = 5000.; // 終了時間
const double N = int((t_end - t_start) / dt); // 計算回数

// ヒートシンクのサイズ設定（長方形）
const double l = 0.7;          // 横幅
const double w = 0.5;          // 厚さ
const double i_max = int(w/d); // 厚さマス数
const double j_max = int(l/d); // 横幅マス数

// アルミ
const double rho_1    = 2688.; // 密度
const double c_1     = 905.;   // 比熱
const double lambda_1 = 237.;  // 熱伝導率
const double alpha_1  = lambda_1 / rho_1 / c_1; // 熱拡散率

// 空気 T =300K
const double rho_2    = 1.1763;       // 密度
const double c_2     = 1007.;         // 比熱
const double lambda_2 = 26.14 / 1000; // 熱伝導率
const double alpha_2  = lambda_2 / rho_2 / c_2 / 50; // 熱拡散率

class heatsink_data
{
    Eigen::MatrixXd T_new, T_old, 
};