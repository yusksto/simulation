#include <iostream>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <omp.h>
#include "Eigen/Core"

constexpr auto PI = 3.141592653589793238462643383279;

bool TF(double, double, double, double, double);
void simurate(double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d&, Eigen::Vector3d&);
Eigen::Vector3d RungeKutta_r(Eigen::Vector3d, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d);
Eigen::Vector3d RungeKutta_v(Eigen::Vector3d, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d);
Eigen::Vector3d cal_a(Eigen::Vector3d, Eigen::Vector3d&, Eigen::Vector3d&);


//シミュレーション各定数
const double t_0 = 0;    //開始時刻
const double dt = 0.001; //刻み幅
const double N = 100;     //公転回数

const double G = 100.0;   //万有引力定数
const double M = 3000.0; //恒星質量
const double m = 100.0;   //惑星質量
const double R_m = 100.0; //惑星公転半径

const double R_M = R_m * m / M;             //恒星回転半径
const double R = R_m + R_M;                 //恒星-惑星間距離
const double ω = sqrt(G * M / R_m) / R;    //公転角周波数
const double t_1 = t_0 + 2.0 * PI * N / ω; //終了時刻

//gif出力
const double range = R_m * 1.75;         //座標表示範囲
const int n = 90;                        //一周当たりのgif枚数
const double p = 2.0 * PI / ω / dt / n; //出力頻度


int main()
{
    std::cout << p << std::endl;
    //gnuplot設定
    FILE* gp = _popen("gnuplot", "w");
    assert(gp);
    fprintf(gp, "reset\n");
    fprintf(gp, "set terminal gif animate optimize delay %d size %d,%d\n", (int)(5), (int)(1920), (int)(1080));
    fprintf(gp, "set output 'tmp.gif'\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "unset autoscale x\n");
    fprintf(gp, "unset autoscale y\n");
    fprintf(gp, "unset autoscale z\n");
    fprintf(gp, "set xr [%lf:%lf]\n", -range, range);
    fprintf(gp, "set yr [%lf:%lf]\n", -range, range);
    fprintf(gp, "set zr [%lf:%lf]\n", -range, range);
    fprintf(gp, "set size square\n");
    fprintf(gp, "set view 45, 60, 1, 1\n");


    //質点初期化
    const double dr = 1;
    const double dv = 0.1;
    const double pv = 1;
    int k = 0;
    for (double v_ = -pv; v_ < pv; v_ += dv)
    {
        for (double x = -range; x < range; x += dr)
        {
            for (double y = -range; y < range; y += dr)
            {
                if (TF(x, y, 0, 0.9, 0.1))
                {
                    k++;
                }
            }
        }
    }
    std::cout << k << std::endl;
    int a = 0;
    Eigen::MatrixXd r = Eigen::MatrixXd::Zero(3, k);
    Eigen::MatrixXd v = Eigen::MatrixXd::Zero(3, k);
    Eigen::VectorXi tf = Eigen::VectorXi::Ones(k);
    for (double v_ = -pv; v_ < pv; v_ += dv)
    {
        for (double x = -range; x < range; x += dr)
        {
            for (double y = -range; y < range; y += dr)
            {
                if (TF(x, y, 0, 0.9, 0.1))
                {
                    r(0, a) = x;
                    r(1, a) = y;
                    r(2, a) = 0;

                    v(0, a) = -y * ω * sqrt(x * x + y * y) / sqrt(x * x + y * y) * (1 - v_);
                    v(1, a) = x * ω * sqrt(x * x + y * y) / sqrt(x * x + y * y) * (1 - v_);
                    v(2, a) = 0;
                    a++;
                }
            }
        }
    }


    //シミュレーション
    Eigen::Vector3d r_(0, 0, 0);
    Eigen::Vector3d v_(0, 0, 0);
    int b = (int)p; //出力頻度調整用
    for (double t = t_0; t < t_1; t += dt)
    {
        b++;

        //gif出力
        if (b > p)
        {
            std::cout << 100 * (t - t_0) / (t_1 - t_0) << std::endl;


            fprintf(gp, "set multiplot layout 1,2\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Y'\n");
            fprintf(gp, "set title 'X-Y t=%lf'\n", t);
            fprintf(gp, "plot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf\n", -R_M * cos(ω * t), -R_M * sin(ω * t));
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf\n", R_m * cos(ω * t), R_m * sin(ω * t));
            fprintf(gp, "e\n");
#pragma omp parallel for
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf\n", r(0, i), r(1, i));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Z'\n");
            fprintf(gp, "set title 'X-Z t=%lf'\n", t);
            fprintf(gp, "plot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf\n", -R_M, 0.0);
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf\n", R_m, 0.0);
            fprintf(gp, "e\n");
#pragma omp parallel for
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf\n", r(0, i) * cos(ω * t) + r(1, i) * sin(ω * t), r(1, i) * cos(ω * t) - r(0, i) * sin(ω * t));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "unset multiplot\n");
            fflush(gp);


            b = 0;
        }

#pragma omp parallel for
        for (int i = 0; i < k; i++)
        {
            if (tf(i))
            {
                Eigen::Vector3d r_;
                Eigen::Vector3d v_;

                r_(0) = r(0, i);
                r_(1) = r(1, i);
                r_(2) = r(2, i);

                v_(0) = v(0, i);
                v_(1) = v(1, i);
                v_(2) = v(2, i);

                simurate(t, r_, v_, r_, v_);

                r(0, i) = r_(0);
                r(1, i) = r_(1);
                r(2, i) = r_(2);

                v(0, i) = v_(0);
                v(1, i) = v_(1);
                v(2, i) = v_(2);

                if (TF(r_(0), r_(1), r_(2), 0.9, 0.1) == false)
                {
                    tf(i) = 0;
                }
            }
        }
    }

    //終了
    fprintf(gp, "set out\n");
    _pclose(gp);
    return 0;
}

bool TF(double x, double y, double z, double dr_1, double dr_2) {
    return R_m * (1 - dr_1) < sqrt(x * x + y * y + z * z) && sqrt(x * x + y * y + z * z) < R_m * (1 + dr_2);
}

void simurate(double t, Eigen::Vector3d r_0, Eigen::Vector3d v_0, Eigen::Vector3d& r_1, Eigen::Vector3d& v_1) {

    Eigen::Vector3d r_M(-R_M * cos(ω * t), -R_M * sin(ω * t), 0);
    Eigen::Vector3d r_m(R_m * cos(ω * t), R_m * sin(ω * t), 0);

    r_1 = RungeKutta_r(r_0, r_M, r_m, v_0);
    v_1 = RungeKutta_v(r_0, r_M, r_m, v_0);
}

Eigen::Vector3d RungeKutta_r(Eigen::Vector3d r, Eigen::Vector3d& r_M, Eigen::Vector3d& r_m, Eigen::Vector3d v) {
    Eigen::Vector3d k_1 = RungeKutta_v(r, r_M, r_m, v);
    Eigen::Vector3d k_2 = RungeKutta_v(r, r_M, r_m, v + dt / 2 * k_1);
    Eigen::Vector3d k_3 = RungeKutta_v(r, r_M, r_m, v + dt / 2 * k_2);
    Eigen::Vector3d k_4 = RungeKutta_v(r, r_M, r_m, v + dt * k_3);
    return r + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
}

Eigen::Vector3d RungeKutta_v(Eigen::Vector3d r, Eigen::Vector3d& r_M, Eigen::Vector3d& r_m, Eigen::Vector3d v) {
    Eigen::Vector3d k_1 = cal_a(r, r_M, r_m);
    Eigen::Vector3d k_2 = cal_a(r + dt / 2 * k_1, r_M, r_m);
    Eigen::Vector3d k_3 = cal_a(r + dt / 2 * k_2, r_M, r_m);
    Eigen::Vector3d k_4 = cal_a(r + dt * k_3, r_M, r_m);
    return v + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
}

Eigen::Vector3d cal_a(Eigen::Vector3d r, Eigen::Vector3d& r_M, Eigen::Vector3d& r_m) {
    return G * M * (r_M - r) * pow(pow(r(0) - r_M(0), 2) + pow(r(1) - r_M(1), 2) + pow(r(2) - r_M(2), 2), -1.5) + G * m * (r_m - r) * pow(pow(r(0) - r_m(0), 2) + pow(r(1) - r_m(1), 2) + pow(r(2) - r_m(2), 2), -1.5);
}

/*
 //質点初期化
    const double dr = 1;
    int k = 0;
    for (double x = -range; x < range; x += dr)
    {
        for (double y = -range; y < range; y += dr)
        {
            for (double z = -range; z < range; z += dr)
            {
                if (TF(x, y, z, 0.05, 0.05))
                {
                    k++;
                }
            }
        }
    }
    std::cout << k << std::endl;
    int a = 0;
    Eigen::MatrixXd r  = Eigen::MatrixXd::Zero(3, k);
    Eigen::MatrixXd v  = Eigen::MatrixXd::Zero(3, k);
    Eigen::VectorXi tf = Eigen::VectorXi::Ones(k);
    for (double x = -range; x < range; x += dr)
    {
        for (double y = -range; y < range; y += dr)
        {
            for (double z = -range; z < range; z += dr)
            {
                if (TF(x, y, z, 0.05, 0.05))
                {
                    r(0, a) = x;
                    r(1, a) = y;
                    r(2, a) = z;

                    v(0, a) = -y * ω * sqrt(x * x + y * y + z * z) / sqrt(x * x + y * y);
                    v(1, a) = x * ω * sqrt(x * x + y * y + z * z) / sqrt(x * x + y * y);
                    v(2, a) = 0;
                    a++;
                }
            }
        }
    }
*/

/*
//質点初期化
    const double dr = 5;
    const double dv = 0.5;
    const double pv = 1;
    int k = 0;
    for (double v_ = -pv; v_ < pv; v_ += dv)
    {
        for (double x = -range; x < range; x += dr)
        {
            for (double y = -range; y < range; y += dr)
            {
                if (TF(x, y, 0, 0.8, 0.2))
                {
                    k++;
                }
            }
        }
    }
    std::cout << k << std::endl;
    int a = 0;
    Eigen::MatrixXd r  = Eigen::MatrixXd::Zero(3, k);
    Eigen::MatrixXd v  = Eigen::MatrixXd::Zero(3, k);
    Eigen::VectorXi tf = Eigen::VectorXi::Ones(k);
    for (double v_ = -pv; v_ < pv; v_ += dv)
    {
        for (double x = -range; x < range; x += dr)
        {
            for (double y = -range; y < range; y += dr)
            {
                if (TF(x, y, 0, 0.8, 0.2))
                {
                    r(0, a) = x;
                    r(1, a) = y;
                    r(2, a) = 0;

                    v(0, a) = -y * ω * sqrt(x * x + y * y) / sqrt(x * x + y * y) * (1 - v_);
                    v(1, a) = x * ω * sqrt(x * x + y * y) / sqrt(x * x + y * y) * (1 - v_);
                    v(2, a) = 0;
                    a++;
                }
            }
        }
    }
*/

/*

fprintf(gp, "set multiplot layout 1,2\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Y'\n");
            fprintf(gp, "set zlabel 'Z'\n");
            fprintf(gp, "set title 'X-Y-Z t=%lf'\n", t);
            fprintf(gp, "splot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf, %lf\n", -R_M * cos(ω * t), -R_M * sin(ω * t), 0.0);
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf, %lf\n", R_m * cos(ω * t), R_m * sin(ω * t), 0.0);
            fprintf(gp, "e\n");
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf, %lf\n", r(0, i), r(1, i), r(2, i));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Y'\n");
            fprintf(gp, "set zlabel 'Z'\n");
            fprintf(gp, "set title 'X-Y-Z t=%lf'\n", t);
            fprintf(gp, "splot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf, %lf\n", -R_M, 0.0, 0.0);
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf, %lf\n", R_m, 0.0, 0.0);
            fprintf(gp, "e\n");
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf, %lf\n", r(0, i) * cos(ω * t) + r(1, i) * sin(ω * t), -r(0, i) * sin(ω * t) + r(1, i) * cos(ω * t), r(2, i));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "unset multiplot\n");
            fflush(gp);

*/

/*

fprintf(gp, "set multiplot layout 1,2\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Y'\n");
            fprintf(gp, "set title 'X-Y t=%lf'\n", t);
            fprintf(gp, "plot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf\n", -R_M * cos(ω * t), -R_M * sin(ω * t));
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf\n", R_m * cos(ω * t), R_m * sin(ω * t));
            fprintf(gp, "e\n");
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf\n", r(0, i), r(1, i));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Z'\n");
            fprintf(gp, "set title 'X-Z t=%lf'\n", t);
            fprintf(gp, "plot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf\n", -R_M * cos(ω * t), 0.0);
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf\n", R_m * cos(ω * t), 0.0);
            fprintf(gp, "e\n");
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf\n", r(0, i), r(2, i));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "unset multiplot\n");
            fflush(gp);

*/

/*
fprintf(gp, "set multiplot layout 1,2\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Y'\n");
            fprintf(gp, "set title 'X-Y t=%lf'\n", t);
            fprintf(gp, "plot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf\n", -R_M, 0.0);
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf\n", R_m, 0.0);
            fprintf(gp, "e\n");
#pragma omp parallel for
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf\n", r(0, i) * cos(ω * t) + r(1, i) * sin(ω * t), -r(0, i) * sin(ω * t) + r(1, i) * cos(ω * t));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Z'\n");
            fprintf(gp, "set title 'X-Z t=%lf'\n", t);
            fprintf(gp, "plot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf\n", -R_M, 0.0);
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf\n", R_m, 0.0);
            fprintf(gp, "e\n");
#pragma omp parallel for
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf\n", r(0, i) * cos(ω * t) + r(1, i) * sin(ω * t), r(2, i));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "unset multiplot\n");
            fflush(gp);
*/

/*
fprintf(gp, "set multiplot layout 1,2\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Y'\n");
            fprintf(gp, "set title 'X-Y t=%lf'\n", t);
            fprintf(gp, "plot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf\n", -R_M * cos(ω * t), -R_M * sin(ω * t));
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf\n", R_m * cos(ω * t), R_m * sin(ω * t));
            fprintf(gp, "e\n");
#pragma omp parallel for
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf\n", r(0, i), r(1, i));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "set xlabel 'X'\n");
            fprintf(gp, "set ylabel 'Z'\n");
            fprintf(gp, "set title 'X-Z t=%lf'\n", t);
            fprintf(gp, "plot '-' pt 7 ps 4 lc \'dark-magenta\', '-' pt 7 ps 2 lc \'web-green\', '-' pt 7 ps 0.6 lc \'web-blue\'\n");
            fprintf(gp, "%lf, %lf\n", -R_M, 0.0);
            fprintf(gp, "e\n");
            fprintf(gp, "%lf, %lf\n", R_m, 0.0);
            fprintf(gp, "e\n");
#pragma omp parallel for
            for (int i = 0; i < k; i++)
            {
                if (tf(i))
                {
                    fprintf(gp, "%lf, %lf\n", r(0, i) * cos(ω * t) + r(1, i) * sin(ω * t), r(1, i) * cos(ω * t) - r(0, i) * sin(ω * t));
                }
            }
            fprintf(gp, "e\n");

            fprintf(gp, "unset multiplot\n");
            fflush(gp);
*/