#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <cstdio>
#include <cassert>


int cal();
float cal_a(float, float, float, float, float, float);
int simulate(float, float, float, float);


constexpr auto PI = 3.141592653589793238462643383279502884197F;


//シミュレーションについて
const float dt = 0.001F; //ステップ
const float n = 10.0F;   //何周するか

const float M = 10000.0F; //中心星の質量
const float m = 200.0F;   //惑星の質量
const float G = 200.0F;   //重力定数
const float R_0 = 100.0F; //惑星の回転半径

const float R_1 = R_0 * m / M;          //中心星の回転半径
const float R = R_0 + R_1;              // 中心星と惑星の距離
const float ω = sqrt(G * M / R_0) / R; //惑星の角速度


//gif出力について
const float xy_M = 1.6F * R;         //gif表示範囲
const int dθ = 360;                 //一周当たりの枚数
const int tround = 10;               //gif一周の時間（秒）
const int times = (int)(n * dθ);    //gif枚数

const float dt_1 = 2.0F * PI * n / ω / times; //何秒単位か
const float t_total = 2.0F * PI * n / ω;      //かかる時間
float t_1;                                     //経過時間



//シミュレーションに使用する変数
float t;                      //対象の時間経過
float r;                      //対象の回転半径
float xp_0, yp_0, xp_1, yp_1; //惑星の座標
float xq_0, yq_0, xq_1, yq_1; //中心星の座標
float x_0, y_0, x_1, y_1;     //対象の座標
float vx_0, vy_0, vx_1, vy_1; //対象の速度

float ωt, ωtdt;//あまり関係ない




int main()
{
    FILE* gp = _popen("gnuplot", "w");
    assert(gp);

    fprintf(gp, "set terminal gif animate optimize delay %d size %d,%d\n", (int)(100 * tround / dθ) , (int)(xy_M * 10) , (int)(xy_M * 10));
    fprintf(gp, "set output 'tmp.gif'\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "set xr [%d:%d]\n", (int)(-xy_M), (int)(xy_M));
    fprintf(gp, "set yr [%d:%d]\n", (int)(-xy_M), (int)(xy_M));

    const int p = 360 * 7;

    float* x = new float[p];
    float* y = new float[p];
    float* vx = new float[p];
    float* vy = new float[p];

    for (int i = 0; i <= p - 1; i += 7)
    {
        x[i] = 0.985F * R_0 * cos(2 * PI * i / p);
        y[i] = 0.985F * R_0 * sin(2 * PI * i / p);
        vx[i] = -y[i] * ω;
        vy[i] = x[i] * ω;

        x[i + 1] = 0.99F * R_0 * cos(2 * PI * i / p);
        y[i + 1] = 0.99F * R_0 * sin(2 * PI * i / p);
        vx[i + 1] = -y[i + 1] * ω;
        vy[i + 1] = x[i + 1] * ω;

        x[i + 2] = 0.995F * R_0 * cos(2 * PI * i / p);
        y[i + 2] = 0.995F * R_0 * sin(2 * PI * i / p);
        vx[i + 2] = -y[i + 2] * ω;
        vy[i + 2] = x[i + 2] * ω;

        x[i + 3] = R_0 * cos(2 * PI * i / p);
        y[i + 3] = R_0 * sin(2 * PI * i / p);
        vx[i + 3] = -y[i + 3] * ω;
        vy[i + 3] = x[i + 3] * ω;

        x[i + 4] = 1.005F * R_0 * cos(2 * PI * i / p);
        y[i + 4] = 1.005F * R_0 * sin(2 * PI * i / p);
        vx[i + 4] = -y[i + 4] * ω;
        vy[i + 4] = x[i + 4] * ω;

        x[i + 5] = 1.01F * R_0 * cos(2 * PI * i / p);
        y[i + 5] = 1.01F * R_0 * sin(2 * PI * i / p);
        vx[i + 5] = -y[i + 5] * ω;
        vy[i + 5] = x[i + 5] * ω;

        x[i + 6] = 1.015F * R_0 * cos(2 * PI * i / p);
        y[i + 6] = 1.015F * R_0 * sin(2 * PI * i / p);
        vx[i + 6] = -y[i + 6] * ω;
        vy[i + 6] = x[i + 6] * ω;
    }

    for (int i = 0; i < times; i++)
    {
        t_1 = i * dt_1;

        std::cout << (int)(100.0F * t_1 / t_total) << std::endl;

        fprintf(gp, "plot '-' pt 7 ps 4, '-' pt 7 ps 2, '-' pt 7 ps 0.6\n");
        fprintf(gp, "%d, %d\n", (int)(-R_1 * cos(ω * t_1)), (int)(-R_1 * sin(ω * t_1)));
        fprintf(gp, "e\n");
        fprintf(gp, "%d, %d\n", (int)(R_0 * cos(ω * t_1)), (int)(R_0 * sin(ω * t_1)));
        fprintf(gp, "e\n");

        for (int j = 0; j < p - 1; j++)
        {
            simulate(x[j], y[j], vx[j], vy[j]);
            x[j] = x_0;
            y[j] = y_0;
            vx[j] = vx_0;
            vy[j] = vy_0;
        }

        for (int j = 0; j < p - 1; j++)
        {
            fprintf(gp, "%d, %d\n", (int)x[j], (int)y[j]);
        }
        fprintf(gp, "e\n");
    }
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;

    _pclose(gp);

    return 0;
}

int simulate(float x, float y, float vx, float vy)
{

    t = 0;

    x_0 = x;
    y_0 = y;

    vx_0 = vx;
    vy_0 = vy;

    for (; t < dt_1; t += dt)
    {
        ωt = ω * (t + t_1);
        ωtdt = ω * (t + +t_1 + dt);

        xp_0 = R_0 * cos(ωt);
        yp_0 = R_0 * sin(ωt);
        xp_1 = R_0 * cos(ωtdt);
        yp_1 = R_0 * sin(ωtdt);

        xq_0 = -R_1 * cos(ωt);
        yq_0 = -R_1 * sin(ωt);
        xq_1 = -R_1 * cos(ωtdt);
        yq_1 = -R_1 * sin(ωtdt);

        cal();
    }
    return 0;
}

int cal()
{
    x_1 = x_0 + (vx_0 + cal_a(x_0, y_0, xp_0, yp_0, xq_0, yq_0) * dt) * dt;
    y_1 = y_0 + (vy_0 + cal_a(y_0, x_0, yp_0, xp_0, yq_0, xq_0) * dt) * dt;

    vx_1 = vx_0 + (cal_a(x_0, y_0, xp_0, yp_0, xq_0, yq_0) + cal_a(x_1, y_1, xp_1, yp_1, xq_1, yq_1)) / 2.0F * dt;
    vy_1 = vy_0 + (cal_a(y_0, x_0, yp_0, xp_0, yq_0, xq_0) + cal_a(y_1, x_1, yp_1, xp_1, yq_1, xq_1)) / 2.0F * dt;

    x_1 = x_0 + (vx_1 + vx_0) / 2.0F * dt;
    y_1 = y_0 + (vy_1 + vy_0) / 2.0F * dt;

    vx_0 = vx_1;
    vy_0 = vy_1;

    x_0 = x_1;
    y_0 = y_1;

    return 0;
}


float cal_a(float x1, float y1, float x2, float y2, float x3, float y3)
{
    return -G * M * (x1 - x3) / ((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3)) / sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3)) - G * m * (x1 - x2) / ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}
