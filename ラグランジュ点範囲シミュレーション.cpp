#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>


int cul();
float cul_a(float, float, float, float, float, float);
bool simulate(float, float);


const float PI = 3.141589793238462643383279502884197F;


const float dt = 0.01F; //ステップ
const float dxy = 1.0F; //シミュレーションする座標の細かさ
const float n = 100.0F; //周回数


const float M = 10000.0F; //中心星の質量
const float m = 200.0F; //惑星の質量
const float G = 200.0F; //重力定数
const float R_0 = 100.0F; //惑星の回転半径
const float dR = 0.1F; //探索範囲


const float R_1 = R_0 * m / M; //中心星の回転半径
const float R = R_0 + R_1; // 中心星と惑星の距離
const float ω = sqrt(G * M / R_0) / R; //惑星の角速度


const float xy_m = 0.70710628F * dR * R;
const float xy_M = (1.0F + dR) * R;
const int p = (int)(xy_M * xy_M * 4 / dxy / dxy);



float t_0, t_1;
float r;
float xp_0, yp_0, xp_1, yp_1; //惑星の座標
float xq_0, yq_0, xq_1, yq_1; //中心星の座標
float x_0, y_0, x_1, y_1;
float vx_0, vy_0, vx_1, vy_1;


int main()
{
    std::ofstream ofs;
    ofs.open("data.txt", std::ios::out);
    ofs << 0 << " " << 0 << std::endl;
    ofs << R_0 << " " << 0 << std::endl;
    ofs << -R_1 << " " << 0 << std::endl;
    std::cout << 0 << " " << 0 << std::endl;
    std::cout << R_0 << " " << 0 << std::endl;
    std::cout << -R_1 << " " << 0 << std::endl;

    int q = 0;
    int r = 0;

    for (float y = -xy_M; y < xy_M; y += dxy)
    {
        for (float x = -xy_M; x < xy_M; x += dxy)
        {
            if (abs(sqrt(x * x + y * y) - R) < dR * R)
            {
                if (simulate(x, y) == true) //シミュレーション
                {
                    ofs << x << " " << y << std::endl;
                }
            }

            q++; //進行度表示
            if (r + 1 < 100 * q / p)
            {
                r = (int)(100 * q / p) - 1;
                std::cout << r << std::endl;
            }
        }
    }
    return 0;
}

bool simulate(float x, float y)
{

    t_0 = 0;

    r = sqrt(x * x + y * y);

    x_0 = x;
    y_0 = y;

    vx_0 = -y / r * sqrt(G * (M + m) / r);
    vy_0 = x / r * sqrt(G * (M + m) / r);

    for (; ω * t_0 < n * 2 * PI; t_0 += dt)
    {
        t_1 = t_0 + dt;
        xp_0 = R_0 * cos(ω * t_0);
        yp_0 = R_0 * sin(ω * t_0);
        xp_1 = R_0 * cos(ω * t_1);
        yp_1 = R_0 * sin(ω * t_1);

        xq_0 = -R_1 * cos(ω * t_0);
        yq_0 = -R_1 * sin(ω * t_0);
        xq_1 = -R_1 * cos(ω * t_1);
        yq_1 = -R_1 * sin(ω * t_1);

        cul();

        if (sqrt((x_0 - (x * cos(ω * t_1) - y * sin(ω * t_1))) * (x_0 - (x * cos(ω * t_1) - y * sin(ω * t_1))) + (y_0 - (x * sin(ω * t_1) + y * cos(ω * t_1))) * (y_0 - (x * sin(ω * t_1) + y * cos(ω * t_1)))) > R / 2.0F)
        {
            return false;
        }
    }
    return true;
}

int cul()
{
    x_1 = x_0 + (vx_0 + cul_a(x_0, y_0, xp_0, yp_0, xq_0, yq_0) * dt) * dt;
    y_1 = y_0 + (vy_0 + cul_a(y_0, x_0, yp_0, xp_0, yq_0, xq_0) * dt) * dt;

    vx_1 = vx_0 + (cul_a(x_0, y_0, xp_0, yp_0, xq_0, yq_0) + cul_a(x_1, y_1, xp_1, yp_1, xq_1, yq_1)) / 2 * dt;
    vy_1 = vy_0 + (cul_a(y_0, x_0, yp_0, xp_0, yq_0, xq_0) + cul_a(y_1, x_1, yp_1, xp_1, yq_1, xq_1)) / 2 * dt;

    x_1 = x_0 + (vx_1 + vx_0) / 2 * dt;
    y_1 = y_0 + (vy_1 + vy_0) / 2 * dt;

    vx_0 = vx_1;
    vy_0 = vy_1;

    x_0 = x_1;
    y_0 = y_1;

    return 0;
}


float cul_a(float x1, float y1, float x2, float y2, float x3, float y3)
{
    float a = -G * M * (x1 - x3) / ((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3)) / sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3)) - G * m * (x1 - x2) / ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    return a;
}


