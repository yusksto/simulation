#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>

const double t_0 = 0;
const double t_1 = 25000;
const double dt = 0.001;

const double GM = 20;
const double a = 0;

const double x_0 = 20;
const double y_0 = 0;
const double vx_0 = 0;
const double vy_0 = sqrt(GM * 2 / x_0);

long double x[2], y[2], vx[2], vy[2];
long double ax[2], ay[2];

const int n = (int)(0.1 / dt + 1);

int cul();

int main()
{
    std::ofstream ofs;
    ofs.open("data.txt", std::ios::out);
    x[0] = x_0;
    y[0] = y_0;
    x[1] = x_0;
    y[1] = y_0;
    vx[0] = vx_0;
    vy[0] = vy_0;
    int m = n;
    for (double t = t_0; t < t_1; t += dt) {
        if (m == n)
        {
            ofs << t << " " << x[0] << " " << y[0] << std::endl;
            m = 1;
        }
        m++;
        cul();
    }
    ofs.close();
    return 0;
}

int cul() {
    ax[0] = -GM * ((x[0] - a) / ((x[0] - a) * (x[0] - a) + y[0] * y[0]) / sqrt((x[0] - a) * (x[0] - a) + y[0] * y[0]) + (x[0] + a) / ((x[0] + a) * (x[0] + a) + y[0] * y[0]) / sqrt((x[0] + a) * (x[0] + a) + y[0] * y[0]));
    ay[0] = -GM * (y[0] / ((x[0] - a) * (x[0] - a) + y[0] * y[0]) / sqrt((x[0] - a) * (x[0] - a) + y[0] * y[0]) + y[0] / ((x[0] + a) * (x[0] + a) + y[0] * y[0]) / sqrt((x[0] + a) * (x[0] + a) + y[0] * y[0]));

    x[1] = x[0] + (vx[0] + ax[0] * dt) * dt;
    y[1] = y[0] + (vy[0] + ay[0] * dt) * dt;


    ax[0] = -GM * ((x[0] - a) / ((x[0] - a) * (x[0] - a) + y[0] * y[0]) / sqrt((x[0] - a) * (x[0] - a) + y[0] * y[0]) + (x[0] + a) / ((x[0] + a) * (x[0] + a) + y[0] * y[0]) / sqrt((x[0] + a) * (x[0] + a) + y[0] * y[0]));
    ay[0] = -GM * (y[0] / ((x[0] - a) * (x[0] - a) + y[0] * y[0]) / sqrt((x[0] - a) * (x[0] - a) + y[0] * y[0]) + y[0] / ((x[0] + a) * (x[0] + a) + y[0] * y[0]) / sqrt((x[0] + a) * (x[0] + a) + y[0] * y[0]));
    ax[1] = -GM * ((x[1] - a) / ((x[1] - a) * (x[1] - a) + y[1] * y[1]) / sqrt((x[1] - a) * (x[1] - a) + y[1] * y[1]) + (x[1] + a) / ((x[1] + a) * (x[1] + a) + y[1] * y[1]) / sqrt((x[1] + a) * (x[1] + a) + y[1] * y[1]));
    ay[1] = -GM * (y[1] / ((x[1] - a) * (x[1] - a) + y[1] * y[1]) / sqrt((x[1] - a) * (x[1] - a) + y[1] * y[1]) + y[1] / ((x[1] + a) * (x[1] + a) + y[1] * y[1]) / sqrt((x[1] + a) * (x[1] + a) + y[1] * y[1]));

    vx[1] = vx[0] + (ax[0] + ax[1]) / 2 * dt;
    vy[1] = vy[0] + (ay[0] + ay[1]) / 2 * dt;

    x[1] = x[0] + (vx[1] + vx[0]) / 2 * dt;
    y[1] = y[0] + (vy[1] + vy[0]) / 2 * dt;

    x[0] = x[1];
    y[0] = y[1];

    vx[0] = vx[1];
    vy[0] = vy[1];
    return 0;
}
