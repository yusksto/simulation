#include <iostream>
#include <fstream>
#include "Eigen/Core"
using namespace Eigen;

Vector3f RungeKutta(float t, Vector3f X);
Vector3f LorentzEquation(float t, Vector3f X);

const float p = 10;
const float r = 28;
const float b = 8 / 3;
const float dt = 0.001;
const float t_0 = 0;
const float t_1 = 50;
const Vector3f X_0(1, 1, 1);

int main() {
    float t = t_0;
    Vector3f X = X_0;

    std::ofstream ofs;
    ofs.open("Lorentz_attractor.dat", std::ios::out);
    ofs << t << " " << X.transpose() << std::endl;

    do {
        X = RungeKutta(t, X);
        t += dt;
        ofs << t << " " << X.transpose() << std::endl;
    } while (t < t_1);

    return 0;
}

Vector3f RungeKutta(float t, Vector3f X) {
    Vector3f k_1 = LorentzEquation(t, X);
    Vector3f k_2 = LorentzEquation(t + dt / 2, X + dt / 2 * k_1);
    Vector3f k_3 = LorentzEquation(t + dt / 2, X + dt / 2 * k_2);
    Vector3f k_4 = LorentzEquation(t + dt, X + dt * k_3);
    Vector3f X_next = X + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);

    return X_next;
}
Vector3f LorentzEquation(float t, Vector3f X) {

    Vector3f val(-p * X[0] + p * X[1], -X[0] * X[2] + r * X[0] - X[1], X[0] * X[1] - b * X[2]);

    return val;
}