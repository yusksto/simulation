#include <iostream>

#include <cmath>
#include <math.h>

#include <cstdio>
#include <cassert>

#include "Eigen/Core"

const float dt = 0.05f;
const float t_0 = 0;
const float t_1 = 100000;

const float dr = 3;
const int xrange = 100;
const int yrange = 3;

const float k = 1000;
const float b = 2;
const float m = 10;

inline float cal_a_x(float, float, float, float);
inline float cal_a_y(float, float, float, float);

int main()
{
	FILE* gp = _popen("gnuplot", "w");
	assert(gp);
	fprintf(gp, "unset key\n");
	fprintf(gp, "set xr [%d:%d]\n", (int)(-yrange * dr), (int)(yrange * dr));
	fprintf(gp, "set yr [%d:%d]\n", (int)(-yrange * dr * 0.5f), (int)(yrange * dr * 1.5f));

	Eigen::MatrixXf y_0 = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf y_1 = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf x_0 = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf x_1 = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf vy = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf vx = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf ay = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf ax = Eigen::MatrixXf::Zero(xrange, yrange);

	for (int i = 0; i < xrange; i++)
	{
		for (int n = 0; n < yrange; n++)
		{
			if (n % 2 == 0)
			{
				x_0(i, n) = i * dr;
				x_1(i, n) = i * dr;
				y_0(i, n) = n * dr * sqrt(3.0f) / 2.0f;
				y_1(i, n) = n * dr * sqrt(3.0f) / 2.0f;
			}
			else
			{
				x_0(i, n) = i * dr + 0.5f * dr;
				x_1(i, n) = i * dr + 0.5f * dr;
				y_0(i, n) = n * dr * sqrt(3.0f) / 2.0f;
				y_1(i, n) = n * dr * sqrt(3.0f) / 2.0f;
			}
		}
	}



	for (float t = t_0; t < t_1; t += dt)
	{
		for (int i = 0; i < xrange; i++)
		{
			for (int n = 0; n < yrange; n++)
			{
				ax(i, n) = 0;
				ay(i, n) = -10;
				if (n % 2 == 0)
				{
					if (n != yrange - 1)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i, n + 1), y_0(i, n + 1));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i, n + 1), y_0(i, n + 1));
					}

					if (i != xrange - 1)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i + 1, n), y_0(i + 1, n));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i + 1, n), y_0(i + 1, n));
					}

					if (n != 0)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i, n - 1), y_0(i, n - 1));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i, n - 1), y_0(i, n - 1));
					}

					if (i != 0 && n != 0)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i - 1, n - 1), y_0(i - 1, n - 1));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i - 1, n - 1), y_0(i - 1, n - 1));
					}

					if (i != 0)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i - 1, n), y_0(i - 1, n));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i - 1, n), y_0(i - 1, n));
					}

					if (i != 0 && n != yrange - 1)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i - 1, n + 1), y_0(i - 1, n + 1));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i - 1, n + 1), y_0(i - 1, n + 1));
					}
				}
				else
				{
					if (i != xrange - 1 && n != yrange - 1)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i + 1, n + 1), y_0(i + 1, n + 1));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i + 1, n + 1), y_0(i + 1, n + 1));
					}

					if (i != xrange - 1)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i + 1, n), y_0(i + 1, n));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i + 1, n), y_0(i + 1, n));
					}

					if (i != xrange - 1 && n != 0)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i + 1, n - 1), y_0(i + 1, n - 1));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i + 1, n - 1), y_0(i + 1, n - 1));
					}

					if (n != 0)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i, n - 1), y_0(i, n - 1));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i, n - 1), y_0(i, n - 1));
					}

					if (i != 0)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i - 1, n), y_0(i - 1, n));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i - 1, n), y_0(i - 1, n));
					}

					if (n != yrange - 1)
					{
						ax(i, n) += cal_a_x(x_0(i, n), y_0(i, n), x_0(i, n + 1), y_0(i, n + 1));
						ay(i, n) += cal_a_y(x_0(i, n), y_0(i, n), x_0(i, n + 1), y_0(i, n + 1));
					}
				}

				ax(i, n) -= b * vx(i, n);
				ay(i, n) -= b * vy(i, n);

				vx(i, n) += ax(i, n) * dt;
				vy(i, n) += ay(i, n) * dt;

				x_1(i, n) += vx(i, n) * dt;
				y_1(i, n) += vy(i, n) * dt;
			}
		}

		fprintf(gp, "plot '-' pt 7 ps 0.6 \n");
		for (int i = 1; i < xrange - 1; i++)
		{
			for (int n = 0; n < yrange; n++)
			{

				x_0(i, n) = x_1(i, n);
				y_0(i, n) = y_1(i, n);
				fprintf(gp, "%lf, %lf\n", (float)x_0(i, n), (float)y_0(i, n));
			}
		}
		fprintf(gp, "e\n");
		fflush(gp);
	}






	

	_pclose(gp);

	return 0;
}

inline float cal_a_x(float xa, float ya, float xb, float yb)
{
	return k * (dr - sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb))) * (xa - xb) / sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb)) / m;
}

inline float cal_a_y(float xa, float ya, float xb, float yb)
{
	return k * (dr - sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb))) * (ya - yb) / sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb)) / m;
}