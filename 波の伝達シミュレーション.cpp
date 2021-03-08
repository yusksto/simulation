#include <iostream>
#include <cmath>
#include <math.h>
#include <cstdio>
#include <cassert>
#include "Eigen/Core"

const float dt = 0.05f;
const float t_0 = 0;
const float t_1 = 1000;

const int xrange = 201;
const int yrange = 201;
const float dxy = 1;
const int N = xrange * yrange;

const float k = 10;
const float b = 0.01f;
const float m = 1;


inline float cal_a(float);

int main()
{
	FILE* gp = _popen("gnuplot", "w");
	assert(gp);
	fprintf(gp, "unset key\n");
	fprintf(gp, "set xr [%d:%d]\n", (int)(-xrange * 0.5f), (int)(xrange * 1.5f));
	fprintf(gp, "set yr [%d:%d]\n", (int)(-yrange * 0.5f), (int)(yrange * 1.5f));
	fprintf(gp, "set zr [-30:30]\n");
	const float n = 9 * 100 * m / (9 - N);

	Eigen::MatrixXf y_0 = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf y_1 = Eigen::MatrixXf::Zero(xrange, yrange);
	Eigen::MatrixXf vy = Eigen::MatrixXf::Constant(xrange, yrange, n);
	Eigen::MatrixXf ay = Eigen::MatrixXf::Zero(xrange, yrange);


	const int midx = (int)(xrange / 2);
	const int midy = (int)(yrange / 2);
	

	vy(midx, midy) = 100;
	vy(midx + 1, midy) = 100;
	vy(midx - 1, midy) = 100;
	vy(midx, midy + 1) = 100;
	vy(midx, midy + 1) = 100;
	vy(midx + 1, midy + 1) = 100;
	vy(midx - 1, midy - 1) = 100;
	vy(midx + 1, midy - 1) = 100;
	vy(midx - 1, midy + 1) = 100;

	


	//シミュレーション開始
	for (float t = t_0; t < t_1; t += dt)
	{
		fprintf(gp, "splot \"-\" with lines\n");

		for (int i = 0; i < xrange; i++)
		{
			for (int n = 0; n < yrange; n++)
			{
				ay(i, n) = 0;

				if (i - 1 >= 0)
				{
					ay(i, n) += cal_a(y_0(i, n) - y_0(i - 1, n));
				}

				if (i + 1 < xrange)
				{
					ay(i, n) += cal_a(y_0(i, n) - y_0(i + 1, n));
				}

				if (n - 1 >= 0)
				{
					ay(i, n) += cal_a(y_0(i, n) - y_0(i, n - 1));
				}

				if (n + 1 < yrange)
				{
					ay(i, n) += cal_a(y_0(i, n) - y_0(i, n + 1));
				}

				ay(i, n) -= b * vy(i, n);

				vy(i, n) += ay(i, n) * dt;
				y_1(i, n) += vy(i, n) * dt;

				fprintf(gp, "%lf, %lf, %lf\n", (float)i, (float)n, (float)y_1(i, n));
			}

			if (i == xrange - 1)
			{
				break;
			}

			i++;

			for (int n = yrange - 1; n >= 0; n--)
			{
				ay(i, n) = 0;

				if (i - 1 >= 0)
				{
					ay(i, n) += cal_a(y_0(i, n) - y_0(i - 1, n));
				}

				if (i + 1 < xrange)
				{
					ay(i, n) += cal_a(y_0(i, n) - y_0(i + 1, n));
				}

				if (n - 1 >= 0)
				{
					ay(i, n) += cal_a(y_0(i, n) - y_0(i, n - 1));
				}

				if (n + 1 < yrange)
				{
					ay(i, n) += cal_a(y_0(i, n) - y_0(i, n + 1));
				}

				ay(i, n) -= b * vy(i, n);

				vy(i, n) += ay(i, n) * dt;
				y_1(i, n) += vy(i, n) * dt;

				fprintf(gp, "%lf, %lf, %lf\n", (float)i, (float)n, (float)y_1(i, n));
			}

			

		}
		fprintf(gp, "e\n");
		for (int i = 0; i < xrange; i++)
		{
			for (int n = 0; n < yrange; n++)
			{
				y_0(i, n) = y_1(i, n);
			}
		}
		fflush(gp);
	}

	_pclose(gp);
	return 0;
}

inline float cal_a(float dy)
{
	return - k * (sqrt(dy * dy + dxy * dxy) - dxy) * dy / sqrt(dy * dy + dxy * dxy) / m;
}






