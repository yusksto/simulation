#include <iostream>
#include <cstddef>
#include <memory>
#include <new>
#include <random>
#include <time.h>
#include <cmath>
#include <omp.h>
#include "Eigen/Core"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

struct RGBA {
	unsigned char r, g, b, a; //赤, 緑, 青, 透過
	RGBA() = default;
	constexpr RGBA(const unsigned char r_, const unsigned char g_, const unsigned char b_, const unsigned char a_) :r(r_), g(g_), b(b_), a(a_) {}
};



int main()
{
	constexpr auto width_0 = 1920;//画面幅 //x軸
	constexpr auto height_0 = 1080;//画面高さ //y軸
	const int N = width_0 * height_0;

	const float P = 0;
	const int d = 9;

	std::mt19937 mt;
	mt.seed((unsigned int)time(NULL));
	std::uniform_int_distribution<> uid(0, d);
	uid(mt); uid(mt); uid(mt); uid(mt); uid(mt);

	Eigen::MatrixXi D = Eigen::MatrixXi::Zero(width_0, height_0);//0-上, 1-右, 2-下, 3-左
	Eigen::MatrixXi Q = Eigen::MatrixXi::Zero(width_0, height_0);//0-なし, 1-途中, 2-分岐, 3-末端
	Eigen::MatrixXi J = Eigen::MatrixXi::Zero(width_0, height_0);//0-なし, 1以上- 分岐判定 本数
	Eigen::MatrixXi I = Eigen::MatrixXi::Zero(width_0, height_0);//電流量

	const int midx = (int)(width_0 / 2);
	const int midy = (int)(height_0 / 2);

	Q(midx, midy) = 1;



	//経路作成	
	std::cout << "経路作成 ... " << std::flush;
	for (bool i = true; i;)
	{
		i = false;
#pragma omp parallel for
		for (int s = 0; s < width_0; s++)
		{
			for (int t = 0; t < height_0; t++)
			{
				if (Q(s, t) == 1)
				{
					i = true;
					if (s - 1 >= 0 && Q(s - 1, t) == 0)
					{
						if (uid(mt) == 0)
						{
							Q(s - 1, t) = 1;
							D(s - 1, t) = 1;
							J(s, t)++;
						}
					}
					else
					{
						if (t - 1 >= 0 && Q(s, t - 1) == 0)
						{
							if (uid(mt) == 0)
							{
								Q(s, t - 1) = 1;
								D(s, t - 1) = 0;
								J(s, t)++;
							}
						}
						else
						{
							if (t + 1 < height_0 && Q(s, t + 1) == 0)
							{
								if (uid(mt) == 0)
								{
									Q(s, t + 1) = 1;
									D(s, t + 1) = 2;
									J(s, t)++;
								}
							}
							else
							{
								if (s + 1 < width_0 && Q(s + 1, t) == 0)
								{
									if (uid(mt) == 0)
									{
										Q(s + 1, t) = 1;
										D(s + 1, t) = 3;
										J(s, t)++;
									}
								}
								else
								{
									if (J(s, t) == 0)
									{
										Q(s, t) = 3;//末端										
									}
									else
									{
										Q(s, t) = 2;//分岐
									}
								}
							}
						}
					}
				}
			}
		}
	}
	std::cout << "done" << std::endl;


	//経路トレース
	std::cout << "経路トレース ... " << std::flush;
	for (int s = 0; s < width_0; s++)
	{
		for (int t = 0; t < height_0; t++)
		{
			if (Q(s, t) == 3)
			{
				I(s, t) = 1;
			}
		}
	}


	for (bool i = true; i;)
	{
		i = false;
#pragma omp parallel for
		for (int s = 0; s < width_0; s++)
		{
			for (int t = 0; t < height_0; t++)
			{
				if (Q(s, t) == 3)
				{
					Q(s, t) = 0;
					switch (D(s, t))
					{
					case 0:
						if (t + 1 < height_0)
						{
							i = true;
							I(s, t + 1) += I(s, t);
							J(s, t + 1)--;
							if (J(s, t + 1) == 0)
							{
								Q(s, t + 1) = 3;
							}
						}
						break;
					case 1:
						if (s + 1 < width_0)
						{
							i = true;
							I(s + 1, t) += I(s, t);
							J(s + 1, t)--;
							if (J(s + 1, t) == 0)
							{
								Q(s + 1, t) = 3;
							}
						}
						break;
					case 2:
						if (t - 1 >= 0)
						{
							i = true;
							I(s, t - 1) += I(s, t);
							J(s, t - 1)--;
							if (J(s, t - 1) == 0)
							{
								Q(s, t - 1) = 3;
							}
						}
						break;
					case 3:
						if (s - 1 >= 0)
						{
							i = true;
							I(s - 1, t) += I(s, t);
							J(s - 1, t)--;
							if (J(s - 1, t) == 0)
							{
								Q(s - 1, t) = 3;
							}
						}
						break;
					}
				}
			}
		}
	}
	std::cout << "done" << std::endl;


	//配色処理
	std::cout << "配色処理 ... " << std::flush;
	for (int s = 0; s < width_0; s++)
	{
		for (int t = 0; t < height_0; t++)
		{
			I(s, t) = (int)log(1 + I(s, t));
		}
	}

	int max = 0;
	for (int s = 0; s < width_0; s++)
	{
		for (int t = 0; t < height_0; t++)
		{
			if (I(s, t) > max)
			{
				max = I(s, t);
			}
		}
	}

	for (int s = 0; s < width_0; s++)
	{
		for (int t = 0; t < height_0; t++)
		{
			I(s, t) = (int)(255 * I(s, t) / max);
		}
	}
	std::cout << "done" << std::endl;

	//画像出力
	std::cout << "画像出力 ... " << std::flush;
	constexpr std::size_t width{ width_0 }, height{ height_0 }; //幅と高さ
	std::unique_ptr<RGBA[][width]> rgba(new(std::nothrow) RGBA[height][width]);
	if (!rgba) return -1;

	for (std::size_t row{}; row < height; ++row)
		for (std::size_t col{}; col < width; ++col) {
			rgba[row][col].r = I(col, row);
			rgba[row][col].g = I(col, row);
			rgba[row][col].b = I(col, row);
			rgba[row][col].a = 255; //不透過
		}
	stbi_write_png("picture_1.png", static_cast<int>(width), static_cast<int>(height), static_cast<int>(sizeof(RGBA)), rgba.get(), 0);
	std::cout << "done" << std::endl;
}
