#include <iostream>
#include <cstddef>
#include <memory>
#include <new>
#include <random>
#include <time.h>
#include <cmath>
#include <math.h>
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
	constexpr auto width_0 = 1080;//画面幅 //x軸
	constexpr auto height_0 = 1080;//画面高さ //y軸
	constexpr auto times = 10000;//凝集個数

	std::mt19937 mt;
	mt.seed((unsigned int)time(NULL));
	std::uniform_int_distribution<> uid_0(0, 10);//初期密度
	std::uniform_int_distribution<> uid_1(0, 5);//0～3 - 上下左右, 4以上 - そのまま
	std::uniform_int_distribution<> uid_width(0, width_0 - 1);
	std::uniform_int_distribution<> uid_height(0, height_0 - 1);
	mt; mt; mt; mt; mt;

	Eigen::MatrixXi N = Eigen::MatrixXi::Zero(width_0, height_0);

	for (int i = 0; i < width_0; i++)
	{
		for (int n = 0; n < height_0; n++)
		{
			if (uid_0(mt) == 0)
			{
				N(i, n) = 1;
			}
		}
	}

	const int midx = (int)(width_0 / 2);
	const int midy = (int)(height_0 / 2);

	N(midx, midy) = 2;

	//画像生成
	std::cout << "画像生成 ... " << std::flush;
	int width_1, height_1;
	for (int t = 0; t < times;)
	{
		for (int i = 0; i < width_0; i++)
		{
			for (int n = 0; n < height_0; n++)
			{
				if (N(i, n) == 1)
				{
					width_1 = uid_width(mt);
					height_1 = uid_height(mt);
					if (i != 0 && N(i - 1, n) == 2)
					{
						N(i, n) = 2;
						t++;
						std::cout << t << std::endl;
					}
					else
					{
						if (i != width_0 - 1 && N(i + 1, n) == 2)
						{
							N(i, n) = 2;
							t++;
							std::cout << t << std::endl;
						}
						else
						{
							if (n != 0 && N(i, n - 1) == 2)
							{
								N(i, n) = 2;
								t++;
								std::cout << t << std::endl;
							}
							else
							{
								if (n != height_0 - 1 && N(i, n + 1) == 2)
								{
									N(i, n) = 2;
									t++;
									std::cout << t << std::endl;
								}
								else
								{
									switch (uid_1(mt))
									{
									case 4:
										if (i != 0 && N(i - 1, n) == 0)
										{
											N(i - 1, n) = 1;
											N(i, n) = 0;
										}
										else
										{
											if (i == 0 && N(width_0 - 1, n) == 0)
											{
												N(width_0 - 1, n) = 1;
												N(i, n) = 0;
											}
										}
										break;
									case 1:
										if (n != 0 && N(i, n - 1) == 0)
										{
											N(i, n - 1) = 1;
											N(i, n) = 0;
										}
										else
										{
											if (n == 0 && N(i, height_0 - 1) == 0)
											{
												N(i, height_0 - 1) = 1;
												N(i, n) = 0;
											}
										}
										break;
									case 2:
										if (i != width_0 - 1 && N(i + 1, n) == 0)
										{
											N(i + 1, n) = 1;
											N(i, n) = 0;
										}
										else
										{
											if (i == width_0 - 1 && N(0, n) == 0)
											{
												N(0, n) = 1;
												N(i, n) = 0;
											}
										}
										break;
									case 3:
										if (n != height_0 - 1 && N(i, n + 1) == 0)
										{
											N(i, n + 1) = 1;
											N(i, n) = 0;
										}
										else
										{
											if (n == height_0 - 1 && N(i, 0) == 0)
											{
												N(i, 0) = 1;
												N(i, n) = 0;
											}
										}
										break;
									default:
										break;
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



	//画像出力
	std::cout << "画像出力 ... " << std::flush;
	constexpr std::size_t width{ width_0 }, height{ height_0 }; //幅と高さ
	std::unique_ptr<RGBA[][width]> rgba(new(std::nothrow) RGBA[height][width]);
	if (!rgba) return -1;

	for (std::size_t row{}; row < height; ++row)
		for (std::size_t col{}; col < width; ++col) {
			if (N(col, row) == 2)
			{
				rgba[row][col].r = 0;
				rgba[row][col].g = 0;
				rgba[row][col].b = 255;
				rgba[row][col].a = 255;
			}
			else
			{
				if (N(col, row) == 1)
				{
					rgba[row][col].r = 255;
					rgba[row][col].g = 0;
					rgba[row][col].b = 0;
					rgba[row][col].a = 255;
				}
				else
				{
					rgba[row][col].r = 255;
					rgba[row][col].g = 255;
					rgba[row][col].b = 255;
					rgba[row][col].a = 255;
				}
			}
		}
	stbi_write_png("picture_1.png", static_cast<int>(width), static_cast<int>(height), static_cast<int>(sizeof(RGBA)), rgba.get(), 0);
	std::cout << "done" << std::endl;
}
