#include <iostream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <cstddef>
#include <memory>
#include <new>

struct RGBA {
	unsigned char r, g, b, a; //赤, 緑, 青, 透過
	RGBA() = default;
	constexpr RGBA(const unsigned char r_, const unsigned char g_, const unsigned char b_, const unsigned char a_) :r(r_), g(g_), b(b_), a(a_) {}
};

int main()
{
	constexpr auto width_0 = 4096;//画面幅
	constexpr auto heigth_0 = 2048;//画面高さ
	const long double coordinateX = 0; //中心X座標
	const long double coordinateY = 0; //中心Y座標
	const long double Magnification = 1; //拡大倍率
	const int calculations = 500;

	//定数 a + bi
	const long double a = -0.74543;
	const long double b = 0.11301;


	constexpr std::size_t width{ width_0 }, height{ heigth_0 }; //幅と高さ
	std::unique_ptr<RGBA[][width]> rgba(new(std::nothrow) RGBA[height][width]);
	if (!rgba) return -1;

	long double x, y;
	long double x1, x2, y1, y2;
	int n = 0;
	for (std::size_t row{}; row < height; ++row)
		for (std::size_t col{}; col < width; ++col) {
			x = ((long double)col - (long double)width_0 / 2) / (long double)heigth_0 * 2 / Magnification + coordinateX;
			y = ((long double)heigth_0 / 2 - (long double)row) / (long double)heigth_0 * 2 / Magnification + coordinateY;
			x1 = x;
			y1 = y;
			n = 0;
			for (; n < calculations; n++) {
				x2 = x1 * x1 - y1 * y1 + a;
				y2 = 2 * x1 * y1 + b;
				x1 = x2;
				y1 = y2;
				if (x2 * x2 + y2 * y2 > 2)
				{
					if (n < calculations / 2) {
						rgba[row][col].r = (int)(255 * (2 * (float)n / (float)calculations));
						rgba[row][col].g = (int)(255 * (2 * (float)n / (float)calculations));
						rgba[row][col].b = 255;
						rgba[row][col].a = 255;
					}
					else {
						rgba[row][col].r = (int)(255 * (2 - 2 * (float)n / (float)calculations));
						rgba[row][col].g = (int)(255 * (2 - 2 * (float)n / (float)calculations));
						rgba[row][col].b = (int)(255 * (2 - 2 * (float)n / (float)calculations));
						rgba[row][col].a = 255;
					}
					break;
				}
			}
			if (n == calculations)
			{
				rgba[row][col].r = 0;
				rgba[row][col].g = 0;
				rgba[row][col].b = 0;
				rgba[row][col].a = 255;
			}
		}
	stbi_write_png("picture_1.png", static_cast<int>(width), static_cast<int>(height), static_cast<int>(sizeof(RGBA)), rgba.get(), 0);
}

/*
青 > 白 > 黒

if (n < calculations / 2) {
						rgba[row][col].r = (int)(255 * (2 * (float)n / (float)calculations));
						rgba[row][col].g = (int)(255 * (2 * (float)n / (float)calculations));
						rgba[row][col].b = 255;
						rgba[row][col].a = 255;
					}
					else {
						rgba[row][col].r = (int)(255 * (2 - 2 * (float)n / (float)calculations));
						rgba[row][col].g = (int)(255 * (2 - 2 * (float)n / (float)calculations));
						rgba[row][col].b = (int)(255 * (2 - 2 * (float)n / (float)calculations));
						rgba[row][col].a = 255;
					}

*/
