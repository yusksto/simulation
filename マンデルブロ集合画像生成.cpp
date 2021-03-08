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
	constexpr auto width_0 = 4096;//画面幅16384
	constexpr auto height_0 = 4096;//画面高さ
	const long double coordinateX = -0.7; //中心X座標
	const long double coordinateY = 0; //中心Y座標
	const long double Magnification = 0.7; //拡大倍率
	const int times = 100;

	

	constexpr std::size_t width{ width_0 }, height{ height_0 }; //幅と高さ
	std::unique_ptr<RGBA[][width]> rgba(new(std::nothrow) RGBA[height][width]);
	if (!rgba) return -1;

	long double x, y;
	long double x1, y1;
	long double x2, y2;
	int n = 0;
	for (std::size_t row{}; row < height; ++row)
		for (std::size_t col{}; col < width; ++col) {
			x = ((long double)col - (long double)width_0 / 2) / (long double)height_0 * 2 / Magnification + coordinateX;
			y = ((long double)height_0 / 2 - (long double)row) / (long double)height_0 * 2 / Magnification + coordinateY;
			x1 = x;
			y1 = y;
			n = 0;
			for (; n < times; n++) {
				x2 = x1 * x1 - y1 * y1 + x;
				y2 = 2 * x1 * y1 + y;
				x1 = x2;
				y1 = y2;
				if (x2 * x2 + y2 * y2 > 4)
				{
					if (n < times / 2) {
						rgba[row][col].r = (int)(255 * (2 * (float)n / (float)times));
						rgba[row][col].g = (int)(255 * (2 * (float)n / (float)times));
						rgba[row][col].b = 255;
						rgba[row][col].a = 255;
					}
					else {
						rgba[row][col].r = (int)(255 * (2 - 2 * (float)n / (float)times));
						rgba[row][col].g = (int)(255 * (2 - 2 * (float)n / (float)times));
						rgba[row][col].b = (int)(255 * (2 - 2 * (float)n / (float)times));
						rgba[row][col].a = 255;
					}
					break;
				}
			}
			if (n == times)
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

if (n < times / 2) {
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

/*
虹

rgba[row][col].r = 127 + ((n - 6) % 17) * 8;
					rgba[row][col].g = 127 + (n % 17) * 8;
					rgba[row][col].b = 127 + ((n + 6) % 17) * 8;
					rgba[row][col].a = 255;
*/

/*
constexpr auto width_0 = 1960;//画面幅16384
	constexpr auto heigth_0 = 1080;//画面高さ
	const long double coordinateX = -0.75013883465; //中心X座標
	const long double coordinateY = 0.0099859; //中心Y座標
	const long double Magnification = 1000000000000; //拡大倍率
	const int calculations = 15000;
*/