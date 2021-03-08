#include <iostream>
#include <cmath>
#include <math.h>
#include <random>
#include <cstdio>
#include <cassert>

template<typename A, size_t N, typename T>
void Fill(A(&array)[N], const T& val) {
	std::fill((T*)array, (T*)(array + N), val);
}

//シグモイド関数
inline float func(float);

constexpr float PI = 3.14159265358979;

int main()
{
	//gnuplotパイプ
	FILE* gp = _popen("gnuplot", "w");
	assert(gp);

	//乱数
	std::random_device rd;
	rd; rd; rd; rd; rd;
	std::mt19937 mt;
	mt.seed(rd());
	int width = 10;  //乱数範囲[-width , width]
	int digit = 10000;//乱数小数点精度（digit分の１まで）
	std::uniform_int_distribution<> uid_0(0, 2 * width * digit);
	std::uniform_int_distribution<> uid_1(0, 500);

	//使用配列
	float weight_0[9][10][10], weight_1[9][10][10];
	float weight_00[2][10], weight_01[10], weight_10[2][10], weight_11[10];
	float bias_0[10][10], bias_1[10][10];
	float network[10][10];

	//配列初期化
	Fill(weight_0, 0);
	Fill(weight_1, 0);
	Fill(weight_00, 0);
	Fill(weight_01, 0);
	Fill(weight_10, 0);
	Fill(weight_11, 0);
	Fill(bias_0, 0);
	Fill(bias_1, 0);
	Fill(network, 0);

	//乱数によるパラメーターの初期化
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 10; j++) {
			for (int k = 0; k < 10; k++) {
				weight_0[i][j][k] = (float)uid_0(mt) / (float)digit - (float)width;
			}
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 10; j++) {
			weight_00[i][j] = (float)uid_0(mt) / (float)digit - (float)width;
		}
	}
	for (int i = 0; i < 10; i++) {
		weight_01[i] = (float)uid_0(mt) / (float)digit - (float)width;
	}
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			bias_0[i][j] = (float)uid_0(mt) / (float)digit - (float)width;
		}
	}

	//パラメーラーコピー
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 10; j++) {
			for (int k = 0; k < 10; k++) {
				weight_1[i][j][k] = weight_0[i][j][k];
			}
		}
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 10; j++) {
			weight_10[i][j] = weight_00[i][j];
		}
	}
	for (int i = 0; i < 10; i++) {
		weight_11[i] = weight_01[i];
	}
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			bias_1[i][j] = bias_0[i][j];
		}
	}

	//開始
	const int n_0 = 100;                 //最大変異回数
	int n = 0;                            //変異回数カウント
	const int R = 1000;                 //リセット回数
	int r = 0;                            //リセット回数カウント
	int t = 1;                            //総リセット回数カウント
	float δ = 99999999999999, δ_0 = 99999999999999, δ_1 = 0;//誤差
	const int F = 1;                      //画像出力頻度
	int f = 0;                            //画像出力カウント用
	bool inA = 0;
	bool inB = 0;
	float outR = 0;
	bool outT = 0;
	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set xrange [-1:4]\n");
	fprintf(gp, "set yrange [-1:2]\n");
	while (n < n_0) {
		r++;
		//変数初期化
		δ_1 = 0;
		outR = 0;
		Fill(network, 0);

		//誤差計算
		for (int m = 0; m < 4; m++) {
			switch (m)
			{
			case(0):
				inA = 0;
				inB = 0;
				break;
			case(1):
				inA = 1;
				inB = 0;
				break;
			case(2):
				inA = 0;
				inB = 1;
				break;
			case(3):
				inA = 1;
				inB = 1;
				break;
			}

			//ニューラルネットワーク計算
			for (int i = 0; i < 10; i++) {
				network[0][i] += (float)inA * weight_10[0][i];
			}
			for (int i = 0; i < 10; i++) {
				network[0][i] += (float)inB * weight_10[1][i];
			}
			for (int i = 0; i < 10; i++) {
				network[0][i] = func(network[0][i] + bias_1[0][i]);
			}
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 10; j++) {
					for (int k = 0; k < 10; k++) {
						network[i + 1][j] += network[i][k] * weight_1[i][j][k];
					}
					network[i + 1][j] = func(network[i + 1][j] + bias_1[i + 1][j]);
				}
			}
			for (int i = 0; i < 10; i++) {
				outR += network[9][i] * weight_11[i];
			}
			outR = func(outR);

			if ((int)inA - (int)inB == 0){
				outT = 0;
			} else {
				outT = 1;
			}

			δ_1 += pow(abs(outR - outT) + 1, 2);
		}

		//比較、表示、上書き
		if (δ_1 < δ_0 && δ_1 < δ) {
			δ_0 = δ_1;
			n++;
			r = 0;

			//表示（F回に１回）
			f++;
			if (F == f) {
				f = 0;
				std::cout << n << " " << δ_0 << " " << t << std::endl;
				fprintf(gp, "set output '%d.png'\n", t * 1000000 + n);
				fprintf(gp, "plot '-' with lines linetype rgbcolor '#0000FF', '-' with lines linetype rgbcolor '#FF0000'\n");

				for (int m = 0; m < 4; m++) {
					switch (m)
					{
					case(0):
						inA = 0;
						inB = 0;
						break;
					case(1):
						inA = 1;
						inB = 0;
						break;
					case(2):
						inA = 0;
						inB = 1;
						break;
					case(3):
						inA = 1;
						inB = 1;
						break;
					}

					//ニューラルネットワーク計算
					for (int i = 0; i < 10; i++) {
						network[0][i] += (float)inA * weight_10[0][i];
					}
					for (int i = 0; i < 10; i++) {
						network[0][i] += (float)inB * weight_10[1][i];
					}
					for (int i = 0; i < 10; i++) {
						network[0][i] = func(network[0][i] + bias_1[0][i]);
					}
					for (int i = 0; i < 9; i++) {
						for (int j = 0; j < 10; j++) {
							for (int k = 0; k < 10; k++) {
								network[i + 1][j] += network[i][k] * weight_1[i][j][k];
							}
							network[i + 1][j] = func(network[i + 1][j] + bias_1[i + 1][j]);
						}
					}
					for (int i = 0; i < 10; i++) {
						outR += network[9][i] * weight_11[i];
					}
					outR = func(outR);

					if ((int)inA - (int)inB == 0) {
						outT = 0;
					}
					else {
						outT = 1;
					}

					fprintf(gp, "%d, %lf\n", m, outR);
				}
				fprintf(gp, "e\n");
				fflush(gp);
				for (int m = 0; m < 4; m++) {
					switch (m)
					{
					case(0):
						inA = 0;
						inB = 0;
						break;
					case(1):
						inA = 1;
						inB = 0;
						break;
					case(2):
						inA = 0;
						inB = 1;
						break;
					case(3):
						inA = 1;
						inB = 1;
						break;
					}
					if ((int)inA - (int)inB == 0) {
						outT = 0;
					}
					else {
						outT = 1;
					}
					fprintf(gp, "%d, %d\n", m, outT);
				}
				fprintf(gp, "e\n");
				fflush(gp);
			}

			//上書き
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 10; j++) {
					for (int k = 0; k < 10; k++) {
						weight_0[i][j][k] = weight_1[i][j][k];
					}
				}
			}
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 10; j++) {
					weight_00[i][j] = weight_10[i][j];
				}
			}
			for (int i = 0; i < 10; i++) {
				weight_01[i] = weight_11[i];
			}
			for (int i = 0; i < 10; i++) {
				for (int j = 0; j < 10; j++) {
					bias_0[i][j] = bias_1[i][j];
				}
			}

		}

		//リセット判定
		if (r > R) {
			//リセット
			r = 0;
			n = 0;
			if (δ > δ_0) {
				δ = δ_0;
			}
			
			δ_0 = 99999999999999;
			t++;
			//乱数によるパラメーターの初期化
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 10; j++) {
					for (int k = 0; k < 10; k++) {
						weight_0[i][j][k] = (float)uid_0(mt) / (float)digit - (float)width;
					}
				}
			}
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 10; j++) {
					weight_00[i][j] = (float)uid_0(mt) / (float)digit - (float)width;
				}
			}
			for (int i = 0; i < 10; i++) {
				weight_01[i] = (float)uid_0(mt) / (float)digit - (float)width;
			}
			for (int i = 0; i < 10; i++) {
				for (int j = 0; j < 10; j++) {
					bias_0[i][j] = (float)uid_0(mt) / (float)digit - (float)width;
				}
			}
		}

		//変異
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 10; j++) {
				for (int k = 0; k < 10; k++) {
					if (uid_1(mt) == 0) {
						weight_1[i][j][k] = weight_0[i][j][k] + (float)uid_0(mt) / (float)digit - (float)width;
					}
					else {
						if (uid_1(mt) == 0) {
							weight_1[i][j][k] = (float)uid_0(mt) / (float)digit - (float)width;
						}
					}
				}
			}
		}
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 10; j++) {
				if (uid_1(mt) == 0) {
					weight_10[i][j] = weight_00[i][j] + (float)uid_0(mt) / (float)digit - (float)width;
				}
				else {
					if (uid_1(mt) == 0) {
						weight_10[i][j] = (float)uid_0(mt) / (float)digit - (float)width;
					}
				}
			}
		}
		for (int i = 0; i < 10; i++) {
			if (uid_1(mt) == 0) {
				weight_11[i] = weight_01[i] + (float)uid_0(mt) / (float)digit - (float)width;
			}
			else {
				if (uid_1(mt) == 0) {
					weight_11[i] = (float)uid_0(mt) / (float)digit - (float)width;
				}
			}
		}
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				if (uid_1(mt) == 0) {
					bias_1[i][j] = bias_0[i][j] + (float)uid_0(mt) / (float)digit - (float)width;
				}
				else {
					if (uid_1(mt) == 0) {
						bias_1[i][j] = (float)uid_0(mt) / (float)digit - (float)width;
					}
				}
			}
		}

	}

	_pclose(gp);
	return 0;
}


//シグモイド関数（活性化関数）
inline float func(float x)
{
	return 1 / (1 + exp(-x));
}
