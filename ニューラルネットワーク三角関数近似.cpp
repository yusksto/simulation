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
	int width = 5;  //乱数範囲[-width , width]
	int digit = 1000;//乱数小数点精度（digit分の１まで）
	std::uniform_int_distribution<> uid_0(0, 2 * width * digit);
	std::uniform_int_distribution<> uid_1(0,500);

	//使用配列
	float weight_0[9][10][10], weight_1[9][10][10];
	float weight_00[10], weight_01[10], weight_10[10], weight_11[10];
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
	for (int i = 0; i < 10; i++) {
		weight_00[i] = (float)uid_0(mt) / (float)digit - (float)width;
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
	for (int i = 0; i < 10; i++) {
		weight_10[i] = weight_00[i];
		weight_11[i] = weight_00[i];
	}
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			bias_1[i][j] = bias_0[i][j];
		}
	}

	//開始
	const int n_0 = 1000;                 //最大変異回数
	int n = 0;                            //変異回数カウント
	const int R = 1000000;                 //リセット回数
	int r = 0;                            //リセット回数カウント
	int t = 1;                            //総リセット回数カウント
	const int N = 100;                    //誤差計算回数
	float δ = 99999999999999, δ_0 = 99999999999999, δ_1 = 0;//誤差
	const int F = 1;                      //画像出力頻度
	int f = 0;                            //画像出力カウント用
	float in = 0;
	float out = 0;
	fprintf(gp, "set terminal png\n");
	while (n < n_0) {
		r++;
		//変数初期化
		δ_1 = 0;
		out = 0;
		Fill(network, 0);

		//誤差計算
		for (int m = 0; m < N; m++) {
			in = 2 * PI * (float)m / (float)N;

			//ニューラルネットワーク計算
			for (int i = 0; i < 10; i++) {
				network[0][i] = func(in * weight_10[i] + bias_1[0][i]);
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
				out += network[9][i] * weight_11[i];
			}
			out = 4 * func(out) - 2;

			δ_1 += pow(abs(out - sin(in)) + 1, 2);
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
				std::cout << n << " " << δ_0 << " " << t <<std::endl;
				fprintf(gp, "set output '%d.png'\n", t * 1000000 + n);
				fprintf(gp, "plot '-' with lines linetype rgbcolor '#0000FF', '-' with lines linetype rgbcolor '#FF0000'\n");
				for (int i = 0; i < N; i++) {
					in = 2 * PI * (float)i / (float)N;

					//ニューラルネットワーク計算
					for (int i = 0; i < 10; i++) {
						network[0][i] = func(in * weight_10[i] + bias_1[0][i]);
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
						out += network[9][i] * weight_11[i];
					}
					out = 4 * func(out) - 2;

					fprintf(gp, "%lf, %lf\n", in, out);
				}
				fprintf(gp, "e\n");
				fflush(gp);
				for (int i = 0; i < N; i++) {
					in = 2 * PI * (float)i / (float)N;
					fprintf(gp, "%lf, %lf\n", in, sin(in));
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
			for (int i = 0; i < 10; i++) {
				weight_00[i] = weight_10[i];
				weight_01[i] = weight_10[i];
			}
			for (int i = 0; i < 10; i++) {
				for (int j = 0; j < 10; j++) {
					bias_0[i][j] = bias_1[i][j];
				}
			}
		}

		//リセット判定
		if (r > R){
			//リセット
			r = 0;
			n = 0;
			if (δ > δ_0){
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
			for (int i = 0; i < 10; i++) {
				weight_00[i] = (float)uid_0(mt) / (float)digit - (float)width;
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
					if (uid_1(mt) == 0){
						weight_1[i][j][k] = weight_0[i][j][k] + (float)uid_0(mt) / (float)digit - (float)width;
					} else {
						if (uid_1(mt) == 0) {
							weight_1[i][j][k] = (float)uid_0(mt) / (float)digit - (float)width;
						}
					}
				}
			}
		}
		for (int i = 0; i < 10; i++) {
			if (uid_1(mt) == 0){
				weight_10[i] = weight_00[i] + (float)uid_0(mt) / (float)digit - (float)width;
			} else {
				if (uid_1(mt) == 0) {
					weight_10[i] =(float)uid_0(mt) / (float)digit - (float)width;
				}
			}
			if (uid_1(mt) == 0){
				weight_11[i] = weight_00[i] + (float)uid_0(mt) / (float)digit - (float)width;
			} else {
				if (uid_1(mt) == 0) {
					weight_11[i] = (float)uid_0(mt) / (float)digit - (float)width;
				}
			}
		}
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				if (uid_1(mt) == 0){
					bias_1[i][j] = bias_0[i][j] + (float)uid_0(mt) / (float)digit - (float)width;
				} else {
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
