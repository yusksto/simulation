#include <iostream>
#include <cstdio>
#include <cassert>
#include <conio.h>

int main()
{
    FILE* gp = _popen("gnuplot", "w");
    assert(gp);
    fprintf(gp, "unset key\n");

    const int N = 100000000;//人口
    const float k_0 = 1.3f;  //基本再生産数
    const int h = 14;        //回復所要日数


    int N_0 = 100; //感染者数
    int N_1 = 0;   //回復者数
    int t = 0;     //日数
    int buffa[h];  //回復者バッファ用
    for (int i = 0; i < h; i++){buffa[i] = 0;}//初期化


    //シミュレーション開始
    fprintf(gp, "plot '-' with lines\n");
    N_0 = 100; //感染者初期数
    N_1 = 1000;    //回復者初期数
    t = 0;      //日数初期数
    buffa[h - 1] = N_0;

    while (N_0 != 0)
    {
        t++;
        fprintf(gp, "%d, %d\n", t, N_0);

        N_0 -= buffa[0];
        N_1 += buffa[0];
        for (int s = 0; s < h - 1; s++)
        {
            buffa[s] = buffa[s + 1];
        }

        buffa[h - 1] = (int)((float)N_0 * (float)k_0 * (float)(N - N_0 - N_1) / (float)N / h);
        N_0 += buffa[h - 1];
    }
    fprintf(gp, "e\n");
    fflush(gp);
    
    float p_ = (float)N_1 / (float)N;
    std::cout << p_ << " " << t << std::endl;

    while (1) {
        if ('\r' == _getch()) break;
    }
}
