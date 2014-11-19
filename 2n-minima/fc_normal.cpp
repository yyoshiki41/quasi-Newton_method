#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define N 2      // N次元
#define K 1000   // イテレーション
#define EPS 0.01 // 収束条件
#define c 0.01   // ステップ幅
#define M 1000   // 試行回数

double __2n_minima (double x1, double x2, int *call);

int main()
{
	char filepath[256];
	sprintf(filepath, "function_call/normal.txt");
	ofstream fout;
	fout.open(filepath);
	fout << "fc数\t関数値" << endl;

	for (int m = 1; m <= M; m++) {
		fout << "試行回数\t" << m << endl;
		int call = 0; // function call数のリセット

		double x[N], result;
		for (int n = 0; n < N; n++) {
			x[n] = ((double)rand() / ((double)RAND_MAX + 1)) * 10 - 5;
		}
		double y[N], f[N], z[N];
		double p[N], q[N];
		double dx1, dx2, norm;

		// 近似行列の初期値は、単位行列
		double j11 = 1, j12 = 0, j21 = 0, j22 = 1;
		double g11 = 0, g12 = 0, g21 = 0, g22 = 0;

		for (int k = 0; k < K; k++) {
			f[0] = 4 * pow(x[0], 3) - 32 * x[0] + 5;
			f[1] = 4 * pow(x[1], 3) - 32 * x[1] + 5;

			if (k != 0) {
				norm = sqrt(dx1 * dx1 + dx2 * dx2);
				if(norm < EPS) break;

				for (int n = 0; n < N; n++) {
					p[n] = x[n] - y[n];
					q[n] = f[n] - z[n];
				}
				g11 = p[0]*p[0] / (p[0]*q[0]+p[1]*q[1]) - ((j11*q[0]+j21*q[1]) * (j11*q[0]+j12*q[1])) / (q[0]*(q[0]*j11+q[1]*j21) + q[1]*(q[0]*j12+q[1]*j22));
				g12 = p[0]*p[1] / (p[0]*q[0]+p[1]*q[1]) - ((j12*q[0]+j22*q[1]) * (j11*q[0]+j12*q[1])) / (q[0]*(q[0]*j11+q[1]*j21) + q[1]*(q[0]*j12+q[1]*j22));
				g21 = p[1]*p[0] / (p[0]*q[0]+p[1]*q[1]) - ((j11*q[0]+j21*q[1]) * (j21*q[0]+j22*q[1])) / (q[0]*(q[0]*j11+q[1]*j21) + q[1]*(q[0]*j12+q[1]*j22));
				g22 = p[1]*p[1] / (p[0]*q[0]+p[1]*q[1]) - ((j12*q[0]+j22*q[1]) * (j21*q[0]+j22*q[1])) / (q[0]*(q[0]*j11+q[1]*j21) + q[1]*(q[0]*j12+q[1]*j22));
			}

			// Hesseの近似行列の各成分
			j11 = j11 + g11;
			j12 = j12 + g12;
			j21 = j21 + g21;
			j22 = j22 + g22;

			dx1 = j11 * f[0] + j12 * f[1];
			dx2 = j21 * f[0] + j22 * f[1];

			for (int n = 0; n < N; n++) {
				y[n] = x[n];
				z[n] = f[n];
			}

			x[0] = x[0] - c * dx1;
			x[1] = x[1] - c * dx2;
			result = __2n_minima (x[0], x[1], &call);

			if (call <= 200 && call % 5 == 0) {
				fout << call << "\t" << result << endl;
			}
		}
	}

	return 0;
}

// 目的関数
double __2n_minima (double x1, double x2, int *call)
{
	(*call)++;
	double f;
	f = pow(x1, 4) - 16 * pow(x1, 2) + 5 * x1 + pow(x2, 4) - 16 * pow(x2, 2) + 5 * x2;
	return f;
}
