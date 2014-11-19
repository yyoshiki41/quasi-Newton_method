#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define N 2
#define K 1000
#define EPS 0.01
#define c 0.01

int main()
{
	double x[N];
	for (int n = 0; n < N; n++) {
		cin >> x[n];
	}
	double y[N], f[N], z[N];
	double p[N], q[N];
	double dx1, dx2, norm;

	// 近似行列の初期値は、単位行列
	double j11 = 1, j12 = 0, j21 = 0, j22 = 1;
	double g11 = 0, g12 = 0, g21 = 0, g22 = 0;

	char filepath[256];
	sprintf(filepath, "normal/%2.1f,%2.1f.txt", x[0], x[1]);
	ofstream fout; // file出力の為の定義
	fout.open(filepath); // fileを開く
	fout << "k\tx1\tx2" << endl;

	for (int k = 0; k < K; k++) {
		fout << k << "\t";
		fout << x[0] << "\t";
		fout << x[1] << endl;

		for (int n = 0; n < N; n++) {
			z[n] = f[n];
		}
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
		}

		x[0] = x[0] - c * dx1;
		x[1] = x[1] - c * dx2;
	}

	return 0;
}
