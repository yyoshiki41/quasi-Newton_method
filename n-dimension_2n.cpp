#include <iostream>
#include <fstream>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

#define K 1000
#define EPS 0.01
#define c 0.01

int main()
{
	// n次元設定
	cout << "Dimention: ";
	int N;
	cin >> N;
	MatrixXd x(N, 1);
	for (int n = 0; n < N; n++) {
		cin >> x(n, 0);
	}

	// file 出力
	char filepath[256];
	sprintf(filepath, "normal_2n/%d-dimension.txt", N);
	ofstream fout;
	fout.open(filepath);
	// 見出し
	fout << "k";
	for (int n = 0; n < N; n++) {
		fout << "\tx" << n+1;
	}
	fout << endl;

	// 行列の初期化
	MatrixXd y(N, 1);
	MatrixXd f(N, 1);
	MatrixXd z(N, 1);
	MatrixXd p(N, 1);
	MatrixXd q(N, 1);
	MatrixXd dx(N, 1);
	MatrixXd J = MatrixXd::Identity(N, N); // Hesse行列の近似行列の初期値は、単位行列
	MatrixXd G = MatrixXd::Zero(N, N); // 更新行列の初期化

	for (int k = 0; k < K; k++) {
		fout << k;
		for (int n = 0; n < N; n++) {
			fout << "\t" << x(n, 0);
		}
		fout << endl;

		for (int n = 0; n < N; n++) {
			f(n, 0) = 4 * pow(x(n, 0), 3) - 32 * x(n, 0) + 5;
		}

		if (k != 0) {
			// 終了条件
			if (dx.norm() < EPS) break;
			// DFP 公式
			p = x - y;
			q = f - z;
			G = (p * p.transpose() / (p.transpose() * q)(0, 0)) - (J * q * q.transpose() * J / (q.transpose() * J * q)(0, 0));
		}

		y = x;
		z = f;

		// Hesse行列の近似行列
		J = J + G;
		dx = J * f;
		x = x - c * dx;
	}

	return 0;
}
