# 準ニュートン法
> quasi-Newton method

## ベンチマーク関数
- 2n-minima
<pre>
// 2変数
f(x, y) = (x^4 - 16*x^2 + 5*x) + (y^4 - 16*y^2 + 5*y)
// 多変数
f(x) = Σ (x^4 - 16*x^2 + 5*x)
</pre>
- rosenbrock
<pre>
// 2変数
f(x, y) = (1-x)^2 + 100*(y - x^2)^2
// 多変数
f(x_n) = Σ {(1-x_n)^2 + 100*(x_n+1 - x_n^2)^2}
</pre>

## 直線探索
黄金分割法


## ディレクトリ
<pre>
.
├── 2n-minima             2n-minima 関数
│   ├── normal            普通の準ニュートン法(ステップ幅 c=0.01)
│   ├── gold              黄金分割法
│   ├── n-dimensions      n次元への拡張
│   ├── run               試行回数 m回での収束回数
│   └── function_call　　　function call数と目的関数値をカウント
|
└── rosenbrock            rosenbrock 関数
    ├── normal            普通の準ニュートン法(ステップ幅 c=0.01)
    ├── gold              黄金分割法
    ├── n-dimensions      n次元への拡張
    ├── run               試行回数 m回での収束回数
    └── function_call　　　function call数と目的関数値をカウント
</pre>
