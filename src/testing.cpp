#include <iostream>
#include <functional>
#include <optional>
#include <vector>
#include <cmath>

#include "..\FFT\FFT.h"
#include "..\FFT\FFT_With_Inverse.h"
#include "..\Modular\Modular.h"
#include "Polynomial.h"

using namespace std;

signed main() {
if (true) {
	Polynomial<long double> a;
	cin >> a;

	int n;
	cin >> n;
	cout << a.inverse_series(n) << '\n';
	cout << a * a.inverse_series(n) << '\n';
	cout << a.log(n) << '\n';
	cout << (a.log(n)).exp(n) << '\n';

	cout << a.power(2, n) << '\n';
	cout << "SQRT: " << (a.power(2, n)).sqrt(n).value() << '\n';

	int u;
	cin >> u;
	cout << (a << u) << '\n' << (a >> u) << '\n';

	int k, M;
	cin >> k >> M;
	cout << a.power(k, M) << '\n';

	int m;
	cin >> m;
	vector <long double> Points(m);
	for (auto &i : Points) cin >> i;

	auto c = a.evaluate(Points);
	for (int i = 0; i < m; ++i) {
		cout << "P(" << Points[i] << ") = " << c[i] << '\n';
	}
}
if (false) {
	const int mod = 998244353;
	Polynomial<ModInteger<mod>> a;
	cin >> a;

	cout << a.derivative() << '\n';
	cout << a.integral() << '\n';

	int n;
	cin >> n;
	cout << a.inverse_series(n) << '\n';
	cout << a * a.inverse_series(n) << '\n';
	cout << a.log(n) << '\n';
	cout << (a.log(n)).exp(n) << '\n';

	cout << a.power(2, n) << '\n';
	cout << "SQRT: " << (a.power(2, n)).sqrt(n).value() << '\n';

	int u;
	cin >> u;
	cout << (a << u) << '\n' << (a >> u) << '\n';

	int k, M;
	cin >> k >> M;
	cout << a.power(k, M) << '\n';

	int m;
	cin >> m;
	vector <ModInteger<mod>> Points(m);
	for (auto &i : Points) cin >> i;

	auto c = a.evaluate(Points);
	for (int i = 0; i < m; ++i) {
		cout << "P(" << Points[i] << ") = " << c[i] << '\n';
	}
}
}
