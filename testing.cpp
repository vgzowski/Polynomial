#include <bits/stdc++.h>
#include "Polynomial.h"

using namespace std;

signed main() {
{
	Polynomial<long double> a;
	cin >> a;

	int n;
	cin >> n;
	cout << a.inverse_series(n) << '\n';
	cout << a * a.inverse_series(n) << '\n';
	cout << a.log(n) << '\n';
	cout << (a.log(n)).exp(n) << '\n';

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
{
	Polynomial<ModInteger<13>> a;
	cin >> a;

	int n;
	cin >> n;
	cout << a.inverse_series(n) << '\n';
	cout << a * a.inverse_series(n) << '\n';
	cout << a.log(n) << '\n';
	cout << (a.log(n)).exp(n) << '\n';

	int u;
	cin >> u;
	cout << (a << u) << '\n' << (a >> u) << '\n';

	int k, M;
	cin >> k >> M;
	cout << a.power(k, M) << '\n';

	int m;
	cin >> m;
	vector <ModInteger<13>> Points(m);
	for (auto &i : Points) cin >> i;

	auto c = a.evaluate(Points);
	for (int i = 0; i < m; ++i) {
		cout << "P(" << Points[i] << ") = " << c[i] << '\n';
	}
}
}
