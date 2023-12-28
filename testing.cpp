#include <iostream>
#include <vector>
#include "..\FFT\FFT.h"

template <class T>
class Polynomial {
private:
	int deg;
	std::vector <T> poly;
	static constexpr int _borderFFT = 32;
public:
	int _deg() const { return deg; }
	std::vector <T> _poly() const { return poly; }

	template <class U> friend std::istream & operator >> (std::istream &, Polynomial<U> &);
	template <class U> friend std::ostream & operator << (std::ostream &, const Polynomial<U> &);

	Polynomial() : deg(0) {
		poly = std::vector <T> (1, T(0));
	}
	Polynomial(int deg) {
		if (deg < 0) {
			throw std::runtime_error("Degree must be non-negative");
		}
		this->deg = deg;
		poly = std::vector <T> (deg + 1, T(0));
	}
	Polynomial(std::vector <T> ar) {
		if (ar.size() == 0) {
			throw std::runtime_error("Creating a polynomial of an empty vector");
		}
		deg = ar.size() - 1;
		poly = ar;
	}

	void changeDegree(int __deg) {
		poly.resize(__deg + 1);
		deg = __deg;
	}

	bool operator == (const Polynomial &);
	bool operator != (const Polynomial &);

	T operator()(const T &) const;

	template <class U> friend Polynomial<U> reverse(const Polynomial<U> &);
	template <class U> friend Polynomial<U> operator - (const Polynomial<U> &);

	Polynomial& operator += (const Polynomial &);
	Polynomial& operator -= (const Polynomial &);
	Polynomial& operator *= (const Polynomial &);

	Polynomial operator + (const Polynomial &) const;
	Polynomial operator - (const Polynomial &) const;
	Polynomial operator * (const Polynomial &) const;

	Polynomial operator % (const Polynomial &) const;
	Polynomial operator / (const Polynomial &) const;

	T operator[](const size_t &);
	const T operator[](const size_t &) const;

	template <class U> Polynomial<T>& operator *= (const U &);
	template <class U> Polynomial<T> operator * (const U &) const;

	Polynomial inverse_series(size_t n) const;
};

/* INPUT / OUTPUT */

template <class U>
std::istream & operator >> (std::istream& in, Polynomial<U>& p) {
	int D;
	in >> D;
	if (D < 0) {
		throw std::runtime_error("Degree must be non-negative");
	}
	p.deg = D;
	p.poly = std::vector <U> (D + 1);
	for (int i = 0; i <= D; ++i) {
		in >> p.poly[i];
	}
	return in;
}
template <class U>
std::ostream & operator << (std::ostream& out, const Polynomial<U>& p) {
	for (int i = 0; i <= p.deg; ++i) {
		out << (true ? p.poly[i] : abs(p.poly[i]));
		if (i) out << "x^{" << i << "}";

		if (i + 1 <= p.deg) out << (p.poly[i + 1] < 0 ? " - " : " + ");;
	}
	return out;
}

/* EQUALITY */

template <class T>
bool Polynomial<T>::operator == (const Polynomial<T>& rhs) {
	if (deg != rhs.deg) return false;
	for (size_t i = 0; i <= deg; ++i) {
		if (poly[i] != rhs.poly[i]) {
			return false;
		}
	}
	return true;
}
template <class T>
bool Polynomial<T>::operator != (const Polynomial<T>& rhs) {
	return !(*this == rhs);
}

/* NEGATION */
template <class U>
Polynomial<U> operator - (const Polynomial<U>& P) {
	Polynomial<U> result = P;
	for (int i = 0; i <= P.deg; ++i) {
		result.poly[i] = -result.poly[i];
	}
	return result;
}

template <class U>
Polynomial<U> reverse(const Polynomial<U>& P) {
	Polynomial<U> result = P;
	for (int i = 0; i <= P.deg - i; ++i) {
		std::swap( result.poly[i], result.poly[P.deg - i] );
	}
	return result;
}

template <class T>
T Polynomial<T>::operator()(const T& x) const {
	T result = T(0), current_power = T(1);
	for (int i = 0; i <= deg; ++i) {
		result += current_power * poly[i];
		current_power *= x;
	}
	return result;
}
template <class T>
T Polynomial<T>::operator[](const size_t& id) {
	if (id > deg) {
		throw std::runtime_error("Index out of bounds (id >= deg)");
	}
	return poly[id];
}
template <class T>
const T Polynomial<T>::operator[](const size_t& id) const {
	if (id > deg) {
		throw std::runtime_error("Index out of bounds (id >= deg)");
	}
	return poly[id];
}

template <class T>
template <class U> Polynomial<T>& Polynomial<T>::operator *= (const U& scalar) {
	for (int i = 0; i <= deg; ++i) poly[i] *= scalar;
	return *this;
}
template <class T>
template <class U> Polynomial<T> Polynomial<T>::operator * (const U& scalar) const {
	Polynomial<T> result = *this;
	result *= scalar;
	return result;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator += (const Polynomial<T>& rhs) {
	deg = std::max(deg, rhs.deg);
	poly.resize(deg + 1);

	for (size_t i = 0; i <= rhs.deg; ++i) poly[i] += rhs.poly[i];
	return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator -= (const Polynomial<T>& rhs) {
	deg = std::max(deg, rhs.deg);
	poly.resize(deg + 1);

	for (size_t i = 0; i <= rhs.deg; ++i) poly[i] -= rhs.poly[i];
	return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator *= (const Polynomial<T>& rhs) {
	if (deg <= _borderFFT && rhs.deg <= _borderFFT) {
		int n = deg;
		int m = rhs.deg;
		std::vector <T> result(n + m + 1);

		for (int i = 0; i <= n; ++i) {
			for (int j = 0; j <= m; ++j) {
				result[i + j] += poly[i] * rhs.poly[j];
			}
		}
		while (result.size() > 1 && result.back() == 0) result.pop_back();

		poly = result;
		deg = poly.size() - 1;

		return *this;
	}
	poly = FFT::multiply<T,T>(poly, rhs.poly);
	deg = poly.size() - 1;
	return *this;
}
template <class T>
Polynomial<T> Polynomial<T>::operator + (const Polynomial<T>& rhs) const {
	Polynomial<T> result = *this;
	result += rhs;
	return result;
}

template <class T>
Polynomial<T> Polynomial<T>::operator - (const Polynomial<T>& rhs) const {
	Polynomial<T> result = *this;
	result -= rhs;
	return result;
}

template <class T>
Polynomial<T> Polynomial<T>::operator * (const Polynomial<T>& rhs) const {
	Polynomial<T> result = *this;
	result *= rhs;
	return result;
}

template <class T>
Polynomial<T> Polynomial<T>::inverse_series(size_t n) const {
	const Polynomial<T> _two( std::vector <T> { T(2) } );

	int curDegree = 1;
	Polynomial<T> current( std::vector <T> { T(1) / poly[0] } );
	while (curDegree - 1 < n) {
		curDegree *= 2;
		current = current * ( _two - (*this) * current );
		current.changeDegree(curDegree);
	}
	current.changeDegree(n);
	return current;
}

template <class T>
Polynomial<T> Polynomial<T>::operator % (const Polynomial<T>& rhs) const {
	int N = deg, M = rhs.deg;
	if (N < M) {
		return *this;
	}
	Polynomial<T> result = reverse(*this) * (reverse(rhs)).inverse_series(N - M);
	result.changeDegree(N - M);
	return reverse(result);
}
template <class T>
Polynomial<T> Polynomial<T>::operator / (const Polynomial<T>& rhs) const {
	return *this - rhs * (*this % rhs);
}

#include <bits/stdc++.h>
using namespace std;

signed main() {
	FFT::init();
	Polynomial<double> a, b;
	cin >> a;
	cin >> b;

	cout << a / b << '\n' << a % b;
}
