#pragma once

#include <iostream>
#include <optional>
#include <functional>
#include <vector>
#include <cmath>

#include "..\FFT\FFT.h"
#include "..\FFT\FFT_With_Inverse.h"
#include "..\Modular\Modular.h"

template < int _value > struct Int2Type { enum { value = _value }; };

template <class T>
class Polynomial {
private:
	int deg;
	std::vector <T> poly;
	static constexpr int _borderFFT = 32;

	inline void multiply( const Polynomial<T>& rhs, Int2Type <false> ) {
		if (deg <= _borderFFT && rhs.deg <= _borderFFT) {
			int n = deg;
			int m = rhs.deg;
			std::vector <T> result(n + m + 1);

			for (int i = 0; i <= n; ++i) {
				for (int j = 0; j <= m; ++j) {
					result[i + j] += poly[i] * rhs.poly[j];
				}
			}
			poly = result;
			removeZeros();
			deg = poly.size() - 1;
			return;
		}
		poly = FFT::multiply(poly, rhs.poly);
		removeZeros();
		deg = poly.size() - 1;
	}
	inline void multiply( const Polynomial<T>& rhs, Int2Type <true> ) {
		if (deg <= _borderFFT && rhs.deg <= _borderFFT) {
			int n = deg;
			int m = rhs.deg;
			std::vector <T> result(n + m + 1);

			for (int i = 0; i <= n; ++i) {
				for (int j = 0; j <= m; ++j) {
					result[i + j] += poly[i] * rhs.poly[j];
				}
			}
			poly = result;
			removeZeros();
			deg = poly.size() - 1;
			return;
		}

		auto Result = FFT_Inverse::multiply_mod(std::vector <int> (begin(poly), end(poly)),
						std::vector <int> (begin(rhs.poly), end(rhs.poly)),
						T::__mod_value);
		poly = std::vector <T> (begin(Result), end(Result));
		removeZeros();
		deg = poly.size() - 1;
	}

	inline std::optional < Polynomial > __sqrt( size_t n, Int2Type <true> ) const {
		size_t first_not_zero = 0;
		while (poly[first_not_zero] == T(0)) ++first_not_zero;

		if (first_not_zero & 1) {
			return std::nullopt;
		}

		Polynomial<T> current = *this >> first_not_zero;
		T coefficient = current[0];

		std::optional<T> sqrt_value = ::sqrt(coefficient);

		if (!sqrt_value.has_value()) {
			return std::nullopt;
		}

		current /= coefficient;

		return ( ((current.log(n) / 2).exp(n) << (std::min(first_not_zero / 2, n)))).upToN(n) * sqrt_value.value();
	}
	inline std::optional < Polynomial > __sqrt( size_t n, Int2Type <false> ) const {
		size_t first_not_zero = 0;
		while (poly[first_not_zero] == T(0)) ++first_not_zero;

		if (first_not_zero & 1) {
			return std::nullopt;
		}

		Polynomial<T> current = *this >> first_not_zero;
		T coefficient = current[0];

		if (coefficient < T(0)) {
			return std::nullopt;
		}

		T sqrt_value = ::sqrt(coefficient);

		current /= coefficient;

		return ( ((current.log(n) / 2).exp(n) << (std::min(first_not_zero / 2, n)))).upToN(n) * sqrt_value;
	}
public:
	int _deg() const { return deg; }
	std::vector <T> _poly() const { return poly; }

	Polynomial();
	Polynomial(int);
	Polynomial(std::vector <T>);

	template <class U> friend std::istream & operator >> (std::istream &, Polynomial<U> &);
	template <class U> friend std::ostream & operator << (std::ostream &, const Polynomial<U> &);

	void removeZeros();
	void changeDegree(int);
	void cutDegree(int);
	Polynomial upToN(int) const;

	template <class U> friend Polynomial<U> reverse(const Polynomial<U> &);
	template <class U> friend Polynomial<U> operator - (const Polynomial<U> &);

	T operator()(const T &) const;

	T& operator[](const size_t &);
	const T operator[](const size_t &) const;

	bool operator == (const Polynomial &) const;
	bool operator != (const Polynomial &) const;

	template <class U> Polynomial<T>& operator *= (const U &);
	template <class U> Polynomial<T>& operator /= (const U &);

	template <class U> Polynomial<T> operator * (const U &) const;
	template <class U> Polynomial<T> operator / (const U &) const;

	Polynomial& operator += (const Polynomial &);
	Polynomial& operator -= (const Polynomial &);
	Polynomial& operator *= (const Polynomial &);
	Polynomial& operator <<= (const size_t &);
	Polynomial& operator >>= (const size_t &);

	Polynomial operator + (const Polynomial &) const;
	Polynomial operator - (const Polynomial &) const;
	Polynomial operator * (const Polynomial &) const;
	Polynomial operator << (const size_t &) const;
	Polynomial operator >> (const size_t &) const;

	Polynomial operator % (const Polynomial &) const;
	Polynomial operator / (const Polynomial &) const;

	Polynomial integral() const;
	Polynomial derivative() const;

	Polynomial inverse_series(size_t) const;
	Polynomial log(size_t) const;
	Polynomial exp(size_t) const;

	template <typename U> Polynomial<T> power(U, size_t) const;

	std::optional < Polynomial <T> > sqrt(size_t) const;

	std::vector <T> evaluate( const std::vector < T > & ) const;
	template <class U> friend Polynomial<U> interpolate( const std::vector < std::pair <U, U> > & );
	template <class U> friend Polynomial<U> interpolate( const std::vector <U> &, const std::vector <U> & );
};

/*
CONSTRUCTOPRS
*/

template <class T>
Polynomial<T>::Polynomial() {
	this->deg = 0;
	this->poly = std::vector <T> (1, T(0));
}
template <class T>
Polynomial<T>::Polynomial(int _deg) {
	if (_deg < 0) {
		throw std::runtime_error("Negative degree");
	}
	this->deg = _deg;
	this->poly = std::vector <T> (_deg + 1, T(0));
}
template <class T>
Polynomial<T>::Polynomial(std::vector <T> _poly) {
	if (_poly.size() == 0) {
		throw std::runtime_error("Creating polynomial out of an empty vector");
	}
	this->deg = _poly.size() - 1;
	this->poly = _poly;
}

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
		out << (i == 0 ? p.poly[i] : abs(p.poly[i]));
		if (i) out << "x^{" << i << "}";

		if (i + 1 <= p.deg) out << (p.poly[i + 1] < 0 ? " - " : " + ");;
	}
	return out;
}

/*
MANIPULATORS
*/

template <class T>
void Polynomial<T>::removeZeros() {
	while (poly.size() > 1 && poly.back() == 0) poly.pop_back();
	deg = poly.size() - 1;
}

template <class T>
void Polynomial<T>::changeDegree(int __deg) {
	deg = __deg;
	poly.resize(__deg + 1);
}

template <class T>
void Polynomial<T>::cutDegree(int __deg) {
	poly.resize( std::min(__deg + 1, deg + 1) );
	deg = poly.size() - 1;
}

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
	for (int i = 0; i < P.deg - i; ++i) {
		std::swap( result.poly[i], result.poly[P.deg - i] );
	}
	return result;
}

template <class T>
Polynomial<T> Polynomial<T>::upToN(int n) const {
	Polynomial <T> result = *this;
	result.changeDegree(n);
	return result;
}

/*
ACCESSING AND EVALUATING
*/

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
T& Polynomial<T>::operator[](const size_t& id) {
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

/*
OPERATORS
*/

template <class T>
bool Polynomial<T>::operator == (const Polynomial<T>& rhs) const {
	if (deg != rhs.deg) return false;
	for (size_t i = 0; i <= deg; ++i) {
		if (poly[i] != rhs.poly[i]) {
			return false;
		}
	}
	return true;
}
template <class T>
bool Polynomial<T>::operator != (const Polynomial<T>& rhs) const {
	return !(*this == rhs);
}

template <class T>
template <class U> Polynomial<T>& Polynomial<T>::operator *= (const U& scalar) {
	for (int i = 0; i <= deg; ++i) poly[i] *= scalar;
	return *this;
}

template <class T>
template <class U> Polynomial<T> Polynomial<T>::operator * (const U& scalar) const {
	Polynomial<T> result = *this;
	return result *= scalar;
}

template <class T>
template <class U> Polynomial<T>& Polynomial<T>::operator /= (const U& scalar) {
	for (int i = 0; i <= deg; ++i) poly[i] /= scalar;
	return *this;
}

template <class T>
template <class U> Polynomial<T> Polynomial<T>::operator / (const U& scalar) const {
	Polynomial<T> result = *this;
	return result /= scalar;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator += (const Polynomial<T>& rhs) {
	deg = std::max(deg, rhs.deg);
	poly.resize(deg + 1);

	for (size_t i = 0; i <= rhs.deg; ++i) poly[i] += rhs.poly[i];
	removeZeros();
	return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator -= (const Polynomial<T>& rhs) {
	deg = std::max(deg, rhs.deg);
	poly.resize(deg + 1);

	for (size_t i = 0; i <= rhs.deg; ++i) poly[i] -= rhs.poly[i];
	removeZeros();
	return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator *= (const Polynomial<T>& rhs) {
	this->multiply(rhs, Int2Type<IsModular<T>::value>());
	return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator <<= (const size_t& _shift) {
	poly.insert(poly.begin(), _shift, T(0));
	deg = poly.size() - 1;
	return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator >>= (const size_t& _shift) {
	if (_shift > deg) {
		return (*this = Polynomial<T>());
	}
	else {
		poly.erase(poly.begin(), poly.begin() + _shift);
		deg = poly.size() - 1;
		return *this;
	}
}

template <class T>
Polynomial<T> Polynomial<T>::operator + (const Polynomial<T>& rhs) const {
	Polynomial<T> result = *this;
	return result += rhs;
}

template <class T>
Polynomial<T> Polynomial<T>::operator - (const Polynomial<T>& rhs) const {
	Polynomial<T> result = *this;
	return result -= rhs;
}

template <class T>
Polynomial<T> Polynomial<T>::operator * (const Polynomial<T>& rhs) const {
	Polynomial<T> result = *this;
	return result *= rhs;
}

template <class T>
Polynomial<T> Polynomial<T>::operator << (const size_t& k) const {
	Polynomial<T> result = *this;
	return result <<= k;
}

template <class T>
Polynomial<T> Polynomial<T>::operator >> (const size_t& k) const {
	Polynomial<T> result = *this;
	return result >>= k;
}

template <class T>
Polynomial<T> Polynomial<T>::operator / (const Polynomial<T>& rhs) const {
	int N = deg, M = rhs.deg;
	if (N < M) {
		return Polynomial<T>();
	}
	Polynomial<T> result = reverse(*this) * (reverse(rhs)).inverse_series(N - M);
	result.changeDegree(N - M);
	result = reverse(result);
	result.removeZeros();

	return result;
}

template <class T>
Polynomial<T> Polynomial<T>::operator % (const Polynomial<T>& rhs) const {
	return *this - rhs * (*this / rhs);
}

/*
FUNCTIONS ON POLYNOMIALS
*/

template <class T>
Polynomial<T> Polynomial<T>::inverse_series(size_t n) const {
	const Polynomial<T> _two( std::vector <T> { T(2) } );

	int curDegree = 1;

	Polynomial<T> current( std::vector <T> { T(1) / poly[0] } );

	while (curDegree - 1 < n) {
		curDegree *= 2;
		current = current * ( _two - (*this) * current );
		current.changeDegree(curDegree - 1);
	}
	current.changeDegree(n);
	return current;
}

template <class T>
Polynomial<T> Polynomial<T>::derivative() const {
	if (deg == 0) return Polynomial<T>();

	Polynomial<T> result(deg - 1);
	for (size_t i = 0; i + 1 <= deg; ++i) {
		result.poly[i] = poly[i + 1] * T(i + 1);
	}
	return result;
}

template <class T>
Polynomial<T> Polynomial<T>::integral() const {
	Polynomial<T> result(deg + 1);
	for (size_t i = 0; i <= deg; ++i) {
		result.poly[i + 1] = poly[i] / T(i + 1);
	}
	return result;
}

template <class T>
Polynomial<T> Polynomial<T>::log(size_t n) const {
	return (derivative().upToN(n) * inverse_series(n)).integral().upToN(n);
}

template <class T>
Polynomial<T> Polynomial<T>::exp(size_t n) const {
	const Polynomial<T> _one( std::vector <T> { T(1) } );

	int curDegree = 1;
	Polynomial<T> current( std::vector <T> { T(1) } );

	while (curDegree < n + 1) {
		curDegree *= 2;
		current = current * ( _one + upToN(curDegree - 1) - current.log(curDegree - 1) );
		current.changeDegree(curDegree - 1);
	}
	current.changeDegree(n);
	return current;
}

template <class T>
template <typename U>
Polynomial<T> Polynomial<T>::power(U k, size_t n) const {
	size_t first_not_zero = 0;
	while (poly[first_not_zero] == T(0)) ++first_not_zero;

	Polynomial<T> current = *this >> first_not_zero;
	
	T coefficient = current[0];

	current /= coefficient;
	return ( ((current.log(n) * k).exp(n) << (std::min(k * first_not_zero, n))) * pow(coefficient, k)).upToN(n);
}

template <class T>
std::optional < Polynomial <T> > Polynomial<T>::sqrt(size_t n) const {
	return __sqrt( n, Int2Type < IsModular<T>::value >() );
}

/*
FUNCTIONS ON POLYNOMIALS
*/

template <class T>
std::vector <T> Polynomial<T>::evaluate(const std::vector <T>& Pts) const {
	const int N = Pts.size();
	std::vector < Polynomial <T> > seg_tree(4 * N);

	std::function < void(int, int, int) > buildSegmentTree = [&]( int _v, int _vl, int _vr ) {
		if (_vl == _vr) {
			seg_tree[_v] = Polynomial <T> ( std::vector <T> { -Pts[_vl], T(1) } );
			return;
		}
		int _m = (_vl + _vr) >> 1;
		buildSegmentTree( _v << 1, _vl, _m );
		buildSegmentTree( _v << 1 | 1, _m + 1, _vr );
		seg_tree[_v] = seg_tree[_v << 1] * seg_tree[_v << 1 | 1];
	};

	buildSegmentTree(1, 0, N - 1);

	std::vector <T> results(N);

	std::function < void(int, int, int, const Polynomial<T> &) > traverseSegmentTree = [&]( int _v, int _vl, int _vr, const Polynomial <T>& P ) {
		if (_vl == _vr) {
			results[_vl] = P[0];
			return;
		}
		int _m = (_vl + _vr) >> 1;
		traverseSegmentTree( _v << 1, _vl, _m, P % seg_tree[_v << 1] );
		traverseSegmentTree( _v << 1 | 1, _m + 1, _vr, P % seg_tree[_v << 1 | 1] );
	};
	traverseSegmentTree(1, 0, N - 1, *this % seg_tree[1]);
	return results;
}

template <class U>
Polynomial<U> interpolate(const std::vector < std::pair <U, U> >& Pts) {
	const int N = Pts.size();
	std::vector <U> X(N), Y(N);
	for (size_t i = 0; i < N; ++i) {
		X[i] = Pts[i].first;
		Y[i] = Pts[i].first;
	}
	return interpolate(X, Y);
}

template <class U>
Polynomial<U> interpolate(const std::vector <U>& X, const std::vector <U>& Y) {
	const int N = X.size();
	std::vector < Polynomial <U> > seg_tree(4 * N);

	std::function < void(int, int, int) > buildSegmentTreeProduct = [&]( int _v, int _vl, int _vr ) {
		if (_vl == _vr) {
			seg_tree[_v] = Polynomial <U> ( std::vector <U> { -X[_vl], U(1) } );
			return;
		}
		int _m = (_vl + _vr) >> 1;
		buildSegmentTreeProduct( _v << 1, _vl, _m );
		buildSegmentTreeProduct( _v << 1 | 1, _m + 1, _vr );
		seg_tree[_v] = seg_tree[_v << 1] * seg_tree[_v << 1 | 1];
	};

	buildSegmentTreeProduct(1, 0, N - 1);

	Polynomial <U> P = seg_tree[1];
	Polynomial <U> dP = P.derivative();

	std::vector <U> Results = dP.evaluate( X );

	std::vector < Polynomial <U> > seg_tree_inter(4 * N);

	std::function < void(int, int, int) > buildSegmentTreeInterpolation = [&]( int _v, int _vl, int _vr ) {
		if (_vl == _vr) {
			seg_tree_inter[_v] = Polynomial <U> ( std::vector <U> { Y[_vl] / Results[_vl] } );
			return;
		}
		int _m = (_vl + _vr) >> 1;
		buildSegmentTreeInterpolation( _v << 1, _vl, _m );
		buildSegmentTreeInterpolation( _v << 1 | 1, _m + 1, _vr );
		seg_tree_inter[_v] = seg_tree_inter[_v << 1] * seg_tree[_v << 1 | 1] + seg_tree_inter[_v << 1 | 1] * seg_tree[_v << 1];
	};

	buildSegmentTreeInterpolation(1, 0, N - 1);
	return seg_tree_inter[1];
}
