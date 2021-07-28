/*
 * TupleMath.h
 *
 *  Created on: May 25, 2018
 *      Author: miguel
 */

#ifndef SRC_INTEGRATORS_QUADRATURE_TUPLEMATH_H_
#define SRC_INTEGRATORS_QUADRATURE_TUPLEMATH_H_

#include <mitsuba/core/spectrum.h>

// Type that stores info
//typedef double TUPLE_T;


//#define VECTORIZED_TUPLEMATH

#ifdef VECTORIZED_TUPLEMATH
	#if !__has_include(<x86intrin.h>)
		#error No se puede importar
	#endif
	//#include <emmintrin.h>
	#include <x86intrin.h>
	//#include <immintrin.h>

	//typedef double TUPLE_T __attribute__ ((mode(V4DF)));
	typedef __m256d  TUPLE_SYSTEM;
	typedef double TUPLE_T;

	#ifndef __AVX2__  //
		#error AVX2 not supported
	#endif
#else
	typedef double TUPLE_T;
	typedef double TUPLE_SYSTEM;
#endif

/*
 * 	Aux class for use to avoid numeric precisision errors
 */

class TupleMath {
public:

	static const std::size_t TUPLE_NVALUES = 3;
	/*
	 * Constructors
	 */

	/*friend TupleMath abs (const TupleMath& tuple) {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res[i] = std::abs(tuple.values[i]);
		}
		return res;
	}*/

	/*
	 * Add custom support to add arithmetic data types 
	 */
	/*template <typename Type,
             std::enable_if_t<std::is_arithmetic<
                std::remove_reference_t<Type>>::value> * = nullptr>
	TupleMath operator+ (const Type&& rhs) {
		TupleMath res;
		for(int d=0; d<TUPLE_NVALUES; d++) {
			res = values[d] + rhs;
		}
		return res;
	}*/

	class Iterator {
	public:
		Iterator(const TUPLE_T * _ptr): ptr(_ptr) {}
		Iterator operator++() { ++ptr; return *this; }
		bool operator!=(const Iterator & other) const { return ptr != other.ptr; }
		TUPLE_T operator*() { return *ptr; }
		const TUPLE_T operator*() const { return *ptr; }
	private:
		const TUPLE_T* ptr;
	};

	Iterator begin() const {
		//return Iterator(values);
		return Iterator(&values[0]);
	}

	Iterator end() const {
		return Iterator(&values[TUPLE_NVALUES-1]);
	}

	TupleMath operator+(double x){
		TupleMath res;
		for(int d=0; d<TUPLE_NVALUES; d++) {
			res = values[d] + x;
		}
		return res;
	}

	TupleMath operator+(double x) const {
		TupleMath res;
		for(int d=0; d<TUPLE_NVALUES; d++) {
			res = values[d] + x;
		}
		return res;
	}

#ifndef VECTORIZED_TUPLEMATH
	/*TupleMath(mitsuba::Spectrum spec) {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] = spec[i];
		}
	}*/

	TupleMath(mitsuba::Spectrum& spec) {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] = spec[i];
		}
	}

	TupleMath(TUPLE_T val) {
		values[0] = values[1] = values[2] = val;
	}

	TupleMath() {
		values[0] = 0;
		values[1] = 0;
		values[2] = 0;
	}

	/*TupleMath(TupleMath& other) {
		values[0] = other.values[0];
		values[1] = other.values[1];
		values[2] = other.values[2];
	}*/
#else
	TupleMath(TUPLE_T val) {
		//values[0] = values[1] = values[2] = val;
		values = _mm256_set1_pd(val);
	}

	TupleMath(mitsuba::Spectrum spec) {
		/*for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] = spec[i];
		}*/
		_mm256_set_pd(0, spec[2], spec[1], spec[0]);
	}

	TupleMath() {
		values = _mm256_set1_pd(-1);
	}
#endif

#ifndef VECTORIZED_TUPLEMATH
	/*
	 * Copy overload
	 */
	TupleMath & operator=(const mitsuba::Spectrum & spec)
	{
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] = spec[i];
		}
		return *this;
	}

	TupleMath & operator=(mitsuba::Spectrum & spec)
	{
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] = spec[i];
		}
		return *this;
	}

	TupleMath & operator=(int val)
	{
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] = val;
		}
		return *this;
	}

#else
	TupleMath & operator=(const mitsuba::Spectrum & spec)
	{
		_mm256_set_pd(0, spec[2], spec[1], spec[0]);
		return *this;
	}

	TupleMath & operator=(mitsuba::Spectrum & spec)
	{
		_mm256_set_pd(0, spec[2], spec[1], spec[0]);
		return *this;
	}

	TupleMath & operator=(int val)
	{
		values = _mm256_set1_pd(val);
		return *this;
	}

#endif

	/*
	 * Cast overloading
	 */
	operator mitsuba::Spectrum() {
		mitsuba::Spectrum res;
		res[0] = values[0];
		res[1] = values[1];
		res[2] = values[2];
		return res;
	}

	/*
	 * Return a spectrum from this data
	 */
	mitsuba::Spectrum spectrum() {
		mitsuba::Spectrum res;
		res[0] = values[0];
		res[1] = values[1];
		res[2] = values[2];
		return res;
	}

	double getLuminance() {
		return values[0] * 0.212671f + values[1] * 0.715160f + values[2] * 0.072169f;
	}

	void nonNegative() {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] = (values[i] >= 0.0) ? values[i] : 0.0;
		}
	}

	bool isNaN() {
		for (int i=0; i<TUPLE_NVALUES; i++) {
            if (std::isnan(values[i])) { 
                return true;
			}
		}
        return false;
	}

	/*
	 * Brackets overloading
	 */
#ifdef VECTORIZED_TUPLEMATH
	double& operator[](int index) {
		return ((double *)&values)[index];
	}

	const double& operator[](int index) const {
		return ((double *)&values)[index];
	}
#else
	TUPLE_T& operator[](int index) {
		return values[index];
	}

	const TUPLE_T& operator[](int index) const {
		return values[index];
	}
#endif

#ifdef VECTORIZED_TUPLEMATH
	double max() const {
		TUPLE_T res = values[0];
		for(int i=1; i<TUPLE_NVALUES; i++) {
			if (values[i] > res) res = values[i];
		}
		return res;
	}

	double min() const {
		TUPLE_T res = values[0];
		for(int i=1; i<TUPLE_NVALUES; i++) {
			if (values[i] < res) res = values[i];
		}
		return res;
	}
#else
	TUPLE_T max() const {
		TUPLE_T res = values[0];
		for(int i=1; i<TUPLE_NVALUES; i++) {
			if (values[i] > res) res = values[i];
		}
		return res;
	}

	TUPLE_T min() const {
		TUPLE_T res = values[0];
		for(int i=1; i<TUPLE_NVALUES; i++) {
			if (values[i] < res) res = values[i];
		}
		return res;
	}
#endif

	TupleMath sqrt() const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res[i] = std::sqrt(values[i]);
		}
		return res;
	}

	TupleMath pow(int a) const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res[i] = std::pow(values[i], a);
		}
		return res;
	}

	TupleMath abs() const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			//res[i] = fabs(values[i]);
			res[i] = std::abs(values[i]);
		}
		return res;
	}

	// Friend function to overload abs in generic context (use ADL trick using std::abs to resolve to this function or the generic one depending of the type)
	friend TupleMath abs(const TupleMath& t) {
		return t.abs();
	}

	TUPLE_T sum() const {
		TUPLE_T res = 0.0;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res += values[i];
		}
		return res;
	}

	std::size_t size() const {
		return TUPLE_NVALUES;
	}


	/*
	 * Math operators
	 */
#ifdef VECTORIZED_TUPLEMATH
	inline TupleMath operator+(const TupleMath& other) {
		TupleMath res;
		res.values = _mm256_add_pd(this->values, other.values);
		return res;
	}

	inline TupleMath operator+=(const TupleMath& other) {
		values = _mm256_add_pd(this->values, other.values);
		return *this;
	}

	inline TupleMath operator-(const TupleMath& other) {
		TupleMath res;
		res.values = _mm256_sub_pd(this->values, other.values);
		return res;
	}

	inline TupleMath operator-=(const TupleMath& other) {
		values = _mm256_sub_pd(this->values, other.values);
		return *this;
	}

	inline TupleMath operator*(const TupleMath& other) {
		TupleMath res;
		res.values = _mm256_mul_pd(this->values, other.values);
		return res;
	}

	inline TupleMath operator*=(const TupleMath& other) {
		values = _mm256_mul_pd(this->values, other.values);
		return *this;
	}

	inline TupleMath operator*(double other) {
		const TUPLE_SYSTEM scalar = _mm256_set1_pd(other);
		TupleMath res;
		res.values = _mm256_mul_pd(this->values, scalar);
		return res;
	}

	inline TupleMath operator*=(double other) {
		const TUPLE_SYSTEM scalar = _mm256_set1_pd(other);
		values = _mm256_mul_pd(this->values, scalar);
		return *this;
	}

	inline TupleMath operator/(const TupleMath& other) {
		TupleMath res;
		res.values = _mm256_div_pd(this->values, other.values);
		return res;
	}

	inline TupleMath operator/=(const TupleMath& other) {
		values = _mm256_div_pd(this->values, other.values);
		return *this;
	}

	inline TupleMath operator/(double other) {
		TupleMath res;
		const TUPLE_SYSTEM scalar = _mm256_set1_pd(other);
		res.values = _mm256_div_pd(this->values, scalar);
		return res;
	}

	inline TupleMath operator/=(double other) {
		const TUPLE_SYSTEM scalar = _mm256_set1_pd(other);
		values = _mm256_div_pd(this->values, scalar);
		return *this;
	}
#else
	inline TupleMath operator+(const TupleMath& other) const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res.values[i] = this->values[i] + other.values[i];
		}

		return res;
	}

	inline TupleMath operator+=(const TupleMath& other) {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] += other.values[i];
		}

		return *this;
	}

	inline TupleMath operator-(const TupleMath& other) const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res.values[i] = this->values[i] - other.values[i];
		}

		return res;
	}

	inline TupleMath operator-=(const TupleMath& other) {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] -= other.values[i];
		}

		return *this;
	}

	inline TupleMath operator*(const TupleMath& other) const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res.values[i] = this->values[i] * other.values[i];
		}

		return res;
	}

	inline TupleMath operator*=(const TupleMath& other) {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] *= other.values[i];
		}

		return *this;
	}

	inline TupleMath operator*(double other) const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res.values[i] = this->values[i] * other;
		}

		return res;
	}

	inline TupleMath operator*=(double other) {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] *= other;
		}

		return *this;
	}

	inline TupleMath operator/(const TupleMath& other) const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res.values[i] = this->values[i] / other.values[i];
		}

		return res;
	}

	inline TupleMath operator/=(const TupleMath& other) {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] /= other.values[i];
		}

		return *this;
	}

	inline TupleMath operator/(double other) const {
		TupleMath res;
		for(int i=0; i<TUPLE_NVALUES; i++) {
			res.values[i] = this->values[i] / other;
		}

		return res;
	}

	inline TupleMath operator/=(double other) {
		for(int i=0; i<TUPLE_NVALUES; i++) {
			values[i] /= other;
		}

		return *this;
	}

	inline bool operator<(const TupleMath& other) {
		return values[0] < other.values[0] && values[1] < other.values[1] && values[2] < other.values[2];
	}

#endif

	inline bool isValid() {
		for (int i=0; i<TUPLE_NVALUES; i++)
			if (!std::isfinite(values[i]) || values[i] < 0.0f)
				return false;
		return true;
	}

	inline bool isZero() {
		for (int i=0; i<TUPLE_NVALUES; i++)
			if (values[i] != 0.0)
				return false;
		return true;
	}

	inline bool anyZero() {
		for (int i=0; i<TUPLE_NVALUES; i++)
			if (values[i] == 0.0)
				return true;
		return false;
	}

	inline TupleMath setEpsilon(TUPLE_T _epsilon=Epsilon) {
		TupleMath res;
		for (int i=0; i<TUPLE_NVALUES; i++) {
			if (values[i] == 0.0) {
				res[i] = _epsilon;
			} else {
				res[i] = values[i];
			}
		}
		return res;
	}

	inline bool operator == (const TupleMath& other) {
		return values[0] == other.values[0] && values[1] == other.values[1] && values[2] == other.values[2];
	}

	inline bool operator != (const TupleMath& other) {
		return values[0] != other.values[0] || values[1] != other.values[1] || values[2] != other.values[2];
	}


#ifdef VECTORIZED_TUPLEMATH
	TUPLE_SYSTEM values;
#else
	//TUPLE_T values[TUPLE_NVALUES];
	TUPLE_T values[TUPLE_NVALUES];
#endif
};

#ifdef VECTORIZED_TUPLEMATH
template <typename A>
TupleMath operator* (A val, TupleMath tuple) {
	const TUPLE_SYSTEM scalar = _mm256_set1_pd(val);
	TupleMath res;
	res.values = _mm256_mul_pd(scalar, tuple.values);
	return res;
}

#else

template <typename A>
TupleMath operator* (A val, TupleMath tuple) {
	TupleMath res;

	res.values[0] = val * tuple.values[0];
	res.values[1] = val * tuple.values[1];
	res.values[2] = val * tuple.values[2];

	return res;
}

#endif



#endif /* SRC_INTEGRATORS_QUADRATURE_TUPLEMATH_H_ */
