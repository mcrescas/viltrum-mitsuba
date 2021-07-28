/*
 * sample.h
 *
 *  Created on: Sep 2, 2019
 *      Author: miguel
 */

#ifndef SRC_INTEGRATORS_MDAS_SAMPLE_H_
#define SRC_INTEGRATORS_MDAS_SAMPLE_H_

/*
 * Sample struct - so bad , we need i guess a sampler specific for this
 */
struct Sample {
	float imageX;
	float imageY;
	float time;
	float lensU;
	float lensV;
	//float twoD[1][2];
	std::array<std::array<float, 2>, 1> twoD;

	Sample(float x, float y) : imageX(x) , imageY(y) {}
	Sample() {}
};


inline int Ceil2Int(double val) {
#ifdef FAST_INT
	return Round2Int(val + _doublemagicroundeps);
#else
	return (int)ceil(val);
#endif
}

inline int Floor2Int(double val) {
#ifdef FAST_INT
	return Round2Int(val - _doublemagicroundeps);
#else
	return (int)floor(val);
#endif
}

inline int Clamp(int val, int low, int high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}

inline float Clamp(float val, float low, float high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}


#endif /* SRC_INTEGRATORS_MDAS_SAMPLE_H_ */
