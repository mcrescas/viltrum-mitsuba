/*
 * tmp_sampler.h
 *
 *  Created on: Sep 3, 2019
 *      Author: miguel
 */

#ifndef SRC_INTEGRATORS_MDAS_TMP_SAMPLER_H_
#define SRC_INTEGRATORS_MDAS_TMP_SAMPLER_H_


#include <mitsuba/render/sampler.h>
#include <mitsuba/core/properties.h>

/*
 * Ad Hoc sampler class
 * should work without problems
 */
class TMP_SAMPLER : public mitsuba::Sampler {
public:
	TMP_SAMPLER(std::array<float, 2> data) : _data(data) , mitsuba::Sampler(mitsuba::Properties()) {}
	mitsuba::Point2 next2D() override {
		if (index < 2) {index += 2; return mitsuba::Point2(_data[0], _data[1]);}
		else return mitsuba::Point2(0.5, 0.5);
	}


	mitsuba::Float next1D() override {
		if (index < 2) { return _data[index++]; }
		else return 0.5;
	}
private:
	std::array<float, 2> _data;
	int index = 0;
};


#endif /* SRC_INTEGRATORS_MDAS_TMP_SAMPLER_H_ */
