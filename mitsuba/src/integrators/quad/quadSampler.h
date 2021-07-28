
#if !defined(__QUADRATURE_SAMPLER_H)
#define __QUADRATURE_SAMPLER_H

#include <mitsuba/render/sampler.h>

//#include "constants.h"
#include <random>

#define GAP_ERROR 1e-7

//#define DIRECT

//#define ERROR_EXTRA_SAMPLE

MTS_NAMESPACE_BEGIN

//template <int C_DIM_N, int TOTAL_SAMPLES_VALID>
template <int TOTAL_SAMPLES_VALID, typename RNG>
class QuadratureSampler : public Sampler {

	using DATA_TYPE = double;
	using ARRAY_TYPE = std::array<DATA_TYPE, TOTAL_SAMPLES_VALID>;
	const DATA_TYPE CONSTANT_DATA = 0.5;

public:
	QuadratureSampler(const ARRAY_TYPE& _values, RNG& _rng, bool useCustomSampler = false)
		: Sampler(Properties()),  m_useCustomSampler(useCustomSampler), rng(_rng), RETURN_CONSTANT(CONSTANT_DATA),uniform_dist(0.0 , 1.0){

			values = _values;
			index = 0;
			indexEffects = 0;
			if (m_useCustomSampler) this->generate(Point2i());

			for (DATA_TYPE& v : values) {
				if (v == 0.0) {
					v += Epsilon;
				} else if (v == 1.0) {
					v -= Epsilon;
				}
			}
	}
    QuadratureSampler(QuadratureSampler *sampler) : Sampler(Properties()), RETURN_CONSTANT(CONSTANT_DATA), m_useCustomSampler(false) {}

    QuadratureSampler(Stream *stream, InstanceManager *manager) : Sampler(Properties()), RETURN_CONSTANT(CONSTANT_DATA), m_useCustomSampler(false) {}

    void serialize(Stream *stream, InstanceManager *manager) const {}

	virtual void advance() {
		//subSampler->advance();
	}

	virtual void generate(const Point2i& offset) {
		//subSampler->generate(offset);
		index = 0;
	}

	virtual void setSampleIndex(size_t sampleIndex) {
		//subSampler->setSampleIndex(sampleIndex);
		index = sampleIndex;
	}

	virtual mitsuba::Float next1D() {
		if (index < TOTAL_SAMPLES_VALID) {
			double res = values[index];
			index++;
			return res;
		} else {
			if (m_useCustomSampler) {
				//return subSampler->next1D();
				return uniform_dist(rng);
			}
			else {
				#ifdef ERROR_EXTRA_SAMPLE
				Log(EInfo, "Using sample outside samples given [1D]");
				#endif
				return RETURN_CONSTANT;
			}
		}
	}

	virtual Point2 next2D() {
		if (index < TOTAL_SAMPLES_VALID) {
			Point2 res;
			res.x = values[index];
			index++;
			if (index < TOTAL_SAMPLES_VALID) {
				res.y = values[index];
				index++;
			} else {
				if (m_useCustomSampler) {
					res.y = uniform_dist(rng);
				}
				else {
					#ifdef ERROR_EXTRA_SAMPLE
					Log(EInfo, "Using sample outside samples given [2D]");
					#endif
					res.y = RETURN_CONSTANT;
				}
			}
			return res;
		} else {
			if (m_useCustomSampler) {
				return Point2(uniform_dist(rng), uniform_dist(rng));
			}
			else {
				#ifdef ERROR_EXTRA_SAMPLE
				Log(EInfo, "Using sample outside samples given [2D]");
				#endif
				return Point2(RETURN_CONSTANT, RETURN_CONSTANT);
			}
		}
	}

	inline size_t getSampleCount() const {
		return TOTAL_SAMPLES_VALID;
	}

	inline size_t getSampleIndex() const {
		return index;
	}

	/*mitsuba::Float motionBlurSample() {
		if (indexEffects < (C_DIM_N - TOTAL_SAMPLES_VALID)) {
			DATA_TYPE res = values[TOTAL_SAMPLES_VALID + indexEffects];
			indexEffects++;
			return res;
		} else {
			return RETURN_CONSTANT;
		}
	}

	Point2 dofSample() {
		if (indexEffects < (C_DIM_N - TOTAL_SAMPLES_VALID)) {
			Point2 res;
			res.x = values[TOTAL_SAMPLES_VALID + indexEffects];
			indexEffects++;
			if (indexEffects < (C_DIM_N - TOTAL_SAMPLES_VALID)) {
				res.y = values[TOTAL_SAMPLES_VALID + indexEffects];
				indexEffects++;
			} else {
				res.y = RETURN_CONSTANT;
			}
			return res;
		} else {
			return Point2(RETURN_CONSTANT, RETURN_CONSTANT);
		}
	}*/

	//MTS_DECLARE_CLASS()
private:
	ARRAY_TYPE values;
	std::size_t index;
	std::size_t indexEffects;

	mitsuba::Float timeSample;
	Point2 dofSample;

	bool m_useCustomSampler;
	RNG rng;
	const double RETURN_CONSTANT;

	std::uniform_real_distribution<mitsuba::Float> uniform_dist;
};

MTS_NAMESPACE_END

#endif /* __QUADRATURE_SAMPLER_H */
