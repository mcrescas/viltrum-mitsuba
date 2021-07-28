
#ifndef SRC_INTEGRATORS_SINGLESCATTERING_SAMPLING_H_
#define SRC_INTEGRATORS_SINGLESCATTERING_SAMPLING_H_


#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN
class SingleScatteringSampling{
public:

    /*
    * Mean Free Path Sampling
    */

    //First value : s, second value factor (1/probability)
	static std::tuple<Float,Float> sample_distance_mfp(const Medium& medium, Float sample) {
        sample = std::min(1.0 - 1.e-20, std::max(1.e-20,1.0-sample));
		Float sigma_ref = medium.getSigmaT().max();
        Float d = - std::log(sample)/sigma_ref;
        return std::tuple<Float,Float>(d, 1.0f/(sigma_ref*std::exp(-d*sigma_ref)));
    }
    static Float pdf_mfp(const Medium& medium, Float t) {
		Float sigma_ref = medium.getSigmaT().max();
        Float sample = std::exp(-t*sigma_ref);

        return sample*sigma_ref;
    }
    static Float cdf_mfp(const Medium& medium, Float t) {
		Float sigma_ref = medium.getSigmaT().max();
        return 1.0 - std::exp(-t*sigma_ref);

    }



    /*
    * Equiangular Sampling
    */

    //First value : s, second value factor (1/probability)
	static std::tuple<Float,Float> sample_distance_equiangular(const Point lightPos, const Ray &ray, Float sample) {
        sample = std::min(1.0 - 1.e-20, std::max(1.e-20,sample));

		// distance from ray to light closest point in ray
        Float t_light = dot(lightPos - ray.o, ray.d);

        // distance of projection point from light
        float distance = std::max((lightPos - (ray.o + t_light*ray.d)).length(),1.e-20);

        Float thetaA = std::atan(-t_light/distance);
        Float thetaB = M_PI/2;

        Float r = math::lerp(sample, thetaA, thetaB);
        Float tan_r = std::tan(r);
        //Float r = std::atan(r);

        Float t = tan_r* distance;
        Float pdf = distance/((thetaB - thetaA)*(distance*distance + t*t));
        return std::tuple<Float,Float>(t_light+t, 1/pdf);
    }

    static Float pdf_equiangular(const Point lightPos, const Ray &ray, const Float dist){
        // distance from ray to light closest point in ray
        Float t_light = dot(lightPos - ray.o, ray.d);
        Float t = dist - t_light;


        // distance of point from light
        float distance = std::max((lightPos - (ray.o + t_light*ray.d)).length(),1.e-20);

        Float thetaA = std::atan(-t_light/distance);
        Float thetaB = M_PI/2;

        Float pdf = distance/((thetaB - thetaA)*(distance*distance + t*t));
        return pdf;
    }

    static Float cdf_equiangular(const Point lightPos, const Ray &ray, const Float dist){
        // distance from ray to light closest point in ray
        Float t_light = dot(lightPos - ray.o, ray.d);
        Float t = dist - t_light;

        // distance of point from light
        float distance = std::max((lightPos - (ray.o + t_light*ray.d)).length(),1.e-20);

        Float thetaA = std::atan(-t_light/distance);
        Float thetaB = M_PI/2;

        //return (std::atan(t/distance)/(thetaB - thetaA) -thetaA/(thetaB - thetaA)) * (thetaB - thetaA)/(-thetaA + M_PI/2);
        return 2*(thetaA - std::atan(t/distance)) / (2*thetaA - M_PI);
    }

};

MTS_NAMESPACE_END
#endif /* SRC_INTEGRATORS_SINGLESCATTERING_SAMPLING_H_ */
