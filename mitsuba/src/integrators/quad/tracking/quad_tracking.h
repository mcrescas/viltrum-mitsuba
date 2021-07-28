#pragma once

#include "../../../../ext/viltrum/quadrature/integrate-bins-adaptive.h"

using namespace viltrum;

template<typename Nested, typename Error, typename RNG>
class AdaptiveQuadTracking {
    mutable RNG rng;

    StepperAdaptive<Nested,Error> cv_stepper;
    unsigned long adaptive_iterations;

    template<typename R, typename Result>
    struct Data {
        std::vector<R> regions;

        Result Tc;
        Result sumatory;
        unsigned long counter;
        unsigned long counter_queries;

        Data( std::vector<R>&& rs, Result _Tc ) :
            regions(std::forward<std::vector<R>>(rs)),
            sumatory(0),counter(0),counter_queries(0), Tc(_Tc) { }
    };
public:

    template<typename F, typename Float>
    auto init(const F& f, const Range<Float,1>& range) const {

        auto regions = cv_stepper.init(f,range);
        for(unsigned int i=0; i<adaptive_iterations; ++i)
                cv_stepper.step(f,range,regions);

        using VECTOR_TYPE = typename decltype(regions)::value_type;

        return Data<VECTOR_TYPE, decltype(f(std::array<Float,1>{{range.min(0)}}))> (std::move(regions), exp(-cv_stepper.integral(f,range,regions)));
    }

    template<typename F, typename Float, typename R>
    void step(const F& f, const Range<Float,1>& range, Data<R,Float>& data, bool verbose = false) const
    {

    }

    template<typename F, typename Float, typename R>
    Float integral(const F& f, const Data<R,Float>& data) const {
        return data.Tc;
    }


    template<typename F, typename Float, typename R>
    unsigned long queries(const F& f, const Data<R,Float>& data) const {
        return (data.counter==0)?decltype(data.counter_queries)(0):(data.counter_queries);
    }


    AdaptiveQuadTracking(
        Nested&& nested, Error&& error, RNG&& r, unsigned long ai) :
            cv_stepper(std::forward<Nested>(nested), std::forward<Error>(error)), adaptive_iterations(ai),
            rng(std::forward<RNG>(r)) { }


    ~AdaptiveQuadTracking() {}
};

template<typename Nested, typename Error, typename RNG>
auto adaptive_quad_tracking(Nested&& nested, Error&& error, unsigned long adaptive_iterations,
                                      RNG&& rng ) {
    return AdaptiveQuadTracking<std::decay_t<Nested>,std::decay_t<Error>, std::decay_t<RNG>>(
        std::forward<Nested>(nested),std::forward<Error>(error), std::forward<RNG>(rng), adaptive_iterations);
}

template<typename Nested, typename RNG>
auto adaptive_quad_tracking(Nested&& nested, unsigned long adaptive_iterations, RNG&& rng)
{
    return adaptive_quad_tracking(std::forward<Nested>(nested), error_single_dimension_standard(), adaptive_iterations, std::forward<RNG>(rng));
}
