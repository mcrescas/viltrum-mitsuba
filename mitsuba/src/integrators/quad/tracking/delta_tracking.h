#pragma once

#include <array>
#include <random>

#include "../../../../ext/viltrum/quadrature/vector-dimensions.h"
#include "../../../../ext/viltrum/quadrature/range.h"

using namespace viltrum;

template<typename RNG>
class DeltaTracking {
    mutable RNG rng;

    double maj_sigma_t;
    FILE *m_f;

    template<typename Result>
    struct Samples {
        Result sumatory;
        unsigned long counter;
        unsigned long counter_queries;
        Samples() : sumatory(0),counter(0),counter_queries(0) { }
    };
public:
    void set_file(const char* name_file)
    {
        if(m_f)
            fclose(m_f);

        m_f = fopen(name_file, "w");
    }

    template<typename F, typename Float>
    auto init(const F& f, const Range<Float,1>& range) const {
        //First element of tuple is sum, second element is number of samples
        return Samples<decltype(f(range.min(0)))>();
    }

    template<typename F, typename Float, typename Result>
    void step(const F& f, const Range<Float,1>& range, Samples<Result>& samples, bool verbose = false) const {

        std::exponential_distribution<Float> dis(maj_sigma_t);
        std::uniform_real_distribution<Float> udis01(0,1);

        Float Tr = 1;
        Float t = range.min(0);

        if( m_f )
            fprintf(m_f,"%f %f", t, 1);

        while(true)
        {
            Float dt = dis(rng);
            Float t0 = t;
            t += dt;

            if( t>range.max(0) )
                break;


            ++ samples.counter_queries;

            if( verbose )
                printf("Step: %f [%f, %f] - mu_t(t)= %f / %f: T(t) = %f \n", t, range.min(0), range.max(0), f(t), maj_sigma_t, Tr);

            if( udis01(rng) < f(t)/maj_sigma_t)
            {
                Tr = 0;

                if( m_f )
                    fprintf(m_f,"%f %f\n%f %f\n", t, 1, t, 0);

                break;
            }

        }


        if( m_f )
        {
            fprintf(m_f,"%f %f\n", range.max(0), Tr);
            fprintf(m_f,"s\n");

        }

        samples.sumatory += Tr;
        ++samples.counter;
    }

    template<typename F, typename Result>
    Result integral(const F& f, const Samples<Result>& samples) const {
        return (samples.counter==0)?decltype(samples.sumatory)(0):(samples.sumatory/samples.counter);
    }

    template<typename F, typename Result>
    unsigned long queries(const F& f, const Samples<Result>& samples) const {
        return (samples.counter==0)?decltype(samples.counter_queries)(0):(samples.counter_queries);
    }

    DeltaTracking(RNG&& r, double majorant) :
        rng(std::forward<RNG>(r)), maj_sigma_t(majorant), m_f(0) { }

    ~DeltaTracking()
    {
        if(m_f)
            fclose(m_f);

        m_f=0;
    }

};


template<typename RNG>
auto delta_tracking(RNG&& rng, double majorant) {
    return DeltaTracking<RNG>(std::forward<RNG>(rng), majorant);
}

