#pragma once

#include <array>
#include <random>

#include "../../../../ext/viltrum/quadrature/vector-dimensions.h"
#include "../../../../ext/viltrum/quadrature/range.h"

using namespace viltrum;

template<typename RNG>
class ResidualRatioTracking {
    mutable RNG rng;

    FILE *m_f;

    double maj_sigma_t;
    double c_sigma_t;

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

        if(!samples.counter)
        {

            if(m_f)
            {
                Float dt = (range.max(0)-range.min(0))/(999.);
                for(Float t=range.min(0); t<= range.max(0); t+= dt)
                    fprintf(m_f,"%f %f\n", t, f(t));

                fprintf(m_f,"n\n");


                fprintf(m_f, "%f %f\n%f %f\n", range.min(0), c_sigma_t, range.max(0), c_sigma_t);
                fprintf(m_f,"n\n");

                fprintf(m_f,"%f %f\n%f %f\n", range.min(0), maj_sigma_t, range.max(0), maj_sigma_t);
                fprintf(m_f,"n\n");


                // Transmittance computed with ray marching (biased)
                Float delta_t = (range.max(0)-range.min(0))/(Float)10000;

                Float Tr = 1, Tc = 1, T = 1;
                auto t = range.min(0)+delta_t/2.;

                while(true)
                {

                    if( t>range.max(0) )
                        break;

                    T *= exp(-f(t)*delta_t);
                    Tc *= exp(-c_sigma_t*delta_t);
                    Tr *= exp(-(f(t)-c_sigma_t)*delta_t);

                    fprintf(m_f, "%f %f %f %f\n", t, T, Tc, Tr);

                    t += delta_t;
                }
                fprintf(m_f,"s\n");


            }
        }

        std::exponential_distribution<Float> dis(maj_sigma_t);

        Float Tc = exp(-(range.max(0)-range.min(0))*c_sigma_t);

        Float Tr = 1;
        Float t = range.min(0);

        if( m_f )
        {
            fprintf(m_f,"%f %f\n", t, Tr);
        }

        while(true)
        {
            Float dt = dis(rng);
            t += dt;

            if( t>range.max(0) )
                break;

            if( m_f )
            {
                fprintf(m_f,"%f %f\n", t, Tr);
            }

            Tr *= (1 - (f(t) - c_sigma_t)/maj_sigma_t);

            if( m_f )
            {
                fprintf(m_f,"%f %f\n", t, Tr);
            }


            if( verbose )
                printf("Step: %f [%f, %f] - mu_t(t)= %f : T(t) = %f; Tc(t) = %f; T_c = %f\n", t, range.min(0), range.max(0), (f(t)-c_sigma_t)/maj_sigma_t, Tr, exp(-c_sigma_t*t), Tc);

            ++ samples.counter_queries;
        }

        if( m_f )
        {
            fprintf(m_f,"%f %f\n", range.max(0), Tr);
            fprintf(m_f,"s\n");

        }

        samples.sumatory += Tr*Tc;
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

    ResidualRatioTracking(RNG&& r, double majorant, double control) :
        rng(std::forward<RNG>(r)), maj_sigma_t(majorant), c_sigma_t(control), m_f(0) { }

    ~ResidualRatioTracking()
    {
        if(m_f)
            fclose(m_f);

        m_f=0;
    }
};


template<typename RNG>
auto ratio_tracking(RNG&& rng, double majorant) {
    return ResidualRatioTracking<RNG>(std::forward<RNG>(rng), majorant, 0.);
}

template<typename RNG>
auto residual_ratio_tracking(RNG&& rng, double majorant, double control) {
    return ResidualRatioTracking<RNG>(std::forward<RNG>(rng), majorant, control);
}
